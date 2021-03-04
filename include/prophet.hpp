#ifndef PROPHET_HPP
#define PROPHET_HPP

#include "table.hpp"
#include "prophet_model.hpp"
#include "prophet_var_context.hpp"
#include "vecops.hpp"
#include "time_delta.hpp"
#include "funcs.hpp"

#include <stan/services/optimize/lbfgs.hpp>
#include <stan/callbacks/interrupt.hpp>
#include <stan/callbacks/stream_logger.hpp>
#include <stan/callbacks/stream_writer.hpp>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/random_access_index.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/member.hpp>

#define pred(expr) [](double e) { return expr; }

std::vector<size_t> read_dims(std::string line) {

    std::vector<size_t> dims{};
    std::string delimiter = ",";

    size_t pos = 0;
    std::string token;
    while ((pos = line.find(delimiter)) != std::string::npos) {
        token = line.substr(0, pos);
        dims.push_back(std::stoi(token));
        line.erase(0, pos + delimiter.length());
    }

    return dims;
}

template <typename T>
std::pair<std::vector<T>, std::vector<size_t> > load(std::string fname) {

    std::cout << "== " << fname << ":" << std::endl;

    std::ifstream infile("../data/" + fname + ".txt");
    std::string line;
    std::getline(infile, line);
    auto dims = read_dims(line);

    std::vector<T> values{};

    T val;
    while (infile >> val) {
        std::cout << val << std::endl;
        values.push_back(val);
    }

    return std::make_pair(values, dims);
}

struct prophet_interrupt : public stan::callbacks::interrupt {
    void operator()() {
    }
};

class value : public stan::callbacks::writer {
    std::vector<double> x_;

public:
    value() {}

    using stan::callbacks::writer::operator();

    void operator()(const std::vector<double>& x) {
        x_ = x;
    }

    const std::vector<double> x() const {
        return x_;
    }
};

namespace prophet {
    struct seasonality {
        std::string name;
        double period;
        int fourier_order;
        double prior_scale;
        std::string mode;
        std::string condition_name;
    };

    struct fit_params {
        double k;
        double m;
        std::vector<double> delta;
        std::vector<double> beta;
    };

    typedef boost::multi_index::multi_index_container<
        seasonality,
        boost::multi_index::indexed_by<
            boost::multi_index::random_access<>,
            boost::multi_index::hashed_unique<boost::multi_index::member<seasonality, std::string, &seasonality::name>>
        >
    > seasonalities_container;

    class prophet {

        tbl::table history;
        std::vector<double> history_dates;
        std::string growth{"linear"};
        double start = 0.0;
        double y_scale = 0.0;
        bool logistic_floor = false;
        double t_scale = 0.0;
        seasonalities_container seasonalities;
        double seasonality_prior_scale;
        std::string seasonality_mode;

        std::string yearly_seasonality;
        std::string weekly_seasonality;
        std::string daily_seasonality;

        std::string train_holiday_names;

        tbl::table train_component_cols;
        std::map<std::string, std::vector<std::string>> component_modes;

        double changepoint_prior_scale = 0.05;
        std::vector<double> changepoints;
        std::vector<double> changepoints_t;
        double changepoint_range;
        int n_changepoints;

        fit_params params;

        int uncertainty_samples = 1000;
        double interval_width = 0.80;


        public:
            explicit prophet()
            : seasonality_prior_scale(10.0),
              seasonality_mode("additive"),
              yearly_seasonality("auto"),
              weekly_seasonality("auto"),
              daily_seasonality("auto"),
              train_holiday_names(""),
              changepoint_range(0.8),
              n_changepoints(25)
            {};
            prophet(prophet const&) = default;
            prophet(prophet&&) = default;
            ~prophet() = default;

            tbl::table setup_table(tbl::table& tbl, bool scales = false) {
                // TODO: Checks on time columns and values (e.g. nan, missing, etc)
                // TODO: extra regressors
                // TODO: seasonalities
                // TODO: initialize_scales

                std::cout << "==> setup_table:" << std::endl;
                auto [tbl_rows, tbl_cols] = tbl.shape();
                std::cout << "rows: " << tbl_rows << ", cols: " << tbl_cols << std::endl;

                initialize_scales(tbl, scales);

                auto rows = tbl.shape().first;
                if (rows == 0) {
                    rows = tbl.get_times().size();
                }
                // TODO: logistic floor
                tbl.push_col("floor", std::vector<double>(rows, 0.0));

                // df['t'] = (df['ds'] - self.start) / self.t_scale
                auto ds = tbl.get_times();
                auto t = (ds - start) / t_scale;
                tbl.push_col("t", t);

                // if 'y' in df:
                //     df['y_scaled'] = (df['y'] - df['floor']) / self.y_scale
                if (tbl.exists("y")) {
                    auto y = tbl.get_col("y");
                    auto floor = tbl.get_col("floor");
                    auto y_scaled = (y - floor) / y_scale;
                    tbl.push_col("y_scaled", y_scaled);
                }

                // TODO: Handle extra regressors

                std::cout << "<== setup_table" << std::endl;

                return tbl;
            }

            void initialize_scales(tbl::table& tbl, bool scales) {
                if (!scales) {
                    return;
                }

                // TODO: Handle logistic growth
                auto floor = 0.0;
                y_scale = (tbl.get_col("y") - floor) >> vec::abs() >> vec::max();
                std::cout << "y_scale: " << y_scale << std::endl;
                if (y_scale == 0) {
                    y_scale = 1;
                }

                auto times = tbl.get_times();
                start = times >> vec::min();
                t_scale = (times >> vec::max()) - start;

                // TODO: handle extra regressors

            }

            void set_auto_seasonalities() {
                auto first = history.get_times() >> vec::min();
                auto last = history.get_times() >> vec::max();
                auto yearly_disable = last - first < time_delta<days>(730);

                int fourier_order;

                fourier_order = parse_seasonality_args("yearly", yearly_seasonality, yearly_disable, 10);
                if (fourier_order > 0) {
                    seasonality s;
                    s.name = "yearly";
                    s.period = 365.25;
                    s.fourier_order = fourier_order;
                    s.prior_scale = seasonality_prior_scale;
                    s.mode = seasonality_mode;
                    s.condition_name = "";
                    seasonalities.push_back(s);
                }

                auto min_dt = history.get_times()
                    >> vec::diff()
                    >> vec::filter(pred( !std::isnan(e) ))
                    >> vec::min();

                auto weekly_disable = (
                    (last - first < time_delta<weeks>(2)) ||
                    (min_dt >= time_delta<weeks>(1)));
                fourier_order = parse_seasonality_args("weekly", weekly_seasonality, weekly_disable, 3);
                if (fourier_order > 0) {
                    seasonality s;
                    s.name = "weekly";
                    s.period = 7;
                    s.fourier_order = fourier_order;
                    s.prior_scale = seasonality_prior_scale;
                    s.mode = seasonality_mode;
                    s.condition_name = "";
                    seasonalities.push_back(s);
                }

                auto daily_disable = (
                    (last - first < time_delta<days>(2)) ||
                    (min_dt >= time_delta<days>(1)));
                fourier_order = parse_seasonality_args("daily", daily_seasonality, daily_disable, 4);
                if (fourier_order > 0) {
                    seasonality s;
                    s.name = "daily";
                    s.period = 1;
                    s.fourier_order = fourier_order;
                    s.prior_scale = seasonality_prior_scale;
                    s.mode = seasonality_mode;
                    s.condition_name = "";
                    seasonalities.push_back(s);
                }
            }

            int parse_seasonality_args(const std::string& name, const std::string& arg, bool auto_disable, int default_order) {
                int fourier_order = 0;
                auto &hash_index = seasonalities.get<1>();

                if (arg == "auto") {
                    if (hash_index.find(name) != hash_index.end()) {
                        std::cout << "Found custom seasonality named '"
                            << name << "', disabling built-in "
                            << name <<" seasonality." << "\n";
                    } else if (auto_disable) {
                        std::cout << "Disabling '"
                            << name << "' seasonality. Run prophet with "
                            << name << "_seasonality=True to override this.\n";
                    } else {
                        fourier_order = default_order;
                    }
                    
                } else if (arg == "true") {
                    fourier_order = default_order;
                } else if (arg == "false") {
                    fourier_order = 0;
                } else {
                    fourier_order = std::stoi(arg);
                }

                return fourier_order;
            }

            std::tuple<
                tbl::table,
                std::vector<double>,
                tbl::table,
                std::map<std::string, std::vector<std::string>>
            > make_all_seasonality_features(const tbl::table& tbl) {

                std::vector<tbl::table> seasonal_features_vec;
                std::vector<double> prior_scales;
                std::map<std::string, std::vector<std::string>> modes = {
                    {"additive", std::vector<std::string>{}},
                    {"multiplicative", std::vector<std::string>{}}
                };

                for (auto& props: seasonalities) {
                    auto features = make_seasonality_features(
                        tbl.get_times(),
                        props.period,
                        props.fourier_order,
                        props.name);
                    if (props.condition_name != "") {
                        // TODO
                        std::cerr << "make_all_seasonality_features: condition_name is not implemented" << std::endl;
                        std::abort();
                    }
                    seasonal_features_vec.push_back(features);
                    // prior_scales.extend(
                    //     [props['prior_scale']] * features.shape[1])
                    for (auto i = 0; i < features.shape().second; ++i) {
                        prior_scales.push_back(props.prior_scale);
                    }
                    // modes[props['mode']].append(name)
                    modes[props.mode].push_back(props.name);
                }

                std::cout << "==> prior_scales:\n";
                for (auto& v: prior_scales) {
                    std::cout << v << std::endl;
                }
                std::cout << "<== prior_scales" << std::endl;

                std::cout << "==> modes\n";
                for (auto& [mode, vec]: modes) {
                    std::cout << "\t" << mode << "\n";
                    for (auto& v: vec) {
                        std::cout << "\t\t" << v << "\n";
                    }
                }
                std::cout << "<== modes" << std::endl;

                // TODO Holiday features
                // TODO Addtional regressors
                // TODO Dummy to prevent empty X

                // seasonal_features = pd.concat(seasonal_features, axis=1)
                tbl::table seasonal_features;
                for (auto& t: seasonal_features_vec) {
                    seasonal_features += t;
                }
                std::cout << "==> make_all_seasonality_features, seasonal_features.columns" << std::endl;
                for (auto& name: seasonal_features.get_names()) {
                    std::cout << name << std::endl;
                }
                std::cout << "<== make_all_seasonality_features, seasonal_features.columns" << std::endl;
                // component_cols, modes = self.regressor_column_matrix(
                //     seasonal_features, modes
                // )
                auto [component_cols, enhanced_modes] = regressor_column_matrix(seasonal_features, modes);
                // return seasonal_features, prior_scales, component_cols, modes

                return std::tuple(seasonal_features, prior_scales, component_cols, enhanced_modes);
            }

            std::pair<tbl::table, std::map<std::string, std::vector<std::string>>> regressor_column_matrix(const tbl::table& seasonal_features, std::map<std::string, std::vector<std::string>> modes) {

                std::vector<std::pair<int, std::string>> components;

                int i = 0;
                for (auto& name: seasonal_features.get_names()) {
                    components.push_back(std::make_pair(i++, first(name, "_delim_")));
                }

                if (train_holiday_names != "") {
                    // TODO
                    std::cerr << "Holiday names are not yet supported" << std::endl;
                    std::abort();
                }

                std::vector<std::string> modes_{"additive", "multiplicative"};
                for (auto& mode: modes_) {
                    // components = self.add_group_component(
                    //     components, mode + '_terms', modes[mode]
                    // )

                    std::vector<std::string> groups = modes[mode];
                    components = add_group_component(components, mode + "_terms", groups);

                    // # Add combination components to modes
                    // modes[mode].append(mode + '_terms')
                    // modes[mode].append('extra_regressors_' + mode)
                    auto& m = modes[mode];
                    m.push_back(mode + "_terms");
                    m.push_back("extra_regressors_" + mode);
                }

                // modes[self.seasonality_mode].append('holidays')
                auto& _m = modes[seasonality_mode];
                _m.push_back("holidays");

                for (auto& [col, component]: components) {
                    std::cout << "col: " << col << ", component: " << component << "\n";
                }

                // Create cross tab cols to component
                auto max = std::max_element(begin(components), end(components),
                    [](const std::pair<int, std::string>& left, const std::pair<int, std::string>& right){
                    return left.first <  right.first;
                }) -> first;
                std::cout << "max col: " << max << std::endl;

                // Get the unique components
                std::vector<std::string> individual_components;
                for (auto& [col, component]: components) {
                    individual_components.push_back(component);
                }

                std::sort(individual_components.begin(), individual_components.end());
                std::vector<std::string>::iterator it;
                it = std::unique(individual_components.begin(), individual_components.end());
                individual_components.resize(std::distance(individual_components.begin(),it));

                std::cout << "individual_components\n";
                for (auto& component: individual_components) {
                    std::cout << component << "\n";
                }
                std::cout << "individual_components done" << std::endl;

                std::map<std::string, std::vector<double>> m;
                for (auto& name: individual_components) {
                    m[name] = std::vector<double>(max + 1, 0.0);
                }

                for (auto& [col, component]: components) {
                    auto ccol = m[component];
                    ccol[col] = 1;
                    m[component] = ccol;
                }

                tbl::table component_cols;
                for (auto& [name, vec]: m) {
                    component_cols.push_col(name, vec);
                }

                // end crosstab

                for (auto& mode: modes_) {
                    if (!component_cols.exists(mode + "_terms")) {
                        component_cols.push_col(mode + "_terms", std::vector<double>(max + 1, 0.0));
                    }
                }

                std::cout << "component_cols\n";
                component_cols.print();
                std::cout << "component_cols done.\n";
                std::cout << "modes\n";
                for (auto& [m, ms]: modes) {
                    std::cout << "\t" << m << "\n";
                    for (auto& blah: ms) {
                        std::cout << "\t\t" << blah << "\n";
                    }
                }
                std::cout << "modes done.\n";

                return std::make_pair(component_cols, modes);
            }

            std::vector<std::pair<int, std::string>> add_group_component(
                std::vector<std::pair<int, std::string>> components,
                std::string name,
                std::vector<std::string> group) {

                for (auto& elem: group) {
                    std::cout << elem << "\n";
                }

                std::vector<int> group_cols;
                for (auto& [col, component]: components) {
                    if(std::find(group.begin(), group.end(), component) != group.end()) {
                        group_cols.push_back(col);
                    }
                }
                // Get unique cols
                std::sort(group_cols.begin(), group_cols.end());
                std::vector<int>::iterator it;
                it = std::unique(group_cols.begin(), group_cols.end());
                group_cols.resize(std::distance(group_cols.begin(), it));
                if (group_cols.size() > 0) {
                    for (auto& col: group_cols) {
                        components.push_back({col, name});
                    }
                }

                return components;
            }

            void set_changepoints() {

                if (false) {
                    // TODO Handle self.changepoints is not None
                } else {
                    int hist_size = int(floor(history.shape().first * changepoint_range));

                    if (n_changepoints + 1 > hist_size) {
                        std::cerr << "Not implemented yet" << std::endl;
                        std::abort();
                    }

                    std::cout << "hist_size: " << hist_size << ", n_changepoints: " << n_changepoints << "\n";

                    if (n_changepoints > 0) {
                        auto cp_indexes = linspace(0, hist_size - 1, n_changepoints + 1)
                            >> vec::round()
                            >> vec::to_int();
                        std::cout << "cp_indexes:\n";
                        for (auto& cp: cp_indexes) {
                            std::cout << cp << std::endl;
                        }

                        for (auto& index: cp_indexes) {
                            changepoints.push_back(history.get_times()[index]);
                        }
                        changepoints.erase(changepoints.begin());
                    }
                }

                if (changepoints.size() > 0) {
                    changepoints_t = (changepoints - start) / t_scale;
                    std::sort(changepoints_t.begin(), changepoints_t.end());
                    std::cout << "changepoints_t:\n";
                    for (auto& v: changepoints_t) {
                        std::cout << v << "\n";
                    }
                    std::cout << "changepoints_t done." << std::endl;
                }
            }

            auto linear_growth_init(const tbl::table& t) {
                size_t i0 = 0;
                size_t i1 = t.get_times().size() - 1;
                auto T = t.get_col("t")[i1] - t.get_col("t")[i0];
                auto k = (t.get_col("y_scaled")[i1] - t.get_col("y_scaled")[i0]) / T;
                auto m = t.get_col("y_scaled")[i0] - k * t.get_col("t")[i0];

                return std::pair(k, m);
            }

            prophet fit(const tbl::table& tbl) {
                // TODO: Re-enable const on tbl
                // TODO: See if self.history exists, bail if it does
                // TODO: See if ds, y exist in tbl, bail if they do not
                // TODO: Copy tbl to history so we can modify history without worry about tbl
                history = tbl;

                history_dates = tbl.get_times();
                std::sort(history_dates.begin(), history_dates.end());
                std::vector<double>::iterator it;
                it = std::unique(history_dates.begin(), history_dates.end());
                history_dates.resize(std::distance(history_dates.begin(),it));

                history = setup_table(history, true);
                set_auto_seasonalities();
                auto [seasonal_features, prior_scales, component_cols, modes] = make_all_seasonality_features(history);

                train_component_cols = component_cols;
                component_modes = modes;

                std::cout << "seasonal_features shape, rows: " << seasonal_features.shape().first << ", cols: " << seasonal_features.shape().second << std::endl;

                set_changepoints();

                std::map<std::string, int> trend_indicator{{"linear", 0}, {"logistic", 1}, {"flat", 2}};

                std::map<std::string,
                    std::pair<std::vector<double>, std::vector<size_t> > > vars_r{};

                auto [history_rows, history_cols] = history.shape();

                std::cout << "vars_r" << std::endl;
                std::vector<double> X;
                for (auto& name: seasonal_features.get_names()) {
                    auto const &col = seasonal_features.get_col(name);
                    for (auto& v: col) {
                        X.push_back(v);
                    }
                }
                auto [seasonal_rows, seasonal_cols] = seasonal_features.shape();
                vars_r["X"] = std::pair(X, std::vector<size_t>{seasonal_rows, seasonal_cols});
                vars_r["sigmas"] = std::pair(prior_scales, std::vector<size_t>{prior_scales.size()});
                auto t = history.get_col("t");
                vars_r["t"] = std::pair(t, std::vector<size_t>{t.size()});
                vars_r["t_change"] = std::pair(changepoints_t, std::vector<size_t>{changepoints_t.size()});
                vars_r["tau"] = std::pair(std::vector<double>{changepoint_prior_scale}, std::vector<size_t>{0});
                auto y = history.get_col("y_scaled");
                vars_r["y"] = std::pair(y, std::vector<size_t>{y.size()});

                double k, m;
                if (growth == "linear") {
                    std::cout << "history_rows: " << history_rows << std::endl;
                    auto cap = std::vector<double>(history_rows, 0);
                    vars_r["cap"] = std::pair(cap, std::vector<size_t>{cap.size()});
                    std::cout << "cap size: " << cap.size() << std::endl;
                    auto [ ktmp, mtmp ] = linear_growth_init(history);
                    k = ktmp;
                    m = mtmp;
                } else if (growth == "flat") {

                } else {

                }

                std::cout << "k: " << k << ", m: " << m << std::endl;

                std::map<std::string,
                    std::pair<std::vector<int>, std::vector<size_t> > > vars_i{};

                auto [tmp1, K] = seasonal_features.shape();
                vars_i["K"] = std::pair(std::vector<int>{static_cast<int>(K)}, std::vector<size_t>{0});
                auto len = changepoints_t.size();
                vars_i["S"] = std::pair(std::vector<int>{static_cast<int>(len)}, std::vector<size_t>{0});
                vars_i["T"] = std::pair(std::vector<int>{static_cast<int>(history_rows)}, std::vector<size_t>{0});
                auto s_a = component_cols.get_col("additive_terms") >> vec::to_int();
                vars_i["s_a"] = std::pair(s_a, std::vector<size_t>{s_a.size()});
                auto s_m = component_cols.get_col("multiplicative_terms") >> vec::to_int();
                vars_i["s_m"] = std::pair(s_m, std::vector<size_t>{s_m.size()});
                vars_i["trend_indicator"] = std::pair(std::vector<int>{trend_indicator[growth]}, std::vector<size_t>{0});

                stan::io::prophet_var_context context(vars_r, vars_i);

                prophet_model_namespace::prophet_model prophet_model(context);
                std::vector<std::string> constrained_param_names;
                prophet_model.constrained_param_names(constrained_param_names);

                std::map<std::string,
                    std::pair<std::vector<double>, std::vector<size_t> > > vars_r_optimize{};

                std::vector<double> beta(K, 0);
                vars_r_optimize["beta"] = std::pair(beta, std::vector<size_t>{beta.size()});
                std::vector<double> delta(changepoints_t.size(), 0);
                vars_r_optimize["delta"] = std::pair(delta, std::vector<size_t>{delta.size()});
                vars_r_optimize["k"] = std::pair(std::vector<double>{k}, std::vector<size_t>{0});
                vars_r_optimize["m"] = std::pair(std::vector<double>{m}, std::vector<size_t>{0});

                std::map<std::string,
                    std::pair<std::vector<int>, std::vector<size_t> > > vars_i_optimize{};

                vars_i_optimize["sigma_obs"] = std::pair(std::vector<int>{1}, std::vector<size_t>{0});

                stan::io::prophet_var_context init_context(vars_r_optimize, vars_i_optimize);

                auto history_size = 5;
                auto init_alpha = 0.001;
                auto tol_obj = 1e-12;
                auto tol_rel_obj = 10000;
                auto tol_grad = 1e-08;
                auto tol_rel_grad = 1e+07;
                auto tol_param = 1e-08;
                auto refresh = 100;
                auto random_seed = 1854378400;
                auto id = 1;
                auto init_radius = 2;
                auto num_iterations = 10000;
                auto save_iterations = 1;

                prophet_interrupt interrupt{};
                stan::callbacks::stream_logger logger(std::cout, std::cout, std::cout,
                                                    std::cerr, std::cerr);
                stan::callbacks::writer init_writer;
                // stan::callbacks::writer sample_writer;
                value sample_writer;

                std::cout << "OPTIMIZING" << std::endl;

                stan::services::optimize::lbfgs(prophet_model,
                                                init_context,
                                                random_seed, id, init_radius,
                                                history_size,
                                                init_alpha,
                                                tol_obj,
                                                tol_rel_obj,
                                                tol_grad,
                                                tol_rel_grad,
                                                tol_param,
                                                num_iterations,
                                                save_iterations,
                                                refresh,
                                                interrupt, logger,
                                                init_writer, sample_writer);

                // double k;
                // double m;
                // std::vector<double> delta;
                delta.clear();
                double sigma_obs;
                // std::vector<double> beta;
                beta.clear();
                std::vector<double> trend;

                std::vector<std::string> param_names;
                std::vector<std::vector<size_t> > param_dims;
                prophet_model.get_param_names(param_names);
                prophet_model.get_dims(param_dims);

                std::vector<double> params = sample_writer.x();
                double lp = params.front();
                params.erase(params.begin());

                auto iter = params.begin();

                auto dims_iter = param_dims.begin();
                for (auto param_name : param_names) {
                    std::cout << "param_name: " << param_name << ", dims: ";
                    size_t param_dim = 1;
                    for (auto dim : *dims_iter) {
                        std::cout << param_dim << ", ";
                        param_dim = dim;
                    }
                    if (param_name == "k") {
                        k = *iter;
                        ++iter;
                    } else if (param_name == "m") {
                        m = *iter;
                        ++iter;
                    } else if (param_name == "delta") {
                        for (auto i = 0; i < param_dim; ++i) {
                            delta.push_back(*iter);
                            ++iter;
                        }
                    } else if (param_name == "sigma_obs") {
                        sigma_obs = *iter;
                        ++iter;
                    } else if (param_name == "beta") {
                        for (auto i = 0; i < param_dim; ++i) {
                            beta.push_back(*iter);
                            ++iter;
                        }
                    } else if (param_name == "trend") {
                        for (auto i = 0; i < param_dim; ++i) {
                            trend.push_back(*iter);
                            ++iter;
                        }
                    }
                    ++dims_iter;
                    std::cout << std::endl;
                }

                std::cout << "k: " << k << "\n";
                std::cout << "m: " << m << "\n";
                std::cout << "delta:\n";
                for (auto v: delta) {
                    std::cout << "\t" << v << "\n";
                }
                std::cout << "sigma_obs: " << sigma_obs << "\n";
                std::cout << "beta:\n";
                for (auto v: beta) {
                    std::cout << "\t" << v << "\n";
                }
                std::cout << "trend:\n";
                for (auto v: trend) {
                    std::cout << "\t" << v << "\n";
                }

                this->params.delta = delta;
                this->params.beta = beta;
                this->params.k = k;
                this->params.m = m;

                return *this;
            }

            tbl::table make_future_dataframe(int periods, const std::string& freq="D", bool include_history=true) {

                if (!(history_dates.size() > 0)) {
                    std::cerr << "Model has not been fit." << std::endl;
                    // TODO: Determine error handling, use exceptions, but need definition
                    std::abort();
                }

                auto last_date = history_dates >> vec::max();

                auto dates = date_range(last_date, periods + 1, freq);

                dates = dates >> vec::filter([last_date](double e) { return e > last_date; });

                tbl::table tbl;

                if (include_history) {
                    for (auto& d: history_dates) {
                        tbl.push_time(d);
                    }
                }

                for (auto& d: dates) {
                    tbl.push_time(d);
                }

                return tbl;
            }

            std::vector<double> predict_trend(const tbl::table& tbl) {
                const auto& t = tbl.get_col("t");
                std::vector<double> trend;
                if (growth == "linear") {
                    trend = piecewise_linear(t, params.delta, params.k, params.m, changepoints_t);
                } else if (growth == "logistic") {
                    std::cerr << "logistic growth not yet supported" << std::endl;
                    std::abort();
                } else if (growth == "flat") {
                    std::cerr << "flat growth not yet supported" << std::endl;
                    std::abort();
                }

                return trend * y_scale + tbl.get_col("floor");
            }

            tbl::table predict_seasonal_components(const tbl::table& tbl) {
                std::cout << "==> predict_seasonal_components" << std::endl;
                auto [seasonal_features, prior_scales, component_cols, modes] = make_all_seasonality_features(tbl);

                std::cout << "seasonal_features shape, rows: " << seasonal_features.shape().first << ", cols: " << seasonal_features.shape().second << std::endl;

                double lower_p = std::numeric_limits<double>::max();
                double upper_p = std::numeric_limits<double>::max();
                if (uncertainty_samples != std::numeric_limits<int>::max()) {
                    lower_p = 100 * (1.0 - interval_width) / 2;
                    upper_p = 100 * (1.0 + interval_width) / 2;

                    std::cout << "lower_p: " << lower_p << std::endl;
                    std::cout << "upper_p: " << upper_p << std::endl;
                }

                std::cout << "seasonal_features.to_eigen()" << std::endl;
                auto X = seasonal_features.to_eigen();
                std::cout << "X.rows(): " << X.rows() << ", X.cols(): " << X.cols() << std::endl;
                std::cout << "seasonal_features.to_eigen() done." << std::endl;
                tbl::table data;
                for (auto& component: component_cols.get_names()) {
                    auto beta_c = params.beta * component_cols.get_col(component);
                    Eigen::VectorXd v2 = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(beta_c.data(), beta_c.size());
                    std::cout << "v2.rows(): " << v2.rows() << ", cols: " << v2.cols() << std::endl;
                    auto comp = X * v2;
                    std::cout << "beta_c:\n";
                    for (auto& v: beta_c) {
                        std::cout << v << std::endl;
                    }
                    std::cout << "==> comp:" << std::endl;
                    std::cout << comp << std::endl;
                    std::cout << "<== comp:" << std::endl;
                    std::cout << "comp rows: " << comp.rows() << ", cols: " << comp.cols() << std::endl;

                    std::vector<double> comp_v;

                    for (auto i = 0; i < comp.rows(); ++i) {
                        for (auto j = 0; j < comp.cols(); ++j) {
                            comp_v.push_back(comp(i, j));
                        }
                    }

                    auto& component_mode = component_modes["additive"];
                    if (std::find(component_mode.begin(), component_mode.end(), component) != component_mode.end()) {
                        comp_v *= y_scale;
                    }

                    // TODO prophet is doing something dubious here. comp is an array
                    // of 1 element arrays and it takes the mean of each, effectively
                    // returning a single dimensional array of what's in each array.
                    // we're skipping that for now.

                    data.push_col(component, comp_v);
                    if (uncertainty_samples != std::numeric_limits<int>::max()) {
                        // TODO Same thing applies here as above, the python code
                        // is taking the percentile of 1 element arrays
                        data.push_col(component + "_lower", comp_v);
                        data.push_col(component + "_upper", comp_v);
                    }
                }
                std::cout << "<== predict_seasonal_components" << std::endl;

                return data;
            }

            std::vector<double> sample_predictive_trend(tbl::table tbl, int iteration) {
                return std::vector<double>{};
            }

            tbl::table sample_model(tbl::table tbl, tbl::table seasonal_features, int iteration, std::vector<double> s_a, std::vector<double> s_m) {
                return tbl::table{};
            }

            std::map<std::string, std::vector<double>> sample_posterior_predictive(tbl::table tbl) {

                // TODO In Python, the 'k' param is a 1x1 shape. But in the
                // stan model it's declared as a real scalar. So we're just going
                // to declare n_iterations = 1 here
                // n_iterations = self.params['k'].shape[0]

                std::cout << "==> sample_posterior_predictive\n";

                int n_iterations = 1;
                int samp_per_iter = fmax(1, ceil(uncertainty_samples / double(n_iterations)));

                std::cout << "samp_per_iter: " << samp_per_iter << "\n";

                auto [seasonal_features, prior_scales, component_cols, modes] = make_all_seasonality_features(tbl);

                std::map<std::string, std::vector<double>> sim_values = {{"yhat", {}}, {"trend", {}}};
                for (auto i = 0; i < n_iterations; ++i) {
                    for (auto j = 0; j < samp_per_iter; ++j) {
                        auto sim = sample_model(
                            tbl,
                            seasonal_features,
                            i,
                            component_cols.get_col("additive_terms"),
                            component_cols.get_col("multiplicative_terms"));
                    }
                }

                std::cout << "<== sample_posterior_predictive" << std::endl;

                return std::map<std::string, std::vector<double>>{};
            }

            tbl::table predict_uncertainty(tbl::table tbl) {
                auto sim_values = sample_posterior_predictive(tbl);

                return tbl::table{};
            }

            tbl::table predict(tbl::table tbl=tbl::table{}) {

                if (!(history.shape().second > 0)) {
                    std::cerr << "Model has not been fit." << std::endl;
                    // TODO: Determine error handling, use exceptions, but need definition
                    std::abort();
                }

                std::cout << "tbl.get_times().size(): " << tbl.get_times().size() << std::endl;

                if (!(tbl.shape().second > 0) && !(tbl.get_times().size() > 0)) {
                    std::cout << "predict, tbl = history" << std::endl;
                    tbl.print();
                    tbl = history;
                } else {
                    if (tbl.shape().first == 0 && tbl.get_times().size() == 0) {
                        std::cerr << "table has no rows" << std::endl;
                        std::abort();
                    }
                    tbl = setup_table(tbl);
                }

                auto [rows, cols] = tbl.shape();

                std::cout << "==> predict, tbl.shape():\n";
                std::cout << "rows: " << rows << ", cols: " << cols << "\n";
                std::cout << "<== predict, tbl.shape()" << std::endl;

                tbl.push_col("trend", predict_trend(tbl));
                std::cout << "==> trend\n";
                for (auto& v: tbl.get_col("trend")) {
                    std::cout << v << "\n";
                }
                std::cout << "<== trend" << std::endl;

                auto seasonal_components = predict_seasonal_components(tbl);
                std::cout << "==> seasonal_components:\n";
                seasonal_components.print();
                std::cout << "<== seasonal_components" << std::endl;

                if (uncertainty_samples != std::numeric_limits<int>::max()) {
                    auto intervals = predict_uncertainty(tbl);
                }

                return tbl;
            }
    };
}

#endif // PROPHET_HPP
