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

    std::ifstream infile("data/" + fname + ".txt");
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
    class model {
        public:
            explicit model() {};
            model(model const&) = default;
            model(model&&) = default;
            ~model() = default;

    };

    struct seasonality {
        std::string name;
        double period;
        int fourier_order;
        double prior_scale;
        std::string mode;
        std::string condition_name;
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
        double start = 0.0;
        double y_scale = 0.0;
        bool logistic_floor = false;
        double t_scale = 0.0;
        seasonalities_container seasonalities;
        // std::map<std::string, seasonality> seasonalities;
        double seasonality_prior_scale;
        std::string seasonality_mode;

        std::string yearly_seasonality;
        std::string weekly_seasonality;
        std::string daily_seasonality;

        std::string train_holiday_names;


        public:
            explicit prophet()
            : seasonality_prior_scale(10.0),
              seasonality_mode("additive"),
              yearly_seasonality("auto"),
              weekly_seasonality("auto"),
              daily_seasonality("auto"),
              train_holiday_names("")
            {};
            prophet(prophet const&) = default;
            prophet(prophet&&) = default;
            ~prophet() = default;

            tbl::table setup_table(tbl::table& tbl, bool scales = false) {
                // TODO: Checks on time columns and values (e.g. nan, missing, etc)
                // TODO: extra regressors
                // TODO: seasonalities
                // TODO: initialize_scales

                initialize_scales(tbl, scales);

                auto rows = tbl.shape().first;
                // TODO: logistic floor
                tbl.push_col("floor", std::vector<double>(rows, 0.0));

                // df['t'] = (df['ds'] - self.start) / self.t_scale
                auto ds = tbl.get_times();
                auto t = (ds - start) / t_scale;
                tbl.push_col("t", t);

                // if 'y' in df:
                //     df['y_scaled'] = (df['y'] - df['floor']) / self.y_scale
                auto y = tbl.get_col("y");
                auto floor = tbl.get_col("floor");
                auto y_scaled = (y - floor) / y_scale;
                tbl.push_col("y_scaled", y_scaled);

                // TODO: Handle extra regressors

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

            auto make_all_seasonality_features(const tbl::table& tbl) {

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

            model fit(const tbl::table& tbl) {
                // TODO: Re-enable const on tbl
                // TODO: See if self.history exists, bail if it does
                // TODO: See if ds, y exist in tbl, bail if they do not
                // TODO: Copy tbl to history so we can modify history without worry about tbl
                history = tbl;
                history = setup_table(history, true);
                set_auto_seasonalities();
                auto [seasonal_features, prior_scales, component_cols, modes] = make_all_seasonality_features(history);

                std::map<std::string,
                    std::pair<std::vector<double>, std::vector<size_t> > > vars_r{};

                std::map<std::string,
                    std::pair<std::vector<int>, std::vector<size_t> > > vars_i{};

                // X:
                // cap:
                // sigmas:
                // t:
                // t_change:
                // tau:
                // y:
                std::map<std::string, std::string> m1{
                    {"X", "X"},
                    {"cap", "cap"},
                    {"sigmas", "sigmas"},
                    {"t", "t"},
                    {"t_change", "t_change"},
                    {"tau", "tau"},
                    {"y", "y"}
                };

                for (auto const& [key, val] : m1) {
                    vars_r[key] = load<double>(val);
                }

                // K
                // S
                // T
                // s_a
                // s_m
                // trend_indicator
                std::map<std::string, std::string> m2{
                    {"K", "K"},
                    {"S", "S"},
                    {"T", "T1"},
                    {"s_a", "s_a"},
                    {"s_m", "s_m"},
                    {"trend_indicator", "trend_indicator"}
                };

                for (auto const& [key, val] : m2) {
                    vars_i[key] = load<int>(val);
                }

                stan::io::prophet_var_context context(vars_r, vars_i);

                prophet_model_namespace::prophet_model prophet_model(context);
                std::vector<std::string> constrained_param_names;
                prophet_model.constrained_param_names(constrained_param_names);

                std::map<std::string,
                    std::pair<std::vector<double>, std::vector<size_t> > > vars_r_optimize{};

                std::map<std::string,
                    std::pair<std::vector<int>, std::vector<size_t> > > vars_i_optimize{};

                std::map<std::string, std::string> m3{
                    {"beta", "beta"},
                    {"delta", "delta"},
                    {"k", "k1"},
                    {"m", "m"}
                };

                for (auto const& [key, val] : m3) {
                    vars_r_optimize[key] = load<double>(val);
                }

                std::map<std::string, std::string> m4{
                    {"sigma_obs", "sigma_obs"}
                };

                for (auto const& [key, val] : m4) {
                    vars_i_optimize[key] = load<int>(val);
                }

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

                double k;
                double m;
                std::vector<double> delta;
                double sigma_obs;
                std::vector<double> beta;
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
                return model{};
            }

    };
}

#endif // PROPHET_HPP
