#ifndef PROPHET_HPP
#define PROPHET_HPP

#include "table.hpp"
#include "prophet_model.hpp"
#include "prophet_var_context.hpp"
#include "vecops.hpp"
#include "time_delta.hpp"

#include <stan/services/optimize/lbfgs.hpp>
#include <stan/callbacks/interrupt.hpp>
#include <stan/callbacks/stream_logger.hpp>
#include <stan/callbacks/stream_writer.hpp>

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
        double period;
        int fourier_order;
        float prior_scale;
        std::string mode;
        std::string condition_name;
    };

    class prophet {


        tbl::table history;
        double start = 0.0;
        double y_scale = 0.0;
        bool logistic_floor = false;
        double t_scale = 0.0;
        std::map<std::string, seasonality> seasonalities;
        double seasonality_prior_scale;
        std::string seasonality_mode;

        std::string yearly_seasonality;
        std::string weekly_seasonality;
        std::string daily_seasonality;


        public:
            explicit prophet()
            : seasonality_prior_scale(10.0),
              seasonality_mode("additive"),
              yearly_seasonality("auto"),
              weekly_seasonality("auto"),
              daily_seasonality("auto")
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
                    s.period = 365.25;
                    s.fourier_order = fourier_order;
                    s.prior_scale = seasonality_prior_scale;
                    s.mode = seasonality_mode;
                    s.condition_name = "";
                    seasonalities["yearly"] = s;
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
                    s.period = 7;
                    s.fourier_order = fourier_order;
                    s.prior_scale = seasonality_prior_scale;
                    s.mode = seasonality_mode;
                    s.condition_name = "";
                    seasonalities["weekly"] = s;
                }

                auto daily_disable = (
                    (last - first < time_delta<days>(2)) ||
                    (min_dt >= time_delta<days>(1)));
                fourier_order = parse_seasonality_args("daily", daily_seasonality, daily_disable, 4);
                if (fourier_order > 0) {
                    seasonality s;
                    s.period = 1;
                    s.fourier_order = fourier_order;
                    s.prior_scale = seasonality_prior_scale;
                    s.mode = seasonality_mode;
                    s.condition_name = "";
                    seasonalities["daily"] = s;
                }
            }

            int parse_seasonality_args(const std::string& name, const std::string& arg, bool auto_disable, int default_order) {
                int fourier_order = std::stoi(arg);
                if (arg == "auto") {
                    if (seasonalities.find(name) != seasonalities.end()) {
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
                }

                return fourier_order;
            }

            void make_all_seasonality_features(const tbl::table& tbl) {

            }

            model fit(const tbl::table& tbl) {
                // TODO: Re-enable const on tbl
                // TODO: See if self.history exists, bail if it does
                // TODO: See if ds, y exist in tbl, bail if they do not
                // TODO: Copy tbl to history so we can modify history without worry about tbl
                history = tbl;
                history = setup_table(history, true);
                set_auto_seasonalities();
                make_all_seasonality_features(history);

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
