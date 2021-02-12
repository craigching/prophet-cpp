#ifndef FUNCS_HPP
#define FUNCS_HPP

#include "table.hpp"

#include <chrono>
#include <iomanip>
#include <cmath>

int date_to_day(const std::string& date){

    using std::chrono::time_point;
    using std::chrono::time_point_cast;
    using std::chrono::system_clock;
    typedef std::chrono::duration<int, std::ratio<60 * 60 * 24>> days_type;

    std::tm t{};
    std::istringstream ss(date);
    ss >> std::get_time(&t, "%Y-%m-%d %HH:%MM%SS");

    time_t tt =  mktime(&t);

    system_clock::time_point tp = system_clock::from_time_t(tt);

    time_point<system_clock, days_type> day = time_point_cast<days_type>(tp);

    return day.time_since_epoch().count();
}

tbl::table fourier_series(std::vector<std::string> dates, double period, int series_order) {

    std::vector<int> t;

    for (auto& date: dates) {
        t.push_back(date_to_day(date));
    }

    tbl::table tbl{series_order * 2};

    for (auto& v: t) {
        std::vector<double> row;
        for (int i = 0; i < series_order; ++i) {
            row.push_back(sin(2.0 * (i + 1) * M_PI * v / period));
            row.push_back(cos(2.0 * (i + 1) * M_PI * v / period));
        }
        tbl.push_row(row);
    }

    // TODO figure out move semantics here
    return tbl;
}

tbl::table make_seasonality_features(std::vector<std::string> dates, double period, int series_order, const std::string& prefix) {

    auto tbl = fourier_series(dates, period, series_order);

    std::vector<std::string> names;
    for (auto i = 0; i < tbl.shape().second; ++i) {
        std::ostringstream name;
        name << prefix << "_delim_" << i + 1;
        names.push_back(name.str());
    }

    tbl.set_names(names);

    // TODO figure out move semantics here
    return tbl;
}

std::vector<double> piecewise_linear(std::vector<double> t, std::vector<double> deltas, double k, double m, std::vector<double> changepoint_ts) {
    std::vector<double> gammas;

    auto cp_iter = changepoint_ts.begin();
    auto cp_end = changepoint_ts.end();
    auto deltas_iter = deltas.begin();
    auto deltas_end = deltas.end();

    for (; cp_iter != cp_end && deltas_iter != deltas_end; ++cp_iter, ++deltas_iter) {
        gammas.push_back(-*cp_iter * *deltas_iter);
    }

    std::vector<double> k_t(t.size(), k);
    std::vector<double> m_t(t.size(), m);

    std::vector<size_t> indx(t.size());
    size_t s = 0;
    for (auto t_s: changepoint_ts) {
        size_t i = 0;
        for (auto t_v: t) {
            if (t_v >= t_s) {
                k_t[i] += deltas[s];
                m_t[i] += gammas[s];
            }
            ++i;
        }
        ++s;
    }

    std::vector<double> result;

    auto k_t_iter = k_t.begin();
    auto t_iter = t.begin();
    auto m_t_iter = m_t.begin();

    for (; k_t_iter != k_t.end(); ++k_t_iter, ++t_iter, ++m_t_iter) {
        result.push_back(*k_t_iter * *t_iter + *m_t_iter);
    }

    return result;
}

std::vector<double> predict_trend(const tbl::table& tbl) {
    double y_scale_scalar = 12.846746888829;
    const auto& t = tbl.get_col("t");
    std::vector<double> deltas {-1.0630456222095487e-07, 5.8925437529570826e-09, 0.3508431622560911, 0.4639196731154721, 8.759015979732021e-08, -0.011336290924911384, -0.23865392437947305, -0.24272073447767095, 1.3398610238885688e-07, 1.329058431681059e-07, 4.270714421451454e-09, 0.2962454425156414, 0.20966993721139043, 0.0066949364852613045, 0.0007410300204416828, -0.8572701207517232, -0.00010040459440185556, 9.646895212082761e-08, 4.503592908572229e-08, 0.462488748011516, 0.013741215328200695, -6.009798744590757e-08, -0.3345212489885134, 4.7709248246728985e-08, 2.7460283310057605e-08};
    double k = -0.3565853370802134;
    double m = 0.626009448344069;
    std::vector<double> changepoint_ts {0.0330745865676679, 0.0651366857914276, 0.10327370907863652, 0.13533580830239622, 0.16672291596355046, 0.1981100236247047, 0.23152210597367534, 0.2642591967600405, 0.2963212959838002, 0.3300708741140736, 0.3614579817752278, 0.392845089436382, 0.4242321970975363, 0.4556193047586905, 0.48768140398245025, 0.5197435032062099, 0.5514681066486669, 0.5828552143098211, 0.6139048261896727, 0.6452919338508268, 0.6766790415119811, 0.7084036449544381, 0.7397907526155924, 0.7715153560580493, 0.8029024637192035};

    std::vector<double> trend = piecewise_linear(t, deltas, k, m, changepoint_ts);

    std::vector<double> y_scale(trend.size(), y_scale_scalar);
    std::vector<double> floor(trend.size(), 0.0);

    std::vector<double> result;

    //         return trend * self.y_scale + df['floor']

    auto trend_iter = trend.begin();
    auto y_scale_iter = y_scale.begin();
    auto floor_iter = floor.begin();

    for ( ; trend_iter != trend.end(); ++trend_iter, ++y_scale_iter, ++floor_iter) {
        result.push_back(*trend_iter * *y_scale_iter + *floor_iter);
    }

    return result;
}


#endif // FUNCS_HPP
