#ifndef FUNCS_HPP
#define FUNCS_HPP

#include "table.hpp"

#include <chrono>
#include <iomanip>
#include <cmath>

inline bool logically_equal(double a, double b, double error_factor=1.0)
{
  return a==b ||
    std::abs(a-b)<std::abs(std::min(a,b))*std::numeric_limits<double>::epsilon()*
                  error_factor;
}

std::string first(const std::string& str, const std::string& delim) {
    return str.substr(0, str.find(delim));
}

std::vector<double> date_range(double start, int periods, const std::string& freq) {
    std::vector<double> range;
    if (freq == "D") {
        for (auto i = 0; i < periods; ++i) {
            range.push_back(start + (i * 24 * 60 * 60));
        }
    }
    return range;
}

int to_day(time_t tt) {
    using std::chrono::system_clock;
    using std::chrono::time_point_cast;
    typedef std::chrono::duration<int, std::ratio<60 * 60 * 24>> days_type;
    using std::chrono::time_point;

    system_clock::time_point tp = system_clock::from_time_t(tt);
    time_point<system_clock, days_type> day = time_point_cast<days_type>(tp);
    return day.time_since_epoch().count();
}

int date_to_day(const std::string& date){
    std::tm t{};
    std::istringstream ss(date);
    ss >> std::get_time(&t, "%Y-%m-%d %HH:%MM%SS");

    time_t tt =  mktime(&t);

    return to_day(tt);
}

std::vector<double> linspace(int start_in, int end_in, int num_in) {

    std::vector<double> linspaced;

    double start = static_cast<double>(start_in);
    double end = static_cast<double>(end_in);
    double num = static_cast<double>(num_in);

    if (num == 0) {
        return linspaced;
    }

    if (num == 1) {
        linspaced.push_back(start);
        return linspaced;
    }

    double delta = (end - start) / (num - 1);

    for(int i=0; i < num-1; ++i) {
        linspaced.push_back(start + delta * i);
    }

    linspaced.push_back(end);

    return linspaced;
}

tbl::table fourier_series(std::vector<double> dates, double period, int series_order) {

    std::vector<int> t;

    for (auto& date: dates) {
        t.push_back(to_day(date));
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

tbl::table make_seasonality_features(std::vector<double> dates, double period, int series_order, const std::string& prefix) {

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

#endif // FUNCS_HPP
