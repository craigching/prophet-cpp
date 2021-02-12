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

#endif // FUNCS_HPP
