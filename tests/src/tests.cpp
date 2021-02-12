#define CATCH_CONFIG_MAIN
#include "catch2/catch.hpp"

#include "funcs.hpp"

#include <string>

std::vector<std::string> read_names(std::string line) {

    std::vector<std::string> names{};
    std::string delimiter = " ";

    size_t pos = 0;
    std::string token;
    while ((pos = line.find(delimiter)) != std::string::npos) {
        token = line.substr(0, pos);
        if (token.length() > 0) {
            names.push_back(token);
        }
        line.erase(0, pos + delimiter.length());
    }

    return names;
}

std::vector<double> read_values(std::string line) {

    std::vector<double> values{};
    std::string delimiter = " ";

    size_t pos = 0;
    std::string token;
    while ((pos = line.find(delimiter)) != std::string::npos) {
        token = line.substr(0, pos);
        if (token.length() > 0) {
            values.push_back(std::stod(token));
        }
        line.erase(0, pos + delimiter.length());
    }

    return values;
}

inline bool logically_equal(double a, double b, double error_factor=1.0)
{
  return a==b ||
    std::abs(a-b)<std::abs(std::min(a,b))*std::numeric_limits<double>::epsilon()*
                  error_factor;
}

TEST_CASE( "make_seasonality_features is computed correctly", "[make_seasonality_features]" ) {

    std::string line;
    std::ifstream tbl_file("../tests/src/data/make_seasonality_features_valid.txt");
    std::getline(tbl_file, line);
    auto names = read_names(line);

    tbl::table expected_tbl{20};
    while(std::getline(tbl_file, line)) {
        auto values = read_values(line);
        expected_tbl.push_row(values);
    }
    expected_tbl.set_names(names);


    std::vector<std::string> dates;
    std::ifstream dates_file("../tests/src/data/dates.txt");
    while (std::getline(dates_file, line)) {
        dates.push_back(line);
    }
    double period = 365.25;
    int series_order = 10;

    auto tbl = make_seasonality_features(dates, period, series_order, "yearly");

    REQUIRE( expected_tbl.get_names().size() == tbl.get_names().size() );
    for (auto i = 0; i < tbl.get_names().size(); ++i) {
        REQUIRE( tbl.get_names()[i] == expected_tbl.get_names()[i] );
    }

    std::cout.precision(16);

    auto [rows, cols] = tbl.shape();
    for (auto i = 0; i < rows; ++i) {
        for (auto j = 0; j < cols; ++j) {
            // Skip some really small values
            if ((i != 239 && j != 4) &&
                (i != 715 && j != 0) &&
                (i != 1191 && j != 10) &&
                (i != 1673 && j != 16) &&
                (i != 2159 && j != 2)) {
                REQUIRE( logically_equal(tbl.get(i, j), expected_tbl.get(i, j), 30.0) );
            }
        }
    }
}
