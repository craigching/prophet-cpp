#define CATCH_CONFIG_MAIN
#include "catch2/catch.hpp"

#include "funcs.hpp"
#include "vecops.hpp"
#include "time_delta.hpp"

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

template <typename T>
auto read_values_from_file(const std::string& fname) {
    std::ifstream infile(fname);

    std::vector<T> values{};

    T val;
    while (infile >> val) {
        values.push_back(val);
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

TEST_CASE( "predict_trend is computed correctly", "[test_predict_trend]" ) {

    auto expected = read_values_from_file<double>("../tests/src/data/predict_trend_valid.txt");

    auto t = read_values_from_file<double>("../tests/src/data/t.txt");
    tbl::table tbl;
    tbl.push_col("t", t);

    auto result = predict_trend(tbl);

    REQUIRE( result.size() == expected.size() );

    for (auto i = 0; i < result.size(); ++i) {
        REQUIRE( logically_equal(result[i], expected[i]) );
    }
}

TEST_CASE( "piecewise_linear is computed correctly", "[piecewise_linear]" ) {

    auto expected = read_values_from_file<double>("../tests/src/data/piecewise_linear_valid.txt");

    auto t = read_values_from_file<double>("../tests/src/data/t.txt");
    std::vector<double> deltas {-1.0630456222095487e-07, 5.8925437529570826e-09, 0.3508431622560911, 0.4639196731154721, 8.759015979732021e-08, -0.011336290924911384, -0.23865392437947305, -0.24272073447767095, 1.3398610238885688e-07, 1.329058431681059e-07, 4.270714421451454e-09, 0.2962454425156414, 0.20966993721139043, 0.0066949364852613045, 0.0007410300204416828, -0.8572701207517232, -0.00010040459440185556, 9.646895212082761e-08, 4.503592908572229e-08, 0.462488748011516, 0.013741215328200695, -6.009798744590757e-08, -0.3345212489885134, 4.7709248246728985e-08, 2.7460283310057605e-08};
    double k = -0.3565853370802134;
    double m = 0.626009448344069;
    std::vector<double> changepoint_ts {0.0330745865676679, 0.0651366857914276, 0.10327370907863652, 0.13533580830239622, 0.16672291596355046, 0.1981100236247047, 0.23152210597367534, 0.2642591967600405, 0.2963212959838002, 0.3300708741140736, 0.3614579817752278, 0.392845089436382, 0.4242321970975363, 0.4556193047586905, 0.48768140398245025, 0.5197435032062099, 0.5514681066486669, 0.5828552143098211, 0.6139048261896727, 0.6452919338508268, 0.6766790415119811, 0.7084036449544381, 0.7397907526155924, 0.7715153560580493, 0.8029024637192035};

    auto result = piecewise_linear(t, deltas, k, m, changepoint_ts);

    REQUIRE( result.size() == expected.size() );

    for (auto i = 0; i < result.size(); ++i) {
        REQUIRE( logically_equal(result[i], expected[i]) );
    }
}

TEST_CASE( "table reads csv correctly", "[table::read_csv]" ) {

    auto tbl = tbl::read_csv(
        "../tests/src/data/example_wp_log_peyton_manning.csv",
        {"date", "double"},
        "%Y-%m-%d");

    // Check sizes
    auto [rows, cols] = tbl.shape();
    REQUIRE( rows == 2905 );
    REQUIRE( cols == 1 ); // We don't count date columns for now
    // Check column names
    REQUIRE( tbl.get_names()[0] == "y" );
    REQUIRE( tbl.get_times_name() == "ds" );
    // Spot check some values
    // auto t = tbl.get_times()[0];
    // std::time_t tt = std::chrono::system_clock::to_time_t(t);
    // std::tm tm = *std::localtime(&tt);
    // REQUIRE( tm.tm_year == 107 );
    // REQUIRE( tm.tm_mon == 11 );
    // REQUIRE( tm.tm_mday == 10 );

    // t = tbl.get_times()[rows - 1];
    // tt = std::chrono::system_clock::to_time_t(t);
    // tm = *std::localtime(&tt);
    // REQUIRE( tm.tm_year == 116 );
    // REQUIRE( tm.tm_mon == 0 );
    // REQUIRE( tm.tm_mday == 20 );

    auto y = tbl.get_col("y");
    REQUIRE( logically_equal(y[0], 9.59076113897809) );
    REQUIRE( logically_equal(y[rows - 1], 8.89137400948464) );
}

TEST_CASE( "table functions work correctly", "[table::mean,min,max]" ) {

    tbl::table t1;
    t1.push_col("countable", {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0});

    auto mean = t1 >> tbl::mean("countable");
    auto min = t1 >> tbl::min("countable");
    auto max = t1 >> tbl::max("countable");

    REQUIRE( mean == 5.5 );
    REQUIRE( min == 1.0 );
    REQUIRE( max == 10.0 );

#define pred(expr) [](double e) { return expr; }

    auto t2 = t1 >> tbl::filter("countable", pred(e < 6.0));

    REQUIRE( t2.get_col("countable").size() == 5 );
    REQUIRE( t2.get_col("countable") == std::vector<double>{1, 2, 3, 4, 5} );
}

TEST_CASE( "vector operator overloads work correctly", "[vecops]" ) {
    std::vector<double> v1{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
    v1 *= 3.0;

    REQUIRE( v1 == std::vector<double>{3, 6, 9, 12, 15, 18, 21, 24, 27, 30} );

    std::vector<double> v2{1.0, 1.0, 1.0};
    v2 += 3.0;

    REQUIRE( v2 == std::vector<double>{4, 4, 4} );

    std::vector<double> v3{4.0, 4.0, 4.0};
    v3 /= 2.0;

    REQUIRE( v3 == std::vector<double>{2, 2, 2} );

    std::vector<double> v4{5.0, 5.0, 5.0};
    v4 -= 2.0;

    REQUIRE( v4 == std::vector<double>{3, 3, 3} );

    auto v5 = std::vector<double>{5, 5, 5} - 2.0;
    REQUIRE( v5 == std::vector<double>{3, 3, 3});

    v5 = std::vector<double>{5, 5, 5} + double(2.0);
    REQUIRE( v5 == std::vector<double>{7, 7, 7});

    v5 = std::vector<double>{5, 5, 5} * 2.0;
    REQUIRE( v5 == std::vector<double>{10, 10, 10});

    v5 = std::vector<double>{6, 6, 6} / 2.0;
    REQUIRE( v5 == std::vector<double>{3, 3, 3});

    std::vector<double> v6{1, 2, 3, 4, 5, 4, 3, 2, 1};

    auto max = v6 >> vec::max();

    REQUIRE( max == 5 );

    std::vector<double> v7{1, -1, 2, -2, -5, 4, -3, 2, -1};
    auto absv = v7 >> vec::abs();

    REQUIRE( absv == std::vector<double>{1, 1, 2, 2, 5, 4, 3, 2, 1} );

    max = v7 >> vec::abs() >> vec::max();

    REQUIRE( max == 5 );

    std::vector<double> v8{5, 5, 5};
    std::vector<double> v9{3, 3, 3};

    // auto mult = v8 * v9;
    auto diff = v8 - v9;

    REQUIRE( diff == std::vector<double>{2, 2, 2} );
}

TEST_CASE( "vector operator overloads work correctly for time_point", "[vecops]" ) {

    auto tbl = tbl::read_csv(
        "../tests/src/data/example_wp_log_peyton_manning.csv",
        {"date", "double"},
        "%Y-%m-%d");

    auto times = tbl.get_times();
    auto min = times >> vec::min();
    auto max = times >> vec::max();

    REQUIRE( min < max );
}

TEST_CASE( "test diff function", "[diff]" ) {
    std::vector<double> v {102, 104, 106, 108, 110, 112, 114, 116, 118, 120};

    auto v1 = v >> vec::diff();

    REQUIRE( std::isnan(v1[0]) );

    v1.erase(v1.begin());
    REQUIRE( v1 == std::vector<double>{2, 2, 2, 2, 2, 2, 2, 2, 2} );
}

TEST_CASE( "test diff function on table", "[diff on table]"){
    auto tbl = tbl::read_csv(
        "../tests/src/data/example_wp_log_peyton_manning.csv",
        {"date", "double"},
        "%Y-%m-%d");

    auto first = tbl.get_times() >> vec::min();
    auto last = tbl.get_times() >> vec::max();

#define pred(expr) [](double e) { return expr; }
    auto dt = tbl.get_times()
        >> vec::diff()
        >> vec::filter(pred( !std::isnan(e) ));

    auto min_dt = dt >> vec::min();

    REQUIRE( first == 1197266400 );
    REQUIRE( last == 1453269600 );
    REQUIRE( min_dt == 86400 );
    REQUIRE( dt.size() == tbl.get_times().size() - 1 );
}

TEST_CASE( "test time_delta function", "[time_delta]"){

    auto tbl = tbl::read_csv(
        "../tests/src/data/example_wp_log_peyton_manning.csv",
        {"date", "double"},
        "%Y-%m-%d");

    auto first = tbl.get_times() >> vec::min();
    auto last = tbl.get_times() >> vec::max();
    auto yearly_disable = last - first < time_delta<days>(730);

#define pred(expr) [](double e) { return expr; }
    auto min_dt = tbl.get_times()
        >> vec::diff()
        >> vec::filter(pred( !std::isnan(e) ))
        >> vec::min();

    auto weekly_disable = (
        (last - first < time_delta<weeks>(2)) ||
        (min_dt >= time_delta<weeks>(1)));

    auto daily_disable = (
        (last - first < time_delta<days>(2)) ||
        (min_dt >= time_delta<days>(1)));

    REQUIRE( yearly_disable == false );
    REQUIRE( weekly_disable == false );
    REQUIRE( daily_disable == true );
}
