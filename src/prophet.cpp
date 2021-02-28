#include "prophet.hpp"
#include "table.hpp"

int main(int argc, char* argv[]) {
    prophet::prophet m;
    auto tbl = tbl::read_csv(
        "../tests/src/data/example_wp_log_peyton_manning.csv",
        {"date", "double"},
        "%Y-%m-%d");

    m.fit(tbl);
    auto future = m.make_future_dataframe(365);
}
