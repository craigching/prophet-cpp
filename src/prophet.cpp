#include "prophet.hpp"
#include "table.hpp"

int main(int argc, char* argv[]) {
    prophet::prophet prophet;
    auto tbl = tbl::read_csv(
        "../tests/src/data/example_wp_log_peyton_manning.csv",
        {"date", "double"},
        "%Y-%m-%d");

    auto fit = prophet.fit(tbl);
}
