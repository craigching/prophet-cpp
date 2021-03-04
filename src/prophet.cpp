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
    std::cout << "future size: " << future.get_times().size() << std::endl;
    // std::cout << "future:" << std::endl;
    // auto [ cols, rows ] = future.shape();
    // std::cout << "future rows: " << rows << ", cols: " << cols << std::endl;
    // std::cout << "future done." << std::endl;
    m.predict(future);
}
