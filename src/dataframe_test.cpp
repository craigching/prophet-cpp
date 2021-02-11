#include "dataframe.hpp"

int main(int argc, char *argv[]) {
    dataframe df{10};
    df.push_row(std::vector<double>{0, 1, 2, 3, 4, 5, 6, 7, 8, 9});
    df.push_row(std::vector<double>{10, 11, 12, 13, 14, 15, 16, 17, 18, 19});
    df.print();
}