#include <iostream>
#include <fstream>
#include <vector>
// #include <numbers>
#include "prophet.hpp"
#include "table.hpp"
#include "funcs.hpp"

void make_all_seasonality_features(tbl::table tbl) {

}

int main(int argc, char* argv[]) {

    prophet::prophet prophet;
    tbl::table tbl;
    auto fit = prophet.fit(tbl);

}
