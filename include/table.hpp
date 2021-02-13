#ifndef TABLE_H
#define TABLE_H

#include <vector>
#include <string>
#include <iostream>
#include <chrono>
#include <iomanip>
#include <sstream>

#include "csv.h"

namespace tbl {
    class table {

    std::vector<std::string> names;
    std::vector<std::vector<double>> columns;
    std::string times_name;
    std::vector<std::chrono::system_clock::time_point> times;

    public:
        explicit table() {};
        table(table const&) = default;
        table(table&&) = default;
        ~table() = default;

        table(const size_t& cols) {
            for (auto i = 0; i < cols; ++i) {
                columns.push_back(std::vector<double>{});
            }
        }

        void push_col(const std::string& name, std::vector<double> column) {
            names.push_back(name); // TODO Check for duplicate names
            columns.push_back(column);
        }

        void push_time(const std::chrono::system_clock::time_point& t) {
            times.push_back(t);
        }

        void push_row(const std::vector<double>& row) {
            // TODO check that row len == num cols

            std::vector<double>::const_iterator v;
            std::vector<std::vector<double>>::iterator col;

            for (v = row.begin(), col = columns.begin(); v != row.end(); ++v, ++col) {
                (*col).push_back(*v);
            }
        }

        auto get(const size_t row, const size_t col) const {
            return columns[col][row];
        }

        auto get_col(const std::string& name) const {
            for (auto i = 0; i < names.size(); ++i) {
                if (names[i] == name) {
                    return columns[i];
                }
            }

            return std::vector<double>{};
        }

        auto get_times() const {
            return times;
        }

        void print() const {

            // TODO Lots of assumptions here. What if there are 0 columns? Need to assert all columns have same number of rows
            auto cols = columns.size();
            auto rows = columns[0].size();

            std::cout.precision(17);

            for (auto name: names) {
                std::cout << name << "  ";
            }
            std::cout << "\n";

            for (auto i = 0; i < rows; ++i) {
                for (auto j = 0; j < cols; ++j) {
                    std::cout << std::fixed << get(i, j) << "  ";
                }
                std::cout << "\n";
            }
        }

        auto shape() const {
            return std::make_pair(columns[0].size(), columns.size());
        }

        void set_names(std::vector<std::string> new_names) {
            // TODO number of new_names == number of columns
            names = new_names;
        }

        void set_times_name(const std::string& name) {
            times_name = name;
        }

        auto get_names() const {
            return names;
        }

        auto get_times_name() const {
            return times_name;
        }

        static table read_csv(
            const std::string& fname,
            const std::vector<std::string>& types,
            const std::string& date_format = "%Y-%m-%d %HH:%MM%SS",
            const std::string& delim = ",",
            const bool header = true
        ) {
            tbl::table tbl;
            std::ifstream file(fname);
            auto i = 0;
            for(auto& row: CSVRange(file)) {

                if (i == 0 && header) {
                    for (auto j = 0; j < row.size(); ++j) {
                        auto name = std::string(row[j]);
                        // Strip double quotes
                        name.erase(remove(name.begin(), name.end(), '"'), name.end());
                        if (types[j] == "date") {
                            tbl.set_times_name(name);
                        } else {
                            tbl.push_col(name, std::vector<double>{});
                        }
                    }
                    ++i;
                    continue;
                }

                std::vector<double> values_row;
                for (auto j = 0; j < row.size(); ++j) {
                    auto type = types[j];
                    auto value = std::string{row[j]};
                    if (type == "date") {
                        std::tm t{};
                        // Strip double quotes
                        value.erase(remove(value.begin(), value.end(), '"'), value.end());
                        std::istringstream ss{value};
                        ss >> std::get_time(&t, date_format.c_str());

                        time_t tt =  mktime(&t);

                        tbl.push_time(std::chrono::system_clock::from_time_t(tt));
                    } else if (type == "double") {
                        values_row.push_back(std::stod(value));
                    }
                }

                tbl.push_row(values_row);
                ++i;
            }

            return tbl;
        }
    };
}

#endif // TABLE_H
