#ifndef TABLE_H
#define TABLE_H

#include <vector>
#include <string>
#include <iostream>

namespace tbl {
    class table {

    std::vector<std::string> names;
    std::vector<std::vector<double>> columns;

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

        void push_col(std::string& name, std::vector<double> column) {
            names.push_back(name); // TODO Check for duplicate names
            columns.push_back(column);
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

        auto get_names() const {
            return names;
        }
    };
}

#endif // TABLE_H
