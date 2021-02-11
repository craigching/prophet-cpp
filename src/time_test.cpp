#include <iostream>
#include <chrono>
#include <time.h>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <functional>
#include <math.h>


int date_to_day(const std::string& date){

    using std::chrono::time_point;
    using std::chrono::time_point_cast;
    using std::chrono::system_clock;
    typedef std::chrono::duration<int, std::ratio<60 * 60 * 24>> days_type;

    std::tm t{};
    std::istringstream ss(date);
    ss >> std::get_time(&t, "%Y-%m-%d");

    time_t tt =  mktime(&t);

    system_clock::time_point tp = system_clock::from_time_t(tt);

    time_point<system_clock, days_type> day = time_point_cast<days_type>(tp);

    return day.time_since_epoch().count();
}

double foo(double d) {
    return d;
}

int main(){
    // int d = 20160718;
    // std::string d("2007-12-10");
    // std::cout << d << ": " << date_to_day(d) << std::endl; //17000
    std::vector<std::function<double(double)>> funs;

    funs.push_back(foo);
    funs.push_back([](double v) -> double { return sin(v); });

    std::cout << sin(90.0) << std::endl;

    return 0;
}
