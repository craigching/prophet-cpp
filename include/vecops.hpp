#ifndef VECOPS_HPP
#define VECOPS_HPP

#include <cmath>

template <typename T>
void operator*=(std::vector<T>& vec, T v) {
    std::transform(vec.begin(), vec.end(), vec.begin(),
        std::bind(std::multiplies<T>(), std::placeholders::_1, v));    
}

template <typename T>
std::vector<T> operator*(const std::vector<T>& vec, T v) {
    std::vector<T> out;
    std::transform(vec.begin(), vec.end(), std::back_inserter(out),
        std::bind(std::multiplies<T>(), std::placeholders::_1, v));
    return out;
}

template <typename T>
void operator+=(std::vector<T>& vec, T v) {
    std::transform(vec.begin(), vec.end(), vec.begin(),
        std::bind(std::plus<T>(), std::placeholders::_1, v));    
}

template <typename T>
std::vector<T> operator+(const std::vector<T>& vec, T v) {
    std::vector<T> out;
    std::transform(vec.begin(), vec.end(), std::back_inserter(out),
        std::bind(std::plus<T>(), std::placeholders::_1, v));
    return out;
}

template <typename T>
void operator/=(std::vector<T>& vec, T v) {
    std::transform(vec.begin(), vec.end(), vec.begin(),
        std::bind(std::divides<T>(), std::placeholders::_1, v));    
}

template <typename T>
std::vector<T> operator/(const std::vector<T>& vec, T v) {
    std::vector<T> out;
    std::transform(vec.begin(), vec.end(), std::back_inserter(out),
        std::bind(std::divides<T>(), std::placeholders::_1, v));
    return out;
}

template <typename T>
void operator-=(std::vector<T>& vec, T v) {
    std::transform(vec.begin(), vec.end(), vec.begin(),
        std::bind(std::minus<T>(), std::placeholders::_1, v));    
}

// vector - scalar
template <typename T>
std::vector<T> operator-(const std::vector<T>& vec, T v) {
    std::vector<T> out;
    std::transform(vec.begin(), vec.end(), std::back_inserter(out),
        std::bind(std::minus<T>(), std::placeholders::_1, v));
    return out;
}

// vector - vector
template <typename T>
std::vector<T> operator-(const std::vector<T>& v1, const std::vector<T>& v2) {
    std::vector<T> out;
    std::transform(v1.begin(), v1.end(), v2.begin(), std::back_inserter(out), std::minus<T>());
    return out;
}

auto operator>>(const std::vector<double>& v, std::function<double (const std::vector<double>&)> f) {
    return f(v);
}

auto operator>>(const std::vector<double>& v, std::function<std::vector<double> (const std::vector<double>&)> f) {
    return f(v);
}

using TimePoint = std::chrono::system_clock::time_point;
auto operator>>(std::vector<TimePoint>& v, std::function<TimePoint (std::vector<TimePoint>&)> f) {
    return f(v);
}

namespace vec {
    std::function<double (const std::vector<double>&)> max() {
        return [](const std::vector<double>& v) {
            return *std::max_element(v.begin(),v.end());
        };
    }

    std::function<double (const std::vector<double>&)> min() {
        return [](const std::vector<double>& v) {
            return *std::min_element(v.begin(),v.end());
        };
    }

    std::function<std::vector<double> (const std::vector<double>&)> abs() {
        return [](const std::vector<double>&v ) {
            std::vector<double> out;
            std::transform(v.begin(), v.end(), std::back_inserter(out), [](double v) { return std::abs(v); });
            return out;
        };
    }

    std::function<std::vector<double> (const std::vector<double>&)> diff() {
        return [](const std::vector<double>&v ) {
            std::vector<double> out;
            std::adjacent_difference(v.begin(), v.end(), std::back_inserter(out));
            out[0] = nan("");
            return out;
        };
    }

    std::function<TimePoint (std::vector<TimePoint>&)> min_time_point() {
        return [](std::vector<TimePoint>& v) {
            return *std::min_element(v.begin(),v.end());
        };
    }

    std::function<TimePoint (std::vector<TimePoint>&)> max_time_point() {
        return [](std::vector<TimePoint>& v) {
            return *std::max_element(v.begin(),v.end());
        };
    }

    std::function<std::vector<double> (const std::vector<double>&)> filter(std::function<bool (double)> filter_expr) {
        return [=](const std::vector<double>& v) {
            std::vector<double> out;
            std::copy_if(v.begin(), v.end(), std::back_inserter(out), filter_expr);
            return out;
        };
    }

    std::vector<int> range(size_t n) {
        std::vector<int> result;
        for (auto i = 0; i < n; ++i) {
            result.push_back(i);
        }
        return result;
    }
}

#endif // VECOPS_HPP
