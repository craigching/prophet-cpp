#ifndef TIME_DELTA_HPP
#define TIME_DELTA_HPP

enum TimeEnum { days, weeks };

template <TimeEnum E>
double time_delta(int n) {
    switch(E) {
        case days:
            return n * 24 * 60 * 60;
        case weeks:
            return n * 7 * 24 * 60 * 60;
    }
    return 0;
}

#endif // TIME_DELTA_HPP
