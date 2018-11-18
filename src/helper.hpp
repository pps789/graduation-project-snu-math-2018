#ifndef SRC_HELPER_HPP_
#define SRC_HELPER_HPP_

#include<algorithm>
const double DOUBLE_EPS = 1e-6;

double nonzero(double x) {
    return std::abs(x) > DOUBLE_EPS;
}

#endif  // SRC_HELPER_HPP_
