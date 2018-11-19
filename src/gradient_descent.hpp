#ifndef SRC_GRADIENT_DESCENT_HPP_
#define SRC_GRADIENT_DESCENT_HPP_

#include<vector>
#include"helper.hpp"

// Assume that A is square.
template<typename T>
std::vector<T> gradient_descent(
        const std::vector<std::vector<T>>& A,
        std::vector<T> b,
        T TOL) {
    const int N = b.size();
    std::vector<T> x(N), r = b, rs = 0;
    for (int i = 0; i < N; i++) rs += r[i]*r[i];

    while (1) {
        std::vector<T> Ar;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                Ar[i] += A[i][j]*r[j];
            }
        }
        T rtAr = 0;
        for (int i = 0; i < N; i++) rtAr += r[i]*Ar[i];
        T alpha = rs/rtAr;

        for (int i = 0; i < N; i++) x[i] += alpha*r[i];
        for (int i = 0; i < N; i++) r[i] -= alpha*Ar[i];

        rs = 0;
        for (int i = 0; i < N; i++) rs += r[i]*r[i];

        if (rs < TOL*TOL) return x;
    }
}

#endif  // SRC_GRADIENT_DESCENT_HPP_
