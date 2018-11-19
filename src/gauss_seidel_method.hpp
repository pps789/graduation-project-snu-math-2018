#ifndef SRC_GAUSS_SEIDEL_METHOD_HPP_
#define SRC_GAUSS_SEIDEL_METHOD_HPP_

#include<vector>
#include"helper.hpp"

// Assume that A is square.
template<typename T>
std::vector<T> gauss_seidel_method(
        const std::vector<std::vector<T>>& A,
        std::vector<T> b,
        T TOL) {
    const int N = b.size();

    // Initial guessing: zero vector.
    std::vector<T> x(N);
    while (1) {
        std::vector<T> s = b;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < i; j++) {
                s[i] -= A[i][j]*s[j];
            }
            for (int j = i+1; j < N; j++) {
                s[i] -= A[i][j]*x[j];
            }
        }

        for (int i = 0; i < N; i++) s[i] /= A[i][i];
        T delta = 0;
        for (int i = 0; i < N; i++) delta += (x[i]-s[i])*(x[i]-s[i]);
        if (delta < TOL * TOL) return s;
        x = s;
    }
}

#endif  // SRC_GAUSS_SEIDEL_METHOD_HPP_
