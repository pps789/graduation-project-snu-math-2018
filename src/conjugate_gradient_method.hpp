#ifndef SRC_CONJUGATE_GRADIENT_METHOD_HPP_
#define SRC_CONJUGATE_GRADIENT_METHOD_HPP_

#include<vector>
#include"helper.hpp"

// Assume that A is square.
template<typename T>
std::vector<T> conjugate_gradient_method(
        const std::vector<std::vector<T>>& A,
        std::vector<T> b,
        T TOL) {
    const int N = b.size();
    std::vector<T> x(N), r = b, p = b;
    T rsold = 0;
    for (int i = 0; i < N; i++) rsold += r[i]*r[i];

    while (1) {
        std::vector<T> Ap(N);
        for (int i = 0; i < N; i++) {
            for(int j = 0; j < N; j++) {
                Ap[i] += A[i][j]*p[j];
            }
        }
        T ptAp = 0;
        for (int i = 0; i < N; i++) ptAp += p[i]*Ap[i];
        T alpha = rsold / ptAp;

        for (int i = 0; i < N; i++) x[i] += alpha*p[i];
        for (int i = 0; i < N; i++) r[i] -= alpha*Ap[i];
        T rsnew = 0;
        for (int i = 0; i < N; i++) rsnew += r[i]*r[i];
        if (rsnew < TOL*TOL) break;

        for (int i = 0; i < N; i++) p[i] = r[i] + rsnew/rsold*p[i];
        rsold = rsnew;
    }
}

#endif  // SRC_CONJUGATE_GRADIENT_METHOD_HPP_
