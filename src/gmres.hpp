#ifndef SRC_GMRES_HPP_
#define SRC_GMRES_HPP_

#include<cmath>
#include<vector>
#include"helper.hpp"
#include"qr_decomposition.hpp"

// Assume that A is square.
template<typename T>
std::vector<T> gmres(
        const std::vector<std::vector<T>>& A,
        std::vector<T> b,
        std::vector<T> x0,
        T TOL) {
    const int N = A.size();
    std::vector<T> v = b;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            v[i] -= A[i][j]*x0[j];
        }
    }
    T vnorm = 0;
    for (int i = 0; i < N; i++) vnorm += v[i]*v[i];
    vnorm = std::sqrt(vnorm);
    for (int i = 0; i < N; i++) v[i] /= vnorm;

    std::vector<std::vector<T>> V(N);
    std::vector<std::vector<T>> H;

    for (int M = 1; ; M++) {
        for (int i = 0; i < N; i++) V[i].push_back(v[i]);
        // Now, V is N by M.
        std::vector<T> w(N);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                w[j] += A[i][j]*v[j];
            }
        }

        // Update H's m-th column; H is M+1 by M from now.
        H.resize(M+1);
        for (int i = 0; i < M+1; i++) H[i].resize(M);
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                H[i].back() += w[j]*V[i][j];
            }
            for (int j = 0; j < N; j++) {
                w[j] -= H[i].back()*V[i][j];
            }
        }
        for (int i = 0; i < N; i++) H[M].back() += w[i]*w[i];
        H[M].back() = std::sqrt(H[M]);
        for (int i = 0; i < N; i++) v[i] = w[i]/H[M].back();

        // Least-square phase: use QR decomposition.
        auto QR = QR_decomposition(H);
        std::vector<T> e(M+1); e[0] = vnorm;
        auto y = QR_solve(QR, e);

        std::vector<T> x = x0;
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < M; j++) {
                x[j] += V[i][j]*y[j];
            }
        }

        auto r = b;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < M; j++) {
                r[j] -= A[i][j]*x[j];
            }
        }
        T res = 0;
        for (int i = 0; i < N; i++) res += r[i]*r[i];
        if (res < TOL*TOL) return x;
    }
}

#endif  // SRC_GMRES_HPP_
