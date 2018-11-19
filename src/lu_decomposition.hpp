#ifndef SRC_LU_DECOMPOSITION_HPP_
#define SRC_LU_DECOMPOSITION_HPP_

#include<algorithm>
#include<utility>
#include<vector>
#include"helper.hpp"

// Assume that Matrix is square.
// Returns (L-I+U, P).
template<typename T>
std::pair<std::vector<std::vector<T>>, std::vector<int>>
LU_decomposition(std::vector<std::vector<T>> A) {
    const int N = A.size();
    std::vector<int> P(N);

    // Initialization.
    for (int i = 0; i < N; i++) P[i] = i;

    for (int i = 0; i < N; i++) {
        // Pivoting.
        int t = i;
        T maxval = 0;
        for (int j = i; j < N; j++) {
            auto curval = A[j][i];
            for (int k = 0; k < i; k++) {
                curval -= A[j][k]*A[k][i];
            }
            if (std::abs(curval) > std::abs(maxval)) {
                maxval = std::abs(curval);
                t = j;
            }
        }

        std::swap(A[i], A[t]); std::swap(P[i], P[t]);
        if (!nonzero(A[i][i]))
            return {std::vector<std::vector<T>>(), std::vector<int>()};

        // Update U.
        for (int j = i; j < N; j++) {
            for (int k = 0; k < i; k++) {
                A[i][j] -= A[i][k]*A[k][j];
            }
        }

        // Update L.
        for (int j = i+1; j < N; j++) {
            for (int k = 0; k < i; k++) {
                A[j][i] -= A[j][k]*A[k][i];
            }
            A[j][i] /= A[i][i];
        }
    }

    return {A, P};
}

template<typename T>
std::vector<T> LU_solve(
        const std::pair<std::vector<std::vector<T>>, std::vector<int>>& LUP,
        const std::vector<T>& b) {
    const auto& LU = LUP.first;
    const auto& P = LUP.second;
    const int N = P.size();
    std::vector<T> Pb(N);
    for (int i = 0; i < N; i++) {
        Pb[P[i]] = b[i];
    }

    std::vector<T> y(N);
    for (int i = 0; i < N; i++) {
        y[i] = Pb[i];
        for (int j = 0; j < i; j++) {
            y[i] -= LU[i][j]*y[j];
        }
    }

    std::vector<T> x(N);
    for (int i = N-1; i >= 0; i--) {
        x[i] = y[i];
        for (int j = i+1; j < N; j++) {
            x[i] -= LU[i][j]*x[j];
        }
        x[i] /= LU[i][i];
    }

    return x;
}

#endif  // SRC_LU_DECOMPOSITION_HPP_
