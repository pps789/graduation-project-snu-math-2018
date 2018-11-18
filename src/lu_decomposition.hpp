#ifndef SRC_LU_DECOMPOSITION_HPP_
#define SRC_LU_DECOMPOSITION_HPP_

#include<algorithm>
#include<utility>
#include<vector>
#include"helper.hpp"

// Assume that Matrix is square.
// Returns (L-I+U, P).
template<typename Matrix, typename Vector>
std::pair<Matrix, std::vector<int>>
LU_decomposition(Matrix A) {
    const int N = A.size();
    std::vector<int> P(N);

    // Initialization.
    for (int i = 0; i < N; i++) P[i] = i;

    for (int i = 0; i < N; i++) {
        // Pivoting.
        int t = i;
        auto maxval = A[i][i];
        for (int k = 0; k < i; k++) {
            maxval -= A[i][k]*A[k][i];
        }
        maxval = std::abs(maxval);

        for (int j = i+1; j < N; j++) {
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
        if (!nonzero(A[i][i])) return {Matrix(), Vector()};

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

template<typename Matrix, typename Vector>
Vector LU_solve(
        const std::pair<Matrix, std::vector<int>>& LUP,
        const Vector& b) {
    const auto& LU = LUP.first;
    const auto& P = LUP.second;
    const int N = P.size();
    Vector Pb(N);
    for (int i = 0; i < N; i++) {
        Pb[P[i]] = b[i];
    }

    Vector y(N);
    for (int i = 0; i < N; i++) {
        y[i] = Pb[i];
        for (int j = 0; j < i; j++) {
            y[i] -= LU[i][j]*y[j];
        }
    }

    Vector x(N);
    for (int i = N-1; i >= 0; i--) {
        x[i] = y[i]/LU[i][i];
        for (int j = i+1; j < N; j++) {
            x[i] -= LU[i][j]*y[j];
        }
    }

    return x;
}

#endif  // SRC_LU_DECOMPOSITION_HPP_
