#ifndef SRC_CHOLESKY_DECOMPOSITION_HPP_
#define SRC_CHOLESKY_DECOMPOSITION_HPP_

#include<algorithm>
#include<cmath>
#include<vector>
#include"helper.hpp"

// Assume that Matrix is square, SPD.
// Returns L.
template<typename Matrix, typename Vector>
Matrix cholesky_decomposition(Matrix A) {
    const int N = A.size();
    Matrix L = A;

    // Initialize.
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) L[i][j] = 0;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j <= i; j++) {
            auto val = A[i][j];
            for (int k = 0; k < j; k++) {
                val -= L[i][k]*L[j][k];
            }
            if (i == j) {
                L[i][j] = std::sqrt(val);
            } else {  // so that j < i.
                L[i][j] = val / L[j][j];
            }
        }
    }

    return L;
}

template<typename Matrix, typename Vector>
Vector cholesky_solve(Matrix L, Vector b) {
    const int N = L.size();
    Vector y(N);

    for (int i = 0; i < N; i++) {
        y[i] = b[i];
        for (int j = 0; j < i; j++) {
            y[i] -= L[i][j]*y[j];
        }
        y[i] /= L[i][i];
    }

    Vector x(N);
    for (int i = N-1; i >= 0; i--) {
        x[i] = y[i];
        for (int j = i+1; j < N; j++) {
            x[i] -= L[j][i]*x[j];
        }
        x[i] /= L[i][i];
    }

    return x;
}

#endif  // SRC_CHOLESKY_DECOMPOSITION_HPP_
