#ifndef SRC_GAUSS_ELIMINATION_HPP_
#define SRC_GAUSS_ELIMINATION_HPP_

#include<algorithm>
#include<utility>
#include"helper.hpp"

// Assume that Matrix is square.
template<typename Matrix, typename Vector>
Vector gauss_elimination(Matrix A, Vector b) {
    const int N = A.size();
    for (int i = 0; i < N; i++) {
        int t = i;
        auto maxval = std::abs(A[i][i]);
        for (int j = i+1; j < N; j++) {
            if (std::abs(A[j][i]) > maxval) {
                t = j;
                maxval = std::abs(A[j][i]);
            }
        }
        std::swap(A[i], A[t]); std::swap(b[i], b[t]);
        if (nonzero(A[i][i])) {
            for (int j = i+1; j < N; j++) {
                auto rate = A[j][i] / A[i][i];
                for (int k = i; k < N; k++) A[j][k] -= A[i][k]*rate;
                b[j] -= b[i]*rate;
            }
        }
    }

    for (int i = N-1; i >= 0; i--) {
        if (nonzero(A[i][i])) {
            b[i] /= A[i][i];
            for (int j = i-1; j >= 0; j--) b[j] -= A[j][i]*b[i];
        } else if (nonzero(b[i])) {
            // Singular case.
            return Vector();
        }
    }
    return b;
}

#endif  // SRC_GAUSS_ELIMINATION_HPP_
