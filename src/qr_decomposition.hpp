#ifndef SRC_QR_DECOMPOSITION_HPP_
#define SRC_QR_DECOMPOSITION_HPP_

#include<cmath>
#include<numeric>
#include<utility>
#include<vector>
#include<cstdio>

// Input matrix can be ractangle.
template<typename T>
std::pair<std::vector<std::vector<T>>, std::vector<std::vector<T>>>
QR_decomposition(std::vector<std::vector<T>> A) {
    // M: #(row), N: #(col)
    const int M = A.size();
    const int N = A.back().size();
    // Note that M >= N.

    // Q is not a matrix; it is a set of column vector.
    std::vector<std::vector<T>> Q;
    for (int i = 0; i < N; i++) {
        std::vector<T> v(M);
        for (int j = i; j < M; j++) {
            v[j] = A[j][i];
        }
        T vsq = 0;
        for (int j = i; j < M; j++) {
            vsq += v[j]*v[j];
        }
        T wsq = vsq - v[i]*v[i];
        v[i] += (v[i] > 0 ? 1 : -1) * std::sqrt(vsq);
        wsq += v[i]*v[i];
        T wnorm = std::sqrt(wsq);

        for (int j = i; j < M; j++) {
            v[j] /= wnorm;
        }

        Q.push_back(v);
        // Update j-th column of A
        for (int j = i; j < N; j++) {
            T va = 0;
            for (int k = i; k < M; k++) {
                va += v[k]*A[k][j];
            }
            for (int k = i; k < M; k++) {
                A[k][j] -= 2*v[k]*va;
            }
        }
    }

    // drop the "0" rows.
    A.resize(N);
    return {Q, A};
}

template<typename T>
std::vector<T> QR_solve(
        const std::pair<
            std::vector<std::vector<T>>,
            std::vector<std::vector<T>>>& QR,
        std::vector<T> b) {
    const auto& Q = QR.first;
    // Q: M by N; but column-wise
    const int N = Q.size();
    const int M = b.size();
    for (int i = 0; i < N; i++) {
        T c = 0;
        for (int j = i; j < M; j++) {
            c += Q[i][j]*b[j];
        }
        for (int j = i; j < M; j++) {
            b[j] -= 2*c*Q[i][j];
        }
    }

    const auto& R = QR.second;
    std::vector<T> x(N);
    for (int i = N-1; i >= 0; i--) {
        x[i] = b[i];
        for (int j = i+1; j < N; j++) {
            x[i] -= R[i][j]*x[j];
        }
        x[i] /= R[i][i];
    }

    return x;
}

#endif  // SRC_QR_DECOMPOSITION_HPP_
