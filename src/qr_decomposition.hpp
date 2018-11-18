#ifndef SRC_QR_DECOMPOSITION_HPP_
#define SRC_QR_DECOMPOSITION_HPP_

#include<cmath>
#include<numeric>
#include<utility>
#include<vector>

template<typename Matrix, typename Vector>
std::pair<std::vector<Vector>, Matrix>
QR_decomposition(Matrix A) {
    const int N = A.size();
    std::vector<Vector> Q;
    for (int i = 0; i < N; i++) {
        Vector v(N);
        for (int j = i; j < N; j++) {
            v[j] = A[j][i];
        }
        auto vsq = v[i]*v[i];
        for (int j = i+1; j < N; j++) {
            vsq += v[j]*v[j];
        }
        auto wsq = vsq - v[i]*v[i];
        v[i] += (v[i] > 0 ? 1 : -1) * std::sqrt(vsq);
        wsq += v[i]*v[i];
        auto wnorm = std::sqrt(wsq);

        for (int j = i; j < N; j++) {
            v[j] /= wnorm;
        }

        Q.push_back(v);
        // Update j-th column of A
        for (int j = i; j < N; j++) {
            auto va = v[i]*A[i][j];
            for (int k = i+1; k < N; k++) {
                va += v[k]*A[k][j];
            }
            for (int k = i; k < N; k++) {
                A[k][j] -= 2*v[k]*va;
            }
        }
    }

    return {Q, A};
}

template<typename Matrix, typename Vector>
Vector QR_solve(
        const std::pair<std::vector<Vector>, Matrix>& QR,
        Vector b) {
    const int N = b.size();
    const auto& Q = QR.first;
    for (int i = 0; i < N; i++) {
        auto c = Q[i][i]*b[i];
        for (int j = i+1; j < N; j++) {
            c += Q[i][j]*b[j];
        }
        for (int j = i; j < N; j++) {
            b[j] -= 2*c*b[j];
        }
    }

    const auto& R = QR.second;
    Vector x(N);
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
