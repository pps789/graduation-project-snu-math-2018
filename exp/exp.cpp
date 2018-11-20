#include<cstdio>
#include<ctime>
#include<vector>
#include"../src/cholesky_decomposition.hpp"
#include"../src/conjugate_gradient_method.hpp"
#include"../src/gauss_elimination.hpp"
#include"../src/gauss_seidel_method.hpp"
#include"../src/gmres.hpp"
#include"../src/gradient_descent.hpp"
#include"../src/jacobi_method.hpp"
#include"../src/lu_decomposition.hpp"
#include"../src/qr_decomposition.hpp"
#include"../src/sor.hpp"
using namespace std;

const int N = 16348;  // Size of matrix!

int P[N + 1];
int Find(int u) { return P[u] = (P[u] == u ? u : Find(P[u])); }
void Merge(int u, int v) {
    u = Find(u), v = Find(v);
    P[u] = v;
}

int main() {
    for (int i = 0; i <= N; i++) P[i] = i;
    vector<vector<double>> A(N + 1, vector<double>(N + 1));
    int merge_cnt = 0;

    for (int i = 0; i < N*N; i++) {
        int u = rand()%(N + 1), v = rand()%(N + 1);
        if (u != v) {
            A[u][u] += 1; A[v][v] += 1;
            A[u][v] -= 1; A[v][u] -= 1;
            if (Find(u) != Find(v)) {
                merge_cnt++;
                Merge(u, v);
            }
        }
    }
    while (merge_cnt < N) {
        int u = rand()%(N + 1), v = rand()%(N + 1);
        if (u != v) {
            A[u][u] += 1; A[v][v] += 1;
            A[u][v] -= 1; A[v][u] -= 1;
            if (Find(u) != Find(v)) {
                merge_cnt++;
                Merge(u, v);
            }
        }
    }

    A.pop_back();
    for (auto& v : A) v.pop_back();

    vector<double> b(N);
    for (int i = 0; i < 5; i++) b[i] = rand();

    const double TOL = 1e-6;

    auto t = clock();

    // cholesky
    auto cholesky_L = cholesky_decomposition(A);
    auto cholesky = cholesky_solve(cholesky_L, b);

    printf("CHOLESKY: %.6f\n", ((double)(clock()-t)) / CLOCKS_PER_SEC);
    t = clock();

    // conjgrad
    auto conj = conjugate_gradient_method(A, b, TOL);

    printf("CONJGRAD: %.6f\n", ((double)(clock()-t)) / CLOCKS_PER_SEC);
    t = clock();

    // gauss-elimination
    auto GE = gauss_elimination(A, b);

    printf("GAUSS: %.6f\n", ((double)(clock()-t)) / CLOCKS_PER_SEC);
    t = clock();

    // gauss-seidel
    auto GS = gauss_seidel_method(A, b, TOL);

    printf("GAUSS-SEIDEL: %.6f\n", ((double)(clock()-t)) / CLOCKS_PER_SEC);
    t = clock();

    // gmres
    auto GMRES = gmres(A, b, vector<double>(5), TOL);

    printf("GMRES: %.6f\n", ((double)(clock()-t)) / CLOCKS_PER_SEC);
    t = clock();

    // gradient descent
    auto GRD = gradient_descent(A, b, TOL);

    printf("GRADDESC: %.6f\n", ((double)(clock()-t)) / CLOCKS_PER_SEC);
    t = clock();

    // jacobi
    auto J = jacobi_method(A, b, TOL);

    printf("JACOBI: %.6f\n", ((double)(clock()-t)) / CLOCKS_PER_SEC);
    t = clock();

    // LU
    auto LU = LU_decomposition(A);
    auto LU_sol = LU_solve(LU, b);

    printf("LU: %.6f\n", ((double)(clock()-t)) / CLOCKS_PER_SEC);
    t = clock();

    // qr
    auto QR = QR_decomposition(A);
    auto QR_sol = QR_solve(QR, b);

    printf("QR: %.6f\n", ((double)(clock()-t)) / CLOCKS_PER_SEC);
    t = clock();

    // SOR
    auto SOR = sor(A, b, 0.5, TOL);

    printf("SOR: %.6f\n", ((double)(clock()-t)) / CLOCKS_PER_SEC);
    t = clock();

    return 0;
}
