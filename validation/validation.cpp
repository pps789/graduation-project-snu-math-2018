#include<cstdio>
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

void print_vector(vector<double> v) {
    for (double x : v) printf("%.6f ", x);
}

int main() {
    // Example matrix of https://en.wikipedia.org/wiki/Laplacian_matrix
    // Drop the last row/col ->  non-singular!
    vector<vector<double>> A(6, vector<double>(6));
    vector<pair<int, int>> edges = {{1, 2}, {1, 5}, {2, 3}, {2, 5}, {3, 4}, {4, 5}, {4, 6}};
    for (auto edge : edges) {
        int u = edge.first, v = edge.second;
        u--; v--;
        A[u][u] += 1; A[v][v] += 1;
        A[u][v] -= 1; A[v][u] -= 1;
    }

    A.pop_back();
    for (auto& v : A) v.pop_back();

    vector<double> b(5);
    for (int i = 0; i < 5; i++) b[i] = i+1;  // {1,2,3,4,5}

    const double TOL = 1e-6;

    // Solution is : 21.909091 21.636364 19.818182 15.000000 21.181818

    // cholesky
    auto cholesky_L = cholesky_decomposition(A);
    auto cholesky = cholesky_solve(cholesky_L, b);
    printf("CHOLESKY: "); print_vector(cholesky); printf("\n");

    // conjgrad
    auto conj = conjugate_gradient_method(A, b, TOL);
    printf("CONJGRAD: "); print_vector(conj); printf("\n");

    // gauss-elimination
    auto GE = gauss_elimination(A, b);
    printf("GE: "); print_vector(GE); printf("\n");

    // gauss-seidel
    auto GS = gauss_seidel_method(A, b, TOL);
    printf("GS: "); print_vector(GS); printf("\n");

    // gmres
    auto GMRES = gmres(A, b, vector<double>(5), TOL);
    printf("GMRES: "); print_vector(GMRES); printf("\n");

    // gradient descent
    auto GRD = gradient_descent(A, b, TOL);
    printf("GRD: "); print_vector(GRD); printf("\n");

    // jacobi
    auto J = jacobi_method(A, b, TOL);
    printf("J: "); print_vector(J); printf("\n");

    // LU
    auto LU = LU_decomposition(A);
    auto LU_sol = LU_solve(LU, b);
    printf("LU: "); print_vector(LU_sol); printf("\n");

    // qr
    auto QR = QR_decomposition(A);
    auto QR_sol = QR_solve(QR, b);
    printf("QR: "); print_vector(QR_sol); printf("\n");

    // SOR
    auto SOR = sor(A, b, 0.5, TOL);
    printf("SOR: "); print_vector(SOR); printf("\n");

    return 0;
}
