#include "algo.h"

using namespace Eigen;
using namespace std;

void time_step_solution(VectorXcd& u, const MatrixXcd& M)
{
    u = M * u;
    ecrit("dat/solution_temporelle_shrodinger.dat", u);
}

void solver(VectorXcd& u, const MatrixXcd& H_m, const MatrixXcd& H_p, int Nt)
{
    ecrit_new("dat/solution_temporelle_shrodinger.dat", u);
    MatrixXcd M = H_p.inverse() * H_m;
    for (int n = 0; n < Nt; n++) time_step_solution(u, M);
}

void set_gauss_conditions(Eigen::VectorXcd& u0)
{
    for (int k = 0; k <= gridParams.N; k++) {
        bool is_boundary = (k == 0 || k == gridParams.N);
        if (is_boundary) {
            u0[k] = 0.0;
        } else {
            double x = k * gridParams.dx;
            double re = gauss_with_velocity_real(x, gaussParams.x0, gaussParams.sigma, gaussParams.k);
            double im = gauss_with_velocity_imag(x, gaussParams.x0, gaussParams.sigma, gaussParams.k);
            u0[k] = std::complex<double>(re, im);
        }
    }
}

Eigen::MatrixXcd laplace_1d(int n)
{
    Eigen::MatrixXcd mat = Eigen::MatrixXcd::Zero(n, n);
    for (int k = 0; k < n - 1; k++) {
        mat(k, k) = -2.0;
        mat(k, k + 1) = 1.0;
        mat(k + 1, k) = 1.0;
    }
    mat(n - 1, n - 1) = -2.0;
    return mat;
}

void cranck_nickolson_solver() {
    double dt = gridParams.dt; // discretization du temps
    double dx = gridParams.dx; // discretisation de l'espace
    int N = gridParams.N; // nombre de divisions

    VectorXcd u0 = VectorXcd(N + 1); // configuration initiale, obtenue en diagonalisant le Hamiltonien
    set_gauss_conditions(u0);

    MatrixXcd V = Eigen::MatrixXcd::Zero(N + 1, N + 1);
    for(int j = 0; j <= N; j++) V(j,j) = potential(j * dx);
    
    VectorXcd unew = u0;
    MatrixXcd I = MatrixXcd::Identity(N + 1, N + 1);
    MatrixXcd D2 = laplace_1d(N + 1) / (dx * dx);
    MatrixXcd H  = - 1 / (2.0 * gridParams.m) * D2 + V;
    MatrixXcd A = I + 0.5i * dt * H; // Matrice du membre de gauche du systeme lineaire
    MatrixXcd B = I - 0.5i * dt * H; // Matrice du membre de droite du systeme lineaire

    A.row(0).setZero(); A.row(N).setZero();
    B.row(0).setZero(); B.row(N).setZero();
    A(0, 0) = 1.0; A(N, N) = 1.0;
    B(0, 0) = 1.0; B(N, N) = 1.0;

    solver(unew, A, B, gridParams.Nt);
}
