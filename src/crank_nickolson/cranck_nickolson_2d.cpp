#include "algo.h"

using namespace Eigen;
using namespace std;
using namespace std::complex_literals;

MatrixXcd laplace_2d(int n);

void set_gauss_conditions_2d(Eigen::VectorXcd& u0)
{
    int d = gridParams2d.N + 1;
    for (int i = 0; i < d; i++) {
        for (int j = 0; j < d; j++) {
            if (i == 0 || j == 0 || i == gridParams2d.N || j == gridParams2d.N) u0(i + j * d) = 0.0;
            else {
                double x = i * gridParams2d.dx;
                double y = j * gridParams2d.dx;
                u0(i + j * d) = gauss(x, gaussParams.x0, gaussParams.sigma) * gauss(y, gaussParams.x0, gaussParams.sigma) + 0i;
                // u0(i + j * d) = x + 1i*y;
            }
        }
    }
}

MatrixXcd laplace_2d(int n){
    int d = pow(n,2); // dimension de la matrice aplattie
    MatrixXcd laplace_x = MatrixXcd::Zero(d,d);
    MatrixXcd laplace_y = MatrixXcd::Zero(d,d);
    MatrixXcd id = MatrixXcd::Identity(n,n);
    MatrixXcd tridiag = MatrixXcd::Zero(n,n);
    // boucle principale
    for (int k = 0; k < n-1; k++){
        tridiag(k,k) = -2;
        tridiag(k,k+1) = 1;
        tridiag(k+1,k) = 1;
    }
    tridiag(n-1,n-1) = -2;
    // matrice laplace_x est simplement matrice diagonale par blocs de matrices tridiagonales
    for (int k = 0; k <= d-n; k=k+n){
        laplace_x.block(k,k,n,n) = 1*tridiag;
    }
    // matrice laplace_y est une matrice tridiagonale par blocs
    for (int k = 0; k < n-1; k++){
        laplace_y.block(n*k,n*k,n,n) = -2*id;
        laplace_y.block(n*k,n*k+n,n,n) = 1*id;
        laplace_y.block(n*k+n,n*k,n,n) = 1*id;
    }
    laplace_y.block(d-n,d-n,n,n) = -2*id;
    return laplace_x + laplace_y;
}

void time_step_solution_2d(VectorXcd& u, const MatrixXcd& M, bool withWrite)
{
    u = M * u;
    if (withWrite) {
        ecrit("dat/cranck_nickolson_2d.dat", u);
    }
}

void solver_2d(VectorXcd& u, const MatrixXcd& H_m, const MatrixXcd& H_p, const MatrixXcd& H, int N_steps)
{
    ecrit_new("dat/cranck_nickolson_2d.dat", u);
    //ecrit_energy("dat/shrodinger_energy.dat", u.adjoint()*H*u, true);
    (void) H;
    cout << H_p.size() << endl;
    MatrixXcd M = H_p.inverse() * H_m;
    for (int n = 0; n < N_steps; n++){
        time_step_solution_2d(u, M, true);
        if(n%50==0){
            cout << n << endl;
        }
    }
}

void cranck_nickolson_solver_2d(){
    VectorXcd u0 = VectorXcd((gridParams2d.N + 1) * (gridParams2d.N + 1));
    MatrixXcd l2d = laplace_2d(gridParams2d.N);
    cout << gridParams2d.N << endl;
    set_gauss_conditions_2d(u0);
    ecrit_new("dat/cranck_nickolson_2d.dat", u0);
    //cout << u << endl;
    //cout << u.adjoint()*l2d*u << endl;
    //cout << laplace_2d(gridParams2d.N) << endl;
    VectorXcd unew = u0;
    MatrixXcd I = MatrixXcd::Identity((gridParams2d.N + 1)*(gridParams2d.N + 1), (gridParams2d.N + 1)*(gridParams2d.N + 1));
    MatrixXcd D2 = laplace_2d(gridParams2d.N + 1) / (gridParams2d.dx * gridParams2d.dx);
    MatrixXcd H  = - (1 / (2.0 * gridParams.m)) * D2; // + V à définir
    MatrixXcd A = I + 0.5i * gridParams2d.dt * H; // Matrice du membre de gauche du systeme lineaire
    MatrixXcd B = I - 0.5i * gridParams2d.dt * H; // Matrice du membre de droite du systeme lineaire

    //A.row(0).setZero(); A.row(N).setZero();
    //B.row(0).setZero(); B.row(N).setZero();
    //A(0, 0) = 1.0; A(N, N) = 1.0;
    //B(0, 0) = 1.0; B(N, N) = 1.0;

    solver_2d(unew, A, B, H, gridParams.N_steps);
    cout << unew(0) << endl;
}

