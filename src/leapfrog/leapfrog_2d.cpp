#include "algo.h"

using namespace Eigen;

void set_gauss_conditions_2D(Eigen::MatrixXd& u0)
{
    for (int i = 0; i <= gridParams.N_x; i++) {
        for (int j = 0; j <= gridParams.N_x; j++) {
            if (i == 0 || j == 0 || i == gridParams.N_x || j == gridParams.N_x) u0(i, j) = 0.0;
            else {
                double x = i * gridParams.dx;
                double y = j * gridParams.dx;
                u0(i, j) = gauss(x, gaussParams.x0, gaussParams.sigma) * gauss(y, gaussParams.x0, gaussParams.sigma);
            }
        }
    }
}

double laplace_2D(const MatrixXd& mat, int i, int j) {
    if (i == 0 || i == gridParams.N_x || j == 0 || j == gridParams.N_x) return 0.0;
    double grad = mat(i + 1, j) + mat(i - 1, j) + mat(i, j + 1) + mat(i, j - 1) - 4.0 * mat(i, j);
    return grad / (gridParams.dx * gridParams.dx);
}


void init_ci_half(MatrixXd& ci_half, const MatrixXd& ci0, const MatrixXd& rho)
{
    for (int i = 0; i <= gridParams.N_x; i++) {
        for (int j = 0; j <= gridParams.N_x; j++) {
            double lap0 = laplace_2D(rho, i, j);
            ci_half(i, j) = ci0(i, j) + 0.5 * gridParams.dt * lap0 * 0.5 / gridParams.m;
        }
    }
}

void solution_step(MatrixXd& rho, MatrixXd& ci_half)
{
    for (int i = 0; i <= gridParams.N_x; i++) {
        for (int j = 0; j <= gridParams.N_x; j++) {
            double lapCi = laplace_2D(ci_half, i, j);
            rho(i, j) -= gridParams.dt * lapCi * 0.5 / gridParams.m;
        }
    }

    for (int i = 0; i <= gridParams.N_x; i++) {
        for (int j = 0; j <= gridParams.N_x; j++) {
            double lapRho = laplace_2D(rho, i, j);
            ci_half(i, j) += gridParams.dt * lapRho * 0.5 / gridParams.m;
        }
    }
}

void write_data_step(const MatrixXd& rho, const MatrixXd& ci_half, std::ofstream& file)
{
    for (int i = 0; i <= gridParams.N_x; i++) {
        for (int j = 0; j <= gridParams.N_x; j++) {
            double rho_n = rho(i, j);
            double ci_n  = ci_half(i, j); // TODO: adjust for time step
            double psi2  = rho_n * rho_n + ci_n * ci_n;
            file << i * gridParams.dx << " " << j * gridParams.dx << " " << psi2 << "\n";
        }
    }
    file << "\n";
}

void leapfrog_solver2()
{
    const int nt = gridParams.N_steps;
    const int n = gridParams.N_x + 1;
    MatrixXd rhoMat = MatrixXd::Zero(n, n);
    MatrixXd ciHalfMat = MatrixXd::Zero(n, n);
    MatrixXd ci0Mat = MatrixXd::Zero(n, n);

    set_gauss_conditions_2D(rhoMat);
    init_ci_half(ciHalfMat, ci0Mat, rhoMat);

    std::ofstream file("dat/solution_shrodinger_leapfrog_2D.dat");
    write_data_step(rhoMat, ciHalfMat, file);
    for (int t = 1; t < nt; t++) {
        std::cout << "Time step: " << t << "/" << nt << "\r" << std::flush;
        solution_step(rhoMat, ciHalfMat);
        //write_data_step(rhoMat, ciHalfMat, file);
    }
    file.close();
    //init_ci_half(ciHalfMat, ci0Mat, rhoMat);
}