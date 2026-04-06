#include "algo.h"

static int idx(int i, int n) { return n * gridParams.N + i; }

double laplace_phi(const std::vector<double>& phi, int i, int n) {
    if (i == 0 || i == gridParams.N - 1) return 0.0;
    return (phi[idx(i + 1, n)] - 2.0 * phi[idx(i, n)] + phi[idx(i - 1, n)]) / (gridParams.dx * gridParams.dx);
}

void set_pi_half(std::vector<double>& phi, std::vector<double>& pi_half, const std::vector<double>& pi0, double dt)
{
    for (int i = 0; i < gridParams.N; i++) {
        double lap0 = laplace_phi(phi, i, 0);
        pi_half[idx(i, 0)] = pi0[i] + 0.5 * dt * (lap0 - gridParams.m * gridParams.m * phi[idx(i, 0)]);
    }
}

void set_gauss_conditions(std::vector<double>& phi)
{
    for(int i = 0; i < gridParams.N; i++) {
        double x_i = i * gridParams.dx;
        phi[idx(i, 0)] = gauss(x_i, gaussParams.x0, gaussParams.sigma);
    }
}

void solution_step(std::vector<double>& phi, std::vector<double>& pi_half, int n)
{
    for (int i = 0; i < gridParams.N; i++) {
        double lap = laplace_phi(phi, i, n);
        double pi_prev_half = (n == 0) ? pi_half[idx(i, 0)] : pi_half[idx(i, n - 1)];
        double pi_new_half  = pi_prev_half + gridParams.dt * (lap - gridParams.m * gridParams.m * phi[idx(i, n)]);
        pi_half[idx(i, n)] = pi_new_half;
        phi[idx(i, n + 1)] = phi[idx(i, n)] + gridParams.dt * pi_new_half;
    }
}

void write_data(const std::vector<double>& phi, const std::vector<double>& pi_half)
{
    std::ofstream file("dat/klein_gordon_output.dat");
    for (int n = 0; n < gridParams.Nt; n++) {
        for (int i = 0; i < gridParams.N; i++) {
            double phi_sq = phi[idx(i, n)] * phi[idx(i, n)];
            double pi_sq  = pi_half[idx(i, n)] * pi_half[idx(i, n)];
            file << phi_sq + pi_sq << " ";
        }
        file << "\n";
    }
    file.close();
}

void klein_gordon_1D() {
    std::vector<double> phi(gridParams.N * (gridParams.Nt + 1), 0.0);
    std::vector<double> pi_half(gridParams.N * gridParams.Nt, 0.0);
    std::vector<double> pi0(gridParams.N, 0.0);
    
    set_gauss_conditions(phi);
    set_pi_half(phi, pi_half, pi0, gridParams.dt);

    for (int n = 0; n < gridParams.Nt; ++n) solution_step(phi, pi_half, n);
        
    write_data(phi, pi_half);
}
