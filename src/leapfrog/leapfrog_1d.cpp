#include "algo.h"

int idx(int i, int n) { return n * (gridParams.N_x + 1) + i; }

double laplace_1d(const std::vector<double>& vec, int i, int n) {
    if (i == 0 || i == gridParams.N_x) return 0.0;
    double grad = vec[idx(i + 1, n)] - 2.0 * vec[idx(i, n)] + vec[idx(i - 1, n)];
    return grad / (gridParams.dx * gridParams.dx);
}

void init_ci_half(std::vector<double>& ci_half, const std::vector<double>& ci0, const std::vector<double>& rho)
{
    for (int i = 0; i <= gridParams.N_x; i++) {
        double lap0 = laplace_1d(rho, i, 0);
        ci_half[idx(i, 0)] = ci0[i] + 0.5 * gridParams.dt * lap0 * 0.5 / gridParams.m;
    }
}

void solution_time_step(std::vector<double>& rho, std::vector<double>& ci_half, int n)
{
    for (int i = 0; i <= gridParams.N_x; i++) {
        double lapCi = laplace_1d(ci_half, i, n);
        double v = potential(i * gridParams.dx);
        rho[idx(i, n + 1)] = rho[idx(i, n)]
            - gridParams.dt * lapCi * 0.5 / gridParams.m
            + gridParams.dt * v * ci_half[idx(i, n)];
    }

    if (n + 1 < gridParams.N_steps) {
        for (int i = 0; i <= gridParams.N_x; i++) {
            double lapRho = laplace_1d(rho, i, n + 1);
            double v = potential(i * gridParams.dx);
            ci_half[idx(i, n + 1)] = ci_half[idx(i, n)]
                + gridParams.dt * lapRho * 0.5 / gridParams.m
                - gridParams.dt * v * rho[idx(i, n + 1)];
        }
    }
}

void set_gauss_conditions(std::vector<double>& u0)
{
    for (int k = 0; k <= gridParams.N_x; k++) {
        if (k == 0 || k == gridParams.N_x) u0[idx(k, 0)] = 0;
        else u0[idx(k, 0)] = gauss(k * gridParams.dx, gaussParams.x0, gaussParams.sigma);
    }
}

void write_data_psi(const std::vector<double>& rho, const std::vector<double>& ci_half)
{
    std::ofstream file("dat/solution_shrodinger_leapfrog.dat");

    for (int n = 1; n < gridParams.N_steps; n++) {
        for (int i = 0; i <= gridParams.N_x; i++) {
            double rho_n = rho[idx(i, n)];
            double ci_n  = 0.5 * (ci_half[idx(i, n - 1)] + ci_half[idx(i, n)]);
            double psi2  = rho_n * rho_n + ci_n * ci_n;
            file << psi2 << " ";
        }
        file << "\n";
    }
    file.close();
}

void leapfrog_solver()
{
    double cfl = gridParams.dt / (gridParams.dx * gridParams.dx);
    if (cfl > 0.5 / gridParams.m) {
        std::cerr << "WARNING: CFL violated! cfl=" << cfl 
                << " limit=" << 0.5 / gridParams.m 
                << " -> reduce dt or increase dx\n";
        return;
    }

    std::vector<double> rho((gridParams.N_x + 1) * (gridParams.N_steps + 1), 0.0);
    std::vector<double> ci_half((gridParams.N_x + 1) * gridParams.N_steps, 0.0);
    std::vector<double> ci0(gridParams.N_x + 1, 0.0);

    set_gauss_conditions(rho);
    init_ci_half(ci_half, ci0, rho);

    for (int n = 0; n < gridParams.N_steps; ++n)
        solution_time_step(rho, ci_half, n);

    write_data_psi(rho, ci_half);
}