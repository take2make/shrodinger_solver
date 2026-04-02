#include "algo.h"

using namespace std;

int idx(int i, int n) { return n * (gridParams.N + 1) + i; }

double laplace_1D(const vector<double>& vec, int i, int n) {
    if (i == 0 || i == gridParams.N) return 0.0;
    double dvec = vec[idx(i + 1, n)] - 2.0 * vec[idx(i, n)] + vec[idx(i - 1, n)];
    return dvec / (gridParams.dx * gridParams.dx);
}

void init_ci_half(vector<double>& ci_half, const vector<double>& ci0, const vector<double>& rho)
{
    for (int i = 0; i <= gridParams.N; i++) {
        double lap0 = laplace_1D(rho, i, 0);
        ci_half[idx(i, 0)] = ci0[i] + 0.5 * gridParams.dt * lap0 * 0.5 / gridParams.m;
    }
}

void solution_time_step(vector<double>& rho, vector<double>& ci_half, int n)
{
    for (int i = 0; i <= gridParams.N; i++) {
        double lapCi = laplace_1D(ci_half, i, n);
        double v = potential(i * gridParams.dx);
        rho[idx(i, n + 1)] = rho[idx(i, n)]
            - gridParams.dt * lapCi * 0.5 / gridParams.m
            + gridParams.dt * v * ci_half[idx(i, n)];
    }

    if (n + 1 < gridParams.Nt) {
        for (int i = 0; i <= gridParams.N; i++) {
            double lapRho = laplace_1D(rho, i, n + 1);
            double v = potential(i * gridParams.dx);
            ci_half[idx(i, n + 1)] = ci_half[idx(i, n)]
                + gridParams.dt * lapRho * 0.5 / gridParams.m
                - gridParams.dt * v * rho[idx(i, n + 1)];
        }
    }
}

void set_gauss_conditions(vector<double>& u0, vector<double>& ci0, double k)
{
    for (int i = 0; i <= gridParams.N; i++) {
        bool is_boundary = (i == 0 || i == gridParams.N);
        if (is_boundary) { u0[idx(i, 0)] = 0; ci0[i] = 0; }
        else {
            u0[idx(i, 0)] = gauss_with_velocity_real(i * gridParams.dx, gaussParams.x0, gaussParams.sigma, k);
            ci0[i] = gauss_with_velocity_imag(i * gridParams.dx, gaussParams.x0, gaussParams.sigma, k);
        }
    }
}

void write_data_psi(const vector<double>& rho, const vector<double>& ci_half)
{
    ofstream file("dat/solution_shrodinger_leapfrog.dat");

    for (int n = 1; n < gridParams.Nt; n++) {
        for (int i = 0; i <= gridParams.N; i++) {
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
    int N = gridParams.N + 1, nt = gridParams.Nt;
    double k = gaussParams.k;
    if (cfl > 0.5 / gridParams.m) {
        std::cerr << "WARNING: CFL violated! cfl=" << cfl << " limit=" << 0.5 / gridParams.m  << " -> reduce dt or increase dx\n";
        return;
    }

    vector<double> rho(N * (nt + 1), 0.0), ci_half(N * nt, 0.0), ci0(N, 0.0);

    set_gauss_conditions(rho, ci0, k);
    init_ci_half(ci_half, ci0, rho);

    for (int n = 0; n < nt; ++n)
        solution_time_step(rho, ci_half, n);

    write_data_psi(rho, ci_half);
}