
#include "algo.h"

using namespace std;

int idx_2D(int i, int j, int N) { return i * N + j; }

void set_gauss_conditions_2D(vector<double>& u0, int N)
{
    vector<double> gx(N, 0.0);
    for (int i = 0; i < N; i++) gx[i] = gauss(i * gridParams.dx, gaussParams.x0, gaussParams.sigma);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            bool is_boundary = (i == 0 || j == 0 || i == N - 1 || j == N - 1);
            if (is_boundary) u0[idx_2D(i, j, N)] = 0.0;
            else u0[idx_2D(i, j, N)] = gx[i] * gx[j];
        }
    }
}

void compute_laplace_2D(const vector<double>& mat, vector<double>& lap, int N)
{
    const double inv_dx2 = 1.0 / (gridParams.dx * gridParams.dx);
    std::fill(lap.begin(), lap.end(), 0.0);

    for (int i = 1; i < N - 1; i++) {
        const int row = i * N;
        for (int j = 1; j < N - 1; j++) {
            const int k = row + j;
            const double grad = mat[k + N] + mat[k - N] + mat[k + 1] + mat[k - 1] - 4.0 * mat[k];
            lap[k] = grad * inv_dx2;
        }
    }
}


void init_ci_half(vector<double>& ci_half, const vector<double>& ci0, const vector<double>& rho, vector<double>& lap, int N)
{
    const double coeff = 0.5 * gridParams.dt * 0.5 / gridParams.m;
    compute_laplace_2D(rho, lap, N);
    for (int k = 0; k < N * N; k++) ci_half[k] = ci0[k] + coeff * lap[k];
}

void solution_step(vector<double>& rho, vector<double>& ci_half, vector<double>& lap, int N)
{
    const double coeff = gridParams.dt * 0.5 / gridParams.m;
    compute_laplace_2D(ci_half, lap, N);
    for (int k = 0; k < N * N; k++) rho[k] -= coeff * lap[k];
    compute_laplace_2D(rho, lap, N);
    for (int k = 0; k < N * N; k++) ci_half[k] += coeff * lap[k];
}

void write_data_step(const vector<double>& rho, const vector<double>& ci_half, ofstream& file, int N)
{
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            const int k = idx_2D(i, j, N);
            const double rho_n = rho[k];
            const double ci_n  = ci_half[k];
            const double psi2  = rho_n * rho_n + ci_n * ci_n;
            file << i * gridParams.dx << " " << j * gridParams.dx << " " << psi2 << "\n";
        }
    }
    file << "\n";
}

void solver(vector<double>& rhoMat, vector<double>& ciHalfMat, vector<double>& lapMat, int N, int nt) {
    ofstream file("dat/solution_shrodinger_leapfrog_2D.dat");
    write_data_step(rhoMat, ciHalfMat, file, N);
    for (int t = 1; t < nt; t++) {
        cout << "Time step: " << t << "/" << nt << "\r" << std::flush;
        solution_step(rhoMat, ciHalfMat, lapMat, N);
        if (t % 50 == 0) write_data_step(rhoMat, ciHalfMat, file, N);
    }
    file.close();
}

void leapfrog_solver_2D()
{
    const int nt = gridParams.Nt, N = gridParams.N + 1, size = N * N;
    vector<double> rhoMat(size, 0.0), ciHalfMat(size, 0.0), ci0Mat(size, 0.0), lapMat(size, 0.0);

    set_gauss_conditions_2D(rhoMat, N);
    init_ci_half(ciHalfMat, ci0Mat, rhoMat, lapMat, N);

    solver(rhoMat, ciHalfMat, lapMat, N, nt);
}