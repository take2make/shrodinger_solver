import matplotlib.pyplot as plt
import numpy as np

dat = "dat/solution_temporelle_shrodinger.dat"
#dat = "dat/solution_shrodinger_leapfrog.dat"
L = 10.0
M = 0.6
DT = 0.001
X0 = 5.0
SIGMA = 1.0
K = 0.0
N_MAX = 200
STRIDE = 1
T_IDX = 3000

def load_data(path: str) -> np.ndarray:
    return np.loadtxt(path)

def gauss(x, x0, sigma):
    coeff = 1.0 / (np.sqrt(2.0 * np.pi) * sigma)
    return coeff * np.exp(-((x - x0) ** 2) / (2.0 * sigma ** 2))

def get_analytic_solution(x, t, num_norm):
    dx = x[1] - x[0]
    psi0 = gauss(x, X0, SIGMA) * np.exp(1j * K * x)
    ana_norm = np.sum(np.abs(psi0) ** 2) * dx
    if ana_norm > 0.0:
        psi0 *= np.sqrt(num_norm / ana_norm)

    n = np.arange(1, N_MAX + 1, dtype=float)
    basis = np.sqrt(2.0 / L) * np.sin(np.outer(n, np.pi * x / L))
    Cn = np.trapezoid(psi0 * basis, x, axis=1)
    En = (n * np.pi / L) ** 2 / (2.0 * M)
    phase = np.exp(-1j * En * t)
    psi_t = (Cn * phase) @ basis
    return psi_t

def main():
    data = load_data(dat)
    nt, nx = data.shape
    x = np.linspace(0.0, L, nx)
    dx = x[1] - x[0]
    num_norm = np.sum(data[0]) * dx
    t_idx = np.arange(0, nt, STRIDE)
    times = t_idx * DT
    psi2_analytic = np.empty((t_idx.size, nx), dtype=float)
    for i, t in enumerate(times):
        psi_t = get_analytic_solution(x, t, num_norm)
        psi2_analytic[i] = np.abs(psi_t) ** 2

    data_slice = data[t_idx]
    max_err = np.max(np.abs(data_slice - psi2_analytic), axis=1)
    l2_err = np.sqrt(np.sum((data_slice - psi2_analytic) ** 2, axis=1) * dx)

    prob_num = np.sum(data_slice, axis=1) * dx
    prob_ana = np.sum(psi2_analytic, axis=1) * dx

    profile_step = t_idx[T_IDX]
    profile_time = times[T_IDX]

    fig, axs = plt.subplots(1, 3, figsize=(15, 4))

    axs[0].plot(x, data_slice[T_IDX], label=f"numeric t={profile_time:.3g}")
    axs[0].plot(
        x, psi2_analytic[T_IDX], "--", label=f"analytic t={profile_time:.3g}"
    )
    axs[0].set_xlabel("x")
    axs[0].set_ylabel(r"$|\psi|^2$")
    axs[0].set_title(f"Profile at step {profile_step}")
    axs[0].legend()
    axs[0].grid()

    axs[1].plot(times, max_err, label="max error")
    axs[1].set_xlabel("time")
    axs[1].set_ylabel("error")
    axs[1].set_title("Analytic vs numeric")
    axs[1].legend()
    axs[1].grid()

    axs[2].plot(times, prob_num, label="numeric")
    axs[2].plot(times, prob_ana, label="analytic")
    axs[2].set_xlabel("time")
    axs[2].set_ylabel("probability")
    axs[2].set_title("Norm conservation")
    axs[2].legend()
    axs[2].grid()
    axs[2].set_ylim(0.9 * num_norm, 1.1 * num_norm)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
