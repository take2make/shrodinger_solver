import matplotlib.pyplot as plt
import numpy as np

dat = "dat/solution_shrodinger_leapfrog_2D.dat"
L = 10.0
M = 0.6
DT = 0.0006
STRIDE = 50
SIGMA = 0.5
N_MAX = 60
T_IDX = 50


def load_data_2d(path):
    blocks = []
    current = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                if current:
                    blocks.append(np.array(current, dtype=float))
                    current = []
                continue
            parts = line.split()
            if len(parts) >= 3:
                current.append([float(parts[0]), float(parts[1]), float(parts[2])])
    if current:
        blocks.append(np.array(current, dtype=float))

    frames = []
    xs = ys = None
    for block in blocks:
        xs = np.unique(block[:, 0])
        ys = np.unique(block[:, 1])
        nx = xs.size
        ny = ys.size
        frame = block[:, 2].reshape(ny, nx)
        frames.append(frame)
    return np.stack(frames, axis=0), xs, ys


def gauss2d(x, y, x0, y0, sigma):
    coeff = 1.0 / (2.0 * np.pi * sigma * sigma)
    return coeff * np.exp(-((x - x0) ** 2 + (y - y0) ** 2) / (2.0 * sigma ** 2))


def get_analytic_solution_2d(t, basis_x, basis_y, cnm, enm):
    phase = np.exp(-1j * enm * t)
    psi_t = basis_x @ (cnm * phase) @ basis_y.T
    return psi_t

def plot_numeric_analytic(data, psi2_analytic, xs, ys, t_idx):
    num_slice = data[t_idx]
    ana_slice = psi2_analytic[t_idx]
    num_max = np.max(num_slice)
    ana_max = np.max(ana_slice)
    num_norm = num_slice / num_max if num_max != 0 else num_slice
    ana_norm = ana_slice / ana_max if ana_max != 0 else ana_slice
    fig, axs = plt.subplots(1, 2, figsize=(10, 4))
    im0 = axs[0].imshow(
        num_norm,
        origin="lower",
        extent=(xs[0], xs[-1], ys[0], ys[-1]),
        cmap="inferno",
        aspect="auto",
        vmin=0,
        vmax=1,
    )
    axs[0].set_title(f"Numeric |psi|^2")
    axs[0].set_xlabel("x")
    axs[0].set_ylabel("y")
    fig.colorbar(im0, ax=axs[0])

    im1 = axs[1].imshow(
        ana_norm,
        origin="lower",
        extent=(xs[0], xs[-1], ys[0], ys[-1]),
        cmap="inferno",
        aspect="auto",
        vmin=0,
        vmax=1,
    )
    axs[1].set_title(f"Analytic |psi|^2")
    axs[1].set_xlabel("x")
    axs[1].set_ylabel("y")
    fig.colorbar(im1, ax=axs[1])

    plt.tight_layout()
    plt.savefig("compare_numeric_analytic.png", dpi=300)
    plt.show()

def plot_error_and_prob(times, max_err, prob_num, prob_ana):
    _, axs = plt.subplots(1, 2, figsize=(10, 4))

    axs[0].plot(times, max_err, label="Max error")
    axs[0].set_xlabel("Time")
    axs[0].set_ylabel("Numeric-analytic error")
    axs[0].set_title("Analytic vs numeric")
    axs[0].legend()
    axs[0].grid()
    axs[0].set_ylim(-0.005, 0.005)

    axs[1].plot(times, prob_num, label="Numeric probability")
    axs[1].plot(times, prob_ana, label="Analytic probability", linestyle="--", color="red")
    axs[1].set_xlabel("Time")
    axs[1].set_ylabel("Probability")
    axs[1].set_title("Norm conservation")
    axs[1].legend()
    axs[1].grid()
    axs[1].set_ylim(0.997 * prob_ana[0], 1.003 * prob_ana[0])
    plt.tight_layout()
    plt.savefig("error_and_prob.png", dpi=300)
    plt.show()

def main():
    data, xs, ys = load_data_2d(dat)
    nt, ny, nx = data.shape
    dx = xs[1] - xs[0]
    dy = ys[1] - ys[0]
    x0 = L * 0.5; y0 = L * 0.5

    num_norm = np.sum(data[0]) * dx * dy
    X, Y = np.meshgrid(xs, ys, indexing="xy")
    psi0 = gauss2d(X, Y, x0, y0, SIGMA)
    ana_norm = np.sum(np.abs(psi0) ** 2) * dx * dy
    if ana_norm > 0.0:
        psi0 *= np.sqrt(num_norm / ana_norm)

    n_max_x = min(N_MAX, nx - 2)
    n_max_y = min(N_MAX, ny - 2)
    nxm = np.arange(1, n_max_x + 1, dtype=float)
    nym = np.arange(1, n_max_y + 1, dtype=float)

    basis_x = np.sqrt(2.0 / L) * np.sin(np.outer(xs, np.pi * nxm / L))
    basis_y = np.sqrt(2.0 / L) * np.sin(np.outer(ys, np.pi * nym / L))

    wx = np.ones_like(xs)
    wy = np.ones_like(ys)
    wx[0] = wx[-1] = 0.5
    wy[0] = wy[-1] = 0.5
    w2d = wy[:, None] * wx[None, :]

    cnm = basis_x.T @ (psi0 * w2d) @ basis_y * dx * dy

    ex = (nxm * np.pi / L) ** 2 / (2.0 * M)
    ey = (nym * np.pi / L) ** 2 / (2.0 * M)
    enm = ex[:, None] + ey[None, :]

    t_idx = np.arange(0, nt, 1)
    times = t_idx * DT * STRIDE

    psi2_analytic = np.empty_like(data, dtype=float)
    for i, t in enumerate(times):
        psi_t = get_analytic_solution_2d(t, basis_x, basis_y, cnm, enm)
        psi2_analytic[i] = np.abs(psi_t) ** 2
    plot_numeric_analytic(data, psi2_analytic, xs, ys, T_IDX)
    max_err = np.max(np.abs(data - psi2_analytic), axis=(1, 2))
    proba = np.sum(data, axis=(1, 2)) * dx * dy
    proba_ana = np.sum(psi2_analytic, axis=(1, 2)) * dx * dy
    plot_error_and_prob(times, max_err, proba, proba_ana)

if __name__ == "__main__":
    main()
