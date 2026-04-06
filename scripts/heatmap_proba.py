import numpy as np
import matplotlib.pyplot as plt

dat_cn = "dat/solution_temporelle_shrodinger.dat"
dat_lp = "dat/solution_shrodinger_leapfrog.dat"

def plot_2_heatmaps(dat1, dat2, title1, title2):
    data1 = np.loadtxt(dat1)
    data2 = np.loadtxt(dat2)
    
    data2_resized = np.resize(data2, data1.shape)
    max1 = np.max(data1)
    max2 = np.max(data2_resized)
    data1_norm = data1 / max1 if max1 != 0 else data1
    data2_norm = data2_resized / max2 if max2 != 0 else data2_resized

    fig, axs = plt.subplots(1, 2, figsize=(10, 4))
    im1 = axs[0].imshow(data1_norm, aspect='auto',extent=[0, 10, 0, 10], cmap='inferno', origin='lower', vmin=0, vmax=1)
    axs[0].set_title(title1)
    axs[0].set_xlabel("$x$")
    axs[0].set_ylabel("$t$")
    fig.colorbar(im1, ax=axs[0])

    im2 = axs[1].imshow(data2_norm, aspect='auto',extent=[0, 10, 0, 10], cmap='inferno', origin='lower', vmin=0, vmax=1)
    axs[1].set_title(title2)
    axs[1].set_xlabel("$x$")
    axs[1].set_ylabel("$t$")
    fig.colorbar(im2, ax=axs[1])

    plt.tight_layout()
    plt.savefig("compare_heatmap.png", dpi=300)
    plt.show()
    
def compare_heatmaps(dat1, dat2):
    data1 = np.loadtxt(dat1)
    data2 = np.loadtxt(dat2)

    data2_resized = np.resize(data2, data1.shape)

    _, axs = plt.subplots(1, 3, figsize=(12, 4))
    axs[0].imshow(data1, aspect='auto', cmap='viridis', origin='lower')
    axs[0].set_title('Heatmap 1')
    axs[1].imshow(data2_resized, aspect='auto', cmap='viridis', origin='lower')
    axs[1].set_title('Heatmap 2 (resized)')
    diff = np.abs(data1 - data2_resized)
    axs[2].imshow(diff, aspect='auto', cmap='viridis', origin='lower')
    axs[2].set_title(f'Difference\nMSE={np.mean(diff**2):.4e}')
    plt.show()

if __name__ == "__main__":
    plot_2_heatmaps(dat_cn, dat_lp, "Crank-Nicolson |psi|^2", "Leapfrog |psi|^2")
