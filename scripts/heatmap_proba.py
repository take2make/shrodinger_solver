import numpy as np
import matplotlib.pyplot as plt
import glob

dat_cn = "dat/solution_temporelle_shrodinger.dat"
dat_lp = "dat/solution_shrodinger_leapfrog.dat"

def plot_heatmap(dat: str, title: str) -> None:
    data = np.loadtxt(dat)
    psi_norm = data

    plt.figure(figsize=(10, 6))
    plt.imshow(psi_norm, aspect='auto',extent=[0, 10, 0, 600], cmap='viridis', origin='lower')
    plt.colorbar(label="psi^2(x, t)")
    plt.xlabel("$x$")
    plt.ylabel("$t$")
    plt.title("Shrodinger evolution of psi^2(x, t)")
    plt.show()
    
def plot_2_heatmaps(dat1: str, dat2: str, title1: str, title2: str) -> None:
    data1 = np.loadtxt(dat1)
    data2 = np.loadtxt(dat2)
    
    data2_resized = np.resize(data2, data1.shape)

    fig, axs = plt.subplots(1, 2, figsize=(10, 4))
    im1 = axs[0].imshow(data1, aspect='auto',extent=[0, 10, 0, 10], cmap='viridis', origin='lower')
    axs[0].set_title(title1)
    axs[0].set_xlabel("$x$")
    axs[0].set_ylabel("$t$")
    fig.colorbar(im1, ax=axs[0], label="psi^2(x, t)")

    im2 = axs[1].imshow(data2_resized, aspect='auto',extent=[0, 10, 0, 10], cmap='viridis', origin='lower')
    axs[1].set_title(title2)
    axs[1].set_xlabel("$x$")
    axs[1].set_ylabel("$t$")
    fig.colorbar(im2, ax=axs[1], label="psi^2(x, t)")

    plt.tight_layout()
    plt.show()
    
def compare_heatmaps(dat1: str, dat2: str) -> None:
    data1 = np.loadtxt(dat1)
    data2 = np.loadtxt(dat2)

    # Приводим к одному размеру (например, к размеру первого массива)
    data2_resized = np.resize(data2, data1.shape)

    # Визуализация
    fig, axs = plt.subplots(1, 3, figsize=(12, 4))
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
    #compare_heatmaps(dat_cn, dat_lp)
