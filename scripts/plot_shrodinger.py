import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

dat_cn = "dat/solution_temporelle_shrodinger.dat"
dat_lp = "dat/solution_shrodinger_leapfrog.dat"
dat_kg = "dat/klein_gordon_output.dat"

def animate_2_evolutions(dat1, dat2):
    data1 = np.loadtxt(dat1)
    data2 = np.loadtxt(dat2)

    fig, axs = plt.subplots(1, 2, figsize=(10, 4))
    x = np.linspace(0,10,data1.shape[1])
    line1, = axs[0].plot(x, data1[0])
    line2, = axs[1].plot(x, data2[0])
    axs[0].set(xlim=[0, 10], ylim=[0, 1], xlabel='Position [m]', ylabel=r'$|\Psi|^{2}$', title="Crank-Nicolson")
    axs[0].axvspan(0, 5.5, facecolor="red", alpha=0.12, hatch="///", edgecolor="red")
    #axs[0].bar(x=4.5,height=1,width=0.02,color="red")
    axs[0].bar(x=5.5,height=1,width=0.02,color="red")
    axs[1].set(xlim=[0, 10], ylim=[0, 1], xlabel='Position [m]', ylabel=r'$|\Psi|^{2}$', title="Leapfrog")
    axs[1].axvspan(0, 5.5, facecolor="red", alpha=0.12, hatch="///", edgecolor="red")
    #axs[1].bar(x=4.5,height=1,width=0.02,color="red")
    axs[1].bar(x=5.5,height=1,width=0.02,color="red")

    def update(frame):
        line1.set_ydata(data1[frame])
        line2.set_ydata(data2[frame])
        return (line1, line2)

    ani = FuncAnimation(fig=fig, func=update, frames=len(data1), interval=10, blit=True)
    plt.show()

def animate_evolution_for(dat):
    data = np.loadtxt(dat)

    fig = plt.figure(figsize=(10, 6))
    x = np.linspace(0,10,data.shape[1])
    plt.plot(x, data[0], label="t=0", color="gray", alpha=0.5, linestyle="--")
    line, = plt.plot(x, data[0])
    plt.xlim(0, 10)
    plt.ylim(0, 1)
    plt.xlabel('Position [m]')
    plt.ylabel(r'$|\Psi|^{2}$')
    plt.title('Evolution of $|\Psi|^2$ over time')
    plt.axvspan(0, 5.5, facecolor="red", alpha=0.05, hatch="///", edgecolor="red")
    #plt.bar(x=4.5,height=1,width=0.02,color="red")
    plt.bar(x=5.5,height=1,width=0.02,color="red")

    def update(frame):
        line.set_ydata(data[frame])
        return (line,)

    ani = FuncAnimation(fig=fig, func=update, frames=len(data), interval=10, blit=True)
    plt.show()

def plot_evolution_at(dat, time_idx):
    data = np.loadtxt(dat)
    x = np.linspace(0,10,data.shape[1])
    plt.plot(x, data[0], label="t=0", color="gray", alpha=0.5, linestyle="--")
    plt.plot(x, data[time_idx], label=f"t={time_idx}")
    plt.xlim(0, 10)
    plt.ylim(0, 1)
    plt.xlabel('Position [m]')
    plt.ylabel(r'$|\Psi|^{2}$')
    plt.title(f'Evolution of $|\Psi|^2$')
    plt.axvspan(0, 5.5, facecolor="red", alpha=0.05, hatch="///", edgecolor="red")
    #plt.bar(x=4.5,height=1,width=0.02,color="red")
    plt.bar(x=5.5,height=1,width=0.02,color="red")
    plt.legend()
    plt.savefig(f"evolution_at_time_{time_idx}.png", dpi=300)
    plt.show()

if __name__ == "__main__":
    #animate_2_evolutions(dat_cn, dat_lp)
    #animate_evolution_for(dat_lp)
    #plot_evolution_at(dat_lp, time_idx=800)
    animate_evolution_for(dat_kg)