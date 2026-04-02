import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

dat_cn = "dat/solution_temporelle_shrodinger.dat"
dat_lp = "dat/solution_shrodinger_leapfrog.dat"

def animate_2_evolutions(dat1: str, dat2: str):
    data1 = np.loadtxt(dat1)
    data2 = np.loadtxt(dat2)

    fig, axs = plt.subplots(1, 2, figsize=(10, 4))
    x = np.linspace(0,10,data1.shape[1])
    line1, = axs[0].plot(x, data1[0])
    line2, = axs[1].plot(x, data2[0])
    axs[0].set(xlim=[0, 10], ylim=[0, 1], xlabel='Position [m]', ylabel=r'$|\Psi|^{2}$', title="Crank-Nicolson")
    axs[0].bar(x=4.5,height=1,width=0.02,color="red")
    axs[0].bar(x=5.5,height=1,width=0.02,color="red")
    axs[1].set(xlim=[0, 10], ylim=[0, 1], xlabel='Position [m]', ylabel=r'$|\Psi|^{2}$', title="Leapfrog")
    axs[1].bar(x=4.5,height=1,width=0.02,color="red")
    axs[1].bar(x=5.5,height=1,width=0.02,color="red")

    def update(frame):
        line1.set_ydata(data1[frame])
        line2.set_ydata(data2[frame])
        return (line1, line2)

    ani = FuncAnimation(fig=fig, func=update, frames=len(data1), interval=10, blit=True)
    plt.show()

def animate_evolution_for(dat: str):
    data = np.loadtxt(dat)

    fig = plt.figure(figsize=(10, 6))
    x = np.linspace(0,10,data.shape[1])
    line, = plt.plot(x, data[0])
    plt.xlim(0, 10)
    plt.ylim(0, 1)
    plt.xlabel('Position [m]')
    plt.ylabel(r'$|\Psi|^{2}$')
    plt.title('Evolution of $|\Psi|^2$ over time')
    plt.bar(x=4.5,height=1,width=0.02,color="red")
    plt.bar(x=5.5,height=1,width=0.02,color="red")

    def update(frame):
        line.set_ydata(data[frame])
        return (line,)

    ani = FuncAnimation(fig=fig, func=update, frames=len(data), interval=10, blit=True)
    plt.show()

if __name__ == "__main__":
    animate_2_evolutions(dat_cn, dat_lp)
    animate_evolution_for(dat_lp)