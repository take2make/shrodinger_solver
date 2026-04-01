import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import glob

dat = "dat/solution_temporelle_shrodinger.dat"

file = glob.glob(dat)[0]
data = np.loadtxt(file)
psi_2 = data

fig, ax = plt.subplots()

x = np.linspace(0,10,201)
line, = ax.plot(x, psi_2[0])
ax.set(xlim=[0, 10], ylim=[0, 1], xlabel='Position [m]', ylabel=r'$|\Psi|^{2}$')
ax.bar(x=4.5,height=1,width=0.02,color="red")
ax.bar(x=5.5,height=1,width=0.02,color="red")

def update(frame):
    # update the line plot
    line.set_ydata(psi_2[frame])
    return (line,)


ani = FuncAnimation(fig=fig, func=update, frames=len(psi_2), interval=10,blit=True)
plt.show()

ani.save("schrodinger1.mp4")