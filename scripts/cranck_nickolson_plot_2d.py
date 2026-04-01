import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation


DATA_PATH = "dat/cranck_nickolson_2d.dat"

def parse(entry: int):
    entry = entry[1:-1].split(",")
    entry = [float(s) for s in entry]
    return entry


def update(frame):
    # update the line plot
    im.set_data(data[frame,:].reshape(11,11))
    ax.set_title(f"t={frame}")
    return (im,)

data = np.loadtxt("dat/cranck_nickolson_2d.dat")

fig, ax = plt.subplots()

im = ax.imshow(data[0,:].reshape(11,11))
ani = FuncAnimation(fig=fig, func=update, frames=500, interval=10,blit=True)
ani.save("cranck_nickolson_2d.mp4")
plt.show()



""" print(data.shape)
fig, ax = plt.subplots()
ax.imshow(data[100,:].reshape(11,11))
plt.savefig("cn_2d.png") """