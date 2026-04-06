import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

dat = "dat/solution_shrodinger_leapfrog_2D.dat"
L = 10.0

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

def animate_heatmap(data, xs, ys):
    fig, ax = plt.subplots(figsize=(6, 5))
    vmin = np.percentile(data, 2)
    vmax = np.percentile(data, 98)
    im = ax.imshow(
        data[0],
        origin="lower",
        extent=(xs[0], xs[-1], ys[0], ys[-1]),
        aspect="equal",
        cmap="inferno",
        vmin=vmin,
        vmax=vmax,
    )
    fig.colorbar(im, ax=ax, label="|psi|^2")
    ax.set_xlabel("x")
    ax.set_ylabel("y")

    def update(frame: int):
        im.set_data(data[frame])
        ax.set_title(f"t={frame}")
        return (im,)
    
    global anim
    anim = FuncAnimation(fig, update, frames=data.shape[0], interval=100, blit=False)
    anim.save("heatmap_animation.gif", writer="pillow", fps=10)
    plt.show()

def heatmap_for_step(data, xs, ys, step):
    fig, ax = plt.subplots(figsize=(6, 5))
    vmin = np.percentile(data, 2)
    vmax = np.percentile(data, 98)
    im = ax.imshow(
        data[step],
        origin="lower",
        extent=(xs[0], xs[-1], ys[0], ys[-1]),
        aspect="equal",
        cmap="inferno",
        vmin=vmin,
        vmax=vmax,
    )
    fig.colorbar(im, ax=ax, label="|psi|^2")
    plt.bar(x=4,height=4,width=0.1,color="white")
    plt.bar(x=4,height=1,width=0.1,color="white", bottom=4.5)
    plt.bar(x=4,height=4,width=0.1,color="white", bottom=6)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(f"t={step}")
    plt.savefig(f"heatmap_step_{step}.png", dpi=300)
    plt.show()


if __name__ == "__main__":
    data, xs, ys = load_data_2d(dat)
    animate_heatmap(data, xs, ys)
    heatmap_for_step(data, xs, ys, step=10)
