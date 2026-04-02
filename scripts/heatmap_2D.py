import glob

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

DATA_PATH = "dat/solution_shrodinger_leapfrog_2D.dat"
FPS = 20


def load_data(path: str) -> np.ndarray:
    if path.endswith(".npy"):
        data = np.load(path)
        if data.ndim != 3:
            raise ValueError(f"Expected (t, y, x) array, got {data.shape}")
        return data

    if path.endswith(".dat"):
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
        for block in blocks:
            xs = np.unique(block[:, 0])
            ys = np.unique(block[:, 1])
            nx = xs.size
            ny = ys.size
            if nx * ny != block.shape[0]:
                raise ValueError("Inconsistent grid in data block.")
            frame = block[:, 2].reshape(ny, nx)
            frames.append(frame)
        return np.stack(frames, axis=0)

    files = sorted(glob.glob(path))
    if not files:
        raise FileNotFoundError(f"No files matched pattern: {path}")
    frames = [np.loadtxt(f) for f in files]
    return np.stack(frames, axis=0)


def main() -> None:
    data = load_data(DATA_PATH)

    fig, ax = plt.subplots(figsize=(6, 5))
    vmin = np.percentile(data, 2)
    vmax = np.percentile(data, 98)
    im = ax.imshow(
        data[0],
        origin="lower",
        aspect="auto",
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

    interval_ms = int(1000 / FPS)
    global anim
    anim = FuncAnimation(fig, update, frames=data.shape[0], interval=interval_ms, blit=False)
    plt.show()


if __name__ == "__main__":
    main()
