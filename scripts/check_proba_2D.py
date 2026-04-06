import matplotlib.pyplot as plt
import numpy as np

dat = "dat/solution_shrodinger_leapfrog_2D.dat"

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
    return np.stack(frames, axis=0)


def check_proba(dat: str) -> None:
	data = load_data_2d(dat)

	proba = np.sum(data, axis=(1, 2))
	proba = proba / proba[0]

	fig, ax = plt.subplots(figsize=(7, 4))
	ax.plot(proba, label="Leapfrog 2D")
	ax.set_xlabel("Time step")
	ax.set_ylabel("Total probability")
	ax.set_title("Leapfrog 2D Total Probability")
	ax.legend()
	ax.grid()
	plt.tight_layout()
	plt.show()


if __name__ == "__main__":
	check_proba(dat)
