import argparse
import glob
import numpy as np
import matplotlib.pyplot as plt

DEFAULT_FILES = [
    "dat/solution_temporelle_shrodinger.dat",
    "dat/solution_shrodinger_leapfrog.dat",
]
    
def check_proba(dat1: str, dat2: str) -> None:
    data1 = np.loadtxt(dat1)
    data2 = np.loadtxt(dat2)

    proba1 = np.sum(data1, axis=1)
    proba2 = np.sum(data2, axis=1)

    proba1 = proba1/proba1[0]
    proba2 = proba2/proba2[0]

    fig, axs = plt.subplots(1, 2, figsize=(12, 5))
    axs[0].plot(proba1, label="Crank-Nicolson")
    axs[0].set_xlabel("Time step")
    axs[0].set_ylabel("Total probability")
    axs[0].set_title("Crank-Nicolson Total Probability")
    #axs[0].set_ylim(0.9, 1.1)
    axs[0].legend()
    axs[0].grid()
    axs[1].plot(proba2, label="Leapfrog")
    axs[1].set_xlabel("Time step")
    axs[1].set_ylabel("Total probability")
    axs[1].set_title("Leapfrog Total Probability")
    #axs[1].set_ylim(0.9, 1.1)
    axs[1].legend()
    axs[1].grid()
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    check_proba(DEFAULT_FILES[0], DEFAULT_FILES[1])
