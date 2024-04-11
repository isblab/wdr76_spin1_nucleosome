import os
import sys
import glob
import numpy as np
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt


def get_binarized_distance_matrix(matrix: np.ndarray, threshold: float) -> np.ndarray:
    return np.where(matrix <= threshold, 0, 1)


THRESHOLD = 25.0
all_files = glob.glob("./*Distance-matrix.csv")
if not "binarized_dmaps" in os.listdir():
    os.mkdir("binarized_dmaps")

for FNAME in tqdm(all_files):
    distance_matrix = np.loadtxt(FNAME, delimiter=",", dtype=np.float64)
    distance_matrix = np.delete(np.delete(distance_matrix, 0, axis=0), 0, axis=1)

    contact_matrix = get_binarized_distance_matrix(distance_matrix, THRESHOLD)
    df = pd.DataFrame(contact_matrix, index=None, columns=None)

    FNAME = FNAME.split("/")[-1]
    plt.imshow(contact_matrix, cmap="hot", vmin=0, vmax=1)
    plt.xlabel(FNAME.split("_", maxsplit=1)[0].split("-")[1])
    plt.ylabel(FNAME.split("_", maxsplit=1)[0].split("-")[0])
    plt.xticks(np.arange(0, contact_matrix.shape[1], 40))
    plt.yticks(np.arange(0, contact_matrix.shape[0], 40))
    plt.savefig(f"binarized_dmaps/binarized_{FNAME[:-10]}map.png")
