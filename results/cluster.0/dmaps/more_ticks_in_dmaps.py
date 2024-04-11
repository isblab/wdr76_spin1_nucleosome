import os
import glob
import numpy as np
import matplotlib.pyplot as plt

dmat_fnames = glob.glob(os.path.join(os.getcwd(), "*csv"))
for dmat_fname in dmat_fnames:
    dmat = np.loadtxt(dmat_fname, delimiter=",")

    dmat_fname = dmat_fname.split("/")[-1]

    plt.imshow(dmat, cmap="hot", vmax=30)
    plt.xticks(np.arange(0, dmat.shape[1], 40))
    plt.yticks(np.arange(0, dmat.shape[0], 40))
    plt.xlabel(dmat_fname.split("_")[0].split("-")[1])
    plt.ylabel(dmat_fname.split("_")[0].split("-")[0])
    plt.colorbar()
    plt.savefig(f"{dmat_fname.split('_')[0]}.png")
    plt.close()
