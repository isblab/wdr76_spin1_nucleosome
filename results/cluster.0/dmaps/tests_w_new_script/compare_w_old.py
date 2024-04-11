import sys
import numpy as np
import matplotlib.pyplot as plt

path_old = "test/WDR76.0-SPIN1.0_Distance-matrix.csv"
path_new = "test/distance_matrices/WDR76.0-SPIN1.0_mean_distances.csv"

old_dmat = np.loadtxt(path_old, delimiter=",")
old_dmat = np.round(old_dmat, 2)
old_dmat = np.delete(old_dmat, 0, 0)
old_dmat = np.delete(old_dmat, 0, 1)

new_dmat = np.loadtxt(path_new)
new_dmat = np.round(new_dmat, 2)

matches = np.where(old_dmat == new_dmat, 1, 0)
# plt.imshow(old_dmat, cmap="hot", vmax=30)
# plt.savefig("old.png")
# plt.close()
# plt.imshow(new_dmat, cmap="hot", vmax=30)
# plt.savefig("new.png")
# plt.close()
plt.imshow(matches, vmin=0, vmax=1, cmap="Greys")
plt.colorbar()
plt.savefig("matches.png")
plt.close()

print(old_dmat[300][1])
print(new_dmat[300][1])
