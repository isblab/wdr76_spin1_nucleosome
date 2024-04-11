import os
import numpy as np
import matplotlib.pyplot as plt

d_mat = np.loadtxt("WDR76.0-SPIN1.0_Distance-matrix.csv", delimiter=",")
print(d_mat.shape)

plt.imshow(d_mat, cmap="hot", vmax=20.0)
plt.colorbar()
plt.xlabel("SPIN1.0")
plt.ylabel("WDR76.0")

plt.show()
