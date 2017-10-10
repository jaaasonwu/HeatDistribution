import numpy as np
import matplotlib.pyplot as plt

a = np.loadtxt("output.dat", skiprows=1)
plt.imshow(a, cmap="hot")
plt.colorbar()
plt.show()