import numpy as np
import matplotlib.pyplot as plt
import sys

a = np.loadtxt(sys.argv[1], skiprows=1)
plt.imshow(a, cmap="hot")
plt.colorbar()
plt.show()