import numpy as np
from matplotlib import pyplot as plt, tri
import matplotlib.colors
from matplotlib.cm import ScalarMappable
from matplotlib.tri import TriAnalyzer

fig, ax = plt.subplots()
plt.gca().set_aspect('equal')

x = np.genfromtxt(r"../output/x_displaced_coordinates.txt", delimiter=' ')
y = np.genfromtxt(r"../output/y_displaced_coordinates.txt", delimiter=' ')

triang = tri.Triangulation(x, y)
mask = TriAnalyzer(triang).get_flat_tri_mask(min_circle_ratio=0.1, rescale=True)
triang.set_mask(mask)

z = np.genfromtxt(r"../output/stress_xy.txt", delimiter=' ')
contour = ax.tricontourf(triang, z, alpha=1, origin='upper', cmap="turbo", levels=100)
cbar = fig.colorbar(contour)
plt.show()
