from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt, tri
from matplotlib.tri import TriAnalyzer


def visualize(output, scale):
    fig, ax = plt.subplots()
    plt.gca().set_aspect('equal')

    path = Path(__file__).parent

    x = (np.genfromtxt(path / "../output/x_coordinates.txt", delimiter=' ') +
         scale * np.genfromtxt(path / "../output/displacement_x.txt", delimiter=' '))
    y = (np.genfromtxt(path / "../output/y_coordinates.txt", delimiter=' ') +
         scale * np.genfromtxt(path / "../output/displacement_y.txt", delimiter=' '))

    triang = tri.Triangulation(x, y)
    mask = TriAnalyzer(triang).get_flat_tri_mask(min_circle_ratio=0.1, rescale=True)
    triang.set_mask(mask)

    z = np.genfromtxt(path / f"../output/{output.name}.txt", delimiter=' ')
    contour = ax.tricontourf(triang, z, alpha=1, origin='upper', cmap="turbo", levels=100)
    cbar = fig.colorbar(contour)
    plt.show()
