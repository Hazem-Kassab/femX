import os

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.cm import ScalarMappable

from femx.viewer import Output

from pathlib import Path

fig, ax = plt.subplots()
plt.gca().set_aspect('equal')


def visualize(output: Output, scale):
    path = Path(__file__).parent
    x = (np.genfromtxt(path / "../output/x_coordinates.txt", delimiter=' ') +
         scale * np.genfromtxt(path / "../output/displacement_x.txt", delimiter=' '))
    y = (np.genfromtxt(path / "../output/y_coordinates.txt", delimiter=' ') +
         scale * np.genfromtxt(path / "../output/displacement_y.txt", delimiter=' '))
    z = np.genfromtxt(path / f"../output/{output.name}.txt", delimiter=' ')

    file = open(Path(__file__).parent / "../output/number_of_elements.txt")
    n = int(file.read())
    entries = len(x)
    entries_per_element = int(entries/n)

    z_min = min(z)
    z_max = max(z)

    cmap = plt.get_cmap('turbo')
    m = ScalarMappable(cmap=cmap)
    cbar = fig.colorbar(mappable=m, cmap=cmap, ax=ax, values=np.linspace(z_min, z_max, 5))
    cbar.set_ticks(np.linspace(z_min, z_max, 20))

    for i in range(n):
        xe = x[i*entries_per_element:(i+1)*entries_per_element]
        ye = y[i*entries_per_element:(i+1)*entries_per_element]
        ze = z[i*entries_per_element:(i+1)*entries_per_element]
        contour = ax.tricontourf(xe, ye, ze, alpha=1, cmap=cmap, origin='upper', vmin=z_min, vmax=z_max, levels=20)

    plt.show()
