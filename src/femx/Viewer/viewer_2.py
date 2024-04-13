import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors
from matplotlib.cm import ScalarMappable

fig, ax = plt.subplots()
plt.gca().set_aspect('equal')

scalar = 1e-13

x = (np.genfromtxt(r"output/x_coordinates.txt", delimiter=' ') +
     scalar * np.genfromtxt(r"output/displacement_x.txt", delimiter=' '))
y = (np.genfromtxt(r"output/y_coordinates.txt", delimiter=' ') +
     scalar * np.genfromtxt(r"output/displacement_y.txt", delimiter=' '))
z = np.genfromtxt(r"output/von_mises_stress.txt", delimiter=' ')
file = open(r"output/number_of_elements.txt")
n = int(file.read())
print(n)
entries = len(x)
print(entries)
entries_per_element = int(entries/n)
print(entries_per_element)

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
