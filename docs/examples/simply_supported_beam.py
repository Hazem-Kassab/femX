# simply_supported_beam.py

import matplotlib.pyplot as plt

from femx.elements.plane_elements.Q8_element import Q8Element
from femx.materials.pressure_independent_materials.steel import Steel
from femx.node import Node
from femx.solvers.nonlinear_solver import NonlinearSolver, Solver
from femx.strucutre import Structure

l = 4000
d = 1000
node_spacing = 100

nodes = []

for i in range(0, l + node_spacing, node_spacing):
    row = []
    for j in range(0, d + node_spacing, node_spacing):
        row.append(Node(i, j))
    nodes.append(row)

steel = Steel(200e3, 0.3, 250, 360)

elements = []

i = 0
while i < int(len(nodes))-2:
    j = 0
    while j < int(len(nodes[i]))-2:
        elements.append(
            Q8Element(nodes[i][j], nodes[i + 1][j], nodes[i + 2][j], nodes[i + 2][j + 1],
                      nodes[i + 2][j + 2], nodes[i + 1][j + 2], nodes[i][j + 2],
                      nodes[i][j + 1], 20, steel))
        j += 2
    i += 2

a = int(l/(node_spacing*2))
b = int(d/(node_spacing*2))

nodes[0][b].y_dof.restrained = True
nodes[-1][b].y_dof.restrained = True

force = -1.5e6
node = nodes[a][b]
node.y_dof.force = force


structure = Structure(elements)
# s = Solver(structure)
# s.solve()
# s.output()

s = NonlinearSolver(structure)
loads, disp = s.solve_incrementally(node.y_dof, node.y_dof, 0, 50)
plt.plot(loads, disp)
plt.show()
s.output()
