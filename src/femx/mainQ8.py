import matplotlib.pyplot as plt

from Elements.PlaneElements.Q8_element import Q8Element
from Materials.PressureIndependentMaterials.steel import Steel
from Solvers.nonlinear_solver import NonlinearSolver
from node import Node
from Solvers.solver import Solver
from strucutre import Structure

l = 6000
d = 500
node_spacing = 250

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
    # print(f"*********** i = {i}")
    j = 0
    while j < int(len(nodes[i]))-2:
        # print(f"********* j ={j}")
        elements.append(
            Q8Element(nodes[i][j], nodes[i + 1][j], nodes[i + 2][j], nodes[i + 2][j + 1],
                      nodes[i + 2][j + 2], nodes[i + 1][j + 2], nodes[i][j + 2],
                      nodes[i][j + 1], 20, steel))
        j += 2
    i += 2

# for node in nodes[0]:
#     node.x_dof.restrained = True
#     node.y_dof.restrained = True

a = int(l/(node_spacing*2))
b = int(d/(node_spacing*2))
print(a)
print(b)

nodes[0][b].y_dof.restrained = True
nodes[-1][b].y_dof.restrained = True

force = -2e5
node = nodes[a][b]
node.y_dof.force = force
# nodes[3][2].y_dof.force = -30e4

structure = Structure(elements)
s = Solver(structure)
s.solve()
node.y_dof.force = -3e5
s = NonlinearSolver(structure)
disp, force = s.solve_incrementally(node.y_dof, node.y_dof, force, steps=15)
plt.plot(disp, force)
plt.xlabel("displacement (mm)")
plt.ylabel("Force (N)")
plt.show()
s.output()
