from Materials.PressureIndependentMaterials.steel import Steel
from Elements.PlaneElements.Q4_element import Q4Element
from node import Node
from Solvers.solver import Solver
from strucutre import Structure

l = 6000
d = 1000
element_size = 125

nodes = []

for i in range(0, l+element_size, element_size):
    row = []
    for j in range(0, d+element_size, element_size):
        row.append(Node(i, j))
    nodes.append(row)

steel = Steel(200e3, 250, 0.3)

elements = []

i = 0
while i <= len(nodes)-2:
    j = 0
    while j <= len(nodes[i])-2:
        elements.append(Q4Element(nodes[i][j], nodes[i+1][j], nodes[i+1][j+1], nodes[i][j+1], 20, steel))
        j += 1
    i += 1

# for i in nodes:
#     print(i)

for node in nodes[0]:
    node.x_dof.restrained = True
    node.y_dof.restrained = True
# nodes[0][0].y_dof.restrained = True

# print(nodes[-1][-1])

nodes[-1][1].y_dof.force = -30e5

# print(elements[0].get_coordinates(0, 0))

structure = Structure(elements)
s = Solver(structure)
s.solve()
# s.average_nodal_values()
s.output()
