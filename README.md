# femx

A package for finite element analysis of structures

## Installation

```bash
$ pip install femx
```

## Usage

```bash
# cantilever.py
# analysis of a cantilever beam using 4-node plane stress elements

from femx.materials.pressure_independent_materials.steel import Steel
from femx.elements.plane_elements.Q4_element import Q4Element
from femx.node import Node
from femx.solvers.solver import Solver
from femx.strucutre import Structure

# define the length and depth of the beam
l = 4000
d = 1000
element_size = 50

# initialize empty list to hold nodes
nodes = []

# create node objects at equal spacing
for i in range(0, l+element_size, element_size):
    row = []
    for j in range(0, d+element_size, element_size):
        row.append(Node(i, j))
    nodes.append(row)

# create material object
steel = Steel(200e3, 0.3, 250, 360)

# initialize empty list to hold elements
elements = []

# create  element objects connecting nodes
i = 0
while i <= len(nodes)-2:
    j = 0
    while j <= len(nodes[i])-2:
        elements.append(Q4Element(nodes[i][j], nodes[i+1][j], nodes[i+1][j+1], nodes[i][j+1], 20, steel))
        j += 1
    i += 1


# assign boundary conditions
for node in nodes[0]:
    node.x_dof.restrained = True
    node.y_dof.restrained = True

# assign forces
nodes[-1][1].y_dof.force = -30e5

# assemble a structure
structure = Structure(elements)

# solve for displacements and reactions
s = Solver(structure)
s.solve()

# view results graphically using viewers
s.output()
```

Upon running the above code, output files will be generated in the output folder. To graphically view
the model, you should run the "viewer.py" python file:

![alt text]()

## Non-linear analysis of simply supported beam using 8-node plane stress elements

![alt text]()

load-displacement curve

![alt text]()

## Contributing

Interested in contributing? Check out the contributing guidelines. Please note that this project is released with a Code of Conduct. By contributing to this project, you agree to abide by its terms.

## License

`femx` was created by Hazem Kassab. It is licensed under the terms of the MIT license.
