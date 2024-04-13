from Elements.element import Element


class Structure:

    def __init__(self, elements: list[Element]):
        self.elements = elements

    @property
    def nodes(self):
        nodes = set()
        for element in self.elements:
            for node in element.nodes:
                nodes.add(node)
        return nodes

    def degrees_of_freedom(self):
        dofs = set()
        for node in self.nodes:
            for dof in node.degrees_of_freedom:
                dofs.add(dof)

        return sorted(dofs, key=lambda dof: dof.id)

    @property
    def free_degrees_of_freedom(self):
        return sorted(set(dof for dof in self.degrees_of_freedom() if not dof.restrained), key=lambda dof: dof.id)

    @property
    def restrained_degrees_of_freedom(self):
        return sorted(set(dof for dof in self.degrees_of_freedom() if dof.restrained), key=lambda dof: dof.id)


