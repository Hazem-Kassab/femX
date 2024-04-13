class DOF:
    counter = 0

    def __init__(self):
        self.restrained = False
        self.displacement = 0
        self.force = 0
        DOF.counter += 1
        self.id = DOF.counter

    @property
    def displacement(self):
        return self._displacement

    @displacement.setter
    def displacement(self, value):
        self._displacement = value

    @property
    def force(self):
        return self._force

    @force.setter
    def force(self, value):
        self._force = value

    @property
    def restrained(self):
        return self._restrained

    @restrained.setter
    def restrained(self, value: bool):
        self._restrained = value

    def __repr__(self):
        return f"degree of freedom {self.id}"
