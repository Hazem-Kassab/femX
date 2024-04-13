import csv

import numpy as np
from matplotlib import pyplot as plt

from Elements.element import Element
from strucutre import Structure


class Visualizer:

    def __init__(self, structure: Structure):
        self.structure = structure

    def visualize(self):
        fig, ax = plt.subplots()
        plt.gca().set_aspect('equal')
        x = np.genfromtxt(r"../output/x_coordinates.txt", delimiter=' ')
        y = np.genfromtxt(r"../output/y_coordinates.txt", delimiter=' ')
        z = np.genfromtxt(r"../output/displacement_x.txt", delimiter=' ')

        contour = ax.tricontourf(x, y, z, alpha=1, cmap='turbo', origin='upper', levels=20)
        fig.colorbar(contour)
        for element in self.structure.elements:
            Visualizer.plot_element(element)
            Visualizer.plot_deformed_element(element)

        plt.show()

    @staticmethod
    def element_coordinates_array(element: Element):
        return np.append(element.node_coordinates(), np.resize(element.node_coordinates()[:, 0],
                                                                             (2, 1)), axis=1)

    @staticmethod
    def plot_element(element: Element):
        coordinates_array = Visualizer.element_coordinates_array(element)
        x = coordinates_array[0, :]
        y = coordinates_array[1, :]
        plt.plot(x, y, color="black")

    @staticmethod
    def plot_deformed_element(element: Element):
        displacements = []
        for node in element.nodes:
            displacements.append([dof.displacement for dof in node.degrees_of_freedom])

        displacements = np.append(np.array(displacements).T, np.resize(np.array([displacements])[:, 0],
                                                                       (2, 1)), axis=1)

        deformed_coordinates_array = displacements + Visualizer.element_coordinates_array(element)
        x = deformed_coordinates_array[0, :]
        y = deformed_coordinates_array[1, :]
        plt.plot(x, y, linestyle="dashed", color="grey")
