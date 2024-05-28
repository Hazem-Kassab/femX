from pathlib import Path

import numpy as np
from PyQt6.QtWidgets import QWidget, QApplication, QVBoxLayout, QComboBox
from matplotlib import pyplot as plt, tri
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.cm import ScalarMappable
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import sys

from femx.viewer import Output


class MyWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.fig, self.ax = plt.subplots()
        plt.gca().set_aspect('equal')
        cmap = plt.get_cmap('turbo')
        m = ScalarMappable(cmap=cmap)
        self.cbar = self.fig.colorbar(mappable=m, cmap=cmap, ax=self.ax)
        layout = QVBoxLayout()
        self.setLayout(layout)
        self.list_box = QComboBox()
        self.list_box.addItems([e.name for e in Output])
        self.list_box.setCurrentText("displacement_x")
        layout.addWidget(self.list_box)
        self.list_box.currentIndexChanged.connect(self.display_results)
        self.canvas = FigureCanvas(self.fig)
        self.tool_bar = NavigationToolbar(self.canvas)
        layout.addWidget(self.tool_bar)
        layout.addWidget(self.canvas)
        self.display_results()

    @staticmethod
    def load_coordinates():
        x = np.genfromtxt(
            r"/output/x_coordinates.txt",
            delimiter=' ')
        y = np.genfromtxt(
            r"/output/y_coordinates.txt",
            delimiter=' ')
        return x, y

    @staticmethod
    def load_results(file_name):
        return np.genfromtxt(
            fr"D:\LENOVO\Desktop\ASU MS\4TH TERM\THEORY OF PLASTICITY CES610\Project\output\{file_name}.txt",
            delimiter=' ')

    def display_results(self):
        file_name = str(self.list_box.currentText())
        path = Path(__file__).parent
        file = open( path / "../output/number_of_elements.txt")
        n = int(file.read())

        xc = np.genfromtxt(path / "../output/x_coordinates.txt", delimiter=' ')
        xd = np.genfromtxt(path / "../output/displacement_x.txt", delimiter=' ')
        yc = np.genfromtxt(path / "../output/y_coordinates.txt", delimiter=' ')
        yd = np.genfromtxt(path / "../output/displacement_y.txt", delimiter=' ')

        scale = 0.02*max(max(abs(xc)), max(abs(yc))) / max(max(abs(xd)), max(abs(yd)))
        x = xc + scale * xd
        y = yc + scale * yd
        z = np.genfromtxt(path / f"../output/{file_name}.txt", delimiter=' ')

        z_min = min(z)
        z_max = max(z)

        cmap = plt.get_cmap('turbo')
        m = ScalarMappable(cmap=cmap)
        # cbar = self.fig.colorbar(mappable=m, cmap=cmap, ax=self.ax, values=np.linspace(z_min, z_max, 5))
        # cbar.set_ticks(np.linspace(z_min, z_max, 20))

        file = open(Path(__file__).parent / "../output/number_of_elements.txt")
        n = int(file.read())
        entries = len(x)
        entries_per_element = int(entries / n)

        for i in range(n):
            xe = x[i * entries_per_element:(i + 1) * entries_per_element]
            ye = y[i * entries_per_element:(i + 1) * entries_per_element]
            ze = z[i * entries_per_element:(i + 1) * entries_per_element]
            contour = self.ax.tricontourf(xe, ye, ze, alpha=1, cmap=cmap, origin='upper', vmin=z_min,
                                          vmax=z_max, levels=20, antialiased=False)
        self.cbar.update_normal(contour)
        self.canvas.draw()


app = QApplication(sys.argv)
window = MyWindow()
window.show()
sys.exit(app.exec())