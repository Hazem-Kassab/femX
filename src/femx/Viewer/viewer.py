import numpy as np
from PyQt6.QtWidgets import QWidget, QApplication, QVBoxLayout, QComboBox
from matplotlib import pyplot as plt
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import sys


class MyWindow(QWidget):
    def __init__(self):
        super().__init__()
        layout = QVBoxLayout()
        self.setLayout(layout)
        self.list_box = QComboBox()
        self.list_box.addItems(["displacement_x", "displacement_y", "stress_x"])
        self.list_box.setCurrentText("displacement_x")
        # x, y = MyWindow.load_coordinates()
        # z = MyWindow.load_results(self.list_box.currentText())
        layout.addWidget(self.list_box)
        self.list_box.currentIndexChanged.connect(self.display_results())
        self.canvas = FigureCanvas(fig)  # create canvas
        layout.addWidget(self.canvas)  # add canvas to layout

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
        file_name = self.list_box.currentText
        # x = np.genfromtxt(
        #     r"D:\LENOVO\Desktop\ASU MS\4TH TERM\THEORY OF PLASTICITY CES610\Project\output\x_coordinates.txt",
        #     delimiter=' ')
        # y = np.genfromtxt(
        #     r"D:\LENOVO\Desktop\ASU MS\4TH TERM\THEORY OF PLASTICITY CES610\Project\output\y_coordinates.txt",
        #     delimiter=' ')
        z = np.genfromtxt(
            fr"/output/displacement_x.txt",
            delimiter=' ')
        ax.remove()
        c = ax.tricontourf(x, y, z, alpha=1, cmap='turbo', origin='upper', levels=20)
        # fig.colorbar(c)
        # self.canvas.draw()

    def display(self):
        print(self.list_box.currentText())


fig, ax = plt.subplots()
plt.gca().set_aspect('equal')

x = np.genfromtxt(r"/output/x_coordinates.txt", delimiter=' ')
y = np.genfromtxt(r"/output/y_coordinates.txt", delimiter=' ')
z = np.genfromtxt(r"/output/von_mises_stress.txt", delimiter=' ')

contour = ax.tricontourf(x, y, z, alpha=1, cmap='turbo', origin='upper', levels=20)
fig.colorbar(contour)

app = QApplication(sys.argv)
window = MyWindow()
window.show()
sys.exit(app.exec())