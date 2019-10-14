import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from scipy import interpolate


def equatorialMagneticField(r, phi):
    B = 1.030e6 * r ** (-3.756 - 0.12 * np.cos(phi - 3.562)) + \
        (3.797 - 4.612 * np.cos(phi - 0.825) + 0.606 * np.cos(2 * (phi - 0.473)) +
         0.847 * np.cos(3 * (phi - 0.913))) * np.exp(-1 * r / 150)
    return B


file = open('Test.txt', 'w')

x = np.arange(20, 101, 1)
y = np.arange(0, 2 * np.pi, 0.1)

xi, yi = np.meshgrid(x, y)
z = np.array(equatorialMagneticField(xi, yi))

xi = x * np.cos(yi)
yi = x * np.sin(yi)

heatmap = plt.contourf(xi, yi, z, 100, cmap=plt.cm.get_cmap('gist_rainbow'), alpha=0.4)
lines = plt.contour(xi, yi, z, 5, colors='k')
plt.clabel(lines, fontsize=18, inline=1, colors='k')
plt.colorbar(heatmap)
plt.show()
