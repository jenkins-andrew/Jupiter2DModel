import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from matplotlib import ticker, cm

def equatorialMagneticField(r, phi):
    B = 1.030e6 * r ** (-3.756 - 0.12 * math.cos(phi - 3.562)) + \
(3.797 - 4.612 * math.cos(phi - 0.825) + 0.606 * math.cos(2 * (phi - 0.473)) +
         0.847 * math.cos(3 * (phi - 0.913))) * math.exp(-1 * r / 150)
    return B


x = []
y = []
z = []

file = open('Test.txt', 'w')

for r in np.arange(1, 51, 0.5):
    for phi in np.arange(0, 2 * math.pi, 0.05):
        x.append(r * math.cos(phi))
        y.append(r * math.sin(phi))
        z.append(equatorialMagneticField(r, phi))


x = np.array(x)
y = np.array(y)
z = np.array(z)

xmin, xmax = -50, 50
ymin, ymax = -50, 50
xi = np.arange(xmin, xmax, 1)
yi = np.arange(ymin, ymax, 1)

xi, yi = np.meshgrid(xi,yi)

zi =griddata((x, y), z, (xi, yi))

heatmap= plt.contourf(xi,yi,zi,100, locator=ticker.LogLocator(), cmap=plt.cm.get_cmap('gist_rainbow'),alpha=0.4)
lines= plt.contour(xi,yi,zi,5,colors='k')
plt.clabel(lines,fontsize=18, inline=1, colors='k')
plt.colorbar(heatmap)
plt.show()