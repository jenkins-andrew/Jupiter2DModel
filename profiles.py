import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata


def equatorialMagneticField(r, phi):
    B = 1030000.0 * r ** (-3.756 - 0.12 * math.cos(phi - 3.562)) + \
(3.797 - 4.612 * math.cos(phi - 0.825) + 0.606 * math.cos(2 * (phi - 0.473)) +
         0.847 * math.cos(3 * (phi - 0.913))) * math.exp(-1 * r / 150)
    return B


x = []
y = []
z = []

for r in np.arange(-150, 150):
    for phi in np.arange(0, 2 * math.pi):
        x.append(r * math.cos(phi))
        y.append(r * math.sin(phi))
        z.append(equatorialMagneticField(r, phi))

x=np.array(x)
y=np.array(y)
z=np.array(z)
xmin, xmax = -100, 100
ymin, ymax = -100, 100
xi = np.arange(xmin, xmax, 1)
yi = np.arange(ymin, ymax, 1)

xi, yi = np.meshgrid(xi,yi)

zi =griddata((x, y), z, (xi, yi), method='nearest')

fig, ax = plt.subplots()

CS = ax.contour(x,y,zi)
