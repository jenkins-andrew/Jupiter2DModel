import math
import numpy as np
import matplotlib.pyplot as plt


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


X, Y =np.meshgrid(x, y)

fig, ax = plt.subplots()
CS = ax.contour(X,Y,z)
