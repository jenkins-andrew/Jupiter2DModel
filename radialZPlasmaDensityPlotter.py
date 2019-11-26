import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker, cm
from scipy.interpolate import griddata


r, z, plasmaDensity = np.loadtxt('zPlasmaDensity.txt', delimiter='\t', unpack=True)
x, y, b, p, alfvenVelocity, corotation, corotationcheck = np.loadtxt('alfvenCheck.txt', delimiter='\t', unpack=True)
xB, yB, zB, BB = np.loadtxt('plotmagfieldlines.txt', delimiter='\t', unpack=True)


maxR = 30
xtest = np.arange(6, maxR+1, 0.5)
ytest = np.arange(np.amin(z), np.amax(z), 0.1)
xtest, ytest = np.meshgrid(xtest, ytest)

theta = -10*np.pi/180
xB = xB*np.cos(theta) - zB*np.sin(theta)
zB = xB*np.sin(theta) + zB*np.cos(theta)

mask = (plasmaDensity > 10**(-5))
maskedFieldLines = (np.abs(xB) > 6)

DensityGrid = griddata((r[mask], z[mask]), plasmaDensity[mask], (xtest, ytest))
FieldGrid = griddata((xB, zB), BB, (xtest, ytest))



plt.figure()
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
heatmap = plt.contourf(xtest, ytest, DensityGrid, cmap=plt.cm.get_cmap('gist_rainbow'), locator=ticker.LogLocator(), alpha=0.4)
plt.plot(xB[maskedFieldLines], zB[maskedFieldLines], 'k')
# lines = plt.contour(xtest, ytest, FieldGrid, 10, colors='k')
# plt.clabel(lines, fontsize=18, inline=1, colors='k')
clb = plt.colorbar(heatmap)
clb.ax.set_title(r'(cm$^{-3}$)', fontsize=18)
plt.xlabel('R $(R_J)$', fontsize=18)
plt.ylabel('Z $(R_J)$', fontsize=18)
plt.xticks(size=18)
plt.yticks(size=18)
plt.ylim(np.amin(z), np.amax(z))

plt.tight_layout()

plt.show()