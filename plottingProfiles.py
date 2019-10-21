import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

# Loading the data
x, y, b, p = np.loadtxt('npout.txt', delimiter='\t', unpack=True)

# Creating grid
xtest = np.arange(-100, 101, 1)
ytest = np.arange(-100, 101, 1)
xtest, ytest = np.meshgrid(xtest, ytest)

# Masking a circle of radius 20 R_J
mask = (np.sqrt(xtest ** 2 + ytest ** 2) < 20)

# Making the 3D grid for the magnetic field
BGrid = griddata((x, y), b, (xtest, ytest), method='linear')
BGrid[mask] = np.nan

# Making the 3D grid for the plasma density
PGrid = griddata((x, y), p, (xtest, ytest), method='linear')
PGrid[mask] = np.nan

# Plotting
plt.subplot(121)
heatmap = plt.contourf(xtest, ytest, BGrid, cmap=plt.cm.get_cmap('gist_rainbow'), alpha=0.4)
lines = plt.contour(xtest, ytest, BGrid, 5, colors='k')
plt.clabel(lines, fontsize=18, inline=1, colors='k')
clb = plt.colorbar(heatmap)
clb.ax.set_title('B$_n$ (nT)', fontsize=18)
plt.title('Magnetic field of Jupiter', fontsize=18, wrap=True)
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
plt.xlabel('x $(R_J)$', fontsize=18)
plt.ylabel('y $(R_J)$', fontsize=18)
plt.axis(xlim=(np.amin(xtest), np.amax(xtest)), ylim=(np.amin(ytest), np.amax(ytest)))
plt.xticks(size=18)
plt.yticks(size=18)
plt.text(10, -70, r'$\rightarrow$ To the Sun')

plt.subplot(122)
heatmap = plt.contourf(xtest, ytest, PGrid, cmap=plt.cm.get_cmap('gist_rainbow'), alpha=0.4)
lines = plt.contour(xtest, ytest, PGrid, 5, colors='k')
plt.clabel(lines, fontsize=18, inline=1, colors='k')
clb = plt.colorbar(heatmap)
clb.ax.set_title(r'$\rho$ (cm$^{-3}$)', fontsize=18)
plt.title('Plasma density at Jupiter', fontsize=18, wrap=True)
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
plt.xlabel('x $(R_J)$', fontsize=18)
plt.ylabel('y $(R_J)$', fontsize=18)
plt.axis(xlim=(np.amin(xtest), np.amax(xtest)), ylim=(np.amin(ytest), np.amax(ytest)))
plt.xticks(size=18)
plt.yticks(size=18)
plt.text(10, -70, r'$\rightarrow$ To the Sun')

plt.tight_layout(pad=1)
plt.savefig('contourwithmask.pdf', bbox_inches='tight')

x, y, alfvenVelocity, corotationcheck = np.loadtxt('alfvenCheck.txt', delimiter='\t', unpack=True)

CheckGrid = griddata((x, y), corotationcheck, (xtest, ytest), method='linear')
CheckGrid[mask] = np.nan

plt.contourf(xtest, ytest, CheckGrid)
plt.show()