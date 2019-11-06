import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from matplotlib import colors

# Loading the data
x, y, b, p, alfvenVelocity, corotation, corotationcheck = np.loadtxt('alfvenCheck.txt', delimiter='\t', unpack=True)
r, scaleHeight = np.loadtxt('scaleheighttest.txt', delimiter='\t', unpack=True)

# Creating grid
xtest = np.arange(-100, 101, 1)
ytest = np.arange(-100, 101, 1)
xtest, ytest = np.meshgrid(xtest, ytest)

# Masking a circle of radius 20 R_J
mask = (np.sqrt(xtest ** 2 + ytest ** 2) <= 6) | (np.sqrt(xtest ** 2 + ytest ** 2) >= 100)

# Making the 3D grid for the magnetic field
BGrid = griddata((x, y), b, (xtest, ytest), method='linear')
BGrid[mask] = np.nan

# Making the 3D grid for the plasma density
PGrid = griddata((x, y), p, (xtest, ytest), method='linear')
PGrid[mask] = np.nan

# Making the 3D grid for the Alfven Velocity
AVGrid = griddata((x, y), alfvenVelocity, (xtest, ytest), method='linear')
AVGrid[mask] = np.nan

# Plotting
# plt.subplot(121)
# heatmap = plt.contourf(xtest, ytest, BGrid, cmap=plt.cm.get_cmap('gist_rainbow'), alpha=0.4)
# lines = plt.contour(xtest, ytest, BGrid, 5, colors='k')
# plt.clabel(lines, fontsize=18, inline=1, colors='k')
# clb = plt.colorbar(heatmap)
# clb.ax.set_title('B$_n$ (nT)', fontsize=18)
# plt.title('Magnetic field of Jupiter', fontsize=18, wrap=True)
# plt.rcParams['xtick.labelsize'] = 18
# plt.rcParams['ytick.labelsize'] = 18
# plt.xlabel('x $(R_J)$', fontsize=18)
# plt.ylabel('y $(R_J)$', fontsize=18)
# #plt.axis(xlim=(np.amin(xtest), np.amax(xtest)), ylim=(np.amin(ytest), np.amax(ytest)))
# plt.xticks(size=18)
# plt.yticks(size=18)
# plt.text(10, -70, r'$\rightarrow$ To the Sun')
#
# plt.subplot(122)
# heatmap = plt.contourf(xtest, ytest, PGrid, cmap=plt.cm.get_cmap('gist_rainbow'), alpha=0.4)
# lines = plt.contour(xtest, ytest, PGrid, 5, colors='k')
# plt.clabel(lines, fontsize=18, inline=1, colors='k')
# clb = plt.colorbar(heatmap)
# clb.ax.set_title(r'$\rho$ (cm$^{-3}$)', fontsize=18)
# plt.title('Plasma density at Jupiter', fontsize=18, wrap=True)
# plt.rcParams['xtick.labelsize'] = 18
# plt.rcParams['ytick.labelsize'] = 18
# plt.xlabel('x $(R_J)$', fontsize=18)
# plt.ylabel('y $(R_J)$', fontsize=18)
# #plt.axis(xlim=(np.amin(xtest), np.amax(xtest)), ylim=(np.amin(ytest), np.amax(ytest)))
# plt.xticks(size=18)
# plt.yticks(size=18)
# plt.text(10, -70, r'$\rightarrow$ To the Sun')
#
# plt.tight_layout(pad=1)
# plt.savefig('contourwithmask.pdf', bbox_inches='tight')

# plt.figure()
# cmap = colors.ListedColormap(['red', '#ffffff'])
# boundaries = [0, 1]
# norm = colors.BoundaryNorm(boundaries, cmap.N, clip=True)
# CheckGrid = griddata((x, y), corotationcheck, (xtest, ytest), method='linear')
# CheckGrid[mask] = np.nan
#
# #heatmap = plt.contourf(xtest, ytest, CheckGrid, cmap=cmap)
# heatmap = plt.contourf(xtest, ytest, CheckGrid, cmap=plt.cm.get_cmap('gist_rainbow'), alpha=0.4)
# lines = plt.contour(xtest, ytest, CheckGrid, 1, colors='k')
# plt.xlabel('x $(R_J)$', fontsize=18)
# plt.ylabel('y $(R_J)$', fontsize=18)
# plt.xticks(size=18)
# plt.yticks(size=18)
# plt.text(10, -70, r'$\rightarrow$ To the Sun')
# plt.tight_layout()

plt.figure()
alfvenmask = (alfvenVelocity > 0.95*corotation) & (alfvenVelocity < 1.05*corotation)
calculatedRadius = np.sqrt(x ** 2 + y ** 2)
phiwrong = np.arctan2(x, y)
phi = np.mod(phiwrong, 2*np.pi) * 180 / np.pi
fit = np.poly1d(np.polyfit(phi[alfvenmask], calculatedRadius[alfvenmask], 3))
plt.scatter(phi[alfvenmask], calculatedRadius[alfvenmask], s=0.1, color='k')
fitrange = np.arange(0, 360, 1)
plt.plot(fit(fitrange))
plt.xlabel('Angle (Degrees)', fontsize=18)
plt.ylabel('Radius (R$_J)$', fontsize=18)
plt.xticks(size=18)
plt.yticks(size=18)
plt.tight_layout()

# plt.figure()
# heatmap = plt.contourf(xtest, ytest, AVGrid, cmap=plt.cm.get_cmap('gist_rainbow'), alpha=0.4)
# lines = plt.contour(xtest, ytest, AVGrid, 5, colors='k')
# plt.clabel(lines, inline=1, colors='k')
# clb = plt.colorbar(heatmap)
# clb.ax.set_title(r'V$_A$ (ms$^{-1}$)', fontsize=18)
# plt.title('Alfven Velocity at Jupiter', fontsize=18, wrap=True)
# plt.rcParams['xtick.labelsize'] = 18
# plt.rcParams['ytick.labelsize'] = 18
# plt.xlabel('x $(R_J)$', fontsize=18)
# plt.ylabel('y $(R_J)$', fontsize=18)
# #plt.axis(xlim=(np.amin(xtest), np.amax(xtest)), ylim=(np.amin(ytest), np.amax(ytest)))
# plt.xticks(size=18)
# plt.yticks(size=18)
# plt.text(10, -70, r'$\rightarrow$ To the Sun')
# plt.tight_layout()

# plt.figure()
# plt.plot(r, scaleHeight)
# plt.axis(xmax=100, ymax=scaleHeight[100]+0.2)
# plt.xlabel('Radius (R$_J)$', fontsize=18)
# plt.ylabel('Scale Height (R$_J)$', fontsize=18)
# plt.xticks(size=18)
# plt.yticks(size=18)
# plt.tight_layout()

plt.show()