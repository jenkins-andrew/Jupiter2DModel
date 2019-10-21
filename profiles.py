import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker, cm
from scipy.ndimage.filters import gaussian_filter
from scipy.interpolate import griddata
from scipy import interpolate


def equatorialMagneticField(r, phi):
    B = 1.030e6 * r ** (-3.756 - 0.12 * np.cos(phi - 3.562)) + \
        (3.797 - 4.612 * np.cos(phi - 0.825) + 0.606 * np.cos(2 * (phi - 0.473)) +
         0.847 * np.cos(3 * (phi - 0.913))) * np.exp(-1 * r / 150)
    return B


def plasmaDensity(r, phi):
    N = 1987 * (r / 6) ** (-8.2) + 14 * (r / 6) ** (-3.2) + 0.05 * (r / 6) ** (-0.65)
    return N


#file = open('Test.txt', 'w')

xMagnetic = np.arange(20, 150, 1)
xPlasma = np.arange(6, 100, 1)
phi = np.arange(0, 2 * np.pi + 0.03, 0.05)

xp, phip = np.meshgrid(xPlasma, phi)
xi, phii = np.meshgrid(xMagnetic, phi)

B = np.array(equatorialMagneticField(xi, phii))
N = np.array(plasmaDensity(xp, phip))
# z = gaussian_filter(z,4,mode='nearest')

xi = xMagnetic * np.cos(phii)
yi = xMagnetic * np.sin(phii)

xp = xPlasma * np.cos(phip)
yp = xPlasma * np.sin(phip)

np.savetxt('test.txt', (xi))

print(xi.shape)

# Plotting
'''plt.subplot(211)
heatmap = plt.contourf(xi, yi, B, cmap=plt.cm.get_cmap('gist_rainbow'), alpha=0.4)
lines = plt.contour(xi, yi, B, 5, colors='k')
plt.clabel(lines, fontsize=18, inline=1, colors='k')
clb = plt.colorbar(heatmap)
clb.ax.set_title('B$_n$ (nT)', fontsize=18)
plt.title('Magnetic field of Jupiter', fontsize=18)
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
plt.xlabel('x $(R_J)$', fontsize=18)
plt.ylabel('y $(R_J)$', fontsize=18)
plt.xticks(size=18)
plt.yticks(size=18)
plt.text(60, -70, r'$\rightarrow$ To the Sun')
'''
#plt.subplot(212)
heatmap = plt.contourf(xp, yp, N, locator=ticker.LogLocator(), cmap=plt.cm.get_cmap('gist_rainbow'), alpha=0.4)
lines = plt.contour(xp, yp, N, 5, colors='k')
plt.clabel(lines, fontsize=18, inline=1, colors='k')
#plt.plot(xp,yp,'k.')
clb = plt.colorbar(heatmap)
plt.tight_layout()
clb.ax.set_title(r'$\rho$ (cm$^{-3}$)', fontsize=18)
plt.title('Plasma density of Plasma sheet', fontsize=18)
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
plt.xlabel('x $(R_J)$', fontsize=18)
plt.ylabel('y $(R_J)$', fontsize=18)
plt.xticks(size=18)
plt.yticks(size=18)
plt.text(10, -10, r'$\rightarrow$ To the Sun')
plt.savefig('plasmaDensity.pdf', bbox_inches='tight')
