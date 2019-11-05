import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker, cm
from scipy.ndimage.filters import gaussian_filter
from scipy.interpolate import griddata
from scipy import interpolate

def equatorialMagneticField(r, phi):
    """
    Finds the equatorial magnetic field strength using Vogt et. all 2011 method
    :param r: The radius
    :param phi: The angle
    :return: The equatorial magnetic field
    """
    B = 1.030e6 * r ** (-3.756 - 0.12 * np.cos(phi - 3.562)) + \
        (3.797 - 4.612 * np.cos(phi - 0.825) + 0.606 * np.cos(2 * (phi - 0.473)) +
         0.847 * np.cos(3 * (phi - 0.913))) * np.exp((-1 * r) / 150)
    return B


def plasmaDensity(r, phi):
    N = 1987 * (r / 6) ** (-8.2) + 14 * (r / 6) ** (-3.2) + 0.05 * (r / 6) ** (-0.65)
    return N


def alfvenVelocityFunc(b, rho):
    Va = b * 1e-9 / np.sqrt(1.25663706212e-6 * rho * 1e6 * 1.67e-27)
    return Va


def corotationVelocityFunc(x, y):
    x = x*69911e3
    y = y*69911e3
    v = (2*np.pi/36000) * np.sqrt(x**2 + y**2)
    return v

def radialScaleHeight(r):
    h = -0.116 + 2.14*np.log10(r) - 2.05*np.log10(r)**2 + 0.491*np.log10(r)**3 + 0.126*np.log10(r)**4
    H = 10 ** h
    return H


x = []
y = []
b = []
p = []
scaleHeight = []
radius = []

for r in np.arange(1, 150, 1):
    radius.append(r)
    scaleHeight.append(radialScaleHeight(r))
    for phi in np.arange(0, 2 * np.pi + 0.03, 0.05):
        x.append(r * np.cos(phi))
        y.append(r * np.sin(phi))
        b.append(equatorialMagneticField(r, phi))
        p.append(plasmaDensity(r, phi))

alfvenVelocity = []
corotationVelocity = []

for i in range(len(x)):
    alfvenVelocity.append(alfvenVelocityFunc(b[i], p[i]))
    corotationVelocity.append(corotationVelocityFunc(x[i], y[i]))

corotationcheck = []

for i in range(len(alfvenVelocity)):
    if alfvenVelocity[i] > corotationVelocity[i]:
        corotationcheck.append(0)
    else:
        corotationcheck.append(1)

np.savetxt('alfvenCheck.txt', np.c_[x, y,  b, p, alfvenVelocity, corotationVelocity, corotationcheck], delimiter='\t', header='x\ty\tb\tp\tAlfven\tCorotation\tCheck')
np.savetxt('scaleheighttest.txt', np.c_[radius, scaleHeight], delimiter='\t', header='r\tscaleHeight')