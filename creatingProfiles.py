import numpy as np

def equatorialMagneticField(r, phi):
    """
    Finds the equatorial magnetic field strength using Vogt et. all 2011 method
    :param r: The radius in R_J
    :param phi: The angle in radians, 0 at the Sun, anti-clockwise
    :return: The equatorial magnetic field in nT
    """
    B = 1.030e6 * r ** (-3.756 - 0.12 * np.cos(phi - 3.562)) + \
        (3.797 - 4.612 * np.cos(phi - 0.825) + 0.606 * np.cos(2 * (phi - 0.473)) +
         0.847 * np.cos(3 * (phi - 0.913))) * np.exp((-1 * r) / 150)
    return B


def plasmaNumberDensity(r, phi):
    """
    Calculates the plasma density at the equator using the method from Bagenal 2011
    :param r: The radius in R_J
    :param phi: The angle in radians, 0 at the Sun, anti-clockwise
    :return: The plasma number density at the equator in cm^-3
    """
    N = 1987 * (r / 6) ** (-8.2) + 14 * (r / 6) ** (-3.2) + 0.05 * (r / 6) ** (-0.65)
    return N


def alfvenVelocityFunc(b, n):
    """
    Calculates the Alfven velocity for a given magnetic field strength and number density
    :param b: The equatorial magnetic field in nT
    :param n: The plasma number density at the equator in cm^-3
    :return: The Alfven Velocity in m/s
    """
    Va = b * 1e-9 / np.sqrt(1.25663706212e-6 * n * 1e6 * 20 * 1.67e-27)
    return Va


def corotationVelocityFunc(x, y):
    """
    Calculate the corotational velocity at x and y
    :param x: In R_J
    :param y: In R_J
    :return: The corotation velocity assuming 10 hr rotation period in m/s
    """
    v = (2*np.pi/36000) * np.sqrt(x**2 + y**2) * 71492e3
    return v


def radialScaleHeight(r):
    """
    Finds the scale height at a radius
    :param r: Radius in R_J
    :return: Scale height in R_J
    """
    h = -0.116 + 2.14*np.log10(r/6) - 2.05*np.log10(r/6)**2 + 0.491*np.log10(r/6)**3 + 0.126*np.log10(r/6)**4
    H = 10 ** h
    return H


# Create a series of arrays to hold values
x = []
y = []
equatorialMagField = []
numberDensity = []
scaleHeight = []
radius = []
alfvenVelocity = []
corotationVelocity = []
corotationcheck = []

# Calculate radius, scale height, x, y, equatorial magnetic field and number density by iterating over radius and angle
for r in np.arange(1, 150, 0.5):
    radius.append(r)
    scaleHeight.append(radialScaleHeight(r))
    for phi in np.arange(0, 2 * np.pi + 0.03, 0.05):
        x.append(r * np.cos(phi))
        y.append(r * np.sin(phi))
        equatorialMagField.append(equatorialMagneticField(r, phi))
        numberDensity.append(plasmaNumberDensity(r, phi))

# Calculate Alfven and corotational velocity
for i in range(len(x)):
    alfvenVelocity.append(alfvenVelocityFunc(equatorialMagField[i], numberDensity[i]))
    corotationVelocity.append(corotationVelocityFunc(x[i], y[i]))

# Check if Alfven velocity is greater than corotational, if so set a binary choice to 0
# will be used to create a plot later
for i in range(len(alfvenVelocity)):
    if alfvenVelocity[i] > corotationVelocity[i]:
        corotationcheck.append(0)
    else:
        corotationcheck.append(1)

# Save outputs
np.savetxt('alfvenCheck.txt', np.c_[x, y, equatorialMagField, numberDensity, alfvenVelocity, corotationVelocity,
                                    corotationcheck], delimiter='\t', header='x\ty\tb\tp\tAlfven\tCorotation\tCheck')
np.savetxt('scaleheighttest.txt', np.c_[radius, scaleHeight], delimiter='\t', header='r\tscaleHeight')