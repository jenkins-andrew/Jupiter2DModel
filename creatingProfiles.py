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


def equatorialPlasmaNumberDensity(r, speciesValues=None):
    """
    Calculates the plasma density at the equator using the method from Bagenal 2011
    :param r: The radius in R_J
    :param speciesValues:
    :return: The plasma number density at the equator in cm^-3
    """
    b2011 = 1987 * (r / 6) ** (-8.2) + 14 * (r / 6) ** (-3.2) + 0.05 * (r / 6) ** (-0.65)

    try:
        percentage, a, b, c = speciesValues
        if r <= 15.2:
            n = a * (r / 6) ** b
        else:
            n = (c * b2011)
    except:
        n = b2011
    return n


def equatorialAveragePlasmaNumberDensity(r, species):
    """
    Calculates the plasma density at the equator using the method from Bagenal 2011
    :param r: The radius in R_J
    :param species:
    :return: The plasma number density at the equator in cm^-3
    """
    b2011 = 1987 * (r / 6) ** (-8.2) + 14 * (r / 6) ** (-3.2) + 0.05 * (r / 6) ** (-0.65)
    n = []
    for i in species:
        percentage, a, b, c = species[i]
        if r <= 15.2:
            n.append(a * (r / 6) ** b)
        else:
            n.append(c * b2011)

    averageN = np.mean(n)

    return averageN


def averageAmu(r, species, massAmuArray):
    """

    :param r:
    :param species:
    :param massAmuArray:
    :return:
    """
    sumofmasses = 0
    N = []
    for i in massAmuArray:
        mass = massAmuArray[i]
        try:
            n = equatorialPlasmaNumberDensity(r, species[i])
            N.append(n)
        except:
            n = 0
            print('Species do not match')
        sumofmasses += n * mass

    amu = sumofmasses/sum(N)
    return amu


def alfvenVelocityFunc(b, n, amu=20):
    """
    Calculates the Alfven velocity for a given magnetic field strength and number density
    :param b: The equatorial magnetic field in nT
    :param n: The plasma number density at the equator in cm^-3
    :param amu:
    :return: The Alfven Velocity in m/s
    """
    Va = b * 1e-9 / np.sqrt(1.25663706212e-6 * n * 1e6 * amu * 1.67e-27)
    return Va


def corotationVelocityFunc(x, y):
    """
    Calculate the corotational velocity at x and y
    :param x: In R_J
    :param y: In R_J
    :return: The corotation velocity assuming 10 hr rotation period in m/s
    """
    v = (2*np.pi/(3600*9.9250)) * np.sqrt(x**2 + y**2) * 71492e3
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


def densityAtZFromEquator(z, r, species):
    """

    :param z:
    :param r:
    :param species:
    :return:
    """

    nZ = equatorialAveragePlasmaNumberDensity(r, species) * np.exp(-1 * (z / radialScaleHeight(r))**2)
    return nZ


def radialVelocityFunc(r, species, massArray):
    """

    :param r:
    :param species:
    :param massArray:
    :return:
    """
    vr = 1000/(2 * equatorialAveragePlasmaNumberDensity(r, species)*1e6 * averageAmu(r, species, massArray)*1.67e-27 *
               radialScaleHeight(r) * np.pi * r * 71492e3**2)
    return vr


# Create a series of arrays to hold values
xInRJ = []
yInRJ = []
zInRJ = []
equatorialMagField = []
numberDensity = []
scaleHeight = []
radius = []
alfvenVelocity = []
corotationVelocity = []
corotationcheck = []
plasmaZDensity = []
radiusForZDensity = []
amuAtR = []
radialVelocity = []
radialVelocityAtPi = []
alfvenVelocityATPi = []

# No longer have to be in the same order
speciesList = {'e-': [1, 2451, -6.27, 4.21],
               'o+': [0.24, 592, -7.36, 0.368],
               'o++': [0.03, 76.3, -6.73, 0.086],
               's+': [0.07, 163, -6.81, 0.169],
               's++': [0.22, 538, -6.74, 0.598],
               's+++': [0.004, 90.7, -6.21, 0.165],
               'h+': [0.02, 50.6, -5.31, 0.212],
               'na+': [0.04, 97.2, -6.75, 0.106],
               'hoto+': [0.06, 134, -4.63, 1.057]}
ME = 0.00054858
speciesMass = {'e-': 0.00054858,
               'o++': 15.999 - (ME * 2),
               's+': 32.065 - ME,
               's++': 32.065 - (ME * 2),
               's+++': 32.065 - (ME * 3),
               'h+': 1.00784 - ME,
               'na+': 22.989769 - ME,
               'hoto+': 15.999 - (ME * 2),
               'o+': 15.999 - ME
               }

# Calculate radius, scale height, x, y, equatorial magnetic field and number density by iterating over radius and angle
for r in np.arange(6, 100, 0.5):
    radius.append(r)
    scaleHeight.append(radialScaleHeight(r))
    radialVelocityAtPi.append(radialVelocityFunc(r, speciesList, speciesMass))
    alfvenVelocityATPi.append(alfvenVelocityFunc(equatorialMagneticField(r, 0), equatorialAveragePlasmaNumberDensity(r, speciesList), averageAmu(r, speciesList, speciesMass)))
    for phi in np.arange(0, 2 * np.pi + 0.03, 0.05):
        xInRJ.append(r * np.cos(phi))
        yInRJ.append(r * np.sin(phi))
        equatorialMagField.append(equatorialMagneticField(r, phi))
        numberDensity.append(equatorialAveragePlasmaNumberDensity(r, speciesList))
        amuAtR.append(averageAmu(r, speciesList, speciesMass))
        radialVelocity.append(radialVelocityFunc(r, speciesList, speciesMass))

# Calculate Alfven and corotational velocity
for i in range(len(xInRJ)):
    alfvenVelocity.append(alfvenVelocityFunc(equatorialMagField[i], numberDensity[i], amuAtR[i]))
    corotationVelocity.append(corotationVelocityFunc(xInRJ[i], yInRJ[i]))

# Check if Alfven velocity is greater than corotational, if so set a binary choice to 0
# will be used to create a plot later
for i in range(len(alfvenVelocity)):
    if alfvenVelocity[i] > radialVelocity[i]:
        corotationcheck.append(0)
    else:
        corotationcheck.append(1)

for r in np.arange(6, 100, 0.5):
    for z in np.arange(-12, 12, 0.1):
        radiusForZDensity.append(r)
        zInRJ.append(z)
        plasmaZDensity.append(densityAtZFromEquator(z, r, speciesList))

# Save outputs
np.savetxt('alfvenCheck.txt', np.c_[xInRJ, yInRJ, equatorialMagField, numberDensity, alfvenVelocity, radialVelocity,
                                    corotationcheck], delimiter='\t', header='x\ty\tb\tp\tAlfven\tCorotation\tCheck')
np.savetxt('scaleheighttest.txt', np.c_[radius, scaleHeight], delimiter='\t', header='r\tscaleHeight')

np.savetxt('zPlasmaDensity.txt', np.c_[radiusForZDensity, zInRJ, plasmaZDensity], delimiter='\t', header='r\tz\tplasmaZDensity')

np.savetxt('alfvenradial.txt', np.c_[radius, alfvenVelocityATPi, radialVelocityAtPi], delimiter='\t', header='r\tscaleHeight')