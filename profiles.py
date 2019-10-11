import math
import matplotlib as plt


def equatorialMagneticField(r, phi):
    B = 1030000.0 * r ** (-3.756 - 0.12 * math.cos(phi - 3.562)) + \
        (3.797 - 4.612 * math.cos(phi - 0.825) + 0.606 * math.cos(2 * (phi - 0.473)) +
         0.847 * math.cos(3 * (phi - 0.913))) * math.exp(-1* r / 150)
    return B


print(equatorialMagneticField(10, 3))

