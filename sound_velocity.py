# -*- coding: utf-8 -*-
"""
Formula for the calculation of the sound velocity.
The pressure is always taken in bar, the temperature in °C, the conductivity in
mS/cm and the salinity in PSU.

@author: Joshua Marks
"""

# %% Imports

# Third party imports
import numpy as np

# %% Formulas for Sound Velocity Calculations


def Belogolskii(p, T):
    """
    Sound velocity after Belogol'skii et al. (1999)

    Input:
        pressure [bar], temperature [°C]

    Returns:
        sound velocity [m/s]
    """
    p_MPa = p*0.1  # bar in MPa
    # a[i] corresponds to a_{i 0} in paper notation
    a = np.array([1402.38744, 5.03836171, -5.81172916e-2, 3.34638117e-4,
                  -1.48259672e-6, 3.16585020e-9])
    # aij[i, j] corresponds to a_{i j+1} in paper notation
    # i = line   (0-3 in array is equal to first number in paper notation)
    # j = column (0-2 in array corresponds to 1-3 as second number in paper
    #             notaion)
    # array in paper notation:
    #              [[ a01,             a02,             a03],
    #               [ a11,             a12,             a13],
    #               [ a21,             a22,             a23],
    #               [ a31,             a32,             a33]]
    aij = np.array([[1.49043589,      4.31532833e-3,  -1.852993525e-5],
                    [1.077850609e-2, -2.938590293e-4,  1.481844713e-6],
                    [-2.232794656e-4, 6.822485943e-6, -3.940994021e-8],
                    [2.718246452e-6, -6.674551162e-8,  3.939902307e-10]])
    M1 = aij[0, 0] + aij[1, 0]*T + aij[2, 0]*T**2 + aij[3, 0]*T**3
    M2 = aij[0, 1] + aij[1, 1]*T + aij[2, 1]*T**2 + aij[3, 1]*T**3
    M3 = aij[0, 2] + aij[1, 2]*T + aij[2, 2]*T**2 + aij[3, 2]*T**3
    W0 = a[0] + a[1]*T + a[2]*T**2 + a[3]*T**3 + a[4]*T**4 + a[5]*T**5
    W1 = M1*(p_MPa - 0.101325)
    W2 = M2*(p_MPa - 0.101325)**2
    W3 = M3*(p_MPa - 0.101325)**3
    return W0 + W1 + W2 + W3
