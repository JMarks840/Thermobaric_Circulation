# -*- coding: utf-8 -*-
"""
Collection of formulas for basic calulations.
The pressure is always taken in bar, the temperature in °C, the conductivity in
mS/cm and the salinity in PSU.

Created on Wed Aug 21 10:47:38 2024

@author: marksj
"""

# %% Imports

# Third party imports
import gsw
import numpy as np
import matplotlib.pyplot as plt

# Local imports
import sound_velocity as sv
import density as dens
# %% Formulas for Basic Calculations


def T_md_appr(p):
    """
    Temperature of maximum density (approximation)

    Input:
        pressure [bar]

    Returns:
        temperature of maximum density [°C]
    """
    p_dbar = p*1e1
    Tmd = 4 - p_dbar*0.0021
    return Tmd


def T_md_CM(p, S):
    """
    Temperature of maximum density after Chen and Millero (1986)

    Input:
        pressure [bar], salinity [PSU]

    Returns:
        temperature of maximum density [°C]
    """
    Tmd = 3.9839 - 1.9911e-2*p - 5.822e-6*p**2 - (0.2219 + 1.106e-4*p)*S
    return Tmd


def T_md_teos10(p, T_low, T_up, T_range_num):
    """
    Temperature of maximum density after TEOS-10

    This function uses the GSW-Python package available at
    https://github.com/TEOS-10/GSW-Python
    Reference:
        McDougall, T.J. and P.M. Barker, 2011: Getting started with TEOS-10 and
        the Gibbs Seawater (GSW) Oceanographic Toolbox, 28pp., SCOR/IAPSO
        WG127, ISBN 978-0-646-55621-5

    The temperature range has to include all possible values for the
    temperature of maximum density for the whole depth range (compare with
    T_md_CM if needed).

    Note:
        This is only a numerical approxcimation of the correct calculation
        del(rho_is)/del(T)=0.
        The resolution depends on the input values for the temperature range.
        The conductivity is assumed to be zero for simplicity.

    Input:
        pressure [bar], lowest possible temperature of maximum density [°C],
        highest possible temperature of maximum density [°C], number of steps
        between T_low and T_up

    Returns:
        array of temperature of maximum density for the whole profile [°C]
    """
    T_md_teos10 = np.array([])
    T_range = np.linspace(T_low, T_up, T_range_num)
    for i, p_i in enumerate(p):
        rho_teos10 = dens.rho_teos10(p_i, T_range)
        T_md_teos10_i = T_range[np.argmax(rho_teos10)]
        T_md_teos10 = np.append(T_md_teos10, T_md_teos10_i)
    return T_md_teos10


def T_md(p, lambda0, lambda1, T_low, T_up, T_range_num):
    """
    Temperature of maximum density based on the used in-situ density

    The temperature range has to include all possible values for the
    temperature of maximum density for the whole depth range (compare with
    T_md_CM if needed).

    Note:
        This is only a numerical approxcimation of the correct calculation
        del(rho_is)/del(T)=0.
        The resolution depends on the input values for the temperature range.
        The conductivity is assumed to be zero for simplicity.

    Input:
        pressure [bar], lambda_0 [kg*cm/(m^3*mS)], lambda_1 [kg*cm/(m^3*mS*K)],
        lowest possible temperature of maximum density [°C], highest possible
        temperature of maximum density [°C], number of steps between T_low and
        T_up

    Returns:
        array of temperature of maximum density for the whole profile [°C]
    """
    T_md = np.array([])
    T_range = np.linspace(T_low, T_up, T_range_num)
    k25 = 0
    for i, p_i in enumerate(p):
        c = sv.Belogolskii(p_i, T_range)
        rho_is = dens.rho_insitu(p_i, T_range, k25, c, lambda0, lambda1)
        T_md_i = T_range[np.argmax(rho_is)]
        T_md = np.append(T_md, T_md_i)
    return T_md


def depth(p):
    """
    Transformation of pressure into depth (approximation)

    Input:
        pressure [bar]

    Returns:
        depth [m]
    """
    p_dbar = p*1e1
    return -p_dbar/0.98


def k25(T, C, alpha):
    """
    Calculation of the conductivity at 25 °C

    Input:
        temperature [°C], conductivity [mS/cm], alpha as part of the lake
        parameters (LakeParameters) [1/K]

    Returns:
        conductivity at 25 °C [mS/cm]
    """
    return C/(alpha*(T - 25) + 1)


def C_from_k25(T, k25, alpha):
    """
    Calculation of the conductivity from the conductivity at 25 °C (reverse of
    k25(T, C, alpha))

    Input:
        temperature [°C], conductivity at 25 °C [mS/cm], alpha as part of the
        lake parameters (LakeParameters) [1/K]

    Returns:
        conductivity [mS/cm]
    """
    return k25*(alpha*(T - 25) + 1)


def UNESCO(p, T, C):
    """
    UNESCO formula to transfrom the measured conductivity into practival
    salinity units (PSU) (implemented after Boeherer and Schultze (2008))

    Input:
        pressure [bar], temperature [°C], conductivity [mS/cm]

    Returns:
        salinity [PSU]
    """
    a = np.array([0.0080, -0.1692, 25.3851, 14.0941, -7.0261, 2.7081])
    b = np.array([0.0005, -0.0056, -0.0066, -0.0375, 0.0636, -0.0144])
    c = np.array([0.6766097, 2.00564e-2, 1.104259e-4, -6.9698e-7,
                  1.0031e-9])
    d = np.array([1, 3.426e-2, 4.464e-4, 0.4215, -3.107e-3])
    e = np.array([0, 2.070e-5, -6.370e-10, 3.989e-15])
    k = 0.0162
    R = C/42.914
    r_t = c[0] + c[1]*T + c[2]*T**2 + c[3]*T**3 + c[4]*T**4
    R_p = (1 + ((e[0] + e[1]*p + e[2]*p**2 + e[3]*p**3)/(d[0] + d[1]*T +
           d[2]*T**2 + (d[3] + d[4]*T)*R)))
    R_t = R/(R_p*r_t)
    S_oc1 = (a[0] + a[1]*R_t**(1/2) + a[2]*R_t**(2/2) + a[3]*R_t**(3/2) +
             a[4]*R_t**(4/2) + a[5]*R_t**(5/2))
    S_oc2 = (b[0] + b[1]*R_t**(1/2) + b[2]*R_t**(2/2) + b[3]*R_t**(3/2) +
             b[4]*R_t**(4/2) + b[5]*R_t**(5/2))
    return S_oc1 + ((T - 15)/(1 + k*(T - 15)))*S_oc2


def N2_teos10(depth, p, T):
    """
    Brunt-Väisälä-Frequency based on the in-situ density after TEOS-10

    Input:
        pressure [bar], temperature [°C]

    Returns:
        Brunt-Väisälä-Frequency [1/s^2] as array with the last element value
        -9999 to ensure the same size as the input arrays

    Note:
        The input for depth, pressure, and temperature has to be a numpy array
        with at least 2 elements each and the same length.
        The last value will be -9999 to ensure the same length as the input
        since the output will be always one values shorter due to the
        calculations.
    """
    N2_column = []
    index = np.linspace(1, (len(depth) - 1), (len(depth) - 1)
                        ).astype("int64")
    for i in index:
        rho_1 = dens.rho_teos10(p[i], T[i-1])
        rho_2 = dens.rho_teos10(p[i], T[i])
        z_1 = depth[i-1]
        z_2 = depth[i]
        N2 = -(9.81/rho_1)*((rho_2 - rho_1)/(z_2 - z_1))
        N2_column.append(N2)
    N2_column = np.array(N2_column)
    N2_column = np.append(N2_column, -9999)
    return N2_column

def N2_pot(depth, T, k25, lambda0, lambda1):
    """
    Brunt-Väisälä-Frequency using the potential density

    Input:
        depth [m], temperature [°C], conductivity at 25 °C [mS/cm],
        lambda_0 [kg*cm/(m^3*mS)], lambda_1 [kg*cm/(m^3*mS*K)]

    Returns:
        Brunt-Väisälä-Frequency [1/s^2] as array with the last element value
        -9999 to ensure the same size as the input arrays

    Note:
        The input for depth, temperature and conductivity has to be a numpy
        array with at least 2 elements each and the same length.
        The last value will be -9999 to ensure the same length as the input
        since the output will be always one values shorter due to the
        calculations.
    """
    N2_column = []
    index = np.linspace(1, (len(depth) - 1), (len(depth) - 1)
                        ).astype("int64")
    for i in index:
        rho_1 = dens.Moreira(T[i-1], k25[i-1], lambda0, lambda1)
        rho_2 = dens.Moreira(T[i], k25[i], lambda0, lambda1)
        z_1 = depth[i-1]
        z_2 = depth[i]
        N2 = -(9.81/rho_1)*((rho_2 - rho_1)/(z_2 - z_1))
        N2_column.append(N2)
    N2_column = np.array(N2_column)
    N2_column = np.append(N2_column, -9999)
    return N2_column

def N2_is(depth, p, T, k25, c, lambda0, lambda1):
    """
    Brunt-Väisälä-Frequency using the in-situ density

    Calculating the stability comparing water parcels at the same pressure.

    Input:
        depth [m], pressure [bar], temperature [°C], conductivity at 25 °C
        [mS/cm], sound velocity [m/s], lambda_0 [kg*cm/(m^3*mS)],
        lambda_1 [kg*cm/(m^3*mS*K)]

    Returns:
        Brunt-Väisälä-Frequency [1/s^2] as array with the last element value
        -9999 to ensure the same size as the input arrays

    Note:
        The input for depth, pressure, temperature and conductivity has to be
        a numpy array with at least 2 elements each and the same length.
        The last value will be -9999 to ensure the same length as the input
        since the output will be always one values shorter due to the
        calculations.
    """
    N2_column = []
    index = np.linspace(1, (len(depth) - 1), (len(depth) - 1)
                        ).astype("int64")
    for i in index:
        rho_1 = dens.rho_insitu(p[i], T[i-1], k25[i-1],
                                sv.Belogolskii(p[i], T[i-1]), lambda0, lambda1)
        rho_2 = dens.rho_insitu(p[i], T[i], k25[i], c[i], lambda0, lambda1)
        z_1 = depth[i-1]
        z_2 = depth[i]
        N2 = -(9.81/rho_1)*((rho_2 - rho_1)/(z_2 - z_1))
        N2_column.append(N2)
    N2_column = np.array(N2_column)
    N2_column = np.append(N2_column, -9999)
    return N2_column


def plot_T_md(p):
    """
    Plotting the approximated temperature of maximum density, see T_md_appr()
    """
    Tmd = T_md_appr(p)
    d = depth(p)
    plt.plot(Tmd, d)
    plt.xlabel("Temperature of maximum density [°C]")
    plt.ylabel("Depth [m]")


def plot_T_md_CM(p, S):
    """
    Plotting temperature of maximum density after Chen and Millero (1986),
    see T_md()
    """
    Tmd = T_md_CM(p, S)
    d = depth(p)
    plt.plot(Tmd, d)
    plt.xlabel("Temperature of maximum density [°C]")
    plt.ylabel("Depth [m]")
