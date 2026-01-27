# -*- coding: utf-8 -*-
"""
Collection of formulas for the calculation of the density based on different
approaches.
The pressure is always taken in bar, the temperature in °C, the conductivity in
mS/cm and the salinity in PSU.

Created on Mon Aug 26 14:19:20 2024

@author: marksj
"""

# %% Imports

# Third party imports
import gsw
import numpy as np

# Local imports
import calculations as calcs
import sound_velocity as sv

# %% Formulas for Denisty Calculations

# First the formulas for the potential density


def ChenAndMillero(p, T, k25, alpha):
    """
    Potential density after Chen and Millero (1977, 1986)
    (also in Boehrer and Schultze (2008))

    Input:
        temperature [°C], salinity [PSU]

    Returns:
        potential density [kg/m^3]
    """
    C = calcs.C_from_k25(T, k25, alpha)
    S = calcs.UNESCO(p, T, C)
    a = np.array([999.8395, 6.7914e-2, -9.0894e-3, 1.0171e-4, -1.2846e-6,
                  1.1592e-8, -5.0125e-11])
    b = np.array([0.8181,  -3.85e-3,    4.96e-5])
    rho_CM = (a[0] + a[1]*T + a[2]*T**2 + a[3]*T**3 + a[4]*T**4 + a[5]*T**5 +
              S*(b[0] + b[1]*T + b[2]*T**2))
    return rho_CM


def Tanaka(T):
    """
    Potential density after Tanaka et al. (2001)

    Input:
        temperature [°C]

    Returns:
        potential density [kg/m^3]
    """
    a1 = -3.983035  # [°C]
    a2 = 301.797  # [°C]
    a3 = 522528.9  # [°C^2]
    a4 = 69.34881  # [°C]
    a5 = 999.974950  # [kg/m^3]
    return a5*(1 - (((T + a1)**2*(T + a2))/(a3*(T + a4))))


def Moreira(T, k25, lambda0, lambda1):
    """
    Potential density after Moreira et al. (2016)

    Input:
        temperature [°C], conductivity at 25 °C [mS/cm],
        lambda_0 [kg*cm/(m^3*mS)], lambda_1 [kg*cm/(m^3*mS*K)]

    Returns:
        potential density [kg/m^3]
    """
    pot_density = Tanaka(T)
    return pot_density + k25*(lambda0 + lambda1*(T - 25))


# Second the formulas for the in-situ density


def rho_teos10(p, T):
    """
    In-situ density after TEOS-10

    This function uses the GSW-Python package available at
    https://github.com/TEOS-10/GSW-Python
    Reference:
        McDougall, T.J. and P.M. Barker, 2011: Getting started with TEOS-10 and
        the Gibbs Seawater (GSW) Oceanographic Toolbox, 28pp., SCOR/IAPSO
        WG127, ISBN 978-0-646-55621-5

    Input:
        pressure [bar], temperature [°C]

    Returns:
        in-situ density [kg/m^3]

    Note:
        This is so far only implemented for (absolute) salinity equals 0.
    """
    p_dbar = p*10
    SA = 0
    CT = gsw.conversions.CT_from_pt(SA, T)
    rho = gsw.density.rho(SA, CT, p_dbar)
    return rho


def rho_CM(p, T, k25, alpha):
    """
    In-situ density after Chen and Millero (1986)
    (also in Boehrer and Schultze (2008))

    Input:
        pressure [bar], temperature [°C], salinity [PSU]

    Returns:
        potential density [kg/m^3]
    """
    C = calcs.C_from_k25(T, k25, alpha)
    S = calcs.UNESCO(p, T, C)
    c = np.array([19652.17, 148.113, -2.293,    1.256e-2, -4.18e-5])
    d = np.array([3.2726,  -2.147e-4, 1.128e-4])
    e = np.array([53.238,  -0.313])
    f = 5.728e-3
    rho = ChenAndMillero(p, T, k25, alpha)
    K = (c[0] + c[1]*T + c[2]*T**2 + c[3]*T**3 + c[4]*T**4 + p*(d[0] + d[1]*T +
         d[2]*T**2) + S*(e[0] + e[1]*T) + f*p*S)
    return rho/(1-(p/K))


def rho_insitu(p, T, k25, c, lambda0, lambda1):
    """
    In-situ density

    Own formula:
        rho_insitu = rho_pot + (p^2/2)*(del(1/c^2)/del(p)) + p/(c(p=0)^2)

    Input:
        pressure [bar], temperature [°C], conductivity at 25 °C [mS/cm],
        sound velocity [m/s], lambda_0 [kg*cm/(m^3*mS)],
        lambda_1 [kg*cm/(m^3*mS*K)]

    Returns:
        in-situ density [kg/m^3]

    Note:
        The formula is only valid for lakes up to 400 m  due to the used linear
        approximation, but this can be adapted.
    """
    p_Pa = p*1e5  # bar in Pa
    rho_pot = Moreira(T, k25, lambda0, lambda1)
    # gradient of del(1/c^2)/del(p) [s^2/m^2/bar]
    del1c2_delp = (((1/sv.Belogolskii(36, T)**2) -
                    (1/sv.Belogolskii(0, T)**2))/(36 - 0))
    # 1/c^2 at surface pressure [s^2/m^2]
    OneOverc2_0 = (1/c**2) - (del1c2_delp*p)
    return rho_pot + ((p_Pa**2/2)*del1c2_delp*1e-5) + (p_Pa*OneOverc2_0)
