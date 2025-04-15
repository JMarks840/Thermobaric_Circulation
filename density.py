# -*- coding: utf-8 -*-
"""
Collection of formulas for the calculation of the density based on different
approaches.
The pressure is always taken in bar, the temperature in °C, the conductivity in
mS/cm and the salinity in PSU.

@author: Joshua Marks
"""

# %% Imports

# Local imports
import sound_velocity as sv

# %% Formulas for Denisty Calculations

# First the formulas for the potential density


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


# Second the formula for the in-situ density


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
    """
    p_Pa = p*1e5  # bar in Pa
    rho_pot = Moreira(T, k25, lambda0, lambda1)
    # gradient of del(1/c^2)/del(p) [s^2/m^2/bar]
    del1c2_delp = (((1/sv.Belogolskii(36, T)**2) -
                    (1/sv.Belogolskii(0, T)**2))/(36 - 0))
    # 1/c^2 at surface pressure [s^2/m^2]
    OneOverc2_0 = (1/c**2) - (del1c2_delp*p)
    return rho_pot + ((p_Pa**2/2)*del1c2_delp*1e-5) + (p_Pa*OneOverc2_0)
