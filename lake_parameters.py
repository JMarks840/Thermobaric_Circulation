# -*- coding: utf-8 -*-
"""
Lake specific parameters are needed for different calculations.
For the calculation of the conductivity at 25 °C the parameter alpha is needed,
for the calculation of the potential density with Moreira et al. (2016) the
lambda coefficients and for the calculation of the sound velocity with Selyakov
(to be published) the sigma coefficients are needed.
They all depend on the composition of the solutes of the lake water itself.

Created on Mon Aug 26 16:39:03 2024

@author: marksj
"""

# %% Lake Parameters Definition


class LakeParameters:

    def __init__(self, alpha, lambda0, lambda1, sigma0, sigma1):
        """
        Setting the lake specific parameters.
        These are:

            alpha for the calculation of the conductivity at 25 °C

            the lambda values for the calculation of the potential density
            after Moreira et al. (2016)

            and the sigma values for the calculation of the sound velocity
            after Selyakov (to be published)

        Input:
            alpha [1/K], lambda_0 [kg*cm/(m^3*mS)],
            lambda_1 [kg*cm/(m^3*mS*K)], sigma_0 [??], sigma_1 [??]
        """
        self.alpha = alpha
        self.lambda0 = lambda0
        self.lambda1 = lambda1
        self.sigma0 = sigma0
        self.sigma1 = sigma1


# Lake specific parameters for Lake Shikotsu
# ([1/K], [kg*cm/(m^3*mS)], [kg*cm/(m^3*mS*K)], [??], [??])
SHIKOTSU_PARAMETERS = LakeParameters(0.0194, 0.55, -0.0012, 0, 0)
