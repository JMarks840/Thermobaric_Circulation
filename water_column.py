# -*- coding: utf-8 -*-
"""
Generating the main object of the simulation: the water column. The profiles of
the different quantities are saved as arrays as its attributes.
All calculations are done on these attributes and every change is saved there
as well.

Note: The conductivity profile, and therefore the k25 profile, are considered
being constant and do not change.

Created on Tue Aug 20 11:04:32 2024

@author: marksj
"""

# %% Imports

# Third party imports
import os
import matplotlib.pyplot as plt
import numpy as np

# Local imports
import calculations as calcs
import sound_velocity as sv
import density as dens
import lake_parameters as lps

# %% Water Column Definition
# Attributes of the water column are the profiles of the corresponding quantity
# in form of numpy arrays. All data of the profiles is stored in its attributes
# and all calculations are done on them. All manipulation is done on this
# object.


class WaterColumn:
    """
    Creating a water column with all needed quantities. The attributes are
    numpy arrays of the corresponding quantities. Take care of the units of the
    quantities. The expected units are:
        pressure: [bar] or [dbar],
        temperature: [°C],
        conductivity [mS/cm],
        sound velocity: [m/s]

    Input:
        lake name (str),

        lake specific parameters (LakeParameters),

        measured lake data as numpy array (name as str; for more information
        see the staticmethod 'load_profiles'),

        pressure unit ('dbar' or 'bar' as str)

        indices of the columns of the corresponding quantities in the array
        depending on the probe used (the first column (index 0) is usually
        pressure):
            pressure, temperature, conductivity and sound velocity (optional)

        chosen in-situ density formulation (see density.py):
            (1) specific in-situ formulation (DEFAULT)
            (2) TEOS-10 formulation (only available with salinity=0) (enter
                'teos10' as string for this option)
            (3) formulation by Chen and Millero (1986) (enter 'CM' as string
                for this option)
            (4) potential density using Moreira et al. (2016); this is not
                intended to be used and the descriptions of other functions do
                not include this option explicitly, it only exists for
                demonstration purposes; if this option is used the stability
                check is based on potential density instead and in-situ density
                is simply set to 0 everywhere (enter 'potdens' as string for
                this option)

    The created object has the following attributes:
            lake name, lake parameters and the following profiles: pressure
            [bar], depth, temperature, k25, sound velocity, potential denstiy
            and in situ density

    Note:
        The temperature input is expected to be potential temperature.
        The sound velocity is only optional for the default in-situ density
        option. For the other in-situ density options the a given sound
        velocity will be ignored and directly calculated instead.
    """

    def __init__(self, lake: str, parameters: lps.LakeParameters,
                 profiles_file_name: str, pressure_unit: str, pressure_index,
                 temperature_index, conductivity_index,
                 sound_velocity_index=None,
                 is_density_type: str | None = None):
        print("Creating water column.")
        profiles = self.load_profiles(profiles_file_name, pressure_unit,
                                      pressure_index)
        self.input_profile = profiles_file_name
        self.lake = lake
        self.parameters = parameters
        if is_density_type == "teos10":
            self.is_density_type = "teos10"
        elif is_density_type == "CM":
            self.is_density_type = "CM"
        elif is_density_type == "potdens":
            self.is_density_type = "potdens"
        else:
            self.is_density_type = "default"
        self.pressure = profiles[:, pressure_index]
        self.depth = calcs.depth(self.pressure)
        self.temperature = profiles[:, temperature_index]
        self.conductivity = profiles[:, conductivity_index]
        self.k25 = calcs.k25(self.temperature, self.conductivity,
                             self.parameters.alpha)
        if self.is_density_type == "default":
            if sound_velocity_index:
                self.sound_velocity = profiles[:, sound_velocity_index]
                self.pot_density = dens.Moreira(self.temperature, self.k25,
                                                self.parameters.lambda0,
                                                self.parameters.lambda1)
                self.insitu_density = dens.rho_insitu(self.pressure,
                                                      self.temperature,
                                                      self.k25,
                                                      self.sound_velocity,
                                                      self.parameters.lambda0,
                                                      self.parameters.lambda1)
            if not sound_velocity_index:
                self.calculate_sv_dens()
        else:
            self.calculate_sv_dens()
        self.calculate_n2()
        self.conv_top_idx = np.zeros(len(self.depth), dtype=np.int32)
        self.conv_bot_idx = np.zeros(len(self.depth), dtype=np.int32)
        self.mix_idx = np.zeros(len(self.depth), dtype=np.int32)
        print("Water Column created.")
        print("Profiles of the created water column of Lake " +
              self.lake + ": pressure, depth, temperature, k25, sound" +
              "velocity, potential and in-situ density.")
        print()

    @staticmethod
    def load_profiles(file_name: str, pressure_unit: str, pressure_index):
        """
        Load data in form of a prepared numpy array. Data has to be stored in
        the folder 'input_data'. The array is suspected to have the data only
        without any headers. To identify the columns and the units, please use
        the original data file. The pressure is usually in the first column
        (start counting of the columns at the pressure with 0) and is converted
        into bar automatically if it is measured in dbar. Other units have to
        be converted to bar or mbar before loading the data.

        Input:
                file name with ending .npy (str), pressure unit ('bar' or
                'dbar' as string);
                expected units for important quantities are:
                    pressure: [bar] or [dbar],
                    temperature: [°C],
                    conductivity [mS/cm],
                    sound velocity: [m/s]
        """
        data = np.load("input_data/" + file_name)
        print("Loading successful.")
        print("Loaded file: " + file_name)
        if pressure_unit == "dbar":
            data[:, pressure_index] = data[:, pressure_index]*0.1
        return data

    @staticmethod
    def create_input(max_depth, depth_steps, temperature_value,
                     conductivity_value):
        """
        Creates an input .txt file with needed columns to create a water column
        in the input_data folder. The temperature and conductivity are set on
        one value for the whole water column, respectively.

        Input:
            maximum depth [m] (positive downward), depth resolution/steps [m],
            temperature [°C] and conductivity [mS/cm] of the whole water column

        Output:
            .txt file and .npy file for creation of a water column

        Returns:
            file name (str)
        """
        depth = np.arange(0, max_depth, depth_steps)
        pressure = depth*0.98/10  # [bar]
        temperature = np.full(len(depth), temperature_value)
        conductivity = np.full(len(depth), conductivity_value)
        combined_profiles = np.column_stack((pressure, temperature,
                                             conductivity))
        file_name = (f"depth{max_depth}steps{depth_steps}" +
                     f"T{temperature_value}C{conductivity_value}")
        os.makedirs("input_data/", exist_ok=True)
        np.savetxt("input_data/" + file_name + ".txt", combined_profiles,
                   delimiter=",")
        np.save("input_data/" + file_name, combined_profiles)
        print(f"{file_name} was created as .txt and as .npy in the folder " +
              "input_data.")
        print(f"The profiles range over a depth of {max_depth} m with values" +
              f" every {depth_steps} m. The assigned temperature is " +
              f"{temperature_value} °C and the conductivity is " +
              f"{conductivity_value} mS/cm.")
        print("To create the water column use")
        print(f"'{file_name}.npy' as 'profiles_file_name'")
        print("'bar' as 'pressure_unit'")
        print("'0' as 'pressure_index'")
        print("'1' as 'temperature_index'")
        print("'2' as 'conductivity_index'")
        print()
        return file_name

    def plot_temperature(self):
        """
        Plotting the temperature profile.

        Note:
            The approximation of the temperature of maximum density is used for
            this plot. Hence, there can occur some discrepancies.
            For the correct temperature of maximum density the
            corresponding function from calculations.py should be used.
        """
        plt.plot(self.temperature, self.depth)
        plt.plot(calcs.T_md_appr(self.pressure), self.depth,
                 color='black',
                 label="Temperature of maximum density (assuming no salinity)")
        plt.xlim(3, 5)
        plt.title("Temperature Profile")
        plt.xlabel("Temperature [°C]")
        plt.ylabel("Depth [m]")
        plt.legend()
        plt.show()

    def plot_k25(self):
        """
        Plotting the k25 profile.
        """
        plt.plot(self.k25, self.depth)
        plt.title("Conductivity Profile")
        plt.xlabel("Conductivity at 25 °C [mS/cm]")
        plt.ylabel("Depth [m]")
        plt.show()

    def plot_sound_velocity(self):
        """
        Plotting the sound velocity profile.
        """
        plt.plot(self.sound_velocity, self.depth)
        plt.xlim(1420, 1430)
        plt.title("Sound Velocity Profile")
        plt.xlabel("Sound velocity [m/s]")
        plt.ylabel("Depth [m]")
        plt.show()

    def plot_pot_density(self):
        """
        Plotting the potential density profile.
        """
        plt.plot(self.pot_density, self.depth)
        plt.title("Potential Density Profile")
        plt.xlabel(r"In-situ density [kg/m^3]")
        plt.ylabel("Depth [m]")
        plt.show()

    def plot_insitu_density(self):
        """
        Plotting the in-situ density profile.
        """
        plt.plot(self.insitu_density, self.depth)
        plt.title("In-situ Density Profile")
        plt.xlabel(r"In-situ density [kg/m^3]")
        plt.ylabel("Depth [m]")
        plt.show()

    def plot_n2(self):
        """
        Plotting the Brunt-Väisälä-Frequency profile.
        """
        plt.plot(self.n2[:-1], self.depth[:-1])
        plt.title("Brunt-Väisälä-Frequency Profile")
        plt.xlabel(r"Brunt-Väisälä-Frequency [1/s^2]")
        plt.ylabel("Depth [m]")
        plt.show()

    def plot_all(self):
        """
        Plotting the temperature, k25, sound velocity, potential density and
        in-situ density profile at once.
        """
        self.plot_temperature()
        self.plot_k25()
        self.plot_sound_velocity()
        self.plot_pot_density()
        if self.is_density_type == "potdens":
            pass
        else:
            self.plot_insitu_density()
        self.plot_n2()

    def calculate_sv_dens(self):
        """
        Calculating the sound velocity, potential and the in-situ density based
        on the assigned pressure, temperature and salinity of the water column,
        since the last two can be changed and are the controlling quantities.
        The pressure, k25 and the depth stay always the same for the arrays.
        This method should always be used after a change of the temperature at
        some point in the column to keep the other quantities up to date.
        """
        self.sound_velocity = sv.Belogolskii(self.pressure, self.temperature)
        self.pot_density = dens.Moreira(self.temperature, self.k25,
                                        self.parameters.lambda0,
                                        self.parameters.lambda1)
        if self.is_density_type == "default":
            self.insitu_density = dens.rho_insitu(self.pressure,
                                                  self.temperature,
                                                  self.k25,
                                                  self.sound_velocity,
                                                  self.parameters.lambda0,
                                                  self.parameters.lambda1)
        elif self.is_density_type == "teos10":
            self.insitu_density = dens.rho_teos10(self.pressure,
                                                  self.temperature)
        elif self.is_density_type == "CM":
            self.insitu_density = dens.rho_CM(self.pressure, self.temperature,
                                              self.k25, self.parameters.alpha)
        elif self.is_density_type == "potdens":
            self.insitu_density = np.zeros(len(self.depth))

    def calculate_n2(self):
        """
        Calculating the Brunt-Väisälä-Frequency.

        Note:
            If the in-situ density formulation of Chen and Millero (1986) is
            used, the Brunt-Väisälä-Frequency is still calculated with the
            default in-situ density formulation, which can lead to
            discrepancies between the stabilisation and the calculated
            Brunt-Väisälä-Frequency.
        """
        if self.is_density_type == "teos10":
            self.n2 = calcs.N2_teos10(self.depth, self.pressure,
                                      self.temperature)
        elif self.is_density_type == "potdens":
            self.n2 = calcs.N2_pot(self.depth, self.temperature, self.k25,
                                   self.parameters.lambda0,
                                   self.parameters.lambda1)
        else:
            self.n2 = calcs.N2_is(self.depth, self.pressure, self.temperature,
                                  self.k25, self.sound_velocity,
                                  self.parameters.lambda0,
                                  self.parameters.lambda1)

    def create_grid(self, grid_size):
        """
        Creating a grid (an array) of a given length for the simulation, which
        means to convert the water column into an array with a certain length:
        The original input array is sliced into a given number of pieces of the
        same size and the last data point of each piece is used in the new
        grid (array). Only the uppermost water parcel is kept as the surface
        water parcel and only represents itself.
        (the indices are not rounded but the decimal places are simply omitted)

        Note: All values except pressure, temperature and k25 are
        recalculated. The inital water column will be overwritten with the new
        grid format and data.

        Input:
            grid size

        Output:
            water column with new grid size as array length
        """
        indices = np.linspace(0, len(self.depth)-1, grid_size)
        indices = indices.astype(np.int64)
        self.pressure = self.pressure[indices]
        self.depth = self.depth[indices]
        self.temperature = self.temperature[indices]
        self.conductivity = self.conductivity[indices]
        self.k25 = self.k25[indices]
        self.sound_velocity = self.sound_velocity[indices]
        self.pot_density = self.pot_density[indices]
        self.insitu_density = self.insitu_density[indices]
        self.n2 = self.n2[indices]
        self.conv_top_idx = self.conv_top_idx[indices]
        self.conv_bot_idx = self.conv_bot_idx[indices]
        self.mix_idx = self.mix_idx[indices]

    def set_surface_temperature(self, surface_temperature):
        """
        Control the surface temperature and set it to a certain level.

        Input:
            surface temperature [°C]
        """
        self.temperature[0] = surface_temperature
        self.calculate_sv_dens()

    def set_bottom_temperature(self, bottom_temperature):
        """
        Control the bottom temperature and set it to a certain level.

        Input:
            bottom temperature [°C]
        """
        self.temperature[-1] = bottom_temperature
        self.calculate_sv_dens()

    def set_k25(self, value):
        """
        Set the whole k25 profile to one vlaue.

        Input:
            k25 value [mS/cm]

        Output:
            k25 profile with the input value over the whole profile
        """
        self.k25 = np.full(self.k25.shape, value)
        self.conductivity = calcs.C_from_k25(self.temperature,
                                             self.k25,
                                             self.parameters.alpha)
        self.calculate_sv_dens()

    def diffuse(self, exchange_volume):
        """
        Diffuses the water column recalculating the temperature and k25 of each
        water parcel of the water column. For this, the value of each water
        parcel is averaged using half of the exchange volume of the water
        parcels above and below and the rest of the old values of the water
        parcel itself.
        The surface values are not changed and the bottom water parcel is
        averaged over half of the exchange volume of the water parcel above
        and the rest of the old value of the water parcel itself.
        Then the sound velocity and the density are calculated based on the new
        temperature and k25 profiles.

        Input:
            exchange volume (0 < exchange volume < 1)
        """
        index = np.linspace(1, (len(self.depth) - 2), (len(self.depth) - 2)
                            ).astype("int64")
        new_temperature = np.zeros(len(self.temperature))
        new_k25 = np.zeros(len(self.k25))
        for i in index:
            new_temperature[i] = ((exchange_volume/2)*self.temperature[i-1] +
                                  (1 - exchange_volume)*self.temperature[i] +
                                  (exchange_volume/2)*self.temperature[i+1])
            new_k25[i] = ((exchange_volume/2)*self.k25[i-1] +
                          (1 - exchange_volume)*self.k25[i] +
                          (exchange_volume/2)*self.k25[i+1])
        new_temperature[-1] = ((exchange_volume/2)*self.temperature[-2] +
                               (1 - (exchange_volume/2))*self.temperature[-1])
        new_k25[0] = ((exchange_volume/2)*self.k25[1] +
                      (1 - (exchange_volume/2))*self.k25[0])
        new_k25[-1] = ((exchange_volume/2)*self.k25[-2] +
                       (1 - (exchange_volume/2))*self.k25[-1])
        self.temperature[1:] = new_temperature[1:]
        self.k25 = new_k25
        self.conductivity = calcs.C_from_k25(self.temperature,
                                             self.k25,
                                             self.parameters.alpha)
        self.calculate_sv_dens()
        self.calculate_n2()

    def convection_top(self):
        """
        Checking the stability of the water column top down based on the
        in-situ density:
            Comparison of in-situ density at pressure of comparison of surface
            temperature, k25 and sound velocity and in-situ temperature, k25
            and sound velocity.
        This stops at the first stable layer.
        Then the column is stabilised by setting the temperature equal to the
        surface temperature over the depth of the instable layers.
        The indices, where instabilities occured, are saved for this execution
        and will be overwritten by the next execution of this method.

        Note:
            The k25 value is averaged for the instable layers already while
            checking the stability. This is done because mixing without further
            input is assumed for the k25, while the temperature is held
            constant at surface temperature due to surface heat control.
            The Brunt-Väisälä-Frequency is not automatically calculated in this
            function, but in the complete stabilisation function stabilise().
        """
        self.conv_top_idx[:] = 0
        T_surface = self.temperature[0]
        k25_average = self.k25[0]
        index = 1
        if self.is_density_type == "default":
            while (dens.rho_insitu(self.pressure[index], T_surface,
                                   k25_average,
                                   sv.Belogolskii(self.pressure[index],
                                                  T_surface),
                                   self.parameters.lambda0,
                                   self.parameters.lambda1
                                   )
                   > (self.insitu_density[index] + 1e-9)):
                k25_average = np.mean(self.k25[:index])
                index += 1
                if index > (len(self.depth) - 1):
                    break
        elif self.is_density_type == "teos10":
            while (dens.rho_teos10(self.pressure[index], T_surface)
                   > (self.insitu_density[index] + 1e-9)):
                k25_average = np.mean(self.k25[:index])
                index += 1
                if index > (len(self.depth) - 1):
                    break
        elif self.is_density_type == "CM":
            while (dens.rho_CM(self.pressure[index], T_surface, k25_average,
                               self.parameters.alpha)
                   > (self.insitu_density[index] + 1e-9)):
                k25_average = np.mean(self.k25[:index])
                index += 1
                if index > (len(self.depth) - 1):
                    break
        elif self.is_density_type == "potdens":
            while (dens.Moreira(T_surface, k25_average,
                                self.parameters.lambda0,
                                self.parameters.lambda1)
                   > self.pot_density[index]):
                k25_average = np.mean(self.k25[:index])
                index += 1
                if index > (len(self.depth) - 1):
                    break
        self.conv_top_idx[:index] = 1
        self.temperature[:index] = T_surface
        self.k25[:index] = k25_average
        self.conductivity = calcs.C_from_k25(self.temperature,
                                             self.k25,
                                             self.parameters.alpha)
        self.calculate_sv_dens()

    def convection_bottom(self):
        """
        Checking the stability of the water column by scanning through all
        layers bottom up and checking their stability compared to the layers
        below based on the in-situ density:
            Comparison of in-situ density at pressure of comparison of 'chosen
            layer' temperature, k25 and sound velocity and 'next lower layer'
            temperature, k25 and sound velocity.
            Then averaging the temperature for instable layers and using this
            average as the new 'chosen layer' and comparing this with the next
            'next lower layer' as new 'next lower layer'. The pressure of
            comparison is always at the 'next lower layer'.
            This is done repeatedly by scanning trhough the water column bottom
            up to select the 'chosen layer', excluding the uppermost surface
            layer, because it is temperature controlled, and the lowermost
            bottom layer.
        The indices, where instabilities occured, are saved for this execution
        and will be overwritten by the next execution of this method.

        Note:
            The k25 value is averaged for the instable layers similar to the
            temperature.
            The Brunt-Väisälä-Frequency is not automatically calculated in this
            function, but in the complete stabilisation function stabilise().
        """
        self.conv_bot_idx[:] = 0
        index = np.linspace((len(self.depth) - 2), 1, (len(self.depth) - 2)
                            ).astype("int64")
        for i in index:
            T_average = self.temperature[i]
            k25_average = self.k25[i]
            j = 1
            if self.is_density_type == "default":
                while (dens.rho_insitu(self.pressure[i+j], T_average,
                                       k25_average,
                                       sv.Belogolskii(self.pressure[i+j],
                                                      T_average),
                                       self.parameters.lambda0,
                                       self.parameters.lambda1)
                       > (self.insitu_density[i+j] + 1e-9)):
                    T_average = np.mean(self.temperature[i:(i+j+1)])
                    k25_average = np.mean(self.k25[i:(i+j+1)])
                    j += 1
                    if (i+j) > (len(self.depth) - 1):
                        break
            elif self.is_density_type == "teos10":
                while (dens.rho_teos10(self.pressure[i+j], T_average)
                       > (self.insitu_density[i+j] + 1e-9)):
                    T_average = np.mean(self.temperature[i:(i+j+1)])
                    k25_average = np.mean(self.k25[i:(i+j+1)])
                    j += 1
                    if (i+j) > (len(self.depth) - 1):
                        break
            elif self.is_density_type == "CM":
                while (dens.rho_CM(self.pressure[i+j], T_average, k25_average,
                                   self.parameters.alpha)
                       > (self.insitu_density[i+j] + 1e-9)):
                    T_average = np.mean(self.temperature[i:(i+j+1)])
                    k25_average = np.mean(self.k25[i:(i+j+1)])
                    j += 1
                    if (i+j) > (len(self.depth) - 1):
                        break
            elif self.is_density_type == "potdens":
                while (dens.Moreira(T_average, k25_average,
                                    self.parameters.lambda0,
                                    self.parameters.lambda1)
                       > self.pot_density[i+j]):
                    T_average = np.mean(self.temperature[i:(i+j+1)])
                    k25_average = np.mean(self.k25[i:(i+j+1)])
                    j += 1
                    if (i+j) > (len(self.depth) - 1):
                        break
            if j > 1:
                self.conv_bot_idx[i:(i+j)] = 1
            self.temperature[i:(i+j)] = T_average
            self.k25[i:(i+j)] = k25_average
            self.conductivity = calcs.C_from_k25(self.temperature,
                                                 self.k25,
                                                 self.parameters.alpha)
            self.calculate_sv_dens()

    def stabilise(self):
        """
        Running convection_top() and convection_bottom() successively.
        """
        self.convection_top()
        self.convection_bottom()
        self.mix_idx = self.conv_top_idx | self.conv_bot_idx
        self.calculate_n2()

    def stepwise_simulation(self, surface_temperature_change, exchange_volume):
        """
        Changing the surface temperature by a given amount and running a
        simulation (diffusion and stabilisation) for this change.

        Input:
            surface temperature change

        Output:
            stabilised water column after the change
        """
        self.set_surface_temperature(self.temperature[0] +
                                     surface_temperature_change)
        self.diffuse(exchange_volume)
        self.stabilise()
        self.calculate_n2()
        self.plot_temperature()
