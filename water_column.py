# -*- coding: utf-8 -*-
"""
Generating the main object of the simulation: the water column. The profiles of
the different quantities are saved as arrays as its attributes.
All calculations are done on these attributes and every change is saved there
as well.

Note: The conductivity profile, and therefore the k25 profile, are considered
being constant and do not change.

@author: Joshua Marks
"""

# %% Imports

# Third party imports
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

    The created object has the following attributes:
            lake name, lake parameters and the following profiles: pressure,
            depth, temperature, k25, sound velocity and in situ density
    """

    def __init__(self, lake: str, parameters: lps.LakeParameters,
                 profiles_file_name: str, pressure_unit: str, pressure_index,
                 temperature_index, conductivity_index,
                 sound_velocity_index=None):
        print("Creating water column.")
        profiles = self.load_profiles(profiles_file_name, pressure_unit,
                                      pressure_index)
        self.lake = lake
        self.parameters = parameters
        self.pressure = profiles[:, pressure_index]
        self.depth = calcs.depth(self.pressure)
        self.temperature = profiles[:, temperature_index]
        conductivity = profiles[:, conductivity_index]
        self.k25 = calcs.k25(self.temperature, conductivity,
                             self.parameters.alpha)
        if sound_velocity_index:
            self.sound_velocity = profiles[:, sound_velocity_index]
            self.insitu_density = dens.rho_insitu(self.pressure,
                                                  self.temperature, self.k25,
                                                  self.sound_velocity,
                                                  self.parameters.lambda0,
                                                  self.parameters.lambda1)
            print("Water Column created.")
            print("Profiles of the created water column of Lake " + self.lake +
                  ": pressure, depth, temperature, k25, " +
                  "sound velocity, and in-situ density.")
        if not sound_velocity_index:
            self.calculate_sv_dens()
            print("Water Column created.")
            print("Profiles of the created water column of Lake " + self.lake +
                  ": pressure, depth, temperature and k25. " +
                  "The sound velocity and in-situ density were calculated.")
        print()
        self.n2 = None
        self.conv_top_idx = np.zeros(len(self.depth), dtype=np.int32)
        self.conv_bot_idx = np.zeros(len(self.depth), dtype=np.int32)
        self.mix_idx = np.zeros(len(self.depth), dtype=np.int32)

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
        np.savetxt("input_data/" + file_name + ".txt", combined_profiles,
                   delimiter=',')
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
        """
        plt.plot(self.temperature, self.depth)
        plt.plot(calcs.T_md_CM(self.pressure, 0), self.depth,
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
        plt.plot(self.n2, self.depth[:-1])
        plt.title("Brunt-Väisälä-Frequency Profile")
        plt.xlabel(r"Brunt-Väisälä-Frequency [1/s^2]")
        plt.ylabel("Depth [m]")
        plt.show()

    def plot_all(self):
        """
        Plotting the temperature, k25, sound velocity and in-situ density
        profile at once.
        """
        self.plot_temperature()
        self.plot_k25()
        self.plot_sound_velocity()
        self.plot_insitu_density()
        if self.n2 is not None:
            self.plot_n2()

    def calculate_sv_dens(self):
        """
        Calculating the sound velocity and the in-situ density, as well as the
        Brunt-Väisälä-Frequency, based on the assigned pressure, temperature
        and salinity of the water column, since the last two can be changed and
        are the controlling quantities.
        The pressure, k25 and the depth stay always the same for the arrays.
        This method should always be used after a change of the temperature at
        some point in the column to keep the other quantities up to date.
        """
        self.sound_velocity = sv.Belogolskii(self.pressure, self.temperature)
        self.insitu_density = dens.rho_insitu(self.pressure, self.temperature,
                                              self.k25, self.sound_velocity,
                                              self.parameters.lambda0,
                                              self.parameters.lambda1)

    def calculate_n2(self):
        self.n2 = calcs.N2_is(self.depth, self.pressure, self.temperature,
                              self.k25, self.sound_velocity,
                              self.parameters.lambda0, self.parameters.lambda1)

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
        self.k25 = self.k25[indices]
        self.sound_velocity = self.sound_velocity[indices]
        self.insitu_density = self.insitu_density[indices]
        self.conv_top_idx = self.conv_top_idx[indices]
        self.conv_bot_idx = self.conv_bot_idx[indices]
        self.mix_idx = self.mix_idx[indices]

    def set_surface_temperature(self, surface_temperature):
        """
        Control the surface temperature and set it to a certain value.

        Input:
            surface temperature [°C]
        """
        self.temperature[0] = surface_temperature
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
        self.calculate_sv_dens()

    def diffuse(self, exchange_volume):
        """
        Diffuses the water column recalculating the temperature of each water
        parcel of the water column. For this, the values are averaged using
        half of the exchange volume of the water parcels above and below and
        the rest of the old value of the water parcel itself and.
        The surface values are not changed and the bottom water parcel is
        averaged over half of the exchange volume of the water parcels above
        and the rest of the old values of the water parcels itself.
        Then the sound velocity and the density are calculated based on the new
        temperature profile.

        Input:
            exchange volume (0 < exchange volume < 1)
        """
        index = np.linspace(1, (len(self.depth) - 2), (len(self.depth) - 2)
                            ).astype("int64")
        new_temperature = np.zeros(len(self.temperature))
        for i in index:
            new_temperature[i] = ((exchange_volume/2)*self.temperature[i-1] +
                                  (1 - exchange_volume)*self.temperature[i] +
                                  (exchange_volume/2)*self.temperature[i+1])
        new_temperature[-1] = ((exchange_volume/2)*self.temperature[-2] +
                               (1 - (exchange_volume/2))*self.temperature[-1])
        self.temperature[1:] = new_temperature[1:]
        self.calculate_sv_dens()

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
        """
        self.conv_top_idx[:] = 0
        T_surface = self.temperature[0]
        index = 1
        while (dens.rho_insitu(self.pressure[index], T_surface,
                               self.k25[index],
                               sv.Belogolskii(self.pressure[index], T_surface),
                               self.parameters.lambda0, self.parameters.lambda1
                               )
               > self.insitu_density[index]):
            index += 1
            if index > (len(self.depth) - 1):
                break
        self.conv_top_idx[:index] = 1
        self.temperature[:index] = T_surface
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
        """
        self.conv_bot_idx[:] = 0
        index = np.linspace((len(self.depth) - 2), 1, (len(self.depth) - 2)
                            ).astype("int64")
        for i in index:
            T_average = self.temperature[i]
            j = 1
            while (dens.rho_insitu(self.pressure[i+j], T_average,
                                   self.k25[i],
                                   sv.Belogolskii(self.pressure[i+j],
                                                  T_average),
                                   self.parameters.lambda0,
                                   self.parameters.lambda1)
                   > self.insitu_density[i+j]):
                T_average = np.mean(self.temperature[i:(i+j+1)])
                j += 1
                if (i+j) > (len(self.depth) - 1):
                    break
            if j > 1:
                self.conv_bot_idx[i:(i+j)] = 1
            self.temperature[i:(i+j)] = T_average
            self.calculate_sv_dens()

    def stabilise(self):
        self.convection_top()
        self.convection_bottom()
        self.mix_idx = self.conv_top_idx | self.conv_bot_idx

    def stepwise_simulation(self, surface_temperature_change):
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
        self.diffuse()
        self.stabilise()
        self.plot_temperature()
