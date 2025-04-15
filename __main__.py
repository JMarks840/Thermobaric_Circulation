# -*- coding: utf-8 -*-
"""
File for creating input and simulation and to run the simulation.

@author: Joshua Marks
"""

# %% Imports

# Local imports
import water_column as wc
import lake_parameters as lps
import surface_temperature as st
import simulation as sim

# %% Example Simulation

"""
Example simulation using Lake Shikotsu as motivation.
The input is created as an isothermal temperature profile without any salinity.
The surface temperature is from Lake Shikotsu and the file is named
'Shikotsu_surfacetemperature.txt', but it will be extrapolated and repeated
over 7 years.
"""

# Setting parameters
# for manual input profile:
depth = 360  # [m]
depth_step_size = 2  # [m]
T = 4.2  # [°C] (whole profile)
C = 0  # [mS/cm] (whole profile)

# for surface temperature:
data_frequency = "min"
repetitions_years = 7

# for simulation:
grid_size = int(depth/depth_step_size)
exchange_volume = 0.5
time_step_size = "hourly"

# Create manual input profile
input_file_name = wc.WaterColumn.create_input(depth, depth_step_size, T, C)

# Create water column
water_column = wc.WaterColumn("Shikotsu_example",
                              lps.SHIKOTSU_PARAMETERS,
                              input_file_name + ".npy", "bar", 0, 1, 2)

# Create surface temperature
surf_temp = st.SurfaceTemperature("Shikotsu", "231590",
                                  "231590_Shikotsu_surfacetemperature.txt",
                                  "Winter2023/2024")
# prolong and repeat surface temperature input
surf_temp.repeat_time_period_annually(data_frequency, repetitions_years)
# applying moving average
surf_temp.moving_average(time_step_size)

# Creating the simulation
simulation_Shikotsu_manualProfiles = sim.Simulation(water_column, grid_size,
                                                    exchange_volume, surf_temp,
                                                    time_step_size)

# Run simulation
simulation_Shikotsu_manualProfiles.run()
