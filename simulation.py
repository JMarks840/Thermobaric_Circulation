# -*- coding: utf-8 -*-
"""
Run a simulation given a water column as starting object and as object to
conduct the simulation on and the surface temperature as driver of the
simulation.

Created on Tue Sep 10 16:16:31 2024

@author: marksj
"""

# %% Imports

# Third party imports
import time
import os
import numpy as np
import xarray as xr

# Local imports
import water_column as wc
import surface_temperature as surft

# %% Simulation


class Simulation:
    """
    Creating a simple 1D simulation by setting the starting and framework
    conditions. The actual simulation can be operated on this object using the
    included water column and surface temperature series.

    Input:
        water column (WaterColumn),

        grid size (as in create_grid() in WaterColumn),

        surface temperature time series (SurfaceTemperature),

        time step size (as in create_time_series() in SurfaceTemperature)

    Note:
        The used water column and surface temperature will be altered during
        the simulation. For restoring the original objects, please recreate and
        reload the resprective objects manually.
    """

    def __init__(self, water_column: wc.WaterColumn, grid_size,
                 exchange_volume,
                 surface_temperature_time_series: surft.SurfaceTemperature,
                 time_step_size: str):
        self.water_column = water_column
        self.grid_size = grid_size
        self.exchange_volume = exchange_volume
        self.surface_temp_ts = surface_temperature_time_series
        self.time_step_size = time_step_size

    def run(self):
        """
        Given the water column as starting values and as objects to perform the
        simulation on, the surface temperature as input variable and the grid
        size, the exchange volume and time step size as simulation parameters
        this runs the simulation for the given intput.

        Output:
            netCDF (*.nc) file containing temperature, pressure, k25, sound
            velocity, potential and in situ density, n2 and mixing log for the
            whole simulation time in hourly timesteps (for smaller time steps
            the last of the corresponding hour and for larger time steps the
            corresponding days are saved) for every depth. The file is named
            after the lake name in the water column input file.
            Plotting the input data and the final output data from the last
            time step in the given grid format.

        Note:
            The surface temperature input can be manipulated by surface
            temperature operations (see surface_temperature.py) before creating
            the simulation with this input.
        """
        print("Running simulation")
        os.makedirs("output/", exist_ok=True)
        print("The input profiles are shown in the gridded format first.")
        start_time = time.time()
        print("Simulation started at", time.strftime("%H:%M:%S",
                                                     time.localtime()))
        # Simulation set-up
        self.water_column.create_grid(self.grid_size)
        self.water_column.plot_all()
        self.surface_temp_ts.create_time_series(self.time_step_size)
        self.surface_temp_ts.plot_surface_temperature()
        # Initializing lists for output
        times = []
        depth = self.water_column.depth
        pressure_profiles = []
        temperature_profiles_diff = []
        temperature_profiles = []
        conductivity_profiles = []
        k25_profiles_diff = []
        k25_profiles = []
        sound_velocity_profiles = []
        pot_density_profiles_diff = []
        pot_density_profiles = []
        is_density_profiles_diff = []
        is_density_profiles = []
        n2_profiles_diff = []
        n2_profiles = []
        mix_log_profiles = []
        # Simulation
        for data_point in self.surface_temp_ts.data.itertuples():
            self.water_column.set_surface_temperature(data_point.Temperature)
            self.water_column.diffuse(self.exchange_volume)
            # Saving profiles after diffusion and before stabilization
            temperature_profiles_diff.append(self.water_column.temperature.copy())
            k25_profiles_diff.append(self.water_column.k25.copy())
            pot_density_profiles_diff.append(self.water_column.pot_density.copy())
            is_density_profiles_diff.append(self.water_column.insitu_density.copy())
            n2_profiles_diff.append(self.water_column.n2.copy())
            # Stabilizing water column
            self.water_column.stabilise()
            # Saving current profiles after stabilization
            times.append(data_point.Index)
            pressure_profiles.append(self.water_column.pressure.copy())
            temperature_profiles.append(self.water_column.temperature.copy())
            conductivity_profiles.append(self.water_column.conductivity.copy())
            k25_profiles.append(self.water_column.k25.copy())
            sound_velocity_profiles.append(self.water_column.sound_velocity.copy())
            pot_density_profiles.append(self.water_column.pot_density.copy())
            is_density_profiles.append(self.water_column.insitu_density.copy())
            n2_profiles.append(self.water_column.n2.copy())
            mix_log_profiles.append(self.water_column.mix_idx.copy())
        # Saving output
        times = np.array(times)
        output = xr.Dataset(
            {"pressure": xr.DataArray(
                np.vstack(pressure_profiles),
                dims=["time", "depth"],
                attrs={"long_name": "pressure", "unit": "bar"}
                ),
             "temperature": xr.DataArray(
                np.vstack(temperature_profiles),
                dims=["time", "depth"],
                attrs={"long_name": "water temperature", "unit": "degC"}
                ),
             "temperature_diff": xr.DataArray(
                 np.vstack(temperature_profiles_diff),
                 dims=["time", "depth"],
                 attrs={"long_name": "water temperature after diffusion and before stabilization",
                        "unit": "degC"}
                 ),
             "conductivity": xr.DataArray(
                 np.vstack(conductivity_profiles),
                 dims=["time", "depth"],
                 attrs={"long_name": "conductivity", "unit": "mS/cm"}
                 ),
             "k25": xr.DataArray(
                 np.vstack(k25_profiles),
                 dims=["time", "depth"],
                 attrs={"long_name": "conductivity at 25 °C", "unit": "mS/cm"}
                 ),
             "k25_diff": xr.DataArray(
                 np.vstack(k25_profiles_diff),
                 dims=["time", "depth"],
                 attrs={"long_name": "conductivity at 25 °C after diffusion and before stabilization",
                        "unit": "mS/cm"}
                 ),
             "sound_velocity": xr.DataArray(
                 np.vstack(sound_velocity_profiles),
                 dims=["time", "depth"],
                 attrs={"long_name": "sound velocity", "unit": "m/s"}
                 ),
             "pot_density": xr.DataArray(
                 np.vstack(pot_density_profiles),
                 dims=["time", "depth"],
                 attrs={"long_name": "potential density", "unit": "kg/m^3"}
                 ),
             "pot_density_diff": xr.DataArray(
                 np.vstack(pot_density_profiles_diff),
                 dims=["time", "depth"],
                 attrs={"long_name": "potential density after diffusion and before stabilization",
                        "unit": "kg/m^3"}
                 ),
             "is_density": xr.DataArray(
                 np.vstack(is_density_profiles),
                 dims=["time", "depth"],
                 attrs={"long_name": "in situ density", "unit": "kg/m^3"}
                 ),
             "is_density_diff": xr.DataArray(
                 np.vstack(is_density_profiles_diff),
                 dims=["time", "depth"],
                 attrs={"long_name": "in situ density after diffusion and before stabilization",
                        "unit": "kg/m^3"}
                 ),
             "stability": xr.DataArray(
                 np.vstack(n2_profiles),
                 dims=["time", "depth"],
                 attrs={"long_name": "stability frequency", "unit": "1/s^2"}
                 ),
             "stability_diff": xr.DataArray(
                 np.vstack(n2_profiles),
                 dims=["time", "depth"],
                 attrs={"long_name": "stability frequencyafter diffusion and before stabilization",
                        "unit": "1/s^2"}
                 ),
             "mix_log": xr.DataArray(
                 np.vstack(mix_log_profiles),
                 dims=["time", "depth"],
                 attrs={"long_name": "mixing log of stabilization process",
                        "unit": "1 for mixed, 0 for not mixed"}
                 )
             },
            coords={
                "time": times,
                "depth": xr.DataArray(
                    depth,
                    dims=["depth"],
                    attrs={"long_name": "depth below water surface",
                           "unit": "m"}
                    )
                },
            attrs={
                "title": "model data output",
                "lake_name": self.water_column.lake,
                "input_profile": self.water_column.input_profile,
                "input_surface_temperature": self.surface_temp_ts.file_name,
                "is_density_type": self.water_column.is_density_type,
                "grid_size": str(self.grid_size),
                "exchange_volume": str(self.exchange_volume),
                "time_step": self.time_step_size
                }
            )
        output.to_netcdf(f"output/{self.water_column.lake}.nc")
        # Last output
        self.water_column.plot_all()
        end_time = time.time()
        runtime = end_time - start_time
        struct_runtime = time.gmtime(runtime)
        print("Simulation successful.")
        print("Simulation ended at", time.strftime("%H:%M:%S",
                                                   time.localtime()))
        print(f"The simulation took {struct_runtime.tm_hour} hours, " +
              f"{struct_runtime.tm_min} minutes and {struct_runtime.tm_sec} " +
              "seconds.")
        print("Cached data can be found in the folder 'output_data'.")
        print("The final profiles are plotted last, after the stability log.")
        print("")
