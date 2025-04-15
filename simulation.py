# -*- coding: utf-8 -*-
"""
Run a simulation given a water column as starting object and as object to
conduct the simulation on with the surface temperature as driver of the
simulation but without altering the conductivity profile if given.

@author: Joshua Marks
"""

# %% Imports

# Third party imports
import time
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap
from matplotlib.colors import SymLogNorm

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
        self.t_logs = []
        self.mix_log = []
        self.n2_log = []

    def run(self):
        """
        Given the water column as starting values and as objects to perform the
        simulation on, the surface temperature as input variable and the grid
        size and time step size as simulation parameters this runs the
        simulation for the given intput.

        Output:
            .txt files with the depth, temperature, k25, sound veloctiy and in-
            situ density data for each hour (for smaller time steps the last of
            the corresponding hour and for larger time steps the corresponding
            days are saved). They are named after the lake, the day and end
            with '_sim.txt'.
            Plotting the input data and the final output data from the last
            time step in the given grid format.
            Saving the temperature over the time period for the uppermost not
            controlled layer, the middle layer and the lowermost layer in the
            file 't_logs.txt'.
            Saving a data frame with the indication of occurred mixing as
            'mix_log.txt'.
            Saving a data frame with the Brunt-Väisälä-Frequency [1/s^2] as
            'n2_log.txt'.

        Note:
            The folders 'output_data' has to exist.
        """
        print("Running simulation")
        os.makedirs("output_data/" + self.water_column.lake, exist_ok=True)
        print(f"Folder {self.water_column.lake} in output_data created.")
        print("The output files can be found there.")
        print("The input profiles are shown in the gridded format first.")
        start_time = time.time()
        print("Simulation started at", time.strftime("%H:%M:%S",
                                                     time.localtime()))
        # Simulation set-up
        self.water_column.create_grid(self.grid_size)
        middle = round(0.5*len(self.water_column.depth))
        self.water_column.plot_all()
        self.surface_temp_ts.create_time_series(self.time_step_size)
        self.surface_temp_ts.plot_surface_temperature()
        # Simulation
        for data_point in self.surface_temp_ts.data.itertuples():
            self.water_column.set_surface_temperature(data_point.Temperature)
            self.water_column.diffuse(self.exchange_volume)
            self.water_column.stabilise()
            self.water_column.calculate_n2()
            # Saving current profiles:
            current_data = np.column_stack((self.water_column.depth,
                                            self.water_column.pressure,
                                            self.water_column.temperature,
                                            self.water_column.k25,
                                            self.water_column.sound_velocity,
                                            self.water_column.insitu_density))
            header = ("Depth, Pressure,Temperature,k25,Sound Velocity," +
                      "In-situ Density")
            np.savetxt("output_data/" + self.water_column.lake + "/"
                       + self.water_column.lake + "_" +
                       data_point.Index.strftime("%Y-%m-%d_%H") + "_sim.txt",
                       current_data, delimiter=",", header=header,
                       comments="")
            # Saving the current T-logs
            t_log_entry = pd.DataFrame({"Time": data_point.Index,
                                        "Temperature top":
                                        [self.water_column.temperature[1]],
                                        "Temperature middle":
                                        [self.water_column.temperature[middle]
                                         ],
                                        "Temperature bottom":
                                        [self.water_column.temperature[-1]]})
            self.t_logs.append(t_log_entry)
            # Saving the current mixing log
            mix_log_entry = pd.DataFrame({"Time": data_point.Index,
                                          "Mixing log":
                                          [self.water_column.mix_idx]})
            self.mix_log.append(mix_log_entry)
            # Saving the current stability profile
            n2_log_entry = pd.DataFrame({"Time": data_point.Index,
                                         "N2": [self.water_column.n2]})
            self.n2_log.append(n2_log_entry)
        # Saving t_logs
        self.t_logs = pd.concat(self.t_logs, ignore_index=True)
        self.t_logs.set_index("Time", inplace=True)
        self.t_logs.to_csv("output_data/" + self.water_column.lake + "/"
                           + self.water_column.lake + "_t_logs.txt")
        # Saving mix_log
        self.mix_log = pd.concat(self.mix_log, ignore_index=True)
        self.mix_log.set_index("Time", inplace=True)
        self.mix_log.to_csv("output_data/" + self.water_column.lake + "/"
                            + self.water_column.lake + "_mix_log.txt")
        # Saving n2_log
        self.n2_log = pd.concat(self.n2_log, ignore_index=True)
        self.n2_log.set_index("Time", inplace=True)
        self.n2_log.to_csv("output_data/" + self.water_column.lake + "/"
                           + self.water_column.lake + "_n2_log.txt")
        # Last output
        self.plot_t_logs()
        self.plot_mix_log()
        self.plot_n2_log()
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

    def plot_t_logs(self):
        """
        Plotting the temperature logs.
        """
        self.t_logs.plot(title="Temperature-Logs", xlabel="Time",
                         ylabel="Temperature [°C]", legend=True)
        plt.show()

    def plot_mix_log(self):
        """
        Plotting the mixing log as 2D colored image with two colors indicating
        mixing events and no-mixing events.
        """
        data = np.array(self.mix_log["Mixing log"].tolist()).transpose()
        time = np.array(self.mix_log.index)
        cmap = ListedColormap(["lightsteelblue", "steelblue"])
        im = plt.imshow(data, cmap=cmap, aspect="auto", interpolation="none",
                        extent=[0, len(time), (self.water_column.depth[-1] +
                                               self.water_column.depth[1]),
                                self.water_column.depth[0]])
        plt.xticks(rotation=20)
        legend_elements = [mpatches.Patch(color=im.cmap(im.norm(0)),
                                          label="no mixing"),
                           mpatches.Patch(color=im.cmap(im.norm(1)),
                                          label="mixing")]
        plt.legend(handles=legend_elements, bbox_to_anchor=(1.01, 1), loc=2,
                   borderaxespad=0)
        plt.title("Mixing")
        plt.xlabel("Time Step")
        plt.ylabel("Depth [m]")
        plt.show()

    def plot_n2_log(self):
        """
        Plotting the N^2 (Brunt-Väisälä-Frequency) log as 2D colored image with
        positive values colored as shown by the colorbar and all negative
        values, regardless of their absolute value, are shown in white.
        """
        data = np.array(self.n2_log["N2"].tolist()).transpose()
        time = np.array(self.n2_log.index)
        cmap = "plasma"
        im = plt.imshow(data, cmap=cmap, norm=SymLogNorm(linthresh=1e-10,
                                                         vmin=0),
                        aspect="auto", interpolation="none",
                        extent=[0, len(time), (self.water_column.depth[-2] +
                                               self.water_column.depth[1]),
                                self.water_column.depth[0]])
        im.cmap.set_under("white")
        plt.colorbar(im).set_label(r"N$^2$ [1/$\text{s}^2$]")
        plt.xticks(rotation=20)
        plt.title("Stability")
        plt.xlabel("Time Step")
        plt.ylabel("Depth [m]")
        plt.show()
