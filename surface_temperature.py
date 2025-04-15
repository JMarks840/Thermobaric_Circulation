# -*- coding: utf-8 -*-
"""
Generating the time series of the surface temperature based on the measured
surface temperature with temperature loggers over a certain period of time.

@author: Joshua Marks
"""

# %% Imports

# Third party imports
import pandas as pd
import matplotlib.pyplot as plt

# %% Surface Temperature Time Series Definition


class SurfaceTemperature:
    """
    Creating a surface temperature input object.

    Input:
        sensor number (str), file name with .txt (str), lake name (str) and the
        time period as e.g. 'Winter2023/2024' or 'Summer2024' (str)
    """

    def __init__(self, lake: str, sensor: str, file_name: str, time_period: str
                 ):
        print("Creating surface temperature series.")
        self.lake = lake
        self.sensor = sensor
        self.file_name = file_name
        self.time_period = time_period
        self.data = self.load_file(file_name)
        print("Surface temperature series created.")
        print()

    @staticmethod
    def load_file(file_name):
        """
        Loading ".txt" files, located in the folder input_data.

        Input:
            file name with ending .txt

        Note:
            The necessary columns (time and temperature) should be called
            'Time' and 'Temperature'.
        """
        file = pd.read_csv("input_data/" + file_name,
                           index_col=0, parse_dates=True)
        print("Loading successful.")
        print("Loaded file: " + file_name)
        if "Temperature" not in file.columns:
            print("Note:")
            print("Even though the surface temperature data was loaded and " +
                  "the surface temperature series will be created, the " +
                  "column with the temperature data is not named " +
                  "'Temperature'. This will leed to problems in the " +
                  "simulation. "
                  "Please rename the column with the temperature data to " +
                  "'Temperature' and reload the file.")
        return file

    def create_time_series(self, time_step_size: str):
        """
        Creating a time series of the surface temperature with a given time
        step size.
        This is done by slicing the input data into several pieces of the same
        size defined by the time step size and using the average of the
        temperature of the slices as new values for the new time series.

        Note: The initial time series will be overwritten with the new time
        series with the corresponding step size and the averaged temperature
        data.

        Input:
            time step size (str) as being one of the following: 'every minute',
            'hourly', 'daily', 'weekly' or 'monthly'

        Output:
            surface temperature as time series with a certain time step size
        """
        time_step_dictionary = {"every minute": "min",
                                "hourly": "h",
                                "daily": "D",
                                "weekly": "W",
                                "monthly": "ME"}
        if time_step_size in time_step_dictionary:
            frequency = time_step_dictionary[time_step_size]
            time_series = self.data.resample(frequency).mean()
            time_series.dropna(subset="Temperature", inplace=True)
            self.data = time_series
        else:
            print("The time step size is not callable. Either change the " +
                  "code or choose one of the following time step sizes:" +
                  "'hourly', 'daily', 'weekly' or 'monthly'." +
                  "The time series was not changed yet.")

    def moving_average(self, time_step_size):
        """
        Transforming the time series into a moving average.
        Depending on the time step size, the moving average will use the time
        window of 2 days for 'every minute', 'hourly', 7 days for 'daily' and
        60 days for 'weekly' and 'monthly'.

        Note: The initial time series will be overwritten with the new time
        series with the corresponding moving average of the temperature for the
        corresponding step size.

        Input:
            time step size (str) as being one of the following: 'every minute',
            'hourly', 'daily', 'weekly' or 'monthly'

        Output:
            moving average of the surface temperature as time series with a
            certain time step size
        """
        if time_step_size in ["every minute", "hourly"]:
            mov_av = self.data["Temperature"].rolling(window="2D",
                                                      center=True).mean()
        elif time_step_size in ["daily"]:
            mov_av = self.data["Temperature"].rolling(window="7D",
                                                      center=True).mean()
        elif time_step_size in ["weekly", "monthly"]:
            mov_av = self.data["Temperature"].rolling(window="60D",
                                                      center=True).mean()
        self.data["Temperature"] = mov_av

    def repeat_time_period_annually(self, data_frequency: str,
                                    repetitions_years):
        """
        Prolongs the surface temperature input by repeating it over the given
        number of years. If the time series does not cover a whole year, the
        missing time will be linearly interpolated between the two ends of the
        measurement.
        The frequency of the measured data has to be known: every minute
        ('min'), hourly ('h'), daily ('D'), weekly ('W').

        Input:
            time steps of data input (str), number of years the series should
            be expanded to

        Output:
            expanded surface temperature time series

        Note:
            February 29th will never be used. If it is in the input the data
            from that day will be deleted and in all leap years this day will
            not occur and will therefore not contribute any surface temperature
            input.
            The input has to be shorter than a year, at maximum exactly one
            year.
        """
        start_date = self.data.index[0]
        end_date_1y = start_date + pd.DateOffset(years=1)
        end_date_1y = end_date_1y - pd.Timedelta(1, data_frequency)
        full_year = pd.date_range(start_date, end_date_1y, freq=data_frequency)
        missing_dates = full_year.difference(self.data.index)
        if len(missing_dates) > 0:
            self.data = self.data.reindex(full_year)
            self.data.iloc[-1, 0] = self.data.iloc[0, 0]
            self.data.loc[:, "Temperature"] = (self.data["Temperature"]
                                               .interpolate())
            print("The surface temperature was extended to one year with " +
                  "linear interpolation.")
        else:
            print("The surface temperature time series was already one year.")
        # Delete entries of February 29th in data if present
        if self.data.loc[((self.data.index.month == 2) &
                          (self.data.index.day == 29))] is not None:
            drop_idx = self.data.index[((self.data.index.month == 2) &
                                        (self.data.index.day == 29))].tolist()
            self.data.drop(drop_idx, inplace=True)
        # Repeat data
        self.data = pd.concat([self.data] * repetitions_years)
        # Create new indices and delete February 29th there as well
        end_date = start_date + pd.DateOffset(years=repetitions_years)
        end_date = end_date - pd.Timedelta(1, data_frequency)
        new_index = pd.date_range(start_date, end_date, freq=data_frequency)
        new_index = new_index[((new_index.month != 2) | (new_index.day != 29))]
        # Set correct index
        self.data.index = new_index
        print(f"The surface time series was extended to {repetitions_years} " +
              "years by repetition of the measurement.")
        print()

    def plot_surface_temperature(self):
        """
        Plotting the surface temperature over the time period.
        """
        self.data.plot(y="Temperature", xlabel="Time",
                       ylabel="Temperature [°C]", title="Surface Temperature")
        plt.show()
