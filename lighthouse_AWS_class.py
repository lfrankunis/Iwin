# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 11:34:12 2021

@author: unismet
"""


import numpy as np
import datetime
from netCDF4 import Dataset
import sys
import os
import pandas as pd
from collections import deque
import yaml
from io import StringIO


class lighthouse_AWS():

    def __init__(self,  station="Bohemanneset",                        # station number: 1885
                        resolution="1min",                # temporal resolution of the data
                        starttime="202104090800",         # INPUT Define start- and end-points in time for retrieving the data
                        endtime="202104091800",           # Format: YYYYmmddHHMM
                        variables=['temperature', 'air_pressure', 'relative_humidity', 'wind_speed', 'wind_direction'],
                        file_type="nc", total_time=datetime.timedelta(days=1), path=False
                        ):

        self.starttime = datetime.datetime.strptime(starttime, "%Y%m%d%H%M").replace(tzinfo=datetime.timezone.utc)
        self.endtime = datetime.datetime.strptime(endtime, "%Y%m%d%H%M").replace(tzinfo=datetime.timezone.utc)
        day_vector = list(np.arange(self.starttime, self.endtime+datetime.timedelta(hours=1), datetime.timedelta(days=1), dtype=datetime.datetime))

        self.total_time = total_time
        self.file_type = file_type

        self.variables = variables
        self.resolution = resolution
        self.station = station

        if not path:
            with open("./config_paths.yaml", "r", encoding='utf-8') as f:
                paths = yaml.safe_load(f)
            self.path = paths['local_data']
        else:
            self.path = path

        # Load all files needed
        self.data = self.read_lighthouse_AWS(day_vector.pop(0), file_type)
        if file_type == "nc":
            for day in day_vector:
                new_data = self.read_lighthouse_AWS(day, file_type)
                for vari in self.variables:
                    self.data[vari] = np.concatenate((self.data[vari], new_data[vari]), axis=0)
                self.data['time'] = np.append(self.data['time'], new_data['time'])

        self.resolution_min = (self.data['time'][-1] - self.data['time'][-2]).seconds / 60.


        ind = np.where((self.data['time'] >= self.starttime) & (self.data['time'] <= self.endtime))[0]
        for vari in self.variables:
            self.data[vari] = self.data[vari][ind]
        self.data['time'] = self.data['time'][ind]
        self.data['local_time'] = self.data['local_time'][ind]





    def read_lighthouse_AWS(self, day, file_type='nc'):
        """
        Method to read CARRA data for the given day
        """

        data = {}

        if self.file_type == 'nc':
            data_file = '{a}lighthouse_AWS_{s}/{b}{c:02d}{d:02d}/nc/{b}{c:02d}{d:02d}_lighthouse_AWS_{s}_Table_{r}.nc'.format(a=self.path, s=self.station, b=day.year, c=day.month, d=day.day, r=self.resolution)
            print('reading {a}'.format(a=os.path.basename(data_file)))

            with Dataset(data_file, 'r') as f:

                for vari in self.variables:
                    data[vari] = np.array(f.variables[vari][:])

                data['time'] = np.array([datetime.datetime.utcfromtimestamp(int(i)) for i in f.variables['time'][:]], dtype=datetime.datetime)
                data['local_time'] = np.array([datetime.datetime.fromtimestamp(int(i)) for i in f.variables['time'][:]], dtype=datetime.datetime)

        elif self.file_type == 'raw':
            
            with open("./config_data_files.yaml", "r", encoding='utf-8') as f:
                data_infiles = yaml.safe_load(f)

            infiles_station = data_infiles["lighthouse_stations"][self.station]
            for d, f in infiles_station.items():
                d1 = datetime.datetime.strptime(str(d)[:8], "%Y%m%d").replace(tzinfo=datetime.timezone.utc)
                d2 = datetime.datetime.strptime(str(d)[8:], "%Y%m%d").replace(tzinfo=datetime.timezone.utc)
                if ((day >= d1) & (day <= d2)):
                    infile = f
                    
            replace_settings = {"path": self.path, "resolution": self.resolution}
            for placeholder, value in replace_settings.items():
                infile = infile.replace(placeholder, value)
            
            
            print(f'reading {infile}')

            col_names = pd.read_csv(infile, header=1, sep=",", nrows=1).to_dict('records')[0]
            
            with open(infile, 'r') as f:
                q = deque(f, 5000)
                
            df_data = pd.read_csv(StringIO(''.join(q)), header=4, sep=",", na_values="NAN", names=list(col_names.keys()))

            # transfer timestamps into Python datetime objects
            df_data["time"] = pd.to_datetime(df_data["TIMESTAMP"]).dt.to_pydatetime()

            df_data = df_data.rename({"temperature_Avg": 'temperature',
                                "air_pressure_Avg": 'air_pressure',
                                "relative_humidity_Avg": 'relative_humidity',
                                "dewpoint_temperature_Avg": 'dewpoint',
                                "wind_speed_Avg": 'wind_speed',
                                "wind_speed_Max": 'wind_speed_max',
                                "wind_direction_Avg": 'wind_direction',
                                "wind_direction_Std": 'wind_direction_std',
                                "BattV_Min": "battery"}, axis='columns')

            df_data = df_data.dropna()

            data_all = df_data.to_dict('list')

            for vari in self.variables:
                data[vari] = np.array(data_all[vari])

            data['time'] = np.array(data_all['time'])
            for i in range(len(data['time'])):
                data['time'][i] = data['time'][i].replace(tzinfo=datetime.timezone.utc)
            data['local_time'] = np.array([i.tz_convert("Europe/Berlin") for i in data["time"]])

        return data



    def only_latest_data(self, period):
        """
        Method to use only the latest data, starting from 'period' before the last value.
        """

        ind = np.where((self.data['time'] >= self.data['time'][-1]-period))[0]
        for vari in self.variables:
            self.data[vari] = self.data[vari][ind]
        self.data['time'] = self.data['time'][ind]
        self.data['local_time'] = self.data['local_time'][ind]

        return



    def calculate_windvector_components(self):
        """
        Method to calculate wind speed and wind direction from the vectorial components u and v.
        """

        self.data['u'] = -np.abs(self.data['wind_speed']) * np.sin(np.deg2rad(self.data['wind_direction']))
        self.data['v'] = -np.abs(self.data['wind_speed']) * np.cos(np.deg2rad(self.data['wind_direction']))
        self.variables.extend(['u', 'v'])

        return


    def calculate_specific_humidity(self):
        """
        Method to calulate the specific humidity (g/kg) from the specific humidity, temperature and pressure measurements.
        """

        e = 0.01*self.data['relative_humidity']*(6.112 * np.exp((17.62*self.data['temperature'])/(243.12+self.data['temperature'])))

        self.data['specific_humidity'] = 1000.*(0.622*e)/(self.data['air_pressure']-0.378*e)
        self.variables.append('specific_humidity')

        return



    def calculate_wind_sector(self):
        """
        Method to calculate the wind direction sector (N, NE, E, SE, ...).
        """

        self.data['wind_sector'] = (np.array((self.data['wind_direction']/45.)+.5, dtype=int) % 8) + 1
        self.variables.append('wind_sector')

        return



    def calculate_wind_in_knots(self):
        """
        Method to calculate the wind speeds (u, v, wspeed) in knots.
        Must be called after the components are calculated.
        """

        self.data['u_knts'] = 1.94384 * self.data['u']
        self.data['v_knts'] = 1.94384 * self.data['v']
        self.data['wind_speed_knts'] = np.sqrt(self.data['u_knts']**2. + self.data['v_knts']**2.)
        self.variables.extend(['u_knts', 'v_knts', 'wind_speed_knts'])

        return