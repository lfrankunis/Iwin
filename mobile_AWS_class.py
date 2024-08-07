# -*- coding: utf-8 -*-
"""
Created on Sat Apr 10 12:21:23 2021
@author: unismet
"""


import numpy as np
import datetime
from netCDF4 import Dataset
import sys
import os
import utm
import copy
import pandas as pd
from collections import deque
from io import StringIO
from pyproj import Proj
import yaml

import matplotlib.pyplot as plt

class mobile_AWS():

    def __init__(self,  station="MSBerg",                        # station number: 1883 or 1872 or 1924
                        resolution="1min",                # temporal resolution of the data
                        starttime="202104090800",         # INPUT Define start- and end-points in time for retrieving the data
                        endtime="202104091800",           # Format: YYYYmmddHHMM
                        variables=['temperature', 'air_pressure', 'relative_humidity', 'wind_speed', 'wind_speed_raw', 'wind_direction_raw', 'wind_direction', 'latitude', 'longitude'],
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
        self.data = self.read_mobile_AWS(day_vector.pop(0), file_type)
        if file_type == "nc":
            for day in day_vector:
                new_data = self.read_mobile_AWS(day, file_type)
                for vari in self.variables:
                    self.data[vari] = np.concatenate((self.data[vari], new_data[vari]), axis=0)
                self.data['time'] = np.append(self.data['time'], new_data['time'])
                self.data['local_time'] = np.append(self.data['local_time'], new_data['local_time'])

        self.resolution_min = (self.data['time'][-1] - self.data['time'][-2]).seconds / 60.


        ind = np.where((self.data['time'] >= self.starttime) & (self.data['time'] <= self.endtime))[0]
        for vari in self.variables:
            self.data[vari] = self.data[vari][ind]
        self.data['time'] = self.data['time'][ind]
        self.data['local_time'] = self.data['local_time'][ind]




    def read_mobile_AWS(self, day, file_type='nc'):
        """
        Method to read CARRA data for the given day
        """

        data = {}

        if self.file_type == 'nc':
            data_file = '{a}mobile_AWS_{s}/{b}{c:02d}{d:02d}/nc/{b}{c:02d}{d:02d}_mobile_AWS_{s}_Table_{r}.nc'.format(a=self.path, s=self.station, b=day.year, c=day.month, d=day.day, r=self.resolution)
            print('reading {a}'.format(a=os.path.basename(data_file)))

            with Dataset(data_file, 'r') as f:

                for vari in self.variables:
                    data[vari] = np.array(f.variables[vari][:])

                data['time'] = np.array([datetime.datetime.utcfromtimestamp(int(i)).replace(tzinfo=datetime.timezone.utc) for i in f.variables['time'][:]], dtype=datetime.datetime)
                data['local_time'] = np.array([datetime.datetime.fromtimestamp(int(i)) for i in f.variables['time'][:]], dtype=datetime.datetime)

        elif self.file_type == 'raw':
            
            with open("./config_data_files.yaml", "r", encoding='utf-8') as f:
                data_infiles = yaml.safe_load(f)

            infiles_station = data_infiles["mobile_stations"][self.station]
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
                
            df_data = pd.read_csv(StringIO(''.join(q)), header=0, skiprows=4, sep=",", na_values="NAN", names=list(col_names.keys()))

            # transfer timestamps into Python datetime objects
            df_data["time"] = pd.to_datetime(df_data["TIMESTAMP"]).dt.to_pydatetime()

            if "v3" in infile:
                df_data["GPRMC_speed_kn"] /= 1.94384

                df_data["GPRMC_longitude"] /= 100.
                df_data["GPRMC_longitude"] = (df_data["GPRMC_longitude"] // 1.) + (((df_data["GPRMC_longitude"] % 1.)*100.)/60.)
                df_data["GPRMC_latitude"] /= 100.
                df_data["GPRMC_latitude"] = (df_data["GPRMC_latitude"] // 1.) + (((df_data["GPRMC_latitude"] % 1.)*100.)/60.)

                df_data[((df_data["GPRMC_longitude"] < 13.38286) | (df_data["GPRMC_longitude"] > 17.46182) | (df_data["GPRMC_latitude"] < 77.95926) | (df_data["GPRMC_latitude"] > 78.85822))] = np.nan
                df_data.drop(columns=["GPS_heading", "GPS_speed", "GPS_location", "wind_speed_corrected_Avg", "wind_speed_corrected", "wind_speed_corrected_Max"],inplace=True)
                df_data = df_data.rename({"temperature": 'temperature',
                                    "temperature_Avg": 'temperature_Avg',
                                    "air_pressure_Avg": 'air_pressure',
                                    "relative_humidity": 'relative_humidity',
                                    "relative_humidity_Avg": 'relative_humidity_Avg',
                                    "dewpoint_temperature_Avg": 'dewpoint',
                                    "sea_surface_temperature": "sea_surface_temperature",
                                    "sea_surface_temperature_Avg": "sea_surface_temperature_Avg",
                                    "wind_speed_corrected_Avg": 'wind_speed',
                                    "wind_speed_corrected": 'wind_speed_smp',
                                    "wind_speed_corrected_Max": 'wind_speed_Max',
                                    "wind_speed_raw_Avg": 'wind_speed_raw',
                                    "wind_speed_raw": 'wind_speed_raw_Smp',
                                    "wind_speed_raw_Max": 'wind_speed_max_raw',
                                    "wind_direction_corrected_Avg": 'wind_direction',
                                    "wind_direction_corrected": 'wind_direction_Smp',
                                    "wind_direction_corrected_Std": 'wind_direction_std',
                                    "wind_direction_raw_Avg": 'wind_direction_raw',
                                    "wind_direction_raw": 'wind_direction_raw_Smp',
                                    "wind_direction_raw_Std": 'wind_direction_std_raw',
                                    "compass_heading": "compass_heading",
                                    "GPRMC_speed_kn": "GPS_speed",
                                    "HEHDT_heading": "GPS_heading",
                                    "GPRMC_latitude": "latitude",
                                    "GPRMC_longitude": "longitude"}, axis='columns')
                
                df_data["wind_speed"] = df_data["wind_speed_raw"]
                df_data["wind_direction"] = df_data["wind_direction_raw"]

            else:
                # extract lat, lon and lat from GPS Location
                latitude = np.ones((len(df_data["TIMESTAMP"]))) * np.nan
                longitude = np.ones((len(df_data["TIMESTAMP"]))) * np.nan
                altitude = np.ones((len(df_data["TIMESTAMP"]))) * np.nan
                for c, gps in enumerate(df_data["GPS_location"]):
                    if type(gps) == float:
                        latitude[c] = np.nan
                        longitude[c] = np.nan
                        altitude[c] = np.nan
                    else:
                        latitude[c] = float(gps.split(":")[0])
                        longitude[c] = float(gps.split(":")[1])
                        altitude[c] = float(gps.split(":")[2])
    
                        if ((latitude[c] < 77.95926) or (latitude[c] > 78.85822) or (longitude[c] < 13.38286) or (longitude[c] > 17.46182)):
                            latitude[c] = np.nan
                            longitude[c] = np.nan
                            altitude[c] = np.nan
    
                df_data['latitude'] = latitude
                df_data['longitude'] = longitude
                df_data['altitude'] = altitude
    
                df_data = df_data.rename({"temperature": 'temperature',
                                    "temperature_Avg": 'temperature_Avg',
                                    "air_pressure_Avg": 'air_pressure',
                                    "relative_humidity": 'relative_humidity',
                                    "relative_humidity_Avg": 'relative_humidity_Avg',
                                    "dewpoint_temperature_Avg": 'dewpoint',
                                    "sea_surface_temperature": "sea_surface_temperature",
                                    "sea_surface_temperature_Avg": "sea_surface_temperature_Avg",
                                    "wind_speed_corrected_Avg": 'wind_speed',
                                    "wind_speed_corrected": 'wind_speed_smp',
                                    "wind_speed_corrected_Max": 'wind_speed_Max',
                                    "wind_speed_raw_Avg": 'wind_speed_raw',
                                    "wind_speed_raw": 'wind_speed_raw_Smp',
                                    "wind_speed_raw_Max": 'wind_speed_max_raw',
                                    "wind_direction_corrected_Avg": 'wind_direction',
                                    "wind_direction_corrected": 'wind_direction_Smp',
                                    "wind_direction_corrected_Std": 'wind_direction_std',
                                    "wind_direction_raw_Avg": 'wind_direction_raw',
                                    "wind_direction_raw": 'wind_direction_raw_Smp',
                                    "wind_direction_raw_Std": 'wind_direction_std_raw',
                                    "GPS_heading": "GPS_heading",
                                    "GSP_speed": "GPS_speed",
                                    "compass_heading": "compass_heading"}, axis='columns')

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


    def filter_GPScoverage(self):
        """
        Method to delete time steps with missing GPS coverage
        """

        ind = np.where((~np.isnan(self.data['latitude'])) & (~np.isnan(self.data['longitude'])))[0]

        for vari in self.variables:
            self.data[vari] = self.data[vari][ind]
        self.data['time'] = self.data['time'][ind]
        self.data['local_time'] = self.data['local_time'][ind]

        return



    def delete_harbors(self):
        """
        Method to delete time steps when the station was in one of the 3 harbors (LYR, BB, PYR)
        """

        # LYR
        ind = list(np.where((self.data['latitude'] > 78.22745) & (self.data['latitude'] < 78.22878) & (self.data['longitude'] > 15.60521) & (self.data['longitude'] < 15.61387))[0])

        for vari in self.variables:
            self.data[vari] = np.delete(self.data[vari], ind)
        self.data['time'] = np.delete(self.data['time'], ind)

        # BB
        ind = list(np.where((self.data['latitude'] > 78.06336) & (self.data['latitude'] < 78.06433) & (self.data['longitude'] > 14.19790) & (self.data['longitude'] < 14.20329))[0])

        for vari in self.variables:
            self.data[vari] = np.delete(self.data[vari], ind)
        self.data['time'] = np.delete(self.data['time'], ind)

        # PYR
        ind = list(np.where((self.data['latitude'] > 78.65447) & (self.data['latitude'] < 78.65518) & (self.data['longitude'] > 16.37723) & (self.data['longitude'] < 16.38635))[0])

        for vari in self.variables:
            self.data[vari] = np.delete(self.data[vari], ind)
        self.data['time'] = np.delete(self.data['time'], ind)

        return



    def masks_for_harbors(self):
        """
        Method to get bool arrays marking the time steps when the station was in one of the 3 harbors (LYR, BB, PYR)
        """

        # LYR
        self.data["mask_LYR"] = np.zeros_like(self.data['time'], dtype=bool)
        ind = np.where((self.data['latitude'] > 78.22745) & (self.data['latitude'] < 78.22878) & (self.data['longitude'] > 15.60521) & (self.data['longitude'] < 15.61387))[0]
        self.data["mask_LYR"][ind] = True

        # BB
        self.data["mask_BB"] = np.zeros_like(self.data['time'], dtype=bool)
        ind = np.where((self.data['latitude'] > 78.06336) & (self.data['latitude'] < 78.06433) & (self.data['longitude'] > 14.19790) & (self.data['longitude'] < 14.20329))[0]
        self.data["mask_BB"][ind] = True

        # PYR
        self.data["mask_PYR"] = np.zeros_like(self.data['time'], dtype=bool)
        ind = np.where((self.data['latitude'] > 78.65447) & (self.data['latitude'] < 78.65518) & (self.data['longitude'] > 16.37723) & (self.data['longitude'] < 16.38635))[0]
        self.data["mask_PYR"][ind] = True

        # sailing: True, when the boat is out
        self.data["mask_sailing"] = ~(self.data['mask_LYR'] | self.data['mask_PYR'] | self.data['mask_BB'])

        self.variables.extend(['mask_LYR', 'mask_PYR', 'mask_BB', 'mask_sailing'])

        return
    
    
    
    def correct_winds(self, threshold=0.25, sat_comp=False):
        """
        Method to correct the wind speed and correction for the motion of the boats.
        """
        
        
        if sat_comp:
            
            u_raw = -np.abs(self.data["wind_speed_raw"]) * np.sin(np.deg2rad(self.data["wind_direction_raw"]))
            v_raw = -np.abs(self.data["wind_speed_raw"]) * np.cos(np.deg2rad(self.data["wind_direction_raw"]))
            u_georef = u_raw * np.cos(np.deg2rad(self.data["GPS_heading"])) + v_raw * np.sin(np.deg2rad(self.data["GPS_heading"]))
            v_georef = -u_raw * np.sin(np.deg2rad(self.data["GPS_heading"])) + v_raw * np.cos(np.deg2rad(self.data["GPS_heading"]))

            boat_dir = (self.data["GPS_heading"] + 180.) % 360.
            boat_u = -np.abs(self.data["GPS_speed"]) * np.sin(np.deg2rad(boat_dir))
            boat_v = -np.abs(self.data["GPS_speed"]) * np.cos(np.deg2rad(boat_dir))

            u_true = u_georef + boat_u
            v_true = v_georef + boat_v
            
        else:
        
            myProj = Proj("+proj=utm +zone=33 +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
            x, y = myProj(self.data["longitude"], self.data["latitude"])
            boat_u = np.gradient(x)/np.asarray(np.gradient(self.data["time"].astype("datetime64[s]")), dtype=float)
            boat_v = np.gradient(y)/np.asarray(np.gradient(self.data["time"].astype("datetime64[s]")), dtype=float)
            boat_speed = self.data["GPS_speed"] 
            boat_heading = self.data["GPS_heading"]                                                                         #OR boat_heading = self.data["compass_heading"]
            
            u = -np.abs(self.data["wind_speed"]) * np.sin(np.deg2rad(self.data["wind_direction"]))
            v = -np.abs(self.data["wind_speed"]) * np.cos(np.deg2rad(self.data["wind_direction"]))
            u_raw = -np.abs(self.data["wind_speed_raw"]) * np.sin(np.deg2rad(self.data["wind_direction_raw"]))
            v_raw = -np.abs(self.data["wind_speed_raw"]) * np.cos(np.deg2rad(self.data["wind_direction_raw"]))
            
            u_georef = u_raw * np.cos(np.deg2rad(boat_heading)) + v_raw * np.sin(np.deg2rad(boat_heading))
            v_georef = -u_raw * np.sin(np.deg2rad(boat_heading)) + v_raw * np.cos(np.deg2rad(boat_heading))
            
            u_shipcorrected = u_georef + boat_u
            v_shipcorrected = v_georef + boat_v
            
            u_true = copy.deepcopy(u)
            v_true = copy.deepcopy(v)
            u_true[boat_speed > threshold] = u_shipcorrected[boat_speed > threshold]
            v_true[boat_speed > threshold] = v_shipcorrected[boat_speed > threshold]
            
            
            
            
        self.data["wind_speed"] = np.sqrt(u_true**2. + v_true**2.)
        self.data["wind_direction"] = (np.rad2deg(np.arctan2(-u_true, -v_true)) + 360.) % 360.
        
        self.data["wind_speed"][~np.isfinite(self.data["wind_speed"])] = np.nan
        self.data["wind_direction"][~np.isfinite(self.data["wind_direction"])] = np.nan
        
        return



    def calculate_windvector_components(self, corrected=True):
        """
        Method to calculate wind speed and wind direction from the vectorial components u and v.
        """

        if corrected:
            self.data['u'] = -np.abs(self.data['wind_speed']) * np.sin(np.deg2rad(self.data['wind_direction']))
            self.data['v'] = -np.abs(self.data['wind_speed']) * np.cos(np.deg2rad(self.data['wind_direction']))
            self.variables.extend(['u', 'v'])

        else:
            self.data['u_raw'] = -np.abs(self.data['wind_speed_raw']) * np.sin(np.deg2rad(self.data['wind_direction_raw']))
            self.data['v_raw'] = -np.abs(self.data['wind_speed_raw']) * np.cos(np.deg2rad(self.data['wind_direction_raw']))
            self.variables.extend(['u_raw', 'v_raw'])

        return


    def calculate_specific_humidity(self):
        """
        Method to calulate the specific humidity (g/kg) from the specific humidity, temperature and pressure measurements.
        """

        e = 0.01*self.data['relative_humidity']*(6.112 * np.exp((17.62*self.data['temperature'])/(243.12+self.data['temperature'])))

        self.data['specific_humidity'] = 1000.*(0.622*e)/(self.data['air_pressure']-0.378*e)
        self.variables.append('specific_humidity')

        return



    def calculate_wind_sector(self, corrected=True):
        """
        Method to calculate the wind direction sector (N, NE, E, SE, ...).
        """

        if corrected:
            self.data['wind_sector'] = (np.array((self.data['wind_direction']/45.)+.5, dtype=int) % 8) + 1
            self.variables.append('wind_sector')
        else:
            self.data['wind_sector_raw'] = (np.array((self.data['wind_direction_raw']/45.)+.5, dtype=int) % 8) + 1
            self.variables.append('wind_sector_raw')
        return


    def calculate_wind_in_knots(self, corrected=True):
        """
        Method to calculate the wind speeds (u, v, wspeed) in knots.
        Must be called after the components are calculated.
        """

        if corrected:
            self.data['u_knts'] = 1.94384 * self.data['u']
            self.data['v_knts'] = 1.94384 * self.data['v']
            self.data['wind_speed_knts'] = np.sqrt(self.data['u_knts']**2. + self.data['v_knts']**2.)
            self.variables.extend(['u_knts', 'v_knts', 'wind_speed_knts'])

        else:
            self.data['u_knts_raw'] = 1.94384 * self.data['u_raw']
            self.data['v_knts_raw'] = 1.94384 * self.data['v_raw']
            self.data['wind_speed_knts_raw'] = np.sqrt(self.data['u_knts_raw']**2. + self.data['v_knts_raw']**2.)
            self.variables.extend(['u_knts_raw', 'v_knts_raw', 'wind_speed_knts_raw'])

        return