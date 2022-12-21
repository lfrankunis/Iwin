# -*- coding: utf-8 -*-
"""
Created on Sat Mar 27 14:37:44 2021
@author: unismet
"""

import numpy as np
import pandas as pd
import datetime
import yaml
import shutil
from pyproj import Proj
import copy
import os
import json


def restructure_mobile_AWS(from_time, to_time, station="1883", resolution="10min"):
    """
    Function to restructure the mobile AWS data into e.g. daily files
    Parameters
    ----------
    path : str
        string defining the path where to find the input data and save the output data
    station : str
        string specifying the station ("1883" or "1872" or "1924")
    resolution : str
        string specifying the resolution of the data ("20sec", "1min", or "10min")
    from_time : datetime object
        start time of the restructured file
    to_time : datetime object
        end time of the restructured file
    Returns a netCDF4 file
    -------
    In case the restructured interval is not one day, remember to change the output file names
    """

    threshold = 0.25
    
    # define path to the data folder
    with open("./config_paths.yaml", "r", encoding='utf-8') as f:
        paths = yaml.safe_load(f)

    if from_time < datetime.datetime(2022,5,17):
        infile = f"{paths['local_data']}mobile_AWS_{station}/mobile_AWS_{station}_Table_{resolution}_v1_{from_time.year}.dat"
    elif from_time < datetime.datetime(2023,1,1):
        infile = f"{paths['local_data']}mobile_AWS_{station}/mobile_AWS_{station}_Table_{resolution}_v2_{from_time.year}.dat"
    else:
        infile = f"{paths['local_data']}mobile_AWS_{station}/mobile_AWS_{station}_Table_{resolution}.dat"
    outfile_nc = f"{paths['local_storage']}mobile_AWS_{station}/{resolution}/mobile_AWS_{station}_Table_{resolution}_{from_time.year}{from_time.month:02d}{from_time.day:02d}.nc"
    backup_nc = f"{paths['harddrive']}mobile_AWS_{station}/{resolution}/mobile_AWS_{station}_Table_{resolution}_{from_time.year}{from_time.month:02d}{from_time.day:02d}.nc"

    print("restructuring {r} data from AWS {s}...".format(r=resolution, s=station))

    col_names = pd.read_csv(infile, header=1, sep=",", nrows=1).to_dict('records')[0]

    time_avail = pd.to_datetime(pd.read_csv(infile, header=3, sep=",", usecols=[0], na_values="NAN").iloc[:,0]).dt.to_pydatetime()
    ind = np.where((time_avail >= from_time) & (time_avail < to_time))[0]
    
    if len(ind) == 0:
        return time_avail[-1]
    
    data = pd.read_csv(infile, header=3+ind[0], nrows=len(ind), sep=",", na_values="NAN", names=list(col_names.keys()))
    
    if from_time < datetime.datetime(2022,5,17):
        new_names = {"TIMESTAMP": "time",
                     'Wind_Speed_S_WVT': 'wind_speed_raw',
                     'Wind_Direction_D1_WVT': 'wind_direction_raw',
                     'Wind_Direction_SDI_WVT': 'wind_direction_raw_Std',
                     'Wind_Speed_raw_Max': 'wind_speed_raw_Max',
                     'Wind_Speed_Corrected_S_WVT': 'wind_speed_corrected',
                     'Wind_Direction_Corrected_D1_WVT': 'wind_direction_corrected',
                     'Wind_Direction_Corrected_SDI_WVT': 'wind_direction_corrected_Std',
                     'Wind_Speed_corrected_Max': 'wind_speed_corrected_Max',
                     'Ambient_Temperature': 'temperature',
                     'Relative_Humidity': 'relative_humidity',
                     'Barometric_Pressure': 'air_pressure',
                     'GPS_Location': 'GPS_location',
                     'Sensor_status': "sensor_status",
                     "GPS_speed": "GPS_speed",
                     "GPS_heading": "GPS_heading"}
        cols_to_drop = [i for i in list(data.columns) if i not in list(new_names.keys())]
        data.drop(columns=cols_to_drop, inplace=True)
        data.rename(new_names, axis=1, inplace=True)

        
        variable_attributes = {
        "units": {'time': "seconds since 1970-01-01T00:00:00Z", 'wind_speed_raw': "m s-1",
                  'wind_direction_raw': "degree", 'wind_direction_raw_Std': "degree", 'wind_speed_raw_Max': "m s-1",
                  'wind_speed_corrected': "m s-1", 'wind_direction_corrected': "degree", 'wind_direction_corrected_Std': "degree", 'wind_speed_corrected_Max': "m s-1",
                  'temperature': "degree_C", 'relative_humidity': "percent", 'air_pressure': "hPa",
                  'GPS_heading': "degree", 'GPS_speed': "m s-1", 'latitude': "degree_N", 'longitude': "degree_E", 'altitude': "m"},
        "long_names" : {'time': "UTC time", 'sensor_status': "sensor status",
                        'wind_speed_raw': "raw wind speed averaged over the sampling interval",
                        'wind_direction_raw': "raw wind direction averaged over the sampling interval",
                        'wind_direction_raw_Std': "standard devitation of the raw wind speed during the sampling interval",
                        'wind_speed_raw_Max': "maximum raw wind speed during the sampling interval",
                        'wind_speed_corrected': "wind speed averaged over the sampling interval, corrected for the movement of the boat",
                        'wind_direction_corrected': "wind direction averaged over the sampling interval, corrected for the movement of the boat",
                        'wind_direction_corrected_Std': "standard deviation of the wind direction during the sampling interval, corrected for the movement of the boat",
                        'wind_speed_corrected_Max': "maximum wind speed during the sampling interval, corrected for the movement of the boat",
                        'temperature': "air temperature sampled at the respective timestamp",
                        'relative_humidity': "air relative humidity sampled at the respective timestamp",
                        'air_pressure': "air pressure sampled at the respective timestamp",
                        'GPS_heading': "heading of the boat, retrieved from the GPS", 'GPS_speed': "speed of the boat, retrieved from the GPS",
                        'latitude': "latitude", 'longitude': "longitude", 'altitude': "height of the sensor over ground, retrieved from the GPS",
                        "exhaust_plume_influence": "flag indicating a possile contamination of the measurements by the exhaust plume"},
        "standard_names": {'time': "time", 'sensor_status': "status_flag", 'wind_speed_raw': "wind_speed", 'wind_direction_raw': "wind_from_direction",
                           'wind_speed_raw_Max': "wind_speed", 
                           'wind_speed_corrected': "wind_speed", 'wind_direction_corrected': "wind_from_direction",
                           'wind_speed_corrected_Max': "wind_speed",
                           'temperature': "air_temperature", 'relative_humidity': 'relative_humidity', 'air_pressure': "air_pressure",
                           'GPS_heading': "platform_azimuth_angle", 'GPS_speed': "platform_speed_wrt_ground",
                           'latitude': "latitude", 'longitude': "longitude", 'altitude': "altitude",
                           "exhaust_plume_influence": "status_flag"},
        "flag_values": {"sensor_status": 'OK, Unknown Fault, Wind Fault', "exhaust_plume_influence": "0, 1"},
        "flag_meanings": {"exhaust_plume_influence": "measurements_not_impacted_by_exhaust_plume, measurements_impacted_by_exhaust_plume"}}
        
    else:
        new_names = {"TIMESTAMP": "time",
                     'wind_speed_raw_Avg': 'wind_speed_raw',
                     'wind_direction_raw_Avg': 'wind_direction_raw',
                     'wind_direction_raw_Std': 'wind_direction_raw_Std',
                     'wind_speed_raw_Max': 'wind_speed_raw_Max',
                     'wind_speed_corrected_Avg': 'wind_speed_corrected',
                     'wind_direction_corrected_Avg': 'wind_direction_corrected',
                     'wind_direction_corrected_Std': 'wind_direction_corrected_Std',
                     'wind_speed_corrected_Max': 'wind_speed_corrected_Max',
                     'temperature': 'temperature',
                     'relative_humidity': 'relative_humidity',
                     'air_pressure': 'air_pressure',
                     'GPS_location': 'GPS_location',
                     'sensor_status': "sensor_status",
                     'compass_heading': 'compass_heading',
                     "GPS_speed": "GPS_speed",
                     "GPS_heading": "GPS_heading"}
        cols_to_drop = [i for i in list(data.columns) if i not in list(new_names.keys())]
        data.drop(columns=cols_to_drop, inplace=True)
        data.rename(new_names, axis=1, inplace=True)

        
        variable_attributes = {
        "units": {'time': "seconds since 1970-01-01T00:00:00Z", 'wind_speed_raw': "m s-1",
                  'wind_direction_raw': "degree", 'wind_direction_raw_Std': "degree", 'wind_speed_raw_Max': "m s-1",
                  'wind_speed_corrected': "m s-1", 'wind_direction_corrected': "degree", 'wind_direction_corrected_Std': "degree",
                  'wind_speed_corrected_Max': "m s-1",
                  'temperature': "degree_C", 'relative_humidity': "percent", 'air_pressure': "hPa",
                  'GPS_heading': "degree", 'GPS_speed': "m s-1", 'compass_heading': "degree", 'latitude': "degree_N", 'longitude': "degree_E", 'altitude': "m"},
        
        "long_name" : {'time': "UTC time", 'sensor_status': "sensor status",
                        'wind_speed_raw': "raw wind speed averaged over the sampling interval",
                        'wind_direction_raw': "raw wind direction averaged over the sampling interval",
                        'wind_direction_raw_Std': "standard devitation of the raw wind speed during the sampling interval",
                        'wind_speed_raw_Max': "maximum raw wind speed during the sampling interval",
                        'wind_speed_corrected': "wind speed averaged over the sampling interval, corrected for the movement of the boat",
                        'wind_direction_corrected': "wind direction averaged over the sampling interval, corrected for the movement of the boat",
                        'wind_direction_corrected_Std': "standard deviation of the wind direction during the sampling interval, corrected for the movement of the boat",
                        'wind_speed_corrected_Max': "maximum wind speed during the sampling interval, corrected for the movement of the boat",
                        'temperature': "air temperature averaged over the sampling interval",
                        'relative_humidity': "air relative humidity averaged over the sampling interval",
                        'air_pressure': "air pressure averaged over the sampling interval",
                        'GPS_heading': "heading of the boat, retrieved from the GPS", 'GPS_speed': "speed of the boat, retrieved from the GPS",
                        'compass_heading': "heading of the boat, retrieved from the compass",
                        'latitude': "latitude", 'longitude': "longitude", 'altitude': "height of the sensor over ground, retrieved from the GPS",
                        "exhaust_plume_influence": "flag indicating a possile contamination of the measurements by the exhaust plume"},
        
        "standard_name": {'time': "time", 'sensor_status': "status_flag",
                           'wind_speed_raw_Max': "wind_speed", 'wind_speed_raw': "wind_speed",
                           'wind_direction_raw': "wind_from_direction",
                           'wind_speed_corrected_Max': "wind_speed", 'wind_speed_corrected': "wind_speed",
                           'wind_direction_corrected': "wind_from_direction",
                           'temperature': "air_temperature",
                           'relative_humidity': 'relative_humidity',
                           'air_pressure': "air_pressure", 'GPS_heading': "platform_azimuth_angle", 'GPS_speed': "platform_speed_wrt_ground",
                           'compass_heading': "platform_azimuth_angle", 'latitude': "latitude", 'longitude': "longitude", 'altitude': "altitude",
                           "exhaust_plume_influence": "status_flag"},
        "flag_values": {"sensor_status": 'OK, Unknown Fault, Wind Fault', "exhaust_plume_influence": "0, 1"},
        "flag_meanings": {"exhaust_plume_influence": "measurements_not_impacted_by_exhaust_plume, measurements_impacted_by_exhaust_plume"}}
        

    # transfer timestamps into Python datetime objects
    data["time"] = [dt.replace(tzinfo=datetime.timezone.utc).timestamp() for dt in pd.to_datetime(data["time"]).dt.to_pydatetime()]

    # extract lat, lon and lat from GPS Location
    latitude = np.ones((len(data["time"]))) * np.nan
    longitude = np.ones((len(data["time"]))) * np.nan
    altitude = np.ones((len(data["time"]))) * np.nan
    for c, gps in enumerate(data["GPS_location"]):
        if type(gps) == float:
            latitude[c] = np.nan
            longitude[c] = np.nan
            altitude[c] = np.nan
        else:
            latitude[c] = float(gps.split(":")[0])
            longitude[c] = float(gps.split(":")[1])
            altitude[c] = float(gps.split(":")[2])
            
            # filter data with erroneous GPS positions
            if ((latitude[c] < 77.95926) or (latitude[c] > 78.85822) or (longitude[c] < 13.38286) or (longitude[c] > 17.46182)):
                latitude[c] = np.nan
                longitude[c] = np.nan
                altitude[c] = np.nan
                
                
    # add lon, lat and alt to dataframe
    data["latitude"] = latitude
    data["longitude"] = longitude
    data["altitude"] = altitude


    # correct wind data for motion of the boat
    myProj = Proj("+proj=utm +zone=33 +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
    x, y = myProj(longitude, latitude)
    boat_u = np.gradient(x)/np.asarray(np.gradient(data["time"].values.astype("datetime64[s]")), dtype=float)
    boat_v = np.gradient(y)/np.asarray(np.gradient(data["time"].values.astype("datetime64[s]")), dtype=float)
    if from_time < datetime.datetime(2022,5,17):
        data["GPS_speed"]  = np.sqrt(boat_u**2. + boat_v**2.)
        data["GPS_heading"] = (((np.rad2deg(np.arctan2(-boat_u, -boat_v)) + 360.) % 360.) + 180.) % 360.

    boat_speed = data["GPS_speed"] 
    boat_heading = data["GPS_heading"] # OR boat_heading = data["compass_heading"]
    
    # average wind 
    u = -np.abs(data["wind_speed_corrected"]) * np.sin(np.deg2rad(data["wind_direction_corrected"]))
    v = -np.abs(data["wind_speed_corrected"]) * np.cos(np.deg2rad(data["wind_direction_corrected"]))
    u_raw = -np.abs(data["wind_speed_raw"]) * np.sin(np.deg2rad(data["wind_direction_raw"]))
    v_raw = -np.abs(data["wind_speed_raw"]) * np.cos(np.deg2rad(data["wind_direction_raw"]))

    u_georef = u_raw * np.cos(np.deg2rad(boat_heading)) + v_raw * np.sin(np.deg2rad(boat_heading))
    v_georef = -u_raw * np.sin(np.deg2rad(boat_heading)) + v_raw * np.cos(np.deg2rad(boat_heading))

    u_shipcorrected = u_georef + boat_u
    v_shipcorrected = v_georef + boat_v

    u_true = copy.deepcopy(u)
    v_true = copy.deepcopy(v)
    u_true[boat_speed > threshold] = u_shipcorrected[boat_speed > threshold]
    v_true[boat_speed > threshold] = v_shipcorrected[boat_speed > threshold]
    

    data["wind_speed_corrected"] = np.sqrt(u_true**2. + v_true**2.)
    data["wind_direction_corrected"] = (np.rad2deg(np.arctan2(-u_true, -v_true)) + 360.) % 360.
    
    data.loc[~np.isfinite(data["wind_speed_corrected"]), "wind_speed_corrected"] = np.nan
    data.loc[~np.isfinite(data["wind_direction_corrected"]), "wind_direction_corrected"] = np.nan
    
    if ((station == "1924") & (from_time < datetime.datetime(2022,8,31))):
        data.loc[boat_speed <= threshold, "wind_direction_corrected"] = np.nan
    
    # maximum wind speed
    u = -np.abs(data["wind_speed_corrected_Max"]) * np.sin(np.deg2rad(data["wind_direction_corrected"]))
    v = -np.abs(data["wind_speed_corrected_Max"]) * np.cos(np.deg2rad(data["wind_direction_corrected"]))
    u_raw = -np.abs(data["wind_speed_raw_Max"]) * np.sin(np.deg2rad(data["wind_direction_raw"])) 
    v_raw = -np.abs(data["wind_speed_raw_Max"]) * np.cos(np.deg2rad(data["wind_direction_raw"]))

    u_georef = u_raw * np.cos(np.deg2rad(boat_heading)) + v_raw * np.sin(np.deg2rad(boat_heading))
    v_georef = -u_raw * np.sin(np.deg2rad(boat_heading)) + v_raw * np.cos(np.deg2rad(boat_heading))

    u_shipcorrected = u_georef + boat_u
    v_shipcorrected = v_georef + boat_v

    u_true = copy.deepcopy(u)
    v_true = copy.deepcopy(v)
    u_true[boat_speed > threshold] = u_shipcorrected[boat_speed > threshold]
    v_true[boat_speed > threshold] = v_shipcorrected[boat_speed > threshold]

    data["wind_speed_corrected_Max"] = np.sqrt(u_true**2. + v_true**2.)
    
    data.loc[~np.isfinite(data["wind_speed_corrected_Max"]), "wind_speed_corrected_Max"] = np.nan

    data.set_index("time", inplace=True)
    data.drop(columns=["GPS_location"], inplace=True)
    
    data = data[~data.index.duplicated(keep='first')]
    
    start_data_coverage = datetime.datetime.fromtimestamp(data.index[0], datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    end_data_coverage = datetime.datetime.fromtimestamp(data.index[-1], datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    
    
    # filter data for unphysical outliers
    for vari in data.columns:
        if "wind_direction" in vari:
            data.loc[data[vari]<0., vari] = np.nan
            data.loc[data[vari]>360., vari] = np.nan
        elif "wind_speed" in vari:
            data.loc[data[vari]<0., vari] = np.nan
            data.loc[data[vari]>80., vari] = np.nan
        elif 'temperature' in vari:
            data.loc[data[vari]<-80., vari] = np.nan
            data.loc[data[vari]>40., vari] = np.nan
        elif "relative_humidity" in vari:
            data.loc[data[vari]<0., vari] = np.nan
            data.loc[data[vari]>100., vari] = np.nan
        elif "pressure" in vari:
            data.loc[data[vari]<800., vari] = np.nan
            data.loc[data[vari]>1100., vari] = np.nan
    
    # determine exhaust plume flag
    exhaust_sectors = {"1872": (215., 235.),
                       "1924": (170., 190.),
                       "1883": (-3., -2.)}      # dummy for MS Bard

    if ((station == "1924") & (from_time > datetime.datetime(2022,8,30)) & (from_time < datetime.datetime(2022,11,20))):
        data["exhaust_plume_influence"] = np.asarray(((data['wind_direction_raw'] < exhaust_sectors["1883"][1]) & 
                                                      (data['wind_direction_raw'] > exhaust_sectors["1883"][0])), dtype=int)
    else:
        data["exhaust_plume_influence"] = np.asarray(((data['wind_direction_raw'] < exhaust_sectors[station][1]) & 
                                                      (data['wind_direction_raw'] > exhaust_sectors[station][0])), dtype=int)
    
    data.fillna(-99999., inplace=True)
    
    ds = data.to_xarray()
    
    for vari in list(ds.variables):
        if vari == "time":
            ds[vari] = ds[vari].astype('float64', keep_attrs=True)
        elif vari in ["exhaust_plume_influence", 'sensor_status']:
            ds[vari] = ds[vari].astype('int32', keep_attrs=True)
        else:
            ds[vari] = ds[vari].astype('float32', keep_attrs=True)
        attris = {}
        for a, d in variable_attributes.items():
            try:
                attris[a] = d[vari]
            except KeyError:
                pass
        ds[vari].attrs = attris
        
    
    # Assign global attributes
    dtnow = datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")
    
    if ((station == "1924") & (from_time > datetime.datetime(2022,8,30)) & (from_time < datetime.datetime(2022,11,20))):
        boat_names = {"1883": {"name": "MS Bard", "installation":  datetime.datetime(2021,5,6,0,0,0,tzinfo=datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")},
                      "1872": {"name": "MS Polargirl", "installation":  datetime.datetime(2021,5,29,0,0,0,tzinfo=datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")},
                      "1924": {"name": "MS Bard", "installation":  datetime.datetime(2022,3,21,0,0,0,tzinfo=datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")}}
    else:
        boat_names = {"1883": {"name": "MS Bard", "installation":  datetime.datetime(2021,5,6,0,0,0,tzinfo=datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")},
                      "1872": {"name": "MS Polargirl", "installation":  datetime.datetime(2021,5,29,0,0,0,tzinfo=datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")},
                      "1924": {"name": "MS Billefjord", "installation":  datetime.datetime(2022,3,21,0,0,0,tzinfo=datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")}}
        
    if from_time < datetime.datetime(2022,1,1):
        creator_insts = "The University Centre in Svalbard; The University Centre in Svalbard; The University Centre in Svalbard; MET Norway"
        creator_names = "L. Frank; F. Schalamon; M. O. Jonassen; T. Remes"
        creator_mails = "lukasf@unis.no; florina.schalamon@uni-graz.at; mariusj@unis.no; teresav@met.no"
        creator_urls = "https://orcid.org/0000-0003-1472-7967; https://orcid.org/0000-0002-2509-4133; https://orcid.org/0000-0002-4745-9009; https://orcid.org/0000-0002-6421-859X"
    elif from_time < datetime.datetime(2022,10,1):
        creator_insts = "The University Centre in Svalbard; The University Centre in Svalbard; The University Centre in Svalbard; The University Centre in Svalbard; MET Norway"
        creator_names = "L. Frank; F. Schalamon; A. Stenlund; M. O. Jonassen; T. Remes"
        creator_mails = "lukasf@unis.no; florina.schalamon@uni-graz.at; ;mariusj@unis.no; teresav@met.no"
        creator_urls = "https://orcid.org/0000-0003-1472-7967; https://orcid.org/0000-0002-2509-4133; https://orcid.org/0000-0003-4241-735X; https://orcid.org/0000-0002-4745-9009; https://orcid.org/0000-0002-6421-859X"
    else:
        creator_insts = "The University Centre in Svalbard; The University Centre in Svalbard; MET Norway"
        creator_names = "L. Frank; M. O. Jonassen; T. Remes"
        creator_mails = "lukasf@unis.no; mariusj@unis.no; teresav@met.no"
        creator_urls = "https://orcid.org/0000-0003-1472-7967; https://orcid.org/0000-0002-4745-9009; https://orcid.org/0000-0002-6421-859X"
    
    
    global_attributes = {"boat": boat_names[station]['name'],
        "title": f"Standard meteorological near-surface observations from {from_time.strftime('%Y-%m-%d')} measured onboard {boat_names[station]['name']} in Isfjorden, Svalbard.",
        "summary": "The file contains time series of the standard meteorological near-surface parameters temperature, humidity, pressure, wind speed and wind direction. The raw data is only filtered for obviously wrong GPS positions, otherwise the data is made available as is. Wind speed and wind direction are corrected for the horizontal movements of the boat using GPS data.",
        "keywords": "EARTH SCIENCE>ATMOSPHERE>ATMOSPHERIC TEMPERATURE>SURFACE TEMPERATURE>AIR TEMPERATURE, EARTH SCIENCE>ATMOSPHERE>ATMOSPHERIC TEMPERATURE>SURFACE TEMPERATURE>DEW POINT TEMPERATURE, EARTH SCIENCE>ATMOSPHERE>ATMOSPHERIC WATER VAPOR>WATER VAPOR INDICATORS>HUMIDITY>RELATIVE HUMIDITY, \
            EARTH SCIENCE>ATMOSPHERE>ATMOSPHERIC WINDS>SURFACE WINDS>WIND DIRECTION, EARTH SCIENCE>ATMOSPHERE>ATMOSPHERIC WINDS>SURFACE WINDS>WIND SPEED, EARTH SCIENCE>ATMOSPHERE>ATMOSPHERIC PRESSURE>SURFACE PRESSURE, \
            ACADEMIC>UNIS, GEOGRAPHIC REGION>POLAR, CONTINENT>EUROPE>NORTHERN EUROPE>SCANDINAVIA>NORWAY, \
            IN SITU/LABORATORY INSTRUMENTS>CURRENT/WIND METERS>SONIC ANEMOMETER, IN SITU/LABORATORY INSTRUMENTS>PRESSURE/HEIGHT METERS>PRESSURE SENSORS, \
            IN SITU/LABORATORY INSTRUMENTS>TEMPERATURE/HUMIDITY SENSORS>TEMPERATURE SENSORS, IN SITU/LABORATORY INSTRUMENTS>TEMPERATURE/HUMIDITY SENSORS>HUMIDITY SENSORS",
        "keywords_vocabulary": 'GCMD Science Keywords',
        "naming_authority": "NMDC",
        "data_type": "netCDF-4",
        "feature_type": "Time Series",
        "geospatial_lat_min": 77.95926,
        "geospatial_lat_max": 78.85822,
        "geospatial_lon_min": 13.38286,
        "geospatial_lon_max": 17.46182,
        "area": "Svalbard, Isfjorden",
        "time_coverage_start": start_data_coverage,
        "time_coverage_end": end_data_coverage,
        "Conventions": "CF-1.6, ACDD-1.3",
        'date_created': dtnow,
        'history': f'File created at {dtnow} using xarray in Python3.',
        "processing_level": "Known bad data has been replaced with null values, ranges applied, data has been scaled using contextual information (wind speed and direction corrected for horizontal movement of boat)",
        "creator_type": "group",
        "creator_institution": creator_insts,
        "creator_name": creator_names,
        "creator_email": creator_mails,
        "creator_url": creator_urls,
        "institution": "The University Centre in Svalbard (UNIS)",
        "project": "IWIN: Isfjorden Weather Information Network",
        "source": "Gill MaxiMet GMX 500",
        "platform": f"WATER-BASED PLATFORMS>VESSELS>SURFACE>{boat_names[station]['name']}",
        "instrument_type": "Gill MaxiMet GMX 500",
        "license": "https://creativecommons.org/licenses/by/4.0/",
        "iso_topic_category": "atmosphere",
        "principal_investigator": "L. Frank",
        "publisher_name": "MET Norway",
        "publisher_url": "https://thredds.met.no/thredds/catalog.html",
        "publisher_email": "thredds@met.no",
        "id": "???"}
    
    for a, value in global_attributes.items():
        ds.attrs[a] = value
    
    myencoding = {v: {'_FillValue': -99999., 'zlib': False} for v in list(ds.variables) if v not in ["time", "sensor_status", "exhaust_plume_influence"]}
    myencoding["exhaust_plume_influence"] = {"_FillValue": -99999}
    myencoding["time"] = {'_FillValue': None}
    
    ds.to_netcdf(outfile_nc, unlimited_dims=["time"], encoding=myencoding)
    
    if os.path.isfile(backup_nc):
        os.remove(backup_nc)
    shutil.copyfile(outfile_nc, backup_nc)

    return time_avail[-1]




def restructure_lighthouse_AWS(from_time, to_time, station="1885", resolution="10min"):
    """
    Function to restructure the lighthouse AWS data into e.g. daily files
    Parameters
    ----------
    path : str
        string defining the path where to find the input data and save the output data
    station : str
        string specifying the station (e.g. "1885")
    resolution : str
        string specifying the resolution of the data ("1min" or "10min")
    from_time : datetime object
        start time of the restructured file
    to_time : datetime object
        end time of the restructured file
    Returns a netCDF4 file
    -------
    In case the restructured interval is not one day, remember to change the output file names
    """
    
    # define path to the data folder
    with open("./config_paths.yaml", "r", encoding='utf-8') as f:
        paths = yaml.safe_load(f)

    if from_time < datetime.datetime(2022,5,8):
        infile = f"{paths['local_data']}lighthouse_AWS_{station}/lighthouse_AWS_{station}_Table_{resolution}_v1_{from_time.year}.dat"
    elif from_time < datetime.datetime(2022,11,18):
        infile = f"{paths['local_data']}lighthouse_AWS_{station}/lighthouse_AWS_{station}_Table_{resolution}_v2_{from_time.year}.dat"
    else:
        infile = f"{paths['local_data']}lighthouse_AWS_{station}/lighthouse_AWS_{station}_Table_{resolution}.dat"
    outfile_nc = f"{paths['local_storage']}lighthouse_AWS_{station}/{resolution}/lighthouse_AWS_{station}_Table_{resolution}_{from_time.year}{from_time.month:02d}{from_time.day:02d}.nc"
    backup_nc = f"{paths['harddrive']}lighthouse_AWS_{station}/{resolution}/lighthouse_AWS_{station}_Table_{resolution}_{from_time.year}{from_time.month:02d}{from_time.day:02d}.nc"

    print("restructuring {r} data from AWS {s}...".format(r=resolution, s=station))

    col_names = pd.read_csv(infile, header=1, sep=",", nrows=1).to_dict('records')[0]

    time_avail = pd.to_datetime(pd.read_csv(infile, header=3, sep=",", usecols=[0], na_values="NAN").iloc[:,0]).dt.to_pydatetime()
    ind = np.where((time_avail >= from_time) & (time_avail < to_time))[0]
    
    if len(ind) == 0:
        return time_avail[-1]
    
    data = pd.read_csv(infile, header=3+ind[0], nrows=len(ind), sep=",", na_values="NAN", names=list(col_names.keys()))
    
    if ((station == "1887") & (from_time < datetime.datetime(2022,11,6))):
        data['wind_direction_Avg'] -= 65.
    
    if from_time < datetime.datetime(2022,5,8):
        new_names = {'TIMESTAMP': "time",
                     'BP_mbar_Avg': "air_pressure",
                     'RH_Avg': "relative_humidity",
                     'AirT_C': "temperature",
                     'WS_ms_Max': "wind_speed_Max",
                     'WS_ms_S_WVT': "wind_speed",
                     'WindDir_D1_WVT': "wind_direction",
                     'WindDir_SD1_WVT': "wind_direction_Std",
                     'MetSENS_Status': "sensor_status"}
        cols_to_drop = [i for i in list(data.columns) if i not in list(new_names.keys())]
        data.drop(columns=cols_to_drop, inplace=True)
        data.rename(new_names, axis=1, inplace=True)

        
        variable_attributes = {
        "units": {'time': "seconds since 1970-01-01T00:00:00Z",
                  'air_pressure': "hPa", 'relative_humidity': "percent",
                  'temperature': "degree_C", 
                  'wind_speed': "m s-1", 'wind_speed_Max': "m s-1", 'wind_direction': "degree",
                  'wind_direction_Std': "degree"},
        
        "long_name": {'time': "UTC time",
                      'air_pressure': "air pressure averaged over the sampling interval",
                      'relative_humidity': "air relative humidity averaged over the sampling interval",
                      'temperature': "air temperature sampled at the respective timestamp",
                      'wind_speed_Max': "maximum wind speed during the sampling interval",
                      'wind_speed': "wind speed averaged over the sampling interval",
                      'wind_direction': "wind direction averaged over the sampling interval",
                      'wind_direction_Std': "standard deviation of the wind direction during the sampling interval",
                      'sensor_status': "sensor status"},
        "standard_name": {'time': "time", 'air_pressure': "air_pressure",
                          'relative_humidity': "relative_humidity",
                          'temperature': "air_temperature",
                          'wind_speed': "wind_speed", 'wind_speed_Max': "wind_speed", 'wind_direction': "wind_from_direction",
                          'sensor_status': "status_flag"}}
        
    else:
        new_names = {"TIMESTAMP": "time",
                     'air_pressure_Avg': "air_pressure",
                     'relative_humidity_Avg': "relative_humidity",
                     'temperature_Avg': "temperature",
                     'wind_speed_Max': "wind_speed_Max",
                     'wind_speed_Avg': "wind_speed",
                     'wind_direction_Avg': "wind_direction",
                     'wind_direction_Std': "wind_direction_Std",
                     'MetSENS_Status': 'sensor_status'}
        cols_to_drop = [i for i in list(data.columns) if i not in list(new_names.keys())]
        data.drop(columns=cols_to_drop, inplace=True)
        data.rename(new_names, axis=1, inplace=True)

        
        variable_attributes = {
        "units": {'time': "seconds since 1970-01-01T00:00:00Z", 'wind_speed': "m s-1", 'wind_direction': "degree", 'wind_direction_Std': "degree",
                  'wind_speed_Max': "m s-1", 'air_pressure': "hPa",
                  'relative_humidity': "percent", 'temperature': "degree_C"},
        
        "long_name": {'time': "UTC time", 'wind_speed': "wind speed averaged over the sampling interval", "sensor_status": "sensor status",
                      'wind_direction': "wind direction averaged over the sampling interval",
                      'wind_direction_Std': "standard deviation of the wind direction during the sampling interval",
                      'wind_speed_Max': "maximum wind speed during the sampling interval",
                      'air_pressure': "air pressure averaged over the sampling interval",
                      'relative_humidity': "air relative humidity averaged over the sampling interval",
                      'temperature': "air temperature averaged over the sampling interval"},
        
        "standard_name": {'time': "time", 'wind_speed': "wind_speed", 'wind_direction': "wind_from_direction",
                  'wind_speed_Max': "wind_speed", 'air_pressure': "air_pressure", "sensor_status": "status_flag",
                  'relative_humidity': "relative_humidity", 'temperature': "air_temperature"}}
        



    # transfer timestamps into Python datetime objects
    data["time"] = [dt.replace(tzinfo=datetime.timezone.utc).timestamp() for dt in pd.to_datetime(data["time"]).dt.to_pydatetime()]

    data.set_index("time", inplace=True)
    
    data = data[~data.index.duplicated(keep='first')]
    
    start_data_coverage = datetime.datetime.fromtimestamp(data.index[0], datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    end_data_coverage = datetime.datetime.fromtimestamp(data.index[-1], datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    
    # filter data for unphysical outliers
    for vari in data.columns:
        if "wind_direction" in vari:
            data.loc[data[vari]<0., vari] = np.nan
            data.loc[data[vari]>360., vari] = np.nan
        elif "wind_speed" in vari:
            data.loc[data[vari]<0., vari] = np.nan
            data.loc[data[vari]>80., vari] = np.nan
        elif vari == 'temperature':
            data.loc[data[vari]<-80., vari] = np.nan
            data.loc[data[vari]>40., vari] = np.nan
        elif vari == "relative_humidity":
            data.loc[data[vari]<0., vari] = np.nan
            data.loc[data[vari]>100., vari] = np.nan
        elif vari == "air_pressure":
            data.loc[data[vari]<800., vari] = np.nan
            data.loc[data[vari]>1100., vari] = np.nan
    
    data.fillna(-99999., inplace=True)
    
    ds = data.to_xarray()
    for vari in list(ds.variables):
        if vari == "time":
            ds[vari] = ds[vari].astype('float64', keep_attrs=True)
        elif vari != 'sensor_status':
            ds[vari] = ds[vari].astype('float32', keep_attrs=True)
        attris = {}
        for a, d in variable_attributes.items():
            try:
                attris[a] = d[vari]
            except KeyError:
                pass
        ds[vari].attrs = attris
        
        
    # Assign global attributes
    dtnow = datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")
    
    lighthouses = {"1884": {"name": "Narveneset", 'lat': 78.56343,'lon': 16.29687, "installation": datetime.datetime(2022,6,18,0,0,0,tzinfo=datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")},
                   "1885": {"name": "Bohemanneset", 'lat': 78.38166, 'lon': 14.75300, "installation": datetime.datetime(2021,8,20,0,0,0,tzinfo=datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")},
                   "1886": {"name": "Daudmannsodden", 'lat': 78.21056,'lon': 12.98685, "installation": datetime.datetime(2022,7,8,0,0,0,tzinfo=datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")},
                   "1887": {"name": "Gasoyane", 'lat': 78.45792,'lon': 16.20082, "installation": datetime.datetime(2022,9,3,0,0,0,tzinfo=datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")}}
    
    
    if from_time < datetime.datetime(2022,1,1):
        creator_insts = "The University Centre in Svalbard; The University Centre in Svalbard; The University Centre in Svalbard; MET Norway"
        creator_names = "L. Frank; F. Schalamon; M. O. Jonassen; T. Remes"
        creator_mails = "lukasf@unis.no; florina.schalamon@uni-graz.at; mariusj@unis.no; teresav@met.no"
        creator_urls = "https://orcid.org/0000-0003-1472-7967; https://orcid.org/0000-0002-2509-4133; https://orcid.org/0000-0002-4745-9009; https://orcid.org/0000-0002-6421-859X"
    elif from_time < datetime.datetime(2022,10,1):
        creator_insts = "The University Centre in Svalbard; The University Centre in Svalbard; The University Centre in Svalbard; The University Centre in Svalbard; MET Norway"
        creator_names = "L. Frank; F. Schalamon; A. Stenlund; M. O. Jonassen; T. Remes"
        creator_mails = "lukasf@unis.no; florina.schalamon@uni-graz.at; ;mariusj@unis.no; teresav@met.no"
        creator_urls = "https://orcid.org/0000-0003-1472-7967; https://orcid.org/0000-0002-2509-4133; https://orcid.org/0000-0003-4241-735X; https://orcid.org/0000-0002-4745-9009; https://orcid.org/0000-0002-6421-859X"
    else:
        creator_insts = "The University Centre in Svalbard; The University Centre in Svalbard; MET Norway"
        creator_names = "L. Frank; M. O. Jonassen; T. Remes"
        creator_mails = "lukasf@unis.no; mariusj@unis.no; teresav@met.no"
        creator_urls = "https://orcid.org/0000-0003-1472-7967; https://orcid.org/0000-0002-4745-9009; https://orcid.org/0000-0002-6421-859X"
    
    
    
    global_attributes = {"location": lighthouses[station]["name"],
        "latitude": lighthouses[station]["lat"],
        "longitude": lighthouses[station]["lon"],
        "title": f"Standard meteorological near-surface observations from {from_time.strftime('%Y-%m-%d')} at {lighthouses[station]['name']} in Isfjorden, Svalbard.",
        "summary": "The file contains time series of the standard meteorological near-surface parameters temperature, humidity, pressure, wind speed and wind direction. The raw data is only filtered for obviously wrong values, otherwise the data is made available as is.",
        "keywords": "EARTH SCIENCE>ATMOSPHERE>ATMOSPHERIC TEMPERATURE>SURFACE TEMPERATURE>AIR TEMPERATURE, EARTH SCIENCE>ATMOSPHERE>ATMOSPHERIC TEMPERATURE>SURFACE TEMPERATURE>DEW POINT TEMPERATURE, EARTH SCIENCE>ATMOSPHERE>ATMOSPHERIC WATER VAPOR>WATER VAPOR INDICATORS>HUMIDITY>RELATIVE HUMIDITY, \
            EARTH SCIENCE>ATMOSPHERE>ATMOSPHERIC WINDS>SURFACE WINDS>WIND DIRECTION, EARTH SCIENCE>ATMOSPHERE>ATMOSPHERIC WINDS>SURFACE WINDS>WIND SPEED, EARTH SCIENCE>ATMOSPHERE>ATMOSPHERIC PRESSURE>SURFACE PRESSURE, \
            ACADEMIC>UNIS, GEOGRAPHIC REGION>POLAR, CONTINENT>EUROPE>NORTHERN EUROPE>SCANDINAVIA>NORWAY, \
            IN SITU/LABORATORY INSTRUMENTS>CURRENT/WIND METERS>SONIC ANEMOMETER, IN SITU/LABORATORY INSTRUMENTS>PRESSURE/HEIGHT METERS>PRESSURE SENSORS, \
            IN SITU/LABORATORY INSTRUMENTS>TEMPERATURE/HUMIDITY SENSORS>TEMPERATURE SENSORS, IN SITU/LABORATORY INSTRUMENTS>TEMPERATURE/HUMIDITY SENSORS>HUMIDITY SENSORS",
        "keywords_vocabulary": 'GCMD Science Keywords',
        "naming_authority": "ADC",
        "data_type": "netCDF-4",
        "feature_type": "Time Series",
        "geospatial_lat_min": lighthouses[station]["lat"],
        "geospatial_lat_max": lighthouses[station]["lat"],
        "geospatial_lon_min": lighthouses[station]["lon"],
        "geospatial_lon_max": lighthouses[station]["lon"],
        "area": "Svalbard, Isfjorden",
        "time_coverage_start": start_data_coverage,
        "time_coverage_end": end_data_coverage,
        "Conventions": "CF-1.6, ACDD-1.3",
        'date_created': dtnow,
        'history': f'File created at {dtnow} using xarray in Python3.',
        "processing_level": "Known bad data has been replaced with null values, ranges applied",
        "creator_type": "group",
        "creator_institution": creator_insts,
        "creator_name": creator_names,
        "creator_email": creator_mails,
        "creator_url": creator_urls,
        "institution": "The University Centre in Svalbard (UNIS)",
        "project": "IWIN: Isfjorden Weather Information Network",
        "source": "Campbell Scientific METSENS 500",
        "platform": "LAND-BASED PLATFORMS>PERMANENT LAND SITES>METEOROLOGICAL STATIONS, LAND-BASED PLATFORMS>PERMANENT LAND SITES>COASTAL STATIONS",
        "platform_vocabulary": 'GCMD Science Keywords',
        "instrument_type": "Campbell Scientific METSENS 500",
        "license": "https://creativecommons.org/licenses/by/4.0/",
        "iso_topic_category": "atmosphere",
        "principal_investigator": "L. Frank",
        "publisher_name": "MET Norway",
        "publisher_url": "https://thredds.met.no/thredds/catalog.html",
        "publisher_email": "thredds@met.no",
        "id": "???"}
    
    for a, value in global_attributes.items():
        ds.attrs[a] = value
    
    myencoding = {v: {'_FillValue': -99999., 'zlib': False} for v in list(ds.variables) if v not in ["time", "sensor_status"]}
    myencoding["time"] = {'_FillValue': None}
    
    ds.to_netcdf(outfile_nc, unlimited_dims=["time"], encoding=myencoding)
    
    if os.path.isfile(backup_nc):
        os.remove(backup_nc)
    shutil.copyfile(outfile_nc, backup_nc)


    return time_avail[-1]


def export_json(df, out_path):

    properties = list(df.drop(['lat', 'lon'], axis=1).columns)
    
    df = df.fillna(value=0)

    geojson = {'type':'FeatureCollection', 'features':[]}

    # loop through each row in the dataframe and convert each row to geojson format
    for _, row in df.iterrows():
        # create a feature template to fill in
        feature = {'type':'Feature',
                   'properties':{},
                   'geometry':{'type':'Point',
                               'coordinates':[row["lon"],row["lat"]]}}


        # for each column, get the value and add it as a new feature property
        for prop in properties:
            feature['properties'][prop] = row[prop]
        
        # add this feature (aka, converted dataframe row) to the list of features inside our dict
        geojson['features'].append(feature)
        
    with open(out_path, 'w', encoding='utf-8') as f:
        f.write(json.dumps(geojson))
    
    return 


def create_GIS_input_file(boats, lighthouses, met_stations, past_hours=3, path="C:/Users/unismet/OneDrive - Universitetssenteret på Svalbard AS/"):
    
    boat_names = {1883: "MS_Bard", 1872: "MS_Polargirl", 1924: "MS_Billefjord"}
    
    boat_temp_res = {1883: 10, 1872: 15, 1924: 15}
    
    lighthouse_names = {1884: {"name": "Narveneset", 'lat': 78.56343,'lon': 16.29687},
                   1885: {"name": "Bohemanneset", 'lat': 78.38166, 'lon': 14.75300},
                   1886: {"name": "Daudmannsodden", 'lat': 78.21056,'lon': 12.98685},
                   1887: {"name": "Gasoyane", 'lat': 78.45792,'lon': 16.20082}}
    
    station_names = {"IR": {"name": "Isfjord Radio", "ID": "SN99790", "height": 7., "lat": 78.0625, "lon": 13.6192},
                     "LYR": {"name": "Longyearbyen Airport", "ID": "SN99840", "height": 28., "lat": 78.2453, "lon": 15.5015},
                     "PYR": {"name": "Pyramiden", "ID": "SN99880", "height": 20., "lat": 78.6557, "lon": 16.3603},
                     "NS": {"name": "Nedre Sassendalen", "ID": "SN99882", "height": 13., "lat": 78.3313, "lon": 16.6818},
                     "PB": {"name": "Plataberget", "ID": "SN99843", "height": 450., "lat": 78.2278, "lon": 15.378},
                     "AD": {"name": "Adventdalen", "ID": "SN99870", "height": 15., "lat":  78.2022, "lon": 15.831},
                     "IH": {"name": "Innerhytta", "ID": "SN99879", "height": 81., "lat":  78.18883, "lon": 16.34423},
                     "AO": {"name": "Akseloya", "ID": "SN99765", "height": 20., "lat":  77.6873, "lon": 14.7578},
                     "KL": {"name": "Klauva", "ID": "SN99884", "height": 480., "lat":  78.3002, "lon": 18.2248},
                     "ID": {"name": "Istjorndalen", "ID": "SN99770", "height": 188., "lat":  78.0092, "lon": 15.2108}, 
                     "JH": {"name": "Janssonhagen", "ID": "SN99874", "height": 250., "lat":  78.18, "lon": 16.41},
                     "RP": {"name": "Reindalspasset", "ID": "SN99763", "height": 181., "lat":  78.0648, "lon":  17.0442},
                     "SV": {"name": "Svea", "ID": "SN99760", "height": 9., "lat":  77.8953, "lon": 16.72}}
    
    rounding_precision = {"temperature": 0, "relative_humidity": 0, "wind_direction": 0, "pressure": 0, "wind_speed": 1, "u": 1, "v": 1, "u_knts": 1, "v_knts": 1}
    
    list_all = []
    for b in boats.keys():
        ind = np.where(np.array([t.minute for t in boats[b].data["local_time"]]) % boat_temp_res[b] == 0)[0]
        df = pd.DataFrame([boat_names[b]]*len(ind), index=boats[b].data['local_time'][ind], columns=["STAT"])
        df["lat"] = boats[b].data["latitude"][ind]
        df["lon"] = boats[b].data["longitude"][ind]
        for v in ['temperature', 'air_pressure', 'relative_humidity', 'wind_speed', 'wind_direction']:
            df[v] = boats[b].data[v][ind]
        list_all.append(df)
    for l in lighthouses.keys():
        df = pd.DataFrame([lighthouse_names[l]["name"]], index=[lighthouses[l].data['local_time'][-1]], columns=["STAT"])
        df["lat"] = lighthouse_names[l]["lat"]
        df["lon"] = lighthouse_names[l]["lon"]
        for v in ['temperature', 'air_pressure', 'relative_humidity', 'wind_speed', 'wind_direction', 'wind_sector']:
            df[v] = [lighthouses[l].data[v][-1]]
        list_all.append(df)
    for m in met_stations.keys():
        df = met_stations[m].iloc[-1:].copy()
        df["STAT"] = station_names[m]["name"]
        df["lat"] = station_names[m]["lat"]
        df["lon"] = station_names[m]["lon"]
        list_all.append(df)
        
    df_total = pd.concat(list_all).drop(["u", "v", "u_knts", "v_knts", "wind_sector"], axis=1)
    
    df_total = df_total.round(rounding_precision)
    df_total = df_total[df_total.index >= (df_total.index.max() - pd.Timedelta(hours=past_hours))]
    
    df_total.index = df_total.index.strftime('%H:%M')
    
    df_total.rename({'temperature': "Temperature", 'air_pressure': "Pressure", 'relative_humidity': "Relative Humidity",
                     "wind_speed": "Wind Speed", "wind_direction": "Wind Direction"})
    
    df_total.to_csv(f"{path}IWIN/GIS/latest_isfjorden_data.csv")
    
    export_json(df_total, f"{path}IWIN/GIS/latest_isfjorden_data.geojson")
    
    return
    
