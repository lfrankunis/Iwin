# -*- coding: utf-8 -*-
"""
Created on Sat Mar 27 14:37:44 2021
@author: unismet
"""

import numpy as np
import pandas as pd
import datetime
from netCDF4 import Dataset


def restructure_mobile_AWS(from_time, to_time, station="1883", resolution="10min", path="C:/Data/"):
    """
    Function to restructure the mobile AWS data into e.g. daily files
    Parameters
    ----------
    path : str
        string defining the path where to find the input data and save the output data
    station : str
        string specifying the station ("1883" or "1872" or "1924")
    resolution : str
        string specifying the resolution of the data ("20sec", "1min", "5min", "10min" or "hour")
    from_time : datetime object
        start time of the restructured file
    to_time : datetime object
        end time of the restructured file
    Returns both an ascii file similar to the direct output files from Loggerent as well as a netCDF4 file
    -------
    In case the restructured interval is not one day, remember to change the output file names
    """
    ## data path for before september
    # infile = "{p}mobile_AWS_{s}/raw_backups/mobile_AWS_{s}_Table_{r}.dat".format(p=path, s=station, r=resolution)
    
    infile = "{p}mobile_AWS_{s}/mobile_AWS_{s}_Table_{r}.dat".format(p=path, s=station, r=resolution)
    outfile_ascii = "{p}mobile_AWS_{s}/{a}{b:02d}{c:02d}/ascii/{a}{b:02d}{c:02d}_mobile_AWS_{s}_Table_{r}.dat".format(p=path, s=station, r=resolution, a=from_time.year, b=from_time.month, c=from_time.day)
    outfile_nc = "{p}mobile_AWS_{s}/{a}{b:02d}{c:02d}/nc/{a}{b:02d}{c:02d}_mobile_AWS_{s}_Table_{r}.nc".format(p=path, s=station, r=resolution, a=from_time.year, b=from_time.month, c=from_time.day)

    print("restructuring {r} data from AWS {s}...".format(r=resolution, s=station))

    col_names = pd.read_csv(infile, header=1, sep=",", nrows=1).to_dict('records')[0]

    time_avail = pd.to_datetime(pd.read_csv(infile, header=3, sep=",", usecols=[0], na_values="NAN").iloc[:,0]).dt.to_pydatetime()

    ind = np.where((time_avail >= from_time) & (time_avail < to_time))[0]

    data = pd.read_csv(infile, header=3+ind[0], nrows=len(ind), sep=",", na_values="NAN", names=list(col_names.keys()))

    # transfer timestamps into Python datetime objects
    data["TIMESTAMP"] = pd.to_datetime(data["TIMESTAMP"]).dt.to_pydatetime()
    data["Wind_Speed_corrected_TMx"] = pd.to_datetime(data["Wind_Speed_corrected_TMx"]).dt.to_pydatetime()
    data["Wind_Speed_raw_TMx"] = pd.to_datetime(data["Wind_Speed_raw_TMx"]).dt.to_pydatetime()

    # extract lat, lon and lat from GPS Location
    latitude = np.ones((len(data["TIMESTAMP"]))) * np.nan
    longitude = np.ones((len(data["TIMESTAMP"]))) * np.nan
    altitude = np.ones((len(data["TIMESTAMP"]))) * np.nan
    for c, gps in enumerate(data["GPS_Location"]):
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


    # write data to new text file containing only that one day
    # read header lines from complete data file
    with open(infile, 'r') as f:
        header = [next(f) for x in range(4)]


    # header lines
    with open(outfile_ascii, 'w') as f:
        for line in header:
            f.write("{a}".format(a=line))
    # data lines
    data.to_csv(outfile_ascii, mode="a", header=False, float_format='%.5f', na_rep="NAN", index=False)

    # write data also into nc file
    with Dataset(outfile_nc, 'w', format="NETCDF4") as f:
        f.Comments = "UNIS {r} mobile AWS {s} data from {a} until {b}.".format(r=resolution, s=station, a=from_time, b=to_time)
        f.createDimension('time', len(ind))

        var = f.createVariable('temperature', 'f8', ('time',))
        var.Longname = "Ambient_Temperature"
        var.Unit = col_names["Ambient_Temperature"]
        var[:] = data["Ambient_Temperature"]

        var = f.createVariable('temperature_max', 'f8', ('time',))
        var.Longname = "Ambient_Temperature_Max"
        var.Unit = col_names["Ambient_Temperature_Max"]
        var[:] = data["Ambient_Temperature_Max"]

        var = f.createVariable('temperature_min', 'f8', ('time',))
        var.Longname = "Ambient_Temperature_Min"
        var.Unit = col_names["Ambient_Temperature_Min"]
        var[:] = data["Ambient_Temperature_Min"]

        var = f.createVariable('pressure', 'f8', ('time',))
        var.Longname = "Barometric_Pressure"
        var.Unit = col_names["Barometric_Pressure"]
        var[:] = data["Barometric_Pressure"]

        var = f.createVariable('relative_humidity', 'f8', ('time',))
        var.Longname = "Relative_Humidity"
        var.Unit = col_names["Relative_Humidity"]
        var[:] = data["Relative_Humidity"]

        var = f.createVariable('relative_humidity_max', 'f8', ('time',))
        var.Longname = "Relative_Humidity_Max"
        var.Unit = col_names["Relative_Humidity_Max"]
        var[:] = data["Relative_Humidity_Max"]

        var = f.createVariable('relative_humidity_min', 'f8', ('time',))
        var.Longname = "Relative_Humidity_Min"
        var.Unit = col_names["Relative_Humidity_Min"]
        var[:] = data["Relative_Humidity_Min"]

        var = f.createVariable('dewpoint', 'f8', ('time',))
        var.Longname = "Sensor_Dewpoint"
        var.Unit = col_names["Sensor_Dewpoint"]
        var[:] = data["Sensor_Dewpoint"]

        var = f.createVariable('wind_speed', 'f8', ('time',))
        var.Longname = "corrected Wind Speed using GPS data"
        var.Unit = col_names["Wind_Speed_Corrected_S_WVT"]
        var[:] = data["Wind_Speed_Corrected_S_WVT"]

        var = f.createVariable('wind_speed_max', 'f8', ('time',))
        var.Longname = "corrected Wind Speed Maximum using GPS data"
        var.Unit = col_names["Wind_Speed_corrected_Max"]
        var[:] = data["Wind_Speed_corrected_Max"]

        var = f.createVariable('wind_speed_max_timestamp', 'f8', ('time',))
        var.Longname = "timestamp of occurrence of the corrected Wind Speed Maximum"
        var.Unit = "Seconds since 1.1.1970"
        var[:] = np.array([dt.replace(tzinfo=datetime.timezone.utc).timestamp() for dt in data["Wind_Speed_corrected_TMx"]])

        var = f.createVariable('wind_speed_raw', 'f8', ('time',))
        var.Longname = "uncorrected Wind Speed"
        var.Unit = col_names["Wind_Speed_S_WVT"]
        var[:] = data["Wind_Speed_S_WVT"]

        var = f.createVariable('wind_speed_max_raw', 'f8', ('time',))
        var.Longname = "uncorrected Wind Speed Maximum"
        var.Unit = col_names["Wind_Speed_raw_Max"]
        var[:] = data["Wind_Speed_raw_Max"]

        var = f.createVariable('wind_speed_max_raw_timestamp', 'f8', ('time',))
        var.Longname = "timestamp of occurrence of the uncorrected Wind Speed Maximum"
        var.Unit = "Seconds since 1.1.1970"
        var[:] = np.array([dt.replace(tzinfo=datetime.timezone.utc).timestamp() for dt in data["Wind_Speed_raw_TMx"]])

        var = f.createVariable('wind_direction', 'f8', ('time',))
        var.Longname = "corrected Wind Direction using compass data"
        var.Unit = col_names["Wind_Direction_Corrected_D1_WVT"]
        var[:] = data["Wind_Direction_Corrected_D1_WVT"]

        var = f.createVariable('wind_direction_std', 'f8', ('time',))
        var.Longname = "corrected Wind Direction Standard Deviation"
        var.Unit = col_names["Wind_Direction_Corrected_SDI_WVT"]
        var[:] = data["Wind_Direction_Corrected_SDI_WVT"]

        var = f.createVariable('wind_direction_raw', 'f8', ('time',))
        var.Longname = "uncorrected Wind Direction"
        var.Unit = col_names["Wind_Direction_D1_WVT"]
        var[:] = data["Wind_Direction_D1_WVT"]

        var = f.createVariable('wind_direction_std_raw', 'f8', ('time',))
        var.Longname = "uncorrected Wind Direction Standard Deviation"
        var.Unit = col_names["Wind_Direction_SDI_WVT"]
        var[:] = data["Wind_Direction_SDI_WVT"]

        var = f.createVariable('latitude', 'f8', ('time',))
        var.Longname = "GPS Latitude"
        var.Unit = "degN"
        var[:] = latitude

        var = f.createVariable('longitude', 'f8', ('time',))
        var.Longname = "GPS Longitude"
        var.Unit = "degE"
        var[:] = longitude

        var = f.createVariable('altitude', 'f8', ('time',))
        var.Longname = "GPS Altitude"
        var.Unit = "m"
        var[:] = altitude

        var = f.createVariable('time', 'f8', ('time',))
        var.Longname = "time"
        var.Unit = "Seconds since 1.1.1970"
        var[:] = np.array([dt.replace(tzinfo=datetime.timezone.utc).timestamp() for dt in data["TIMESTAMP"]])

        if resolution == "hour":

            data["Sensor_Supply_volt_TMn"] = pd.to_datetime(data["Sensor_Supply_volt_TMn"]).dt.to_pydatetime()

            var = f.createVariable('BattV', 'f8', ('time',))
            var.Longname = "Battery_Voltage"
            var.Unit = col_names["BattV"]
            var[:] = data["BattV"]

            var = f.createVariable('PTemp_C', 'f8', ('time',))
            var.Longname = "Internal Logger_Panel_Temperature"
            var.Unit = col_names["PTemp_C"]
            var[:] = data["PTemp_C"]

            var = f.createVariable('sonic_status', 'i4', ('time',))
            var.Longname = "Sonic_Status_Code"
            var.Unit = col_names["Sensor_status"]
            var[:] = data["Sensor_status"]

            var = f.createVariable('SonicV_min', 'f8', ('time',))
            var.Longname = "minimal_Sonic_Supply_Voltage"
            var.Unit = col_names["Sensor_Supply_volt_Min"]
            var[:] = data["Sensor_Supply_volt_Min"]

            var = f.createVariable('SonicV_min_timestamp', 'f8', ('time',))
            var.Longname = "timestamp_of_minimal_Sonic_Supply_Voltage"
            var.Unit = col_names["Sensor_Supply_volt_TMn"]
            var[:] = data["Sensor_Supply_volt_TMn"]

    return




def restructure_lighthouse_AWS(from_time, to_time, station="1885", resolution="10min", path="C:/Data/"):
    """
    Function to restructure the lighthouse AWS data into e.g. daily files
    Parameters
    ----------
    path : str
        string defining the path where to find the input data and save the output data
    station : str
        string specifying the station (e.g. "1885")
    resolution : str
        string specifying the resolution of the data ("1min", "10min" or "hour")
    from_time : datetime object
        start time of the restructured file
    to_time : datetime object
        end time of the restructured file
    Returns both an ascii file similar to the direct output files from Loggerent as well as a netCDF4 file
    -------
    In case the restructured interval is not one day, remember to change the output file names
    """

    infile = "{p}lighthouse_AWS_{s}/lighthouse_AWS_{s}_Table_{r}.dat".format(p=path, s=station, r=resolution)
    outfile_ascii = "{p}lighthouse_AWS_{s}/{a}{b:02d}{c:02d}/ascii/{a}{b:02d}{c:02d}_lighthouse_AWS_{s}_Table_{r}.dat".format(p=path, s=station, r=resolution, a=from_time.year, b=from_time.month, c=from_time.day)
    outfile_nc = "{p}lighthouse_AWS_{s}/{a}{b:02d}{c:02d}/nc/{a}{b:02d}{c:02d}_lighthouse_AWS_{s}_Table_{r}.nc".format(p=path, s=station, r=resolution, a=from_time.year, b=from_time.month, c=from_time.day)

    print("restructuring {r} data from AWS {s}...".format(r=resolution, s=station))

    col_names = pd.read_csv(infile, header=1, sep=",", nrows=1).to_dict('records')[0]

    time_avail = pd.to_datetime(pd.read_csv(infile, header=3, sep=",", usecols=[0], na_values="NAN").iloc[:,0]).dt.to_pydatetime()

    ind = np.where((time_avail >= from_time) & (time_avail < to_time))[0]

    data = pd.read_csv(infile, header=3+ind[0], nrows=len(ind), sep=",", na_values="NAN", names=list(col_names.keys()))

    # transfer timestamps into Python datetime objects
    data["TIMESTAMP"] = pd.to_datetime(data["TIMESTAMP"]).dt.to_pydatetime()

    # write data to new text file containing only that one day
    # read header lines from complete data file
    with open(infile, 'r') as f:
        header = [next(f) for x in range(4)]


    # header lines
    with open(outfile_ascii, 'w') as f:
        for line in header:
            f.write("{a}".format(a=line))
    # data lines
    data.to_csv(outfile_ascii, mode="a", header=False, float_format='%.5f', na_rep="NAN", index=False)

    # write data also into nc file
    with Dataset(outfile_nc, 'w', format="NETCDF4") as f:
        f.Comments = "UNIS {r} lighthouse AWS {s} data from {a} until {b}.".format(r=resolution, s=station, a=from_time, b=to_time)
        f.createDimension('time', len(ind))

        var = f.createVariable('temperature', 'f8', ('time',))
        var.Longname = "Ambient_Temperature"
        var.Unit = col_names["AirT_C"]
        var[:] = data["AirT_C"]

        var = f.createVariable('temperature_max', 'f8', ('time',))
        var.Longname = "Maximum_Ambient_Temperature"
        var.Unit = col_names["AirT_C_Max"]
        var[:] = data["AirT_C_Max"]

        var = f.createVariable('temperature_min', 'f8', ('time',))
        var.Longname = "Minimum_Ambient_Temperature"
        var.Unit = col_names["AirT_C_Min"]
        var[:] = data["AirT_C_Min"]

        var = f.createVariable('pressure', 'f8', ('time',))
        var.Longname = "Barometric_Pressure"
        var.Unit = col_names["BP_mbar"]
        var[:] = data["BP_mbar"]

        var = f.createVariable('pressure_avg', 'f8', ('time',))
        var.Longname = "Average_Barometric_Pressure"
        var.Unit = col_names["BP_mbar_Avg"]
        var[:] = data["BP_mbar_Avg"]

        var = f.createVariable('relative_humidity', 'f8', ('time',))
        var.Longname = "Relative_Humidity"
        var.Unit = col_names["RH"]
        var[:] = data["RH"]

        var = f.createVariable('relative_humidity_avg', 'f8', ('time',))
        var.Longname = "Average_Relative_Humidity"
        var.Unit = col_names["RH_Avg"]
        var[:] = data["RH_Avg"]

        var = f.createVariable('dewpoint', 'f8', ('time',))
        var.Longname = "Sensor_Dewpoint"
        var.Unit = col_names["DP_C"]
        var[:] = data["DP_C"]

        var = f.createVariable('wind_speed', 'f8', ('time',))
        var.Longname = "Wind_Speed"
        var.Unit = col_names["WS_ms"]
        var[:] = data["WS_ms"]

        var = f.createVariable('wind_speed_avg', 'f8', ('time',))
        var.Longname = "Average_Wind_Speed"
        var.Unit = col_names["WS_ms_Avg"]
        var[:] = data["WS_ms_Avg"]

        var = f.createVariable('wind_speed_min', 'f8', ('time',))
        var.Longname = "Minimum_Wind_Speed"
        var.Unit = col_names["WS_ms_Min"]
        var[:] = data["WS_ms_Min"]

        var = f.createVariable('wind_speed_max', 'f8', ('time',))
        var.Longname = "Maximum_Wind_Speed"
        var.Unit = col_names["WS_ms_Max"]
        var[:] = data["WS_ms_Max"]

        var = f.createVariable('wind_speed_vector_avg', 'f8', ('time',))
        var.Longname = "Vector_Average_Wind_Speed"
        var.Unit = col_names["WS_ms_S_WVT"]
        var[:] = data["WS_ms_S_WVT"]

        var = f.createVariable('wind_direction', 'f8', ('time',))
        var.Longname = "Wind_Direction"
        var.Unit = col_names["WindDir"]
        var[:] = data["WindDir"]

        var = f.createVariable('wind_direction_avg', 'f8', ('time',))
        var.Longname = "Average_Wind_Direction"
        var.Unit = col_names["WindDir_D1_WVT"]
        var[:] = data["WindDir_D1_WVT"]

        var = f.createVariable('wind_direction_std', 'f8', ('time',))
        var.Longname = "Standard_Deviation_Wind_Direction"
        var.Unit = col_names["WindDir_SD1_WVT"]
        var[:] = data["WindDir_SD1_WVT"]

        var = f.createVariable('time', 'f8', ('time',))
        var.Longname = "time"
        var.Unit = "Seconds since 1.1.1970"
        var[:] = np.array([dt.replace(tzinfo=datetime.timezone.utc).timestamp() for dt in data["TIMESTAMP"]])


        var = f.createVariable('BattV', 'f8', ('time',))
        var.Longname = "Battery_Voltage"
        var.Unit = col_names["BattV_Min"]
        var[:] = data["BattV_Min"]

        var = f.createVariable('PTemp_C', 'f8', ('time',))
        var.Longname = "Internal Logger_Panel_Temperature"
        var.Unit = col_names["PTemp_C"]
        var[:] = data["PTemp_C"]

        var = f.createVariable('sensor_status', 'S1', ('time',))
        var.Longname = "Sensor_Status_Code"
        var.Unit = col_names["MetSENS_Status"]
        var[:] = data["MetSENS_Status"]

    return 
