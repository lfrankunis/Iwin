#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 07:27:19 2023

@author: lukasf
"""


import sys
import pandas as pd
import datetime
import xarray as xr


#%% input

path_lighthouse_data = "https://thredds.met.no/thredds/dodsC/met.no/observations/unis/lighthouse_AWS"
path_mobile_data = "https://thredds.met.no/thredds/dodsC/met.no/observations/unis/mobile_AWS"
path_out = "/Users/lukasf/OneDrive - Universitetssenteret på Svalbard AS/IWIN/Paper/defined_dataset/"

#%% read data

data = {}
mobile_stations = ["MSBard", "MSBerg", "MSPolargirl", "MSBillefjord"]
for s in mobile_stations:
    print(s)
    with xr.open_dataset(f"{path_mobile_data}_{s}_1min") as ds:
        data[s] = ds.sortby("time").sel(time=slice("2021-01-01T00:00:00", "2023-06-23T00:00:00"))
    del data[s].attrs["publisher_name"]
    del data[s].attrs["publisher_institution"]
    del data[s].attrs["publisher_url"]
    del data[s].attrs["publisher_email"]
    del data[s].attrs["publisher_type"]
    
    del data[s].attrs["title"]
    data[s].attrs["title"] = f"Standard meteorological near-surface observations measured onboard {s} in Isfjorden, Svalbard."
    
    start_data_coverage = datetime.datetime.fromtimestamp(data[s].time[0].values.item()*1.e-9, datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    del data[s].attrs["time_coverage_start"]
    data[s].attrs["time_coverage_start"] = start_data_coverage
    
    del data[s].attrs["time_coverage_end"]
    data[s].attrs["time_coverage_end"] = "2023-06-23T00:00:00Z"
    
    data[s].to_netcdf(f"{path_out}mobile_AWS_{s}_1min.nc", unlimited_dims=["time"])

lighthouses = ["Daudmannsodden", "Bohemanneset", "Gasoyane", "Narveneset"]

for s in lighthouses:
    print(s)
    with xr.open_dataset(f"{path_lighthouse_data}_{s}_1min") as ds:
        data[s] = ds.sortby("time").sel(time=slice("2021-01-01T00:00:00", "2023-06-23T00:00:00"))
    del data[s].attrs["publisher_name"]
    del data[s].attrs["publisher_institution"]
    del data[s].attrs["publisher_url"]
    del data[s].attrs["publisher_email"]
    del data[s].attrs["publisher_type"]
    
    del data[s].attrs["title"]
    data[s].attrs["title"] = f"Standard meteorological near-surface observations at {s} in Isfjorden, Svalbard."
    
    start_data_coverage = datetime.datetime.fromtimestamp(data[s].time[0].values.item()*1.e-9, datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    del data[s].attrs["time_coverage_start"]
    data[s].attrs["time_coverage_start"] = start_data_coverage
    
    del data[s].attrs["time_coverage_end"]
    data[s].attrs["time_coverage_end"] = "2023-06-23T00:00:00Z"
    
    data[s].to_netcdf(f"{path_out}lighthouse_AWS_{s}_1min.nc", unlimited_dims=["time"])







