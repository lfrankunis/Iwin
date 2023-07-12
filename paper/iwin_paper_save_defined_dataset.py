#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 07:27:19 2023

@author: lukasf
"""


import sys
import pandas as pd
import pandas as pd
import cmocean as cmo
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition
import numpy as np
import geopandas as gpd
from pyproj import Proj
import matplotlib as mpl
import xarray as xr
import rioxarray as rxr
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import latex_helpers
import windrose
from scipy import stats
from scalebar import scale_bar

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
        data[s] = ds.sortby("time").sel(time=slice("2021-01-01T00:00:00", "2023-06-23T00:00:00")).to_netcdf(f"{path_out}mobile_AWS_{s}_1min.nc", unlimited_dims=["time"])

lighthouses = ["Daudmannsodden", "Bohemanneset", "Gasoyane", "Narveneset"]

for s in lighthouses:
    print(s)
    with xr.open_dataset(f"{path_lighthouse_data}_{s}_1min") as ds:
        data[s] = ds.sortby("time").sel(time=slice("2021-01-01T00:00:00", "2023-06-23T00:00:00")).to_netcdf(f"{path_out}lighthouse_AWS_{s}_1min.nc", unlimited_dims=["time"])

