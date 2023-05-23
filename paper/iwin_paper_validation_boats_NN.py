#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 17:05:34 2023

@author: lukasf
"""

import sys
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


lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()

mpl.rcParams.update({
    "pgf.texsystem": "lualatex",
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
    "axes.labelsize": 10,
    "font.size": 10,
    "legend.fontsize": 10,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10})

#%% input

path_iwin_data = "https://thredds.met.no/thredds/dodsC/met.no/observations/unis/mobile_AWS"
path_narveneset = "https://thredds.met.no/thredds/dodsC/met.no/observations/unis/lighthouse_AWS_Narveneset_1min"

path_out = "/Users/lukasf/Desktop/Iwin_paper_figures/iwin_paper_fig_eval_boats_NN.pdf"


#%% constants etc

station_labels = {"MSBard": "MS Bard", "MSPolargirl": "MS Polargirl", "MSBillefjord": "MS Billefjord"}

vari_labels = {"temperature": "temperature [°C]", "wind_speed_corrected": "wind speed [m/s]", "air_pressure": "pressure [hPa]", "relative_humidity": "rel. humidity [\%]", "wind_direction_corrected": "wind direction [deg]"}

nn_lat = 78.56343
nn_lon = 16.29687
d_lat = 0.02
d_lon = 0.08


with xr.open_dataset(path_narveneset) as ds:
    df_narve = ds.to_dataframe()
    df_narve.set_index(ds.time.values, inplace=True)
    df_narve = df_narve["2022-06-17":"2022-08-30"]
    df_narve = df_narve[["temperature", "relative_humidity", "air_pressure", "wind_speed", "wind_direction"]]
    df_narve.rename({'wind_speed': "wind_speed_corrected", "wind_direction": "wind_direction_corrected"}, axis=1, inplace=True)    
    # e = 0.01*df["relative_humidity"]*(611.2 * np.exp((17.62*df["temperature"])/(243.12+df["temperature"])))
    # df['specific_humidity'] = 1000.*(0.622*e)/(100.*df["air_pressure"]-0.378*e)
    # df.drop(["relative_humidity"], axis=1, inplace=True)
    df_narve["air_pressure"] += 7.*(1./8.)



mobile_stations = ["MSBard", "MSPolargirl", "MSBillefjord"]
boat_data = {}
for s in mobile_stations:
    print(s)
    with xr.open_dataset(f"{path_iwin_data}_{s}_1min") as ds:
        boat_data[s] = ds.sel(time=slice("2022-01-01", "2023-01-01")).load()


for b, b_data in boat_data.items():
    mask = ((b_data["latitude"] > nn_lat-d_lat) & (b_data["latitude"] < nn_lat+d_lat) & (b_data["longitude"] > nn_lon-d_lon) & (b_data["longitude"] < nn_lon+d_lon))    # Narveneset + max 3 km
    boat_data[b] = xr.where(mask, b_data, np.nan)

    
for b, b_data in boat_data.items():
    df = b_data.to_dataframe()
    df = df.set_index(b_data.time.values)
    df = df["2022-06-17":"2022-10-30"]
    # e = 0.01*df["relative_humidity"]*(611.2 * np.exp((17.62*df["temperature"])/(243.12+df["temperature"])))
    # df['specific_humidity'] = 1000.*(0.622*e)/(100.*df["air_pressure"]-0.378*e)
    df["air_pressure"] += 18.*(1./8.)
    df["air_pressure"][df["air_pressure"] < 970.] = np.nan
    boat_data[b] = df


#%% error statistics wrt Narveneset

biases = pd.DataFrame(index=list(boat_data.keys()))
maes = pd.DataFrame(index=list(boat_data.keys()))

for v, vari in enumerate(df_narve.columns):
    df = pd.DataFrame(df_narve[vari])
    df.rename({vari: "narveneset"}, axis=1, inplace=True)
    for s, station in enumerate(boat_data.keys()):
        df[station] = boat_data[station][vari]
        if vari == "wind_direction_corrected":
            df[station][df[station] > 180.] -= 360.
            df["narveneset"][df["narveneset"] > 180.] -= 360.
        df[station] -= df["narveneset"]
        
    biases[vari] = df.mean().drop("narveneset")
    maes[vari] = abs(df).mean().drop("narveneset")


biases.rename(station_labels, inplace=True)
maes.rename(station_labels, inplace=True)


#%% plots

fig, ax = plt.subplots(5,1, sharex=True, figsize=latex_helpers.set_size(503.6, whr=0.9))
for v, vari in enumerate(df_narve.columns):
    ax[v].plot(df_narve.index, df_narve[vari], color="k", lw=1., label="NN")
    for s, station in enumerate(boat_data.keys()):
        ax[v].plot(boat_data[station].index, boat_data[station][vari], color=f"C{s}", lw=1., label=station_labels[station])
        
ax[0].legend(loc="lower center", ncols=4, bbox_to_anchor=(0.5, 1.05))
ax[-1].set_xlim((df_narve.index[0], df_narve.index[-1]))

ax[1].set_ylim(top=100.)
ax[-1].set_ylim(bottom=0.)


for v, vari in enumerate(df_narve.columns):
    ax[v].grid()
    ax[v].set_ylabel(vari_labels[vari])

# plt.savefig(path_out, dpi=300)