#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 10:29:32 2023

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
path_met = "/Users/lukasf/OneDrive - Universitetssenteret på Svalbard AS/IWIN/reference_data_MET_stations/airport_summer.csv"
path_out = "/Users/lukasf/Desktop/Iwin_paper_figures/iwin_paper_fig_eval_boats.pdf"


#%% constants etc

station_labels = {"MSBard": "MS Bard", "MSPolargirl": "MS Polargirl", "MSBillefjord": "MS Billefjord"}

vari_labels = {"temperature": "temperature [°C]", "wind_speed_corrected": "wind speed [m/s]", "air_pressure": "pressure [hPa]", "relative_humidity": "rel. humidity [\%]"}


#%% read data

with open(path_met) as f:
    ncols = len(f.readline().split(';'))

df_lyr = pd.read_csv(path_met, dayfirst=True, sep=";", skipfooter=1, header=0, usecols=range(2,ncols+1), parse_dates=[0], decimal=",", na_values=["-"])

df_lyr["Time"] = df_lyr["Time(norwegian mean time)"] - pd.Timedelta("1h")
df_lyr.set_index("Time", inplace=True)
df_lyr.drop(["Time(norwegian mean time)"], axis=1, inplace=True)
df_lyr = df_lyr["2022-06-01":"2022-08-30"]
df_lyr.rename({'Air temperature': "temperature", 'Mean wind speed': "wind_speed_corrected", 'Relative air humidity': "relative_humidity",
       'Air pressure at station level': "air_pressure"}, axis=1, inplace=True)
# e = 0.01*df_lyr["relative_humidity"]*(611.2 * np.exp((17.62*df_lyr["temperature"])/(243.12+df_lyr["temperature"])))
# df_lyr['specific_humidity'] = 1000.*(0.622*e)/(100.*df_lyr["air_pressure"]-0.378*e)
# df_lyr.drop(["relative_humidity"], axis=1, inplace=True)
df_lyr = df_lyr[["temperature", "wind_speed_corrected", "relative_humidity", "air_pressure"]]
df_lyr["air_pressure"] += 28.*(1./8.)

all_data = {"airport": df_lyr}

mobile_stations = ["MSBard", "MSPolargirl", "MSBillefjord"]
boat_data = {}
for s in mobile_stations:
    print(s)
    with xr.open_dataset(f"{path_iwin_data}_{s}_1min") as ds:
        boat_data[s] = ds.sel(time=slice("2022-01-01", "2023-01-01")).load()


for b, b_data in boat_data.items():
    mask = ((b_data["latitude"] > 78.22745) & (b_data["latitude"] < 78.22878) & (b_data["longitude"] > 15.60521) & (b_data["longitude"] < 15.61387))
    boat_data[b] = xr.where(mask, b_data, np.nan)

    
    
for b, b_data in boat_data.items():
    df = b_data.to_dataframe()
    df = df.set_index(b_data.time.values)
    df = df["2022-06-01":"2022-08-30"]
    # e = 0.01*df["relative_humidity"]*(611.2 * np.exp((17.62*df["temperature"])/(243.12+df["temperature"])))
    # df['specific_humidity'] = 1000.*(0.622*e)/(100.*df["air_pressure"]-0.378*e)
    df["air_pressure"] += 18.*(1./8.)
    boat_data[b] = df
    
    all_data[b] = df
    
    
#%% error statistics wrt airport

biases = pd.DataFrame(index=list(boat_data.keys()))
maes = pd.DataFrame(index=list(boat_data.keys()))

for v, vari in enumerate(all_data["airport"].keys()):
    df = pd.DataFrame(df_lyr[vari])
    df.rename({vari: "airport"}, axis=1, inplace=True)
    for s, station in enumerate(boat_data.keys()):
        df[station] = boat_data[station][vari]
        df[station] -= df["airport"]
        
    biases[vari] = df.mean().drop("airport")
    maes[vari] = abs(df).mean().drop("airport")


biases.rename(station_labels, inplace=True)
maes.rename(station_labels, inplace=True)

#%% error statistics wrt to ensemble mean (only three boats)


biases_ens = pd.DataFrame(index=list(boat_data.keys()))
maes_ens = pd.DataFrame(index=list(boat_data.keys()))

dfs_for_hists = {}

for v, vari in enumerate(df_lyr.columns):
    df = pd.DataFrame(index=boat_data["MSBard"].index)
    for s, station in enumerate(boat_data.keys()):
        df[station] = boat_data[station][vari]
    dfs_for_hists[vari] = df.copy()
    df["ensemble_mean"] = df.mean(axis=1)
    for s, station in enumerate(boat_data.keys()):
        df[station] -= df["ensemble_mean"]
        
    biases_ens[vari] = df.mean().drop("ensemble_mean")
    maes_ens[vari] = abs(df).mean().drop("ensemble_mean")

biases_ens.rename(station_labels, inplace=True)
maes_ens.rename(station_labels, inplace=True)
    
  
#%% plots

fig, ax = plt.subplots(4,1, sharex=True, figsize=latex_helpers.set_size(503.6, whr=0.9))
for v, vari in enumerate(all_data["airport"].keys()):
    ax[v].plot(df_lyr.index, df_lyr[vari], color="k", lw=1., label="Svalbard Airport")
    for s, station in enumerate(boat_data.keys()):
        ax[v].plot(all_data[station].index, all_data[station][vari], color=f"C{s}", lw=1., label=station_labels[station])
        
ax[0].legend(loc="lower center", ncols=4, bbox_to_anchor=(0.5, 1.05))
ax[-1].set_xlim((all_data["airport"].index[0], all_data["airport"].index[-1]))

ax[1].set_ylim(bottom=0.)
# ax[-1].set_ylim(bottom=0.)


for v, vari in enumerate(all_data["airport"].keys()):
    ax[v].grid()
    ax[v].set_ylabel(vari_labels[vari])



#

fig, ax = plt.subplots(2,1, figsize=latex_helpers.set_size(503.6, whr=0.5), sharex=True)
ax[0].axhline(0., color="k", ls="--")
for v, vari in enumerate(all_data["airport"].keys()):
    ax[0].plot(np.arange(len(boat_data.keys())), biases[vari], color=f"C{v}", marker="x", label=" ".join(vari_labels[vari].split()[:-1]))

    ax[1].plot(np.arange(len(boat_data.keys())), maes[vari], color=f"C{v}", marker="x")

ax[0].legend(loc="lower center", ncols=4, bbox_to_anchor=(0.5, 1.05))


ax[0].set_ylabel("bias [K; m/s; \%; hPa]")
ax[1].set_ylabel("MAEs [K; m/s; \%; hPa]")

ax[1].set_ylim(bottom=0.)
ax[1].set_xticks(np.arange(len(boat_data.keys())))
ax[1].set_xticklabels(list(biases.index))

for a in ax:
    a.grid()

plt.show()

plt.savefig(f"{path_out}iwin_paper_fig_eval_boats_errorstats.pdf", dpi=300)



#

# fig, ax = plt.subplots(2,1, figsize=latex_helpers.set_size(503.6, whr=0.5), sharex=True)
# ax[0].axhline(0., color="k", ls="--")
# for v, vari in enumerate(all_data["airport"].keys()):
#     ax[0].plot(np.arange(len(boat_data.keys())), biases_ens[vari], color=f"C{v}", marker="x", label=" ".join(vari_labels[vari].split()[:-1]))

#     ax[1].plot(np.arange(len(boat_data.keys())), maes_ens[vari], color=f"C{v}", marker="x")

# ax[0].legend(loc="lower center", ncols=4, bbox_to_anchor=(0.5, 1.05))


# ax[0].set_ylabel("bias [K; m/s; \%; hPa]")
# ax[1].set_ylabel("MAEs [K; m/s; \%; hPa]")

# ax[1].set_ylim(bottom=0.)
# ax[1].set_xticks(np.arange(len(boat_data.keys())))
# ax[1].set_xticklabels(list(biases.index))

# for a in ax:
#     a.grid()

# plt.show()


