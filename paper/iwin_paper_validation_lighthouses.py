#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 10:56:08 2023

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

path_iwin_data = "https://thredds.met.no/thredds/dodsC/met.no/observations/unis/lighthouse_AWS"
path_met = "/Users/lukasf/OneDrive - Universitetssenteret på Svalbard AS/IWIN/reference_data_MET_stations/"
path_out = "/Users/lukasf/Desktop/Iwin_paper_figures/"


#%% constants etc.

station_heights = {"airport": 28., "nedre_sassendalen": 13., "isfjord_radio": 7.,
                   "Narveneset": 7., "Bohemanneset": 12., "Daudmannsodden": 39., "Gasoyane": 19.}

station_colors = {"airport": "grey", "nedre_sassendalen": "grey", "isfjord_radio": "grey",
                   "Narveneset": "r", "Bohemanneset": "r", "Daudmannsodden": "r", "Gasoyane": "r"}

station_labels = {"airport": "Svalbard\nAirport", "nedre_sassendalen": "Nedre\nSassendalen", "isfjord_radio": "Isfjord\nRadio",
                   "Narveneset": "Narveneset", "Bohemanneset": "Bohemanneset", "Daudmannsodden": "Daudmannsodden", "Gasoyane": "Gasoyane"}

vari_labels = {"temperature": "temperature [°C]", "wind_speed": "wind speed [m/s]", "air_pressure": "pressure [hPa]", "relative_humidity": "rel. humidity [\%]"}

#%% read data

met_stations = ["isfjord_radio", "airport", "nedre_sassendalen"]


with open(f"{path_met}airport.csv") as f:
    ncols = len(f.readline().split(';'))

met_data = {}
for m in met_stations:
    df = pd.read_csv(f"{path_met}{m}.csv", dayfirst=True, sep=";", skipfooter=1, header=0, usecols=range(2,ncols+1), parse_dates=[0], decimal=",", na_values=["-"])

    df["Time"] = df["Time(norwegian mean time)"] - pd.Timedelta("1h")
    df.set_index("Time", inplace=True)
    df.drop(["Time(norwegian mean time)"], axis=1, inplace=True)
    df = df["2022-09-03":"2023-05-30"]
    df.rename({'Air temperature': "temperature", 'Mean wind speed': "wind_speed", 'Relative air humidity': "relative_humidity",
           'Air pressure at station level': "air_pressure"}, axis=1, inplace=True)
    # e = 0.01*df["relative_humidity"]*(611.2 * np.exp((17.62*df["temperature"])/(243.12+df["temperature"])))
    # df['specific_humidity'] = 1000.*(0.622*e)/(100.*df["air_pressure"]-0.378*e)
    # df.drop(["relative_humidity"], axis=1, inplace=True)
    df = df[["temperature", "relative_humidity", "air_pressure", "wind_speed"]]
    met_data[m] = df

    

lighthouses = ["Daudmannsodden", "Bohemanneset", "Gasoyane", "Narveneset"]

lighthouse_data = {}
for s in lighthouses:
    print(s)
    with xr.open_dataset(f"{path_iwin_data}_{s}_1min") as ds:
        df = ds.to_dataframe()
        df.set_index(ds.time.values, inplace=True)
        df = df["2022-09-03":"2023-05-30"]
        # e = 0.01*df["relative_humidity"]*(611.2 * np.exp((17.62*df["temperature"])/(243.12+df["temperature"])))
        # df['specific_humidity'] = 1000.*(0.622*e)/(100.*df["air_pressure"]-0.378*e)
        # df.drop(["relative_humidity"], axis=1, inplace=True)
        lighthouse_data[s] = df
        
all_data = met_data | lighthouse_data

for s in all_data.keys():
    all_data[s]["air_pressure"] += station_heights[s]*(1./8.)
        


#%% error statistics wrt ensemble mean

biases = pd.DataFrame(index=list(all_data.keys()))
maes = pd.DataFrame(index=list(all_data.keys()))

for v, vari in enumerate(all_data["airport"].keys()):
    df = pd.DataFrame(index=all_data["airport"].index)
    for s, station in enumerate(all_data.keys()):
        df[station] = all_data[station][vari]
    df["ensemble_mean"] = df.mean(axis=1)
    for s, station in enumerate(all_data.keys()):
        df[station] -= df["ensemble_mean"]
        
    biases[vari] = df.mean().drop("ensemble_mean")
    maes[vari] = abs(df).mean().drop("ensemble_mean")


biases = biases.reindex(['Daudmannsodden', "isfjord_radio", "Bohemanneset", "airport", "Gasoyane", "Narveneset", "nedre_sassendalen"])
maes =  maes.reindex(['Daudmannsodden', "isfjord_radio", "Bohemanneset", "airport", "Gasoyane", "Narveneset", "nedre_sassendalen"])

biases.rename(station_labels, inplace=True)
maes.rename(station_labels, inplace=True)


#%% plots

fig, ax = plt.subplots(4,1, sharex=True, figsize=latex_helpers.set_size(503.6, whr=0.9))
for v, vari in enumerate(all_data["airport"].keys()):
    for s, station in enumerate(all_data.keys()):
        ax[v].plot(all_data[station].index, all_data[station][vari], color=f"C{s}", lw=1., label=station_labels[station])
        
ax[0].legend(loc="lower center", ncols=4, bbox_to_anchor=(0.5, 1.05))
ax[-1].set_xlim((all_data["airport"].index[0], all_data["airport"].index[-1]))

ax[1].set_ylim(top=100.)
ax[-1].set_ylim(bottom=0.)


for v, vari in enumerate(all_data["airport"].keys()):
    ax[v].grid()
    ax[v].set_ylabel(vari_labels[vari])

plt.savefig(f"{path_out}iwin_paper_fig_eval_lighthouses_timeseries.pdf", dpi=300)





#

fig, ax = plt.subplots(2,1, figsize=latex_helpers.set_size(503.6, whr=0.6), sharex=True)
ax[0].axhline(0., color="k", ls="--")
for v, vari in enumerate(all_data["airport"].keys()):
    ax[0].plot(np.arange(len(all_data.keys())), biases[vari], color=f"C{v}", label=" ".join(vari_labels[vari].split()[:-1]))
    ax[0].scatter([1.,3.,6.], biases.loc[[station_labels.get(s) for s in met_stations]][vari], s=20, c=f"C{v}", marker="o")
    ax[0].scatter([0.,2.,4.,5.], biases.loc[[station_labels.get(s) for s in lighthouses]][vari], s=50, c=f"C{v}", marker="x")

    ax[1].plot(np.arange(len(all_data.keys())), maes[vari], color=f"C{v}")
    ax[1].scatter([1.,3.,6.], maes.loc[[station_labels.get(s) for s in met_stations]][vari], s=20, c=f"C{v}", marker="o")
    ax[1].scatter([0.,2.,4.,5.], maes.loc[[station_labels.get(s) for s in lighthouses]][vari], s=50, c=f"C{v}", marker="x")

ax[0].legend(loc="lower center", ncols=4, bbox_to_anchor=(0.5, 1.05))

ax[0].set_ylabel("bias [K; m/s; \%; hPa]")
ax[1].set_ylabel("MAE [K; m/s; \%; hPa]")

ax[1].set_ylim(bottom=0.)
ax[1].set_xticks(np.arange(len(all_data.keys())))
ax[1].set_xticklabels(list(biases.index))

for a in ax:
    a.grid()
    
plt.savefig(f"{path_out}iwin_paper_fig_eval_lighthouses_errorstats.pdf", dpi=300)


plt.show()



#%% histograms

bins = {"temperature": np.arange(-25., 10.1, 5.), "wind_speed": np.arange(0., 30.1, 5.), "relative_humidity": np.arange(30., 100.1, 10.), "air_pressure": np.arange(960., 1040.1, 10.)}

xlabels = {"temperature": "temperature [°C]", "wind_speed": "wind speed [m/s]", "relative_humidity": "rel. humidity [\%]", "air_pressure": "pressure [hPa]"}

colors = ["b", "r", "g", "orange"]

fig, axes = plt.subplots(2,2, sharey=True, figsize=latex_helpers.set_size(503.6, whr=0.6))
ax = axes.flatten()
for v, vari in enumerate(["temperature", "relative_humidity", "wind_speed", "air_pressure"]):
    
    ax[v].hist([lighthouse_data[s][vari] for s in lighthouse_data.keys()], bins=bins[vari], color=colors, label=list(lighthouse_data.keys()))
    ax[v].set_xlabel(xlabels[vari])
    ax[v].grid()
    
ax[0].legend(ncols=4, bbox_to_anchor=(2.3, 1.3))
fig.subplots_adjust(hspace=0.35)
plt.show()

plt.savefig(f"{path_out}iwin_paper_fig_eval_lighthouses_histograms.pdf", dpi=300)
































