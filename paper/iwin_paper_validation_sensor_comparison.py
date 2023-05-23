#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 18:42:58 2023

@author: lukasf
"""

import sys
import pandas as pd
import cmocean as cmo
import matplotlib.pyplot as plt
import xarray as xr
import glob
import numpy as np
import matplotlib as mpl
import latex_helpers

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

path_data = "/Users/lukasf/OneDrive - Universitetssenteret på Svalbard AS/IWIN/sensor_comparison_MSBard/"
path_out = "/Users/lukasf/Desktop/Iwin_paper_figures/iwin_paper_fig_eval_sensors_MSBard.pdf"


sensors = [1883, 1924]

data = {}
for s in sensors:
    print(s)
    with xr.open_mfdataset(sorted(glob.glob(f"{path_data}{s}/mobile_AWS_{s}_Table_1min_2022????.nc"))) as ds:
        mask = ~((ds["latitude"] > 78.22745) & (ds["latitude"] < 78.22878) & (ds["longitude"] > 15.60521) & (ds["longitude"] < 15.61387))
        ds = xr.where(mask, ds, np.nan)
        df = ds.to_dataframe()
    df = df.set_index(ds.time.values)
    df.dropna(axis=0, how="all", inplace=True)
    df["air_pressure"][df["air_pressure"] <= 974.] = np.nan
    data[s] = df

#%%
    
wdir_sectors = pd.concat([data[1883]["wind_direction_raw"], data[1924]["wind_direction_raw"]], axis=1).mean(axis=1)
# wdir_sectors = data[1883]["wind_direction_raw"]
wdir_sectors = (((wdir_sectors/45.)+.5).astype(int) % 8) + 1

occ_count = wdir_sectors.groupby(wdir_sectors).count()
occ_count /= 0.01*occ_count.sum()

#%%
maes = {}
overall_maes = {}
for v in ["wind_speed_corrected", "wind_direction_corrected", "temperature", "relative_humidity", "air_pressure"]:
    maes[v] = []
    data_vari = pd.concat([data[1883][v], data[1924][v]], axis=1, keys=[1883, 1924])
    overall_maes[v] = abs((data_vari[1883]-data_vari[1924])).mean()
    overall_std = (data_vari[1883]-data_vari[1924]).std()
    if v == "wind_direction_corrected":
        data_vari[1883][data_vari[1883] > 180.] -= 360.
        data_vari[1924][data_vari[1924] > 180.] -= 360.
        
    data_vari["wdir_sector"] = wdir_sectors
    grouped = data_vari.groupby('wdir_sector')
    
    for name, group in grouped:
        maes[v].append(abs((group[1883]-group[1924])).mean())
        

#%%

# t = "raw"

# grouped_windspeed = {1883: pd.concat([data[1883][f"wind_speed_{t}"], wdir_sectors], axis=1, keys=["wspeed", "wdir_sector"]).groupby("wdir_sector").mean(),
#                      1924: pd.concat([data[1924][f"wind_speed_{t}"], wdir_sectors], axis=1, keys=["wspeed", "wdir_sector"]).groupby("wdir_sector").mean()}
        
#%% plot

theta = np.arange(0., 2.*np.pi+0.01, np.pi/4.)


vari_labels = {"temperature": "temperature [°C]", "wind_speed_corrected": "wind speed [m/s]", "air_pressure": "pressure [hPa]", "relative_humidity": "rel. humidity [\%]", "wind_direction_corrected": "wind direction [°]"}

vari_rticks = {"temperature": [0.05, 0.15], "wind_speed_corrected": [0.25, 0.75], "air_pressure": [1.0], "relative_humidity": [0.25, 0.75], "wind_direction_corrected": [20., 40.]}

fig, ax = plt.subplots(2,3, subplot_kw={"projection": "polar"}, figsize=latex_helpers.set_size(503.6, whr=0.65))
for a in ax.flatten():
    a.set_theta_direction(-1)
    a.set_theta_zero_location('N')
    
for i, (v, e) in enumerate(maes.items()):
    ax.flatten()[i].plot(np.linspace(0., 2.*np.pi, 10000), np.ones((10000))*overall_maes[v], color="k", ls="--")
    ax.flatten()[i].plot(theta, list(e)+[e[0]], label=v, marker="x", color="b")
    ax.flatten()[i].tick_params(axis='x', pad=.5)
    ax.flatten()[i].set_yticks(vari_rticks[v])
    if i < 3:
        ax.flatten()[i].set_title(vari_labels[v], fontsize=12)
    else:
        ax.flatten()[i].set_xlabel(vari_labels[v], fontsize=12)

    
ax[-1,-1].bar(theta[:-1], occ_count, color="b", width=np.pi/5.)
ax[-1,-1].set_xlabel("occurrence [\%]", fontsize=12)
ax[-1,-1].tick_params(axis='x', pad=.5)
ax[-1,-1].set_yticks([20., 40.])
ax[-1,-1].set_yticklabels(["20", "40"])

for i, a in enumerate(ax[:,0]):
    a.set_xticks(theta[:-1])
    if i == 0:
        a.set_xticklabels(["0°", "45°", "90°\n 270°", "135°", "180°\n0°", "225°", "270°", "315°"])
    else:
        a.set_xticklabels(["", "45°", "90°\n 270°", "135°", "180°", "225°", "270°", "315°"])


for i, a in enumerate(ax[:,1]):
    a.set_xticks(theta[:-1])
    if i == 0:
        a.set_xticklabels(["0°", "45°", "90°\n 270°", "135°", "180°\n0°", "225°", "", "315°"])
    else:
        a.set_xticklabels(["", "45°", "90°\n 270°", "135°", "180°", "225°", "", "315°"])

    
for i, a in enumerate(ax[:,2]):
    a.set_xticks(theta[:-1])
    if i == 0:
        a.set_xticklabels(["0°", "45°", "90°", "135°", "180°\n0°", "225°", "", "315°"])
    else:
        a.set_xticklabels(["", "45°", "90°", "135°", "180°", "225°", "", "315°"])

plt.savefig(path_out, dpi=300)

plt.show()