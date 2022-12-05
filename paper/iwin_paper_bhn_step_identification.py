#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 26 12:11:23 2022

@author: lukasf
"""

import pandas as pd
import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from mpl_point_clicker import clicker
from mpl_interactions import zoom_factory, panhandler
from scipy.stats import pearsonr
import cmocean as cmo


path_iwin_data = "/Users/lukasf/OneDrive - Universitetssenteret på Svalbard AS/IWIN/Storage/"


#%% load data

with xr.open_mfdataset(f"{path_iwin_data}lighthouse_AWS_1885/1min/lighthouse_AWS_1885_Table_1min_*.nc") as ds:
    df_full = ds.to_dataframe()
df_full.set_index(ds.time.values, inplace=True)


#%%

land_sector_limits = [270., 360.]

plt.ion()

fig, ax = plt.subplots(2,1, sharex=True, constrained_layout=True)
ax[0].plot(df_full.index, df_full["temperature"], "k-")                          # .rolling(7).median(
ax[0].plot(df_full.index, df_full["temperature_Avg"], "k-")                          # .rolling(7).median(
ax[-1].plot(df_full.index, df_full["wind_direction_Avg"], "k.")
ax[-1].set_ylim((0., 360.))
ax[-1].fill_between(df_full.index, land_sector_limits[0], land_sector_limits[1], color='orange', alpha=0.2)
for a in ax:
    a.grid()
zoom_factory(ax[0])
zoom_factory(ax[1])
ph = panhandler(fig, button=2)
klicker_T = clicker(
    ax[0],
    ["land", "sea"],
    markers=["o", "x"],
    colors=["r", "b"]
)

# plt.show(block=True)


#%%

shifts_dict = klicker_T.get_positions()

# shift_land_times = shifts_dict["land"][:,0]
# shift_sea_times = shifts_dict["sea"][:,0]

# land_T = shifts_dict["land"][:,1]
# sea_T = shifts_dict["sea"][:,1]

# shifts = pd.DataFrame(columns=["T_land", "T_sea", "RH_land", "RH_sea", "WS_land", "WS_sea", "WD_land", "WD_sea", "q_land", "q_sea"])
# for li, lt in enumerate(shift_land_times):
#     si = (np.abs(shift_sea_times - lt)).argmin()
#     st = shift_sea_times[si]
#     middle_time = pd.Timestamp(np.mean([lt, st]), unit="d").round(freq="T")
#     land_time = pd.Timestamp(lt, unit="d").round(freq="T")
#     sea_time = pd.Timestamp(st, unit="d").round(freq="T")
#     data_row = pd.DataFrame({"T_land": land_T[li],
#                              "T_sea": sea_T[si],
#                              "RH_land": df_full["RH"][land_time],
#                              "RH_sea": df_full["RH"][sea_time],
#                              "WS_land": df_full["WS"][land_time],
#                              "WS_sea": df_full["WS"][sea_time],
#                              "WD_land": df_full["WD"][land_time],
#                              "WD_sea": df_full["WD"][sea_time],
#                              "q_land": df_full["q"][land_time],
#                              "q_sea": df_full["q"][sea_time]}, index=[middle_time])
#     shifts = pd.concat([shifts, data_row], axis=0)
# for si, st in enumerate(shift_sea_times):
#     li = (np.abs(shift_land_times - st)).argmin()
#     lt = shift_land_times[li]
#     middle_time = pd.Timestamp(np.mean([lt, st]), unit="d").round(freq="T")
#     land_time = pd.Timestamp(lt, unit="d").round(freq="T")
#     sea_time = pd.Timestamp(st, unit="d").round(freq="T")
#     data_row = pd.DataFrame({"T_land": land_T[li],
#                              "T_sea": sea_T[si],
#                              "RH_land": df_full["RH"][land_time],
#                              "RH_sea": df_full["RH"][sea_time],
#                              "WS_land": df_full["WS"][land_time],
#                              "WS_sea": df_full["WS"][sea_time],
#                              "WD_land": df_full["WD"][land_time],
#                              "WD_sea": df_full["WD"][sea_time],
#                              "q_land": df_full["q"][land_time],
#                              "q_sea": df_full["q"][sea_time]}, index=[middle_time])
#     shifts = pd.concat([shifts, data_row], axis=0)
    
# shifts.drop_duplicates(inplace=True)
# shifts.sort_index(axis=0, inplace=True)

# shifts["T_diff"] = shifts["T_land"] - shifts["T_sea"]
# shifts["RH_diff"] = shifts["RH_land"] - shifts["RH_sea"]
# shifts["WS_diff"] = shifts["WS_land"] - shifts["WS_sea"]
# shifts["q_diff"] = shifts["q_land"] - shifts["q_sea"]

# shifts.to_csv(f'/Users/lukasf/Desktop/Bohemanneset/bohemanneset_shifts_{shifts.index[0].strftime("%Y%m%d%H%M")}_{shifts.index[-1].strftime("%Y%m%d%H%M")}.csv')




#%% plots

df = pd.read_csv("/Users/lukasf/Desktop/Bohemanneset/steps_analysis_output/bohemanneset_shifts.csv", index_col=0, parse_dates=[0])


fig, ax = plt.subplots(2,1,sharex=True)
df.plot(y="T_land", ax=ax[0], c="r", ls="", marker="x", legend=False)
df.plot(y="T_sea", ax=ax[0], c="b", ls="", marker="x", legend=False)
ax[1].plot(df.index, df.T_land-df.T_sea, "k-", ls="", marker="x")
for a in ax:
    a.grid()   


    
plt.show()









#%% rotate wind direction

# with xr.open_mfdataset(f"{path_iwin_data}lighthouse_AWS_1885/1min/lighthouse_AWS_1885_Table_1min_*.nc") as ds:
#     df_full = ds.to_dataframe()
# df_full.set_index(ds.time.values, inplace=True)

df_full["u"] = -np.abs(df_full["wind_speed_Avg"]) * np.sin(np.deg2rad(df_full["wind_direction_Avg"]))
df_full["v"] = -np.abs(df_full["wind_speed_Avg"]) * np.cos(np.deg2rad(df_full["wind_direction_Avg"]))


# rotate  wind vector into onshore-offshore-system
land_sector_limits = [280., 350.]
sea_sector_limits = [260., 10.]




land_sector_middle = np.nanmean(land_sector_limits)
# pure onshore: x deg (clockwise --> rotation = x deg - 360.
angle = land_sector_middle - 360.

land_sector_limits_rot = np.array(land_sector_limits) + angle
offshore_wind_limit = np.nanmean(np.sin(np.deg2rad(land_sector_limits_rot)))
sea_sector_limits_rot = np.array(sea_sector_limits) + angle
onshore_wind_limit = np.nanmean(np.sin(np.deg2rad(sea_sector_limits_rot)))


coast_perp = df_full["u"] * np.cos(np.deg2rad(angle)) + df_full["v"] * np.sin(np.deg2rad(angle))
coast_para = -df_full["u"] * np.sin(np.deg2rad(angle)) + df_full["v"] * np.cos(np.deg2rad(angle))

df_full["WD_coast"] = (np.rad2deg(np.arctan2(-coast_perp, -coast_para)) + 360.) % 360.
df_full["onshore_wind"] = np.sin(np.deg2rad(df_full["WD_coast"]))

df_full.interpolate("linear", inplace=True)
df_full_1 = df_full["2021-08-19 09:38:00":"2022-02-20 23:59:00"]
df_full_2= df_full["2022-02-22 00:02:00":"2022-06-25 23:59:00"]
df_full = pd.concat([df_full_1, df_full_2])



# df_full = pd.concat([df_full,
#                      df_full[["T","onshore_wind"]].rolling("6H").corr()["T"][1:][::2].reset_index()\
#                          .drop(["level_1"], axis=1).set_index("time").rename({"T": "corr_T_onshore_winds"}, axis=1)], axis=1)


df_land = df_full[df_full["onshore_wind"] <= offshore_wind_limit]
df_sea = df_full[df_full["onshore_wind"] > onshore_wind_limit]





#%% plot

fig, ax = plt.subplots(3,1,sharex=True)
df_sea.plot(y="temperature", ax=ax[0], c="r", ls="", marker=".", ms=1, legend=False, label="T_sea")
df_sea.plot(y="temperature_Avg", ax=ax[0], c="r", ls="", marker=".", ms=1, legend=False)
df_land.plot(y="temperature", ax=ax[0], c="b", ls="", marker=".", ms=1, legend=False, label="T_land")
df_land.plot(y="temperature_Avg", ax=ax[0], c="b", ls="", marker=".", ms=1, legend=False)

df_sea.plot(y="onshore_wind", ax=ax[1], c="r", ls="", marker=".", ms=1, legend=False)
df_land.plot(y="onshore_wind", ax=ax[1], c="b", ls="", marker=".", ms=1, legend=False)

df_full.plot(y="wind_direction_Avg", ax=ax[-1], c="k", ls="", marker=".", ms=1, legend=False)
# df_full.plot(y="corr_T_onshore_winds", ax=ax[2], c="k", ls="", marker=".", ms=2, legend=False)

ax[-1].set_ylim((0., 360.))
ax[-1].fill_between(df_full.index, land_sector_limits[0], land_sector_limits[1], color='orange', alpha=0.2)
for a in ax:
    a.grid()   
    
ax[0].legend()


#%%

df_land_daily = {"min": df_land.resample("1D").min(), "mean": df_land.resample("1D").mean(),
                  "median": df_land.resample("1D").median(), "max": df_land.resample("1D").max()}


df_sea_daily = {"min": df_sea.resample("1D").min(), "mean": df_sea.resample("1D").mean(),
                  "median": df_sea.resample("1D").median(), "max": df_sea.resample("1D").max()}

vari = "temperature"


fig, ax = plt.subplots(2,1, sharex=True)
for a in ax:
    a.grid() 
    
for c, stat_combi in enumerate([["median", "median"], ["mean", "mean"]]):


    df_daily = pd.concat([pd.DataFrame(df_land_daily[stat_combi[0]][vari]).rename({vari: "land"}, axis=1),
                          pd.DataFrame(df_sea_daily[stat_combi[1]][vari]).rename({vari: "sea"}, axis=1)], axis=1)
    df_daily["diff"] = df_daily["sea"]-df_daily["land"]
    
    monthly_mean_diff = df_daily["diff"].resample("1M", label="left").mean()
    monthly_mean_diff.index += pd.tseries.frequencies.to_offset("15D")
    
    
    ax[0].plot(df_daily["land"], c=f"C{c}", ls="-")
    ax[0].plot(df_daily["sea"], c=f"C{c}", ls=":")

    pic = ax[1].scatter(x=df_daily.index, y=df_daily["diff"],
                  c=f"C{c}", s=20, marker="x")
    ax[1].plot(monthly_mean_diff, c=f"C{c}", label=stat_combi)


plt.legend()
plt.show()







