#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 08:37:14 2022

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

path_iwin_data = "/Users/lukasf/OneDrive - Universitetssenteret på Svalbard AS/IWIN/Storage/"
path_map_data = "/Users/lukasf/OneDrive - Universitetssenteret på Svalbard AS/Svalbard_map_data/"
path_out = "/Users/lukasf/Desktop/Iwin_paper_figures/iwin_paper_overview_map.pdf"


lat_lims = [78.0, 78.8]
lon_lims = [12.5, 17.5]



#%% load map data

input_file = f'{path_map_data}NP_S250_SHP/S250_Isbreer_f.shp'
df_layer_glaciers = gpd.read_file(input_file)
df_layer_glaciers = df_layer_glaciers.to_crs(ccrs.Mercator().proj4_init)

dem = rxr.open_rasterio(f"{path_map_data}NP_S0_DTM20/S0_DTM20.tif", masked=True).squeeze()
dem = dem.rio.reproject("EPSG:4326")
dem = dem.rio.clip_box(minx=lon_lims[0], miny=lat_lims[0], maxx=lon_lims[1], maxy=lat_lims[1])
dem = dem.rio.reproject(ccrs.Mercator().proj4_init)
dem = dem.where(dem > 0.)


#%% read data

mobile_stations = [1883, 1872, 1924]
boat_data = {}
for s in mobile_stations:
    print(s)
    with xr.open_mfdataset(f"{path_iwin_data}mobile_AWS_{s}/1min/mobile_AWS_{s}_Table_1min_*.nc") as ds:
        boat_data[s] = ds.load()


for b, b_data in boat_data.items():
    mask = ~((b_data["latitude"] > 78.22745) & (b_data["latitude"] < 78.22878) & (b_data["longitude"] > 15.60521) & (b_data["longitude"] < 15.61387))
    boat_data[b] = xr.where(mask, b_data, np.nan)
    mask = ~((b_data["latitude"] > 78.06336) & (b_data["latitude"] < 78.06433) & (b_data["longitude"] > 14.1979) & (b_data["longitude"] < 14.20329))
    boat_data[b] = xr.where(mask, b_data, np.nan)
    mask = ~((b_data["latitude"] > 78.65447) & (b_data["latitude"] < 78.65518) & (b_data["longitude"] > 16.37723) & (b_data["longitude"] < 16.38635))
    boat_data[b] = xr.where(mask, b_data, np.nan)
    
    
for b, b_data in boat_data.items():
    df = b_data.to_dataframe()
    boat_data[b] = df.set_index(b_data.time.values)
    
all_boat_data = pd.concat([b for b in boat_data.values()], axis=0, ignore_index=True)
all_boat_data = all_boat_data[all_boat_data['latitude'].notna()]


#%% count how often the boats are in a certain location

lat_bins = np.arange(lat_lims[0], lat_lims[1], 0.01)
lon_bins = np.arange(lon_lims[0], lon_lims[1], 0.05)

hist_counts = stats.binned_statistic_2d(all_boat_data.longitude, all_boat_data.latitude, all_boat_data.temperature, bins=[lon_bins, lat_bins], statistic="count")
hist_counts.statistic[hist_counts.statistic == 0] = np.nan


#%%

lighthouses = {1884: {"name": "Narveneset", 'lat': 78.56343,'lon': 16.29687, "abbrev": "\\bf NN", "x_off": -23., "y_off": 0.},
               1885: {"name": "Bohemanneset", 'lat': 78.38166, 'lon': 14.75300, "abbrev": "\\bf BHN", "x_off": -28., "y_off": 7.},
               1886: {"name": "Daudmannsodden", 'lat': 78.21056,'lon': 12.98685, "abbrev": "\\bf DMO", "x_off": -25., "y_off": -15.},
               1887: {"name": "Gasoyane", 'lat': 78.45792,'lon': 16.20082, "abbrev": "\\bf GO", "x_off": 8., "y_off": -5.}}

MET_stations = {"\\bf PYR": {"lat": 78.6557, "lon": 16.3603, "x_off": -30., "y_off": 0.},
                "\\bf BB":  {"lat": 78.0609, "lon": 14.2103, "x_off": 5., "y_off": 0.},
                "\\bf LYR": {"lat": 78.2453, "lon": 15.5015, "x_off": -10., "y_off": -15.},
                "\\bf NS":  {"lat": 78.3313, "lon": 16.6818, "x_off": 5., "y_off": -10.},
                "\\bf IR":  {"lat": 78.0625, "lon": 13.6192, "x_off": 5., "y_off": -10.}}

#%% plot map

fig, ax_main = plt.subplots(1,1, figsize=latex_helpers.set_size(503.6, whr=0.6), subplot_kw={'projection': ccrs.Mercator()})
ax_main.set_xticks([13., 14., 15., 16., 17.], crs=ccrs.PlateCarree())
ax_main.set_yticks([78.1, 78.3, 78.5, 78.7], crs=ccrs.PlateCarree())
ax_main.xaxis.set_major_formatter(lon_formatter)
ax_main.yaxis.set_major_formatter(lat_formatter)
ax_main.set_facecolor("lightgrey")

pic = ax_main.pcolormesh((lon_bins[1:]+lon_bins[:-1])/2., (lat_bins[1:]+lat_bins[:-1])/2., np.transpose(hist_counts.statistic),
                         vmin=0., vmax=1.e3, cmap=cmo.cm.amp, shading="auto", transform=ccrs.PlateCarree())
cbar = plt.colorbar(pic, ax=ax_main, extend="max")
cbar.ax.set_ylabel("\# of observations")

df_coastline = gpd.read_file(f"{path_map_data}NP_S250_SHP/S250_Land_l.shp")
df_coastline = df_coastline.to_crs(ccrs.Mercator().proj4_init)
df_coastline.plot(ax=ax_main, edgecolor="k", facecolor="none", zorder=1300, lw=1.)
dem.plot.imshow(ax=ax_main, cmap=mpl.colors.ListedColormap([cmo.cm.topo(a) for a in np.linspace(0.6,1.,255)]), levels=np.arange(0, 50. * np.ceil(np.nanmax(dem)/50.)+1., 50.),
                  interpolation=None, add_colorbar=False, zorder=100)

df_layer_glaciers.plot(ax=ax_main, edgecolor=None, facecolor="#FFFFFF", zorder=120)


for m, mdata in MET_stations.items():
    ax_main.scatter([mdata["lon"]], [mdata["lat"]], color="k", marker='o', s=30, transform=ccrs.PlateCarree(), zorder=2000)
    ax_main.annotate(m, xy=(mdata["lon"], mdata["lat"]), xytext=(mdata["x_off"], mdata["y_off"]), xycoords=ccrs.PlateCarree()._as_mpl_transform(ax_main), 
                     textcoords="offset points", fontsize=10, color="k", zorder=2000)


for l in lighthouses.values():
    ax_main.scatter([l["lon"]], [l["lat"]], color="c", marker='x', s=70, lw=3, transform=ccrs.PlateCarree(), zorder=2000)
    ax_main.annotate(l["abbrev"], xy=(l["lon"], l["lat"]), xytext=(l["x_off"], l["y_off"]), xycoords=ccrs.PlateCarree()._as_mpl_transform(ax_main), 
                     textcoords="offset points", fontsize=10, color="c", zorder=2000)




ax_main.set_extent(lon_lims+lat_lims, crs=ccrs.PlateCarree())
ax_main.set_title(None)
ax_main.set_xlabel(None)
ax_main.set_ylabel(None)

scale_bar(ax_main, (0.04, 0.03), 10, text_kwargs={"weight": "bold"}, zorder=400)


plt.savefig(path_out, dpi=300)


plt.show()
















