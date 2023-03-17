#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 20 14:45:27 2022

@author: lukasf
"""
import sys
import pandas as pd
import pandas as pd
import cmocean as cmo
import matplotlib.pyplot as plt
import  numpy as np
import geopandas as gpd
from pyproj import Proj
import matplotlib as mpl
import xarray as xr
import rioxarray as rxr
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import latex_helpers
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

path_boat_data = "https://thredds.met.no/thredds/dodsC/met.no/observations/unis/mobile_AWS"
path_bhn_data = "https://thredds.met.no/thredds/dodsC/met.no/observations/unis/lighthouse_AWS_Bohemanneset_1min"
path_map_data = "/Users/lukasf/OneDrive - Universitetssenteret på Svalbard AS/Svalbard_map_data/"
path_out = "/Users/lukasf/Desktop/Iwin_paper_figures/iwin_paper_glaciermicroclimate.pdf"

day_str = "20220805"
day = pd.to_datetime(day_str)

lat_lims = [78.05, 78.55]
lon_lims = [13.5, 16.]


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
mobile_stations = ["MSBard", "MSPolargirl"]
boat_data = {}
for s in mobile_stations:
    print(s)
    with xr.open_dataset(f"{path_boat_data}_{s}_20sec") as ds:
        boat_data[s] = ds.sel(time=slice("2022-08-05T00:00:00", "2022-08-06T00:00:00")).load()



with xr.open_dataset(path_bhn_data) as ds:
    lighthouse_data = ds.sel(time=slice("2022-08-05T00:00:00", "2022-08-06T00:00:00")).load()
lighthouse_data = lighthouse_data.interp(time=pd.date_range("2022-08-05 00:00:00", "2022-08-06 00:00:00", freq="20S"), method="linear")

for b, b_data in boat_data.items():
    boat_data[b]["temperature"] -= lighthouse_data["temperature"]
    mask = ~((b_data["latitude"] > 78.22745) & (b_data["latitude"] < 78.22878) & (b_data["longitude"] > 15.60521) & (b_data["longitude"] < 15.61387))
    boat_data[b] = xr.where(mask, b_data, np.nan)
    mask = ~((b_data["latitude"] > 78.06336) & (b_data["latitude"] < 78.06433) & (b_data["longitude"] > 14.1979) & (b_data["longitude"] < 14.20329))
    boat_data[b] = xr.where(mask, b_data, np.nan)
    mask = ~((b_data["latitude"] > 78.65447) & (b_data["latitude"] < 78.65518) & (b_data["longitude"] > 16.37723) & (b_data["longitude"] < 16.38635))
    boat_data[b] = xr.where(mask, b_data, np.nan)

boat_data["MSBard"] = boat_data["MSBard"].sel(time=slice(pd.Timestamp("2022-08-05 12:30"), pd.Timestamp("2022-08-05 17:00")))
boat_data["MSPolargirl"] = boat_data["MSPolargirl"].sel(time=slice(pd.Timestamp("2022-08-05 06:00"), pd.Timestamp("2022-08-05 17:00")))


#%% times for arrows

times_arrows = {"MSBard": list(np.array([pd.Timestamp("2022-08-05 14:15")])),
                "MSPolargirl": list(np.array([pd.Timestamp("2022-08-05 10:00")], dtype="datetime64"))}

wind_arrow_data = {}
for b, b_data in boat_data.items():
    wind_arrow_data[b] = b_data.sel(time = times_arrows[b])
    
#%% plot

vmin = -3.
vmax = 0.
cmap = cmo.cm.thermal

# vmin = 6.5
# vmax= 9.5
# cmap = cmo.cm.thermal

fig, ax = plt.subplots(1,1, figsize=latex_helpers.set_size(503.6, whr=0.7), subplot_kw={'projection': ccrs.Mercator()})
ax.set_xticks([14., 14.5, 15., 15.5], crs=ccrs.PlateCarree())
ax.set_yticks([78.1, 78.2, 78.3, 78.4, 78.5], crs=ccrs.PlateCarree())
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
# gl = ax.gridlines(draw_labels=False)
ax.set_facecolor("lightblue")
df_coastline = gpd.read_file(f"{path_map_data}NP_S250_SHP/S250_Land_l.shp")
df_coastline = df_coastline.to_crs(ccrs.Mercator().proj4_init)
df_coastline.plot(ax=ax, edgecolor="k", facecolor="none", zorder=20, lw=1.)
dem.plot.imshow(ax=ax, cmap=mpl.colors.ListedColormap([cmo.cm.topo(a) for a in np.linspace(0.6,1.,255)]), levels=np.arange(0, 50. * np.ceil(np.nanmax(dem)/50.)+1., 50.),
                  interpolation=None, add_colorbar=False)

df_layer_glaciers.plot(ax=ax, edgecolor=None, facecolor="#FFFFFF")

# LYR
ax.scatter([15.63083], [78.22433], color="k", marker='d', s=50, transform=ccrs.PlateCarree(), zorder=200)
# BB
ax.scatter([14.21033], [78.06091], color="k", marker='+', lw=2., s=70, transform=ccrs.PlateCarree(), zorder=200)
# BHN
ax.scatter([14.75300], [78.38166], color="k", marker='*', lw=2., s=70, transform=ccrs.PlateCarree(), zorder=200)

for b, b_data in boat_data.items():
    df = pd.DataFrame({'latitude': b_data["latitude"], 'longitude': b_data["longitude"], "color": b_data["temperature"]})
    gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.longitude, df.latitude), crs="EPSG:4326")
    gdf = gdf.to_crs(ccrs.Mercator().proj4_init)
    pic = ax.scatter(x=gdf["longitude"], y=gdf["latitude"], c=gdf["color"], s=1, zorder=100, cmap=cmap, transform=ccrs.PlateCarree(), vmin=vmin, vmax=vmax)
    
    u = -np.abs(wind_arrow_data[b]["wind_speed_corrected"]) * np.sin(np.deg2rad(wind_arrow_data[b]["wind_direction_corrected"]))
    v = -np.abs(wind_arrow_data[b]["wind_speed_corrected"]) * np.cos(np.deg2rad(wind_arrow_data[b]["wind_direction_corrected"]))
    df = pd.DataFrame({'latitude': wind_arrow_data[b]["latitude"], 'longitude': wind_arrow_data[b]["longitude"], "u": 1.94384*u, "v": 1.94384*v})
    gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.longitude, df.latitude), crs="EPSG:4326")
    gdf = gdf.to_crs(ccrs.Mercator().proj4_init)
    ax.barbs(gdf['geometry'].x, gdf['geometry'].y, gdf['u'], gdf['v'], color="deepskyblue", length=8., linewidth=2., zorder=130)

cbar = plt.colorbar(pic, ax=ax, orientation="vertical")
cbar.ax.set_ylabel("Temperature [°C]")

ax.set_extent(lon_lims+lat_lims, crs=ccrs.PlateCarree())
ax.set_title(None)
ax.set_xlabel(None)
ax.set_ylabel(None)

scale_bar(ax, (0.8, 0.03), 10, text_kwargs={"weight": "bold"})

plt.savefig(path_out, dpi=300)


plt.show()





















