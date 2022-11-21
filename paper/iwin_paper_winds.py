#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 20 14:45:27 2022

@author: lukasf
"""
import sys
import pandas as pd
from unisacsi import Meteo as met
import unisacsi
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
path_out = "/Users/lukasf/Desktop/Iwin_paper_figures/iwin_paper_sassenjet.pdf"

day_str = "20220920"
day = pd.to_datetime(day_str)

lat_lims = [78.05, 78.75]
lon_lims = [13.5, 17.5]

lat_lims_FF = [78.15, 78.25]
lon_lims_FF = [15.0, 15.5]


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
boat_data = {"MS_Bard": met.read_IWIN(f"{path_iwin_data}mobile_AWS_1883/20sec/mobile_AWS_1883_Table_20sec_{day_str}.nc"),
             "MS_Polargirl": met.read_IWIN(f"{path_iwin_data}mobile_AWS_1872/20sec/mobile_AWS_1872_Table_20sec_{day_str}.nc")}


lighthouse_data = met.read_IWIN(f"{path_iwin_data}lighthouse_AWS_1887/1min/lighthouse_AWS_1887_Table_1min_{day_str}.nc")
# lighthouse_data = lighthouse_data.interp(time=pd.date_range("2022-08-05 00:00:00", "2022-08-05 23:59:50", freq="20S"), method="linear")

for b, b_data in boat_data.items():
    # boat_data[b]["temperature_Avg"] -= lighthouse_data["temperature_Avg"]
    mask = ~((b_data["latitude"] > 78.22745) & (b_data["latitude"] < 78.22878) & (b_data["longitude"] > 15.60521) & (b_data["longitude"] < 15.61387))
    boat_data[b] = xr.where(mask, b_data, np.nan)
    mask = ~((b_data["latitude"] > 78.06336) & (b_data["latitude"] < 78.06433) & (b_data["longitude"] > 14.1979) & (b_data["longitude"] < 14.20329))
    boat_data[b] = xr.where(mask, b_data, np.nan)
    mask = ~((b_data["latitude"] > 78.65447) & (b_data["latitude"] < 78.65518) & (b_data["longitude"] > 16.37723) & (b_data["longitude"] < 16.38635))
    boat_data[b] = xr.where(mask, b_data, np.nan)


boat_data_hr = boat_data["MS_Polargirl"].sel(time =slice(pd.Timestamp("2022-09-20 12:30"), pd.Timestamp("2022-09-20 17:00")))

    
#%% plot

vmin = 0.
vmax = 20.

cmap = cmo.cm.amp

times_arrows = {"MS_Bard": list(np.array([ pd.Timestamp("2022-09-20 07:05"),
                                    pd.Timestamp("2022-09-20 07:40"),
                                    pd.Timestamp("2022-09-20 09:00"),
                                    pd.Timestamp("2022-09-20 09:58"),
                                    pd.Timestamp("2022-09-20 10:26"),
                                    pd.Timestamp("2022-09-20 10:50"),
                                    pd.Timestamp("2022-09-20 11:30")], dtype="datetime64")),
                "MS_Polargirl": list(np.array([pd.Timestamp("2022-09-20 13:48"),
                                               pd.Timestamp("2022-09-20 13:51"),
                                               pd.Timestamp("2022-09-20 14:00"),
                                               pd.Timestamp("2022-09-20 14:05"),
                                               pd.Timestamp("2022-09-20 14:10")], dtype="datetime64"))}

wind_arrow_data = {}
for b, b_data in boat_data.items():
    wind_arrow_data[b] = b_data.sel(time = times_arrows[b])
    
# lighthouse_arrow_data = lighthouse_data.sel(time = times_arrows["MS_Bard"]+times_arrows["MS_Polargirl"])

fig, ax = plt.subplots(1,1, figsize=latex_helpers.set_size(397.4, whr=0.7), subplot_kw={'projection': ccrs.Mercator()})
ax.set_xticks([14., 15., 16., 17.], crs=ccrs.PlateCarree())
ax.set_yticks([78.1, 78.3], crs=ccrs.PlateCarree())
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
# gl = ax.gridlines(draw_labels=False)
ax.set_facecolor("lightgrey")
met.map_add_coastline(fig, ax, option=1, color="k", lat_limits=lat_lims, lon_limits=lon_lims, path_mapdata=path_map_data)
dem.plot.imshow(ax=ax, cmap=mpl.colors.ListedColormap([cmo.cm.topo(a) for a in np.linspace(0.6,1.,255)]), levels=np.arange(0, 50. * np.ceil(np.nanmax(dem)/50.)+1., 50.),
                  interpolation=None, add_colorbar=False)

df_layer_glaciers.plot(ax=ax, edgecolor=None, facecolor="#FFFFFF")

# LYR
ax.scatter([15.63083], [78.22433], color="k", marker='d', s=50, transform=ccrs.PlateCarree(), zorder=200)
# PYR
ax.scatter([16.33432], [78.65312], color="k", marker='^', s=50, transform=ccrs.PlateCarree(), zorder=200)
# BB
ax.scatter([14.21033], [78.06091], color="k", marker='+', lw=3., s=70, transform=ccrs.PlateCarree(), zorder=200)


for b, b_data in boat_data.items():
    df = pd.DataFrame({'latitude': b_data["latitude"], 'longitude': b_data["longitude"], "color": b_data["wind_speed_corrected_Avg"]})
    gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.longitude, df.latitude), crs="EPSG:4326")
    gdf = gdf.to_crs(ccrs.Mercator().proj4_init)
    pic = ax.scatter(x=gdf["longitude"], y=gdf["latitude"], c=gdf["color"], s=1, zorder=100, cmap=cmap, transform=ccrs.PlateCarree(), vmin=vmin, vmax=vmax)
    
cbar = plt.colorbar(pic, ax=ax, orientation="vertical")
cbar.ax.set_ylabel("Wind Speed [m s-1]")


b = "MS_Bard"
u = -np.abs(wind_arrow_data[b]["wind_speed_corrected_Avg"]) * np.sin(np.deg2rad(wind_arrow_data[b]["wind_direction_corrected_Avg"]))
v = -np.abs(wind_arrow_data[b]["wind_speed_corrected_Avg"]) * np.cos(np.deg2rad(wind_arrow_data[b]["wind_direction_corrected_Avg"]))
df = pd.DataFrame({'latitude': wind_arrow_data[b]["latitude"], 'longitude': wind_arrow_data[b]["longitude"], "u": 1.94384*u, "v": 1.94384*v})
gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.longitude, df.latitude), crs="EPSG:4326")
gdf = gdf.to_crs(ccrs.Mercator().proj4_init)
ax.barbs(gdf['geometry'].x, gdf['geometry'].y, gdf['u'], gdf['v'], color="k", length=6., linewidth=1., zorder=130)



ax.set_extent(lon_lims+lat_lims, crs=ccrs.PlateCarree())
ax.set_title(None)
ax.set_xlabel(None)
ax.set_ylabel(None)

ax.plot(lon_lims_FF, [lat_lims_FF[0], lat_lims_FF[0]], "k-", lw=1.5, transform=ccrs.PlateCarree(), zorder=300)
ax.plot(lon_lims_FF, [lat_lims_FF[1], lat_lims_FF[1]], "k-", lw=1.5, transform=ccrs.PlateCarree(), zorder=300)
ax.plot([lon_lims_FF[0], lon_lims_FF[0]], lat_lims_FF, "k-", lw=1.5, transform=ccrs.PlateCarree(), zorder=300)
ax.plot([lon_lims_FF[1], lon_lims_FF[1]], lat_lims_FF, "k-", lw=1.5, transform=ccrs.PlateCarree(), zorder=300)

ax.plot([13.4, lon_lims_FF[0]], [78.41014, lat_lims_FF[1]], "k--", lw=1.5, transform=ccrs.PlateCarree(), zorder=300)
ax.plot([15.33, lon_lims_FF[1]], [78.41014, lat_lims_FF[1]], "k--", lw=1.5, transform=ccrs.PlateCarree(), zorder=300)

ax.text(0.94, 0.95, "(b)", transform=ax.transAxes, fontsize=10)


pos = ax.get_position()
ax_FF = fig.add_axes([pos.x0-0.02, pos.y0+pos.height/2.0, pos.width/2.0, pos.height/2.0+0.03], projection=ccrs.Mercator())
ax_FF.spines['geo'].set_linewidth(1.5)
ax_FF.set_xticks([15.1, 15.25, 15.4], crs=ccrs.PlateCarree())
ax_FF.set_yticks([78.17, 78.20, 78.23], crs=ccrs.PlateCarree())
ax_FF.xaxis.set_major_formatter(lon_formatter)
ax_FF.yaxis.set_major_formatter(lat_formatter)
ax_FF.set_facecolor("lightgrey")
met.map_add_coastline(fig, ax_FF, option=1, color="k", lat_limits=lat_lims, lon_limits=lon_lims, path_mapdata=path_map_data)
dem.plot.imshow(ax=ax_FF, cmap=mpl.colors.ListedColormap([cmo.cm.topo(a) for a in np.linspace(0.6,1.,255)]), levels=np.arange(0, 50. * np.ceil(np.nanmax(dem)/50.)+1., 50.),
                  interpolation=None, add_colorbar=False)

df_layer_glaciers.plot(ax=ax_FF, edgecolor=None, facecolor="#FFFFFF")

df = pd.DataFrame({'latitude': boat_data_hr["latitude"], 'longitude': boat_data_hr["longitude"], "color": boat_data_hr["wind_speed_corrected_Avg"]})
gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.longitude, df.latitude), crs="EPSG:4326")
gdf = gdf.to_crs(ccrs.Mercator().proj4_init)
ax_FF.scatter(x=gdf["longitude"], y=gdf["latitude"], c=gdf["color"], s=2, zorder=100, cmap=cmap, transform=ccrs.PlateCarree(), vmin=vmin, vmax=vmax)

b = "MS_Polargirl"
u = -np.abs(wind_arrow_data[b]["wind_speed_corrected_Avg"]) * np.sin(np.deg2rad(wind_arrow_data[b]["wind_direction_corrected_Avg"]))
v = -np.abs(wind_arrow_data[b]["wind_speed_corrected_Avg"]) * np.cos(np.deg2rad(wind_arrow_data[b]["wind_direction_corrected_Avg"]))
df = pd.DataFrame({'latitude': wind_arrow_data[b]["latitude"], 'longitude': wind_arrow_data[b]["longitude"], "u": 1.94384*u, "v": 1.94384*v})
gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.longitude, df.latitude), crs="EPSG:4326")
gdf = gdf.to_crs(ccrs.Mercator().proj4_init)
ax_FF.barbs(gdf['geometry'].x, gdf['geometry'].y, gdf['u'], gdf['v'], color="k", length=6., linewidth=1., zorder=130)


ax_FF.set_extent(lon_lims_FF+lat_lims_FF, crs=ccrs.PlateCarree())
ax_FF.set_title(None)
ax_FF.set_xlabel(None)
ax_FF.set_ylabel(None)
ax_FF.xaxis.tick_top()
ax_FF.text(0.03, 0.9, "(a)", transform=ax_FF.transAxes, fontsize=10)


plt.savefig(path_out, dpi=300)


plt.show()




















