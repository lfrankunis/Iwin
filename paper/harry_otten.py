#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 20 14:45:27 2022

@author: lukasf
"""
import sys
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

path_iwin_data = "/Users/lukasf/OneDrive - Universitetssenteret på Svalbard AS/IWIN/Storage/sorted_by_sensor/"
path_map_data = "/Users/lukasf/OneDrive - Universitetssenteret på Svalbard AS/Svalbard_map_data/"
path_out = "/Users/lukasf/Desktop/harry_otten.png"

day_str = "20220920"
day = pd.to_datetime(day_str)

lat_lims = [78.0, 78.8]
lon_lims = [12.5, 17.5]

lat_lims_sv = [76., 81.]
lon_lims_sv = [7., 30.]

#%% load map data

input_file = f'{path_map_data}NP_S250_SHP/S250_Isbreer_f.shp'
df_layer_glaciers = gpd.read_file(input_file)
df_layer_glaciers = df_layer_glaciers.to_crs(ccrs.Mercator().proj4_init)

dem = rxr.open_rasterio(f"{path_map_data}NP_S0_DTM20/S0_DTM20.tif", masked=True).squeeze()
dem = dem.rio.reproject("EPSG:4326")
dem = dem.rio.clip_box(minx=lon_lims[0], miny=lat_lims[0], maxx=lon_lims[1], maxy=lat_lims[1])
dem = dem.rio.reproject(ccrs.Mercator().proj4_init)
dem = dem.where(dem > 0.)

bathy = rxr.open_rasterio(f"{path_map_data}IBCAO/IBCAO_v4_1_200m_t4x1y0.tif", masked=True).squeeze()
bathy.rio.set_crs(3996)
bathy = bathy.rio.reproject("EPSG:4326")
bathy = bathy.rio.clip_box(minx=5., miny=76., maxx=30., maxy=81.)
bathy = bathy.rio.reproject(ccrs.Mercator().proj4_init)

#%% stations

lighthouses = {1884: {"name": "Narveneset", 'lat': 78.56343,'lon': 16.29687, "abbrev": "\\bf NN", "x_off": -23., "y_off": 0.},
               1885: {"name": "Bohemanneset", 'lat': 78.38166, 'lon': 14.75300, "abbrev": "\\bf BHN", "x_off": -28., "y_off": 7.},
               1886: {"name": "Daudmannsodden", 'lat': 78.21056,'lon': 12.98685, "abbrev": "\\bf DMO", "x_off": -25., "y_off": -15.},
               1887: {"name": "Gasoyane", 'lat': 78.45792,'lon': 16.20082, "abbrev": "\\bf GO", "x_off": 8., "y_off": -5.}}

MET_stations = {"\\bf PYR": {"lat": 78.6557, "lon": 16.3603, "x_off": -30., "y_off": 0.},
                "\\bf LYR": {"lat": 78.2453, "lon": 15.5015, "x_off": -10., "y_off": -15.},
                "\\bf NS":  {"lat": 78.3313, "lon": 16.6818, "x_off": 5., "y_off": -10.},
                "\\bf IR":  {"lat": 78.0625, "lon": 13.6192, "x_off": 5., "y_off": -10.}}



#%% read data
mobile_stations = {1883: "MS_Bard", 1872: "MS_Polargirl"}
boat_data = {}
for s, boat in mobile_stations.items():
    with xr.open_dataset(f"{path_iwin_data}mobile_AWS_{s}/20sec/{day_str[:4]}/{day_str[4:6]}/mobile_AWS_{s}_Table_20sec_{day_str}.nc") as ds:
        boat_data[boat] = ds.load()


for b, b_data in boat_data.items():
    mask = ~((b_data["latitude"] > 78.22745) & (b_data["latitude"] < 78.22878) & (b_data["longitude"] > 15.60521) & (b_data["longitude"] < 15.61387))
    boat_data[b] = xr.where(mask, b_data, np.nan)
    mask = ~((b_data["latitude"] > 78.06336) & (b_data["latitude"] < 78.06433) & (b_data["longitude"] > 14.1979) & (b_data["longitude"] < 14.20329))
    boat_data[b] = xr.where(mask, b_data, np.nan)
    mask = ~((b_data["latitude"] > 78.65447) & (b_data["latitude"] < 78.65518) & (b_data["longitude"] > 16.37723) & (b_data["longitude"] < 16.38635))
    boat_data[b] = xr.where(mask, b_data, np.nan)


    
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
                                               pd.Timestamp("2022-09-20 14:00"),
                                               pd.Timestamp("2022-09-20 14:10"),
                                               pd.Timestamp("2022-09-20 07:58"),
                                               pd.Timestamp("2022-09-20 08:50"),
                                               pd.Timestamp("2022-09-20 13:12")], dtype="datetime64"))}

wind_arrow_data = {}
for b, b_data in boat_data.items():
    wind_arrow_data[b] = b_data.sel(time = times_arrows[b])
    

fig, ax = plt.subplots(1,1, figsize=latex_helpers.set_size(503.6, whr=0.6), subplot_kw={'projection': ccrs.Mercator()})
ax.set_xticks([14., 15., 16., 17.], crs=ccrs.PlateCarree())
ax.set_yticks([78.1, 78.3], crs=ccrs.PlateCarree())
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

for m, mdata in MET_stations.items():
    ax.scatter([mdata["lon"]], [mdata["lat"]], color="orange", marker='x', s=70, lw=3, transform=ccrs.PlateCarree(), zorder=2000)
    # ax.annotate(m, xy=(mdata["lon"], mdata["lat"]), xytext=(mdata["x_off"], mdata["y_off"]), xycoords=ccrs.PlateCarree()._as_mpl_transform(ax), 
    #                  textcoords="offset points", fontsize=10, color="darkorange", zorder=2000)


for l in lighthouses.values():
    ax.scatter([l["lon"]], [l["lat"]], color="b", marker='x', s=70, lw=3, transform=ccrs.PlateCarree(), zorder=2000)
    # ax.annotate(l["abbrev"], xy=(l["lon"], l["lat"]), xytext=(l["x_off"], l["y_off"]), xycoords=ccrs.PlateCarree()._as_mpl_transform(ax), 
    #                  textcoords="offset points", fontsize=10, color="m", zorder=2000)
    

for b, b_data in boat_data.items():
    df = pd.DataFrame({'latitude': b_data["latitude"], 'longitude': b_data["longitude"], "color": b_data["wind_speed_corrected"]})
    gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.longitude, df.latitude), crs="EPSG:4326")
    gdf = gdf.to_crs(ccrs.Mercator().proj4_init)
    pic = ax.scatter(x=gdf["longitude"], y=gdf["latitude"], c=gdf["color"], s=1, zorder=100, cmap=cmap, transform=ccrs.PlateCarree(), vmin=vmin, vmax=vmax)
    
cbar = plt.colorbar(pic, ax=ax, orientation="vertical")
cbar.ax.set_ylabel("Wind Speed [m s-1]")


for b in boat_data.keys():
    u = -np.abs(wind_arrow_data[b]["wind_speed_corrected"]) * np.sin(np.deg2rad(wind_arrow_data[b]["wind_direction_corrected"]))
    v = -np.abs(wind_arrow_data[b]["wind_speed_corrected"]) * np.cos(np.deg2rad(wind_arrow_data[b]["wind_direction_corrected"]))
    df = pd.DataFrame({'latitude': wind_arrow_data[b]["latitude"], 'longitude': wind_arrow_data[b]["longitude"], "u": 1.94384*u, "v": 1.94384*v})
    gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.longitude, df.latitude), crs="EPSG:4326")
    gdf = gdf.to_crs(ccrs.Mercator().proj4_init)
    ax.barbs(gdf['geometry'].x, gdf['geometry'].y, gdf['u'], gdf['v'], color="k", length=6., linewidth=1., zorder=130)

ax.set_extent(lon_lims+lat_lims, crs=ccrs.PlateCarree())
ax.set_title(None)
ax.set_xlabel(None)
ax.set_ylabel(None)


scale_bar(ax, (0.87, 0.03), 10, text_kwargs={"weight": "bold"})

ax_sv = fig.add_axes([-0.05, .47, .5, .45], projection=ccrs.Mercator())
ax_sv.set_extent(lon_lims_sv+lat_lims_sv, crs=ccrs.PlateCarree())
ax_sv.set_facecolor("lightblue")
bathy.plot.imshow(ax=ax_sv, cmap=cmo.cm.topo, norm=mpl.colors.TwoSlopeNorm(0., 200. * np.floor(np.nanmin(bathy)/200.), 200. * np.ceil(np.nanmax(bathy)/200.)),
                  interpolation="none", add_colorbar=False)

ax_sv.plot([lon_lims[0], lon_lims[1], lon_lims[1], lon_lims[0], lon_lims[0]],
            [lat_lims[0], lat_lims[0], lat_lims[1], lat_lims[1], lat_lims[0]], "r-", transform=ccrs.PlateCarree())

ax_sv.xaxis.tick_top()
ax_sv.set_title(None)
ax_sv.set_xlabel(None)
ax_sv.set_ylabel(None)
ax_sv.set_xticks([10., 20.], crs=ccrs.PlateCarree())
ax_sv.set_yticks([77., 78., 79., 80.], crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax_sv.xaxis.set_major_formatter(lon_formatter)
ax_sv.yaxis.set_major_formatter(lat_formatter)

plt.savefig(path_out, dpi=300)


plt.show()








