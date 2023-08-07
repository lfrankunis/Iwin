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

path_msbard_data = "https://thredds.met.no/thredds/dodsC/met.no/observations/unis/mobile_AWS_MSBard_20sec"
path_nn_data = "https://thredds.met.no/thredds/dodsC/met.no/observations/unis/lighthouse_AWS_Narveneset/1min/2022/10/lighthouse_AWS_Narveneset_Table_1min_20221020.nc"
path_go_data = "https://thredds.met.no/thredds/dodsC/met.no/observations/unis/lighthouse_AWS_Gasoyane/1min/2022/10/lighthouse_AWS_Gasoyane_Table_1min_20221020.nc"
path_map_data = "/Users/lukasf/OneDrive - Universitetssenteret på Svalbard AS/Svalbard_map_data/"
path_out = "/Users/lukasf/Desktop/Iwin_paper_figures/iwin_paper_fig_11.pdf"
path_out_2 = "/Users/lukasf/Desktop/Iwin_paper_figures/iwin_paper_fig_12.pdf"

day_str = "20221020"
day = pd.to_datetime(day_str)

lat_lims = [78.37, 78.73]
lon_lims = [15.25, 17.25]

station_heights = {"airport": 28., "nedre_sassendalen": 13., "isfjord_radio": 7., "msbard": 20.,
                   "Narveneset": 7., "Bohemanneset": 12., "Daudmannsodden": 39., "Gasoyane": 23.}



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

with xr.open_dataset(path_msbard_data) as ds:
    boat_data = ds.sel(time=slice("2022-10-20T00:00:00", "2022-10-21T00:00:00")).load()

mask = ~((boat_data["latitude"] > 78.22745) & (boat_data["latitude"] < 78.22878) & (boat_data["longitude"] > 15.60521) & (boat_data["longitude"] < 15.61387))
boat_data = xr.where(mask, boat_data, np.nan)
mask = ~((boat_data["latitude"] > 78.06336) & (boat_data["latitude"] < 78.06433) & (boat_data["longitude"] > 14.1979) & (boat_data["longitude"] < 14.20329))
boat_data = xr.where(mask, boat_data, np.nan)
mask = ~((boat_data["latitude"] > 78.65447) & (boat_data["latitude"] < 78.65518) & (boat_data["longitude"] > 16.37723) & (boat_data["longitude"] < 16.38635))
boat_data = xr.where(mask, boat_data, np.nan)

boat_data["air_pressure"] += station_heights["msbard"]*(1./8.)
#%%
print("NN")    
with xr.open_dataset(path_nn_data) as ds:
    narveneset_data = ds.sel(time=slice("2022-10-20 08:50:00", "2022-10-20 09:10:00")).load()

print("GO")
with xr.open_dataset(path_go_data) as ds:
    gasoyane_data = ds.sel(time=slice("2022-10-20 11:25:00", "2022-10-20 11:45:00")).load()    




narveneset_data["air_pressure"] += station_heights["Narveneset"]*(1./8.)
gasoyane_data["air_pressure"] += station_heights["Gasoyane"]*(1./8.)
                

#%%

boat_data_nn = boat_data.sel(time=slice("2022-10-20 08:50:00", "2022-10-20 09:10:00"))
boat_data_go = boat_data.sel(time=slice("2022-10-20 11:25:00", "2022-10-20 11:45:00"))

#%% plot

vmin = -3.8
vmax = -1.5

cmap = cmo.cm.thermal
cmap_barbs = mpl.colors.ListedColormap([cmo.cm.ice_r(a) for a in np.linspace(0.2,.6,255)])

times_drainage_signals = np.array([ pd.Timestamp("2022-10-20 07:57"),
                                    pd.Timestamp("2022-10-20 08:06"),
                                    pd.Timestamp("2022-10-20 08:29"),
                                    pd.Timestamp("2022-10-20 08:56"),
                                    pd.Timestamp("2022-10-20 09:20"),
                                    pd.Timestamp("2022-10-20 10:04"),
                                    pd.Timestamp("2022-10-20 11:07")], dtype="datetime64")

wind_arrow_data = boat_data.sel(time = times_drainage_signals)

fig, ax = plt.subplots(1,1, figsize=latex_helpers.set_size(503.6, whr=0.7), subplot_kw={'projection': ccrs.Mercator()})
ax.set_xticks([15.5, 16., 16.5, 17.], crs=ccrs.PlateCarree())
ax.set_yticks([78.4, 78.5, 78.6, 78.7], crs=ccrs.PlateCarree())
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
# PYR
ax.scatter([16.33432], [78.65312], color="k", marker='^', s=50, transform=ccrs.PlateCarree(), zorder=100)

# location names
ax.annotate("\\bf Kapp Ekholm", (16.58, 78.575), xycoords=ccrs.PlateCarree()._as_mpl_transform(ax), c="k", weight="bold", fontsize=8, rotation=0., zorder=1000.)
ax.annotate("\\bf Skansbukta", (15.65, 78.53), xycoords=ccrs.PlateCarree()._as_mpl_transform(ax), c="k", weight="bold", fontsize=8, rotation=0., zorder=1000.)
ax.annotate("\\bf Rundodden", (15.7, 78.51), xycoords=ccrs.PlateCarree()._as_mpl_transform(ax), c="k", weight="bold", fontsize=8, rotation=0., zorder=1000.)
ax.annotate("\\bf Narveneset", (16., 78.58), xycoords=ccrs.PlateCarree()._as_mpl_transform(ax), c="k", weight="bold", fontsize=8, rotation=0., zorder=1000.)
ax.annotate("\\bf Pyramiden", (15.98, 78.65), xycoords=ccrs.PlateCarree()._as_mpl_transform(ax), c="k", weight="bold", fontsize=8, rotation=0., zorder=1000.)
ax.annotate("\\bf Nordenskiöldbreen", (17., 78.63), xycoords=ccrs.PlateCarree()._as_mpl_transform(ax), c="k", weight="bold", fontsize=8, rotation=70., zorder=1000.)
# ax.annotate("\\bf Gåsøyane", (16.135, 78.435), xycoords=ccrs.PlateCarree()._as_mpl_transform(ax), c="k", weight="bold", fontsize=8, rotation=0., zorder=1000.)


df = pd.DataFrame({'latitude': boat_data["latitude"], 'longitude': boat_data["longitude"], "color": boat_data["temperature"]})
gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.longitude, df.latitude), crs="EPSG:4326")
gdf = gdf.to_crs(ccrs.Mercator().proj4_init)
pic = ax.scatter(x=gdf["longitude"], y=gdf["latitude"], c=gdf["color"], s=1, zorder=100, cmap=cmap, transform=ccrs.PlateCarree(), vmin=vmin, vmax=vmax)
cbar = plt.colorbar(pic, ax=ax, orientation="vertical")
cbar.ax.set_ylabel("Temperature [°C]")

df = pd.DataFrame({'latitude': boat_data_nn["latitude"], 'longitude': boat_data_nn["longitude"], "color": boat_data_nn["temperature"]})
gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.longitude, df.latitude), crs="EPSG:4326")
gdf = gdf.to_crs(ccrs.Mercator().proj4_init)
pic = ax.scatter(x=gdf["longitude"], y=gdf["latitude"], c=gdf["color"], s=8, zorder=100, cmap=cmap, transform=ccrs.PlateCarree(), vmin=vmin, vmax=vmax)

df = pd.DataFrame({'latitude': boat_data_go["latitude"], 'longitude': boat_data_go["longitude"], "color": boat_data_go["temperature"]})
gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.longitude, df.latitude), crs="EPSG:4326")
gdf = gdf.to_crs(ccrs.Mercator().proj4_init)
pic = ax.scatter(x=gdf["longitude"], y=gdf["latitude"], c=gdf["color"], s=8, zorder=100, cmap=cmap, transform=ccrs.PlateCarree(), vmin=vmin, vmax=vmax)


u = -np.abs(wind_arrow_data["wind_speed_corrected"]) * np.sin(np.deg2rad(wind_arrow_data["wind_direction_corrected"]))
v = -np.abs(wind_arrow_data["wind_speed_corrected"]) * np.cos(np.deg2rad(wind_arrow_data["wind_direction_corrected"]))
df = pd.DataFrame({'latitude': wind_arrow_data["latitude"], 'longitude': wind_arrow_data["longitude"], "u": 1.94384*u, "v": 1.94384*v})
gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.longitude, df.latitude), crs="EPSG:4326")
gdf = gdf.to_crs(ccrs.Mercator().proj4_init)
ax.barbs(gdf['geometry'].x, gdf['geometry'].y, gdf['u'], gdf['v'], wind_arrow_data.time.values.astype(float), length=6., linewidth=1.5, color="deepskyblue", zorder=130)

# u = -np.abs(lighthouse_arrow_data["wind_speed_Avg"]) * np.sin(np.deg2rad(lighthouse_arrow_data["wind_direction_Avg"]))
# v = -np.abs(lighthouse_arrow_data["wind_speed_Avg"]) * np.cos(np.deg2rad(lighthouse_arrow_data["wind_direction_Avg"]))
# df = pd.DataFrame({'latitude': np.ones_like(lighthouse_arrow_data["wind_speed_Avg"].values)*78.45792, 'longitude': np.ones_like(lighthouse_arrow_data["wind_speed_Avg"].values)*16.20082, "u": 1.94384*u, "v": 1.94384*v})
# gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.longitude, df.latitude), crs="EPSG:4326")
# gdf = gdf.to_crs(ccrs.Mercator().proj4_init)
# ax.barbs(gdf['geometry'].x, gdf['geometry'].y, gdf['u'], gdf['v'], lighthouse_arrow_data.time.values.astype(float), length=6., linewidth=1., cmap=cmap_barbs, zorder=430)


ax.set_extent(lon_lims+lat_lims, crs=ccrs.PlateCarree())
ax.set_title(None)
ax.set_xlabel(None)
ax.set_ylabel(None)

scale_bar(ax, (0.02, 0.02), 5, text_kwargs={"weight": "bold"})

ax.text(0.65, 0.015, "© Norwegian Polar Institute", fontsize=7, weight='bold', transform=ax.transAxes, zorder=1001)


plt.savefig(path_out, dpi=300)


plt.show()




#%% comparison with lighthouses

vari_labels = {"temperature": "temperature [°C]", "wind_speed": "wind speed [m/s]", "air_pressure": "pressure [hPa]", "relative_humidity": "rel. hum. [\%]", "wind_direction": "wind direction [°]"}


vari_names = {"temperature": "temperature", "wind_speed": "wind_speed_corrected", "air_pressure": "air_pressure", "relative_humidity": "relative_humidity", "wind_direction": "wind_direction_corrected"}


plt.close("all")




fig, ax = plt.subplots(5,2, sharex="col", sharey="row", figsize=latex_helpers.set_size(503.6, whr=0.9))
for v, vari in enumerate(vari_labels.keys()):
    if vari == "wind_direction":
        boat_data_nn[vari_names[vari]].plot(ax=ax[v,0], color="r", ls=" ", marker=".", ms=1.5)
        boat_data_go[vari_names[vari]].plot(ax=ax[v,1], color="r", ls=" ", marker=".", ms=1.5)
        narveneset_data[vari].plot(ax=ax[v,0], color="deepskyblue", ls=" ", marker=".", ms=1.5)
        gasoyane_data[vari].plot(ax=ax[v,1], color="b", ls=" ", marker=".", ms=1.5)
    elif vari == "temperature":
        lns1 = boat_data_nn[vari_names[vari]].plot(ax=ax[v,0], color="r", lw=1., label="MS Bard/Berg")
        boat_data_go[vari_names[vari]].plot(ax=ax[v,1], color="r", lw=1.)
        lns2 = narveneset_data[vari].plot(ax=ax[v,0], color="deepskyblue", lw=1., label="Narveneset")
        lns3 = gasoyane_data[vari].plot(ax=ax[v,1], color="b", lw=1., label="Gåsøyane")
    else:
        boat_data_nn[vari_names[vari]].plot(ax=ax[v,0], color="r", lw=1.)
        boat_data_go[vari_names[vari]].plot(ax=ax[v,1], color="r", lw=1.)
        narveneset_data[vari].plot(ax=ax[v,0], color="deepskyblue", lw=1.)
        gasoyane_data[vari].plot(ax=ax[v,1], color="b", lw=1.)
    ax[v,0].set_xlabel(None)
    ax[v,0].set_ylabel(vari_labels[vari])
    ax[v,1].set_xlabel(None)
    ax[v,1].set_ylabel(None)
    
    ax[v,0].axvspan(pd.Timestamp("2022-10-20 08:55:00"), pd.Timestamp("2022-10-20 09:02:00"), alpha=0.2, color='red')
    ax[v,1].axvspan(pd.Timestamp("2022-10-20 11:33:00"), pd.Timestamp("2022-10-20 11:38:00"), alpha=0.2, color='red')
    
ax[-1,0].set_yticks([0., 90., 180., 270., 360.])
ax[-1,0].set_yticklabels(["N", "E", "S", "W", "N"])
    
for a in ax.flatten():
    a.grid()
    a.xaxis.set_major_formatter(mpl.dates.DateFormatter('%H:%M'))
    for label in a.get_xticklabels(which='major'):
        label.set(rotation=30, horizontalalignment='center')

lns = lns2+lns1+lns3
labs = [l.get_label() for l in lns]
ax[0,0].legend(lns, labs, ncols=3, bbox_to_anchor=(1.1, 1.5), loc='upper center')



plt.savefig(path_out_2, dpi=300)







