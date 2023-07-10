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

path_boat_data = "https://thredds.met.no/thredds/dodsC/met.no/observations/unis/mobile_AWS"
path_nn_data = "https://thredds.met.no/thredds/dodsC/met.no/observations/unis/lighthouse_AWS_Narveneset/1min/2022/09/lighthouse_AWS_Narveneset_Table_1min_20220920.nc"
path_go_data = "https://thredds.met.no/thredds/dodsC/met.no/observations/unis/lighthouse_AWS_Gasoyane/1min/2022/09/lighthouse_AWS_Gasoyane_Table_1min_20220920.nc"
path_map_data = "/Users/lukasf/OneDrive - Universitetssenteret på Svalbard AS/Svalbard_map_data/"
path_out = "/Users/lukasf/Desktop/Iwin_paper_figures/iwin_paper_fig_15.pdf"
path_out_2 = "/Users/lukasf/Desktop/Iwin_paper_figures/iwin_paper_fig_16.pdf"


day_str = "20220920"
day = pd.to_datetime(day_str)

lat_lims = [78.05, 78.75]
lon_lims = [13.5, 17.5]

lat_lims_FF = [78.15, 78.25]
lon_lims_FF = [15.0, 15.5]

station_heights = {"airport": 28., "nedre_sassendalen": 13., "isfjord_radio": 7., "msbard": 25.,
                   "Narveneset": 15., "Bohemanneset": 12., "Daudmannsodden": 39., "Gasoyane": 30.}



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
        boat_data[s] = ds.sel(time=slice("2022-09-20T00:00:00", "2022-09-21T00:00:00")).load()


for b, b_data in boat_data.items():
    # boat_data[b]["temperature"] -= lighthouse_data["temperature"]
    mask = ~((b_data["latitude"] > 78.22745) & (b_data["latitude"] < 78.22878) & (b_data["longitude"] > 15.60521) & (b_data["longitude"] < 15.61387))
    boat_data[b] = xr.where(mask, b_data, np.nan)
    mask = ~((b_data["latitude"] > 78.06336) & (b_data["latitude"] < 78.06433) & (b_data["longitude"] > 14.1979) & (b_data["longitude"] < 14.20329))
    boat_data[b] = xr.where(mask, b_data, np.nan)
    mask = ~((b_data["latitude"] > 78.65447) & (b_data["latitude"] < 78.65518) & (b_data["longitude"] > 16.37723) & (b_data["longitude"] < 16.38635))
    boat_data[b] = xr.where(mask, b_data, np.nan)

    boat_data[b]["air_pressure"] += station_heights["msbard"]*(1./8.)


boat_data_hr = boat_data["MSPolargirl"].sel(time =slice(pd.Timestamp("2022-09-20 12:30"), pd.Timestamp("2022-09-20 17:00")))


#%%
print("NN")    
with xr.open_dataset(path_nn_data) as ds:
    narveneset_data = ds.sel(time=slice("2022-09-20 08:50:00", "2022-09-20 09:25:00")).load()

print("GO")
with xr.open_dataset(path_go_data) as ds:
    gasoyane_data = ds.sel(time=slice("2022-09-20 11:00:00", "2022-09-20 11:20:00")).load()    

narveneset_data["air_pressure"] += station_heights["Narveneset"]*(1./8.)
gasoyane_data["air_pressure"] += station_heights["Gasoyane"]*(1./8.)


#%%

boat_data_nn = boat_data["MSBard"].sel(time=slice("2022-09-20 08:50:00", "2022-09-20 09:25:00"))
boat_data_go = boat_data["MSBard"].sel(time=slice("2022-09-20 11:00:00", "2022-09-20 11:20:00"))

       
#%% plot

vmin = 0.
vmax = 20.

cmap = cmo.cm.amp

times_arrows = {"MSBard": list(np.array([ pd.Timestamp("2022-09-20 07:05"),
                                    pd.Timestamp("2022-09-20 07:40"),
                                    pd.Timestamp("2022-09-20 09:00"),
                                    pd.Timestamp("2022-09-20 09:58"),
                                    pd.Timestamp("2022-09-20 10:26"),
                                    pd.Timestamp("2022-09-20 10:50"),
                                    pd.Timestamp("2022-09-20 11:30")], dtype="datetime64")),
                "MSPolargirl": list(np.array([pd.Timestamp("2022-09-20 13:48"),
                                               pd.Timestamp("2022-09-20 13:51"),
                                               pd.Timestamp("2022-09-20 14:00"),
                                               pd.Timestamp("2022-09-20 14:05"),
                                               pd.Timestamp("2022-09-20 14:10")], dtype="datetime64")),
                "MSPolargirl_2": list(np.array([pd.Timestamp("2022-09-20 07:58"),
                                                 pd.Timestamp("2022-09-20 08:50"),
                                                 pd.Timestamp("2022-09-20 13:12")], dtype="datetime64"))}

wind_arrow_data = {}
for b, b_data in boat_data.items():
    wind_arrow_data[b] = b_data.sel(time = times_arrows[b])
wind_arrow_data["MSPolargirl_2"] = boat_data["MSPolargirl"].sel(time = times_arrows["MSPolargirl_2"])
    

fig, ax = plt.subplots(1,1, figsize=latex_helpers.set_size(503.6, whr=0.7), subplot_kw={'projection': ccrs.Mercator()})
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

# LYR
ax.scatter([15.63083], [78.22433], color="k", marker='d', s=50, transform=ccrs.PlateCarree(), zorder=200)
# PYR
ax.scatter([16.33432], [78.65312], color="k", marker='^', s=50, transform=ccrs.PlateCarree(), zorder=200)
# BB
ax.scatter([14.21033], [78.06091], color="k", marker='+', lw=3., s=70, transform=ccrs.PlateCarree(), zorder=200)


for b, b_data in boat_data.items():
    df = pd.DataFrame({'latitude': b_data["latitude"], 'longitude': b_data["longitude"], "color": b_data["wind_speed_corrected"]})
    gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.longitude, df.latitude), crs="EPSG:4326")
    gdf = gdf.to_crs(ccrs.Mercator().proj4_init)
    pic = ax.scatter(x=gdf["longitude"], y=gdf["latitude"], c=gdf["color"], s=1, zorder=100, cmap=cmap, transform=ccrs.PlateCarree(), vmin=vmin, vmax=vmax)
    
cbar = plt.colorbar(pic, ax=ax, orientation="vertical")
cbar.ax.set_ylabel("Wind Speed [m s-1]")

df = pd.DataFrame({'latitude': boat_data_nn["latitude"], 'longitude': boat_data_nn["longitude"], "color": boat_data_nn["wind_speed_corrected"]})
gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.longitude, df.latitude), crs="EPSG:4326")
gdf = gdf.to_crs(ccrs.Mercator().proj4_init)
pic = ax.scatter(x=gdf["longitude"], y=gdf["latitude"], c=gdf["color"], s=8, zorder=100, cmap=cmap, transform=ccrs.PlateCarree(), vmin=vmin, vmax=vmax)

df = pd.DataFrame({'latitude': boat_data_go["latitude"], 'longitude': boat_data_go["longitude"], "color": boat_data_go["wind_speed_corrected"]})
gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.longitude, df.latitude), crs="EPSG:4326")
gdf = gdf.to_crs(ccrs.Mercator().proj4_init)
pic = ax.scatter(x=gdf["longitude"], y=gdf["latitude"], c=gdf["color"], s=8, zorder=100, cmap=cmap, transform=ccrs.PlateCarree(), vmin=vmin, vmax=vmax)


b = "MSBard"
u = -np.abs(wind_arrow_data[b]["wind_speed_corrected"]) * np.sin(np.deg2rad(wind_arrow_data[b]["wind_direction_corrected"]))
v = -np.abs(wind_arrow_data[b]["wind_speed_corrected"]) * np.cos(np.deg2rad(wind_arrow_data[b]["wind_direction_corrected"]))
df = pd.DataFrame({'latitude': wind_arrow_data[b]["latitude"], 'longitude': wind_arrow_data[b]["longitude"], "u": 1.94384*u, "v": 1.94384*v})
gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.longitude, df.latitude), crs="EPSG:4326")
gdf = gdf.to_crs(ccrs.Mercator().proj4_init)
ax.barbs(gdf['geometry'].x, gdf['geometry'].y, gdf['u'], gdf['v'], color="k", length=6., linewidth=1., zorder=130)

b = "MSPolargirl_2"
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

ax.plot(lon_lims_FF, [lat_lims_FF[0], lat_lims_FF[0]], "k-", lw=1.5, transform=ccrs.PlateCarree(), zorder=300)
ax.plot(lon_lims_FF, [lat_lims_FF[1], lat_lims_FF[1]], "k-", lw=1.5, transform=ccrs.PlateCarree(), zorder=300)
ax.plot([lon_lims_FF[0], lon_lims_FF[0]], lat_lims_FF, "k-", lw=1.5, transform=ccrs.PlateCarree(), zorder=300)
ax.plot([lon_lims_FF[1], lon_lims_FF[1]], lat_lims_FF, "k-", lw=1.5, transform=ccrs.PlateCarree(), zorder=300)

ax.plot([13.4, lon_lims_FF[0]], [78.41014, lat_lims_FF[1]], "k--", lw=1.5, transform=ccrs.PlateCarree(), zorder=300)
ax.plot([15.33, lon_lims_FF[1]], [78.41014, lat_lims_FF[1]], "k--", lw=1.5, transform=ccrs.PlateCarree(), zorder=300)

ax.text(0.94, 0.95, "(a)", transform=ax.transAxes, fontsize=10)

scale_bar(ax, (0.02, 0.02), 10, text_kwargs={"weight": "bold"})



pos = ax.get_position()
ax_FF = fig.add_axes([pos.x0-0.02, pos.y0+pos.height/2.0, pos.width/2.0, pos.height/2.0+0.03], projection=ccrs.Mercator())
ax_FF.spines['geo'].set_linewidth(1.5)
ax_FF.set_xticks([15.1, 15.25, 15.4], crs=ccrs.PlateCarree())
ax_FF.set_yticks([78.17, 78.20, 78.23], crs=ccrs.PlateCarree())
ax_FF.xaxis.set_major_formatter(lon_formatter)
ax_FF.yaxis.set_major_formatter(lat_formatter)
ax_FF.set_facecolor("lightblue")
df_coastline = gpd.read_file(f"{path_map_data}NP_S250_SHP/S250_Land_l.shp")
df_coastline = df_coastline.to_crs(ccrs.Mercator().proj4_init)
df_coastline.plot(ax=ax, edgecolor="k", facecolor="none", zorder=20, lw=1.)
dem.plot.imshow(ax=ax_FF, cmap=mpl.colors.ListedColormap([cmo.cm.topo(a) for a in np.linspace(0.6,1.,255)]), levels=np.arange(0, 50. * np.ceil(np.nanmax(dem)/50.)+1., 50.),
                  interpolation=None, add_colorbar=False)

df_layer_glaciers.plot(ax=ax_FF, edgecolor=None, facecolor="#FFFFFF")

df = pd.DataFrame({'latitude': boat_data_hr["latitude"], 'longitude': boat_data_hr["longitude"], "color": boat_data_hr["wind_speed_corrected"]})
gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.longitude, df.latitude), crs="EPSG:4326")
gdf = gdf.to_crs(ccrs.Mercator().proj4_init)
ax_FF.scatter(x=gdf["longitude"], y=gdf["latitude"], c=gdf["color"], s=2, zorder=100, cmap=cmap, transform=ccrs.PlateCarree(), vmin=vmin, vmax=vmax)

b = "MSPolargirl"
u = -np.abs(wind_arrow_data[b]["wind_speed_corrected"]) * np.sin(np.deg2rad(wind_arrow_data[b]["wind_direction_corrected"]))
v = -np.abs(wind_arrow_data[b]["wind_speed_corrected"]) * np.cos(np.deg2rad(wind_arrow_data[b]["wind_direction_corrected"]))
df = pd.DataFrame({'latitude': wind_arrow_data[b]["latitude"], 'longitude': wind_arrow_data[b]["longitude"], "u": 1.94384*u, "v": 1.94384*v})
gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.longitude, df.latitude), crs="EPSG:4326")
gdf = gdf.to_crs(ccrs.Mercator().proj4_init)
ax_FF.barbs(gdf['geometry'].x, gdf['geometry'].y, gdf['u'], gdf['v'], color="k", length=6., linewidth=1., zorder=130)


ax_FF.set_extent(lon_lims_FF+lat_lims_FF, crs=ccrs.PlateCarree())
ax_FF.set_title(None)
ax_FF.set_xlabel(None)
ax_FF.set_ylabel(None)
ax_FF.xaxis.tick_top()
ax_FF.text(0.03, 0.9, "(b)", transform=ax_FF.transAxes, fontsize=10)

scale_bar(ax_FF, (0.78, 0.03), 2, text_kwargs={"weight": "bold"})

ax.text(0.67, 0.015, "© Norwegian Polar Institute", fontsize=7, weight='bold', transform=ax.transAxes, zorder=1001)


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
        lns1 = boat_data_nn[vari_names[vari]].plot(ax=ax[v,0], color="r", lw=1., label="MS Bard")
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
    
    ax[v,0].axvspan(pd.Timestamp("2022-09-20 08:52:00"), pd.Timestamp("2022-09-20 09:20:00"), alpha=0.2, color='red')
    ax[v,1].axvspan(pd.Timestamp("2022-09-20 11:07:00"), pd.Timestamp("2022-09-20 11:10:00"), alpha=0.2, color='red')
    
    
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















