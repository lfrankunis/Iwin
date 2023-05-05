#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 10:28:25 2022

@author: lukasf
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition
import pandas as pd
import geopandas as gpd
import xarray as xr
import rioxarray as rxr
from rioxarray.merge import merge_arrays
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import matplotlib.pyplot as plt
from scalebar import scale_bar
import latex_helpers
import rasterio
import glob
from pyproj import CRS, Proj
from rasterio import plot
from affine import Affine
from datetime import date, timedelta as td

lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()

mpl.rcParams.update({
    #"pgf.texsystem": "lualatex",
    'font.family': 'serif',
    #'text.usetex': True,
    'pgf.rcfonts': False,
    "axes.labelsize": 10,
    "font.size": 10,
    "legend.fontsize": 10,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10})



#%% input

path_iwin_data = r'C:\Users\admin\Desktop\iwin\\'
path_sentinel_data = r'C:\Users\admin\Desktop\iwin\sentinel\\'
path_map_data = r'C:\Users\admin\Desktop\Exjobb_Svalbard\Svalbard_map_data\\'
path_out = r'C:\Users\admin\Desktop\iwin\satellite_image_map.pdf'


lat_lims = [78.2, 78.6]
lon_lims = [13.5, 16.]


#%% load dem data

resolution_el = 20
input_path = path_map_data + 'NP_S0_DTM{b}\S0_DTM{b}.tif'.format(b=resolution_el)
crs_proj4 = ccrs.Mercator()

dem = rxr.open_rasterio(input_path, masked=True).squeeze()
dem = dem.rio.reproject("EPSG:4326")
dem = dem.rio.clip_box(minx=lon_lims[0], miny=lat_lims[0], maxx=lon_lims[1], maxy=lat_lims[1])
dem = dem.where(dem > 0.)

dem = dem.rio.reproject(crs_proj4)

#%% create mosaic of satellite images and reproject (only done once)

wl = "mosaic_TCI"
sat_files = sorted(glob.glob(f"{path_sentinel_data}{wl}.tif"))
tiles = []
for f in sat_files:
    with rxr.open_rasterio(f) as tile:
        tiles.append(tile)
merged_tiles = merge_arrays(tiles)
merged_tiles = merged_tiles.rio.reproject("EPSG:4326")
merged_tiles = merged_tiles.rio.clip_box(minx=lon_lims[0], miny=lat_lims[0], maxx=lon_lims[1], maxy=lat_lims[1])
merged_tiles = merged_tiles.rio.reproject(ccrs.Mercator().proj4_init)

merged_tiles.rio.to_raster(f"{path_sentinel_data}mosaics/mosaic_{wl}.tif")

#%%

with xr.open_mfdataset(f"{path_iwin_data}1min\lighthouse_AWS_1885_Table_1min_*.nc") as ds:
    df_full = ds.to_dataframe()
df_full.set_index(ds.time.values, inplace=True)



def calculate_specific_humidity(rel_hum, temp, press):
    """
    Method to calulate the specific humidity (g/kg) from the specific humidity, temperature and pressure measurements.
    """
    e = 0.01*rel_hum*(6.112 * np.exp((17.62*temp)/(243.12+temp)))
    spec_hum = 1000.*(0.622*e)/(press-0.378*e)
    return spec_hum

df_full["q"] = calculate_specific_humidity(df_full["relative_humidity_Avg"],df_full['temperature'],df_full['air_pressure'])

df_cut = df_full.loc[(df_full.index >= "2022-04-04 00:00:00") & (df_full.index <= "2022-04-06 00:00:00")]


#%% load mosaic of satellite images

wl = "TCI"

with rxr.open_rasterio(f"{path_sentinel_data}mosaic_{wl}.tif") as f:
    merged_tiles = f.squeeze().load()

#%%

df = pd.read_csv(r'C:\Users\admin\Desktop\iwin\bhn_shifts_1.csv', engine = 'python',
                      dayfirst=True, sep=';', header=0, 
                      parse_dates=[0], decimal=".")
t = df["time"]
df["time"] = pd.to_datetime(df["time"], format = '%Y-%m-%d %H:%M')
df.set_index("time", inplace=True)

#%%

df_apr = pd.read_csv(r'C:\Users\admin\Desktop\iwin\bhn_shifts_april_1.csv', engine = 'python',
                      dayfirst=True, sep=';', header=0,  
                      parse_dates=[0], decimal=".")
t1 = df_apr["time"]
df_apr["time"] = pd.to_datetime(df_apr["time"], format = '%Y-%m-%d %H:%M')
df_apr.set_index("time", inplace=True)

#%% plot map

df = df.loc[(df["WD_land"] >= 270) | (df["WD_land"] <= 20)]#|(mob_df1.index >= "2022-03-01 00:00:00") & (mob_df1.index <= "2022-03-31 00:00:00")]
df = df.loc[(df["WD_sea"] >= 20) & (df["WD_sea"] <= 270)]

fig, [ax_right, ax_main] = plt.subplots(1,2, figsize=(10.5,5), gridspec_kw={'width_ratios': [0.8, 1]}, subplot_kw={'projection': ccrs.Mercator()}, dpi=300)#
ax_main.set_xticks([13.5, 14., 14.5, 15., 15.5, 16.], crs=ccrs.PlateCarree())
ax_main.set_yticks([78.2, 78.3, 78.4, 78.5, 78.6], crs=ccrs.PlateCarree())
ax_main.tick_params(axis="y", labelsize=4, width=0.2)
ax_main.tick_params(axis="x", labelsize=4, width=0.2)
ax_main.xaxis.set_major_formatter(lon_formatter)
ax_main.yaxis.set_major_formatter(lat_formatter)
ax_main.text(0.625, 0.015, "© European Space Agency - ESA", fontsize=3, weight='bold', transform=ax_main.transAxes)#, zorder=1001)


df_coastline = gpd.read_file(f"{path_map_data}NP_S250_SHP/S250_Land_l.shp")
df_coastline = df_coastline.to_crs(ccrs.Mercator().proj4_init)
df_coastline.plot(ax=ax_main, edgecolor="k", facecolor="none",lw=0.3)#, zorder=10, lw=1.)
merged_tiles.plot.imshow(ax=ax_main, rgb="band")

ax_main.set_extent(lon_lims+lat_lims, crs=ccrs.PlateCarree())
ax_main.set_title(None)
ax_main.set_xlabel(None)
ax_main.set_ylabel(None)


inset_size = .4

lon0 = 14.75300
lat0 = 78.38166

# BHN
ax_main.scatter([lon0], [lat0], color="g", marker='*', lw=0.5, s=15, transform=ccrs.PlateCarree(), zorder=200)

# project these to xy map (data) coords
x, y = ax_main.projection.transform_point(lon0, lat0, ccrs.PlateCarree())

# add a polar projection axes to the figure
polar_inset = fig.add_axes((0, 0, 1, 1), projection='polar')

# setup a transformation function from data to axes coords and get the axes coords of the polar origin location
data2axes = (ax_main.transAxes + ax_main.transData.inverted()).inverted()
xp, yp = data2axes.transform((x, y))

# set the position of the polar projection axes
ip = InsetPosition(ax_main, [xp - inset_size / 2, yp - inset_size / 2, inset_size, inset_size])
polar_inset.set_axes_locator(ip)

# set up the polar axes
polar_inset.set_theta_direction(-1)
polar_inset.set_theta_zero_location('N')
polar_inset.set_facecolor('none')
polar_inset.tick_params(labelleft=False, labelbottom=False)
polar_inset.grid(axis="x")
polar_inset.grid(axis="y", lw=0.3)
polar_inset.set_yticks([0., .9])# fontsize=4)
polar_inset.set_ylim([0., 1.])#, fontsize=4)
polar_inset.tick_params(width=0.2)
polar_inset.spines.polar.set_visible(False)


# plot wind directions
polar_inset.scatter(x=np.deg2rad(df_apr["WD_sea"][df_apr["WD_sea"].notnull()]), y=np.ones(len(df_apr["WD_sea"][df_apr["WD_sea"].notnull()]))-0.1, s=0.05 * abs(df_apr["T_diff"][df_apr["WD_sea"].notnull()]), color="r", zorder=15)
polar_inset.scatter(x=np.deg2rad(df_apr["WD_land"][df_apr["WD_land"].notnull()]), y=np.ones(len(df_apr["WD_land"][df_apr["WD_land"].notnull()]))-0.1, s=0.05 * abs(df_apr["T_diff"][df_apr["WD_land"].notnull()]), color='b', zorder=15)

ax_main.set_title('b)', fontsize=8)

for k, spine in ax_main.spines.items():
    spine.set_zorder(14)
    spine.set_linewidth(0.2)


ax_right.remove()

#times for shifts
z = pd.to_datetime(["2022-04-04 03:38:00", "2022-04-04 12:33:00", "2022-04-04 14:14:00",
"2022-04-04 22:31:00", "2022-04-05 02:27:00", "2022-04-05 05:10:00",
"2022-04-05 06:50:00", "2022-04-05 08:25:00", "2022-04-05 19:57"])

#times "before and after" shifts
a = pd.to_datetime([z[0] + td(minutes = 10), z[1], z[2] - td(minutes = 10), z[3] + td(minutes = 10),
z[4]- td(minutes = 10), z[5], z[6] - td(minutes = 10), z[7], z[8] - td(minutes = 10)])

b = pd.to_datetime([z[0], z[1] + td(minutes = 10), z[2], z[3],
z[4], z[5] - td(minutes = 10), z[6], z[7]- td(minutes = 10), z[8]])

#plot times series and markers
ax2 = fig.add_axes([0.05, 
                  0.725,
                  0.38,
                  0.14])
df_cut["temperature"].plot(ax=ax2,color="black", linewidth=0.5, zorder=1)
ax2.scatter(a, df_apr["T_land"], c="b", s = 1,zorder=4)
ax2.scatter(b, df_apr["T_sea"], c="r", s = 1,zorder=5)
ax2.set(xticklabels=[])
ax2.set_ylabel("T [°C]", fontsize=4, labelpad=0.3)
ax2.yaxis.set_label_position("left")
ax2.yaxis.set_ticks_position("left")
ax2.tick_params(axis="y", labelsize=3, width=0.2)
ax2.set_yticks([-18,-14,-10], ["-18","-14","-10"])
ax2.tick_params(axis="x", width=0.2)
ax2.grid(linewidth=0.2, alpha=0.5)

for axisss in ['top', 'bottom', 'left', 'right']:

    ax2.spines[axisss].set_linewidth(0.2)  # change width

ax2.set_title('a)',fontsize=8)

ax3 = fig.add_axes([0.05, 
                  0.53,
                  0.38,
                  0.14])
ax3.plot(df_cut.index, df_cut["wind_speed"], c="black", linewidth = 0.4)
ax3.scatter(a, df_apr["WS_land"], c="b", s = 0.6,zorder=4)
ax3.scatter(b, df_apr["WS_sea"], c="r", s = 0.6,zorder=5)
ax3.set(xticklabels=[])
ax3.set_ylabel("WS [m/s]", fontsize=4, labelpad=0.3)
ax3.yaxis.set_label_position("left")
ax3.yaxis.set_ticks_position("left")
ax3.tick_params(axis="y", labelsize=3, width=0.2)
ax3.set_yticks([0,10,20])
ax3.tick_params(axis="x", width=0.2)
ax3.grid(linewidth=0.2, alpha=0.5)

for axisss in ['top', 'bottom', 'left', 'right']:

    ax3.spines[axisss].set_linewidth(0.2)

ax4 = fig.add_axes([0.05, 
                  0.34,
                  0.38,
                  0.14])
ax4.plot(df_cut.index, df_cut["q"], c="black", linewidth = 0.5)
ax4.scatter(a, df_apr["q_land"], c="b", s = 1,zorder=4)
ax4.scatter(b, df_apr["q_sea"], c="r", s = 1,zorder=5)
ax4.set(xticklabels=[])
ax4.set_ylabel("q [g/kg]", fontsize=4, labelpad=0.3)
ax4.yaxis.set_label_position("left")
ax4.yaxis.set_ticks_position("left")
ax4.tick_params(axis="y", labelsize=3, width=0.2)
ax4.set_yticks([0.5,1,1.5])
ax4.tick_params(axis="x", width=0.2)
ax4.grid(linewidth=0.2, alpha=0.5)

for axisss in ['top', 'bottom', 'left', 'right']:

    ax4.spines[axisss].set_linewidth(0.2)

ax5 = fig.add_axes([0.05, 
                  0.155,
                  0.38,
                  0.14])
ax5.scatter(df_cut.index, df_cut["wind_direction"], c="black", s=0.2)
ax5.scatter(a, df_apr["WD_land"], c="b", s = 2)
ax5.scatter(b, df_apr["WD_sea"], c="r", s = 2)
ax5.set_ylabel("WD [°]", fontsize=4, labelpad=0.3)
ax5.yaxis.set_label_position("left")
ax5.yaxis.set_ticks_position("left")
ax5.set_xticklabels(["4 April","","","","","5 April","","",""], fontsize=4)
ax5.tick_params(axis="y", labelsize=3, width=0.2)
ax5.set_yticks([0,90,180,270,360], ["N","E","S","W","N"])
ax5.tick_params(axis="x", width=0.2)
ax5.grid(linewidth=0.2, alpha=0.5)
for axisss in ['top', 'bottom', 'left', 'right']:

    ax5.spines[axisss].set_linewidth(0.2)

plt.tight_layout()
box = ax_main.get_position()
box.x0 = box.x0 + 0.025
box.x1 = box.x1 + 0.025
ax_main.set_position(box)

plt.show()

fig.savefig(f"{path_iwin_data}iwin_paper_fig_10.pdf")   # save the figure to file
plt.close(fig)    # close the figure window