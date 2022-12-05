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
path_sentinel_data = "/Users/lukasf/Desktop/Sentinel/"
path_map_data = "/Users/lukasf/OneDrive - Universitetssenteret på Svalbard AS/Svalbard_map_data/"
path_out = "/Users/lukasf/Desktop/Iwin_paper_figures/satellite_image_map.pdf"


lat_lims = [78.2, 78.6]
lon_lims = [13.5, 16.]


#%% load dem data

# with rxr.open_rasterio(f"{path_map_data}NP_S0_DTM20/S0_DTM20.tif", masked=True) as f:
#     dem = f.squeeze().load()
# dem = dem.rio.reproject("EPSG:4326")
# dem = dem.rio.clip_box(minx=lon_lims[0], miny=lat_lims[0], maxx=lon_lims[1], maxy=lat_lims[1])
# dem = dem.rio.reproject(ccrs.Mercator().proj4_init)
# dem = dem.where(dem > 0.)

#%% create mosaic of satellite images and reproject (only done once)

# wl = "TCI"
# sat_files = sorted(glob.glob(f"{path_sentinel_data}band_data/{wl}/*.tif"))
# tiles = []
# for f in sat_files:
#     with rxr.open_rasterio(f) as tile:
#         tiles.append(tile)
# merged_tiles = merge_arrays(tiles)
# merged_tiles = merged_tiles.rio.reproject("EPSG:4326")
# merged_tiles = merged_tiles.rio.clip_box(minx=lon_lims[0], miny=lat_lims[0], maxx=lon_lims[1], maxy=lat_lims[1])
# merged_tiles = merged_tiles.rio.reproject(ccrs.Mercator().proj4_init)

# merged_tiles.rio.to_raster(f"{path_sentinel_data}mosaics/mosaic_{wl}.tif")

#%% load mosaic of satellite images

wl = "TCI"

with rxr.open_rasterio(f"{path_sentinel_data}mosaics/mosaic_{wl}.tif") as f:
    merged_tiles = f.squeeze().load()

#%% plot map

fig, ax_main = plt.subplots(1,1, figsize=latex_helpers.set_size(503.6, whr=0.6), subplot_kw={'projection': ccrs.Mercator()})
ax_main.set_xticks([13.5, 14., 14.5, 15., 15.5, 16.], crs=ccrs.PlateCarree())
ax_main.set_yticks([78.2, 78.3, 78.4, 78.5, 78.6], crs=ccrs.PlateCarree())
ax_main.xaxis.set_major_formatter(lon_formatter)
ax_main.yaxis.set_major_formatter(lat_formatter)


df_coastline = gpd.read_file(f"{path_map_data}NP_S250_SHP/S250_Land_l.shp")
df_coastline = df_coastline.to_crs(ccrs.Mercator().proj4_init)
df_coastline.plot(ax=ax_main, edgecolor="k", facecolor="none", zorder=10, lw=1.)
# dem.plot.contour(ax=ax_main, linestyles="-", linewidths=0.5, colors="k", levels=np.arange(0, 100. * np.ceil(np.nanmax_main(dem)/100.)+1., 100.), zorder=100)
merged_tiles.plot.imshow(ax=ax_main, rgb="band")

ax_main.set_extent(lon_lims+lat_lims, crs=ccrs.PlateCarree())
ax_main.set_title(None)
ax_main.set_xlabel(None)
ax_main.set_ylabel(None)



inset_size = .4

lon0 = 14.75300
lat0 = 78.38166

# BHN
ax_main.scatter([lon0], [lat0], color="g", marker='*', lw=2., s=70, transform=ccrs.PlateCarree(), zorder=200)


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
polar_inset.set_yticks([0., .9])
polar_inset.set_ylim([0., 1.])
polar_inset.spines.polar.set_visible(False)


# plot some dummy data data
polar_inset.scatter(x=np.deg2rad(np.arange(20., 200., 20.)), y=np.ones(9)-0.1, s=np.arange(11., 20., 20.), color="r", zorder=15)
polar_inset.scatter(x=np.deg2rad(np.arange(270., 360., 10.)), y=np.ones(9)-0.1, s=np.arange(11., 20., 1.), color='b', zorder=15)


for k, spine in ax_main.spines.items():
    spine.set_zorder(14)

plt.savefig(path_out, dpi=300)


plt.show()