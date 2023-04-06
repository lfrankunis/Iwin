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

path_iwin_data = "https://thredds.met.no/thredds/dodsC/met.no/observations/unis/lighthouse_AWS"
path_map_data = "/Users/lukasf/OneDrive - Universitetssenteret på Svalbard AS/Svalbard_map_data/"
path_out = "/Users/lukasf/Desktop/Iwin_paper_figures/iwin_paper_fig"


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


lighthouses = {"Narveneset": {"name": "Narveneset", 'lat': 78.56343,'lon': 16.29687, "abbrev": "\\bf NN", "x_off": -23., "y_off": 0., "map_lon_lims": [15.5, 17.4], "map_lat_lims": [78.46, 78.7], "lonticks": [15.7, 16.2, 16.7, 17.2], "latticks": [78.46, 78.52, 78.58, 78.64, 78.7], "fig_name": "04"},
                "Bohemanneset": {"name": "Bohemanneset", 'lat': 78.38166, 'lon': 14.75300, "abbrev": "\\bf BHN", "x_off": -28., "y_off": 7., "map_lon_lims": [13.5, 16.2], "map_lat_lims": [78.21, 78.55], "lonticks": [13.6, 14.1, 14.6, 15.1, 15.6, 16.1], "latticks": [78.22, 78.3, 78.38, 78.46, 78.54], "fig_name": "03"},
                "Daudmannsodden": {"name": "Daudmannsodden", 'lat': 78.21056,'lon': 12.98685, "abbrev": "\\bf DMO", "x_off": -25., "y_off": -15., "map_lon_lims": [11.8, 14.5], "map_lat_lims": [78.01, 78.36], "lonticks": [11.9, 12.4, 12.9, 13.4, 13.9, 14.4], "latticks": [78.01, 78.08, 78.15, 78.22, 78.29, 78.36], "fig_name": "05"},
                "Gasoyane": {"name": "Gasoyane", 'lat': 78.45792,'lon': 16.20082, "abbrev": "\\bf GO", "x_off": 8., "y_off": -5., "map_lon_lims": [15.25, 17.45], "map_lat_lims": [78.32, 78.6], "lonticks": [15.35, 15.85, 16.35, 16.85, 17.35], "latticks": [78.34, 78.40, 78.46, 78.52, 78.58], "fig_name": "06"}}



#%%

lighthouse_data = {}
for s in lighthouses.keys():
    print(s)
    with xr.open_dataset(f"{path_iwin_data}_{s}_1min") as ds:
        lighthouse_data[s] = ds.to_dataframe()
        lighthouse_data[s].set_index(ds.time.values, inplace=True)


#%% plot map


for l, ldata in lighthouses.items():
    print(l)

    fig, ax_main = plt.subplots(1,1, figsize=latex_helpers.set_size(503.6, whr=0.6), subplot_kw={'projection': ccrs.Mercator()})
    ax_main.set_xticks(ldata["lonticks"], crs=ccrs.PlateCarree())
    ax_main.set_yticks(ldata["latticks"], crs=ccrs.PlateCarree())
    ax_main.xaxis.set_major_formatter(lon_formatter)
    ax_main.yaxis.set_major_formatter(lat_formatter)
    ax_main.set_facecolor("lightblue")
    
    # df_coastline = gpd.read_file(f"{path_map_data}NP_S100_SHP/S100_Land_l.shp")
    # df_coastline = df_coastline.to_crs(ccrs.Mercator().proj4_init)
    # df_coastline.plot(ax=ax_main, edgecolor="k", facecolor="none", zorder=1300, lw=1.)
    dem.plot.imshow(ax=ax_main, cmap=mpl.colors.ListedColormap([cmo.cm.topo(a) for a in np.linspace(0.6,1.,255)]), levels=np.arange(0, 50. * np.ceil(np.nanmax(dem)/50.)+1., 50.),
                      interpolation=None, add_colorbar=False, zorder=100)
    
    df_layer_glaciers.plot(ax=ax_main, edgecolor=None, facecolor="#FFFFFF", zorder=120)
    
    ax_main.set_extent(ldata["map_lon_lims"]+ldata["map_lat_lims"], crs=ccrs.PlateCarree())
    ax_main.set_title(None)
    ax_main.set_xlabel(None)
    ax_main.set_ylabel(None)
    
    
    
    wspeed_bins = np.arange(0. ,21., 5.)
    handles = []
    colors = cmo.cm.thermal(np.linspace(0,1,len(wspeed_bins)))
    for i, ws in enumerate(wspeed_bins):
        if i < len(wspeed_bins)-1:
            handles.append(mpl.patches.Patch(color=colors[i], label=f'{int(ws)}-{int(wspeed_bins[i+1])} m/s'))
        else:
            handles.append(mpl.patches.Patch(color=colors[i], label=f'$>${int(ws)} m/s'))
    ax_main.legend(handles=handles, ncols=len(wspeed_bins), loc="upper center", bbox_to_anchor=(0.5, 1.12), handlelength=1.)
    
    for k, spine in ax_main.spines.items():  #ax.spines is a dictionary
        spine.set_zorder(500)
    
    
    inset_size = .5
    
    # the geographical coords where the polar origin will be placed
    lon0 = ldata["lon"]
    lat0 = ldata["lat"]

    # project these to xy map (data) coords
    x, y = ax_main.projection.transform_point(lon0, lat0, ccrs.PlateCarree())

    # add a polar projection axes to the figure
    polar_inset = fig.add_axes((0, 0, 1, 1), axes_class=windrose.WindroseAxes)

    # setup a transformation function from data to axes coords and get the axes coords of the polar origin location
    data2axes = (ax_main.transAxes + ax_main.transData.inverted()).inverted()
    xp, yp = data2axes.transform((x, y))

    # set the position of the polar projection axes
    ip = InsetPosition(ax_main, [xp - inset_size / 2, yp - inset_size / 2, inset_size, inset_size])
    polar_inset.set_axes_locator(ip)

    # set up the polar axes
    polar_inset.axis("off")
    polar_inset.set_facecolor('none')
    polar_inset.tick_params(labelleft=False, labelbottom=False)
    polar_inset.grid(False)
    
    # plot data
    polar_inset.bar(lighthouse_data[l]["wind_direction"], lighthouse_data[l]["wind_speed"], bins=wspeed_bins, normed=True, opening=0.8, cmap=cmo.cm.thermal)
    
    scale_bar(ax_main, (0.03, 0.03), 10, text_kwargs={"weight": "bold"}, zorder=400)
    
    ax_main.text(0.72, 0.015, "© Norwegian Polar Institute", fontsize=7, weight='bold', transform=ax_main.transAxes, zorder=1000)

    
    plt.savefig(f"{path_out}_{lighthouses[l]['fig_name']}.pdf", dpi=300)
    
    
plt.show()
















