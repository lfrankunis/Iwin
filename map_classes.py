#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 11:40:35 2021

@author: unismet
"""


import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import math
import sys
import os
import copy
import pandas as pd
import geopandas as gpd
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import rioxarray as rxr
from scalebar import scale_bar
from shapely.geometry import LineString


ocean = copy.copy(mpl.cm.get_cmap('ocean'))
ocean.set_bad('grey')

terrain = copy.copy(mpl.cm.get_cmap('terrain'))
terrain.set_bad('white')

coolwarm = copy.copy(mpl.cm.get_cmap('coolwarm'))
coolwarm.set_bad('grey')

bwr = copy.copy(mpl.cm.get_cmap('bwr'))
bwr.set_bad('grey')

twilight = copy.copy(mpl.cm.get_cmap('twilight'))
twilight.set_bad('grey')

sali = copy.copy(mpl.cm.get_cmap('summer'))
sali.set_bad('grey')

dens = copy.copy(mpl.cm.get_cmap('BuPu'))
dens.set_bad('grey')

plt.rcParams['contour.negative_linestyle'] = 'solid'

labelsize= 16





def get_paths():
    """
    Method to set all paths to the different data files, depending on e.g. the operating system
    """

    if sys.platform == 'linux':
        if os.path.isdir("/media/lukas/ACSI"):
            pre_path = "/media/lukas/ACSI"
        else:
            pre_path = "/media/lukas/ACSI_backup"
    elif sys.platform == "win32":
        pre_path = "D:"
    path_map_data = '{p}/Data/Svalbard_map_data/'.format(p=pre_path)
    path_ocean_data = '{p}/Data/Ocean/'.format(p=pre_path)
    path_carra = '{p}/Data/CARRA/'.format(p=pre_path)

    return [path_map_data, path_ocean_data, path_carra]








# set the colormap and centre the colorbar
class MidpointNormalize(mpl.colors.Normalize):
    """
    Normalise the colorbar so that diverging bars work there way either side from a prescribed midpoint value)
    e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
    """

    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        mpl.colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))



class Shared_Methods:



    def add_crosssection_line(self, lat, lon, c='k', lw=3):
        """
        Method to simply add the model grid point locations as dots.
        lat and lon 1D-arrays
        """

        df = pd.DataFrame({'latitude': lat, 'longitude': lon, 'section': 1})
        gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.longitude, df.latitude))
        gdf = gdf.groupby(['section'])['geometry'].apply(lambda x: LineString(x.tolist()))
        gdf = gpd.GeoDataFrame(gdf, geometry='geometry', crs="EPSG:4326")
        gdf = gdf.to_crs(self.crs_proj4)

        try:
            gdf.plot(ax=self.ax, color=c, linewidth=lw)
        except AttributeError:
            for i in range(len(self.ax)):
                gdf.plot(ax=self.ax[i], color=c, linewidth=lw)

        return



    def add_grid_points(self, lat, lon, c='k', s=50, marker="o"):
        """
        Method to simply add the model grid point locations as dots.
        lat and lon 1D-arrays
        """

        df = pd.DataFrame({'latitude': lat, 'longitude': lon})
        gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.longitude, df.latitude), crs="EPSG:4326")
        gdf = gdf.to_crs(self.crs_proj4)

        try:
            gdf.plot(ax=self.ax, color=c, markersize=s, marker=marker, zorder=10)
        except AttributeError:
            for i in range(len(self.ax)):
                gdf.plot(ax=self.ax[i], color=c, markersize=s, marker=marker, zorder=10)

        return



    def add_grid_points_scalar(self, lat, lon, data, normalize_midpoint=False, min_cbar_range=False, cbar_switch=True, cbar_label="Temperature [°C]", markersize=50, marker="o", cmap=coolwarm):
        """
        Method to add dots at the locations, colored according to a scalar.
        data, lat and lon 1D-arrays
        """

        df = pd.DataFrame({'latitude': lat, 'longitude': lon, 'data': data})
        gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.longitude, df.latitude), crs="EPSG:4326")
        gdf = gdf.to_crs(self.crs_proj4)

        try:
            if normalize_midpoint:
                gdf.plot(column='data', ax=self.ax, norm=MidpointNormalize(midpoint=normalize_midpoint), legend=cbar_switch, legend_kwds={'label': cbar_label}, marker=marker, markersize=markersize, zorder=10, cmap=cmap)
            elif min_cbar_range:
                if len(data) > 0:
                    cmin = np.nanmin(data) - (min_cbar_range - np.abs(np.nanmax(data)-np.nanmin(data)))/2
                    cmax = np.nanmax(data) + (min_cbar_range - np.abs(np.nanmax(data)-np.nanmin(data)))/2
                else:
                    cmin = 5.
                    cmax = cmin + min_cbar_range
                norm = mpl.colors.Normalize(vmin=cmin, vmax=cmax)
                gdf.plot(column='data', ax=self.ax, norm=norm, legend=cbar_switch, legend_kwds={'label': cbar_label}, marker=marker, markersize=markersize, zorder=10, cmap=cmap)
            else:
                gdf.plot(column='data', ax=self.ax, legend=cbar_switch, legend_kwds={'label': cbar_label}, marker=marker, markersize=markersize, zorder=10, cmap=cmap)
        except AttributeError:
            for i in range(len(self.ax)):
                if normalize_midpoint:
                    gdf.plot(column='data', ax=self.ax[i], norm=MidpointNormalize(midpoint=normalize_midpoint), legend=cbar_switch, legend_kwds={'label': cbar_label, 'orientation': 'horizontal'}, marker=marker, markersize=markersize, zorder=10, cmap=cmap)
                elif min_cbar_range:
                    if len(data) > 0:
                        cmin = np.nanmin(data) - (min_cbar_range - np.abs(np.nanmax(data)-np.nanmin(data)))/2
                        cmax = np.nanmax(data) + (min_cbar_range - np.abs(np.nanmax(data)-np.nanmin(data)))/2
                    else:
                        cmin = 5.
                        cmax = cmin + min_cbar_range
                    norm = mpl.colors.Normalize(vmin=cmin, vmax=cmax)
                    gdf.plot(column='data', ax=self.ax[i], norm=norm, legend=cbar_switch, legend_kwds={'label': cbar_label}, marker=marker, markersize=markersize, zorder=10, cmap=cmap)
                else:
                    gdf.plot(column='data', ax=self.ax[i], legend=cbar_switch, legend_kwds={'label': cbar_label, 'orientation': 'horizontal'}, marker=marker, markersize=markersize, zorder=10, cmap=cmap)

        return



    def add_grid_points_arrows(self, lat, lon, u, v):
        """
        Method to add arrows at the locations specified by lat and lon, direction according to u (pos east) and v (pos north).
        u, v, lat and lon 1D-arrays
        """

        df = pd.DataFrame({'latitude': lat, 'longitude': lon, 'u': u, 'v': v})
        gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.longitude, df.latitude), crs="EPSG:4326")
        gdf = gdf.to_crs(self.crs_proj4)

        try:
            q = self.ax.quiver(gdf['geometry'].x, gdf['geometry'].y, gdf['u'], gdf['v'], zorder=10)
            self.ax.quiverkey(q, 0.95, 0.92, 10, r'$10 \frac{m}{s}$', labelpos='N', coordinates='axes', transform=self.crs_proj4)
        except AttributeError:
            for i in range(len(self.ax)):
                q = self.ax[i].quiver(gdf['geometry'].x, gdf['geometry'].y, gdf['u'], gdf['v'], zorder=10)
                self.ax[i].quiverkey(q, 0.95, 0.92, 10, r'$10 \frac{m}{s}$', labelpos='N', coordinates='axes', transform=self.crs_proj4)

        return



    def add_grid_points_meteo_arrows(self, lat, lon, u, v, length=10, lw=1, zorder=10):
        """
        Method to add meteorological wind arrows at the locations specified by lat and lon, direction according to u (pos east) and v (pos north).
        u, v, lat and lon 1D-arrays, u and v in knots!
        """

        df = pd.DataFrame({'latitude': lat, 'longitude': lon, 'u': u, 'v': v})
        gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.longitude, df.latitude), crs="EPSG:4326")
        gdf = gdf.to_crs(self.crs_proj4)

        try:
            self.ax.barbs(gdf['geometry'].x, gdf['geometry'].y, gdf['u'], gdf['v'], zorder=zorder, length=length, linewidth=lw)
        except AttributeError:
            for i in range(len(self.ax)):
                self.ax[i].barbs(gdf['geometry'].x, gdf['geometry'].y, gdf['u'], gdf['v'], zorder=zorder, length=length, lw=lw)

        return










class bathymetry_map(Shared_Methods):
    """
    Contains a figure with a map of the bathymetry.
    """

    def __init__(self, ax, projection, area, mappath, cbar_switch=True, line_contours=False, scale_length_km=10):

        self.area = area
        self.line_contours = line_contours

        self.path_map_data = mappath

        self.ax = ax
        self.projection = projection

        input_path = self.path_map_data + 'IBCAO/IBCAO_v4_1_200m_t4x1y0.tif'

        # limits from file
        limits = pd.read_csv(self.path_map_data + 'area_latlon_limits.txt', delim_whitespace=True, dtype=str)
        i = limits["area"] == self.area
        list_limits = [float(limits["lon_min"][i]), float(limits["lon_max"][i]), float(limits["lat_min"][i]), float(limits["lat_max"][i])]


        # read elevation data
        dbm = rxr.open_rasterio(input_path, masked=True).squeeze()
        dbm.rio.set_crs(3996)
        dbm = dbm.rio.reproject("EPSG:4326")

        dbm = dbm.rio.clip_box(minx=list_limits[0], miny=list_limits[2], maxx=list_limits[1], maxy=list_limits[3])
        dbm = dbm.where(dbm < 0.)

        # plot the data
        self.crs_proj4 = self.projection.proj4_init
        dbm = dbm.rio.reproject(self.crs_proj4)

        if self.line_contours:
            dbm.plot.contour(ax=self.ax, linewidths=0.5, colors='k', levels=20)
        else:
            pic = dbm.plot.imshow(ax=self.ax, cmap=ocean, levels=20, interpolation=None, add_colorbar=False)
            if cbar_switch:
                cbar = plt.colorbar(pic, pad=0.02)
                cbar.ax.tick_params('y', labelsize=16)
                cbar.ax.set_ylabel('Height [m]', fontsize=16)

        self.ax.set_title(None)

        gl = self.ax.gridlines(draw_labels=False)
        gl.left_labels = True
        gl.bottom_labels = True
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER

        self.ax.set_extent(list_limits, crs=ccrs.PlateCarree())

        # scale bar
        scale_bar(self.ax, (0.9, 0.02), scale_length_km)








#############################################################################################################################################
#############################################################################################################################################




class elevation_map(Shared_Methods):
    """
    Contains a figure with a map of the elevation.
    """

    def __init__(self, ax, projection, area, mappath, resolution=50, cbar_switch=True, line_contours=False, scale_length_km=10):

        self.area = area
        self.resolution = resolution
        self.line_contours = line_contours

        self.path_map_data = mappath

        self.ax = ax
        self.projection = projection


        input_path = self.path_map_data + 'NP_S0_DTM{b}/S0_DTM{b}.tif'.format(b=self.resolution)


        # limits from file
        limits = pd.read_csv(self.path_map_data + 'area_latlon_limits.txt', delim_whitespace=True, dtype=str)
        i = limits["area"] == self.area
        list_limits = [float(limits["lon_min"][i]), float(limits["lon_max"][i]), float(limits["lat_min"][i]), float(limits["lat_max"][i])]


        # read elevation data
        dem = rxr.open_rasterio(input_path, masked=True).squeeze()
        dem = dem.rio.reproject("EPSG:4326")
        dem = dem.rio.clip_box(minx=list_limits[0], miny=list_limits[2], maxx=list_limits[1], maxy=list_limits[3])
        dem = dem.where(dem > 0.)


        # plot the data
        self.crs_proj4 = self.projection.proj4_init
        dem = dem.rio.reproject(self.crs_proj4)

        if self.line_contours:
            dem.plot.contour(ax=self.ax, linewidths=0.5, colors='k', levels=np.arange(0., math.ceil(np.nanmax(dem)*1.e-2)/(1.e-2), 2.*float(self.resolution)))
        else:
            pic = dem.plot.imshow(ax=self.ax, cmap=terrain, interpolation=None, add_colorbar=False)
            if cbar_switch:
                cbar = plt.colorbar(pic, pad=0.02)
                cbar.ax.tick_params('y', labelsize=16)
                cbar.ax.set_ylabel('Height [m]', fontsize=16)

        self.ax.set_title(None)

        gl = self.ax.gridlines(draw_labels=False)
        gl.left_labels = True
        gl.bottom_labels = True
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER

        self.ax.set_extent(list_limits, crs=ccrs.PlateCarree())

        # scale bar
        scale_bar(self.ax, (0.9, 0.02), scale_length_km)











class surface_cover_map(Shared_Methods):
    """
    Contains a figure with a map of the surface cover.
    """

    def __init__(self, ax, projection, area, mappath, resolution=250, line_contours=False, scale_length_km=10):

        self.area = area
        self.resolution = resolution
        self.line_contours = line_contours

        self.path_map_data = mappath

        self.ax = ax
        self.projection = projection

        # colors from file
        colors = pd.read_csv(self.path_map_data + 'rgb_colors.txt', dtype=str, comment='%', delim_whitespace=True, index_col=0, header=None).to_dict()[1]

        # limits from file
        limits = pd.read_csv(self.path_map_data + 'area_latlon_limits.txt', delim_whitespace=True, dtype=str)
        i = limits["area"] == self.area
        list_limits = [float(limits["lon_min"][i]), float(limits["lon_max"][i]), float(limits["lat_min"][i]), float(limits["lat_max"][i])]

        if str(self.resolution) == "1000":
            layers = ['Land', 'Vann', 'Isbreer']
        else:
            layers = ['Land', 'Vann', 'Elvesletter', 'Isbreer', 'Morener', 'TekniskSituasjon']

        # set up map
        self.crs_proj4 = self.projection.proj4_init

        #  self.ax.set_facecolor(colors['Hav'])
        for layer in layers:
            if self.line_contours:
                input_file = self.path_map_data + 'NP_S{b}_SHP/S{b}_{c}_l.shp'.format(b=self.resolution, c=layer)
                df_layer = gpd.read_file(input_file)
                df_layer = df_layer.to_crs(self.crs_proj4)
                df_layer.plot(ax=self.ax, edgecolor='k', facecolor=None)
            else:
                input_file = self.path_map_data + 'NP_S{b}_SHP/S{b}_{c}_f.shp'.format(b=self.resolution, c=layer)
                df_layer = gpd.read_file(input_file)
                df_layer = df_layer.to_crs(self.crs_proj4)
                df_layer.plot(ax=self.ax, edgecolor=None, facecolor=colors[layer])

        self.ax.set_title(None)

        gl = self.ax.gridlines(draw_labels=False)
        gl.left_labels = True
        gl.bottom_labels = True
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER

        self.ax.set_extent(list_limits, crs=ccrs.PlateCarree())

        # scale bar
        scale_bar(self.ax, (0.9, 0.02), scale_length_km)





class topo_map(Shared_Methods):
    """
    Contains a figure with a combined map of the surface cover and elevation or alternatively two maps next to each other.
    """

    def __init__(self, ax, projection, area, mappath, resolution_sc=250, resolution_el=50, scale_length_km=10):

        self.area = area
        self.resolution_el = resolution_el
        self.resolution_sc = resolution_sc

        self.path_map_data = mappath

        self.ax = ax
        self.projection = projection


        input_path = self.path_map_data + 'NP_S0_DTM{b}/S0_DTM{b}.tif'.format(b=self.resolution_el)

        # colors from file
        colors = pd.read_csv(self.path_map_data + 'rgb_colors.txt', dtype=str, comment='%', delim_whitespace=True, index_col=0, header=None).to_dict()[1]

        # limits from file
        limits = pd.read_csv(self.path_map_data + 'area_latlon_limits.txt', delim_whitespace=True, dtype=str)
        i = limits["area"] == self.area
        list_limits = [float(limits["lon_min"][i]), float(limits["lon_max"][i]), float(limits["lat_min"][i]), float(limits["lat_max"][i])]


        if str(self.resolution_sc) == "1000":
            layers = ['Land', 'Vann', 'Isbreer']
        else:
            layers = ['Land', 'Vann', 'Elvesletter', 'Isbreer', 'Morener', 'TekniskSituasjon']


        # read elevation data
        dem = rxr.open_rasterio(input_path, masked=True).squeeze()
        dem = dem.rio.reproject("EPSG:4326")
        dem = dem.rio.clip_box(minx=list_limits[0], miny=list_limits[2], maxx=list_limits[1], maxy=list_limits[3])
        dem = dem.where(dem > 0.)

        self.crs_proj4 = self.projection.proj4_init
        dem = dem.rio.reproject(self.crs_proj4)

        # combined map
        #self.ax.set_facecolor(colors['Hav'])

        for layer in layers:
            input_file = self.path_map_data + 'NP_S{b}_SHP/S{b}_{c}_f.shp'.format(b=self.resolution_sc, c=layer)
            df_layer = gpd.read_file(input_file)
            df_layer = df_layer.to_crs(self.crs_proj4)
            df_layer.plot(ax=self.ax, edgecolor=None, facecolor=colors[layer])

        # elevation contour lines
        dem.plot.contour(ax=self.ax, linewidths=0.5, colors='k', levels=np.arange(0., math.ceil(np.nanmax(dem)*1.e-2)/(1.e-2), 2.*float(self.resolution_el)))

        # finish map
        self.ax.set_extent(list_limits, crs=ccrs.PlateCarree())

        self.ax.set_title(None)

        gl = self.ax.gridlines(draw_labels=False)
        gl.left_labels = True
        gl.bottom_labels = True
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER

        # scale bar
        scale_bar(self.ax, (0.9, 0.02), scale_length_km)