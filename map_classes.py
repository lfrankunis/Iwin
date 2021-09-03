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
import datetime
from scipy import interpolate
import pandas as pd
from gsw import freezing, density, conversions
import utm
import sqlite3
import geopandas as gpd
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from geopy.distance import geodesic
import rioxarray as rxr
from scalebar import scale_bar

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



def read_all_stations_and_sections(path):

    all_stations = pd.read_excel(path + 'ctd_stations_sections.ods', sheet_name='Stations', index_col=0)
    all_sections = pd.read_excel(path + 'ctd_stations_sections.ods', sheet_name='Sections')

    return [all_stations, all_sections]



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





def unis_hd_stations_time_query(starttime, endtime, section, path):
    """
    Function to extract those measurements, which were made at the given stations in the specified time period.
    """

    accuracy_m = 100.
    delta_lat = accuracy_m / 111000.
    delta_lon = accuracy_m / (111000.*np.cos(78.*(np.pi/180.)))

    input_path = path + 'unis_hd_Isfjorden.db'

    start = datetime.datetime.strptime(starttime, '%Y%m%d').timestamp()
    end = datetime.datetime.strptime(endtime, '%Y%m%d').timestamp()

    stations = {}
    with sqlite3.connect(input_path) as conn:
        cur = conn.cursor()
        for i in section.index.values:
            exec_statement = "SELECT * FROM ctd_attributes WHERE date BETWEEN ? AND ? \
                                                             AND latitude BETWEEN ? AND ? \
                                                             AND longitude BETWEEN ? AND ?"
            cur.execute(exec_statement, (start, end, \
                                         section.loc[i]['latitude'] - delta_lat, section.loc[i]['latitude'] + delta_lat,
                                         section.loc[i]['longitude'] - delta_lon, section.loc[i]['longitude'] + delta_lon))
            rows = cur.fetchall()

            df = pd.DataFrame(rows, columns=['id','station','latitude','longitude','date','depth'])
            df['station'].str.strip()

            stations[i] = df


    keys_to_remove = []
    for key, item in stations.items():
        if item.empty:
            keys_to_remove.append(key)
    for i in keys_to_remove:
        del stations[i]
        section.drop(i, inplace=True)

    return [stations, section]






def unis_hd_extract_profile(id, path):
    """
    Method to extract the data from one profile from the UNIS HD database
    Takes the ID of the profile within the database as input
    """


    input_path = path + 'unis_hd_Isfjorden.db'

    exec_statement = "SELECT P,T,S FROM measurements_{a}".format(a=id)
    with sqlite3.connect(input_path) as conn:
        df = pd.read_sql_query(exec_statement, conn)

    return df





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
            q = self.ax.barbs(gdf['geometry'].x, gdf['geometry'].y, gdf['u'], gdf['v'], zorder=zorder, length=length, linewidth=lw)
        except AttributeError:
            for i in range(len(self.ax)):
                q = self.ax[i].barbs(gdf['geometry'].x, gdf['geometry'].y, gdf['u'], gdf['v'], zorder=zorder, length=length, lw=lw)

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

        limits_utm_x, limits_utm_y, _, _ = utm.from_latlon(np.array(list_limits[2:]), np.array(list_limits[:2]))

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












class crosssection(Shared_Methods):
    """
    Contains a figure with a crosssection along the given coordinates.
    """

    def __init__(self, area, lat, lon, mappath, oceanpath, plotting=False, ocean_only=True, atm_only=False):

        self.area = area
        self.lat = lat
        self.lon = lon
        self.plotting = plotting

        self.path_map_data = mappath
        self.path_ocean_data = oceanpath

        self.resolution = 200.


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

        LON, LAT = np.meshgrid(dbm.x, dbm.y)

        # determine coordinates of base points of crosssection
        cross_lat = np.array([self.lat[0]])
        cross_lon = np.array([self.lon[0]])
        stations_dist = np.array([0])
        for i in range(len(self.lat)-1):
            delta = geodesic((self.lat[i], self.lon[i]), (self.lat[i+1], self.lon[i+1])).m
            stations_dist = np.append(stations_dist, delta)
            number_of_points = int(np.ceil(delta/self.resolution)+1.)
            new_lat = np.linspace(self.lat[i], self.lat[i+1], number_of_points)[1:]
            new_lon = np.linspace(self.lon[i], self.lon[i+1], number_of_points)[1:]
            cross_lat = np.append(cross_lat, new_lat)
            cross_lon = np.append(cross_lon, new_lon)


        cross_dist = np.array([0.])
        for i in range(len(cross_lat)-1):
            cross_dist = np.append(cross_dist, geodesic((cross_lat[i], cross_lon[i]), (cross_lat[i+1], cross_lon[i+1])).m)

        cross_dist = np.cumsum(cross_dist) / 1.e3
        stations_dist = np.cumsum(stations_dist) / 1.e3

        # interpolate height at points defining the crossection
        z = interpolate.griddata((LON.flatten(), LAT.flatten()), np.array(dbm).flatten(), (cross_lon, cross_lat), method='linear')

        if ocean_only:
            self.top = 0.
            self.bottom = np.nanmin(z)-0.1*abs(np.nanmin(z))
        elif atm_only:
            self.top = np.nanmax(z)+0.1*abs(np.nanmax(z))
            self.bottom = -0.1*abs(np.nanmax(z))
        else:
            self.top = np.nanmax(z)+0.1*abs(np.nanmax(z)-np.nanmin(z))
            self.bottom = np.nanmin(z)-0.1*abs(np.nanmax(z)-np.nanmin(z))

        if self.plotting:
            # plot the crosssection
            fig, ax = plt.subplots(1,1,figsize=(18,12))
            ax.fill_between(cross_dist, z, self.bottom, color='grey', zorder=4)
            if atm_only:
                ax.fill_between(cross_dist, 0., self.bottom, color='grey', zorder=4)
            if ((ocean_only == False) & (atm_only == False)):
                ax.axhline(y=0., color='k', linestyle='--', zorder=3)
            ax.set_xlabel("Distance [km]", fontsize=labelsize)
            ax.set_ylabel("Depth [m]", fontsize=labelsize)
            ax.tick_params('both', labelsize=labelsize)
            ax.set_xlim((cross_dist[0], cross_dist[-1]))
            ax.set_ylim((self.bottom, self.top))
            ax.set_yticklabels(ax.get_yticks())
            labels = [str(int(abs(item.get_position()[1]))) for item in ax.get_yticklabels()]
            ax.set_yticklabels(labels)
            plt.tight_layout()

            self.fig = fig
            self.ax = ax

        self.basepoints_dist = cross_dist
        self.basepoints_z = z
        self.stations_dist = stations_dist




    def grid_unishd(self, stations):
        """
        Grid the measurement data from the UNIS HD according to the bathymetry crosssection.
        Adding it afterwards with add_unishd_in_figure
        stations: List with the station numbers of the section
        time_start, time_end: strings in format YYYYmmdd
        """

        self.stations = stations

        d = int(-self.bottom)

        measurements = {}
        temperature = np.ones((d, len(self.stations))) * np.nan
        salinity = np.ones((d,len(self.stations))) * np.nan
        self.depth = -1.*np.arange(d)
        stat_ind = 0
        for station, attributes in self.stations.items():
            list_of_measurements = [unis_hd_extract_profile(i, self.path_ocean_data) for i in attributes['id'].tolist()]
            measurements[station] = np.ones((d,3,len(list_of_measurements))) * np.nan
            for count, m in enumerate(list_of_measurements):
                ind = [int(i) for i in m['P']]
                values = list_of_measurements[count].to_numpy()
                measurements[station][ind,:,count] = values

            measurements[station] = np.nanmean(measurements[station], axis=2)
            temperature[:,stat_ind] = measurements[station][:,1]
            salinity[:,stat_ind] = measurements[station][:,2]
            stat_ind += 1

        SD, PR = np.meshgrid(self.stations_dist, self.depth)
        BD, PR2 = np.meshgrid(self.basepoints_dist, self.depth)

        temperature = temperature.flatten()
        salinity = salinity.flatten()
        ind = ~np.isnan(temperature)
        temperature = temperature[ind]
        salinity = salinity[ind]
        SD = SD.flatten()[ind]
        PR = PR.flatten()[ind]

#        self.temperature = interpolate.griddata((SD, PR), temperature, (BD, PR2), method='linear')
#        self.salinity = interpolate.griddata((SD, PR), salinity, (BD, PR2), method='linear')


        rbf_mode = 'linear'
        rbf_temperature = interpolate.Rbf(SD, PR, temperature, function=rbf_mode, smooth=0)
        rbf_salinity = interpolate.Rbf(SD, PR, salinity, function=rbf_mode, smooth=0)
        self.temperature = rbf_temperature(BD, PR2)
        self.salinity = rbf_salinity(BD, PR2)

        self.pressure = -1.*PR2

        for i in range(len(self.basepoints_dist)):
            ind = np.where(self.depth < self.basepoints_z[i])[0]
            ind = ind[np.min([len(ind)-1, 40])]
            self.temperature[ind:,i] = np.nan
            self.salinity[ind:,i] = np.nan
            self.pressure[ind:,i] = np.nan

        T_f = freezing.t_freezing(self.salinity, self.pressure, 0.5)
        self.temperature[self.temperature < T_f] = T_f[self.temperature < T_f]

        self.absolute_salinity = conversions.SA_from_SP(self.salinity, self.pressure, 13., 78.)
        self.conservative_temperature = conversions.CT_from_t(self.absolute_salinity, self.temperature, self.pressure)

        self.density = density.rho(self.absolute_salinity, self.conservative_temperature, self.pressure)
        self.sigma0 = density.sigma0(self.absolute_salinity, self.conservative_temperature)

        return


    def water_mass_identification(self):
        """
        Function to identify the water masses along the crosssection
        """

        water_mass_def = pd.read_excel(self.path_ocean_data + 'water_masses.ods', sheet_name='Watermasses')

        self.water_masses = np.ones_like(self.temperature) * np.nan
        for index, row in water_mass_def.iterrows():
            ind = np.all(np.array([self.temperature > row['T_min'],
                                    self.temperature <= row['T_max'],
                                    self.salinity > row['S_psu_min'],
                                    self.salinity <= row['S_psu_max']]), axis=0)
            self.water_masses[ind] = index

        return







    def add_unishd_in_figure(self, vari):
        """
        Function to add the gridded temperature or salinity data into a prepared bathymetry crosssection
        """


        water_mass_def = pd.read_excel(self.path_ocean_data + 'water_masses.ods', sheet_name='Watermasses')


        if vari == 'T':
            pic = self.ax.contourf(self.basepoints_dist, self.depth, self.temperature,
                                    cmap=coolwarm, levels=np.linspace(-2.,7.,19))
            cbar = plt.colorbar(pic, fraction= 0.1)
            cbar.ax.tick_params('y', labelsize=16)
            cbar.ax.set_ylabel('Temperature [°C]', fontsize=16)
            cs = self.ax.contour(self.basepoints_dist, self.depth, self.temperature,
                                 colors='w', linewidths=0.8, levels=np.linspace(-1.,6.,8))
            self.ax.clabel(cs, fontsize=labelsize-6, fmt='%1.1f', zorder=3)
        elif vari == 'S':
            pic = self.ax.contourf(self.basepoints_dist, self.depth, self.salinity,
                                    cmap=sali, levels=np.linspace(28., 36., 17))
            cbar = plt.colorbar(pic, fraction= 0.1)
            cbar.ax.tick_params('y', labelsize=16)
            cbar.ax.set_ylabel('Salinity [PSU]', fontsize=16)
            cs = self.ax.contour(self.basepoints_dist, self.depth, self.salinity,
                                    colors='k', linewidths=0.8, levels=np.linspace(29.,35.,7))
            self.ax.clabel(cs, fontsize=labelsize-6, fmt='%1.1f', zorder=3)
        elif vari == 'sigma0':
            pic = self.ax.contourf(self.basepoints_dist, self.depth, self.sigma0,
                                    cmap=dens, levels=np.linspace(25., 29., 17))
            cbar = plt.colorbar(pic, fraction= 0.1)
            cbar.ax.tick_params('y', labelsize=16)
            cbar.ax.set_ylabel('Sigma0 [kg/m³]', fontsize=16)
            cs = self.ax.contour(self.basepoints_dist, self.depth, self.sigma0,
                                    colors='w', linewidths=0.8, levels=np.linspace(25.5,28.5,7))
            self.ax.clabel(cs, fontsize=labelsize-6, fmt='%1.1f', zorder=3)
        elif vari == 'WM':
            pic = self.ax.contourf(self.basepoints_dist, self.depth, self.water_masses,
                                    colors=water_mass_def['color'].tolist(), levels=np.linspace(-0.5, 6.5, 8))
            legend_elements = [mpl.patches.Patch(facecolor=water_mass_def['color'][i], edgecolor=water_mass_def['color'][i],
                                     label=water_mass_def['Abbr'][i]) for i in water_mass_def.index.values if i != 3]
            self.ax.legend(handles=legend_elements, loc='lower right', fontsize=labelsize)


        ax2 = plt.twiny(self.ax)
        ax2.tick_params('x', labelsize=labelsize, labelrotation=90.)
        ax2.set_xlabel('Station', fontsize=labelsize)
        ax2.set_xticks(self.stations_dist)
        ax2.set_xticklabels(list(self.stations.keys()))

        plt.tight_layout()