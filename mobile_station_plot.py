#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 11:35:11 2021

@author: unismet
"""


import numpy as np
import map_classes
import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import mobile_AWS_class
import lighthouse_AWS_class
from matplotlib.dates import HourLocator, DateFormatter
import datetime
import timeit
import pymannkendall as mk



timer_start = timeit.default_timer()

path_map_data, path_ocean_data, path_carra= map_classes.get_paths()

plt.rcParams.update({'font.size': 12, "timezone": "Europe/Paris"})

# load boat data
boat = mobile_AWS_class.mobile_AWS(station=1, resolution="1min",
                                    starttime="202108192200", endtime="202108202200",
                                    variables=['temperature', 'pressure', 'relative_humidity', 'wind_speed', 'wind_direction', 'latitude', 'longitude'],
                                    file_type="raw")

boat.filter_GPScoverage()
boat.masks_for_harbors()
boat.calculate_windvector_components(corrected=True)
boat.calculate_wind_sector(corrected=True)
boat.calculate_wind_in_knots(corrected=True)

status = 'live'
if status == "live":
    boat.only_latest_data(datetime.timedelta(hours=12))

    # load lighthouse data --> here only for live plots
    bohemanneset = lighthouse_AWS_class.lighthouse_AWS(station=1885, resolution="1min",
                                        starttime="202108192200", endtime="202108202200",
                                        variables=['temperature', 'pressure', 'relative_humidity', 'wind_speed', 'wind_direction'],
                                        file_type="raw")

    bohemanneset.calculate_windvector_components()
    bohemanneset.calculate_wind_sector()
    bohemanneset.calculate_wind_in_knots()
    bohemanneset.only_latest_data(datetime.timedelta(hours=12))

    lighthouses = {"Bohemanneset": {'lat': 78.38166, 'lon': 14.75300}}

wdir_classes = {1:['N', 0.], 2:['NE', 45.], 3:['E', 90.], 4:['SE', 135.], 5:['S', 180.], 6:['SW', 225.], 7:['W', 270.], 8:['NW', 315.]}

# create figure and subplot grid
fig = plt.figure(figsize=(19,6), constrained_layout=True)
gs = fig.add_gridspec(4,2, width_ratios=[3,2], wspace=0.1)

# plot map
ax_map = fig.add_subplot(gs[:,1], projection=ccrs.Mercator())
sc_map = map_classes.surface_cover_map(ax_map, ccrs.Mercator(), area='Isfjorden', mappath=path_map_data,
                                        resolution=250, scale_length_km=10)

# marker for BB and PYR
settlements = {"LYR": {'lat': 78.22702, 'lon': 15.61482},
               "PYR": {'lat': 78.65486, 'lon': 16.32826},
               "BB":  {'lat': 78.06284, 'lon': 14.21341}}

ax_map.scatter(settlements['BB']['lon'], settlements['BB']['lat'], c="k", marker="d", s=120, zorder=20, transform=ccrs.PlateCarree())
ax_map.scatter(settlements['PYR']['lon'], settlements['PYR']['lat'], c="k", marker="^", s=120, zorder=20, transform=ccrs.PlateCarree())

# name annotations for settlements
transform = ccrs.PlateCarree()._as_mpl_transform(ax_map)
ax_map.annotate("Longyearbyen", (settlements['LYR']['lon'], settlements['LYR']['lat']), xytext=(settlements['LYR']['lon'], settlements['LYR']['lat']-0.04), ha="left", xycoords=transform, zorder=21)
ax_map.annotate("Pyramiden", (settlements['PYR']['lon'], settlements['PYR']['lat']), xytext=(settlements['PYR']['lon']-0.08, settlements['PYR']['lat']+0.01), ha="right", xycoords=transform, zorder=21)
ax_map.annotate("Barentburg", (settlements['BB']['lon'], settlements['BB']['lat']), xytext=(settlements['BB']['lon']+0.1, settlements['BB']['lat']-0.02),  ha="left", xycoords=transform, zorder=21)



# plot boat data
sc_map.add_grid_points_scalar(boat.data['latitude'], boat.data['longitude'], boat.data['temperature'], min_cbar_range=3.,
                                cbar_switch=True, cbar_label="Temperature [°C]", cmap="RdYlBu_r", markersize=15, marker="o")

# all sailing winds
try:
    sc_map.add_grid_points_meteo_arrows(boat.data['latitude'][boat.data["mask_sailing"]][::int(10./boat.resolution_min)], boat.data['longitude'][boat.data["mask_sailing"]][::int(10./boat.resolution_min)],
                                        boat.data['u_knts'][boat.data["mask_sailing"]][::int(10./boat.resolution_min)], boat.data['v_knts'][boat.data["mask_sailing"]][::int(10./boat.resolution_min)], length=7, lw=.8)
except IndexError:
    pass

# latest harbour winds from each harbor, only last two, to not overload the map
try:
    sc_map.add_grid_points_meteo_arrows(boat.data['latitude'][boat.data["mask_LYR"]][-2:], boat.data['longitude'][boat.data["mask_LYR"]][-2:],
                                        boat.data['u_knts'][boat.data["mask_LYR"]][-2:], boat.data['v_knts'][boat.data["mask_LYR"]][-2:], length=7, lw=.8)
except IndexError:
    pass

try:
    sc_map.add_grid_points_meteo_arrows(boat.data['latitude'][boat.data["mask_PYR"]][-2:], boat.data['longitude'][boat.data["mask_PYR"]][-2:],
                                        boat.data['u_knts'][boat.data["mask_PYR"]][-2:], boat.data['v_knts'][boat.data["mask_PYR"]][-2:], length=7, lw=.8)
except IndexError:
    pass
try:
    sc_map.add_grid_points_meteo_arrows(boat.data['latitude'][boat.data["mask_BB"]][-2:], boat.data['longitude'][boat.data["mask_BB"]][-2:],
                                        boat.data['u_knts'][boat.data["mask_BB"]][-2:], boat.data['v_knts'][boat.data["mask_BB"]][-2:], length=7, lw=.8)
except IndexError:
    pass


if status == "live":
    # marker for ship's position
    ax_map.scatter(boat.data['longitude'][-1], boat.data['latitude'][-1], c="orangered", marker="X", s=150, zorder=25, transform=ccrs.PlateCarree())

    # legend
    lc = boat.data["local_time"][-1].strftime("%H:%M")
    text_legend = mpl.lines.Line2D([],[],marker="None",linestyle="None",color="k",label=f"(upd. {lc} [local time])")

    ship_legend = mpl.lines.Line2D([],[],marker="X",markersize=12,linestyle="None",color="orangered",label="Last received position")

    ind = np.where((boat.data['time'] >= boat.data['time'][-1]-datetime.timedelta(hours=3)))[0]
    p = np.nanmean(boat.data['pressure'][ind])
    _, ptrend_bool, _, _, _, _, _, ptend, _ = mk.original_test(boat.data['pressure'][ind])
    if ptrend_bool:
        if ptend > 0:
            pressure_legend = mpl.lines.Line2D([],[],marker=10,markersize=12,linestyle="None",color="lime",label=f"Pressure: {p:.1f} hPa")
        else:
            pressure_legend = mpl.lines.Line2D([],[],marker=11,markersize=12,linestyle="None",color="red",label=f"Pressure: {p:.1f} hPa")
    else:
        pressure_legend = mpl.lines.Line2D([],[],linewidth=3,linestyle="-",color="k",label=f"Pressure: {p:.1f} hPa")

    handles = [pressure_legend, ship_legend, text_legend]
    ax_map.legend(handles=handles, loc="upper left")

    ax_map.annotate("\u265C", (lighthouses['Bohemanneset']['lon'], lighthouses['Bohemanneset']['lat']), xytext=(lighthouses['Bohemanneset']['lon'], lighthouses['Bohemanneset']['lat']), ha="center", va="center", xycoords=transform, zorder=21)
    props = dict(boxstyle='round', facecolor='white', alpha=0.6)
    data_str = '{t:<5} {a:.1f}°C\n{r:<3} {b:.1f}%\n{w:<3} {c:.1f}m/s\n{g:<3} {d}\n(upd. {l})'.format(t="T", r="RH", w="WS", g="WD", l=bohemanneset.data["local_time"][-1].strftime("%H:%M"),
                                                                                             a=bohemanneset.data['temperature'][-1],
                                                                                             b=bohemanneset.data['relative_humidity'][-1],
                                                                                             c=bohemanneset.data['wind_speed'][-1],
                                                                                             d=wdir_classes[bohemanneset.data['wind_sector'][-1]][0])

    ax_map.annotate(data_str, xy=(lighthouses['Bohemanneset']['lon'], lighthouses['Bohemanneset']['lat']), bbox=props, ha="left", va="bottom",
                    xytext=(lighthouses['Bohemanneset']['lon']-1.2, lighthouses['Bohemanneset']['lat']+0.05), xycoords=transform, zorder=21)

    sc_map.add_grid_points_meteo_arrows(np.array([lighthouses['Bohemanneset']['lat'], lighthouses['Bohemanneset']['lat']]), np.array([lighthouses['Bohemanneset']['lon'], lighthouses['Bohemanneset']['lon']]),
                                        np.array([-8., -15.]), np.array([15., 8.]), length=7, lw=.8, zorder=40)


# plot time series

future_time_plot = {'live': 4, 'overview': 0}

plot_variables = {0: 'temperature', 1: 'relative_humidity', 2: 'wind_speed', 3: 'wind_direction'}

plot_colors = {'temperature': '#377eb8',
          'relative_humidity': '#4daf4a',
          'wind_speed': '#984ea3',
          'wind_direction': '#e41a1c'}

plot_labels = {'temperature': 'Temp. [°C]',
          'relative_humidity': 'Rel. Humidity [%]',
          'wind_speed': 'Wind Speed [m/s]',
          'wind_direction': 'Wind Direction [°]'}

plot_units = {'temperature': '°C',
          'relative_humidity': '%',
          'wind_speed': 'm/s',
          'wind_direction': '°'}

# initialize complete wind direction plot, to get the x axis for the following subplots
ax_4 = fig.add_subplot(gs[3,0])
ax_4.scatter(boat.data['local_time'], boat.data["wind_direction"], s=10, color=plot_colors['wind_direction'])



ax_1 = fig.add_subplot(gs[0,0], sharex=ax_4)
ax_2 = fig.add_subplot(gs[1,0], sharex=ax_4)
ax_3 = fig.add_subplot(gs[2,0], sharex=ax_4)


# settings for all 4 time series plots
axes = [ax_1, ax_2, ax_3, ax_4]
for i, ax in enumerate(axes):

    # everything that goes into all but the last frame (wind direction)
    if i < len(axes)-1:
        ax.plot(boat.data['local_time'], boat.data[plot_variables[i]], color=plot_colors[plot_variables[i]])
        ax.set(ylabel = plot_labels[plot_variables[i]])
        plt.setp(ax.get_xticklabels(), visible=False)
        ymin = np.nanmin(boat.data[plot_variables[i]]) - 0.1*np.abs(np.nanmin(boat.data[plot_variables[i]]))
        if ymin == 0.:
            ymin = -0.1
        ymax = np.nanmax(boat.data[plot_variables[i]]) + 0.1*np.abs(np.nanmax(boat.data[plot_variables[i]]))
        if ymax == 0.:
            ymax = 0.1
        ax.set_ylim((ymin, ymax))
        ax.fill_between(boat.data['local_time'], ymin, ymax, where = boat.data['mask_LYR'] == True, facecolor=plot_colors[plot_variables[i]], alpha=0.2)
        ax.yaxis.set_major_locator(plt.MaxNLocator(5))
    # everything that goes only into the last subplot (wind direction)
    else:
        ax.set(xlabel = 'Local time', ylabel = plot_labels['wind_direction'])
        ax.set_xlim((boat.data['local_time'][0], boat.data['local_time'][-1]+datetime.timedelta(hours=future_time_plot[status])))
        ax.set_ylim((0., 360.))
        ax.fill_between(boat.data['local_time'], 0., 360., where = boat.data['mask_LYR'] == True, facecolor=plot_colors[plot_variables[i]], alpha=0.2)
        ax.set_yticks([0., 90., 180., 270., 360.])
        ax.set_yticklabels(["N", "E", "S", "W", "N"])
        ax.xaxis.set_major_locator(HourLocator(interval=2))
        ax.xaxis.set_major_formatter(DateFormatter('%H'))



    # everything that goes into the live plots
    if status == 'live':
        if i < len(axes)-1:
            latest_value_label = "{a:.1f} {p}".format(a=boat.data[plot_variables[i]][-1], p=plot_units[plot_variables[i]])
        else:
            latest_value_label = wdir_classes[boat.data['wind_sector'][-1]][0]
        ax.scatter(boat.data['local_time'][-1], boat.data[plot_variables[i]][-1], s=50, color=plot_colors[plot_variables[i]])
        ax.annotate(latest_value_label, (boat.data['local_time'][-1], boat.data[plot_variables[i]][-1]),
                    xytext=(10,0), textcoords="offset points", va="center", bbox=dict(boxstyle="round", alpha=0.1))

    # plot markers for BB and PYR (if possible)
    try:
        ax.scatter([boat.data['local_time'][boat.data['mask_BB']][0], boat.data['local_time'][boat.data['mask_BB']][-1]],
                    [boat.data[plot_variables[i]][boat.data['mask_BB']][0], boat.data[plot_variables[i]][boat.data['mask_BB']][-1]],
                    color=plot_colors[plot_variables[i]], marker='d', s=250)
    except IndexError:
        pass
    try:
        ax.scatter([boat.data['local_time'][boat.data['mask_PYR']][0], boat.data['local_time'][boat.data['mask_PYR']][-1]],
                    [boat.data[plot_variables[i]][boat.data['mask_PYR']][0], boat.data[plot_variables[i]][boat.data['mask_PYR']][-1]],
                    color=plot_colors[plot_variables[i]], marker='^', s=250)
    except IndexError:
        pass

    ax.grid()

timer_stop = timeit.default_timer()
print('Time: {a}s'.format(a=timer_stop - timer_start))

# fig.savefig("C:/Users/lukasf/Desktop/iwin_mobile_3.png" )
plt.show()