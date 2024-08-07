# -*- coding: utf-8 -*-
"""
Created on Wed Sep 22 07:05:14 2021

@author: unismet
"""

import numpy as np
import map_classes
import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from matplotlib.dates import HourLocator, DateFormatter
import datetime
import pymannkendall as mk



wdir_classes = {1:['N', 0.], 2:['NE', 45.], 3:['E', 90.], 4:['SE', 135.], 5:['S', 180.], 6:['SW', 225.], 7:['W', 270.], 8:['NW', 315.]}



def initialize_halfpage_map():
    
    plt.close("all")
    
    plt.rcParams.update({"font.size": 12, "timezone": "Europe/Paris"})
    
    # create figure and subplot grid
    fig = plt.figure(figsize=(19,6), constrained_layout=True)
    gs = fig.add_gridspec(4,2, width_ratios=[3,2], wspace=0.1)

    # plot map
    ax_map = fig.add_subplot(gs[:,1], projection=ccrs.Mercator())
    sc_map = map_classes.surface_cover_map(ax_map, ccrs.Mercator(), area='Isfjorden', mappath="/mnt/c/Svalbard_map_data/",
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
    ax_map.annotate("Pyramiden", (settlements['PYR']['lon'], settlements['PYR']['lat']), xytext=(settlements['PYR']['lon']+0.06, settlements['PYR']['lat']+0.02), ha="left", xycoords=transform, zorder=21)
    ax_map.annotate("Barentsburg", (settlements['BB']['lon'], settlements['BB']['lat']), xytext=(settlements['BB']['lon']+0.1, settlements['BB']['lat']-0.02),  ha="left", xycoords=transform, zorder=21)

    ax_map.text(0.98, 0.97, datetime.date.today().strftime("%d. %b %Y"), transform=ax_map.transAxes, verticalalignment='top', horizontalalignment="right", bbox=dict(boxstyle='round', facecolor='none', alpha=0., edgecolor="k"))


    return([fig, gs, ax_map, sc_map])


def initialize_fullpage_map():
    
    plt.close("all")
    
    plt.rcParams.update({"font.size": 16, "timezone": "Europe/Paris"})
    
    # create figure and subplot grid
    fig, ax_map = plt.subplots(1,1,figsize=(19,12), subplot_kw={"projection": ccrs.Mercator()})

    # plot map
    sc_map = map_classes.surface_cover_map(ax_map, ccrs.Mercator(), area='Isfjorden', mappath="/mnt/c/Svalbard_map_data/",
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
    ax_map.annotate("Pyramiden", (settlements['PYR']['lon'], settlements['PYR']['lat']), xytext=(settlements['PYR']['lon']+0.06, settlements['PYR']['lat']+0.02), ha="left", xycoords=transform, zorder=21)
    ax_map.annotate("Barentsburg", (settlements['BB']['lon'], settlements['BB']['lat']), xytext=(settlements['BB']['lon']+0.1, settlements['BB']['lat']-0.02),  ha="left", xycoords=transform, zorder=21)

    ax_map.text(0.98, 0.97, datetime.date.today().strftime("%d. %b %Y"), transform=ax_map.transAxes, verticalalignment='top', horizontalalignment="right", bbox=dict(boxstyle='round', facecolor='none', alpha=0., edgecolor="k"))

    return([fig, False, ax_map, sc_map])



def get_cbar_range(list_data, min_cbar_range):
    
    if np.any([len(i) for i in list_data]):
        total_min = np.nanmin([np.nanmin(i) for i in list_data])
        total_max = np.nanmax([np.nanmax(i) for i in list_data])
        
        if np.abs(total_max-total_min) >= min_cbar_range:
            cmin = total_min
            cmax = total_max
        else:
            cmin = total_min - (min_cbar_range - np.abs(total_max-total_min))/2
            cmax = total_max + (min_cbar_range - np.abs(total_max-total_min))/2
    else:
        cmin = 5.
        cmax = cmin + min_cbar_range
        
    return [cmin, cmax]



def plot_boat_on_map(ax_map, sc_map, boat, variable="temperature", position_switch=True, legend_switch=True, cbar_switch=True, fixed_cbar_range=[3., 6.]):

    
    cbar_labels = {'temperature': 'Temperature [°C]',
              'relative_humidity': 'Relative Humidity [%]',
              'wind_speed': 'Wind Speed [m/s]',
              'wind_direction': 'Wind Direction [°]',
              "SST": "Sea surface temperature [°C]"}
    
    sc_map.add_grid_points_scalar(boat.data['latitude'], boat.data['longitude'], boat.data[variable], fixed_cbar_range=fixed_cbar_range,
                                    cbar_switch=cbar_switch, cbar_label=cbar_labels[variable], cmap="RdYlBu_r", markersize=15, marker="o")

    # all sailing winds
    try:
        sc_map.add_grid_points_meteo_arrows(boat.data['latitude'][boat.data["mask_sailing"]][::int(10./boat.resolution_min)], boat.data['longitude'][boat.data["mask_sailing"]][::int(10./boat.resolution_min)],
                                            boat.data['u_knts'][boat.data["mask_sailing"]][::int(10./boat.resolution_min)], boat.data['v_knts'][boat.data["mask_sailing"]][::int(10./boat.resolution_min)], length=7, lw=.8)
    except (IndexError, AttributeError, ValueError):
        pass

    # latest harbour winds from each harbor, only last two, to not overload the map
    try:
        sc_map.add_grid_points_meteo_arrows([boat.data['latitude'][boat.data["mask_LYR"]][-1]], [boat.data['longitude'][boat.data["mask_LYR"]][-1]],
                                            [boat.data['u_knts'][boat.data["mask_LYR"]][-1]], [boat.data['v_knts'][boat.data["mask_LYR"]][-1]], length=7, lw=.8)
    except (IndexError, AttributeError, ValueError):
        pass

    try:
        sc_map.add_grid_points_meteo_arrows([boat.data['latitude'][boat.data["mask_PYR"]][-1]], [boat.data['longitude'][boat.data["mask_PYR"]][-1]],
                                            [boat.data['u_knts'][boat.data["mask_PYR"]][-1]], [boat.data['v_knts'][boat.data["mask_PYR"]][-1]], length=7, lw=.8)
    except (IndexError, AttributeError, ValueError):
        pass
    try:
        sc_map.add_grid_points_meteo_arrows([boat.data['latitude'][boat.data["mask_BB"]][-1]], [boat.data['longitude'][boat.data["mask_BB"]][-1]],
                                            [boat.data['u_knts'][boat.data["mask_BB"]][-1]], [boat.data['v_knts'][boat.data["mask_BB"]][-1]], length=7, lw=.8)
    except (IndexError, AttributeError, ValueError):
        pass

    # marker for ship's position
    if position_switch:
        ax_map.scatter(boat.data['longitude'][-1], boat.data['latitude'][-1], c="m", marker="X", s=150, zorder=25, transform=ccrs.PlateCarree())

    # legend
    if legend_switch:
        lc = boat.data["local_time"][-1].strftime("%H:%M")
        text_legend = mpl.lines.Line2D([],[],marker="None",linestyle="None",color="k",label=f"(upd. {lc} [local time])")

        ship_legend = mpl.lines.Line2D([],[],marker="X",markersize=13,linestyle="None",color="m",label="Last received position")

        ind = np.where((boat.data['time'] >= boat.data['time'][-1]-datetime.timedelta(hours=3)))[0]
        p = np.nanmean(boat.data['air_pressure'][ind])
        _, ptrend_bool, _, _, _, _, _, ptend, _ = mk.original_test(boat.data['air_pressure'][ind])
        if ptrend_bool:
            if ptend > 0:
                pressure_legend = mpl.lines.Line2D([],[],marker=10,markersize=12,linestyle="None",color="lime",label=f"Pressure: {p:.1f} hPa")
            else:
                pressure_legend = mpl.lines.Line2D([],[],marker=11,markersize=12,linestyle="None",color="red",label=f"Pressure: {p:.1f} hPa")
        else:
            pressure_legend = mpl.lines.Line2D([],[],linewidth=3,linestyle="-",color="k",label=f"Pressure: {p:.1f} hPa")

        handles = [pressure_legend, ship_legend, text_legend]
        ax_map.legend(handles=handles, loc="upper left")
    
    return



def combined_legend_positions(ax_map, boat, boat_names):
    
    markers = {"MSBard": "X", "MSPolargirl": "*", "MSBillefjord": "P", "MSBerg": "X"}
    colors = {"MSBard": "m", "MSPolargirl": "g", "MSBillefjord": "purple", "MSBerg": "m"}
    sizes = {"MSBard": 150, "MSPolargirl": 200, "MSBillefjord": 150, "MSBerg": 150}
    
    i = list(boat.keys())[0]
    ind = np.where((boat[i].data['time'] >= boat[i].data['time'][-1]-datetime.timedelta(hours=3)))[0]
    p = np.nanmean(boat[i].data['air_pressure'][ind])
    _, ptrend_bool, _, _, _, _, _, ptend, _ = mk.original_test(boat[i].data['air_pressure'][ind])
    if ptrend_bool:
        if ptend > 0:
            pressure_legend = mpl.lines.Line2D([],[],marker=10,markersize=12,linestyle="None",color="lime",label=f"Pressure: {p:.1f} hPa")
        else:
            pressure_legend = mpl.lines.Line2D([],[],marker=11,markersize=12,linestyle="None",color="red",label=f"Pressure: {p:.1f} hPa")
    else:
        pressure_legend = mpl.lines.Line2D([],[],linewidth=3,linestyle="-",color="k",label=f"Pressure: {p:.1f} hPa")
        
    handles = [pressure_legend]
        
    for i in boat.keys():
        # marker for ship's position
        ax_map.scatter(boat[i].data['longitude'][-1], boat[i].data['latitude'][-1], c=colors[i], marker=markers[i], s=sizes[i], zorder=25, transform=ccrs.PlateCarree())

        # legend
        lc = boat[i].data["local_time"][-1].strftime("%H:%M")
        text_legend = mpl.lines.Line2D([],[], marker="None", linestyle="None", color=colors[i], label=f"(upd. {lc} [local time])")
    
        ship_legend = mpl.lines.Line2D([],[], marker=markers[i], markersize=sizes[i]/10., linestyle="None", color=colors[i], label=boat_names[i]["name"])
        
        handles.append(ship_legend)
        handles.append(text_legend)
        
    ax_map.legend(handles=handles, loc="upper left")
    
    return
    
    



def plot_lighthouse_on_map(lighthouse, ax_map, sc_map):
    
    transform = ccrs.PlateCarree()._as_mpl_transform(ax_map)
    
    ax_map.annotate("\u265C", (lighthouse.longitude, lighthouse.latitude), xytext=(lighthouse.longitude, lighthouse.latitude), ha="center", va="center", xycoords=transform, zorder=21)
    props = dict(boxstyle='round', facecolor='white', alpha=0.6)
    data_str = '{t:<5} {a:.1f}°C\n{r:<3} {b:.1f}%\n{w:<3} {c:.1f}m/s\n{g:<3} {d}\n(upd. {l})'.format(t="T", r="RH", w="WS", g="WD", l=lighthouse.data["local_time"][-1].strftime("%H:%M"),
                                                                                             a=lighthouse.data['temperature'][-1],
                                                                                             b=lighthouse.data['relative_humidity'][-1],
                                                                                             c=lighthouse.data['wind_speed'][-1],
                                                                                             d=wdir_classes[lighthouse.data['wind_sector'][-1]][0])

    if lighthouse.station_name == "Bohemanneset":
        x_shift = -0.08
        y_shift = 0.02
        if wdir_classes[lighthouse.data['wind_sector'][-1]][0] in ["N", "NW", "W"]:
            x_shift = -0.3
            y_shift = 0.06
        ax_map.annotate(data_str, xy=(lighthouse.longitude, lighthouse.latitude), bbox=props, ma='left', ha="right", va="bottom",
                        xytext=(lighthouse.longitude+x_shift, lighthouse.latitude+y_shift), xycoords=transform, zorder=21)
    elif lighthouse.station_name == "Gasoyane":
        x_shift = 0.1
        y_shift = 0.00
        if wdir_classes[lighthouse.data['wind_sector'][-1]][0] in ["N", "NW", "W"]:
            x_shift = 0.2
            y_shift = 0.00
        ax_map.annotate(data_str, xy=(lighthouse.longitude, lighthouse.latitude), bbox=props, ma='left', ha="left", va="center",
                        xytext=(lighthouse.longitude+x_shift, lighthouse.latitude+y_shift), xycoords=transform, zorder=21)
    elif lighthouse.station_name == "Narveneset":
        x_shift = -0.15
        y_shift = 0.05
        ax_map.annotate(data_str, xy=(lighthouse.longitude, lighthouse.latitude), bbox=props, ma='left', ha="right", va="bottom",
                        xytext=(lighthouse.longitude+x_shift, lighthouse.latitude), xycoords=transform, zorder=21)
    elif lighthouse.station_name == "Daudmannsodden":
        x_shift = -0.08
        y_shift = 0.02
        if wdir_classes[lighthouse.data['wind_sector'][-1]][0] in ["N", "NW", "W"]:
            x_shift = -0.3
            y_shift = 0.06
        ax_map.annotate(data_str, xy=(lighthouse.longitude, lighthouse.latitude), bbox=props, ma='left', ha="right", va="bottom",
                        xytext=(lighthouse.longitude+x_shift, lighthouse.latitude+y_shift), xycoords=transform, zorder=21)
    elif lighthouse.station_name == "Kapp Thordsen":
        x_shift = -0.08
        y_shift = 0.02
        if wdir_classes[lighthouse.data['wind_sector'][-1]][0] in ["N", "NW", "W"]:
            x_shift = -0.2
            y_shift = 0.04
        ax_map.annotate(data_str, xy=(lighthouse.longitude, lighthouse.latitude), bbox=props, ma='left', ha="right", va="bottom",
                        xytext=(lighthouse.longitude+x_shift, lighthouse.latitude+y_shift), xycoords=transform, zorder=21)
            

    sc_map.add_grid_points_meteo_arrows([lighthouse.latitude], [lighthouse.longitude],
                                        [lighthouse.data['u_knts'][-1]], [lighthouse.data['v_knts'][-1]], length=7, lw=.8, zorder=40)

    return


def plot_MET_station_on_map(station, data, ax_map, sc_map):
    
    station_names = {"IR": {"name": "Isfjord Radio", "ID": "SN99790", "height": 7., "lat": 78.0625, "lon": 13.6192},
                     "LYR": {"name": "Longyearbyen Airport", "ID": "SN99840", "height": 28., "lat": 78.2453, "lon": 15.5015},
                     "PYR": {"name": "Pyramiden", "ID": "SN99880", "height": 20., "lat": 78.6557, "lon": 16.3603},
                     "NS": {"name": "Nedre Sassendalen", "ID": "SN99882", "height": 13., "lat": 78.3313, "lon": 16.6818},
                     "PB": {"name": "Plataberget", "ID": "SN99843", "height": 450., "lat": 78.2278, "lon": 15.378},
                     "AD": {"name": "Adventdalen", "ID": "SN99870", "height": 15., "lat":  78.2022, "lon": 15.831},
                     "IH": {"name": "Innerhytta", "ID": "SN99879", "height": 81., "lat":  78.18883, "lon": 16.34423},
                     "AO": {"name": "Akseloya", "ID": "SN99765", "height": 20., "lat":  77.6873, "lon": 14.7578},
                     "KL": {"name": "Klauva", "ID": "SN99884", "height": 480., "lat":  78.3002, "lon": 18.2248},
                     "ID": {"name": "Istjorndalen", "ID": "SN99770", "height": 188., "lat":  78.0092, "lon": 15.2108}, 
                     "JH": {"name": "Janssonhagen", "ID": "SN99874", "height": 250., "lat":  78.18, "lon": 16.41},
                     "RP": {"name": "Reindalspasset", "ID": "SN99763", "height": 181., "lat":  78.0648, "lon":  17.0442},
                     "SV": {"name": "Svea", "ID": "SN99760", "height": 9., "lat":  77.8953, "lon": 16.72}}
    
    transform = ccrs.PlateCarree()._as_mpl_transform(ax_map)
    
    if station not in ["LYR", "PYR"]:
        ax_map.annotate("\u265C", (station_names[station]["lon"], station_names[station]["lat"]), xytext=(station_names[station]["lon"], station_names[station]["lat"]), ha="center", va="center", xycoords=transform, zorder=21)
    
    
    props = dict(boxstyle='round', facecolor='white', alpha=0.6)
    data_str = '{t:<5} {a:.1f}°C\n{r:<3} {b:.1f}%\n{w:<3} {c:.1f}m/s\n{g:<3} {d}\n(upd. {l})'.format(t="T", r="RH", w="WS", g="WD", l=data.index[-1].strftime("%H:%M"),
                                                                                             a=data['temperature'].iloc[-1],
                                                                                             b=data['relative_humidity'].iloc[-1],
                                                                                             c=data['wind_speed'].iloc[-1],
                                                                                             d=wdir_classes[data['wind_sector'].iloc[-1]][0])


    if station == "LYR":
        x_shift = 0.0
        y_shift = -0.08
        # if wdir_classes[data['wind_sector'].iloc[-1]][0] in ["S", "SW", "W"]:
            # x_shift = -0.15
            # y_shift = -0.03
        ax_map.annotate(data_str, xy=(station_names[station]["lon"], station_names[station]["lat"]), bbox=props, ma='left', ha="center", va="top",
                        xytext=(station_names[station]["lon"]+x_shift, station_names[station]["lat"]+y_shift), xycoords=transform, zorder=21)
    elif station == "IR":
        x_shift = -0.08
        y_shift = -0.0
        if wdir_classes[data['wind_sector'].iloc[-1]][0] in ["S", "SW", "W"]:
            x_shift = -0.3
            # y_shift = -0.06
        ax_map.annotate(data_str, xy=(station_names[station]["lon"], station_names[station]["lat"]), bbox=props, ma='left', ha="right", va="top",
                        xytext=(station_names[station]["lon"]+x_shift, station_names[station]["lat"]+y_shift), xycoords=transform, zorder=21)
    elif station == "NS":
        x_shift = 0.08
        y_shift = -0.02
        if wdir_classes[data['wind_sector'].iloc[-1]][0] in ["S", "SE", "E"]:
            x_shift = 0.3
            y_shift = -0.06
        ax_map.annotate(data_str, xy=(station_names[station]["lon"], station_names[station]["lat"]), bbox=props, ma='left', ha="left", va="top",
                        xytext=(station_names[station]["lon"]+x_shift, station_names[station]["lat"]+y_shift), xycoords=transform, zorder=21)
    elif station == "PYR":
        x_shift = 0.0
        y_shift = 0.06
        # if wdir_classes[data['wind_sector'].iloc[-1]][0] in ["N", "NW", "W"]:
        #     x_shift = -0.3
        #     y_shift = 0.06
        ax_map.annotate(data_str, xy=(station_names[station]["lon"], station_names[station]["lat"]), bbox=props, ma='left', ha="left", va="bottom",
                        xytext=(station_names[station]["lon"]+x_shift, station_names[station]["lat"]+y_shift), xycoords=transform, zorder=21)
            

    sc_map.add_grid_points_meteo_arrows([station_names[station]["lat"]], [station_names[station]["lon"]],
                                        [data['u_knts'].iloc[-1]], [data['v_knts'].iloc[-1]], length=7, lw=.8, zorder=40)

    return





def plot_boat_timeseries(boat, fig, gs, status):
    
    future_time_plot = {'live': 4, 'overview': 0}

    plot_variables = {0: 'temperature', 1: 'relative_humidity', 2: 'wind_speed', 3: 'wind_direction'}

    plot_colors = {'temperature': '#377eb8',
              'relative_humidity': '#4daf4a',
              'wind_speed': '#984ea3',
              'wind_direction': '#e41a1c'}

    plot_labels = {'temperature': 'Temp. [°C]',
              'relative_humidity': 'Rel. Humi. [%]',
              'wind_speed': 'Wind Speed [m/s]',
              'wind_direction': 'Wind Dir. [°]'}

    plot_units = {'temperature': '°C',
              'relative_humidity': '%',
              'wind_speed': 'm/s',
              'wind_direction': '°'}



    # initialize complete wind direction plot, to get the x axis for the following subplots
    ax_4 = fig.add_subplot(gs[3,0])
    ax_4.scatter(boat.data['local_time'], boat.data["wind_direction"], s=3, color=plot_colors['wind_direction'])



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
            
            ax.fill_between(boat.data['local_time'], ymin, ymax, where = boat.data['mask_BB'] == True, facecolor="none", linewidth=2.0, edgecolor=plot_colors[plot_variables[i]], hatch="/", alpha=0.5)
            ax.fill_between(boat.data['local_time'], ymin, ymax, where = boat.data['mask_PYR'] == True, facecolor="none", linewidth=2.0, edgecolor=plot_colors[plot_variables[i]], hatch='\\', alpha=0.5)
            
            ax.yaxis.set_major_locator(plt.MaxNLocator(5))
        # everything that goes only into the last subplot (wind direction)
        else:
            ax.set(xlabel = 'Local time', ylabel = plot_labels['wind_direction'])
            ax.set_xlim((boat.data['local_time'][0], boat.data['local_time'][-1]+datetime.timedelta(hours=future_time_plot[status])))
            ax.set_ylim((0., 360.))
            ax.fill_between(boat.data['local_time'], 0., 360., where = boat.data['mask_LYR'] == True, facecolor=plot_colors[plot_variables[i]], alpha=0.2)
            
            ax.fill_between(boat.data['local_time'], 0., 360., where = boat.data['mask_BB'] == True, facecolor="none", linewidth=2.0, edgecolor=plot_colors[plot_variables[i]], hatch="/", alpha=0.5)
            ax.fill_between(boat.data['local_time'], 0., 360., where = boat.data['mask_PYR'] == True, facecolor="none", linewidth=2.0, edgecolor=plot_colors[plot_variables[i]], hatch='\\', alpha=0.5)

            ax.set_yticks([0., 90., 180., 270., 360.])
            ax.set_yticklabels(["N", "E", "S", "W", "N"])
            ax.xaxis.set_major_locator(HourLocator(interval=2))
            ax.xaxis.set_major_formatter(DateFormatter('%H'))



        if status == "live":
            if i < len(axes)-1:
                latest_value_label = "{a:.1f} {p}".format(a=boat.data[plot_variables[i]][-1], p=plot_units[plot_variables[i]])
            else:
                latest_value_label = wdir_classes[boat.data['wind_sector'][-1]][0]
            ax.scatter(boat.data['local_time'][-1], boat.data[plot_variables[i]][-1], s=50, color=plot_colors[plot_variables[i]])
            ax.annotate(latest_value_label, (boat.data['local_time'][-1], boat.data[plot_variables[i]][-1]),
                        xytext=(10,0), textcoords="offset points", va="center", bbox=dict(boxstyle="round", alpha=0.1))

        # plot markers for BB and PYR (if possible)
        try:
            ax.scatter(boat.data['local_time'][boat.data['mask_BB']][np.sum(boat.data['mask_BB'])//2], boat.data[plot_variables[i]][boat.data['mask_BB']][np.sum(boat.data['mask_BB'])//2],
                        color=plot_colors[plot_variables[i]], marker='d', s=250)
        except IndexError:
            pass
        try:
            ax.scatter(boat.data['local_time'][boat.data['mask_PYR']][np.sum(boat.data['mask_PYR'])//2], boat.data[plot_variables[i]][boat.data['mask_PYR']][np.sum(boat.data['mask_PYR'])//2],
                        color=plot_colors[plot_variables[i]], marker='^', s=250)
        except IndexError:
            pass

        ax.grid()
        
    return



def plot_lighthouse_timeseries(lighthouse, status):
    
    plt.rcParams.update({"font.size": 12, "timezone": "Europe/Paris"})
    
    future_time_plot = {'live': 4, 'overview': 0}

    plot_variables = {0: 'temperature', 1: 'relative_humidity', 2: 'wind_speed', 3: 'wind_direction'}
    # plot_variables = {0: 'battery', 1: 'relative_humidity', 2: 'wind_speed', 3: 'wind_direction'}


    plot_colors = {'temperature': '#377eb8',
              'relative_humidity': '#4daf4a',
              'wind_speed': '#984ea3',
              'wind_direction': '#e41a1c',
              'battery': 'k'}

    plot_labels = {'temperature': 'Temp.\n[°C]',
              'relative_humidity': 'Rel. Humi.\n[%]',
              'wind_speed': 'Wind Speed\n[m/s]',
              'wind_direction': 'Wind Dir.\n[°]',
              'battery': 'k'}

    plot_units = {'temperature': '°C',
              'relative_humidity': '%',
              'wind_speed': 'm/s',
              'wind_direction': '°',
              'battery': "V"}
    

    fig, axes = plt.subplots(nrows=4, ncols=len(lighthouse.keys()), sharex=True, sharey="row", constrained_layout=True)
    for l_ind, l in enumerate(lighthouse.keys()):
        if len(lighthouse.keys()) > 1:
            axes_lighthouse = axes[:,l_ind]
        else:
            axes_lighthouse = axes
            
        axes_lighthouse[0].set_title(lighthouse[l].station_name)
            
        for i, ax in enumerate(axes_lighthouse):
            if i == len(plot_variables.keys())-1:
                ax.scatter(lighthouse[l].data['local_time'], lighthouse[l].data["wind_direction"], s=3, color=plot_colors['wind_direction'])
                ax.set_xlim((lighthouse[l].data['local_time'][0], lighthouse[l].data['local_time'][-1]+datetime.timedelta(hours=future_time_plot[status])))
                ax.set_ylim((0., 360.))
                ax.set_yticks([0., 90., 180., 270., 360.])
                ax.set_yticklabels(["N", "E", "S", "W", "N"])
                ax.xaxis.set_major_locator(HourLocator(interval=2))
                ax.xaxis.set_major_formatter(DateFormatter('%H'))
                if l_ind == 0:
                    ax.set_ylabel(plot_labels['wind_direction'])
                ax.set_xlabel("Local time")
            else:
                ax.plot(lighthouse[l].data['local_time'], lighthouse[l].data[plot_variables[i]], color=plot_colors[plot_variables[i]])
                if l_ind == 0:
                    ax.set_ylabel(plot_labels[plot_variables[i]])
                
            if status == "live":
                if i < len(axes_lighthouse)-1:
                    latest_value_label = "{a:.1f} {p}".format(a=lighthouse[l].data[plot_variables[i]][-1], p=plot_units[plot_variables[i]])
                else:
                    latest_value_label = wdir_classes[lighthouse[l].data['wind_sector'][-1]][0]
                ax.scatter(lighthouse[l].data['local_time'][-1], lighthouse[l].data[plot_variables[i]][-1], s=50, color=plot_colors[plot_variables[i]])
                ax.annotate(latest_value_label, (lighthouse[l].data['local_time'][-1], lighthouse[l].data[plot_variables[i]][-1]),
                            xytext=(10,0), textcoords="offset points", va="center", bbox=dict(boxstyle="round", alpha=0.1))
            ax.grid()

            
    return
                    
    