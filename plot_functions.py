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
    sc_map = map_classes.surface_cover_map(ax_map, ccrs.Mercator(), area='Isfjorden', mappath="C:/Svalbard_map_data/",
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
    ax_map.annotate("Barentsburg", (settlements['BB']['lon'], settlements['BB']['lat']), xytext=(settlements['BB']['lon']+0.1, settlements['BB']['lat']-0.02),  ha="left", xycoords=transform, zorder=21)


    return([fig, gs, ax_map, sc_map])


def initialize_fullpage_map():
    
    plt.close("all")
    
    plt.rcParams.update({"font.size": 16, "timezone": "Europe/Paris"})
    
    # create figure and subplot grid
    fig, ax_map = plt.subplots(1,1,figsize=(19,12), subplot_kw={"projection": ccrs.Mercator()})

    # plot map
    sc_map = map_classes.surface_cover_map(ax_map, ccrs.Mercator(), area='Isfjorden', mappath="C:/Svalbard_map_data/",
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
    ax_map.annotate("Barentsburg", (settlements['BB']['lon'], settlements['BB']['lat']), xytext=(settlements['BB']['lon']+0.1, settlements['BB']['lat']-0.02),  ha="left", xycoords=transform, zorder=21)


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
              'wind_direction': 'Wind Direction [°]'}
    
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
    
    return



def combined_legend_positions(ax_map, boat, boat_names):
    
    markers = {1: "X", 2: "*"}
    colors = {1: "m", 2: "g"}
    sizes = {1: 150, 2: 200}
    
    i = list(boat.keys())[0]
    ind = np.where((boat[i].data['time'] >= boat[i].data['time'][-1]-datetime.timedelta(hours=3)))[0]
    p = np.nanmean(boat[i].data['pressure'][ind])
    _, ptrend_bool, _, _, _, _, _, ptend, _ = mk.original_test(boat[i].data['pressure'][ind])
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
    
        ship_legend = mpl.lines.Line2D([],[], marker=markers[i], markersize=11+2*i, linestyle="None", color=colors[i], label=boat_names[i])
        
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

    x_shift = -0.08
    y_shift = 0.02
    
    if lighthouse.station_name == "Bohemanneset":
        if wdir_classes[lighthouse.data['wind_sector'][-1]][0] in ["N", "NW", "W"]:
            x_shift = -0.3
            y_shift = 0.06
            
            
    ax_map.annotate(data_str, xy=(lighthouse.longitude, lighthouse.latitude), bbox=props, ma='left', ha="right", va="bottom",
                    xytext=(lighthouse.longitude+x_shift, lighthouse.latitude+y_shift), xycoords=transform, zorder=21)

    sc_map.add_grid_points_meteo_arrows([lighthouse.latitude], [lighthouse.longitude],
                                        [lighthouse.data['u_knts'][-1]], [lighthouse.data['v_knts'][-1]], length=7, lw=.8, zorder=40)

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
