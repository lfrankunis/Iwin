# -*- coding: utf-8 -*-
"""
Created on Sat Mar 27 14:37:44 2021

@author: unismet
"""

import os
import datetime
from IWIN_plotting_functions import initialize_halfpage_map, initialize_fullpage_map, \
    plot_boat_on_map, plot_lighthouse_on_map, plot_MET_station_on_map, \
        plot_boat_timeseries, get_cbar_range, combined_legend_positions, plot_lighthouse_timeseries
from MET_stations import download_MET_stations
import matplotlib.pyplot as plt
import mobile_AWS_class
import lighthouse_AWS_class
import ftplib
import yaml
from IWIN_structure_data_functions import create_GIS_input_file

import multiprocessing as mp



def update_all_plots(boats_to_plot, lighthouses_to_plot, MET_stations_to_plot):
    
    # define constants
    with open("./config_paths.yaml", "r", encoding='utf-8') as f:
        paths = yaml.safe_load(f)
        
    update_time = datetime.datetime.now().replace(second=0, microsecond=0)
    
    with open("./config_metadata.yaml", "r", encoding='utf-8') as f:
        station_metadata = yaml.safe_load(f)
        
    lighthouses = station_metadata["lighthouse_stations"]


    boats = station_metadata["mobile_stations"]
    
    
    with open("./config_sensors.yaml", "r", encoding='utf-8') as f:
        sensors = yaml.safe_load(f)
    

    status = "live"

    map_vari = "temperature"
    min_cbar_range = 3.
    
    
    
    
    ####################################################################################
    ####################################################################################
    
    # define what data time period to load
    latest_update_time = datetime.datetime.strftime(update_time, "%Y%m%d%H%M")
    period_to_load = datetime.timedelta(days=1)
    start_time = (update_time - period_to_load).strftime("%Y%m%d%H%M")


    ####################################################################################
    ####################################################################################
    
    # load boat data
    boat = {}
    for b in boats_to_plot:
        boat[b] = mobile_AWS_class.mobile_AWS(station=b, resolution="1min",
                                            starttime=start_time, endtime=latest_update_time,
                                            variables=['temperature', 'air_pressure', 'relative_humidity', 'wind_speed', 'wind_direction',
                                                       "wind_speed_raw", "wind_direction_raw", 'latitude', 'longitude', "GPS_speed", "GPS_heading", "compass_heading"],
                                            file_type="raw", path=paths["local_data"])
        
        boat[b].only_latest_data(datetime.timedelta(hours=12))
        boat[b].filter_GPScoverage()
        boat[b].masks_for_harbors()
        boat[b].correct_winds()
        boat[b].calculate_windvector_components(corrected=True)
        boat[b].calculate_wind_sector(corrected=True)
        boat[b].calculate_wind_in_knots(corrected=True)
        
        # add info to the class instance
        boat[b].boat_name = boats[b]["name"]
        
    
    ####################################################################################
    ####################################################################################
    
    # load lighthouse data
    lighthouse = {}
    for l in lighthouses_to_plot:
        lighthouse[l] = lighthouse_AWS_class.lighthouse_AWS(station=l, resolution="1min",
                                            starttime=start_time, endtime=latest_update_time,
                                            variables=['temperature', 'air_pressure', 'relative_humidity', 'wind_speed', 'wind_direction'],
                                            file_type="raw", path=paths["local_data"])
    
        lighthouse[l].calculate_windvector_components()
        lighthouse[l].calculate_wind_sector()
        lighthouse[l].calculate_wind_in_knots()
        lighthouse[l].only_latest_data(datetime.timedelta(hours=12))
        
        # add info to the class instance
        lighthouse[l].station_name = lighthouses[l]["name"]
        lighthouse[l].latitude = lighthouses[l]["lat"]
        lighthouse[l].longitude = lighthouses[l]["lon"]

    ####################################################################################
    ####################################################################################
    
    # load MET station data
    
    met_stations = download_MET_stations(update_time)
    
    ####################################################################################
    ####################################################################################
    

    # create and upload individual boat plots
    for b in boats_to_plot:
        
        fig, gs, ax_map, sc_map = initialize_halfpage_map()
        cbar_range = get_cbar_range([boat[b].data[map_vari]], min_cbar_range)
        plot_boat_on_map(ax_map, sc_map, boat[b], variable=map_vari, position_switch=True, legend_switch=True, cbar_switch=True, fixed_cbar_range=cbar_range)
        for l in lighthouses_to_plot:
            plot_lighthouse_on_map(lighthouse[l], ax_map, sc_map)
        for m in MET_stations_to_plot:
            plot_MET_station_on_map(m, met_stations[m], ax_map, sc_map)
        plot_boat_timeseries(boat[b], fig, gs, status)
        
        local_output_path = f"{paths['local_desktop']}liveplot_{b}.png"
                
        plt.savefig(local_output_path)
        plt.close("all")
        
        upload_picture(local_output_path, os.path.basename(local_output_path))
        
        
    # create and upload overview plot
    fig, gs, ax_map, sc_map = initialize_fullpage_map()
    cbar_range = get_cbar_range([boat[i].data[map_vari] for i in boats_to_plot], min_cbar_range)
    for i, b in enumerate(boats_to_plot):
        if i == 0:
            plot_boat_on_map(ax_map, sc_map, boat[b], variable=map_vari, position_switch=False, legend_switch=False, cbar_switch=True, fixed_cbar_range=cbar_range)
        else:
            plot_boat_on_map(ax_map, sc_map, boat[b], variable=map_vari, position_switch=False, legend_switch=False, cbar_switch=False, fixed_cbar_range=cbar_range)
    for l in lighthouses_to_plot:
        plot_lighthouse_on_map(lighthouse[l], ax_map, sc_map)
    for m in MET_stations_to_plot:
        plot_MET_station_on_map(m, met_stations[m], ax_map, sc_map)
    if len(boats_to_plot) > 0:
        combined_legend_positions(ax_map, boat, boats) # combined legend
    
    local_output_path = f"{paths['local_desktop']}liveplot_overview_map.png"
    
    plt.savefig(local_output_path)
    plt.close("all")
        
    upload_picture(local_output_path, os.path.basename(local_output_path))
    # create_GIS_input_file(boat, lighthouse, met_stations, past_hours=3, path=paths['onedrive'])


    # create and upload lighthouse plot
    plot_lighthouse_timeseries(lighthouse, status)

    local_output_path = f"{paths['local_desktop']}liveplot_lighthouses.png"
    
    plt.savefig(local_output_path)
    plt.close("all")
        
    upload_picture(local_output_path, os.path.basename(local_output_path))
        
    return



def upload_picture(local_file, online_file_name):
    
    with ftplib.FTP("unisacsi.com", "u834143596", "unis12345") as session:
        with open(local_file, "rb") as f:
            session.storbinary(f"STOR {online_file_name}", f)

    return


####################################################################################
####################################################################################



if __name__ == '__main__':
    
    # read config file
    with open("./config_paths.yaml", "r", encoding='utf-8') as f:
        paths = yaml.safe_load(f)
    
    with open("./config_processing_plotting.yaml", "r", encoding='utf-8') as f:
        stations_to_plot = yaml.safe_load(f)
    
    
    for k in stations_to_plot.keys():
        if stations_to_plot[k] is None:
            stations_to_plot[k] = []

    
    MET_switches = {"LYR": True,
                    "IR":  True,
                    "PYR": True,
                    "NS":  True}
    
    
    ###########################################################################
    ###########################################################################
    

    MET_stations = [s for s, sw in MET_switches.items() if sw]

    proc=mp.Process(target=update_all_plots, args=[stations_to_plot["mobile_stations"], stations_to_plot["lighthouse_stations"], MET_stations])
    proc.start()
    proc.join()