# -*- coding: utf-8 -*-
"""
Created on Sat Mar 27 14:37:44 2021

@author: unismet
"""

import os
import time
import datetime
from plot_functions import initialize_halfpage_map, initialize_fullpage_map, plot_boat_on_map, plot_lighthouse_on_map, plot_MET_station_on_map, plot_boat_timeseries, get_cbar_range, combined_legend_positions
from MET_stations import download_MET_stations
import matplotlib.pyplot as plt
import mobile_AWS_class
import lighthouse_AWS_class
import ftplib
from AWS_structure_data_functions import create_GIS_input_file

import multiprocessing as mp


def next_wakeup():
    
    # refresh period
    dt_minutes = 2
    dt_hours = 0
    time_delta=datetime.timedelta(hours=dt_hours, minutes=dt_minutes)
    round_to = time_delta.total_seconds()                   # 60s

    dt = datetime.datetime.now()
    seconds = (dt - dt.min).seconds                       

    if seconds % round_to == 0 and dt.microsecond == 0:
        rounding = (seconds + round_to / 2) // round_to * round_to
    else:
        rounding = (seconds + dt.microsecond/1000000 + round_to) // round_to * round_to

    next_wakeup_time = dt + datetime.timedelta(0, rounding - seconds, - dt.microsecond) + datetime.timedelta(seconds=30)
    print("The next wakeup is scheduled for: {a}".format(a=next_wakeup_time))

    return next_wakeup_time



def update_all_plots(update_time):
    
    # define constants
    path_data = "C:/Data/"


    lighthouses = {1885: {"name": "Bohemanneset", 'lat': 78.38166, 'lon': 14.75300},
                   1886: {"name": "Gasoyane", 'lat': 78.45792,'lon': 16.20082},
                   #1886: {"name": "Kapp Thordsen", 'lat': 78.45638,'lon': 15.46793},
                   1884: {"name": "Narveneset", 'lat': 78.56343,'lon': 16.29687},
                   1887: {"name": "Daudmannsodden", 'lat': 78.21056,'lon': 12.98685}}
    lighthouses_to_plot = [1884, 1885]

    boat_names = {1883: "MS_Bard", 1872: "MS_Polargirl", 1924: "MS_Billefjord"}
    boats_to_plot = [1883, 1924, 1872]
    
    MET_stations_to_plot = ["LYR", "IR", "PYR", "NS"]

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
                                            variables=['temperature', 'pressure', 'relative_humidity', 'wind_speed', 'wind_direction',
                                                       "wind_speed_raw", "wind_direction_raw", 'latitude', 'longitude', "GPS_speed", "GPS_heading", "compass_heading"],
                                            file_type="raw", path=path_data)
        
        boat[b].only_latest_data(datetime.timedelta(hours=12))
        boat[b].filter_GPScoverage()
        boat[b].masks_for_harbors()
        boat[b].correct_winds()
        boat[b].calculate_windvector_components(corrected=True)
        boat[b].calculate_wind_sector(corrected=True)
        boat[b].calculate_wind_in_knots(corrected=True)
        
        # add info to the class instance
        boat[b].boat_name = boat_names[b]
    
    ####################################################################################
    ####################################################################################
    
    # load lighthouse data
    lighthouse = {}
    for l in lighthouses_to_plot:
        lighthouse[l] = lighthouse_AWS_class.lighthouse_AWS(station=l, resolution="1min",
                                            starttime=start_time, endtime=latest_update_time,
                                            variables=['temperature', 'pressure', 'relative_humidity', 'wind_speed', 'wind_direction'],
                                            file_type="raw", path=path_data)
    
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
    
    met_stations = download_MET_stations(update_time, MET_stations_to_plot)
    
    ####################################################################################
    ####################################################################################
    

    # create and upload plots
    for b in boats_to_plot:
        
        fig, gs, ax_map, sc_map = initialize_halfpage_map()
        cbar_range = get_cbar_range([boat[b].data[map_vari]], min_cbar_range)
        plot_boat_on_map(ax_map, sc_map, boat[b], variable=map_vari, position_switch=True, legend_switch=True, cbar_switch=True, fixed_cbar_range=cbar_range)
        for l in lighthouses_to_plot:
            plot_lighthouse_on_map(lighthouse[l], ax_map, sc_map)
        for m in MET_stations_to_plot:
            plot_MET_station_on_map(m, met_stations[m], ax_map, sc_map)
        plot_boat_timeseries(boat[b], fig, gs, status)
        
        local_output_path = "C:/Users/unismet/Desktop/liveplot_{b}.png".format(b=boat[b].boat_name)
        
        plt.savefig(local_output_path)
        plt.close("all")
        
        upload_picture(local_output_path, os.path.basename(local_output_path))
        
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
        combined_legend_positions(ax_map, boat, boat_names) # combined legend
    
    local_output_path = "C:/Users/unismet/Desktop/liveplot_overview_map.png"
    
    plt.savefig(local_output_path)
    plt.close("all")
        
    upload_picture(local_output_path, os.path.basename(local_output_path))
    
    
    create_GIS_input_file(boat, lighthouse, met_stations, past_hours=3)
    
    
        
    return



def upload_picture(local_file, online_file_name):
    
    with ftplib.FTP("unisacsi.com", "u834143596", "unis12345") as session:
        with open(local_file, "rb") as f:
            session.storbinary(f"STOR {online_file_name}", f)

    return


####################################################################################
####################################################################################



if __name__ == '__main__':

    next_wakeup_time = next_wakeup()    # start first wake-up
    
    # always true, to keep the script running forever
    while True:                 
        while datetime.datetime.now() < next_wakeup_time:       # sleep until the next scheduled wakeup time
            time.sleep(1)
            
        proc=mp.Process(target=update_all_plots, args=[next_wakeup_time])
        proc.daemon=True
        proc.start()
        proc.join()
        
        next_wakeup_time = next_wakeup()           # don't forget to update the next wakeup time!!!