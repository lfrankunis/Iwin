# -*- coding: utf-8 -*-
"""
Created on Sat Mar 27 14:37:44 2021

@author: unismet
"""

import os
import time
import datetime
from plot_functions import initialize_halfpage_map, plot_boat_on_map, plot_lighthouse_on_map, plot_boat_timeseries, get_cbar_range, plot_lighthouse_timeseries
import matplotlib.pyplot as plt
import numpy as np
import mobile_AWS_class
import lighthouse_AWS_class
import ftplib
import yaml

import multiprocessing as mp


def next_wakeup():
    
    # refresh period
    dt_days = 1
    dt_minutes_offset = 15      # to shift the python precessing until the Loggernet data downloading is completed 
    
    time_delta=datetime.timedelta(days=dt_days)
    round_to = time_delta.total_seconds()

    dt = datetime.datetime.now()
    seconds = (dt - dt.min).seconds                       

    if seconds % round_to == 0 and dt.microsecond == 0:
        rounding = (seconds + round_to / 2) // round_to * round_to
    else:
        rounding = (seconds + dt.microsecond/1000000 + round_to) // round_to * round_to

    next_wakeup_time = dt + datetime.timedelta(0, rounding - seconds, - dt.microsecond)
    next_wakeup_time += datetime.timedelta(minutes=dt_minutes_offset)
    print("The next wakeup is scheduled for: {a}".format(a=next_wakeup_time))

    return next_wakeup_time



def update_overview_plot(update_time, file_type="raw", for_website=True):
    

    # define constants
    with open("./path_config.yaml", "r") as f:
        paths = yaml.safe_load(f)
    
    
    lighthouses = {1885: {"name": "Bohemanneset", 'lat': 78.38166, 'lon': 14.75300},
                   1887: {"name": "Gasoyane", 'lat': 78.45792,'lon': 16.20082},
                   1884: {"name": "Narveneset", 'lat': 78.56343,'lon': 16.29687},
                   1886: {"name": "Daudmannsodden", 'lat': 78.21056,'lon': 12.98685}}
    lighthouses_to_plot = [1884, 1885, 1886, 1887]
    
    boat_names = {1883: "MS_Bard", 1872: "MS_Polargirl", 1924 : "MS_Bard_2"}
    boats_to_plot = []
    
    status = "overview"
    
    map_vari = "temperature"
    min_cbar_range = 3.
    
    dt_minutes_offset = 15      # to shift the python precessing until the Loggernet data downloading is completed 

    ####################################################################################
    
    # define what data time period to load
    latest_update_time = datetime.datetime.strftime(update_time - datetime.timedelta(minutes=dt_minutes_offset), "%Y%m%d%H%M")
    period_to_load = datetime.timedelta(days=1)
    start_time = (update_time - datetime.timedelta(minutes=dt_minutes_offset) - period_to_load).strftime("%Y%m%d%H%M")

    ####################################################################################

    # load boat data
    boat = {}
    for b in boats_to_plot:
        boat[b] = mobile_AWS_class.mobile_AWS(station=b, resolution="1min",
                                            starttime=start_time, endtime=latest_update_time,
                                            variables=['temperature', 'air_pressure', 'relative_humidity', 'wind_speed', 'wind_direction', 'wind_speed_raw',
                                                       'wind_direction_raw', 'latitude', 'longitude', "GPS_speed", "GPS_heading", "compass_heading"],
                                            file_type=file_type, path=paths["local_data"])
    
        boat[b].filter_GPScoverage()
        boat[b].masks_for_harbors()
        boat[b].correct_winds()
        boat[b].calculate_windvector_components(corrected=True)
        boat[b].calculate_wind_sector(corrected=True)
        boat[b].calculate_wind_in_knots(corrected=True)
        
        # add info to the class instance
        boat[b].boat_name = boat_names[b]
        
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
        
        # add info to the class instance
        lighthouse[l].station_name = lighthouses[l]["name"]
        lighthouse[l].latitude = lighthouses[l]["lat"]
        lighthouse[l].longitude = lighthouses[l]["lon"]
        
        
    ####################################################################################    
    
    # create plot(s)
    for b in boats_to_plot:
        
        fig, gs, ax_map, sc_map = initialize_halfpage_map()
        cbar_range = get_cbar_range([boat[b].data[map_vari]], min_cbar_range)
        plot_boat_on_map(ax_map, sc_map, boat[b], variable=map_vari, position_switch=False, legend_switch=False, cbar_switch=True, fixed_cbar_range=cbar_range)
        plot_boat_timeseries(boat[b], fig, gs, status)
        
        if for_website:
            local_output_path = f"{paths['local_desktop']}daily_overview_plot_{boat[b].boat_name}.png"
        
            plt.savefig(local_output_path)
            plt.close("all")
            
            upload_picture(local_output_path, os.path.basename(local_output_path))
            
        else:
            local_output_path = f"{paths['local_desktop']}{boat[b].boat_name}_{(update_time-datetime.timedelta(days=1)).strftime('%Y%m%d')}.png"
        
            plt.savefig(local_output_path)
            plt.close("all")
            
        
    plot_lighthouse_timeseries(lighthouse, status)
    
    if for_website:
        local_output_path = f"{paths['local_desktop']}daily_overview_plot_lighthouses.png"
    
        plt.savefig(local_output_path)
        plt.close("all")
        
        upload_picture(local_output_path, os.path.basename(local_output_path))
        
    else:
        local_output_path = f"{paths['local_desktop']}lighthouses_{(update_time-datetime.timedelta(days=1)).strftime('%Y%m%d')}.png"
    
        plt.savefig(local_output_path)
        plt.close("all")
    
    return






def upload_picture(local_file, online_file_name):
    
    with ftplib.FTP("unisacsi.com", "u834143596", "unis12345") as session:
        with open(local_file, "rb") as f:
            session.storbinary(f"STOR {online_file_name}", f)

    return






####################################################################################
####################################################################################







if __name__ == '__main__':

    next_wakeup_time = next_wakeup() # start first wake-up

    # always true, to keep the script running forever
    while True:                 
        while datetime.datetime.now() < next_wakeup_time:       # sleep until the next scheduled wakeup time
            time.sleep(1)
            
        proc=mp.Process(target=update_overview_plot, args=[next_wakeup_time])
        proc.daemon=True
        proc.start()
        proc.join()
        
        next_wakeup_time = next_wakeup()           # don't forget to update the next wakeup time!!!