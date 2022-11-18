# -*- coding: utf-8 -*-
"""
Created on Sat Mar 27 14:37:44 2021

@author: unismet
"""

import os
import time
import datetime
from plot_functions import plot_lighthouse_timeseries
import matplotlib.pyplot as plt
import lighthouse_AWS_class
import ftplib
import yaml

import multiprocessing as mp


def next_wakeup():
    
    # refresh period
    dt_minutes = 5
    dt_hours = 0
    time_delta=datetime.timedelta(hours=dt_hours, minutes=dt_minutes)
    round_to = time_delta.total_seconds()                   # 60s

    dt = datetime.datetime.now()
    seconds = (dt - dt.min).seconds                       

    if seconds % round_to == 0 and dt.microsecond == 0:
        rounding = (seconds + round_to / 2) // round_to * round_to
    else:
        rounding = (seconds + dt.microsecond/1000000 + round_to) // round_to * round_to

    next_wakeup_time = dt + datetime.timedelta(0, rounding - seconds, - dt.microsecond) + datetime.timedelta(minutes=2)
    print("The next wakeup is scheduled for: {a}".format(a=next_wakeup_time))

    return next_wakeup_time



def update_all_plots(update_time):
    
    # define constants
    with open("./path_config.yaml", "r") as f:
        paths = yaml.safe_load(f)


    lighthouses = {1884: {"name": "Narveneset", 'lat': 78.56343,'lon': 16.29687},
                   1885: {"name": "Bohemanneset", 'lat': 78.38166, 'lon': 14.75300},
                   1886: {"name": "Daudmannsodden", 'lat': 78.21056,'lon': 12.98685},
                   1887: {"name": "Gasoyane", 'lat': 78.45792,'lon': 16.20082}}
    
    lighthouses_to_plot = [1884, 1885, 1886, 1887]

    status = "live"
    
    ####################################################################################
    ####################################################################################
    
    # define what data time period to load
    latest_update_time = datetime.datetime.strftime(update_time, "%Y%m%d%H%M")
    period_to_load = datetime.timedelta(days=1)
    start_time = (update_time - period_to_load).strftime("%Y%m%d%H%M")
    
    ####################################################################################
    ####################################################################################
    
    # load lighthouse data
    lighthouse = {}
    for l in lighthouses_to_plot:
        lighthouse[l] = lighthouse_AWS_class.lighthouse_AWS(station=l, resolution="1min",
                                            starttime=start_time, endtime=latest_update_time,
                                            variables=['temperature', 'air_pressure', 'relative_humidity', 'wind_speed', 'wind_direction', "battery"],
                                            file_type="raw", path=paths["local_data"])
    
        lighthouse[l].calculate_windvector_components()
        lighthouse[l].calculate_wind_sector()
        lighthouse[l].calculate_wind_in_knots()
        lighthouse[l].only_latest_data(datetime.timedelta(hours=72))
        
        # add info to the class instance
        lighthouse[l].station_name = lighthouses[l]["name"]
        lighthouse[l].latitude = lighthouses[l]["lat"]
        lighthouse[l].longitude = lighthouses[l]["lon"]

    ####################################################################################
    ####################################################################################


    # create and upload plots

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