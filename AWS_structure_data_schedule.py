# -*- coding: utf-8 -*-
"""
Created on Sat Mar 27 14:37:44 2021
@author: unismet
"""


import time
import datetime
from AWS_structure_data_functions import restructure_mobile_AWS, restructure_lighthouse_AWS
import os
import shutil
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




def restructure_AWS(now, stat, res, latest_day, new_day):
    
    mobile_stations = [1883, 1872, 1924]
    lighthouse_stations = [1884, 1885, 1886, 1887]
    

    midnight = datetime.datetime.combine(now, datetime.datetime.min.time())
    if stat in mobile_stations:
        latest_avail_time = restructure_mobile_AWS(midnight, now, station=str(stat), resolution=res)
    elif stat in lighthouse_stations:
        latest_avail_time = restructure_lighthouse_AWS(midnight, now, station=str(stat), resolution=res)

    if ((latest_day != now.date()) & (latest_avail_time > midnight)):
        print("final restructure yesterday")
        if stat in mobile_stations:
            restructure_mobile_AWS(midnight-datetime.timedelta(days=1), midnight, station=str(stat), resolution=res)
        elif stat in lighthouse_stations:
            restructure_lighthouse_AWS(midnight-datetime.timedelta(days=1), midnight, station=str(stat), resolution=res)
     
        new_day.put(now.date())
        
    else:
        new_day.put(latest_day)


####################################################################################
####################################################################################



if __name__ == '__main__':
    
    # define for which stations the program should run
    mobile_switches = {1883: False, 
                       1872: False,
                       1924: False}

    lighthouse_switches = {1885: True,
                           1884: True,
                           1886: True,
                           1887: True}
    
    stations_to_restructure = [s for s, sw in mobile_switches.items() if sw] + [s for s, sw in lighthouse_switches.items() if sw]
    
    mobile_resolutions = ["10min", "1min", "20sec"]
    lighthouse_resolutions = ["10min", "1min"]


    latest_day_not_fully_restructured = {}
    for station in stations_to_restructure:
        if station in mobile_switches.keys():
            resolutions = mobile_resolutions
        elif station in lighthouse_switches.keys():
            resolutions = lighthouse_resolutions
              
        latest_day_not_fully_restructured[station] = {r:datetime.date.today()-datetime.timedelta(days=1) for r in resolutions}
        
    # start first wake-up
    next_wakeup_time = next_wakeup()   
    
    # always true, to keep the script running forever
    while True:                 
        while datetime.datetime.now() < next_wakeup_time:       # sleep until the next scheduled wakeup time
            time.sleep(1)
                    
        for station in stations_to_restructure:
            
            if station in mobile_switches.keys():
                resolutions = mobile_resolutions
            elif station in lighthouse_switches.keys():
                resolutions = lighthouse_resolutions
                
            for res in resolutions:
                new_latest_day = mp.Queue()
                proc=mp.Process(target=restructure_AWS, args=[next_wakeup_time, station, res, latest_day_not_fully_restructured[station][res], new_latest_day])
                proc.start()
                latest_day_not_fully_restructured[station][res] = new_latest_day.get()
                proc.join()
        
        next_wakeup_time = next_wakeup()           # don't forget to update the next wakeup time!!!
        