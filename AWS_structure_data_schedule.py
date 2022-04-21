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

def next_wakeup():
    global next_wakeup_time
    global round_to
    global dt_minutes_offset


    dt = datetime.datetime.now()
    seconds = (dt - dt.min).seconds

    if seconds % round_to == 0 and dt.microsecond == 0:
        rounding = (seconds + round_to / 2) // round_to * round_to
    else:
        rounding = (seconds + dt.microsecond/1000000 + round_to) // round_to * round_to

    next_wakeup_time = dt + datetime.timedelta(0, rounding - seconds, - dt.microsecond)
    next_wakeup_time += datetime.timedelta(minutes=dt_minutes_offset)
    print("The next wakeup is scheduled for: {a}".format(a=next_wakeup_time))

    return

# define for which stations the program should run
mobile_switches = {1883: False, 
                   1872: False,
                   1924: True}

lighthouse_switches = {1885: True,
                       1884: False,
                       1886: False,
                       1887: False}

# define path to the data folder
path = "C:/Data/"

# define resolution of output files (daily files/hourly files/minute files)
dt_days = 1
dt_hours = 0
dt_minutes = 0
dt_minutes_offset = 15      # to shift the python precessing until the Loggernet data downloading is completed

time_delta=datetime.timedelta(days=dt_days, hours=dt_hours, minutes=dt_minutes)
round_to = time_delta.total_seconds()

next_wakeup()


while True:                 # always true, to keep the script running forever
    while datetime.datetime.now() < next_wakeup_time:       # sleep until the next scheduled wakeup time
        time.sleep(1)

    from_time = next_wakeup_time - time_delta
    to_time = next_wakeup_time

    # take away the offset again
    from_time -= datetime.timedelta(minutes=dt_minutes_offset)
    to_time -= datetime.timedelta(minutes=dt_minutes_offset)

    for station, switch in mobile_switches.items():
        if switch:
        # create directories
            os.mkdir(f"{path}mobile_AWS_{station}/{from_time.year}{from_time.month:02d}{from_time.day:02d}")
            os.mkdir(f"{path}mobile_AWS_{station}/{from_time.year}{from_time.month:02d}{from_time.day:02d}/ascii")
            os.mkdir(f"{path}mobile_AWS_{station}/{from_time.year}{from_time.month:02d}{from_time.day:02d}/nc")
    
            # call the function to restructure
            for res in ["hour", "10min", "5min", "1min", "20sec"]:
                try:
                    restructure_mobile_AWS(from_time, to_time, station=str(station), resolution=res, path=path)
                except FileNotFoundError:
                    pass

            shutil.copytree(f"{path}mobile_AWS_{station}/{from_time.year}{from_time.month:02d}{from_time.day:02d}",
                            f"D:/DATA/mobile_AWS_{station}/{from_time.year}{from_time.month:02d}{from_time.day:02d}")
    



    for station, switch in lighthouse_switches.items():
        if switch:
            # create directories
            os.mkdir(f"{path}lighthouse_AWS_{station}/{from_time.year}{from_time.month:02d}{from_time.day:02d}")
            os.mkdir(f"{path}lighthouse_AWS_{station}/{from_time.year}{from_time.month:02d}{from_time.day:02d}/ascii")
            os.mkdir(f"{path}lighthouse_AWS_{station}/{from_time.year}{from_time.month:02d}{from_time.day:02d}/nc")

        # call the function to restructure
        for res in ["hour", "10min", "1min"]:
            try:
                restructure_lighthouse_AWS(from_time, to_time, station=str(station), resolution=res, path=path)
            except FileNotFoundError:
                pass


        shutil.copytree(f"{path}lighthouse_AWS_{station}/{from_time.year}{from_time.month:02d}{from_time.day:02d}",
                        f"D:/DATA/lighthouse_AWS_{station}/{from_time.year}{from_time.month:02d}{from_time.day:02d}")


    next_wakeup()           # don't forget to update the next wakeup time!!!