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
AWS_1883 = False
AWS_1872 = False
AWS_1924 = False
Bohemanneset = True

# define path to the data folder
path = "C:/Data/"
path_agf350 = "C:/Data/AGF350/"

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

    if AWS_1883:
        # create directories
        os.mkdir(path_agf350 + "mobile_AWS_1883/{a}{b:02d}{c:02d}".format(a=from_time.year, b=from_time.month, c=from_time.day))
        os.mkdir(path_agf350 + "mobile_AWS_1883/{a}{b:02d}{c:02d}/ascii".format(a=from_time.year, b=from_time.month, c=from_time.day))
        os.mkdir(path_agf350 + "mobile_AWS_1883/{a}{b:02d}{c:02d}/nc".format(a=from_time.year, b=from_time.month, c=from_time.day))

        # call the function to restructure
        try:
            restructure_mobile_AWS(from_time, to_time, station="1883", resolution="hour", path=path_agf350)
        except FileNotFoundError:
            pass
        try:
            restructure_mobile_AWS(from_time, to_time, station="1883", resolution="10min", path=path_agf350)
        except FileNotFoundError:
            pass
        try:
            restructure_mobile_AWS(from_time, to_time, station="1883", resolution="5min", path=path_agf350)
        except FileNotFoundError:
            pass
        try:
            restructure_mobile_AWS(from_time, to_time, station="1883", resolution="1min", path=path_agf350)
        except FileNotFoundError:
            pass
        try:
            restructure_mobile_AWS(from_time, to_time, station="1883", resolution="20sec", path=path_agf350)
        except FileNotFoundError:
            pass

        shutil.copytree(path_agf350 + "mobile_AWS_1883/{a}{b:02d}{c:02d}".format(a=from_time.year, b=from_time.month, c=from_time.day),
                        "D:/DATA/AGF350/mobile_AWS_1883/{a}{b:02d}{c:02d}".format(a=from_time.year, b=from_time.month, c=from_time.day))


    if AWS_1872:
        # create directories
        os.mkdir(path_agf350 + "mobile_AWS_1872/{a}{b:02d}{c:02d}".format(a=from_time.year, b=from_time.month, c=from_time.day))
        os.mkdir(path_agf350 + "mobile_AWS_1872/{a}{b:02d}{c:02d}/ascii".format(a=from_time.year, b=from_time.month, c=from_time.day))
        os.mkdir(path_agf350 + "mobile_AWS_1872/{a}{b:02d}{c:02d}/nc".format(a=from_time.year, b=from_time.month, c=from_time.day))

        # call the function to restructure
        try:
            restructure_mobile_AWS(from_time, to_time, station="1872", resolution="hour", path=path_agf350)
        except FileNotFoundError:
            pass
        try:
            restructure_mobile_AWS(from_time, to_time, station="1872", resolution="10min", path=path_agf350)
        except FileNotFoundError:
            pass
        try:
            restructure_mobile_AWS(from_time, to_time, station="1872", resolution="5min", path=path_agf350)
        except FileNotFoundError:
            pass
        try:
            restructure_mobile_AWS(from_time, to_time, station="1872", resolution="1min", path=path_agf350)
        except FileNotFoundError:
            pass
        try:
            restructure_mobile_AWS(from_time, to_time, station="1872", resolution="20sec", path=path_agf350)
        except FileNotFoundError:
            pass

        shutil.copytree(path_agf350 + "mobile_AWS_1872/{a}{b:02d}{c:02d}".format(a=from_time.year, b=from_time.month, c=from_time.day),
                        "D:/DATA/AGF350/mobile_AWS_1872/{a}{b:02d}{c:02d}".format(a=from_time.year, b=from_time.month, c=from_time.day))
    if AWS_1924:
        # create directories
        os.mkdir(path_agf350 + "mobile_AWS_1924/{a}{b:02d}{c:02d}".format(a=from_time.year, b=from_time.month, c=from_time.day))
        os.mkdir(path_agf350 + "mobile_AWS_1924/{a}{b:02d}{c:02d}/ascii".format(a=from_time.year, b=from_time.month, c=from_time.day))
        os.mkdir(path_agf350 + "mobile_AWS_1924/{a}{b:02d}{c:02d}/nc".format(a=from_time.year, b=from_time.month, c=from_time.day))

        # call the function to restructure
        try:
            restructure_mobile_AWS(from_time, to_time, station="1924", resolution="hour", path=path_agf350)
        except FileNotFoundError:
            pass
        try:
            restructure_mobile_AWS(from_time, to_time, station="1924", resolution="10min", path=path_agf350)
        except FileNotFoundError:
            pass
        try:
            restructure_mobile_AWS(from_time, to_time, station="1924", resolution="5min", path=path_agf350)
        except FileNotFoundError:
            pass
        try:
            restructure_mobile_AWS(from_time, to_time, station="1924", resolution="1min", path=path_agf350)
        except FileNotFoundError:
            pass
        try:
            restructure_mobile_AWS(from_time, to_time, station="1924", resolution="20sec", path=path_agf350)
        except FileNotFoundError:
            pass

        shutil.copytree(path_agf350 + "mobile_AWS_1924/{a}{b:02d}{c:02d}".format(a=from_time.year, b=from_time.month, c=from_time.day),
                        "D:/DATA/AGF350/mobile_AWS_1924/{a}{b:02d}{c:02d}".format(a=from_time.year, b=from_time.month, c=from_time.day))


    if Bohemanneset:
        # create directories
        os.mkdir(path + "lighthouse_AWS_1885/{a}{b:02d}{c:02d}".format(a=from_time.year, b=from_time.month, c=from_time.day))
        os.mkdir(path + "lighthouse_AWS_1885/{a}{b:02d}{c:02d}/ascii".format(a=from_time.year, b=from_time.month, c=from_time.day))
        os.mkdir(path + "lighthouse_AWS_1885/{a}{b:02d}{c:02d}/nc".format(a=from_time.year, b=from_time.month, c=from_time.day))

        # call the function to restructure
        try:
            restructure_lighthouse_AWS(from_time, to_time, station="1885", resolution="hour", path=path)
        except FileNotFoundError:
            pass
        try:
            restructure_lighthouse_AWS(from_time, to_time, station="1885", resolution="10min", path=path)
        except FileNotFoundError:
            pass
        try:
            restructure_lighthouse_AWS(from_time, to_time, station="1885", resolution="1min", path=path)
        except FileNotFoundError:
            pass

        shutil.copytree(path + "lighthouse_AWS_1885/{a}{b:02d}{c:02d}".format(a=from_time.year, b=from_time.month, c=from_time.day),
                        "D:/DATA/lighthouse_AWS_1885/{a}{b:02d}{c:02d}".format(a=from_time.year, b=from_time.month, c=from_time.day))


    next_wakeup()           # don't forget to update the next wakeup time!!!