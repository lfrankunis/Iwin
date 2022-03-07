# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 18:07:58 2021

@author: unismet
"""

import datetime
from AWS_structure_data_functions import restructure_mobile_AWS, restructure_lighthouse_AWS
import os
import sys
import shutil



# define for which stations the program should run
AWS_1883 = False
AWS_1872 = False
AWS_1924 = False
Bohemanneset = True




#################################################################################################################
#################################################################################################################


days_to_process = [str(i) for i in sys.argv[1:]]


# define path to the data folder
path = "C:/Data/"
path_agf350 = "C:/Data/AGF350/"



for day in days_to_process:
    
    from_time = datetime.datetime.strptime(day, "%Y%m%d").replace(hour=0, minute=0, second=0)
    to_time = from_time + datetime.timedelta(days=1)
    
    print(f"Manual processing of {from_time}")
    
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
        # try:
        #     restructure_lighthouse_AWS(from_time, to_time, station="1885", resolution="hour", path=path)
        # except FileNotFoundError:
        #     pass
        # try:
        #     restructure_lighthouse_AWS(from_time, to_time, station="1885", resolution="10min", path=path)
        # except FileNotFoundError:
        #     pass
        try:
            restructure_lighthouse_AWS(from_time, to_time, station="1885", resolution="1min", path=path)
        except FileNotFoundError:
            pass
    
        shutil.copytree(path + "lighthouse_AWS_1885/{a}{b:02d}{c:02d}".format(a=from_time.year, b=from_time.month, c=from_time.day),
                        "D:/DATA/lighthouse_AWS_1885/{a}{b:02d}{c:02d}".format(a=from_time.year, b=from_time.month, c=from_time.day))
