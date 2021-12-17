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
AWS_1 = True
AWS_2 = True
AWS_3 = True
Bohemanneset = True




#################################################################################################################
#################################################################################################################


days_to_process = [str(i) for i in sys.argv[1:]]


# define path to the data folder
path = "C:/Data/"
path_rooftest = "C:/Data/UNIS_roof_test_winter/"



for day in days_to_process:
    
    from_time = datetime.datetime.strptime(day, "%Y%m%d").replace(hour=0, minute=0, second=0)
    to_time = from_time + datetime.timedelta(days=1)
    
    print(f"Manual processing of {from_time}")
    
    if AWS_1:
        # create directories
        os.mkdir(path_rooftest + "mobile_AWS_1/{a}{b:02d}{c:02d}".format(a=from_time.year, b=from_time.month, c=from_time.day))
        os.mkdir(path_rooftest + "mobile_AWS_1/{a}{b:02d}{c:02d}/ascii".format(a=from_time.year, b=from_time.month, c=from_time.day))
        os.mkdir(path_rooftest + "mobile_AWS_1/{a}{b:02d}{c:02d}/nc".format(a=from_time.year, b=from_time.month, c=from_time.day))
    
        # call the function to restructure
        try:
            restructure_mobile_AWS(from_time, to_time, station="1", resolution="hour", path=path_rooftest)
        except FileNotFoundError:
            pass
        try:
            restructure_mobile_AWS(from_time, to_time, station="1", resolution="10min", path=path_rooftest)
        except FileNotFoundError:
            pass
        try:
            restructure_mobile_AWS(from_time, to_time, station="1", resolution="5min", path=path_rooftest)
        except FileNotFoundError:
            pass
        try:
            restructure_mobile_AWS(from_time, to_time, station="1", resolution="1min", path=path_rooftest)
        except FileNotFoundError:
            pass
        try:
            restructure_mobile_AWS(from_time, to_time, station="1", resolution="20sec", path=path_rooftest)
        except FileNotFoundError:
            pass
    
        shutil.copytree(path_rooftest + "mobile_AWS_1/{a}{b:02d}{c:02d}".format(a=from_time.year, b=from_time.month, c=from_time.day),
                        "D:/DATA/UNIS_roof_test_winter/mobile_AWS_1/{a}{b:02d}{c:02d}".format(a=from_time.year, b=from_time.month, c=from_time.day))
    
    
    if AWS_2:
        # create directories
        os.mkdir(path_rooftest + "mobile_AWS_2/{a}{b:02d}{c:02d}".format(a=from_time.year, b=from_time.month, c=from_time.day))
        os.mkdir(path_rooftest + "mobile_AWS_2/{a}{b:02d}{c:02d}/ascii".format(a=from_time.year, b=from_time.month, c=from_time.day))
        os.mkdir(path_rooftest + "mobile_AWS_2/{a}{b:02d}{c:02d}/nc".format(a=from_time.year, b=from_time.month, c=from_time.day))
    
        # call the function to restructure
        try:
            restructure_mobile_AWS(from_time, to_time, station="2", resolution="hour", path=path_rooftest)
        except FileNotFoundError:
            pass
        try:
            restructure_mobile_AWS(from_time, to_time, station="2", resolution="10min", path=path_rooftest)
        except FileNotFoundError:
            pass
        try:
            restructure_mobile_AWS(from_time, to_time, station="2", resolution="5min", path=path_rooftest)
        except FileNotFoundError:
            pass
        try:
            restructure_mobile_AWS(from_time, to_time, station="2", resolution="1min", path=path_rooftest)
        except FileNotFoundError:
            pass
        try:
            restructure_mobile_AWS(from_time, to_time, station="2", resolution="20sec", path=path_rooftest)
        except FileNotFoundError:
            pass
    
        shutil.copytree(path_rooftest + "mobile_AWS_2/{a}{b:02d}{c:02d}".format(a=from_time.year, b=from_time.month, c=from_time.day),
                        "D:/DATA/UNIS_roof_test_winter/mobile_AWS_2/{a}{b:02d}{c:02d}".format(a=from_time.year, b=from_time.month, c=from_time.day))
    if AWS_3:
        # create directories
        os.mkdir(path_rooftest + "mobile_AWS_3/{a}{b:02d}{c:02d}".format(a=from_time.year, b=from_time.month, c=from_time.day))
        os.mkdir(path_rooftest + "mobile_AWS_3/{a}{b:02d}{c:02d}/ascii".format(a=from_time.year, b=from_time.month, c=from_time.day))
        os.mkdir(path_rooftest + "mobile_AWS_3/{a}{b:02d}{c:02d}/nc".format(a=from_time.year, b=from_time.month, c=from_time.day))
     
        # call the function to restructure
        try:
            restructure_mobile_AWS(from_time, to_time, station="3", resolution="hour", path=path_rooftest)
        except FileNotFoundError:
            pass
        try:
            restructure_mobile_AWS(from_time, to_time, station="3", resolution="10min", path=path_rooftest)
        except FileNotFoundError:
            pass
        try:
            restructure_mobile_AWS(from_time, to_time, station="3", resolution="5min", path=path_rooftest)
        except FileNotFoundError:
            pass
        try:
            restructure_mobile_AWS(from_time, to_time, station="3", resolution="1min", path=path_rooftest)
        except FileNotFoundError:
            pass
        try:
            restructure_mobile_AWS(from_time, to_time, station="3", resolution="20sec", path=path_rooftest)
        except FileNotFoundError:
            pass
     
        shutil.copytree(path_rooftest + "mobile_AWS_3/{a}{b:02d}{c:02d}".format(a=from_time.year, b=from_time.month, c=from_time.day),
                        "D:/DATA/UNIS_roof_test_winter/mobile_AWS_3/{a}{b:02d}{c:02d}".format(a=from_time.year, b=from_time.month, c=from_time.day))
 
    
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
