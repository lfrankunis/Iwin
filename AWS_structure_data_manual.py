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
import yaml


# define for which stations the program should run
mobile_switches = {1883: False, 
                   1872: False,
                   1924: False}

lighthouse_switches = {1885: False,
                       1884: False,
                       1886: False,
                       1887: False}




#################################################################################################################
#################################################################################################################


days_to_process = [str(i) for i in sys.argv[1:]]


# define path to the data folder
with open("./path_config.yaml", "r") as f:
    paths = yaml.safe_load(f)
    



for day in days_to_process:
    
    from_time = datetime.datetime.strptime(day, "%Y%m%d").replace(hour=0, minute=0, second=0)
    to_time = from_time + datetime.timedelta(days=1)
    
    print(f"Manual processing of {from_time}")
        
        
    for station, switch in mobile_switches.items():
        if switch:
        # create directories
            os.mkdir(f"{paths['local_data']}mobile_AWS_{station}/{from_time.year}{from_time.month:02d}{from_time.day:02d}")
            os.mkdir(f"{paths['local_data']}mobile_AWS_{station}/{from_time.year}{from_time.month:02d}{from_time.day:02d}/ascii")
            os.mkdir(f"{paths['local_data']}mobile_AWS_{station}/{from_time.year}{from_time.month:02d}{from_time.day:02d}/nc")
    
            # call the function to restructure
            for res in ["10min", "1min", "20sec"]:          # for data prior to 20220507, add "5min", "hour" 
                try:
                    restructure_mobile_AWS(from_time, to_time, station=str(station), resolution=res, path=paths['local_data'])
                except FileNotFoundError:
                    pass

            shutil.copytree(f"{paths['local_data']}mobile_AWS_{station}/{from_time.year}{from_time.month:02d}{from_time.day:02d}",
                            f"{paths['harddrive']}mobile_AWS_{station}/{from_time.year}{from_time.month:02d}{from_time.day:02d}")
    

    for station, switch in lighthouse_switches.items():
        if switch:
            # create directories
            os.mkdir(f"{paths['local_data']}lighthouse_AWS_{station}/{from_time.year}{from_time.month:02d}{from_time.day:02d}")
            os.mkdir(f"{paths['local_data']}lighthouse_AWS_{station}/{from_time.year}{from_time.month:02d}{from_time.day:02d}/ascii")
            os.mkdir(f"{paths['local_data']}lighthouse_AWS_{station}/{from_time.year}{from_time.month:02d}{from_time.day:02d}/nc")
    
            # call the function to restructure
            for res in ["10min", "1min"]:                   # for data prior to 20220507, add "hour" 
                try:
                    restructure_lighthouse_AWS(from_time, to_time, station=str(station), resolution=res, path=paths['local_data'])
                except FileNotFoundError:
                    pass
    
            shutil.copytree(f"{paths['local_data']}lighthouse_AWS_{station}/{from_time.year}{from_time.month:02d}{from_time.day:02d}",
                            f"{paths['harddrive']}lighthouse_AWS_{station}/{from_time.year}{from_time.month:02d}{from_time.day:02d}")

