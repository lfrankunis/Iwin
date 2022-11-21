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
import pandas as pd


# define for which stations the program should run
mobile_switches = {1883: False, 
                   1872: False,
                   1924: False}

lighthouse_switches = {1884: False,
                       1885: True,
                       1886: False,
                       1887: False}




#################################################################################################################
#################################################################################################################

if len(sys.argv) == 3:
    days = pd.date_range(sys.argv[1], sys.argv[2], freq="D")
elif len(sys.argv) == 2:
    days = [pd.Timestamp(sys.argv[1])]

days_to_process = [i.strftime("%Y%m%d") for i in days]



# define path to the data folder
with open("./config_paths.yaml", "r", encoding='utf-8') as f:
    paths = yaml.safe_load(f)
    



for day in days_to_process:
    print(day)
    
    from_time = datetime.datetime.strptime(day, "%Y%m%d").replace(hour=0, minute=0, second=0)
    to_time = from_time + datetime.timedelta(days=1)
    
    print(f"Manual processing of {from_time}")
    
    if from_time < datetime.datetime(2022,5,17):
        mobile_resolutions = ["hour", "10min", "5min", "1min", "20sec"]
    else:
        mobile_resolutions = ["10min", "1min", "20sec"]
        
    if ((from_time < datetime.datetime(2021,10,18)) | ((from_time >= datetime.datetime(2022,3,7)) & ((from_time < datetime.datetime(2022,5,8))))):
        lighthouse_resolutions = ["hour", "10min", "1min"]
    elif ((from_time >= datetime.datetime(2021,10,18)) & ((from_time < datetime.datetime(2022,3,7)))):
        lighthouse_resolutions = ["1min"]
    else:
        lighthouse_resolutions = ["10min", "1min"]

        
    for station, switch in mobile_switches.items():
        if ((station == 1883) & (day in ["20210509", "20210510", "20210511", "20210622"])):
            continue
        if switch:
            # call the function to restructure
            for res in mobile_resolutions:
                try:
                    restructure_mobile_AWS(from_time, to_time, station=str(station), resolution=res, path_in=paths['local_data'], path_out=paths["local_storage"])
                    
                    # shutil.copyfile(f"{paths['local_storage']}mobile_AWS_{station}/{res}/mobile_AWS_{station}_Table_{res}_{from_time.year}{from_time.month:02d}{from_time.day:02d}.nc",
                    #                 f"{paths['harddrive']}mobile_AWS_{station}/{res}/mobile_AWS_{station}_Table_{res}_{from_time.year}{from_time.month:02d}{from_time.day:02d}.nc")
                    
                except FileNotFoundError:
                    pass



    for station, switch in lighthouse_switches.items():
        if switch:
            # call the function to restructure
            for res in lighthouse_resolutions:
                try:
                    restructure_lighthouse_AWS(from_time, to_time, station=str(station), resolution=res, path_in=paths['local_data'], path_out=paths["local_storage"])
                    
                    # shutil.copyfile(f"{paths['local_storage']}lighthouse_AWS_{station}/{res}/lighthouse_AWS_{station}_Table_{res}_{from_time.year}{from_time.month:02d}{from_time.day:02d}.nc",
                    #                 f"{paths['harddrive']}lighthouse_AWS_{station}/{res}/lighthouse_AWS_{station}_Table_{res}_{from_time.year}{from_time.month:02d}{from_time.day:02d}.nc")
                    
                except FileNotFoundError:
                    pass
    
