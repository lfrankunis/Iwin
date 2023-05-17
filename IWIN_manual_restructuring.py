# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 18:07:58 2021

@author: unismet
"""

import datetime
from IWIN_structure_data_functions import restructure_mobile_AWS, restructure_lighthouse_AWS
import shutil
import yaml
import pandas as pd
import argparse



p = argparse.ArgumentParser()
p.add_argument('-s', "--stations", dest='stations', nargs='*')
p.add_argument('-t', "--time", dest='time', nargs='*')

args = p.parse_args()


#################################################################################################################
#################################################################################################################

if len(args.time) == 2:
    days = pd.date_range(args.time[0], args.time[1], freq="D")
elif len(args.time) == 1:
    days = [pd.Timestamp(args.time[0])]

days_to_process = [i.strftime("%Y%m%d") for i in days]


mobile_stations = [s for s in args.stations if s in ["MSBard", "MSBillefjord", "MSPolargirl", "RVHannaResvoll"]]
lighthouse_stations = [s for s in args.stations if s in ["Narveneset", "Bohemanneset", "Daudmannsodden", "Gasoyane", "KappThordsen"]]



# define path to the data folder
with open("./config_paths.yaml", "r", encoding='utf-8') as f:
    paths = yaml.safe_load(f)
    

for day in days_to_process:
    
    from_time = datetime.datetime.strptime(day, "%Y%m%d").replace(hour=0, minute=0, second=0, tzinfo=datetime.timezone.utc)
    to_time = from_time + datetime.timedelta(days=1)
    
    print(f"Manual processing of {from_time.date()}")
    
    if from_time < datetime.datetime(2022,5,17, tzinfo=datetime.timezone.utc):
        mobile_resolutions = ["hour", "10min", "5min", "1min", "20sec"]
    else:
        mobile_resolutions = ["10min", "1min", "20sec"]
        
    if ((from_time < datetime.datetime(2021,10,18, tzinfo=datetime.timezone.utc)) | ((from_time >= datetime.datetime(2022,3,7, tzinfo=datetime.timezone.utc)) & ((from_time < datetime.datetime(2022,5,8, tzinfo=datetime.timezone.utc))))):
        lighthouse_resolutions = ["hour", "10min", "1min"]
    elif ((from_time >= datetime.datetime(2021,10,18, tzinfo=datetime.timezone.utc)) & ((from_time < datetime.datetime(2022,3,7, tzinfo=datetime.timezone.utc)))):
        lighthouse_resolutions = ["1min"]
    else:
        lighthouse_resolutions = ["10min", "1min"]

        
    for station in mobile_stations:
        if ((station == "MSBard") & (day in ["20210509", "20210510", "20210511", "20210622"])):
            continue
        for res in mobile_resolutions:
            restructure_mobile_AWS(from_time, to_time, station=station, resolution=res)

    for station in lighthouse_stations:
        for res in lighthouse_resolutions:
            restructure_lighthouse_AWS(from_time, to_time, station=station, resolution=res)

