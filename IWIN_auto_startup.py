# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 13:17:19 2022

@author: unismet
"""

import datetime
import yaml

# read config files
with open("./config_paths.yaml", "r", encoding='utf-8') as f:
    paths = yaml.safe_load(f)
    
    
with open("./config_processing_plotting.yaml", "r", encoding='utf-8') as f:
    stations_to_process = yaml.safe_load(f)


for k in stations_to_process.keys():
    if stations_to_process[k] is None:
        stations_to_process[k] = []
stations_to_restructure = stations_to_process["mobile_stations"] + stations_to_process["lighthouse_stations"]


mobile_resolutions = ["10min", "1min", "20sec"]
lighthouse_resolutions = ["10min", "1min"]


latest_day_not_fully_restructured = {}
for station in stations_to_restructure:
    if station in stations_to_process["mobile_stations"]:
        resolutions = mobile_resolutions
    elif station in stations_to_process["lighthouse_stations"]:
        resolutions = lighthouse_resolutions
          
    latest_day_not_fully_restructured[station] = {r:datetime.date.today()-datetime.timedelta(days=1) for r in resolutions}
    


with open(f"{paths['local_desktop']}iwin_latest_day_not_fully_restructured.yaml", 'w') as f:
    yaml.dump(latest_day_not_fully_restructured, f, default_flow_style=False)