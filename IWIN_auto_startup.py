# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 13:17:19 2022

@author: unismet
"""

import datetime
import yaml

# define for which stations the program should run
mobile_switches = {1883: True, 
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
    

with open("./config_paths.yaml", "r", encoding='utf-8') as f:
    paths = yaml.safe_load(f)

with open(f"{paths['local_desktop']}iwin_latest_day_not_fully_restructured.yaml", 'w') as f:
    yaml.dump(latest_day_not_fully_restructured, f, default_flow_style=False)