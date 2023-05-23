# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 13:27:07 2022

@author: unismet
"""

import time
import yaml
import multiprocessing as mp
import datetime
from IWIN_structure_data_functions import restructure_mobile_AWS, restructure_lighthouse_AWS




def restructure_AWS(stat, res, latest_day, new_day):
    
    mobile_stations = ["MSBard", "MSBerg", "MSBillefjord", "MSPolargirl", "RVHannaResvoll"]
    lighthouse_stations = ["Narveneset", "Bohemanneset", "Daudmannsodden", "Gasoyane", "KappThordsen"]
    
    now = datetime.datetime.now()
    now = now.replace(tzinfo=datetime.timezone.utc)
    
    midnight = datetime.datetime.combine(now, datetime.datetime.min.time())
    midnight = midnight.replace(tzinfo=datetime.timezone.utc)
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





###########################################################################################################################
###########################################################################################################################
###########################################################################################################################



if __name__ == '__main__':

    with open("./config_paths.yaml", "r", encoding='utf-8') as f:
        paths = yaml.safe_load(f)
    
    with open(f"{paths['local_desktop']}iwin_latest_day_not_fully_restructured.yaml", "r", encoding='utf-8') as f:
        latest_day_not_fully_restructured = yaml.safe_load(f)
        
    
    for station, resolutions in latest_day_not_fully_restructured.items():
            
        for res in resolutions.keys():
            new_latest_day = mp.Queue()
            proc=mp.Process(target=restructure_AWS, args=[station, res, latest_day_not_fully_restructured[station][res], new_latest_day])
            proc.start()
            latest_day_not_fully_restructured[station][res] = new_latest_day.get()
            proc.join()
            
    with open(f"{paths['local_desktop']}iwin_latest_day_not_fully_restructured.yaml", 'w') as f:
        yaml.dump(latest_day_not_fully_restructured, f, default_flow_style=False)
        
        
        
        
        
        
        
        