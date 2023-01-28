# -*- coding: utf-8 -*-
"""
Created on Wed Dec 21 12:14:17 2022

@author: unismet
"""

import datetime
import yaml
import os
import argparse
import shutil









###########################################################################################################################
###########################################################################################################################
###########################################################################################################################



if __name__ == '__main__':
    
    p = argparse.ArgumentParser()
    p.add_argument("number_of_days")
    args = p.parse_args()

    number_of_days = int(args.number_of_days)
    
    mobile_stations = [1883, 1872, 1924]

    lighthouse_stations = [1884, 1885, 1886, 1887]
    
    resolutions = ["20sec", "1min", "10min"]
    
    
    with open("./config_paths.yaml", "r", encoding='utf-8') as f:
        paths = yaml.safe_load(f)
        
    days = sorted([datetime.date.today()-datetime.timedelta(days=i) for i in range(number_of_days)])

    # first delete and re-create the transfer folder
    shutil.rmtree(paths["MET_transfer"], ignore_errors=True)
    os.makedirs(paths["MET_transfer"])

    
    # copy all available files from the last x days
    for d in days:
        for resolution in resolutions:
            for station in mobile_stations:
                if os.path.exists(f"{paths['local_storage']}mobile_AWS_{station}/{resolution}/{d.year}/{d.month:02d}/mobile_AWS_{station}_Table_{resolution}_{d.year}{d.month:02d}{d.day:02d}.nc"):
                    if not os.path.exists(f"{paths['MET_transfer']}mobile_AWS_{station}/"):
                        os.makedirs(f"{paths['MET_transfer']}mobile_AWS_{station}/")
                    if not os.path.exists(f"{paths['MET_transfer']}mobile_AWS_{station}/{resolution}/"):
                        os.makedirs(f"{paths['MET_transfer']}mobile_AWS_{station}/{resolution}/")
                    if not os.path.exists(f"{paths['MET_transfer']}mobile_AWS_{station}/{resolution}/{d.year}/"):
                        os.makedirs(f"{paths['MET_transfer']}mobile_AWS_{station}/{resolution}/{d.year}/")
                    if not os.path.exists(f"{paths['MET_transfer']}mobile_AWS_{station}/{resolution}/{d.year}/{d.month:02d}/"):
                        os.makedirs(f"{paths['MET_transfer']}mobile_AWS_{station}/{resolution}/{d.year}/{d.month:02d}/")
                    shutil.copyfile(f"{paths['local_storage']}mobile_AWS_{station}/{resolution}/{d.year}/{d.month:02d}/mobile_AWS_{station}_Table_{resolution}_{d.year}{d.month:02d}{d.day:02d}.nc",
                                    f"{paths['MET_transfer']}mobile_AWS_{station}/{resolution}/{d.year}/{d.month:02d}/mobile_AWS_{station}_Table_{resolution}_{d.year}{d.month:02d}{d.day:02d}.nc")
                    
            for station in lighthouse_stations:
                if os.path.exists(f"{paths['local_storage']}lighthouse_AWS_{station}/{resolution}/{d.year}/{d.month:02d}/lighthouse_AWS_{station}_Table_{resolution}_{d.year}{d.month:02d}{d.day:02d}.nc"):
                    if not os.path.exists(f"{paths['MET_transfer']}lighthouse_AWS_{station}/"):
                        os.makedirs(f"{paths['MET_transfer']}lighthouse_AWS_{station}/")
                    if not os.path.exists(f"{paths['MET_transfer']}lighthouse_AWS_{station}/{resolution}/"):
                        os.makedirs(f"{paths['MET_transfer']}lighthouse_AWS_{station}/{resolution}/")
                    if not os.path.exists(f"{paths['MET_transfer']}lighthouse_AWS_{station}/{resolution}/{d.year}/"):
                        os.makedirs(f"{paths['MET_transfer']}lighthouse_AWS_{station}/{resolution}/{d.year}/")
                    if not os.path.exists(f"{paths['MET_transfer']}lighthouse_AWS_{station}/{resolution}/{d.year}/{d.month:02d}/"):
                        os.makedirs(f"{paths['MET_transfer']}lighthouse_AWS_{station}/{resolution}/{d.year}/{d.month:02d}/")
                    shutil.copyfile(f"{paths['local_storage']}lighthouse_AWS_{station}/{resolution}/{d.year}/{d.month:02d}/lighthouse_AWS_{station}_Table_{resolution}_{d.year}{d.month:02d}{d.day:02d}.nc",
                                    f"{paths['MET_transfer']}lighthouse_AWS_{station}/{resolution}/{d.year}/{d.month:02d}/lighthouse_AWS_{station}_Table_{resolution}_{d.year}{d.month:02d}{d.day:02d}.nc")
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        