# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 07:41:03 2022

@author: unismet
"""

import datetime
import requests
import pandas as pd
import numpy as np


def download_MET_stations(end, variables_to_read=["temperature", "relative_humidity", "wind_speed", "wind_direction", "pressure"]):
    """
    Requests the data from the API, converts the timestamps and returns the data in a pandas dataframe
    """
    
    station_names = {"IR": {"name": "Isfjord Radio", "ID": "SN99790", "height": 7., "lat": 78.0625, "lon": 13.6192},
                     "LYR": {"name": "Longyearbyen Airport", "ID": "SN99840", "height": 28., "lat": 78.2453, "lon": 15.5015},
                     "PYR": {"name": "Pyramiden", "ID": "SN99880", "height": 20., "lat": 78.6557, "lon": 16.3603},
                     "NS": {"name": "Nedre Sassendalen", "ID": "SN99882", "height": 13., "lat": 78.3313, "lon": 16.6818},
                     # "PB": {"name": "Plataberget", "ID": "SN99843", "height": 450., "lat": 78.2278, "lon": 15.378},
                     # "AD": {"name": "Adventdalen", "ID": "SN99870", "height": 15., "lat":  78.2022, "lon": 15.831},
                     # "IH": {"name": "Innerhytta", "ID": "SN99879", "height": 81., "lat":  78.18883, "lon": 16.34423},
                     # "AO": {"name": "Akseloya", "ID": "SN99765", "height": 20., "lat":  77.6873, "lon": 14.7578},
                     # "KL": {"name": "Klauva", "ID": "SN99884", "height": 480., "lat":  78.3002, "lon": 18.2248},
                     # "ID": {"name": "Istjorndalen", "ID": "SN99770", "height": 188., "lat":  78.0092, "lon": 15.2108}, 
                     # "JH": {"name": "Janssonhagen", "ID": "SN99874", "height": 250., "lat":  78.18, "lon": 16.41},
                     # "RP": {"name": "Reindalspasset", "ID": "SN99763", "height": 181., "lat":  78.0648, "lon":  17.0442},
                     # "SV": {"name": "Svea", "ID": "SN99760", "height": 9., "lat":  77.8953, "lon": 16.72}
                     }
    
    stations = ",".join([station_names[s]["ID"] for s in list(station_names.keys())])

    variables_dict = {"temperature": "air_temperature", "relative_humidity": "relative_humidity",
                      "SLP": "air_pressure_at_sea_level", "pressure": "surface_air_pressure",
                      "wind_speed": "wind_speed", "wind_direction": "wind_from_direction"}
    
    variables = ",".join([variables_dict[v] for v in variables_to_read])
    

    client_ID = 'bd5334a6-4be1-4a60-a64b-67bb111bc5cc'
    # client_secret = '21bb21e4-cc59-4dcf-8480-51bd0c448896'


    start = end - datetime.timedelta(hours=12)
    start_str = start.strftime('%Y-%m-%dT%H:%M')
    end_str = end.strftime('%Y-%m-%dT%H:%M')

 

    # Define endpoint and parameters
    endpoint = 'https://frost.met.no/observations/v0.jsonld'
    parameters = {
        'sources': stations,
        'elements': variables,
        'referencetime': start_str + '/' + end_str
    }

    # Issue an HTTP GET request
    r = requests.get(endpoint, parameters, auth=(client_ID, ''))
    # Extract JSON data
    json = r.json()
   
    # Check if the request worked, print out any errors
    if r.status_code == 200:
        data = json['data']
        print('Data retrieved from frost.met.no!')

    else:
        print('Error! Returned status code %s' % r.status_code)
        print('Message: %s' % json['error']['message'])
        print('Reason: %s' % json['error']['reason'])
        assert False, "Download failed!"

 

    # This will return a Dataframe with all of the observations in a table format
    df = pd.DataFrame()
    for i in range(len(data)):
        row = pd.DataFrame(data[i]['observations'])
        row['referenceTime'] = data[i]['referenceTime']
        row['sourceId'] = data[i]['sourceId']
        df = pd.concat([df, row])

    # Convert the time value to something Python understands
    df['TIMESTAMP'] = pd.to_datetime(df['referenceTime'])
    df = df.drop(columns='referenceTime')
    df.set_index("TIMESTAMP", inplace=True)
    df.index = df.index.tz_convert('Europe/Berlin')
    
    
    data_dict = {}
    for s in list(station_names.keys()):
        data_station = []
        for vari in variables_to_read:
            d = df[((df["sourceId"] == f"{station_names[s]['ID']}:0") & (df["elementId"] == variables_dict[vari]))]["value"]
            d.sort_index(inplace=True)
            d.rename(vari, inplace=True)
            d = d[d.index.minute == 0]
            data_station.append(d)
        data_dict[s] = pd.concat(data_station, axis=1).drop_duplicates().dropna(thresh=len(variables_to_read)//2)
        data_dict[s]["u"] = -np.abs(data_dict[s]['wind_speed']) * np.sin(np.deg2rad(data_dict[s]['wind_direction']))
        data_dict[s]["v"] = -np.abs(data_dict[s]['wind_speed']) * np.cos(np.deg2rad(data_dict[s]['wind_direction']))
        data_dict[s]['wind_sector'] = (np.array((data_dict[s]['wind_direction']/45.)+.5, dtype=int) % 8) + 1
        data_dict[s]['u_knts'] = 1.94384 * data_dict[s]['u']
        data_dict[s]['v_knts'] = 1.94384 * data_dict[s]['v']

    return data_dict
