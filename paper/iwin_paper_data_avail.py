# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 13:26:46 2023

@author: Florina
"""
# -*- coding: utf-8 -*-


import matplotlib.pyplot as plt
import netCDF4
import numpy as np
import datetime as dt
from matplotlib.dates import HourLocator, DateFormatter
import matplotlib as mpl
import pandas as pd
import glob
import os
from datetime import datetime, timedelta

import xarray as xr

aws_pos = ['mobile_AWS_MSBard',
           'mobile_AWS_MSBerg',
           'mobile_AWS_MSBillefjord', 
           'mobile_AWS_MSPolargirl', 
           'lighthouse_AWS_Bohemanneset',
           'lighthouse_AWS_Daudmannsodden',
           'lighthouse_AWS_Gasoyane',
           'lighthouse_AWS_Narveneset']

for i,aws in enumerate(aws_pos):
    print(aws)
    ds = xr.open_dataset("https://thredds.met.no/thredds/dodsC/met.no/observations/unis/{aws}_1min".format(aws = aws))
    date = ds.time[:].dt.date
    locals()[aws] = np.unique(date)
   
locals()["mobile_AWS_MSBerg"] = np.concatenate([locals()["mobile_AWS_MSBard"], locals()["mobile_AWS_MSBerg"]])
aws_pos = aws_pos[1:]


startdate = np.min(locals()[aws_pos[0]])
enddate =  np.max(locals()[aws_pos[0]])

for j, aws in enumerate(aws_pos[1:-1]):
    temp_start = np.min(locals()[aws])
    temp_end = np.max(locals()[aws])
    
    if startdate > temp_start:
        startdate = temp_start
        
    if enddate < temp_end:
        enddate = temp_end
        
day_vector = list(np.arange(startdate, enddate, dt.timedelta(days=1), dtype=dt.datetime))
data_avail = pd.DataFrame({'day':day_vector})

for k, aws in enumerate(aws_pos):
    bools = np.zeros_like(day_vector,dtype = bool)
    for l, day in enumerate(day_vector):
        ind = np.where(day_vector[l] == locals()[aws])
        if len(ind[0])>0:
            bools[l] = True
    data_avail[aws] = bools
    

#%% plot 

mpl.rcdefaults()

params = {'legend.fontsize': 'x-large',
                    'figure.figsize': (18, 9),
                    'font.family':'serif',
                    'axes.labelsize': 10,
                    'axes.titlesize':16,
                    'xtick.labelsize':16,
                    'ytick.labelsize':10,
                    'lines.linewidth' : 3,
                    'lines.markersize' : 10,
                    'hatch.linewidth' : 0.1, }
plt.rcParams.update(params)


  
from matplotlib import colors
short = ['MS BD\nMS BG', 'MS BF', 'MS PG','BHN','DMO', 'GO','NN']
cmap1 = colors.ListedColormap(['#FFFFFF','#03009E'])#(0.33725490196078434, 0.7058823529411765, 0.9137254901960784)])
cmap2 = colors.ListedColormap(['#FFFFFF','#34AA00'])#(0.33725490196078434, 0.7058823529411765, 0.9137254901960784)])

col = data_avail.columns.values.tolist() 

gs_top = plt.GridSpec(22, 3, top=0.97, hspace = 0.2, bottom = 0.1, right = 0.985, left = 0.085)

fig = plt.figure(figsize = ([12,5]))

fig.align_ylabels() 

for i, aws in enumerate(aws_pos):
    if i <3:
        j = i*3
        k = (i*3)+3
    else:
        j = (i*3)+1
        k = (i*3)+4
    ax = fig.add_subplot(gs_top[j:k,:])
    
    Z = np.random.rand(1, int(len(data_avail['day'])))
    Z[0] = data_avail[col[i+1]][0:int(len(data_avail['day']))]
    X = np.array(data_avail['day'])
    X = np.append(X, X[-1]+timedelta(days = 1))
    Y = np.arange(Z.shape[0]+1)
    if i <3:
        ax.pcolormesh(X,Y,Z, cmap = cmap1)
    else: 
        ax.pcolormesh(X,Y,Z, cmap = cmap2)
    ax.set_ylabel(short[i],ha='right', va="center", y=0.5, fontsize = 16, rotation=0, labelpad=2)#loc="bottom", rotation="horizontal", fontsize=16, labelpad=40)
    ax.set(yticklabels=[])
    ax.set_yticks([])
    
    if i<6:
        ax.set(xticklabels=[])
        ax.set_xticks([])
    
ax.xaxis.set_major_formatter(DateFormatter('%b-%Y'))
ax.xaxis.set_major_locator(plt.MaxNLocator(10))


plt.show()

plt.savefig("/Users/lukasf/Desktop/Iwin_paper_figures/iwin_paper_fig_08.pdf", dpi=300)


