# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 18:19:30 2021

@author: unismet
"""

import datetime
from mobile_stations_schedule_daily_overviews import update_overview_plot



day = "20211011"



#################################################################################################################
#################################################################################################################


from_time = datetime.datetime.strptime(day, "%Y%m%d").replace(hour=0, minute=15, second=0)
update_time = from_time + datetime.timedelta(days=1)

update_overview_plot(update_time)