# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 18:19:30 2021

@author: unismet
"""

import datetime
import sys
from mobile_stations_schedule_daily_overviews import update_overview_plot



day = str(sys.argv[1])


#################################################################################################################
#################################################################################################################


from_time = datetime.datetime.strptime(day, "%Y%m%d").replace(hour=0, minute=15, second=0)
update_time = from_time + datetime.timedelta(days=1)

# for reprocessing daily overview plot for website:
# update_overview_plot(update_time, file_type="raw", for_website=True)


# for quicklook what happened during one specific day: (comment again after use, and uncomment regular line)
update_overview_plot(update_time, file_type="nc", for_website=False)
