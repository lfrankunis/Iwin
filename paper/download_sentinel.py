#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 08:19:00 2022

@author: lukasf
"""

from collections import OrderedDict
from sentinelsat import SentinelAPI

api = SentinelAPI('lfrank', 'Dicksonland2020')

tiles = {'33XVG': ('2022-04-11T12:00:00.0Z', '2022-04-11T13:00:00.0Z'),
	 '33XVH': ('2022-04-11T12:00:00.0Z', '2022-04-11T13:00:00.0Z'),
	 '33XWG': ('2022-04-05T12:00:00.0Z', '2022-04-05T13:00:00.0Z'),
	 '33XWH': ('2022-04-05T12:00:00.0Z', '2022-04-05T13:00:00.0Z')}

query_kwargs = {
        'platformname': 'Sentinel-2',
	'producttype': 'S2MSI1C',
        'cloudcoverpercentage': (0, 10)}

products = OrderedDict()
for tile, date in tiles.items():
    kw = query_kwargs.copy()
    kw['tileid'] = tile
    kw["date"] = date
    pp = api.query(**kw)
    products.update(pp)

api.to_geodataframe(products)

#api.download_all(products, directory_path="/Users/lukasf/Desktop/Iwin_paper_figures/Sentinel/")