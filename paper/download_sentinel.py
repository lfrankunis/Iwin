#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 08:19:00 2022

@author: lukasf
"""

from collections import OrderedDict
from sentinelsat import SentinelAPI

api = SentinelAPI('lukasf', 'Dicksonland2020', 'https://apihub.copernicus.eu/apihub/')

tiles = ['33XVG', '33XWG', '33XWH', '33XVH']

query_kwargs = {
        'platformname': 'Sentinel-2',
        'producttype': 'S2MSI2A',
        'date': ('20220401', '20220414'),
        'cloudcoverpercentage': (0, 10)}

products = OrderedDict()
for tile in tiles:
    kw = query_kwargs.copy()
    kw['tileid'] = tile
    pp = api.query(**kw)
    products.update(pp)

api.download_all(products, directory_path="/Users/lukasf/Desktop/Iwin_paper_figures/Sentinel/")