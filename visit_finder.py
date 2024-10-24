# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 13:03:08 2023

@author: hannah

script that pulls the time data directly from the fits files for the HST visits
this data is then saved out to a dataframe so it can be used to identify appropriate
solar wind data for each visit
"""

import numpy as np
import pandas as pd
from dateutil.parser import parse
from tqdm import tqdm
import glob
from astropy.io import fits
import datetime as dt
import scipy.constants as c
import pandas as pd

# mac
root_folder = '/Users/hannah/OneDrive - Lancaster University/aurora/python_scripts/dataframes/'


visit_list = ['01','02','03','04','05','08','09','10','11','12','13','15','16','17','18','19','20','21','23','24','25','26','27','28','34','35'] #,'09','10','11','12','13','15','16','18','19','20','21'
#travel_time_df = pd.read_csv('/Users/hannah/OneDrive - Lancaster University/aurora/python_scripts/datatravel_times_sw_df.csv',delimiter=',')

all_times = []
#visit_05_time = []
# we grab all the files we are interested getting times from
for l in visit_list:      
        year= 2016
        time = str(100)
        extra = 'nichols/'
        arch = '*_v'+ l
        ti = str('/*0'+time+'*')
        ab = glob.glob(f'/Users/hannah/OneDrive - Lancaster University/aurora/data/{year}/extract/{extra}'+arch+'/nopolar'+time+ti)
        ab.sort()
        
        num_image = []
        visit_time=[]
        for  n,i in tqdm(enumerate(ab)):
             hdulist = fits.open(i)
             header = hdulist[1].header
             # jupiter times
             exp_time = header['EXPT']
             start_time = parse(header['UDATE'])     # create datetime object
             try:
                 dist_org = header['DIST_ORG']
                 ltime = dist_org*c.au/c.c
                 lighttime = dt.timedelta(seconds=ltime)
             except KeyError:
                    lighttime = dt.timedelta(seconds=2524.42) 
                    exposure = dt.timedelta(seconds=exp_time)
        
             start_time_jup = start_time - lighttime     
             visit_time.append(start_time_jup)
             num_image.append(n+1)
             # visit_times_df = visit_times_df.assign(Times=start_time_jup)
             # visit_times_df = visit_times_df.assign(Visit_Number=l)
        
        all_times.append(visit_time)

        all_times_array = np.array(all_times)
        all_times_array = all_times_array.T
        
        
df = pd.DataFrame(data=all_times_array, index = num_image, columns=visit_list)
df = df.add_prefix('Visit_')

df.to_csv(root_folder+'visit_times.csv',index=False)
