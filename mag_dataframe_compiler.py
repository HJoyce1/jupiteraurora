#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 20 09:54:32 2023

@author: hannah

function to compile all seperate mag data files into one dataframe
"""

import pandas as pd
import spiceypy as sp
import glob

# first we want to glob all the magnetic field data files and sort them into order
files = sorted(glob.glob('/Users/hannah/OneDrive - Lancaster University/aurora/pds/1minavg/*.csv'))

# rearrange time columns so can translate to UTC format 
big_file = pd.read_csv(files[0],delimiter=',')
for i in files:
    csv = pd.read_csv(i,delimiter=',')
    big_file = pd.concat([big_file, csv])

# fill in extra 0s to make strings into time format
big_file['hour'] = big_file['hour'].astype(str).str.pad(2, side='left', fillchar='0')
big_file['min'] = big_file['min'].astype(str).str.pad(2, side='left', fillchar='0')
big_file = big_file.rename(columns={'min':'minute'})

# make column for date-time
big_file_UTC = big_file.assign(UTC = big_file.year.astype(str)+'-' + big_file.doy.astype(str)+'T'+big_file.hour.astype(str)+':'+big_file.minute.astype(str)+':'+big_file.sec.astype(str)+'.'+big_file.msec.astype(str))

# can save out file
big_file_UTC.to_csv('/Users/hannah/OneDrive - Lancaster University/aurora/python_scripts/dataframes/mag_df.csv',index=False)
