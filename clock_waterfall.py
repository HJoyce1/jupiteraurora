#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 25 11:50:47 2026

@author: hannah

waterfall plotting for clock angle data
"""

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
import pandas as pd
from dateutil.parser import parse
import datetime as dt
import matplotlib.patheffects as pe
from matplotlib.path import Path
import spiceypy as spice
import scipy.constants as c
import glob
import datetime
import matplotlib.dates as mdates
import os
from datetime import timedelta
import pandas as pd
import matplotlib.ticker as ticker


root = '/Users/hannah/OneDrive - Lancaster University/aurora/python_scripts/dataframes/'

# data for each visit
visit_02_data = pd.read_csv(root+'visit_02_data_aug.csv', delimiter=',')
visit_03_data = pd.read_csv(root+'visit_03_data_aug.csv', delimiter=',')
visit_04_data = pd.read_csv(root+'visit_04_data_aug.csv', delimiter=',')
visit_05_data = pd.read_csv(root+'visit_05_data_aug.csv', delimiter=',')
visit_08_data = pd.read_csv(root+'visit_08_data_aug.csv', delimiter=',')
visit_09_data = pd.read_csv(root+'visit_09_data_aug.csv', delimiter=',')
visit_10_data = pd.read_csv(root+'visit_10_data_aug.csv', delimiter=',')
visit_11_data = pd.read_csv(root+'visit_11_data_aug.csv', delimiter=',')
visit_16_data = pd.read_csv(root+'visit_16_data_aug.csv', delimiter=',')
visit_17_data = pd.read_csv(root+'visit_17_data_aug.csv', delimiter=',')
visit_18_data = pd.read_csv(root+'visit_18_data_aug.csv', delimiter=',')
visit_19_data = pd.read_csv(root+'visit_19_data_aug.csv', delimiter=',')
visit_20_data = pd.read_csv(root+'visit_20_data_aug.csv', delimiter=',')
visit_21_data = pd.read_csv(root+'visit_21_data_aug.csv', delimiter=',')
visit_24_data = pd.read_csv(root+'visit_24_data_aug.csv', delimiter=',')
visit_25_data = pd.read_csv(root+'visit_25_data_aug.csv', delimiter=',')
visit_27_data = pd.read_csv(root+'visit_27_data_aug.csv', delimiter=',')
visit_28_data = pd.read_csv(root+'visit_28_data_aug.csv', delimiter=',')
visit_34_data = pd.read_csv(root+'visit_34_data_aug.csv', delimiter=',')
visit_35_data = pd.read_csv(root+'visit_35_data_aug.csv', delimiter=',')


# # high cml - 01, 12, 15, 26
visit_01_data = pd.read_csv(root+'visit_01_data_aug.csv', delimiter=',')
visit_12_data = pd.read_csv(root+'visit_12_data_aug.csv', delimiter=',')
visit_15_data = pd.read_csv(root+'visit_15_data_aug.csv', delimiter=',')
visit_26_data = pd.read_csv(root+'visit_26_data_aug.csv', delimiter=',')

# standard times
times_02 = visit_02_data['Ionospheric_Time'].to_numpy()
times_03 = visit_03_data['Ionospheric_Time'].to_numpy()
times_04 = visit_04_data['Ionospheric_Time'].to_numpy()
times_05 = visit_05_data['Ionospheric_Time'].to_numpy()
times_08 = visit_08_data['Ionospheric_Time'].to_numpy()
times_09 = visit_09_data['Ionospheric_Time'].to_numpy()
times_10 = visit_10_data['Ionospheric_Time'].to_numpy()
times_11 = visit_11_data['Ionospheric_Time'].to_numpy()
times_16 = visit_16_data['Ionospheric_Time'].to_numpy()
times_17 = visit_17_data['Ionospheric_Time'].to_numpy()
times_18 = visit_18_data['Ionospheric_Time'].to_numpy()
times_19 = visit_19_data['Ionospheric_Time'].to_numpy()
times_20 = visit_20_data['Ionospheric_Time'].to_numpy()
times_21 = visit_21_data['Ionospheric_Time'].to_numpy()
times_24 = visit_24_data['Ionospheric_Time'].to_numpy()
times_25 = visit_25_data['Ionospheric_Time'].to_numpy()
times_27 = visit_27_data['Ionospheric_Time'].to_numpy()
times_28 = visit_28_data['Ionospheric_Time'].to_numpy()
times_34 = visit_34_data['Ionospheric_Time'].to_numpy()
times_35 = visit_35_data['Ionospheric_Time'].to_numpy()

# high times
times_01 = visit_01_data['Ionospheric_Time'].to_numpy()
times_12 = visit_12_data['Ionospheric_Time'].to_numpy()
times_15 = visit_15_data['Ionospheric_Time'].to_numpy()
times_26 = visit_26_data['Ionospheric_Time'].to_numpy()
'''
# Group datasets into lists
standard_visits = [visit_02_data, visit_03_data, visit_04_data, visit_05_data, 
                   visit_08_data, visit_09_data, visit_10_data, visit_11_data, 
                   visit_16_data, visit_17_data, visit_18_data, visit_19_data, 
                   visit_20_data, visit_21_data, visit_24_data, visit_25_data, 
                   visit_27_data, visit_28_data, visit_34_data, visit_35_data]

high_res_visits = [visit_01_data, visit_12_data, visit_15_data, visit_26_data]

all_visits = standard_visits + high_res_visits

for i, raw_array in enumerate(all_visits):
    # 1. Convert the NumPy array to datetime
    dt_times = pd.to_datetime(raw_array)
    
    # 2. Calculate minutes from the first point of THIS array
    # This forces every visit to start at 0 on the X-axis
    minutes = (dt_times - dt_times[0]) / pd.Timedelta(minutes=1)
    
    # 3. Get your corresponding data values
    # y = all_data_arrays[i] + (i * spacing)
'''

# Standard Visits
m_02 = (pd.to_datetime(times_02) - pd.to_datetime(times_02)[0]) / pd.Timedelta(minutes=1)
m_03 = (pd.to_datetime(times_03) - pd.to_datetime(times_03)[0]) / pd.Timedelta(minutes=1)
m_04 = (pd.to_datetime(times_04) - pd.to_datetime(times_04)[0]) / pd.Timedelta(minutes=1)
m_05 = (pd.to_datetime(times_05) - pd.to_datetime(times_05)[0]) / pd.Timedelta(minutes=1)
m_08 = (pd.to_datetime(times_08) - pd.to_datetime(times_08)[0]) / pd.Timedelta(minutes=1)
m_09 = (pd.to_datetime(times_09) - pd.to_datetime(times_09)[0]) / pd.Timedelta(minutes=1)
m_10 = (pd.to_datetime(times_10) - pd.to_datetime(times_10)[0]) / pd.Timedelta(minutes=1)
m_11 = (pd.to_datetime(times_11) - pd.to_datetime(times_11)[0]) / pd.Timedelta(minutes=1)
m_16 = (pd.to_datetime(times_16) - pd.to_datetime(times_16)[0]) / pd.Timedelta(minutes=1)
m_17 = (pd.to_datetime(times_17) - pd.to_datetime(times_17)[0]) / pd.Timedelta(minutes=1)
m_18 = (pd.to_datetime(times_18) - pd.to_datetime(times_18)[0]) / pd.Timedelta(minutes=1)
m_19 = (pd.to_datetime(times_19) - pd.to_datetime(times_19)[0]) / pd.Timedelta(minutes=1)
m_20 = (pd.to_datetime(times_20) - pd.to_datetime(times_20)[0]) / pd.Timedelta(minutes=1)
m_21 = (pd.to_datetime(times_21) - pd.to_datetime(times_21)[0]) / pd.Timedelta(minutes=1)
m_24 = (pd.to_datetime(times_24) - pd.to_datetime(times_24)[0]) / pd.Timedelta(minutes=1)
m_25 = (pd.to_datetime(times_25) - pd.to_datetime(times_25)[0]) / pd.Timedelta(minutes=1)
m_27 = (pd.to_datetime(times_27) - pd.to_datetime(times_27)[0]) / pd.Timedelta(minutes=1)
m_28 = (pd.to_datetime(times_28) - pd.to_datetime(times_28)[0]) / pd.Timedelta(minutes=1)
m_34 = (pd.to_datetime(times_34) - pd.to_datetime(times_34)[0]) / pd.Timedelta(minutes=1)
m_35 = (pd.to_datetime(times_35) - pd.to_datetime(times_35)[0]) / pd.Timedelta(minutes=1)

# High Res Visits
m_01 = (pd.to_datetime(times_01) - pd.to_datetime(times_01)[0]) / pd.Timedelta(minutes=1)
m_12 = (pd.to_datetime(times_12) - pd.to_datetime(times_12)[0]) / pd.Timedelta(minutes=1)
m_15 = (pd.to_datetime(times_15) - pd.to_datetime(times_15)[0]) / pd.Timedelta(minutes=1)
m_26 = (pd.to_datetime(times_26) - pd.to_datetime(times_26)[0]) / pd.Timedelta(minutes=1)

# extract clock angle for each visit
clock_02 = visit_02_data['Clock_Angle'].to_numpy()
clock_03 = visit_03_data['Clock_Angle'].to_numpy()
clock_04 = visit_04_data['Clock_Angle'].to_numpy()
clock_05 = visit_05_data['Clock_Angle'].to_numpy()
clock_08 = visit_08_data['Clock_Angle'].to_numpy()
clock_09 = visit_09_data['Clock_Angle'].to_numpy()
clock_10 = visit_10_data['Clock_Angle'].to_numpy()
clock_11 = visit_11_data['Clock_Angle'].to_numpy()
clock_16 = visit_16_data['Clock_Angle'].to_numpy()
clock_17 = visit_17_data['Clock_Angle'].to_numpy()
clock_18 = visit_18_data['Clock_Angle'].to_numpy()
clock_19 = visit_19_data['Clock_Angle'].to_numpy()
clock_20 = visit_20_data['Clock_Angle'].to_numpy()
clock_21 = visit_21_data['Clock_Angle'].to_numpy()
clock_24 = visit_24_data['Clock_Angle'].to_numpy()
clock_25 = visit_25_data['Clock_Angle'].to_numpy()
clock_27 = visit_27_data['Clock_Angle'].to_numpy()
clock_28 = visit_28_data['Clock_Angle'].to_numpy()
clock_34 = visit_34_data['Clock_Angle'].to_numpy()
clock_35 = visit_35_data['Clock_Angle'].to_numpy()

# clock_13 = visit_13_data['Clock_Angle'].to_numpy()

# high cml
clock_01 = visit_01_data['Clock_Angle'].to_numpy()
clock_12 = visit_12_data['Clock_Angle'].to_numpy()
clock_15 = visit_15_data['Clock_Angle'].to_numpy()
clock_26 = visit_26_data['Clock_Angle'].to_numpy()

fig = plt.figure(figsize=(13, 25))
plt.subplots_adjust(hspace=0, wspace=0.15)
ax1 = plt.subplot(24,1,1)#hspace=0.20)
ax1.scatter(m_01,clock_01,s=8)
ax1.plot(m_01,clock_01)
ax1.set_ylabel('Clock \nAngle \n($^o$)',fontsize=14)
ax1.set_ylim(-190, 190)
ax1.tick_params(axis='x',which='minor',direction='in')
ax1.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax1.minorticks_on()
ax1.set_yticks([-90, 0, 90])
ax1.tick_params(which='major',direction='in',bottom=False, top=True, left=True, right=True, labelsize=14, length=5)
ax1.tick_params(which='minor', direction='in', bottom=False, top=True, left=True, right=True, length=2)
ax1.text(1.02, 0.5, 'Visit 01', transform=ax1.transAxes, 
         va='center', ha='left', fontsize=14)

ax2 = plt.subplot(24,1,2)
ax2.scatter(m_02,clock_02,s=8)
ax2.plot(m_02,clock_02)
ax2.set_ylabel('Clock \nAngle \n($^o$)',fontsize=14)
ax2.set_ylim(-190, 190)
ax2.tick_params(axis='x',which='minor',direction='in')
ax2.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax2.minorticks_on()
ax2.set_yticks([-90, 0, 90])
ax2.tick_params(which='major',direction='in',bottom=False, top=False, left=True, right=True, labelsize=14, length=5)
ax2.tick_params(which='minor', direction='in', bottom=False, top=False, left=True, right=True, length=2)
ax2.text(1.02, 0.5, 'Visit 02', transform=ax2.transAxes, 
         va='center', ha='left', fontsize=14)

ax3 = plt.subplot(24,1,3)#hspace=0.20)
ax3.scatter(m_03,clock_03,s=8)
ax3.plot(m_03,clock_03)
ax3.set_ylabel('Clock \nAngle \n($^o$)',fontsize=14)
ax3.set_ylim(-190, 190)
ax3.tick_params(axis='x',which='minor',direction='in')
ax3.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax3.minorticks_on()
ax3.set_yticks([-90, 0, 90])
ax3.tick_params(which='major',direction='in',bottom=False, top=False, left=True, right=True, labelsize=14, length=5)
ax3.tick_params(which='minor', direction='in', bottom=False, top=False, left=True, right=True, length=2)
ax3.text(1.02, 0.5, 'Visit 03', transform=ax3.transAxes, 
         va='center', ha='left', fontsize=14)

ax4 = plt.subplot(24,1,4)#hspace=0.20)
ax4.scatter(m_04,clock_04,s=8)
ax4.plot(m_04,clock_04)
ax4.set_ylabel('Clock \nAngle \n($^o$)',fontsize=14)
ax4.set_ylim(-190, 190)
ax4.tick_params(axis='x',which='minor',direction='in')
ax4.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax4.minorticks_on()
ax4.set_yticks([-90, 0, 90])
ax4.tick_params(which='major',direction='in',bottom=False, top=False, left=True, right=True, labelsize=14, length=5)
ax4.tick_params(which='minor', direction='in', bottom=False, top=False, left=True, right=True, length=2)
ax4.text(1.02, 0.5, 'Visit 04', transform=ax4.transAxes, 
         va='center', ha='left', fontsize=14)

ax5 = plt.subplot(24,1,5)#hspace=0.20)
ax5.scatter(m_05,clock_05, s=8)
ax5.plot(m_05,clock_05)
ax5.set_ylabel('Clock\n Angle \n($^o$)',fontsize=14)
ax5.set_ylim(-190, 190)
ax5.tick_params(axis='x',which='minor',direction='in')
ax5.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax5.minorticks_on()
ax5.set_yticks([-90, 0, 90])
ax5.tick_params(which='major',direction='in',bottom=False, top=False, left=True, right=True, labelsize=14, length=5)
ax5.tick_params(which='minor', direction='in', bottom=False, top=False, left=True, right=True, length=2)
ax5.text(1.02, 0.5, 'Visit 05', transform=ax5.transAxes, 
         va='center', ha='left', fontsize=14)

ax6 = plt.subplot(24,1,6)#hspace=0.20)
ax6.scatter(m_08,clock_08,s=8)
ax6.plot(m_08,clock_08)
ax6.set_ylabel('Clock \nAngle \n($^o$)',fontsize=14)
ax6.set_ylim(-190, 190)
ax6.tick_params(axis='x',which='minor',direction='in')
ax6.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax6.minorticks_on()
ax6.set_yticks([-90, 0, 90])
ax6.tick_params(which='major',direction='in',bottom=False, top=False, left=True, right=True, labelsize=14, length=5)
ax6.tick_params(which='minor', direction='in', bottom=False, top=False, left=True, right=True, length=2)
ax6.text(1.02, 0.5, 'Visit 08', transform=ax6.transAxes, 
         va='center', ha='left', fontsize=14)

ax7 = plt.subplot(24,1,7)#hspace=0.20)
ax7.scatter(m_09,clock_09,s=8)
ax7.plot(m_09,clock_09)
ax7.set_ylabel('Clock \nAngle \n($^o$)',fontsize=14)
ax7.set_ylim(-190, 190)
ax7.tick_params(axis='x',which='minor',direction='in')
ax7.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax7.minorticks_on()
ax7.set_yticks([-90, 0, 90])
ax7.tick_params(which='major',direction='in',bottom=False, top=False, left=True, right=True, labelsize=14, length=5)
ax7.tick_params(which='minor', direction='in', bottom=False, top=False, left=True, right=True, length=2)
ax7.text(1.02, 0.5, 'Visit 09', transform=ax7.transAxes, 
         va='center', ha='left', fontsize=14)

ax8 = plt.subplot(24,1,8)#hspace=0.20)
ax8.scatter(m_10,clock_10,s=8)
ax8.plot(m_10,clock_10)
ax8.set_ylabel('Clock \nAngle \n($^o$)',fontsize=14)
ax8.set_ylim(-190, 190)
ax8.tick_params(axis='x',which='minor',direction='in')
ax8.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax8.minorticks_on()
ax8.set_yticks([-90, 0, 90])
ax8.tick_params(which='major',direction='in',bottom=False, top=False, left=True, right=True, labelsize=14, length=5)
ax8.tick_params(which='minor', direction='in', bottom=False, top=False, left=True, right=True, length=2)
ax8.text(1.02, 0.5, 'Visit 10', transform=ax8.transAxes, 
         va='center', ha='left', fontsize=14)

ax9 = plt.subplot(24,1,9)#hspace=0.20)
ax9.scatter(m_11,clock_11,s=8)
ax9.plot(m_11,clock_11)
ax9.set_ylabel('Clock \nAngle \n($^o$)',fontsize=14)
ax9.set_ylim(-190, 190)
ax9.tick_params(axis='x',which='minor',direction='in')
ax9.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax9.minorticks_on()
ax9.set_yticks([-90, 0, 90])
ax9.tick_params(which='major',direction='in',bottom=False, top=False, left=True, right=True, labelsize=14, length=5)
ax9.tick_params(which='minor', direction='in', bottom=False, top=False, left=True, right=True, length=2)
ax9.text(1.02, 0.5, 'Visit 11', transform=ax9.transAxes, 
         va='center', ha='left', fontsize=14)

ax10 = plt.subplot(24,1,10)#hspace=0.20)
ax10.scatter(m_12,clock_12,s=8)
ax10.plot(m_12,clock_12)
ax10.set_ylabel('Clock \nAngle \n($^o$)',fontsize=14)
ax10.set_ylim(-190, 190)
ax10.tick_params(axis='x',which='minor',direction='in')
ax10.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax10.minorticks_on()
ax10.set_yticks([-90, 0, 90])
ax10.tick_params(which='major',direction='in',bottom=False, top=False, left=True, right=True, labelsize=14, length=5)
ax10.tick_params(which='minor', direction='in', bottom=False, top=False, left=True, right=True, length=2)
ax10.text(1.02, 0.5, 'Visit 12', transform=ax10.transAxes, 
         va='center', ha='left', fontsize=14)

ax11 = plt.subplot(24,1,11)#hspace=0.20)
ax11.scatter(m_15,clock_15,s=8)
ax11.plot(m_15,clock_15)
ax11.set_ylabel('Clock \nAngle \n($^o$)',fontsize=14)
ax11.set_ylim(-190, 190)
ax11.tick_params(axis='x',which='minor',direction='in')
ax11.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax11.minorticks_on()
ax11.set_yticks([-90, 0, 90])
ax11.tick_params(which='major',direction='in',bottom=False, top=False, left=True, right=True, labelsize=14, length=5)
ax11.tick_params(which='minor', direction='in', bottom=False, top=False, left=True, right=True, length=2)
ax11.text(1.02, 0.5, 'Visit 15', transform=ax11.transAxes, 
         va='center', ha='left', fontsize=14)

ax12 = plt.subplot(24,1,12)#hspace=0.20)
ax12.scatter(m_16,clock_16,s=8)
ax12.plot(m_16,clock_16)
ax12.set_ylabel('Clock \nAngle \n($^o$)',fontsize=14)
ax12.set_ylim(-190, 190)
ax12.tick_params(axis='x',which='minor',direction='in')
ax12.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax12.minorticks_on()
ax12.set_yticks([-90, 0, 90])
ax12.tick_params(which='major',direction='in',bottom=False, top=False, left=True, right=True, labelsize=14, length=5)
ax12.tick_params(which='minor', direction='in', bottom=False, top=False, left=True, right=True, length=2)
ax12.text(1.02, 0.5, 'Visit 16', transform=ax12.transAxes, 
         va='center', ha='left', fontsize=14)

ax13 = plt.subplot(24,1,13)#hspace=0.20)
ax13.scatter(m_17,clock_17,s=8)
ax13.plot(m_17,clock_17)
ax13.set_ylabel('Clock \nAngle \n($^o$)',fontsize=14)
ax13.set_ylim(-190, 190)
ax13.tick_params(axis='x',which='minor',direction='in')
ax13.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax13.minorticks_on()
ax13.set_yticks([-90, 0, 90])
ax13.tick_params(which='major',direction='in',bottom=False, top=False, left=True, right=True, labelsize=14, length=5)
ax13.tick_params(which='minor', direction='in', bottom=False, top=False, left=True, right=True, length=2)
ax13.text(1.02, 0.5, 'Visit 17', transform=ax13.transAxes, 
         va='center', ha='left', fontsize=14)

ax14 = plt.subplot(24,1,14)#hspace=0.20)
ax14.scatter(m_18,clock_18,s=8)
ax14.plot(m_18,clock_18)
ax14.set_ylabel('Clock \nAngle \n($^o$)',fontsize=14)
ax14.set_ylim(-190, 190)
ax14.tick_params(axis='x',which='minor',direction='in')
ax14.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax14.minorticks_on()
ax14.set_yticks([-90, 0, 90])
ax14.tick_params(which='major',direction='in',bottom=False, top=False, left=True, right=True, labelsize=14, length=5)
ax14.tick_params(which='minor', direction='in', bottom=False, top=False, left=True, right=True, length=2)
ax14.text(1.02, 0.5, 'Visit 18', transform=ax14.transAxes, 
         va='center', ha='left', fontsize=14)

ax15 = plt.subplot(24,1,15)#hspace=0.20)
ax15.scatter(m_19,clock_19,s=8)
ax15.plot(m_19,clock_19)
ax15.set_ylabel('Clock \nAngle \n($^o$)',fontsize=14)
ax15.set_ylim(-190, 190)
ax15.tick_params(axis='x',which='minor',direction='in')
ax15.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax15.minorticks_on()
ax15.set_yticks([-90, 0, 90])
ax15.tick_params(which='major',direction='in',bottom=False, top=False, left=True, right=True, labelsize=14, length=5)
ax15.tick_params(which='minor', direction='in', bottom=False, top=False, left=True, right=True, length=2)
ax15.text(1.02, 0.5, 'Visit 19', transform=ax15.transAxes, 
         va='center', ha='left', fontsize=14)

ax16 = plt.subplot(24,1,16)#hspace=0.20)
ax16.scatter(m_20,clock_20,s=8)
ax16.plot(m_20,clock_20)
ax16.set_ylabel('Clock \nAngle \n($^o$)',fontsize=14)
ax16.set_ylim(-190, 190)
ax16.tick_params(axis='x',which='minor',direction='in')
ax16.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax16.minorticks_on()
ax16.set_yticks([-90, 0, 90])
ax16.tick_params(which='major',direction='in',bottom=False, top=False, left=True, right=True, labelsize=14, length=5)
ax16.tick_params(which='minor', direction='in', bottom=False, top=False, left=True, right=True, length=2)
ax16.text(1.02, 0.5, 'Visit 20', transform=ax16.transAxes, 
         va='center', ha='left', fontsize=14)

ax17 = plt.subplot(24,1,17)#hspace=0.20)
ax17.scatter(m_21,clock_21,s=8)
ax17.plot(m_21,clock_21)
ax17.set_ylabel('Clock \nAngle \n($^o$)',fontsize=14)
ax17.set_ylim(-190, 190)
ax17.tick_params(axis='x',which='minor',direction='in')
ax17.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax17.minorticks_on()
ax17.set_yticks([-90, 0, 90])
ax17.tick_params(which='major',direction='in',bottom=False, top=False, left=True, right=True, labelsize=14, length=5)
ax17.tick_params(which='minor', direction='in', bottom=False, top=False, left=True, right=True, length=2)
ax17.text(1.02, 0.5, 'Visit 21', transform=ax17.transAxes, 
         va='center', ha='left', fontsize=14)

ax18 = plt.subplot(24,1,18)#hspace=0.20)
ax18.scatter(m_24,clock_24,s=8)
ax18.plot(m_24,clock_24)
ax18.set_ylabel('Clock \nAngle \n($^o$)',fontsize=12)
ax18.set_ylim(-190, 190)
ax18.tick_params(axis='x',which='minor',direction='in')
ax18.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax18.minorticks_on()
ax18.set_yticks([-90, 0, 90])
ax18.tick_params(which='major',direction='in',bottom=False, top=False, left=True, right=True, labelsize=14, length=5)
ax18.tick_params(which='minor', direction='in', bottom=False, top=False, left=True, right=True, length=2)
ax18.text(1.02, 0.5, 'Visit 24', transform=ax18.transAxes, 
         va='center', ha='left', fontsize=14)

ax19 = plt.subplot(24,1,19)#hspace=0.20)
ax19.scatter(m_25,clock_25,s=8)
ax19.plot(m_25,clock_25)
ax19.set_ylabel('Clock \nAngle \n($^o$)',fontsize=14)
ax19.set_ylim(-190, 190)
ax19.tick_params(axis='x',which='minor',direction='in')
ax19.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax19.minorticks_on()
ax19.set_yticks([-90, 0, 90])
ax19.tick_params(which='major',direction='in',bottom=False, top=False, left=True, right=True, labelsize=14, length=5)
ax19.tick_params(which='minor', direction='in', bottom=False, top=False, left=True, right=True, length=2)
ax19.text(1.02, 0.5, 'Visit 25', transform=ax19.transAxes, 
         va='center', ha='left', fontsize=14)

ax20 = plt.subplot(24,1,20)#hspace=0.20)
ax20.scatter(m_26,clock_26,s=8)
ax20.plot(m_26,clock_26)
ax20.set_ylabel('Clock \nAngle \n($^o$)',fontsize=14)
ax20.set_ylim(-190, 190)
ax20.tick_params(axis='x',which='minor',direction='in')
ax20.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax20.minorticks_on()
ax20.set_yticks([-90, 0, 90])
ax20.tick_params(which='major',direction='in',bottom=False, top=False, left=True, right=True, labelsize=14, length=5)
ax20.tick_params(which='minor', direction='in', bottom=False, top=False, left=True, right=True, length=2)
ax20.text(1.02, 0.5, 'Visit 26', transform=ax20.transAxes, 
         va='center', ha='left', fontsize=14)

ax21 = plt.subplot(24,1,21)#hspace=0.20)
ax21.scatter(m_27,clock_27,s=8)
ax21.plot(m_27,clock_27)
ax21.set_ylabel('Clock \nAngle \n($^o$)',fontsize=14)
ax21.set_ylim(-190, 190)
ax21.tick_params(axis='x',which='minor',direction='in')
ax21.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax21.minorticks_on()
ax21.set_yticks([-90, 0, 90])
ax21.tick_params(which='major',direction='in',bottom=False, top=False, left=True, right=True, labelsize=14, length=5)
ax21.tick_params(which='minor', direction='in', bottom=False, top=False, left=True, right=True, length=2)
ax21.text(1.02, 0.5, 'Visit 27', transform=ax21.transAxes, 
         va='center', ha='left', fontsize=14)

ax22 = plt.subplot(24,1,22)#hspace=0.20)
ax22.scatter(m_28,clock_28,s=8)
ax22.plot(m_28,clock_28)
ax22.set_ylabel('Clock \nAngle \n($^o$)',fontsize=14)
ax22.set_ylim(-190, 190)
ax22.tick_params(axis='x',which='minor',direction='in')
ax22.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax22.minorticks_on()
ax22.set_yticks([-90, 0, 90])
ax22.tick_params(which='major',direction='in',bottom=False, top=False, left=True, right=True, labelsize=14, length=5)
ax22.tick_params(which='minor', direction='in', bottom=False, top=False, left=True, right=True, length=2)
ax22.text(1.02, 0.5, 'Visit 28', transform=ax22.transAxes, 
         va='center', ha='left', fontsize=14)

ax23 = plt.subplot(24,1,23)#hspace=0.20)
ax23.scatter(m_34,clock_34,s=8)
ax23.plot(m_34,clock_34)
ax23.set_ylabel('Clock \nAngle \n($^o$)',fontsize=14)
ax23.set_ylim(-190, 190)
ax23.tick_params(axis='x',which='minor',direction='in')
ax23.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax23.minorticks_on()
ax23.set_yticks([-90, 0, 90])
ax23.tick_params(which='major',direction='in',bottom=False, top=False, left=True, right=True, labelsize=14, length=5)
ax23.tick_params(which='minor', direction='in', bottom=False, top=False, left=True, right=True, length=2)
ax23.text(1.02, 0.5, 'Visit 34', transform=ax23.transAxes, 
         va='center', ha='left', fontsize=14)

ax24 = plt.subplot(24,1,24)#hspace=0.20)
ax24.scatter(m_35,clock_35,s=8)
ax24.plot(m_35,clock_35)
ax24.set_ylabel('Clock \nAngle \n($^o$)',fontsize=14)
ax24.set_ylim(-190, 190)
ax24.tick_params(axis='x',which='minor',direction='in')
ax24.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax24.minorticks_on()
ax24.set_yticks([-90, 0, 90])
ax24.tick_params(which='major',direction='in',bottom=True, top=False, left=True, right=True, labelsize=14, length=5)
ax24.tick_params(which='minor', direction='in', bottom=True, top=False, left=True, right=True, length=2)
ax24.set_xlabel('Visit Time (Minutes)', fontsize=14)
ax24.text(1.02, 0.5, 'Visit 35', transform=ax24.transAxes, 
         va='center', ha='left', fontsize=14)

'''
# This adds one label for the entire left side of the figure
fig.text(0.04, 0.5, 'Clock Angle ($^o$)', 
         va='center', rotation='vertical', fontsize=26)

# Remove the ylabel from all blocks and add this to ax1
ax1.set_title('Clock Angle ($^\circ$)', loc='left', fontsize=18, pad=20)
'''

saveloc = ('/Users/hannah/OneDrive - Lancaster University/aurora/waterfall_bigger.jpg') 
#saveloc = (f'/Users/hannah/OneDrive - Lancaster University/aurora/median_KH_dusk_V_vs_region_power_{error}.jpg')
plt.savefig(saveloc,bbox_inches='tight',dpi=400)






