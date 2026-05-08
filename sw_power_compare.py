#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 28 17:40:46 2024

@author: hannah

this script produces plots of clock angle/pressure/magnetic field strength/reconnection voltge/kelvin-helmholtz vs power
takes total powers per region and calcualtes medians and errors for the data as well as for solar wind data

applies R^2 analysis to the fitted data
"""

import matplotlib.pyplot as plt
import numpy as np
#import spiceypy as sp
import math
import pandas as pd
import datetime as dt
from datetime import datetime
import matplotlib.dates as mdates
from scipy.optimize import curve_fit
import matplotlib.ticker as plticker
from pyestimate import sin_param_estimate
import matplotlib.ticker
from matplotlib.patches import Patch
import matplotlib.ticker as ticker
from matplotlib.ticker import LogLocator, LogFormatter
from scipy.stats import f
from numpy.linalg import lstsq

root_folder = '/Users/hannah/OneDrive - Lancaster University/aurora/python_scripts/dataframes/'
root_saves = '/Users/hannah/OneDrive - Lancaster University/aurora/corrections_graphs/'

error='20'
error_plot='no_error'
plotting ='clock_median'
key ='tang'

# 'main' region definitions
powers02_df = pd.read_csv(root_folder+'new_normalised_powers_02_CR_adj.csv', delimiter=',') #march / CR_adj
powers03_df = pd.read_csv(root_folder+'new_normalised_powers_03_CR_adj.csv', delimiter=',')
powers04_df = pd.read_csv(root_folder+'new_normalised_powers_04_CR_adj.csv', delimiter=',')
powers05_df = pd.read_csv(root_folder+'new_normalised_powers_05_CR_adj.csv', delimiter=',')
powers08_df = pd.read_csv(root_folder+'new_normalised_powers_08_CR_adj.csv', delimiter=',')
powers09_df = pd.read_csv(root_folder+'new_normalised_powers_09_CR_adj.csv', delimiter=',')
powers10_df = pd.read_csv(root_folder+'new_normalised_powers_10_CR_adj.csv', delimiter=',')
powers11_df = pd.read_csv(root_folder+'new_normalised_powers_11_CR_adj.csv', delimiter=',')
powers16_df = pd.read_csv(root_folder+'new_normalised_powers_16_CR_adj.csv', delimiter=',')
powers17_df = pd.read_csv(root_folder+'new_normalised_powers_17_CR_adj.csv', delimiter=',')
powers18_df = pd.read_csv(root_folder+'new_normalised_powers_18_CR_adj.csv', delimiter=',')
powers19_df = pd.read_csv(root_folder+'new_normalised_powers_19_CR_adj.csv', delimiter=',')
powers20_df = pd.read_csv(root_folder+'new_normalised_powers_20_CR_adj.csv', delimiter=',')
powers21_df = pd.read_csv(root_folder+'new_normalised_powers_21_CR_adj.csv', delimiter=',')
powers24_df = pd.read_csv(root_folder+'new_normalised_powers_24_CR_adj.csv', delimiter=',')
powers25_df = pd.read_csv(root_folder+'new_normalised_powers_25_CR_adj.csv', delimiter=',')
powers27_df = pd.read_csv(root_folder+'new_normalised_powers_27_CR_adj.csv', delimiter=',')
powers28_df = pd.read_csv(root_folder+'new_normalised_powers_28_CR_adj.csv', delimiter=',')
powers34_df = pd.read_csv(root_folder+'new_normalised_powers_34_CR_adj.csv', delimiter=',')
powers35_df = pd.read_csv(root_folder+'new_normalised_powers_35_CR_adj.csv', delimiter=',')

# high cml
powers01_df = pd.read_csv(root_folder+'new_normalised_powers_01_CR_adj.csv', delimiter=',')
powers12_df = pd.read_csv(root_folder+'new_normalised_powers_12_CR_adj.csv', delimiter=',')
powers15_df = pd.read_csv(root_folder+'new_normalised_powers_15_CR_adj.csv', delimiter=',')
powers26_df = pd.read_csv(root_folder+'new_normalised_powers_26_CR_adj_26update.csv', delimiter=',')
powers26_df_old = pd.read_csv(root_folder+'new_normalised_powers_26_CR_adj.csv', delimiter=',')


cml_01 = powers01_df['CML'].to_numpy()
swirl_01 = powers01_df['Corrected_Swirl_Powers'].to_numpy()
swirl_err01 = powers01_df['Swirl_Error'].to_numpy()
dusk_err01 = powers01_df['Dusk_Active_Error'].to_numpy()
#noon_err01 = powers01_df['Noon_Active_Error'].to_numpy()
dusk_01 = powers01_df['Corrected_Dusk_Active_Powers'].to_numpy()
noon_01 = powers01_df['Corrected_Noon_Active_Powers'].to_numpy()

cml_12 = powers12_df['CML'].to_numpy()
swirl_12 = powers12_df['Corrected_Swirl_Powers'].to_numpy()
swirl_err12 = powers12_df['Swirl_Error'].to_numpy()
dusk_err12 = powers12_df['Dusk_Active_Error'].to_numpy()
#noon_err12 = powers12_df['Noon_Active_Error'].to_numpy()
dusk_12 = powers12_df['Corrected_Dusk_Active_Powers'].to_numpy()
noon_12 = powers12_df['Corrected_Noon_Active_Powers'].to_numpy()

cml_15 = powers15_df['CML'].to_numpy()
swirl_15 = powers15_df['Corrected_Swirl_Powers'].to_numpy()
swirl_err15 = powers15_df['Swirl_Error'].to_numpy()
dusk_err15 = powers15_df['Dusk_Active_Error'].to_numpy()
#noon_err15 = powers15_df['Noon_Active_Error'].to_numpy()
dusk_15 = powers15_df['Corrected_Dusk_Active_Powers'].to_numpy()
noon_15 = powers15_df['Corrected_Noon_Active_Powers'].to_numpy()

cml_26 = powers26_df['CML'].to_numpy()
swirl_26 = powers26_df['Corrected_Swirl_Powers'].to_numpy()
dusk_26 = powers26_df['Corrected_Dusk_Active_Powers'].to_numpy()
noon_26 = powers26_df['Corrected_Noon_Active_Powers'].to_numpy()
#noon_err26 = powers26_df['Noon_Active_Error'].to_numpy()]
swirl_26_old = powers26_df_old['Corrected_Swirl_Powers'].to_numpy()
dusk_26_old = powers26_df_old['Corrected_Dusk_Active_Powers'].to_numpy()
noon_26_old = powers26_df_old['Corrected_Noon_Active_Powers'].to_numpy()

swirl_err26 = powers26_df['Swirl_Error'].to_numpy()
dusk_err26 = powers26_df['Dusk_Active_Error'].to_numpy()

cml_02 = powers02_df['CML'].to_numpy()
swirl_02 = powers02_df['Corrected_Swirl_Powers'].to_numpy()
swirl_err02 = powers02_df['Swirl_Error'].to_numpy()
dusk_err02 = powers02_df['Dusk_Active_Error'].to_numpy()
#noon_err02 = powers02_df['Noon_Active_Error'].to_numpy()
dusk_02 = powers02_df['Corrected_Dusk_Active_Powers'].to_numpy()
noon_02 = powers02_df['Corrected_Noon_Active_Powers'].to_numpy()

cml_03 = powers03_df['CML'].to_numpy()
swirl_03 = powers03_df['Corrected_Swirl_Powers'].to_numpy()
swirl_err03 = powers03_df['Swirl_Error'].to_numpy()
dusk_err03 = powers03_df['Dusk_Active_Error'].to_numpy()
#noon_err03 = powers03_df['Noon_Active_Error'].to_numpy()
dusk_03 = powers03_df['Corrected_Dusk_Active_Powers'].to_numpy()
noon_03 = powers03_df['Corrected_Noon_Active_Powers'].to_numpy()

cml_04 = powers04_df['CML'].to_numpy()
swirl_04 = powers04_df['Corrected_Swirl_Powers'].to_numpy()
swirl_err04 = powers04_df['Swirl_Error'].to_numpy()
dusk_err04 = powers04_df['Dusk_Active_Error'].to_numpy()
#noon_err04 = powers04_df['Noon_Active_Error'].to_numpy()
dusk_04 = powers04_df['Corrected_Dusk_Active_Powers'].to_numpy()
noon_04 = powers04_df['Corrected_Noon_Active_Powers'].to_numpy()

cml_05 = powers05_df['CML'].to_numpy()
swirl_05 = powers05_df['Corrected_Swirl_Powers'].to_numpy()
swirl_err05 = powers05_df['Swirl_Error'].to_numpy()
dusk_err05 = powers05_df['Dusk_Active_Error'].to_numpy()
#noon_err05 = powers05_df['Noon_Active_Error'].to_numpy()
dusk_05 = powers05_df['Corrected_Dusk_Active_Powers'].to_numpy()
noon_05 = powers05_df['Corrected_Noon_Active_Powers'].to_numpy()

cml_08 = powers08_df['CML'].to_numpy()
swirl_08 = powers08_df['Corrected_Swirl_Powers'].to_numpy()
swirl_err08 = powers08_df['Swirl_Error'].to_numpy()
dusk_err08 = powers08_df['Dusk_Active_Error'].to_numpy()
#noon_err08 = powers08_df['Noon_Active_Error'].to_numpy()
dusk_08 = powers08_df['Corrected_Dusk_Active_Powers'].to_numpy()
noon_08 = powers08_df['Corrected_Noon_Active_Powers'].to_numpy()

cml_09 = powers09_df['CML'].to_numpy()
swirl_err09 = powers09_df['Swirl_Error'].to_numpy()
dusk_err09 = powers09_df['Dusk_Active_Error'].to_numpy()
#noon_err09 = powers09_df['Noon_Active_Error'].to_numpy()
swirl_09 = powers09_df['Corrected_Swirl_Powers'].to_numpy()
dusk_09 = powers09_df['Corrected_Dusk_Active_Powers'].to_numpy()
noon_09 = powers09_df['Corrected_Noon_Active_Powers'].to_numpy()

cml_10 = powers10_df['CML'].to_numpy()
swirl_err10 = powers10_df['Swirl_Error'].to_numpy()
dusk_err10 = powers10_df['Dusk_Active_Error'].to_numpy()
#noon_err10 = powers10_df['Noon_Active_Error'].to_numpy()
swirl_10 = powers10_df['Corrected_Swirl_Powers'].to_numpy()
dusk_10 = powers10_df['Corrected_Dusk_Active_Powers'].to_numpy()
noon_10 = powers10_df['Corrected_Noon_Active_Powers'].to_numpy()

cml_11 = powers11_df['CML'].to_numpy()
swirl_err11 = powers11_df['Swirl_Error'].to_numpy()
dusk_err11 = powers11_df['Dusk_Active_Error'].to_numpy()
#noon_err11 = powers11_df['Noon_Active_Error'].to_numpy()
swirl_11 = powers11_df['Corrected_Swirl_Powers'].to_numpy()
dusk_11 = powers11_df['Corrected_Dusk_Active_Powers'].to_numpy()
noon_11 = powers11_df['Corrected_Noon_Active_Powers'].to_numpy()

cml_16 = powers16_df['CML'].to_numpy()
swirl_err16 = powers16_df['Swirl_Error'].to_numpy()
dusk_err16 = powers16_df['Dusk_Active_Error'].to_numpy()
#noon_err16 = powers16_df['Noon_Active_Error'].to_numpy()
swirl_16 = powers16_df['Corrected_Swirl_Powers'].to_numpy()
dusk_16 = powers16_df['Corrected_Dusk_Active_Powers'].to_numpy()
noon_16 = powers16_df['Corrected_Noon_Active_Powers'].to_numpy()

cml_17 = powers17_df['CML'].to_numpy()
swirl_17 = powers17_df['Corrected_Swirl_Powers'].to_numpy()
swirl_err17 = powers17_df['Swirl_Error'].to_numpy()
dusk_err17 = powers17_df['Dusk_Active_Error'].to_numpy()
#noon_err17 = powers17_df['Noon_Active_Error'].to_numpy()
dusk_17 = powers17_df['Corrected_Dusk_Active_Powers'].to_numpy()
noon_17 = powers17_df['Corrected_Noon_Active_Powers'].to_numpy()

cml_18 = powers18_df['CML'].to_numpy()
swirl_err18 = powers18_df['Swirl_Error'].to_numpy()
dusk_err18 = powers18_df['Dusk_Active_Error'].to_numpy()
#noon_err18 = powers18_df['Noon_Active_Error'].to_numpy()
swirl_18 = powers18_df['Corrected_Swirl_Powers'].to_numpy()
dusk_18 = powers18_df['Corrected_Dusk_Active_Powers'].to_numpy()
noon_18 = powers18_df['Corrected_Noon_Active_Powers'].to_numpy()

cml_19 = powers19_df['CML'].to_numpy()
swirl_err19 = powers19_df['Swirl_Error'].to_numpy()
dusk_err19 = powers19_df['Dusk_Active_Error'].to_numpy()
#noon_err19 = powers19_df['Noon_Active_Error'].to_numpy()
swirl_19 = powers19_df['Corrected_Swirl_Powers'].to_numpy()
dusk_19 = powers19_df['Corrected_Dusk_Active_Powers'].to_numpy()
noon_19 = powers19_df['Corrected_Noon_Active_Powers'].to_numpy()

cml_20 = powers20_df['CML'].to_numpy()
swirl_err20 = powers20_df['Swirl_Error'].to_numpy()
dusk_err20 = powers20_df['Dusk_Active_Error'].to_numpy()
#noon_err20 = powers20_df['Noon_Active_Error'].to_numpy()
swirl_20 = powers20_df['Corrected_Swirl_Powers'].to_numpy()
dusk_20 = powers20_df['Corrected_Dusk_Active_Powers'].to_numpy()
noon_20 = powers20_df['Corrected_Noon_Active_Powers'].to_numpy()

cml_21 = powers21_df['CML'].to_numpy()
swirl_err21 = powers21_df['Swirl_Error'].to_numpy()
dusk_err21 = powers21_df['Dusk_Active_Error'].to_numpy()
#noon_err21 = powers21_df['Noon_Active_Error'].to_numpy()
swirl_21 = powers21_df['Corrected_Swirl_Powers'].to_numpy()
dusk_21 = powers21_df['Corrected_Dusk_Active_Powers'].to_numpy()
noon_21 = powers21_df['Corrected_Noon_Active_Powers'].to_numpy()

cml_24 = powers24_df['CML'].to_numpy()
swirl_err24 = powers24_df['Swirl_Error'].to_numpy()
dusk_err24 = powers24_df['Dusk_Active_Error'].to_numpy()
#noon_err24 = powers24_df['Noon_Active_Error'].to_numpy()
swirl_24 = powers24_df['Corrected_Swirl_Powers'].to_numpy()
dusk_24 = powers24_df['Corrected_Dusk_Active_Powers'].to_numpy()
noon_24 = powers24_df['Corrected_Noon_Active_Powers'].to_numpy()

cml_25 = powers25_df['CML'].to_numpy()
swirl_err25 = powers25_df['Swirl_Error'].to_numpy()
dusk_err25 = powers25_df['Dusk_Active_Error'].to_numpy()
#noon_err25 = powers25_df['Noon_Active_Error'].to_numpy()
swirl_25 = powers25_df['Corrected_Swirl_Powers'].to_numpy()
dusk_25 = powers25_df['Corrected_Dusk_Active_Powers'].to_numpy()
noon_25 = powers25_df['Corrected_Noon_Active_Powers'].to_numpy()

cml_27 = powers27_df['CML'].to_numpy()
swirl_err27 = powers27_df['Swirl_Error'].to_numpy()
dusk_err27 = powers27_df['Dusk_Active_Error'].to_numpy()
#noon_err27 = powers27_df['Noon_Active_Error'].to_numpy()
swirl_27 = powers27_df['Corrected_Swirl_Powers'].to_numpy()
dusk_27 = powers27_df['Corrected_Dusk_Active_Powers'].to_numpy()
noon_27 = powers27_df['Corrected_Noon_Active_Powers'].to_numpy()

cml_28 = powers28_df['CML'].to_numpy()
swirl_err28 = powers28_df['Swirl_Error'].to_numpy()
dusk_err28 = powers28_df['Dusk_Active_Error'].to_numpy()
#noon_err28 = powers28_df['Noon_Active_Error'].to_numpy()
swirl_28 = powers28_df['Corrected_Swirl_Powers'].to_numpy()
dusk_28 = powers28_df['Corrected_Dusk_Active_Powers'].to_numpy()
noon_28 = powers28_df['Corrected_Noon_Active_Powers'].to_numpy()

cml_34 = powers34_df['CML'].to_numpy()
swirl_err34 = powers34_df['Swirl_Error'].to_numpy()
dusk_err34 = powers34_df['Dusk_Active_Error'].to_numpy()
#noon_err34 = powers34_df['Noon_Active_Error'].to_numpy()
swirl_34 = powers34_df['Corrected_Swirl_Powers'].to_numpy()
dusk_34 = powers34_df['Corrected_Dusk_Active_Powers'].to_numpy()
noon_34 = powers34_df['Corrected_Noon_Active_Powers'].to_numpy()

cml_35 = powers35_df['CML'].to_numpy()
swirl_err35 = powers35_df['Swirl_Error'].to_numpy()
dusk_err35 = powers35_df['Dusk_Active_Error'].to_numpy()
#noon_err35 = powers35_df['Noon_Active_Error'].to_numpy()
swirl_35 = powers35_df['Corrected_Swirl_Powers'].to_numpy()
dusk_35 = powers35_df['Corrected_Dusk_Active_Powers'].to_numpy()
noon_35 = powers35_df['Corrected_Noon_Active_Powers'].to_numpy()




visit_times = pd.read_csv(root_folder+'visit_times.csv')

# data for each visit
visit_02_data = pd.read_csv(root_folder+'visit_02_data_aug.csv', delimiter=',')
visit_03_data = pd.read_csv(root_folder+'visit_03_data_aug.csv', delimiter=',')
visit_04_data = pd.read_csv(root_folder+'visit_04_data_aug.csv', delimiter=',')
visit_05_data = pd.read_csv(root_folder+'visit_05_data_aug.csv', delimiter=',')
visit_08_data = pd.read_csv(root_folder+'visit_08_data_aug.csv', delimiter=',')
visit_09_data = pd.read_csv(root_folder+'visit_09_data_aug.csv', delimiter=',')
visit_10_data = pd.read_csv(root_folder+'visit_10_data_aug.csv', delimiter=',')
visit_11_data = pd.read_csv(root_folder+'visit_11_data_aug.csv', delimiter=',')
visit_16_data = pd.read_csv(root_folder+'visit_16_data_aug.csv', delimiter=',')
visit_17_data = pd.read_csv(root_folder+'visit_17_data_aug.csv', delimiter=',')
visit_18_data = pd.read_csv(root_folder+'visit_18_data_aug.csv', delimiter=',')
visit_19_data = pd.read_csv(root_folder+'visit_19_data_aug.csv', delimiter=',')
visit_20_data = pd.read_csv(root_folder+'visit_20_data_aug.csv', delimiter=',')
visit_21_data = pd.read_csv(root_folder+'visit_21_data_aug.csv', delimiter=',')
visit_24_data = pd.read_csv(root_folder+'visit_24_data_aug.csv', delimiter=',')
visit_25_data = pd.read_csv(root_folder+'visit_25_data_aug.csv', delimiter=',')
visit_27_data = pd.read_csv(root_folder+'visit_27_data_aug.csv', delimiter=',')
visit_28_data = pd.read_csv(root_folder+'visit_28_data_aug.csv', delimiter=',')
visit_34_data = pd.read_csv(root_folder+'visit_34_data_aug.csv', delimiter=',')
visit_35_data = pd.read_csv(root_folder+'visit_35_data_aug.csv', delimiter=',')


# # high cml - 01, 12, 15, 26
visit_01_data = pd.read_csv(root_folder+'visit_01_data_aug.csv', delimiter=',')
visit_12_data = pd.read_csv(root_folder+'visit_12_data_aug.csv', delimiter=',')
visit_15_data = pd.read_csv(root_folder+'visit_15_data_aug.csv', delimiter=',')
visit_26_data = pd.read_csv(root_folder+'visit_26_data_aug.csv', delimiter=',')


# dataframes of +% error 
visit_02_plus = pd.read_csv(f'{root_folder}visit_02_times_plus_error_{error}_aug.csv', delimiter=',')
visit_03_plus = pd.read_csv(f'{root_folder}visit_03_times_plus_error_{error}_aug.csv', delimiter=',')
visit_04_plus = pd.read_csv(f'{root_folder}visit_04_times_plus_error_{error}_aug.csv', delimiter=',')
visit_05_plus = pd.read_csv(f'{root_folder}visit_05_times_plus_error_{error}_aug.csv', delimiter=',')
visit_08_plus = pd.read_csv(f'{root_folder}visit_08_times_plus_error_{error}_aug.csv', delimiter=',')
visit_09_plus = pd.read_csv(f'{root_folder}visit_09_times_plus_error_{error}_aug.csv', delimiter=',')
visit_10_plus = pd.read_csv(f'{root_folder}visit_10_times_plus_error_{error}_aug.csv', delimiter=',')
visit_11_plus = pd.read_csv(f'{root_folder}visit_11_times_plus_error_{error}_aug.csv', delimiter=',')
visit_16_plus = pd.read_csv(f'{root_folder}visit_16_times_plus_error_{error}_aug.csv', delimiter=',')
visit_17_plus = pd.read_csv(f'{root_folder}visit_17_times_plus_error_{error}_aug.csv', delimiter=',')
visit_18_plus = pd.read_csv(f'{root_folder}visit_18_times_plus_error_{error}_aug.csv', delimiter=',')
visit_19_plus = pd.read_csv(f'{root_folder}visit_19_times_plus_error_{error}_aug.csv', delimiter=',')
visit_20_plus = pd.read_csv(f'{root_folder}visit_20_times_plus_error_{error}_aug.csv', delimiter=',')
visit_21_plus = pd.read_csv(f'{root_folder}visit_21_times_plus_error_{error}_aug.csv', delimiter=',')
visit_24_plus = pd.read_csv(f'{root_folder}visit_24_times_plus_error_{error}_aug.csv', delimiter=',')
visit_25_plus = pd.read_csv(f'{root_folder}visit_25_times_plus_error_{error}_aug.csv', delimiter=',')
visit_27_plus = pd.read_csv(f'{root_folder}visit_27_times_plus_error_{error}_aug.csv', delimiter=',')
visit_28_plus = pd.read_csv(f'{root_folder}visit_28_times_plus_error_{error}_aug.csv', delimiter=',')
visit_34_plus = pd.read_csv(f'{root_folder}visit_34_times_plus_error_{error}_aug.csv', delimiter=',')
visit_35_plus = pd.read_csv(f'{root_folder}visit_35_times_plus_error_{error}_aug.csv', delimiter=',')

# high cml
visit_01_plus = pd.read_csv(f'{root_folder}visit_01_times_plus_error_{error}_aug.csv', delimiter=',')
visit_12_plus = pd.read_csv(f'{root_folder}visit_12_times_plus_error_{error}_aug.csv', delimiter=',')
visit_15_plus = pd.read_csv(f'{root_folder}visit_15_times_plus_error_{error}_aug.csv', delimiter=',')
visit_26_plus = pd.read_csv(f'{root_folder}visit_26_times_plus_error_{error}_aug.csv', delimiter=',')

# dataframes of -% error
visit_02_minus = pd.read_csv(f'{root_folder}visit_02_times_minus_error_{error}_aug.csv', delimiter=',')
visit_03_minus = pd.read_csv(f'{root_folder}visit_03_times_minus_error_{error}_aug.csv', delimiter=',')
visit_04_minus = pd.read_csv(f'{root_folder}visit_04_times_minus_error_{error}_aug.csv', delimiter=',')
visit_05_minus = pd.read_csv(f'{root_folder}visit_05_times_minus_error_{error}_aug.csv', delimiter=',')
visit_08_minus = pd.read_csv(f'{root_folder}visit_08_times_minus_error_{error}_aug.csv', delimiter=',')
visit_09_minus = pd.read_csv(f'{root_folder}visit_09_times_minus_error_{error}_aug.csv', delimiter=',')
visit_10_minus = pd.read_csv(f'{root_folder}visit_10_times_minus_error_{error}_aug.csv', delimiter=',')
visit_11_minus = pd.read_csv(f'{root_folder}visit_11_times_minus_error_{error}_aug.csv', delimiter=',')
visit_16_minus = pd.read_csv(f'{root_folder}visit_16_times_minus_error_{error}_aug.csv', delimiter=',')
visit_17_minus = pd.read_csv(f'{root_folder}visit_17_times_minus_error_{error}_aug.csv', delimiter=',')
visit_18_minus = pd.read_csv(f'{root_folder}visit_18_times_minus_error_{error}_aug.csv', delimiter=',')
visit_19_minus = pd.read_csv(f'{root_folder}visit_19_times_minus_error_{error}_aug.csv', delimiter=',')
visit_20_minus = pd.read_csv(f'{root_folder}visit_20_times_minus_error_{error}_aug.csv', delimiter=',')
visit_21_minus = pd.read_csv(f'{root_folder}visit_21_times_minus_error_{error}_aug.csv', delimiter=',')
visit_24_minus = pd.read_csv(f'{root_folder}visit_24_times_minus_error_{error}_aug.csv', delimiter=',')
visit_25_minus = pd.read_csv(f'{root_folder}visit_25_times_minus_error_{error}_aug.csv', delimiter=',')
visit_27_minus = pd.read_csv(f'{root_folder}visit_27_times_minus_error_{error}_aug.csv', delimiter=',')
visit_28_minus = pd.read_csv(f'{root_folder}visit_28_times_minus_error_{error}_aug.csv', delimiter=',')
visit_34_minus = pd.read_csv(f'{root_folder}visit_34_times_minus_error_{error}_aug.csv', delimiter=',')
visit_35_minus = pd.read_csv(f'{root_folder}visit_35_times_minus_error_{error}_aug.csv', delimiter=',')


# high cmls
visit_01_minus = pd.read_csv(f'{root_folder}visit_01_times_minus_error_{error}_aug.csv', delimiter=',')
visit_12_minus = pd.read_csv(f'{root_folder}visit_12_times_minus_error_{error}_aug.csv', delimiter=',')
visit_15_minus = pd.read_csv(f'{root_folder}visit_15_times_minus_error_{error}_aug.csv', delimiter=',')
visit_26_minus = pd.read_csv(f'{root_folder}visit_26_times_minus_error_{error}_aug.csv', delimiter=',')


# ---------- total mag strength ----------

def total_mag(Bx,By,Bz):
    Btot = np.sqrt((Bx**2)+(By**2)+(Bz**2))
    return Btot

bx_01 = visit_01_data['Bx'].to_numpy()
by_01 = visit_01_data['By'].to_numpy()
bz_01 = visit_01_data['Bz'].to_numpy()

btot_01 = total_mag(bx_01,by_01,bz_01)

bx_12 = visit_12_data['Bx'].to_numpy()
by_12 = visit_12_data['By'].to_numpy()
bz_12 = visit_12_data['Bz'].to_numpy()

btot_12 = total_mag(bx_12,by_12,bz_12)

bx_15 = visit_15_data['Bx'].to_numpy()
by_15 = visit_15_data['By'].to_numpy()
bz_15 = visit_15_data['Bz'].to_numpy()

btot_15 = total_mag(bx_15,by_15,bz_15)

bx_26 = visit_26_data['Bx'].to_numpy()
by_26 = visit_26_data['By'].to_numpy()
bz_26 = visit_26_data['Bz'].to_numpy()

btot_26 = total_mag(bx_26,by_26,bz_26)



bx_02 = visit_02_data['Bx'].to_numpy()
by_02 = visit_02_data['By'].to_numpy()
bz_02 = visit_02_data['Bz'].to_numpy()

btot_02 = total_mag(bx_02,by_02,bz_02)


bx_03 = visit_03_data['Bx'].to_numpy()
by_03 = visit_03_data['By'].to_numpy()
bz_03 = visit_03_data['Bz'].to_numpy()

btot_03 = total_mag(bx_03,by_03,bz_03)


bx_04 = visit_04_data['Bx'].to_numpy()
by_04 = visit_04_data['By'].to_numpy()
bz_04 = visit_04_data['Bz'].to_numpy()

btot_04 = total_mag(bx_04,by_04,bz_04)


bx_05 = visit_05_data['Bx'].to_numpy()
by_05 = visit_05_data['By'].to_numpy()
bz_05 = visit_05_data['Bz'].to_numpy()

btot_05 = total_mag(bx_05,by_05,bz_05)


bx_08 = visit_08_data['Bx'].to_numpy()
by_08 = visit_08_data['By'].to_numpy()
bz_08 = visit_08_data['Bz'].to_numpy()

btot_08 = total_mag(bx_08,by_08,bz_08)


bx_09 = visit_09_data['Bx'].to_numpy()
by_09 = visit_09_data['By'].to_numpy()
bz_09 = visit_09_data['Bz'].to_numpy()

btot_09 = total_mag(bx_09,by_09,bz_09)


bx_10 = visit_10_data['Bx'].to_numpy()
by_10 = visit_10_data['By'].to_numpy()
bz_10 = visit_10_data['Bz'].to_numpy()

btot_10 = total_mag(bx_10,by_10,bz_10)


bx_11 = visit_11_data['Bx'].to_numpy()
by_11 = visit_11_data['By'].to_numpy()
bz_11 = visit_11_data['Bz'].to_numpy()

btot_11 = total_mag(bx_11,by_11,bz_11)


bx_16 = visit_16_data['Bx'].to_numpy()
by_16 = visit_16_data['By'].to_numpy()
bz_16 = visit_16_data['Bz'].to_numpy()

btot_16 = total_mag(bx_16,by_16,bz_16)


bx_17 = visit_17_data['Bx'].to_numpy()
by_17 = visit_17_data['By'].to_numpy()
bz_17 = visit_17_data['Bz'].to_numpy()

btot_17 = total_mag(bx_17,by_17,bz_17)


bx_18 = visit_18_data['Bx'].to_numpy()
by_18 = visit_18_data['By'].to_numpy()
bz_18 = visit_18_data['Bz'].to_numpy()

btot_18 = total_mag(bx_18,by_18,bz_18)


bx_19 = visit_19_data['Bx'].to_numpy()
by_19 = visit_19_data['By'].to_numpy()
bz_19 = visit_19_data['Bz'].to_numpy()

btot_19 = total_mag(bx_19,by_19,bz_19)


bx_20 = visit_20_data['Bx'].to_numpy()
by_20 = visit_20_data['By'].to_numpy()
bz_20 = visit_20_data['Bz'].to_numpy()

btot_20 = total_mag(bx_20,by_20,bz_20)


bx_21 = visit_21_data['Bx'].to_numpy()
by_21 = visit_21_data['By'].to_numpy()
bz_21 = visit_21_data['Bz'].to_numpy()

btot_21 = total_mag(bx_21,by_21,bz_21)


bx_24 = visit_24_data['Bx'].to_numpy()
by_24 = visit_24_data['By'].to_numpy()
bz_24 = visit_24_data['Bz'].to_numpy()

btot_24 = total_mag(bx_24,by_24,bz_24)


bx_25 = visit_25_data['Bx'].to_numpy()
by_25 = visit_25_data['By'].to_numpy()
bz_25 = visit_25_data['Bz'].to_numpy()

btot_25 = total_mag(bx_25,by_25,bz_25)


bx_27 = visit_27_data['Bx'].to_numpy()
by_27 = visit_27_data['By'].to_numpy()
bz_27 = visit_27_data['Bz'].to_numpy()

btot_27 = total_mag(bx_27,by_27,bz_27)


bx_28 = visit_28_data['Bx'].to_numpy()
by_28 = visit_28_data['By'].to_numpy()
bz_28 = visit_28_data['Bz'].to_numpy()

btot_28 = total_mag(bx_28,by_28,bz_28)


bx_34 = visit_34_data['Bx'].to_numpy()
by_34 = visit_34_data['By'].to_numpy()
bz_34 = visit_34_data['Bz'].to_numpy()

btot_34 = total_mag(bx_34,by_34,bz_34)


bx_35 = visit_35_data['Bx'].to_numpy()
by_35 = visit_35_data['By'].to_numpy()
bz_35 = visit_35_data['Bz'].to_numpy()

btot_35 = total_mag(bx_35,by_35,bz_35)


# ---- + % ------

bx_01_p = visit_01_plus['PLUS_Bx'].to_numpy()
by_01_p = visit_01_plus['PLUS_By'].to_numpy()
bz_01_p = visit_01_plus['PLUS_Bz'].to_numpy()

btot_01_p = total_mag(bx_01_p,by_01_p,bz_01_p)

bx_12_p = visit_12_plus['PLUS_Bx'].to_numpy()
by_12_p = visit_12_plus['PLUS_By'].to_numpy()
bz_12_p = visit_12_plus['PLUS_Bz'].to_numpy()

btot_12_p = total_mag(bx_12_p,by_12_p,bz_12_p)

bx_15_p = visit_15_plus['PLUS_Bx'].to_numpy()
by_15_p = visit_15_plus['PLUS_By'].to_numpy()
bz_15_p = visit_15_plus['PLUS_Bz'].to_numpy()

btot_15_p = total_mag(bx_15_p,by_15_p,bz_15_p)

bx_26_p = visit_26_plus['PLUS_Bx'].to_numpy()
by_26_p = visit_26_plus['PLUS_By'].to_numpy()
bz_26_p = visit_26_plus['PLUS_Bz'].to_numpy()

btot_26_p = total_mag(bx_26_p,by_26_p,bz_26_p)



bx_02_p = visit_02_plus['PLUS_Bx'].to_numpy()
by_02_p = visit_02_plus['PLUS_By'].to_numpy()
bz_02_p = visit_02_plus['PLUS_Bz'].to_numpy()

btot_02_p = total_mag(bx_02_p,by_02_p,bz_02_p)


bx_03_p = visit_03_plus['PLUS_Bx'].to_numpy()
by_03_p = visit_03_plus['PLUS_By'].to_numpy()
bz_03_p = visit_03_plus['PLUS_Bz'].to_numpy()

btot_03_p = total_mag(bx_03_p,by_03_p,bz_03_p)


bx_04_p = visit_04_plus['PLUS_Bx'].to_numpy()
by_04_p = visit_04_plus['PLUS_By'].to_numpy()
bz_04_p = visit_04_plus['PLUS_Bz'].to_numpy()

btot_04_p = total_mag(bx_04_p,by_04_p,bz_04_p)


bx_05_p = visit_05_plus['PLUS_Bx'].to_numpy()
by_05_p = visit_05_plus['PLUS_By'].to_numpy()
bz_05_p = visit_05_plus['PLUS_Bz'].to_numpy()

btot_05_p = total_mag(bx_05_p,by_05_p,bz_05_p)

bx_08_p = visit_08_plus['PLUS_Bx'].to_numpy()
by_08_p = visit_08_plus['PLUS_By'].to_numpy()
bz_08_p = visit_08_plus['PLUS_Bz'].to_numpy()

btot_08_p = total_mag(bx_08_p,by_08_p,bz_08_p)

bx_09_p = visit_09_plus['PLUS_Bx'].to_numpy()
by_09_p = visit_09_plus['PLUS_By'].to_numpy()
bz_09_p = visit_09_plus['PLUS_Bz'].to_numpy()

btot_09_p = total_mag(bx_09_p,by_09_p,bz_09_p)

bx_10_p = visit_10_plus['PLUS_Bx'].to_numpy()
by_10_p = visit_10_plus['PLUS_By'].to_numpy()
bz_10_p = visit_10_plus['PLUS_Bz'].to_numpy()

btot_10_p = total_mag(bx_10_p,by_10_p,bz_10_p)

bx_11_p = visit_11_plus['PLUS_Bx'].to_numpy()
by_11_p = visit_11_plus['PLUS_By'].to_numpy()
bz_11_p = visit_11_plus['PLUS_Bz'].to_numpy()

btot_11_p = total_mag(bx_11_p,by_11_p,bz_11_p)

bx_16_p = visit_16_plus['PLUS_Bx'].to_numpy()
by_16_p = visit_16_plus['PLUS_By'].to_numpy()
bz_16_p = visit_16_plus['PLUS_Bz'].to_numpy()

btot_16_p = total_mag(bx_16_p,by_16_p,bz_16_p)

bx_17_p = visit_17_plus['PLUS_Bx'].to_numpy()
by_17_p = visit_17_plus['PLUS_By'].to_numpy()
bz_17_p = visit_17_plus['PLUS_Bz'].to_numpy()

btot_17_p = total_mag(bx_17_p,by_17_p,bz_17_p)

bx_18_p = visit_18_plus['PLUS_Bx'].to_numpy()
by_18_p = visit_18_plus['PLUS_By'].to_numpy()
bz_18_p = visit_18_plus['PLUS_Bz'].to_numpy()

btot_18_p = total_mag(bx_18_p,by_18_p,bz_18_p)

bx_19_p = visit_19_plus['PLUS_Bx'].to_numpy()
by_19_p = visit_19_plus['PLUS_By'].to_numpy()
bz_19_p = visit_19_plus['PLUS_Bz'].to_numpy()

btot_19_p = total_mag(bx_19_p,by_19_p,bz_19_p)

bx_20_p = visit_20_plus['PLUS_Bx'].to_numpy()
by_20_p = visit_20_plus['PLUS_By'].to_numpy()
bz_20_p = visit_20_plus['PLUS_Bz'].to_numpy()

btot_20_p = total_mag(bx_20_p,by_20_p,bz_20_p)

bx_21_p = visit_21_plus['PLUS_Bx'].to_numpy()
by_21_p = visit_21_plus['PLUS_By'].to_numpy()
bz_21_p = visit_21_plus['PLUS_Bz'].to_numpy()

btot_21_p = total_mag(bx_21_p,by_21_p,bz_21_p)

bx_24_p = visit_24_plus['PLUS_Bx'].to_numpy()
by_24_p = visit_24_plus['PLUS_By'].to_numpy()
bz_24_p = visit_24_plus['PLUS_Bz'].to_numpy()

btot_24_p = total_mag(bx_24_p,by_24_p,bz_24_p)

bx_25_p = visit_25_plus['PLUS_Bx'].to_numpy()
by_25_p = visit_25_plus['PLUS_By'].to_numpy()
bz_25_p = visit_25_plus['PLUS_Bz'].to_numpy()

btot_25_p = total_mag(bx_25_p,by_25_p,bz_25_p)

bx_27_p = visit_27_plus['PLUS_Bx'].to_numpy()
by_27_p = visit_27_plus['PLUS_By'].to_numpy()
bz_27_p = visit_27_plus['PLUS_Bz'].to_numpy()

btot_27_p = total_mag(bx_27_p,by_27_p,bz_27_p)

bx_28_p = visit_28_plus['PLUS_Bx'].to_numpy()
by_28_p = visit_28_plus['PLUS_By'].to_numpy()
bz_28_p = visit_28_plus['PLUS_Bz'].to_numpy()

btot_28_p = total_mag(bx_28_p,by_28_p,bz_28_p)

bx_34_p = visit_34_plus['PLUS_Bx'].to_numpy()
by_34_p = visit_34_plus['PLUS_By'].to_numpy()
bz_34_p = visit_34_plus['PLUS_Bz'].to_numpy()

btot_34_p = total_mag(bx_34_p,by_34_p,bz_34_p)

bx_35_p = visit_35_plus['PLUS_Bx'].to_numpy()
by_35_p = visit_35_plus['PLUS_By'].to_numpy()
bz_35_p = visit_35_plus['PLUS_Bz'].to_numpy()

btot_35_p = total_mag(bx_35_p,by_35_p,bz_35_p)

# -------- - % ----------

bx_01_m = visit_01_minus['MINUS_Bx'].to_numpy()
by_01_m = visit_01_minus['MINUS_By'].to_numpy()
bz_01_m = visit_01_minus['MINUS_Bz'].to_numpy()

btot_01_m = total_mag(bx_01_m,by_01_m,bz_01_m)

bx_12_m = visit_12_minus['MINUS_Bx'].to_numpy()
by_12_m = visit_12_minus['MINUS_By'].to_numpy()
bz_12_m = visit_12_minus['MINUS_Bz'].to_numpy()

btot_12_m = total_mag(bx_12_m,by_12_m,bz_12_m)

bx_15_m = visit_15_minus['MINUS_Bx'].to_numpy()
by_15_m = visit_15_minus['MINUS_By'].to_numpy()
bz_15_m = visit_15_minus['MINUS_Bz'].to_numpy()

btot_15_m = total_mag(bx_15_m,by_15_m,bz_15_m)

bx_26_m = visit_26_minus['MINUS_Bx'].to_numpy()
by_26_m = visit_26_minus['MINUS_By'].to_numpy()
bz_26_m = visit_26_minus['MINUS_Bz'].to_numpy()

btot_26_m = total_mag(bx_26_m,by_26_m,bz_26_m)



bx_02_m = visit_02_minus['MINUS_Bx'].to_numpy()
by_02_m = visit_02_minus['MINUS_By'].to_numpy()
bz_02_m = visit_02_minus['MINUS_Bz'].to_numpy()

btot_02_m = total_mag(bx_02_m,by_02_m,bz_02_m)

bx_03_m = visit_03_minus['MINUS_Bx'].to_numpy()
by_03_m = visit_03_minus['MINUS_By'].to_numpy()
bz_03_m = visit_03_minus['MINUS_Bz'].to_numpy()

btot_03_m = total_mag(bx_03_m,by_03_m,bz_03_m)

bx_04_m = visit_04_minus['MINUS_Bx'].to_numpy()
by_04_m = visit_04_minus['MINUS_By'].to_numpy()
bz_04_m = visit_04_minus['MINUS_Bz'].to_numpy()

btot_04_m = total_mag(bx_04_m,by_04_m,bz_04_m)

bx_05_m = visit_05_minus['MINUS_Bx'].to_numpy()
by_05_m = visit_05_minus['MINUS_By'].to_numpy()
bz_05_m = visit_05_minus['MINUS_Bz'].to_numpy()

btot_05_m = total_mag(bx_05_m,by_05_m,bz_05_m)

bx_08_m = visit_08_minus['MINUS_Bx'].to_numpy()
by_08_m = visit_08_minus['MINUS_By'].to_numpy()
bz_08_m = visit_08_minus['MINUS_Bz'].to_numpy()

btot_08_m = total_mag(bx_08_m,by_08_m,bz_08_m)

bx_09_m = visit_09_minus['MINUS_Bx'].to_numpy()
by_09_m = visit_09_minus['MINUS_By'].to_numpy()
bz_09_m = visit_09_minus['MINUS_Bz'].to_numpy()

btot_09_m = total_mag(bx_09_m,by_09_m,bz_09_m)

bx_10_m = visit_10_minus['MINUS_Bx'].to_numpy()
by_10_m = visit_10_minus['MINUS_By'].to_numpy()
bz_10_m = visit_10_minus['MINUS_Bz'].to_numpy()

btot_10_m = total_mag(bx_10_m,by_10_m,bz_10_m)

bx_11_m = visit_11_minus['MINUS_Bx'].to_numpy()
by_11_m = visit_11_minus['MINUS_By'].to_numpy()
bz_11_m = visit_11_minus['MINUS_Bz'].to_numpy()

btot_11_m = total_mag(bx_11_m,by_11_m,bz_11_m)

bx_16_m = visit_16_minus['MINUS_Bx'].to_numpy()
by_16_m = visit_16_minus['MINUS_By'].to_numpy()
bz_16_m = visit_16_minus['MINUS_Bz'].to_numpy()

btot_16_m = total_mag(bx_16_m,by_16_m,bz_16_m)

bx_17_m = visit_17_minus['MINUS_Bx'].to_numpy()
by_17_m = visit_17_minus['MINUS_By'].to_numpy()
bz_17_m = visit_17_minus['MINUS_Bz'].to_numpy()

btot_17_m = total_mag(bx_17_m,by_17_m,bz_17_m)

bx_18_m = visit_18_minus['MINUS_Bx'].to_numpy()
by_18_m = visit_18_minus['MINUS_By'].to_numpy()
bz_18_m = visit_18_minus['MINUS_Bz'].to_numpy()

btot_18_m = total_mag(bx_18_m,by_18_m,bz_18_m)

bx_19_m = visit_19_minus['MINUS_Bx'].to_numpy()
by_19_m = visit_19_minus['MINUS_By'].to_numpy()
bz_19_m = visit_19_minus['MINUS_Bz'].to_numpy()

btot_19_m = total_mag(bx_19_m,by_19_m,bz_19_m)

bx_20_m = visit_20_minus['MINUS_Bx'].to_numpy()
by_20_m = visit_20_minus['MINUS_By'].to_numpy()
bz_20_m = visit_20_minus['MINUS_Bz'].to_numpy()

btot_20_m = total_mag(bx_20_m,by_20_m,bz_20_m)

bx_21_m = visit_21_minus['MINUS_Bx'].to_numpy()
by_21_m = visit_21_minus['MINUS_By'].to_numpy()
bz_21_m = visit_21_minus['MINUS_Bz'].to_numpy()

btot_21_m = total_mag(bx_21_m,by_21_m,bz_21_m)

bx_24_m = visit_24_minus['MINUS_Bx'].to_numpy()
by_24_m = visit_24_minus['MINUS_By'].to_numpy()
bz_24_m = visit_24_minus['MINUS_Bz'].to_numpy()

btot_24_m = total_mag(bx_24_m,by_24_m,bz_24_m)

bx_25_m = visit_25_minus['MINUS_Bx'].to_numpy()
by_25_m = visit_25_minus['MINUS_By'].to_numpy()
bz_25_m = visit_25_minus['MINUS_Bz'].to_numpy()

btot_25_m = total_mag(bx_25_m,by_25_m,bz_25_m)

bx_27_m = visit_27_minus['MINUS_Bx'].to_numpy()
by_27_m = visit_27_minus['MINUS_By'].to_numpy()
bz_27_m = visit_27_minus['MINUS_Bz'].to_numpy()

btot_27_m = total_mag(bx_27_m,by_27_m,bz_27_m)

bx_28_m = visit_28_minus['MINUS_Bx'].to_numpy()
by_28_m = visit_28_minus['MINUS_By'].to_numpy()
bz_28_m = visit_28_minus['MINUS_Bz'].to_numpy()

btot_28_m = total_mag(bx_28_m,by_28_m,bz_28_m)

bx_34_m = visit_34_minus['MINUS_Bx'].to_numpy()
by_34_m = visit_34_minus['MINUS_By'].to_numpy()
bz_34_m = visit_34_minus['MINUS_Bz'].to_numpy()

btot_34_m = total_mag(bx_34_m,by_34_m,bz_34_m)

bx_35_m = visit_35_minus['MINUS_Bx'].to_numpy()
by_35_m = visit_35_minus['MINUS_By'].to_numpy()
bz_35_m = visit_35_minus['MINUS_Bz'].to_numpy()

btot_35_m = total_mag(bx_35_m,by_35_m,bz_35_m)


# ---------- clock angles --------------

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


# extract clock angle plus error for each visit
clock_02_plus = visit_02_plus['PLUS_Clock_Angle'].to_numpy()
clock_03_plus = visit_03_plus['PLUS_Clock_Angle'].to_numpy()
clock_04_plus = visit_04_plus['PLUS_Clock_Angle'].to_numpy()
clock_05_plus = visit_05_plus['PLUS_Clock_Angle'].to_numpy()
clock_08_plus = visit_08_plus['PLUS_Clock_Angle'].to_numpy()
clock_09_plus = visit_09_plus['PLUS_Clock_Angle'].to_numpy()
clock_10_plus = visit_10_plus['PLUS_Clock_Angle'].to_numpy()
clock_11_plus = visit_11_plus['PLUS_Clock_Angle'].to_numpy()
clock_16_plus = visit_16_plus['PLUS_Clock_Angle'].to_numpy()
clock_17_plus = visit_17_plus['PLUS_Clock_Angle'].to_numpy()
clock_18_plus = visit_18_plus['PLUS_Clock_Angle'].to_numpy()
clock_19_plus = visit_19_plus['PLUS_Clock_Angle'].to_numpy()
clock_20_plus = visit_20_plus['PLUS_Clock_Angle'].to_numpy()
clock_21_plus = visit_21_plus['PLUS_Clock_Angle'].to_numpy()
clock_24_plus = visit_24_plus['PLUS_Clock_Angle'].to_numpy()
clock_25_plus = visit_25_plus['PLUS_Clock_Angle'].to_numpy()
clock_27_plus = visit_27_plus['PLUS_Clock_Angle'].to_numpy()
clock_28_plus = visit_28_plus['PLUS_Clock_Angle'].to_numpy()
clock_34_plus = visit_34_plus['PLUS_Clock_Angle'].to_numpy()
clock_35_plus = visit_35_plus['PLUS_Clock_Angle'].to_numpy()

# high cmls
clock_01_plus = visit_01_plus['PLUS_Clock_Angle'].to_numpy()
clock_12_plus = visit_12_plus['PLUS_Clock_Angle'].to_numpy()
clock_15_plus = visit_15_plus['PLUS_Clock_Angle'].to_numpy()
clock_26_plus = visit_26_plus['PLUS_Clock_Angle'].to_numpy()

# extract clock angle minus error for each visit
clock_02_minus = visit_02_minus['MINUS_Clock_Angle'].to_numpy()
clock_03_minus = visit_03_minus['MINUS_Clock_Angle'].to_numpy()
clock_04_minus = visit_04_minus['MINUS_Clock_Angle'].to_numpy()
clock_05_minus = visit_05_minus['MINUS_Clock_Angle'].to_numpy()
clock_08_minus = visit_08_minus['MINUS_Clock_Angle'].to_numpy()
clock_09_minus = visit_09_minus['MINUS_Clock_Angle'].to_numpy()
clock_10_minus = visit_10_minus['MINUS_Clock_Angle'].to_numpy()
clock_11_minus = visit_11_minus['MINUS_Clock_Angle'].to_numpy()
clock_16_minus = visit_16_minus['MINUS_Clock_Angle'].to_numpy()
clock_17_minus = visit_17_minus['MINUS_Clock_Angle'].to_numpy()
clock_18_minus = visit_18_minus['MINUS_Clock_Angle'].to_numpy()
clock_19_minus = visit_19_minus['MINUS_Clock_Angle'].to_numpy()
clock_20_minus = visit_20_minus['MINUS_Clock_Angle'].to_numpy()
clock_21_minus = visit_21_minus['MINUS_Clock_Angle'].to_numpy()
clock_24_minus = visit_24_minus['MINUS_Clock_Angle'].to_numpy()
clock_25_minus = visit_25_minus['MINUS_Clock_Angle'].to_numpy()
clock_27_minus = visit_27_minus['MINUS_Clock_Angle'].to_numpy()
clock_28_minus = visit_28_minus['MINUS_Clock_Angle'].to_numpy()
clock_34_minus = visit_34_minus['MINUS_Clock_Angle'].to_numpy()
clock_35_minus = visit_35_minus['MINUS_Clock_Angle'].to_numpy()

# high cmls
clock_01_minus = visit_01_minus['MINUS_Clock_Angle'].to_numpy()
clock_12_minus = visit_12_minus['MINUS_Clock_Angle'].to_numpy()
clock_15_minus = visit_15_minus['MINUS_Clock_Angle'].to_numpy()
clock_26_minus = visit_26_minus['MINUS_Clock_Angle'].to_numpy()


# ------------- pressures ----------------

# extract pressure for each visit
pressure_02 = visit_02_data['Sw_Pressure'].to_numpy()
pressure_03 = visit_03_data['Sw_Pressure'].to_numpy()
pressure_04 = visit_04_data['Sw_Pressure'].to_numpy()
pressure_05 = visit_05_data['Sw_Pressure'].to_numpy()
pressure_08 = visit_08_data['Sw_Pressure'].to_numpy()
pressure_09 = visit_09_data['Sw_Pressure'].to_numpy()
pressure_10 = visit_10_data['Sw_Pressure'].to_numpy()
pressure_11 = visit_11_data['Sw_Pressure'].to_numpy()
pressure_16 = visit_16_data['Sw_Pressure'].to_numpy()
pressure_17 = visit_17_data['Sw_Pressure'].to_numpy()
pressure_18 = visit_18_data['Sw_Pressure'].to_numpy()
pressure_19 = visit_19_data['Sw_Pressure'].to_numpy()
pressure_20 = visit_20_data['Sw_Pressure'].to_numpy()
pressure_21 = visit_21_data['Sw_Pressure'].to_numpy()
pressure_24 = visit_24_data['Sw_Pressure'].to_numpy()
pressure_25 = visit_25_data['Sw_Pressure'].to_numpy()
pressure_27 = visit_27_data['Sw_Pressure'].to_numpy()
pressure_28 = visit_28_data['Sw_Pressure'].to_numpy()
pressure_34 = visit_34_data['Sw_Pressure'].to_numpy()
pressure_35 = visit_35_data['Sw_Pressure'].to_numpy()

#pressure_13 = visit_13_data['Sw_Pressure'].to_numpy()

# high cml
pressure_01 = visit_01_data['Sw_Pressure'].to_numpy()
pressure_12 = visit_12_data['Sw_Pressure'].to_numpy()
pressure_15 = visit_15_data['Sw_Pressure'].to_numpy()
pressure_26 = visit_26_data['Sw_Pressure'].to_numpy()


pressure_02_plus = visit_02_plus['PLUS_Pressure'].to_numpy()
pressure_03_plus = visit_03_plus['PLUS_Pressure'].to_numpy()
pressure_04_plus = visit_04_plus['PLUS_Pressure'].to_numpy()
pressure_05_plus = visit_05_plus['PLUS_Pressure'].to_numpy()
pressure_08_plus = visit_08_plus['PLUS_Pressure'].to_numpy()
pressure_09_plus = visit_09_plus['PLUS_Pressure'].to_numpy()
pressure_10_plus = visit_10_plus['PLUS_Pressure'].to_numpy()
pressure_11_plus = visit_11_plus['PLUS_Pressure'].to_numpy()
pressure_16_plus = visit_16_plus['PLUS_Pressure'].to_numpy()
pressure_17_plus = visit_17_plus['PLUS_Pressure'].to_numpy()
pressure_18_plus = visit_18_plus['PLUS_Pressure'].to_numpy()
pressure_19_plus = visit_19_plus['PLUS_Pressure'].to_numpy()
pressure_20_plus = visit_20_plus['PLUS_Pressure'].to_numpy()
pressure_21_plus = visit_21_plus['PLUS_Pressure'].to_numpy()
pressure_24_plus = visit_24_plus['PLUS_Pressure'].to_numpy()
pressure_25_plus = visit_25_plus['PLUS_Pressure'].to_numpy()
pressure_27_plus = visit_27_plus['PLUS_Pressure'].to_numpy()
pressure_28_plus = visit_28_plus['PLUS_Pressure'].to_numpy()
pressure_34_plus = visit_34_plus['PLUS_Pressure'].to_numpy()
pressure_35_plus = visit_35_plus['PLUS_Pressure'].to_numpy()

# high cmls
pressure_01_plus = visit_01_plus['PLUS_Pressure'].to_numpy()
pressure_12_plus = visit_12_plus['PLUS_Pressure'].to_numpy()
pressure_15_plus = visit_15_plus['PLUS_Pressure'].to_numpy()
pressure_26_plus = visit_26_plus['PLUS_Pressure'].to_numpy()

pressure_02_minus = visit_02_minus['MINUS_Pressure'].to_numpy()
pressure_03_minus = visit_03_minus['MINUS_Pressure'].to_numpy()
pressure_04_minus = visit_04_minus['MINUS_Pressure'].to_numpy()
pressure_05_minus = visit_05_minus['MINUS_Pressure'].to_numpy()
pressure_08_minus = visit_08_minus['MINUS_Pressure'].to_numpy()
pressure_09_minus = visit_09_minus['MINUS_Pressure'].to_numpy()
pressure_10_minus = visit_10_minus['MINUS_Pressure'].to_numpy()
pressure_11_minus = visit_11_minus['MINUS_Pressure'].to_numpy()
pressure_16_minus = visit_16_minus['MINUS_Pressure'].to_numpy()
pressure_17_minus = visit_17_minus['MINUS_Pressure'].to_numpy()
pressure_18_minus = visit_18_minus['MINUS_Pressure'].to_numpy()
pressure_19_minus = visit_19_minus['MINUS_Pressure'].to_numpy()
pressure_20_minus = visit_20_minus['MINUS_Pressure'].to_numpy()
pressure_21_minus = visit_21_minus['MINUS_Pressure'].to_numpy()
pressure_24_minus = visit_24_minus['MINUS_Pressure'].to_numpy()
pressure_25_minus = visit_25_minus['MINUS_Pressure'].to_numpy()
pressure_27_minus = visit_27_minus['MINUS_Pressure'].to_numpy()
pressure_28_minus = visit_28_minus['MINUS_Pressure'].to_numpy()
pressure_34_minus = visit_34_minus['MINUS_Pressure'].to_numpy()
pressure_35_minus = visit_35_minus['MINUS_Pressure'].to_numpy()

# high cmls
pressure_01_minus = visit_01_minus['MINUS_Pressure'].to_numpy()
pressure_12_minus = visit_12_minus['MINUS_Pressure'].to_numpy()
pressure_15_minus = visit_15_minus['MINUS_Pressure'].to_numpy()
pressure_26_minus = visit_26_minus['MINUS_Pressure'].to_numpy()


# ------------ perpendicular magnetic field -------------

b_perp_02 = visit_02_data['B_Perp'].to_numpy()
b_perp_03 = visit_03_data['B_Perp'].to_numpy()
b_perp_04 = visit_04_data['B_Perp'].to_numpy()
b_perp_05 = visit_05_data['B_Perp'].to_numpy()
b_perp_08 = visit_08_data['B_Perp'].to_numpy()
b_perp_09 = visit_09_data['B_Perp'].to_numpy()
b_perp_10 = visit_10_data['B_Perp'].to_numpy()
b_perp_11 = visit_11_data['B_Perp'].to_numpy()
b_perp_16 = visit_16_data['B_Perp'].to_numpy()
b_perp_17 = visit_17_data['B_Perp'].to_numpy()
b_perp_18 = visit_18_data['B_Perp'].to_numpy()
b_perp_19 = visit_19_data['B_Perp'].to_numpy()
b_perp_20 = visit_20_data['B_Perp'].to_numpy()
b_perp_21 = visit_21_data['B_Perp'].to_numpy()
b_perp_24 = visit_24_data['B_Perp'].to_numpy()
b_perp_25 = visit_25_data['B_Perp'].to_numpy()
b_perp_27 = visit_27_data['B_Perp'].to_numpy()
b_perp_28 = visit_28_data['B_Perp'].to_numpy()
b_perp_34 = visit_34_data['B_Perp'].to_numpy()
b_perp_35 = visit_35_data['B_Perp'].to_numpy()

# high cml
b_perp_01 = visit_01_data['B_Perp'].to_numpy()
b_perp_12 = visit_12_data['B_Perp'].to_numpy()
b_perp_15 = visit_15_data['B_Perp'].to_numpy()
b_perp_26 = visit_26_data['B_Perp'].to_numpy()

b_perp_02_plus = visit_02_plus['PLUS_B_Perp'].to_numpy()
b_perp_03_plus = visit_03_plus['PLUS_B_Perp'].to_numpy()
b_perp_04_plus = visit_04_plus['PLUS_B_Perp'].to_numpy()
b_perp_05_plus = visit_05_plus['PLUS_B_Perp'].to_numpy()
b_perp_08_plus = visit_08_plus['PLUS_B_Perp'].to_numpy()
b_perp_09_plus = visit_09_plus['PLUS_B_Perp'].to_numpy()
b_perp_10_plus = visit_10_plus['PLUS_B_Perp'].to_numpy()
b_perp_11_plus = visit_11_plus['PLUS_B_Perp'].to_numpy()
b_perp_16_plus = visit_16_plus['PLUS_B_Perp'].to_numpy()
b_perp_17_plus = visit_17_plus['PLUS_B_Perp'].to_numpy()
b_perp_18_plus = visit_18_plus['PLUS_B_Perp'].to_numpy()
b_perp_19_plus = visit_19_plus['PLUS_B_Perp'].to_numpy()
b_perp_20_plus = visit_20_plus['PLUS_B_Perp'].to_numpy()
b_perp_21_plus = visit_21_plus['PLUS_B_Perp'].to_numpy()
b_perp_24_plus = visit_24_plus['PLUS_B_Perp'].to_numpy()
b_perp_25_plus = visit_25_plus['PLUS_B_Perp'].to_numpy()
b_perp_27_plus = visit_27_plus['PLUS_B_Perp'].to_numpy()
b_perp_28_plus = visit_28_plus['PLUS_B_Perp'].to_numpy()
b_perp_34_plus = visit_34_plus['PLUS_B_Perp'].to_numpy()
b_perp_35_plus = visit_35_plus['PLUS_B_Perp'].to_numpy()

# high cmls
b_perp_01_plus = visit_01_plus['PLUS_B_Perp'].to_numpy()
b_perp_12_plus = visit_12_plus['PLUS_B_Perp'].to_numpy()
b_perp_15_plus = visit_15_plus['PLUS_B_Perp'].to_numpy()
b_perp_26_plus = visit_26_plus['PLUS_B_Perp'].to_numpy()

b_perp_02_minus = visit_02_minus['MINUS_B_Perp'].to_numpy()
b_perp_03_minus = visit_03_minus['MINUS_B_Perp'].to_numpy()
b_perp_04_minus = visit_04_minus['MINUS_B_Perp'].to_numpy()
b_perp_05_minus = visit_05_minus['MINUS_B_Perp'].to_numpy()
b_perp_08_minus = visit_08_minus['MINUS_B_Perp'].to_numpy()
b_perp_09_minus = visit_09_minus['MINUS_B_Perp'].to_numpy()
b_perp_10_minus = visit_10_minus['MINUS_B_Perp'].to_numpy()
b_perp_11_minus = visit_11_minus['MINUS_B_Perp'].to_numpy()
b_perp_16_minus = visit_16_minus['MINUS_B_Perp'].to_numpy()
b_perp_17_minus = visit_17_minus['MINUS_B_Perp'].to_numpy()
b_perp_18_minus = visit_18_minus['MINUS_B_Perp'].to_numpy()
b_perp_19_minus = visit_19_minus['MINUS_B_Perp'].to_numpy()
b_perp_20_minus = visit_20_minus['MINUS_B_Perp'].to_numpy()
b_perp_21_minus = visit_21_minus['MINUS_B_Perp'].to_numpy()
b_perp_24_minus = visit_24_minus['MINUS_B_Perp'].to_numpy()
b_perp_25_minus = visit_25_minus['MINUS_B_Perp'].to_numpy()
b_perp_27_minus = visit_27_minus['MINUS_B_Perp'].to_numpy()
b_perp_28_minus = visit_28_minus['MINUS_B_Perp'].to_numpy()
b_perp_34_minus = visit_34_minus['MINUS_B_Perp'].to_numpy()
b_perp_35_minus = visit_35_minus['MINUS_B_Perp'].to_numpy()

# high cmls
b_perp_01_minus = visit_01_minus['MINUS_B_Perp'].to_numpy()
b_perp_12_minus = visit_12_minus['MINUS_B_Perp'].to_numpy()
b_perp_15_minus = visit_15_minus['MINUS_B_Perp'].to_numpy()
b_perp_26_minus = visit_26_minus['MINUS_B_Perp'].to_numpy()


# ----------- low latitude reconnection voltage ------------

LL_02 = visit_02_data["LL_Recon_V"]
LL_03 = visit_03_data["LL_Recon_V"]
LL_04 = visit_04_data["LL_Recon_V"]
LL_05 = visit_05_data["LL_Recon_V"]
LL_08 = visit_08_data["LL_Recon_V"]
LL_09 = visit_09_data["LL_Recon_V"]
LL_10 = visit_10_data["LL_Recon_V"]
LL_11 = visit_11_data["LL_Recon_V"]
LL_16 = visit_16_data["LL_Recon_V"]
LL_17 = visit_17_data["LL_Recon_V"]
LL_18 = visit_18_data["LL_Recon_V"]
LL_19 = visit_19_data["LL_Recon_V"]
LL_20 = visit_20_data["LL_Recon_V"]
LL_21 = visit_21_data["LL_Recon_V"]
LL_24 = visit_24_data["LL_Recon_V"]
LL_25 = visit_25_data["LL_Recon_V"]
LL_27 = visit_27_data["LL_Recon_V"]
LL_28 = visit_28_data["LL_Recon_V"]
LL_34 = visit_34_data["LL_Recon_V"]
LL_35 = visit_35_data["LL_Recon_V"]

# high cml
LL_01 = visit_01_data["LL_Recon_V"]
LL_12 = visit_12_data["LL_Recon_V"]
LL_15 = visit_15_data["LL_Recon_V"]
LL_26 = visit_26_data["LL_Recon_V"]

LL_02_plus = visit_02_plus["PLUS_LL_rec_V"]
LL_03_plus = visit_03_plus["PLUS_LL_rec_V"]
LL_04_plus = visit_04_plus["PLUS_LL_rec_V"]
LL_05_plus = visit_05_plus["PLUS_LL_rec_V"]
LL_08_plus = visit_08_plus["PLUS_LL_rec_V"]
LL_09_plus = visit_09_plus["PLUS_LL_rec_V"]
LL_10_plus = visit_10_plus["PLUS_LL_rec_V"]
LL_11_plus = visit_11_plus["PLUS_LL_rec_V"]
LL_16_plus = visit_16_plus["PLUS_LL_rec_V"]
LL_17_plus = visit_17_plus["PLUS_LL_rec_V"]
LL_18_plus = visit_18_plus["PLUS_LL_rec_V"]
LL_19_plus = visit_19_plus["PLUS_LL_rec_V"]
LL_20_plus = visit_20_plus["PLUS_LL_rec_V"]
LL_21_plus = visit_21_plus["PLUS_LL_rec_V"]
LL_24_plus = visit_24_plus["PLUS_LL_rec_V"]
LL_25_plus = visit_25_plus["PLUS_LL_rec_V"]
LL_27_plus = visit_27_plus["PLUS_LL_rec_V"]
LL_28_plus = visit_28_plus["PLUS_LL_rec_V"]
LL_34_plus = visit_34_plus["PLUS_LL_rec_V"]
LL_35_plus = visit_35_plus["PLUS_LL_rec_V"]

# high cmls
LL_01_plus = visit_01_plus["PLUS_LL_rec_V"]
LL_12_plus = visit_12_plus["PLUS_LL_rec_V"]
LL_15_plus = visit_15_plus["PLUS_LL_rec_V"]
LL_26_plus = visit_26_plus["PLUS_LL_rec_V"]

LL_02_minus = visit_02_minus["MINUS_LL_rec_V"]
LL_03_minus = visit_03_minus["MINUS_LL_rec_V"]
LL_04_minus = visit_04_minus["MINUS_LL_rec_V"]
LL_05_minus = visit_05_minus["MINUS_LL_rec_V"]
LL_08_minus = visit_08_minus["MINUS_LL_rec_V"]
LL_09_minus = visit_09_minus["MINUS_LL_rec_V"]
LL_10_minus = visit_10_minus["MINUS_LL_rec_V"]
LL_11_minus = visit_11_minus["MINUS_LL_rec_V"]
LL_16_minus = visit_16_minus["MINUS_LL_rec_V"]
LL_17_minus = visit_17_minus["MINUS_LL_rec_V"]
LL_18_minus = visit_18_minus["MINUS_LL_rec_V"]
LL_19_minus = visit_19_minus["MINUS_LL_rec_V"]
LL_20_minus = visit_20_minus["MINUS_LL_rec_V"]
LL_21_minus = visit_21_minus["MINUS_LL_rec_V"]
LL_24_minus = visit_24_minus["MINUS_LL_rec_V"]
LL_25_minus = visit_25_minus["MINUS_LL_rec_V"]
LL_27_minus = visit_27_minus["MINUS_LL_rec_V"]
LL_28_minus = visit_28_minus["MINUS_LL_rec_V"]
LL_34_minus = visit_34_minus["MINUS_LL_rec_V"]
LL_35_minus = visit_35_minus["MINUS_LL_rec_V"]

# high cmls
LL_01_minus = visit_01_minus["MINUS_LL_rec_V"]
LL_12_minus = visit_12_minus["MINUS_LL_rec_V"]
LL_15_minus = visit_15_minus["MINUS_LL_rec_V"]
LL_26_minus = visit_26_minus["MINUS_LL_rec_V"]


# --------- high latitude (+By) reconnection ------------


HL_02_pos = visit_02_data["HL_Recon_V_By_Plus"]
HL_03_pos = visit_03_data["HL_Recon_V_By_Plus"]
HL_04_pos = visit_04_data["HL_Recon_V_By_Plus"]
HL_05_pos = visit_05_data["HL_Recon_V_By_Plus"]
HL_08_pos = visit_08_data["HL_Recon_V_By_Plus"]
HL_09_pos = visit_09_data["HL_Recon_V_By_Plus"]
HL_10_pos = visit_10_data["HL_Recon_V_By_Plus"]
HL_11_pos = visit_11_data["HL_Recon_V_By_Plus"]
HL_16_pos = visit_16_data["HL_Recon_V_By_Plus"]
HL_17_pos = visit_17_data["HL_Recon_V_By_Plus"]
HL_18_pos = visit_18_data["HL_Recon_V_By_Plus"]
HL_19_pos = visit_19_data["HL_Recon_V_By_Plus"]
HL_20_pos = visit_20_data["HL_Recon_V_By_Plus"]
HL_21_pos = visit_21_data["HL_Recon_V_By_Plus"]
HL_24_pos = visit_24_data["HL_Recon_V_By_Plus"]
HL_25_pos = visit_25_data["HL_Recon_V_By_Plus"]
HL_27_pos = visit_27_data["HL_Recon_V_By_Plus"]
HL_28_pos = visit_28_data["HL_Recon_V_By_Plus"]
HL_34_pos = visit_34_data["HL_Recon_V_By_Plus"]
HL_35_pos = visit_35_data["HL_Recon_V_By_Plus"]

# high cml
HL_01_pos = visit_01_data["HL_Recon_V_By_Plus"]
HL_12_pos = visit_12_data["HL_Recon_V_By_Plus"]
HL_15_pos = visit_15_data["HL_Recon_V_By_Plus"]
HL_26_pos = visit_26_data["HL_Recon_V_By_Plus"]


HL_02_plus_pos = visit_02_plus["PLUS_HL_rec_V_pos"]
HL_03_plus_pos = visit_03_plus["PLUS_HL_rec_V_pos"]
HL_04_plus_pos = visit_04_plus["PLUS_HL_rec_V_pos"]
HL_05_plus_pos = visit_05_plus["PLUS_HL_rec_V_pos"]
HL_08_plus_pos = visit_08_plus["PLUS_HL_rec_V_pos"]
HL_09_plus_pos = visit_09_plus["PLUS_HL_rec_V_pos"]
HL_10_plus_pos = visit_10_plus["PLUS_HL_rec_V_pos"]
HL_11_plus_pos = visit_11_plus["PLUS_HL_rec_V_pos"]
HL_16_plus_pos = visit_16_plus["PLUS_HL_rec_V_pos"]
HL_17_plus_pos = visit_17_plus["PLUS_HL_rec_V_pos"]
HL_18_plus_pos = visit_18_plus["PLUS_HL_rec_V_pos"]
HL_19_plus_pos = visit_19_plus["PLUS_HL_rec_V_pos"]
HL_20_plus_pos = visit_20_plus["PLUS_HL_rec_V_pos"]
HL_21_plus_pos = visit_21_plus["PLUS_HL_rec_V_pos"]
HL_24_plus_pos = visit_24_plus["PLUS_HL_rec_V_pos"]
HL_25_plus_pos = visit_25_plus["PLUS_HL_rec_V_pos"]
HL_27_plus_pos = visit_27_plus["PLUS_HL_rec_V_pos"]
HL_28_plus_pos = visit_28_plus["PLUS_HL_rec_V_pos"]
HL_34_plus_pos = visit_34_plus["PLUS_HL_rec_V_pos"]
HL_35_plus_pos = visit_35_plus["PLUS_HL_rec_V_pos"]

# high cmls
HL_01_plus_pos = visit_01_plus["PLUS_HL_rec_V_pos"]
HL_12_plus_pos = visit_12_plus["PLUS_HL_rec_V_pos"]
HL_15_plus_pos = visit_15_plus["PLUS_HL_rec_V_pos"]
HL_26_plus_pos = visit_26_plus["PLUS_HL_rec_V_pos"]


HL_02_minus_pos = visit_02_minus["MINUS_HL_rec_V_pos"]
HL_03_minus_pos = visit_03_minus["MINUS_HL_rec_V_pos"]
HL_04_minus_pos = visit_04_minus["MINUS_HL_rec_V_pos"]
HL_05_minus_pos = visit_05_minus["MINUS_HL_rec_V_pos"]
HL_08_minus_pos = visit_08_minus["MINUS_HL_rec_V_pos"]
HL_09_minus_pos = visit_09_minus["MINUS_HL_rec_V_pos"]
HL_10_minus_pos = visit_10_minus["MINUS_HL_rec_V_pos"]
HL_11_minus_pos = visit_11_minus["MINUS_HL_rec_V_pos"]
HL_16_minus_pos = visit_16_minus["MINUS_HL_rec_V_pos"]
HL_17_minus_pos = visit_17_minus["MINUS_HL_rec_V_pos"]
HL_18_minus_pos = visit_18_minus["MINUS_HL_rec_V_pos"]
HL_19_minus_pos = visit_19_minus["MINUS_HL_rec_V_pos"]
HL_20_minus_pos = visit_20_minus["MINUS_HL_rec_V_pos"]
HL_21_minus_pos = visit_21_minus["MINUS_HL_rec_V_pos"]
HL_24_minus_pos = visit_24_minus["MINUS_HL_rec_V_pos"]
HL_25_minus_pos = visit_25_minus["MINUS_HL_rec_V_pos"]
HL_27_minus_pos = visit_27_minus["MINUS_HL_rec_V_pos"]
HL_28_minus_pos = visit_28_minus["MINUS_HL_rec_V_pos"]
HL_34_minus_pos = visit_34_minus["MINUS_HL_rec_V_pos"]
HL_35_minus_pos = visit_35_minus["MINUS_HL_rec_V_pos"]

# high cmls
HL_01_minus_pos = visit_01_minus["MINUS_HL_rec_V_pos"]
HL_12_minus_pos = visit_12_minus["MINUS_HL_rec_V_pos"]
HL_15_minus_pos = visit_15_minus["MINUS_HL_rec_V_pos"]
HL_26_minus_pos = visit_26_minus["MINUS_HL_rec_V_pos"]


# ------ high latitude (-By) reconnection ----------

HL_02_neg = visit_02_data["HL_Recon_V_By_Neg"]
HL_03_neg = visit_03_data["HL_Recon_V_By_Neg"]
HL_04_neg = visit_04_data["HL_Recon_V_By_Neg"]
HL_05_neg = visit_05_data["HL_Recon_V_By_Neg"]
HL_08_neg = visit_08_data["HL_Recon_V_By_Neg"]
HL_09_neg = visit_09_data["HL_Recon_V_By_Neg"]
HL_10_neg = visit_10_data["HL_Recon_V_By_Neg"]
HL_11_neg = visit_11_data["HL_Recon_V_By_Neg"]
HL_16_neg = visit_16_data["HL_Recon_V_By_Neg"]
HL_17_neg = visit_17_data["HL_Recon_V_By_Neg"]
HL_18_neg = visit_18_data["HL_Recon_V_By_Neg"]
HL_19_neg = visit_19_data["HL_Recon_V_By_Neg"]
HL_20_neg = visit_20_data["HL_Recon_V_By_Neg"]
HL_21_neg = visit_21_data["HL_Recon_V_By_Neg"]
HL_24_neg = visit_24_data["HL_Recon_V_By_Neg"]
HL_25_neg = visit_25_data["HL_Recon_V_By_Neg"]
HL_27_neg = visit_27_data["HL_Recon_V_By_Neg"]
HL_28_neg = visit_28_data["HL_Recon_V_By_Neg"]
HL_34_neg = visit_34_data["HL_Recon_V_By_Neg"]
HL_35_neg = visit_35_data["HL_Recon_V_By_Neg"]

# high cml
HL_01_neg = visit_01_data["HL_Recon_V_By_Neg"]
HL_12_neg = visit_12_data["HL_Recon_V_By_Neg"]
HL_15_neg = visit_15_data["HL_Recon_V_By_Neg"]
HL_26_neg = visit_26_data["HL_Recon_V_By_Neg"]


HL_02_plus_neg = visit_02_plus["PLUS_HL_rec_V_neg"]
HL_03_plus_neg = visit_03_plus["PLUS_HL_rec_V_neg"]
HL_04_plus_neg = visit_04_plus["PLUS_HL_rec_V_neg"]
HL_05_plus_neg = visit_05_plus["PLUS_HL_rec_V_neg"]
HL_08_plus_neg = visit_08_plus["PLUS_HL_rec_V_neg"]
HL_09_plus_neg = visit_09_plus["PLUS_HL_rec_V_neg"]
HL_10_plus_neg = visit_10_plus["PLUS_HL_rec_V_neg"]
HL_11_plus_neg = visit_11_plus["PLUS_HL_rec_V_neg"]
HL_16_plus_neg = visit_16_plus["PLUS_HL_rec_V_neg"]
HL_17_plus_neg = visit_17_plus["PLUS_HL_rec_V_neg"]
HL_18_plus_neg = visit_18_plus["PLUS_HL_rec_V_neg"]
HL_19_plus_neg = visit_19_plus["PLUS_HL_rec_V_neg"]
HL_20_plus_neg = visit_20_plus["PLUS_HL_rec_V_neg"]
HL_21_plus_neg = visit_21_plus["PLUS_HL_rec_V_neg"]
HL_24_plus_neg = visit_24_plus["PLUS_HL_rec_V_neg"]
HL_25_plus_neg = visit_25_plus["PLUS_HL_rec_V_neg"]
HL_27_plus_neg = visit_27_plus["PLUS_HL_rec_V_neg"]
HL_28_plus_neg = visit_28_plus["PLUS_HL_rec_V_neg"]
HL_34_plus_neg = visit_34_plus["PLUS_HL_rec_V_neg"]
HL_35_plus_neg = visit_35_plus["PLUS_HL_rec_V_neg"]

# high cmls
HL_01_plus_neg = visit_01_plus["PLUS_HL_rec_V_neg"]
HL_12_plus_neg = visit_12_plus["PLUS_HL_rec_V_neg"]
HL_15_plus_neg = visit_15_plus["PLUS_HL_rec_V_neg"]
HL_26_plus_neg = visit_26_plus["PLUS_HL_rec_V_neg"]


HL_02_minus_neg = visit_02_minus["MINUS_HL_rec_V_neg"]
HL_03_minus_neg = visit_03_minus["MINUS_HL_rec_V_neg"]
HL_04_minus_neg = visit_04_minus["MINUS_HL_rec_V_neg"]
HL_05_minus_neg = visit_05_minus["MINUS_HL_rec_V_neg"]
HL_08_minus_neg = visit_08_minus["MINUS_HL_rec_V_neg"]
HL_09_minus_neg = visit_09_minus["MINUS_HL_rec_V_neg"]
HL_10_minus_neg = visit_10_minus["MINUS_HL_rec_V_neg"]
HL_11_minus_neg = visit_11_minus["MINUS_HL_rec_V_neg"]
HL_16_minus_neg = visit_16_minus["MINUS_HL_rec_V_neg"]
HL_17_minus_neg = visit_17_minus["MINUS_HL_rec_V_neg"]
HL_18_minus_neg = visit_18_minus["MINUS_HL_rec_V_neg"]
HL_19_minus_neg = visit_19_minus["MINUS_HL_rec_V_neg"]
HL_20_minus_neg = visit_20_minus["MINUS_HL_rec_V_neg"]
HL_21_minus_neg = visit_21_minus["MINUS_HL_rec_V_neg"]
HL_24_minus_neg = visit_24_minus["MINUS_HL_rec_V_neg"]
HL_25_minus_neg = visit_25_minus["MINUS_HL_rec_V_neg"]
HL_27_minus_neg = visit_27_minus["MINUS_HL_rec_V_neg"]
HL_28_minus_neg = visit_28_minus["MINUS_HL_rec_V_neg"]
HL_34_minus_neg = visit_34_minus["MINUS_HL_rec_V_neg"]
HL_35_minus_neg = visit_35_minus["MINUS_HL_rec_V_neg"]

# high cmls
HL_01_minus_neg = visit_01_minus["MINUS_HL_rec_V_neg"]
HL_12_minus_neg = visit_12_minus["MINUS_HL_rec_V_neg"]
HL_15_minus_neg = visit_15_minus["MINUS_HL_rec_V_neg"]
HL_26_minus_neg = visit_26_minus["MINUS_HL_rec_V_neg"]


# -------- Gershman Reconnection Voltage --------------

gersh_02 = visit_02_data["Gershman_Rec"]
gersh_03 = visit_03_data["Gershman_Rec"]
gersh_04 = visit_04_data["Gershman_Rec"]
gersh_05 = visit_05_data["Gershman_Rec"]
gersh_08 = visit_08_data["Gershman_Rec"]
gersh_09 = visit_09_data["Gershman_Rec"]
gersh_10 = visit_10_data["Gershman_Rec"]
gersh_11 = visit_11_data["Gershman_Rec"]
gersh_16 = visit_16_data["Gershman_Rec"]
gersh_17 = visit_17_data["Gershman_Rec_V"]
gersh_18 = visit_18_data["Gershman_Rec"]
gersh_19 = visit_19_data["Gershman_Rec"]
gersh_20 = visit_20_data["Gershman_Rec"]
gersh_21 = visit_21_data["Gershman_Rec"]
gersh_24 = visit_24_data["Gershman_Rec"]
gersh_25 = visit_25_data["Gershman_Rec"]
gersh_27 = visit_27_data["Gershman_Rec"]
gersh_28 = visit_28_data["Gershman_Rec"]
gersh_34 = visit_34_data["Gershman_Rec"]
gersh_35 = visit_35_data["Gershman_Rec"]

gersh_01 = visit_01_data["Gershman_Rec"]
gersh_12 = visit_12_data["Gershman_Rec"]
gersh_15 = visit_15_data["Gershman_Rec"]
gersh_26 = visit_26_data["Gershman_Rec"]

# ------ plus % ------
gersh_02_plus = visit_02_plus["PLUS_Gershman_rec"]
gersh_03_plus = visit_03_plus["PLUS_Gershman_rec"]
gersh_04_plus = visit_04_plus["PLUS_Gershman_rec"]
gersh_05_plus = visit_05_plus["PLUS_Gershman_rec"]
gersh_08_plus = visit_08_plus["PLUS_Gershman_rec"]
gersh_09_plus = visit_09_plus["PLUS_Gershman_rec"]
gersh_10_plus = visit_10_plus["PLUS_Gershman_rec"]
gersh_11_plus = visit_11_plus["PLUS_Gershman_rec"]
gersh_16_plus = visit_16_plus["PLUS_Gershman_rec"]
gersh_17_plus = visit_17_plus["PLUS_Gershman_rec"]
gersh_18_plus = visit_18_plus["PLUS_Gershman_rec"]
gersh_19_plus = visit_19_plus["PLUS_Gershman_rec"]
gersh_20_plus = visit_20_plus["PLUS_Gershman_rec"]
gersh_21_plus = visit_21_plus["PLUS_Gershman_rec"]
gersh_24_plus = visit_24_plus["PLUS_Gershman_rec"]
gersh_25_plus = visit_25_plus["PLUS_Gershman_rec"]
gersh_27_plus = visit_27_plus["PLUS_Gershman_rec"]
gersh_28_plus = visit_28_plus["PLUS_Gershman_rec"]
gersh_34_plus = visit_34_plus["PLUS_Gershman_rec"]
gersh_35_plus = visit_35_plus["PLUS_Gershman_rec"]

gersh_01_plus = visit_01_plus["PLUS_Gershman_rec"]
gersh_12_plus = visit_12_plus["PLUS_Gershman_rec"]
gersh_15_plus = visit_15_plus["PLUS_Gershman_rec"]
gersh_26_plus = visit_26_plus["PLUS_Gershman_rec"]

# ---------- minus % --------

gersh_02_minus = visit_02_minus["MINUS_Gershman_rec"]
gersh_03_minus = visit_03_minus["MINUS_Gershman_rec"]
gersh_04_minus = visit_04_minus["MINUS_Gershman_rec"]
gersh_05_minus = visit_05_minus["MINUS_Gershman_rec"]
gersh_08_minus = visit_08_minus["MINUS_Gershman_rec"]
gersh_09_minus = visit_09_minus["MINUS_Gershman_rec"]
gersh_10_minus = visit_10_minus["MINUS_Gershman_rec"]
gersh_11_minus = visit_11_minus["MINUS_Gershman_rec"]
gersh_16_minus = visit_16_minus["MINUS_Gershman_rec"]
gersh_17_minus = visit_17_minus["MINUS_Gershman_rec"]
gersh_18_minus = visit_18_minus["MINUS_Gershman_rec"]
gersh_19_minus = visit_19_minus["MINUS_Gershman_rec"]
gersh_20_minus = visit_20_minus["MINUS_Gershman_rec"]
gersh_21_minus = visit_21_minus["MINUS_Gershman_rec"]
gersh_24_minus = visit_24_minus["MINUS_Gershman_rec"]
gersh_25_minus = visit_25_minus["MINUS_Gershman_rec"]
gersh_27_minus = visit_27_minus["MINUS_Gershman_rec"]
gersh_28_minus = visit_28_minus["MINUS_Gershman_rec"]
gersh_34_minus = visit_34_minus["MINUS_Gershman_rec"]
gersh_35_minus = visit_35_minus["MINUS_Gershman_rec"]

gersh_01_minus = visit_01_minus["MINUS_Gershman_rec"]
gersh_12_minus = visit_12_minus["MINUS_Gershman_rec"]
gersh_15_minus = visit_15_minus["MINUS_Gershman_rec"]
gersh_26_minus = visit_26_minus["MINUS_Gershman_rec"]


# ---------- Kelvin Helmhotlz ------------

KH_02 = visit_02_data["Kelvin_H"]
KH_03 = visit_03_data["Kelvin_H"]
KH_04 = visit_04_data["Kelvin_H"]
KH_05 = visit_05_data["Kelvin_H"]
KH_08 = visit_08_data["Kelvin_H"]
KH_09 = visit_09_data["Kelvin_H"]
KH_10 = visit_10_data["Kelvin_H"]
KH_11 = visit_11_data["Kelvin_H"]
KH_16 = visit_16_data["Kelvin_H"]
KH_17 = visit_17_data["Kelvin_H"]
KH_18 = visit_18_data["Kelvin_H"]
KH_19 = visit_19_data["Kelvin_H"]
KH_20 = visit_20_data["Kelvin_H"]
KH_21 = visit_21_data["Kelvin_H"]
KH_24 = visit_24_data["Kelvin_H"]
KH_25 = visit_25_data["Kelvin_H"]
KH_27 = visit_27_data["Kelvin_H"]
KH_28 = visit_28_data["Kelvin_H"]
KH_34 = visit_34_data["Kelvin_H"]
KH_35 = visit_35_data["Kelvin_H"]

KH_01 = visit_01_data["Kelvin_H"]
KH_12 = visit_12_data["Kelvin_H"]
KH_15 = visit_15_data["Kelvin_H"]
KH_26 = visit_26_data["Kelvin_H"]

# --------- plus % ---------------

KH_02_plus = visit_02_plus["PLUS_KH"]
KH_03_plus = visit_03_plus["PLUS_KH"]
KH_04_plus = visit_04_plus["PLUS_KH"]
KH_05_plus = visit_05_plus["PLUS_KH"]
KH_08_plus = visit_08_plus["PLUS_KH"]
KH_09_plus = visit_09_plus["PLUS_KH"]
KH_10_plus = visit_10_plus["PLUS_KH"]
KH_11_plus = visit_11_plus["PLUS_KH"]
KH_16_plus = visit_16_plus["PLUS_KH"]
KH_17_plus = visit_17_plus["PLUS_KH"]
KH_18_plus = visit_18_plus["PLUS_KH"]
KH_19_plus = visit_19_plus["PLUS_KH"]
KH_20_plus = visit_20_plus["PLUS_KH"]
KH_21_plus = visit_21_plus["PLUS_KH"]
KH_24_plus = visit_24_plus["PLUS_KH"]
KH_25_plus = visit_25_plus["PLUS_KH"]
KH_27_plus = visit_27_plus["PLUS_KH"]
KH_28_plus = visit_28_plus["PLUS_KH"]
KH_34_plus = visit_34_plus["PLUS_KH"]
KH_35_plus = visit_35_plus["PLUS_KH"]

KH_01_plus = visit_01_plus["PLUS_KH"]
KH_12_plus = visit_12_plus["PLUS_KH"]
KH_15_plus = visit_15_plus["PLUS_KH"]
KH_26_plus = visit_26_plus["PLUS_KH"]

# --------- minus % ---------------

KH_02_minus = visit_02_minus["MINUS_KH"]
KH_03_minus = visit_03_minus["MINUS_KH"]
KH_04_minus = visit_04_minus["MINUS_KH"]
KH_05_minus = visit_05_minus["MINUS_KH"]
KH_08_minus = visit_08_minus["MINUS_KH"]
KH_09_minus = visit_09_minus["MINUS_KH"]
KH_10_minus = visit_10_minus["MINUS_KH"]
KH_11_minus = visit_11_minus["MINUS_KH"]
KH_16_minus = visit_16_minus["MINUS_KH"]
KH_17_minus = visit_17_minus["MINUS_KH"]
KH_18_minus = visit_18_minus["MINUS_KH"]
KH_19_minus = visit_19_minus["MINUS_KH"]
KH_20_minus = visit_20_minus["MINUS_KH"]
KH_21_minus = visit_21_minus["MINUS_KH"]
KH_24_minus = visit_24_minus["MINUS_KH"]
KH_25_minus = visit_25_minus["MINUS_KH"]
KH_27_minus = visit_27_minus["MINUS_KH"]
KH_28_minus = visit_28_minus["MINUS_KH"]
KH_34_minus = visit_34_minus["MINUS_KH"]
KH_35_minus = visit_35_minus["MINUS_KH"]

KH_01_minus = visit_01_minus["MINUS_KH"]
KH_12_minus = visit_12_minus["MINUS_KH"]
KH_15_minus = visit_15_minus["MINUS_KH"]
KH_26_minus = visit_26_minus["MINUS_KH"]

# ---------- Kelvin Helmhotlz - new, dawn ------------

KH_dawn_02 = visit_02_data["Kelvin_H_Dawn"]
KH_dawn_03 = visit_03_data["Kelvin_H_Dawn"]
KH_dawn_04 = visit_04_data["Kelvin_H_Dawn"]
KH_dawn_05 = visit_05_data["Kelvin_H_Dawn"]
KH_dawn_08 = visit_08_data["Kelvin_H_Dawn"]
KH_dawn_09 = visit_09_data["Kelvin_H_Dawn"]
KH_dawn_10 = visit_10_data["Kelvin_H_Dawn"]
KH_dawn_11 = visit_11_data["Kelvin_H_Dawn"]
KH_dawn_16 = visit_16_data["Kelvin_H_Dawn"]
KH_dawn_17 = visit_17_data["Kelvin_H_Dawn"]
KH_dawn_18 = visit_18_data["Kelvin_H_Dawn"]
KH_dawn_19 = visit_19_data["Kelvin_H_Dawn"]
KH_dawn_20 = visit_20_data["Kelvin_H_Dawn"]
KH_dawn_21 = visit_21_data["Kelvin_H_Dawn"]
KH_dawn_24 = visit_24_data["Kelvin_H_Dawn"]
KH_dawn_25 = visit_25_data["Kelvin_H_Dawn"]
KH_dawn_27 = visit_27_data["Kelvin_H_Dawn"]
KH_dawn_28 = visit_28_data["Kelvin_H_Dawn"]
KH_dawn_34 = visit_34_data["Kelvin_H_Dawn"]
KH_dawn_35 = visit_35_data["Kelvin_H_Dawn"]

KH_dawn_01 = visit_01_data["Kelvin_H_Dawn"]
KH_dawn_12 = visit_12_data["Kelvin_H_Dawn"]
KH_dawn_15 = visit_15_data["Kelvin_H_Dawn"]
KH_dawn_26 = visit_26_data["Kelvin_H_Dawn"]

# --------- plus % ---------------

KH_dawn_02_plus = visit_02_plus["PLUS_KH_dawn"]
KH_dawn_03_plus = visit_03_plus["PLUS_KH_dawn"]
KH_dawn_04_plus = visit_04_plus["PLUS_KH_dawn"]
KH_dawn_05_plus = visit_05_plus["PLUS_KH_dawn"]
KH_dawn_08_plus = visit_08_plus["PLUS_KH_dawn"]
KH_dawn_09_plus = visit_09_plus["PLUS_KH_dawn"]
KH_dawn_10_plus = visit_10_plus["PLUS_KH_dawn"]
KH_dawn_11_plus = visit_11_plus["PLUS_KH_dawn"]
KH_dawn_16_plus = visit_16_plus["PLUS_KH_dawn"]
KH_dawn_17_plus = visit_17_plus["PLUS_KH_dawn"]
KH_dawn_18_plus = visit_18_plus["PLUS_KH_dawn"]
KH_dawn_19_plus = visit_19_plus["PLUS_KH_dawn"]
KH_dawn_20_plus = visit_20_plus["PLUS_KH_dawn"]
KH_dawn_21_plus = visit_21_plus["PLUS_KH_dawn"]
KH_dawn_24_plus = visit_24_plus["PLUS_KH_dawn"]
KH_dawn_25_plus = visit_25_plus["PLUS_KH_dawn"]
KH_dawn_27_plus = visit_27_plus["PLUS_KH_dawn"]
KH_dawn_28_plus = visit_28_plus["PLUS_KH_dawn"]
KH_dawn_34_plus = visit_34_plus["PLUS_KH_dawn"]
KH_dawn_35_plus = visit_35_plus["PLUS_KH_dawn"]

KH_dawn_01_plus = visit_01_plus["PLUS_KH_dawn"]
KH_dawn_12_plus = visit_12_plus["PLUS_KH_dawn"]
KH_dawn_15_plus = visit_15_plus["PLUS_KH_dawn"]
KH_dawn_26_plus = visit_26_plus["PLUS_KH_dawn"]

# --------- minus % ---------------

KH_dawn_02_minus = visit_02_minus["MINUS_KH_dawn"]
KH_dawn_03_minus = visit_03_minus["MINUS_KH_dawn"]
KH_dawn_04_minus = visit_04_minus["MINUS_KH_dawn"]
KH_dawn_05_minus = visit_05_minus["MINUS_KH_dawn"]
KH_dawn_08_minus = visit_08_minus["MINUS_KH_dawn"]
KH_dawn_09_minus = visit_09_minus["MINUS_KH_dawn"]
KH_dawn_10_minus = visit_10_minus["MINUS_KH_dawn"]
KH_dawn_11_minus = visit_11_minus["MINUS_KH_dawn"]
KH_dawn_16_minus = visit_16_minus["MINUS_KH_dawn"]
KH_dawn_17_minus = visit_17_minus["MINUS_KH_dawn"]
KH_dawn_18_minus = visit_18_minus["MINUS_KH_dawn"]
KH_dawn_19_minus = visit_19_minus["MINUS_KH_dawn"]
KH_dawn_20_minus = visit_20_minus["MINUS_KH_dawn"]
KH_dawn_21_minus = visit_21_minus["MINUS_KH_dawn"]
KH_dawn_24_minus = visit_24_minus["MINUS_KH_dawn"]
KH_dawn_25_minus = visit_25_minus["MINUS_KH_dawn"]
KH_dawn_27_minus = visit_27_minus["MINUS_KH_dawn"]
KH_dawn_28_minus = visit_28_minus["MINUS_KH_dawn"]
KH_dawn_34_minus = visit_34_minus["MINUS_KH_dawn"]
KH_dawn_35_minus = visit_35_minus["MINUS_KH_dawn"]

KH_dawn_01_minus = visit_01_minus["MINUS_KH_dawn"]
KH_dawn_12_minus = visit_12_minus["MINUS_KH_dawn"]
KH_dawn_15_minus = visit_15_minus["MINUS_KH_dawn"]
KH_dawn_26_minus = visit_26_minus["MINUS_KH_dawn"]

# ---------- Kelvin Helmhotlz - new, dusk ------------

KH_dusk_02 = visit_02_data["Kelvin_H_Dusk"]
KH_dusk_03 = visit_03_data["Kelvin_H_Dusk"]
KH_dusk_04 = visit_04_data["Kelvin_H_Dusk"]
KH_dusk_05 = visit_05_data["Kelvin_H_Dusk"]
KH_dusk_08 = visit_08_data["Kelvin_H_Dusk"]
KH_dusk_09 = visit_09_data["Kelvin_H_Dusk"]
KH_dusk_10 = visit_10_data["Kelvin_H_Dusk"]
KH_dusk_11 = visit_11_data["Kelvin_H_Dusk"]
KH_dusk_16 = visit_16_data["Kelvin_H_Dusk"]
KH_dusk_17 = visit_17_data["Kelvin_H_Dusk"]
KH_dusk_18 = visit_18_data["Kelvin_H_Dusk"]
KH_dusk_19 = visit_19_data["Kelvin_H_Dusk"]
KH_dusk_20 = visit_20_data["Kelvin_H_Dusk"]
KH_dusk_21 = visit_21_data["Kelvin_H_Dusk"]
KH_dusk_24 = visit_24_data["Kelvin_H_Dusk"]
KH_dusk_25 = visit_25_data["Kelvin_H_Dusk"]
KH_dusk_27 = visit_27_data["Kelvin_H_Dusk"]
KH_dusk_28 = visit_28_data["Kelvin_H_Dusk"]
KH_dusk_34 = visit_34_data["Kelvin_H_Dusk"]
KH_dusk_35 = visit_35_data["Kelvin_H_Dusk"]

KH_dusk_01 = visit_01_data["Kelvin_H_Dusk"]
KH_dusk_12 = visit_12_data["Kelvin_H_Dusk"]
KH_dusk_15 = visit_15_data["Kelvin_H_Dusk"]
KH_dusk_26 = visit_26_data["Kelvin_H_Dusk"]

# --------- plus % ---------------

KH_dusk_02_plus = visit_02_plus["PLUS_KH_dusk"]
KH_dusk_03_plus = visit_03_plus["PLUS_KH_dusk"]
KH_dusk_04_plus = visit_04_plus["PLUS_KH_dusk"]
KH_dusk_05_plus = visit_05_plus["PLUS_KH_dusk"]
KH_dusk_08_plus = visit_08_plus["PLUS_KH_dusk"]
KH_dusk_09_plus = visit_09_plus["PLUS_KH_dusk"]
KH_dusk_10_plus = visit_10_plus["PLUS_KH_dusk"]
KH_dusk_11_plus = visit_11_plus["PLUS_KH_dusk"]
KH_dusk_16_plus = visit_16_plus["PLUS_KH_dusk"]
KH_dusk_17_plus = visit_17_plus["PLUS_KH_dusk"]
KH_dusk_18_plus = visit_18_plus["PLUS_KH_dusk"]
KH_dusk_19_plus = visit_19_plus["PLUS_KH_dusk"]
KH_dusk_20_plus = visit_20_plus["PLUS_KH_dusk"]
KH_dusk_21_plus = visit_21_plus["PLUS_KH_dusk"]
KH_dusk_24_plus = visit_24_plus["PLUS_KH_dusk"]
KH_dusk_25_plus = visit_25_plus["PLUS_KH_dusk"]
KH_dusk_27_plus = visit_27_plus["PLUS_KH_dusk"]
KH_dusk_28_plus = visit_28_plus["PLUS_KH_dusk"]
KH_dusk_34_plus = visit_34_plus["PLUS_KH_dusk"]
KH_dusk_35_plus = visit_35_plus["PLUS_KH_dusk"]

KH_dusk_01_plus = visit_01_plus["PLUS_KH_dusk"]
KH_dusk_12_plus = visit_12_plus["PLUS_KH_dusk"]
KH_dusk_15_plus = visit_15_plus["PLUS_KH_dusk"]
KH_dusk_26_plus = visit_26_plus["PLUS_KH_dusk"]

# --------- minus % ---------------

KH_dusk_02_minus = visit_02_minus["MINUS_KH_dusk"]
KH_dusk_03_minus = visit_03_minus["MINUS_KH_dusk"]
KH_dusk_04_minus = visit_04_minus["MINUS_KH_dusk"]
KH_dusk_05_minus = visit_05_minus["MINUS_KH_dusk"]
KH_dusk_08_minus = visit_08_minus["MINUS_KH_dusk"]
KH_dusk_09_minus = visit_09_minus["MINUS_KH_dusk"]
KH_dusk_10_minus = visit_10_minus["MINUS_KH_dusk"]
KH_dusk_11_minus = visit_11_minus["MINUS_KH_dusk"]
KH_dusk_16_minus = visit_16_minus["MINUS_KH_dusk"]
KH_dusk_17_minus = visit_17_minus["MINUS_KH_dusk"]
KH_dusk_18_minus = visit_18_minus["MINUS_KH_dusk"]
KH_dusk_19_minus = visit_19_minus["MINUS_KH_dusk"]
KH_dusk_20_minus = visit_20_minus["MINUS_KH_dusk"]
KH_dusk_21_minus = visit_21_minus["MINUS_KH_dusk"]
KH_dusk_24_minus = visit_24_minus["MINUS_KH_dusk"]
KH_dusk_25_minus = visit_25_minus["MINUS_KH_dusk"]
KH_dusk_27_minus = visit_27_minus["MINUS_KH_dusk"]
KH_dusk_28_minus = visit_28_minus["MINUS_KH_dusk"]
KH_dusk_34_minus = visit_34_minus["MINUS_KH_dusk"]
KH_dusk_35_minus = visit_35_minus["MINUS_KH_dusk"]

KH_dusk_01_minus = visit_01_minus["MINUS_KH_dusk"]
KH_dusk_12_minus = visit_12_minus["MINUS_KH_dusk"]
KH_dusk_15_minus = visit_15_minus["MINUS_KH_dusk"]
KH_dusk_26_minus = visit_26_minus["MINUS_KH_dusk"]


# -------- powers -------

def polar_power(swirl, dusk, noon):
    total_power = []
    for i in range(len(swirl)):
        power = np.sum([swirl[i], dusk[i], noon[i]])
        total_power.append(power)
    return total_power

# error, standard deviation, ie
def standard_deviation(array, mean): # variable ie swirl_08, dusk_08 etc, mean as in swirl_09_mean
    var = []
    for g in range(len(array)):
        diff = (array[g] - mean)**2 
        var.append(diff)

    variance = (np.sum(var))/len(array)
    sd = np.sqrt(variance)
    return sd

polar_02 = np.array(polar_power(swirl_02, dusk_02, noon_02))
polar_03 = np.array(polar_power(swirl_03, dusk_03, noon_03))
polar_04 = np.array(polar_power(swirl_04, dusk_04, noon_04))
polar_05 = np.array(polar_power(swirl_05, dusk_05, noon_05))
polar_08 = np.array(polar_power(swirl_08, dusk_08, noon_08))
polar_09 = np.array(polar_power(swirl_09, dusk_09, noon_09))
polar_10 = np.array(polar_power(swirl_10, dusk_10, noon_10))
polar_11 = np.array(polar_power(swirl_11, dusk_11, noon_11))
polar_16 = np.array(polar_power(swirl_16, dusk_16, noon_16))
polar_17 = np.array(polar_power(swirl_17, dusk_17, noon_17))
polar_18 = np.array(polar_power(swirl_18, dusk_18, noon_18))
polar_19 = np.array(polar_power(swirl_19, dusk_19, noon_19))
polar_20 = np.array(polar_power(swirl_20, dusk_20, noon_20))
polar_21 = np.array(polar_power(swirl_21, dusk_21, noon_21))
polar_24 = np.array(polar_power(swirl_24, dusk_24, noon_24))
polar_25 = np.array(polar_power(swirl_25, dusk_25, noon_25))
polar_27 = np.array(polar_power(swirl_27, dusk_27, noon_27))
polar_28 = np.array(polar_power(swirl_28, dusk_28, noon_28))
polar_34 = np.array(polar_power(swirl_34, dusk_34, noon_34))
polar_35 = np.array(polar_power(swirl_35, dusk_35, noon_35))

# high cml
polar_01 = np.array(polar_power(swirl_01, dusk_01, noon_01))
polar_12 = np.array(polar_power(swirl_12, dusk_12, noon_12))
polar_15 = np.array(polar_power(swirl_15, dusk_15, noon_15))
polar_26 = np.array(polar_power(swirl_26, dusk_26, noon_26))


def error_calc(swirl_err, dusk_err):#, noon_err):
    total_err = []
    for i in range(len(swirl_err)):
        err = swirl_err[i]**2 + dusk_err[i]**2# + noon_err[i]**2
        error = np.sqrt(err)
        total_err.append(error)
    return total_err


polar_err02 = np.array(error_calc(swirl_err02, dusk_err02))#, noon_err02))
polar_err03 = np.array(error_calc(swirl_err03, dusk_err03))#, noon_err03))
polar_err04 = np.array(error_calc(swirl_err04, dusk_err04))#, noon_err04))
polar_err05 = np.array(error_calc(swirl_err05, dusk_err05))#, noon_err05))
polar_err08 = np.array(error_calc(swirl_err08, dusk_err08))#, noon_err08))
polar_err09 = np.array(error_calc(swirl_err09, dusk_err09))#, noon_err09))
polar_err10 = np.array(error_calc(swirl_err10, dusk_err10))#, noon_err10))
polar_err11 = np.array(error_calc(swirl_err11, dusk_err11))#, noon_err11))
polar_err16 = np.array(error_calc(swirl_err16, dusk_err16))#, noon_err16))
polar_err17 = np.array(error_calc(swirl_err17, dusk_err17))#, noon_err17))
polar_err18 = np.array(error_calc(swirl_err18, dusk_err18))##, noon_err18))
polar_err19 = np.array(error_calc(swirl_err19, dusk_err19))#, noon_err19))
polar_err20 = np.array(error_calc(swirl_err20, dusk_err20))#, noon_err20))
polar_err21 = np.array(error_calc(swirl_err21, dusk_err21))#, noon_err21))
polar_err24 = np.array(error_calc(swirl_err24, dusk_err24))#, noon_err24))
polar_err25 = np.array(error_calc(swirl_err25, dusk_err25))#, noon_err25))
polar_err27 = np.array(error_calc(swirl_err27, dusk_err27))#, noon_err27))
polar_err28 = np.array(error_calc(swirl_err28, dusk_err28))#, noon_err28))
polar_err34 = np.array(error_calc(swirl_err34, dusk_err34))#, noon_err34))
polar_err35 = np.array(error_calc(swirl_err35, dusk_err35))#, noon_err35))

# high cml
polar_err01 = np.array(error_calc(swirl_err01, dusk_err01))#, noon_err01))
polar_err12 = np.array(error_calc(swirl_err12, dusk_err12))#, noon_err12))
polar_err15 = np.array(error_calc(swirl_err15, dusk_err15))#, noon_err15))
polar_err26 = np.array(error_calc(swirl_err26, dusk_err26))#, noon_err26))


'''
AVERAGES - POLAR WHOLE - MAIN REGIONS
'''

polar_02_mean = np.mean(polar_02)
polar_03_mean = np.mean(polar_03)
polar_04_mean = np.mean(polar_04)
polar_05_mean = np.mean(polar_05)
polar_08_mean = np.mean(polar_08)
polar_09_mean = np.mean(polar_09)
polar_10_mean = np.mean(polar_10)
polar_11_mean = np.mean(polar_11)
polar_16_mean = np.mean(polar_16)
polar_17_mean = np.mean(polar_17)
polar_18_mean = np.mean(polar_18)
polar_19_mean = np.mean(polar_19)
polar_20_mean = np.mean(polar_20)
polar_21_mean = np.mean(polar_21)
polar_24_mean = np.mean(polar_24)
polar_25_mean = np.mean(polar_25)
polar_27_mean = np.mean(polar_27)
polar_28_mean = np.mean(polar_28)
polar_34_mean = np.mean(polar_34)
polar_35_mean = np.mean(polar_35)

polar_01_mean = np.mean(polar_01)
polar_12_mean = np.mean(polar_12)
polar_15_mean = np.mean(polar_15)
polar_26_mean = np.mean(polar_26)

#polar_means = [polar_02_mean, polar_03_mean, polar_04_mean, polar_05_mean, polar_08_mean, polar_09_mean, polar_10_mean, polar_11_mean, polar_16_mean, polar_17_mean, polar_18_mean, polar_19_mean, polar_20_mean, polar_21_mean, polar_24_mean, polar_25_mean, polar_27_mean, polar_28_mean, polar_34_mean, polar_35_mean]
polar_means = [polar_01_mean, polar_02_mean, polar_03_mean, polar_04_mean, polar_05_mean, polar_08_mean, polar_09_mean, polar_10_mean, polar_11_mean, polar_12_mean, polar_15_mean, polar_16_mean, polar_17_mean, polar_18_mean, polar_19_mean, polar_20_mean, polar_21_mean, polar_24_mean, polar_25_mean, polar_26_mean, polar_27_mean, polar_28_mean, polar_34_mean, polar_35_mean]
polar_means = np.array(polar_means)

sd_polar_02 = standard_deviation(polar_02, polar_02_mean)
sd_polar_03 = standard_deviation(polar_03, polar_03_mean)
sd_polar_04 = standard_deviation(polar_04, polar_04_mean)
sd_polar_05 = standard_deviation(polar_05, polar_05_mean)
sd_polar_08 = standard_deviation(polar_08, polar_08_mean)
sd_polar_09 = standard_deviation(polar_09, polar_09_mean)
sd_polar_10 = standard_deviation(polar_10, polar_10_mean)
sd_polar_11 = standard_deviation(polar_11, polar_11_mean)
sd_polar_16 = standard_deviation(polar_16, polar_16_mean)
sd_polar_17 = standard_deviation(polar_17, polar_17_mean)
sd_polar_18 = standard_deviation(polar_18, polar_18_mean)
sd_polar_19 = standard_deviation(polar_19, polar_19_mean)
sd_polar_20 = standard_deviation(polar_20, polar_20_mean)
sd_polar_21 = standard_deviation(polar_21, polar_21_mean)
sd_polar_24 = standard_deviation(polar_24, polar_24_mean)
sd_polar_25 = standard_deviation(polar_25, polar_25_mean)
sd_polar_27 = standard_deviation(polar_27, polar_27_mean)
sd_polar_28 = standard_deviation(polar_28, polar_28_mean)
sd_polar_34 = standard_deviation(polar_34, polar_34_mean)
sd_polar_35 = standard_deviation(polar_35, polar_35_mean)

sd_polar_01 = standard_deviation(polar_01, polar_01_mean)
sd_polar_12 = standard_deviation(polar_12, polar_12_mean)
sd_polar_15 = standard_deviation(polar_15, polar_15_mean)
sd_polar_26 = standard_deviation(polar_26, polar_26_mean)


polar_02_median = np.median(polar_02)
polar_03_median = np.median(polar_03)
polar_04_median = np.median(polar_04)
polar_05_median = np.median(polar_05)
polar_08_median = np.median(polar_08)
polar_09_median = np.median(polar_09)
polar_10_median = np.median(polar_10)
polar_11_median = np.median(polar_11)
polar_16_median = np.median(polar_16)
polar_17_median = np.median(polar_17)
polar_18_median = np.median(polar_18)
polar_19_median = np.median(polar_19)
polar_20_median = np.median(polar_20)
polar_21_median = np.median(polar_21)
polar_24_median = np.median(polar_24)
polar_25_median = np.median(polar_25)
polar_27_median = np.median(polar_27)
polar_28_median = np.median(polar_28)
polar_34_median = np.median(polar_34)
polar_35_median = np.median(polar_35)

polar_01_median = np.median(polar_01)
polar_12_median = np.median(polar_12)
polar_15_median = np.median(polar_15)
polar_26_median = np.median(polar_26)

#polar_medians = [polar_02_median, polar_03_median, polar_04_median, polar_05_median, polar_08_median, polar_09_median, polar_10_median, polar_11_median, polar_16_median, polar_17_median, polar_18_median, polar_19_median, polar_20_median, polar_21_median, polar_24_median, polar_25_median, polar_27_median, polar_28_median, polar_34_median, polar_35_median]
polar_medians = [polar_01_median, polar_02_median, polar_03_median, polar_04_median, polar_05_median, polar_08_median, polar_09_median, polar_10_median, polar_11_median, polar_12_median, polar_15_median, polar_16_median, polar_17_median, polar_18_median, polar_19_median, polar_20_median, polar_21_median, polar_24_median, polar_25_median, polar_26_median, polar_27_median, polar_28_median, polar_34_median, polar_35_median]
polar_medians = np.array(polar_medians)

med_polar_02_err = 1.2533 * sd_polar_02
med_polar_03_err = 1.2533 * sd_polar_03
med_polar_04_err = 1.2533 * sd_polar_04
med_polar_05_err = 1.2533 * sd_polar_05
med_polar_08_err = 1.2533 * sd_polar_08
med_polar_09_err = 1.2533 * sd_polar_09
med_polar_10_err = 1.2533 * sd_polar_10
med_polar_11_err = 1.2533 * sd_polar_11
med_polar_16_err = 1.2533 * sd_polar_16
med_polar_17_err = 1.2533 * sd_polar_17
med_polar_18_err = 1.2533 * sd_polar_18
med_polar_19_err = 1.2533 * sd_polar_19
med_polar_20_err = 1.2533 * sd_polar_20
med_polar_21_err = 1.2533 * sd_polar_21
med_polar_24_err = 1.2533 * sd_polar_24
med_polar_25_err = 1.2533 * sd_polar_25
med_polar_27_err = 1.2533 * sd_polar_27
med_polar_28_err = 1.2533 * sd_polar_28
med_polar_34_err = 1.2533 * sd_polar_34
med_polar_35_err = 1.2533 * sd_polar_35

med_polar_01_err = 1.2533 * sd_polar_01
med_polar_12_err = 1.2533 * sd_polar_12
med_polar_15_err = 1.2533 * sd_polar_15
med_polar_26_err = 1.2533 * sd_polar_26

#polar_med_errs = [med_polar_02_err, med_polar_03_err, med_polar_04_err, med_polar_05_err, med_polar_08_err, med_polar_09_err, med_polar_10_err, med_polar_11_err, med_polar_16_err, med_polar_17_err, med_polar_18_err, med_polar_19_err, med_polar_20_err, med_polar_21_err, med_polar_24_err, med_polar_25_err, med_polar_27_err, med_polar_28_err, med_polar_34_err, med_polar_35_err]
polar_med_errs = [med_polar_01_err, med_polar_02_err, med_polar_03_err, med_polar_04_err, med_polar_05_err, med_polar_08_err, med_polar_09_err, med_polar_10_err, med_polar_11_err, med_polar_12_err, med_polar_15_err, med_polar_16_err, med_polar_17_err, med_polar_18_err, med_polar_19_err, med_polar_20_err, med_polar_21_err, med_polar_24_err, med_polar_25_err, med_polar_26_err, med_polar_27_err, med_polar_28_err, med_polar_34_err, med_polar_35_err]
polar_med_errs = np.array(polar_med_errs)



'''
AVERAGES -  SWIRL
'''
# means - plot time as midpoint ie time_08[12]
swirl_02_mean = np.mean(swirl_02)
swirl_03_mean = np.mean(swirl_03)
swirl_04_mean = np.mean(swirl_04)
swirl_05_mean = np.mean(swirl_05)
swirl_08_mean = np.mean(swirl_08)
swirl_09_mean = np.mean(swirl_09)
swirl_10_mean = np.mean(swirl_10)
swirl_11_mean = np.mean(swirl_11)
swirl_16_mean = np.mean(swirl_16)
swirl_17_mean = np.mean(swirl_17)
swirl_18_mean = np.mean(swirl_18)
swirl_19_mean = np.mean(swirl_19)
swirl_20_mean = np.mean(swirl_20)
swirl_21_mean = np.mean(swirl_21)
swirl_24_mean = np.mean(swirl_24)
swirl_25_mean = np.mean(swirl_25)
swirl_27_mean = np.mean(swirl_27)
swirl_28_mean = np.mean(swirl_28)
swirl_34_mean = np.mean(swirl_34)
swirl_35_mean = np.mean(swirl_35)

# high cmls
swirl_01_mean = np.mean(swirl_01)
swirl_12_mean = np.mean(swirl_12)
swirl_15_mean = np.mean(swirl_15)
swirl_26_mean = np.mean(swirl_26)

#swirl_means = [swirl_02_mean, swirl_03_mean, swirl_04_mean, swirl_05_mean, swirl_08_mean, swirl_09_mean, swirl_10_mean, swirl_11_mean, swirl_16_mean, swirl_17_mean, swirl_18_mean, swirl_19_mean, swirl_20_mean, swirl_21_mean, swirl_24_mean, swirl_25_mean, swirl_27_mean, swirl_28_mean, swirl_34_mean, swirl_35_mean]
swirl_means = [swirl_01_mean, swirl_02_mean, swirl_03_mean, swirl_04_mean, swirl_05_mean, swirl_08_mean, swirl_09_mean, swirl_10_mean, swirl_11_mean, swirl_12_mean, swirl_15_mean, swirl_16_mean, swirl_17_mean, swirl_18_mean, swirl_19_mean, swirl_20_mean, swirl_21_mean, swirl_24_mean, swirl_25_mean, swirl_26_mean, swirl_27_mean, swirl_28_mean, swirl_34_mean, swirl_35_mean]
swirl_means = np.array(swirl_means)


sd_swirl_02 = standard_deviation(swirl_02, swirl_02_mean)
sd_swirl_03 = standard_deviation(swirl_03, swirl_03_mean)
sd_swirl_04 = standard_deviation(swirl_04, swirl_04_mean)
sd_swirl_05 = standard_deviation(swirl_05, swirl_05_mean)
sd_swirl_08 = standard_deviation(swirl_08, swirl_08_mean)
sd_swirl_09 = standard_deviation(swirl_09, swirl_09_mean)
sd_swirl_10 = standard_deviation(swirl_10, swirl_10_mean)
sd_swirl_11 = standard_deviation(swirl_11, swirl_11_mean)
sd_swirl_16 = standard_deviation(swirl_16, swirl_16_mean)
sd_swirl_17 = standard_deviation(swirl_17, swirl_17_mean)
sd_swirl_18 = standard_deviation(swirl_18, swirl_18_mean)
sd_swirl_19 = standard_deviation(swirl_19, swirl_19_mean)
sd_swirl_20 = standard_deviation(swirl_20, swirl_20_mean)
sd_swirl_21 = standard_deviation(swirl_21, swirl_21_mean)
sd_swirl_24 = standard_deviation(swirl_24, swirl_24_mean)
sd_swirl_27 = standard_deviation(swirl_27, swirl_27_mean)
sd_swirl_28 = standard_deviation(swirl_28, swirl_28_mean)
sd_swirl_25 = standard_deviation(swirl_25, swirl_25_mean)
sd_swirl_34 = standard_deviation(swirl_34, swirl_34_mean)
sd_swirl_35 = standard_deviation(swirl_35, swirl_35_mean)

# high cmls
sd_swirl_01 = standard_deviation(swirl_01, swirl_01_mean)
sd_swirl_12 = standard_deviation(swirl_12, swirl_12_mean)
sd_swirl_15 = standard_deviation(swirl_15, swirl_15_mean)
sd_swirl_26 = standard_deviation(swirl_26, swirl_26_mean)


# medians
swirl_02_median = np.median(swirl_02)
swirl_02_median = np.median(swirl_02)
swirl_03_median = np.median(swirl_03)
swirl_04_median = np.median(swirl_04)
swirl_05_median = np.median(swirl_05)
swirl_08_median = np.median(swirl_08)
swirl_09_median = np.median(swirl_09)
swirl_10_median = np.median(swirl_10)
swirl_11_median = np.median(swirl_11)
swirl_16_median = np.median(swirl_16)
swirl_17_median = np.median(swirl_17)
swirl_18_median = np.median(swirl_18)
swirl_19_median = np.median(swirl_19)
swirl_20_median = np.median(swirl_20)
swirl_21_median = np.median(swirl_21)
swirl_24_median = np.median(swirl_24)
swirl_25_median = np.median(swirl_25)
swirl_27_median = np.median(swirl_27)
swirl_28_median = np.median(swirl_28)
swirl_34_median = np.median(swirl_34)
swirl_35_median = np.median(swirl_35)

# high cmls
swirl_01_median = np.median(swirl_01)
swirl_12_median = np.median(swirl_12)
swirl_15_median = np.median(swirl_15)
swirl_26_median = np.median(swirl_26)

#medians_swirl = [swirl_02_median, swirl_03_median, swirl_04_median, swirl_05_median, swirl_08_median, swirl_09_median, swirl_10_median, swirl_11_median, swirl_16_median, swirl_17_median, swirl_18_median, swirl_19_median, swirl_20_median, swirl_21_median, swirl_24_median, swirl_25_median, swirl_27_median, swirl_28_median, swirl_34_median, swirl_35_median]
medians_swirl = [swirl_01_median, swirl_02_median, swirl_03_median, swirl_04_median, swirl_05_median, swirl_08_median, swirl_09_median, swirl_10_median, swirl_11_median, swirl_12_median, swirl_15_median, swirl_16_median, swirl_17_median, swirl_18_median, swirl_19_median, swirl_20_median, swirl_21_median, swirl_24_median, swirl_25_median, swirl_26_median, swirl_27_median, swirl_28_median, swirl_34_median, swirl_35_median]
medians_swirl = np.array(medians_swirl)


med_swirl_02_err = 1.2533 * sd_swirl_02
med_swirl_03_err = 1.2533 * sd_swirl_03
med_swirl_04_err = 1.2533 * sd_swirl_04
med_swirl_05_err = 1.2533 * sd_swirl_05
med_swirl_08_err = 1.2533 * sd_swirl_08
med_swirl_09_err = 1.2533 * sd_swirl_09
med_swirl_10_err = 1.2533 * sd_swirl_10
med_swirl_11_err = 1.2533 * sd_swirl_11
med_swirl_16_err = 1.2533 * sd_swirl_16
med_swirl_17_err = 1.2533 * sd_swirl_17
med_swirl_18_err = 1.2533 * sd_swirl_18
med_swirl_19_err = 1.2533 * sd_swirl_19
med_swirl_20_err = 1.2533 * sd_swirl_20
med_swirl_21_err = 1.2533 * sd_swirl_21
med_swirl_24_err = 1.2533 * sd_swirl_24
med_swirl_25_err = 1.2533 * sd_swirl_25
med_swirl_27_err = 1.2533 * sd_swirl_27
med_swirl_28_err = 1.2533 * sd_swirl_28
med_swirl_34_err = 1.2533 * sd_swirl_34
med_swirl_35_err = 1.2533 * sd_swirl_35

# high cmls
med_swirl_01_err = 1.2533 * sd_swirl_01
med_swirl_12_err = 1.2533 * sd_swirl_12
med_swirl_15_err = 1.2533 * sd_swirl_15
med_swirl_26_err = 1.2533 * sd_swirl_26

#swirl_med_errs = [med_swirl_02_err, med_swirl_03_err, med_swirl_04_err, med_swirl_05_err, med_swirl_08_err, med_swirl_09_err, med_swirl_10_err, med_swirl_11_err, med_swirl_16_err, med_swirl_17_err, med_swirl_18_err, med_swirl_19_err, med_swirl_20_err, med_swirl_21_err, med_swirl_24_err, med_swirl_25_err, med_swirl_27_err, med_swirl_28_err, med_swirl_34_err, med_swirl_35_err]
swirl_med_errs = [med_swirl_01_err, med_swirl_02_err, med_swirl_03_err, med_swirl_04_err, med_swirl_05_err, med_swirl_08_err, med_swirl_09_err, med_swirl_10_err, med_swirl_11_err, med_swirl_12_err, med_swirl_15_err, med_swirl_16_err, med_swirl_17_err, med_swirl_18_err, med_swirl_19_err, med_swirl_20_err, med_swirl_21_err, med_swirl_24_err, med_swirl_25_err, med_swirl_26_err, med_swirl_27_err, med_swirl_28_err, med_swirl_34_err, med_swirl_35_err]
swirl_med_errs = np.array(swirl_med_errs)



'''
AVERAGES - DUSK

'''

dusk_02_mean = np.mean(dusk_02)
dusk_03_mean = np.mean(dusk_03)
dusk_04_mean = np.mean(dusk_04)
dusk_05_mean = np.mean(dusk_05)
dusk_08_mean = np.mean(dusk_08)
dusk_09_mean = np.mean(dusk_09)
dusk_10_mean = np.mean(dusk_10)
dusk_11_mean = np.mean(dusk_11)
dusk_16_mean = np.mean(dusk_16)
dusk_17_mean = np.mean(dusk_17)
dusk_18_mean = np.mean(dusk_18)
dusk_19_mean = np.mean(dusk_19)
dusk_20_mean = np.mean(dusk_20)
dusk_21_mean = np.mean(dusk_21)
dusk_24_mean = np.mean(dusk_24)
dusk_25_mean = np.mean(dusk_25)
dusk_27_mean = np.mean(dusk_27)
dusk_28_mean = np.mean(dusk_28)
dusk_34_mean = np.mean(dusk_34)
dusk_35_mean = np.mean(dusk_35)

#high cmls
dusk_01_mean = np.mean(dusk_11)
dusk_12_mean = np.mean(dusk_12)
dusk_15_mean = np.mean(dusk_15)
dusk_26_mean = np.mean(dusk_26)

#dusk_means = [dusk_02_mean, dusk_03_mean, dusk_04_mean, dusk_05_mean, dusk_08_mean, dusk_09_mean, dusk_10_mean, dusk_11_mean, dusk_16_mean, dusk_17_mean, dusk_18_mean, dusk_19_mean, dusk_20_mean, dusk_21_mean, dusk_24_mean, dusk_25_mean, dusk_27_mean, dusk_28_mean, dusk_34_mean, dusk_35_mean]
dusk_means = [dusk_01_mean, dusk_02_mean, dusk_03_mean, dusk_04_mean, dusk_05_mean, dusk_08_mean, dusk_09_mean, dusk_10_mean, dusk_11_mean, dusk_12_mean, dusk_15_mean, dusk_16_mean, dusk_17_mean, dusk_18_mean, dusk_19_mean, dusk_20_mean, dusk_21_mean, dusk_24_mean, dusk_25_mean, dusk_26_mean, dusk_27_mean, dusk_28_mean, dusk_34_mean, dusk_35_mean]
dusk_means = np.array(dusk_means)


sd_dusk_02 = standard_deviation(dusk_02, dusk_02_mean)
sd_dusk_03 = standard_deviation(dusk_03, dusk_03_mean)
sd_dusk_04 = standard_deviation(dusk_04, dusk_04_mean)
sd_dusk_05 = standard_deviation(dusk_05, dusk_05_mean)
sd_dusk_08 = standard_deviation(dusk_08, dusk_08_mean)
sd_dusk_09 = standard_deviation(dusk_09, dusk_09_mean)
sd_dusk_10 = standard_deviation(dusk_10, dusk_10_mean)
sd_dusk_11 = standard_deviation(dusk_11, dusk_11_mean)
sd_dusk_16 = standard_deviation(dusk_16, dusk_16_mean)
sd_dusk_17 = standard_deviation(dusk_17, dusk_17_mean)
sd_dusk_18 = standard_deviation(dusk_18, dusk_18_mean)
sd_dusk_19 = standard_deviation(dusk_19, dusk_19_mean)
sd_dusk_20 = standard_deviation(dusk_20, dusk_20_mean)
sd_dusk_21 = standard_deviation(dusk_21, dusk_21_mean)
sd_dusk_24 = standard_deviation(dusk_24, dusk_24_mean)
sd_dusk_25 = standard_deviation(dusk_25, dusk_25_mean)
sd_dusk_27 = standard_deviation(dusk_27, dusk_27_mean)
sd_dusk_28 = standard_deviation(dusk_28, dusk_28_mean)
sd_dusk_34 = standard_deviation(dusk_34, dusk_34_mean)
sd_dusk_35 = standard_deviation(dusk_35, dusk_35_mean)

# high cmls
sd_dusk_01 = standard_deviation(dusk_01, dusk_01_mean)
sd_dusk_12 = standard_deviation(dusk_12, dusk_12_mean)
sd_dusk_15 = standard_deviation(dusk_15, dusk_15_mean)
sd_dusk_26 = standard_deviation(dusk_26, dusk_26_mean)


dusk_02_median = np.median(dusk_02)
dusk_03_median = np.median(dusk_03)
dusk_04_median = np.median(dusk_04)
dusk_05_median = np.median(dusk_05)
dusk_08_median = np.median(dusk_08)
dusk_09_median = np.median(dusk_09)
dusk_10_median = np.median(dusk_10)
dusk_11_median = np.median(dusk_11)
dusk_16_median = np.median(dusk_16)
dusk_17_median = np.median(dusk_17)
dusk_18_median = np.median(dusk_18)
dusk_19_median = np.median(dusk_19)
dusk_20_median = np.median(dusk_20)
dusk_21_median = np.median(dusk_21)
dusk_24_median = np.median(dusk_24)
dusk_25_median = np.median(dusk_25)
dusk_27_median = np.median(dusk_27)
dusk_28_median = np.median(dusk_28)
dusk_34_median = np.median(dusk_34)
dusk_35_median = np.median(dusk_35)

# high cmls
dusk_01_median = np.median(dusk_01)
dusk_12_median = np.median(dusk_12)
dusk_15_median = np.median(dusk_15)
dusk_26_median = np.median(dusk_26)

#dusk_medians = [dusk_02_median, dusk_03_median, dusk_04_median, dusk_05_median, dusk_08_median, dusk_09_median, dusk_10_median, dusk_11_median, dusk_16_median, dusk_17_median, dusk_18_median, dusk_19_median, dusk_20_median, dusk_21_median, dusk_24_median, dusk_25_median, dusk_27_median, dusk_28_median, dusk_34_median, dusk_35_median]
dusk_medians = [dusk_01_median, dusk_02_median, dusk_03_median, dusk_04_median, dusk_05_median, dusk_08_median, dusk_09_median, dusk_10_median, dusk_11_median, dusk_12_median, dusk_15_median, dusk_16_median, dusk_17_median, dusk_18_median, dusk_19_median, dusk_20_median, dusk_21_median, dusk_24_median, dusk_25_median, dusk_26_median, dusk_27_median, dusk_28_median, dusk_34_median, dusk_35_median]
dusk_medians = np.array(dusk_medians)


med_dusk_02_err = 1.2533 * sd_dusk_02
med_dusk_03_err = 1.2533 * sd_dusk_03
med_dusk_04_err = 1.2533 * sd_dusk_04
med_dusk_05_err = 1.2533 * sd_dusk_05
med_dusk_08_err = 1.2533 * sd_dusk_08
med_dusk_09_err = 1.2533 * sd_dusk_09
med_dusk_10_err = 1.2533 * sd_dusk_10
med_dusk_11_err = 1.2533 * sd_dusk_11
med_dusk_16_err = 1.2533 * sd_dusk_16
med_dusk_17_err = 1.2533 * sd_dusk_17
med_dusk_18_err = 1.2533 * sd_dusk_18
med_dusk_19_err = 1.2533 * sd_dusk_19
med_dusk_20_err = 1.2533 * sd_dusk_20
med_dusk_21_err = 1.2533 * sd_dusk_21
med_dusk_24_err = 1.2533 * sd_dusk_24
med_dusk_25_err = 1.2533 * sd_dusk_25
med_dusk_27_err = 1.2533 * sd_dusk_27
med_dusk_28_err = 1.2533 * sd_dusk_28
med_dusk_34_err = 1.2533 * sd_dusk_34
med_dusk_35_err = 1.2533 * sd_dusk_35

# high cmls
med_dusk_01_err = 1.2533 * sd_dusk_01
med_dusk_12_err = 1.2533 * sd_dusk_12
med_dusk_15_err = 1.2533 * sd_dusk_15
med_dusk_26_err = 1.2533 * sd_dusk_26

#dusk_med_errs = [med_dusk_02_err, med_dusk_03_err, med_dusk_04_err, med_dusk_05_err, med_dusk_08_err, med_dusk_09_err, med_dusk_10_err, med_dusk_11_err, med_dusk_16_err, med_dusk_17_err, med_dusk_18_err, med_dusk_19_err, med_dusk_20_err, med_dusk_21_err, med_dusk_24_err, med_dusk_25_err, med_dusk_27_err, med_dusk_28_err, med_dusk_34_err, med_dusk_35_err]
dusk_med_errs = [med_dusk_01_err, med_dusk_02_err, med_dusk_03_err, med_dusk_04_err, med_dusk_05_err, med_dusk_08_err, med_dusk_09_err, med_dusk_10_err, med_dusk_11_err, med_dusk_12_err, med_dusk_15_err, med_dusk_16_err, med_dusk_17_err, med_dusk_18_err, med_dusk_19_err, med_dusk_20_err, med_dusk_21_err, med_dusk_24_err, med_dusk_25_err, med_dusk_26_err, med_dusk_27_err, med_dusk_28_err, med_dusk_34_err, med_dusk_35_err]
dusk_med_errs = np.array(dusk_med_errs)



'''
AVERAGES - NOON
'''

noon_02_mean = np.mean(noon_02)
noon_03_mean = np.mean(noon_03)
noon_04_mean = np.mean(noon_04)
noon_05_mean = np.mean(noon_05)
noon_08_mean = np.mean(noon_08)
noon_09_mean = np.mean(noon_09)
noon_10_mean = np.mean(noon_10)
noon_11_mean = np.mean(noon_11)
noon_16_mean = np.mean(noon_16)
noon_17_mean = np.mean(noon_17)
noon_18_mean = np.mean(noon_18)
noon_19_mean = np.mean(noon_19)
noon_20_mean = np.mean(noon_20)
noon_21_mean = np.mean(noon_21)
noon_24_mean = np.mean(noon_24)
noon_25_mean = np.mean(noon_25)
noon_27_mean = np.mean(noon_27)
noon_28_mean = np.mean(noon_28)
noon_34_mean = np.mean(noon_34)
noon_35_mean = np.mean(noon_35)

# high cmls
noon_01_mean = np.mean(noon_01)
noon_12_mean = np.mean(noon_12)
noon_15_mean = np.mean(noon_15)
noon_26_mean = np.mean(noon_26)

#noon_means = [noon_02_mean, noon_03_mean, noon_04_mean, noon_05_mean, noon_08_mean, noon_09_mean, noon_10_mean, noon_11_mean, noon_16_mean, noon_17_mean, noon_18_mean, noon_19_mean, noon_20_mean, noon_21_mean, noon_24_mean, noon_25_mean, noon_27_mean, noon_28_mean, noon_34_mean, noon_35_mean]
noon_means = [noon_01_mean, noon_02_mean, noon_03_mean, noon_04_mean, noon_05_mean, noon_08_mean, noon_09_mean, noon_10_mean, noon_11_mean, noon_12_mean, noon_15_mean, noon_16_mean, noon_17_mean, noon_18_mean, noon_19_mean, noon_20_mean, noon_21_mean, noon_24_mean, noon_25_mean, noon_26_mean, noon_27_mean, noon_28_mean, noon_34_mean, noon_35_mean]
noon_means = np.array(noon_means)


sd_noon_02 = standard_deviation(noon_02, noon_02_mean)
sd_noon_03 = standard_deviation(noon_03, noon_03_mean)
sd_noon_04 = standard_deviation(noon_04, noon_04_mean)
sd_noon_05 = standard_deviation(noon_05, noon_05_mean)
sd_noon_08 = standard_deviation(noon_08, noon_08_mean)
sd_noon_09 = standard_deviation(noon_09, noon_09_mean)
sd_noon_10 = standard_deviation(noon_10, noon_10_mean)
sd_noon_11 = standard_deviation(noon_11, noon_11_mean)
sd_noon_16 = standard_deviation(noon_16, noon_16_mean)
sd_noon_17 = standard_deviation(noon_17, noon_17_mean)
sd_noon_18 = standard_deviation(noon_18, noon_18_mean)
sd_noon_19 = standard_deviation(noon_19, noon_19_mean)
sd_noon_20 = standard_deviation(noon_20, noon_20_mean)
sd_noon_21 = standard_deviation(noon_21, noon_21_mean)
sd_noon_24 = standard_deviation(noon_24, noon_24_mean)
sd_noon_25 = standard_deviation(noon_25, noon_25_mean)
sd_noon_27 = standard_deviation(noon_27, noon_27_mean)
sd_noon_28 = standard_deviation(noon_28, noon_28_mean)
sd_noon_34 = standard_deviation(noon_34, noon_34_mean)
sd_noon_35 = standard_deviation(noon_35, noon_35_mean)

# high cmls
sd_noon_01 = standard_deviation(noon_01, noon_01_mean)
sd_noon_12 = standard_deviation(noon_12, noon_12_mean)
sd_noon_15 = standard_deviation(noon_15, noon_15_mean)
sd_noon_26 = standard_deviation(noon_26, noon_26_mean)


noon_02_median = np.median(noon_02)
noon_03_median = np.median(noon_03)
noon_04_median = np.median(noon_04)
noon_05_median = np.median(noon_05)
noon_08_median = np.median(noon_08)
noon_09_median = np.median(noon_09)
noon_10_median = np.median(noon_10)
noon_11_median = np.median(noon_11)
noon_16_median = np.median(noon_16)
noon_17_median = np.median(noon_17)
noon_18_median = np.median(noon_18)
noon_19_median = np.median(noon_19)
noon_20_median = np.median(noon_20)
noon_21_median = np.median(noon_21)
noon_24_median = np.median(noon_24)
noon_25_median = np.median(noon_25)
noon_27_median = np.median(noon_27)
noon_28_median = np.median(noon_28)
noon_34_median = np.median(noon_34)
noon_35_median = np.median(noon_35)

# high cmls
noon_01_median = np.median(noon_01)
noon_12_median = np.median(noon_12)
noon_15_median = np.median(noon_15)
noon_26_median = np.median(noon_26)

#noon_medians = [noon_02_mean, noon_03_mean, noon_04_median, noon_05_median, noon_08_median, noon_09_median, noon_10_median, noon_11_median, noon_16_median, noon_17_median, noon_18_median, noon_19_median, noon_20_median, noon_21_median, noon_24_median, noon_25_median, noon_27_median, noon_28_median, noon_34_median, noon_35_median]
noon_medians = [noon_01_mean, noon_02_mean, noon_03_mean, noon_04_median, noon_05_median, noon_08_median, noon_09_median, noon_10_median, noon_11_median, noon_12_median, noon_15_median, noon_16_median, noon_17_median, noon_18_median, noon_19_median, noon_20_median, noon_21_median, noon_24_median, noon_25_median, noon_26_median, noon_27_median, noon_28_median, noon_34_median, noon_35_median]
noon_medians = np.array(noon_medians)


med_noon_02_err = 1.2533 * sd_noon_02
med_noon_03_err = 1.2533 * sd_noon_03
med_noon_04_err = 1.2533 * sd_noon_04
med_noon_05_err = 1.2533 * sd_noon_05
med_noon_08_err = 1.2533 * sd_noon_08
med_noon_09_err = 1.2533 * sd_noon_09
med_noon_10_err = 1.2533 * sd_noon_10
med_noon_11_err = 1.2533 * sd_noon_11
med_noon_16_err = 1.2533 * sd_noon_16
med_noon_17_err = 1.2533 * sd_noon_17
med_noon_18_err = 1.2533 * sd_noon_18
med_noon_19_err = 1.2533 * sd_noon_19
med_noon_20_err = 1.2533 * sd_noon_20
med_noon_21_err = 1.2533 * sd_noon_21
med_noon_24_err = 1.2533 * sd_noon_24
med_noon_25_err = 1.2533 * sd_noon_25
med_noon_27_err = 1.2533 * sd_noon_27
med_noon_28_err = 1.2533 * sd_noon_28
med_noon_34_err = 1.2533 * sd_noon_34
med_noon_35_err = 1.2533 * sd_noon_35

# high cmls
med_noon_01_err = 1.2533 * sd_noon_01
med_noon_12_err = 1.2533 * sd_noon_12
med_noon_15_err = 1.2533 * sd_noon_15
med_noon_26_err = 1.2533 * sd_noon_26

#noon_med_errs = [med_noon_02_err, med_noon_03_err, med_noon_04_err, med_noon_05_err, med_noon_08_err, med_noon_09_err, med_noon_10_err, med_noon_11_err, med_noon_16_err, med_noon_17_err, med_noon_18_err, med_noon_19_err, med_noon_20_err, med_noon_21_err, med_noon_24_err, med_noon_25_err, med_noon_27_err, med_noon_28_err, med_noon_34_err, med_noon_35_err]
noon_med_errs = [med_noon_01_err, med_noon_02_err, med_noon_03_err, med_noon_04_err, med_noon_05_err, med_noon_08_err, med_noon_09_err, med_noon_10_err, med_noon_11_err, med_noon_12_err, med_noon_15_err, med_noon_16_err, med_noon_17_err, med_noon_18_err, med_noon_19_err, med_noon_20_err, med_noon_21_err, med_noon_24_err, med_noon_25_err, med_noon_26_err, med_noon_27_err, med_noon_28_err, med_noon_34_err, med_noon_35_err]
noon_med_errs = np.array(noon_med_errs)


'''
AVERAGE PRESSURE
'''

pressure_02_mean = np.mean(pressure_02)
pressure_03_mean = np.mean(pressure_03)
pressure_04_mean = np.mean(pressure_04)
pressure_05_mean = np.mean(pressure_05)
pressure_08_mean = np.mean(pressure_08)
pressure_09_mean = np.mean(pressure_09)
pressure_10_mean = np.mean(pressure_10)
pressure_11_mean = np.mean(pressure_11)
pressure_16_mean = np.mean(pressure_16)
pressure_17_mean = np.mean(pressure_17)
pressure_18_mean = np.mean(pressure_18)
pressure_18_mean = np.mean(pressure_18)
pressure_19_mean = np.mean(pressure_19)
pressure_20_mean = np.mean(pressure_20)
pressure_21_mean = np.mean(pressure_21)
pressure_24_mean = np.mean(pressure_24)
pressure_25_mean = np.mean(pressure_25)
pressure_27_mean = np.mean(pressure_27)
pressure_28_mean = np.mean(pressure_28)
pressure_34_mean = np.mean(pressure_34)
pressure_35_mean = np.mean(pressure_35)

# #pressure_13_mean = np.mean(pressure_13)

# high cmls
pressure_01_mean = np.mean(pressure_01)
pressure_12_mean = np.mean(pressure_12)
pressure_15_mean = np.mean(pressure_15)
pressure_26_mean = np.mean(pressure_26)

#pressure_means = [pressure_02_mean, pressure_03_mean, pressure_04_mean, pressure_05_mean, pressure_08_mean, pressure_09_mean, pressure_10_mean, pressure_11_mean, pressure_16_mean, pressure_17_mean, pressure_18_mean, pressure_19_mean, pressure_20_mean, pressure_21_mean, pressure_24_mean, pressure_25_mean, pressure_27_mean, pressure_28_mean, pressure_34_mean, pressure_35_mean]
pressure_means = [pressure_01_mean, pressure_02_mean, pressure_03_mean, pressure_04_mean, pressure_05_mean, pressure_08_mean, pressure_09_mean, pressure_10_mean, pressure_11_mean, pressure_12_mean, pressure_15_mean, pressure_16_mean, pressure_17_mean, pressure_18_mean, pressure_19_mean, pressure_20_mean, pressure_21_mean, pressure_24_mean, pressure_25_mean, pressure_26_mean, pressure_27_mean, pressure_28_mean, pressure_34_mean, pressure_35_mean]
pressure_means = np.array(pressure_means)


sd_pressure_02 = standard_deviation(pressure_02, pressure_02_mean)
sd_pressure_03 = standard_deviation(pressure_03, pressure_03_mean)
sd_pressure_04 = standard_deviation(pressure_04, pressure_04_mean)
sd_pressure_05 = standard_deviation(pressure_05, pressure_05_mean)
sd_pressure_08 = standard_deviation(pressure_08, pressure_08_mean)
sd_pressure_09 = standard_deviation(pressure_09, pressure_09_mean)
sd_pressure_10 = standard_deviation(pressure_10, pressure_10_mean)
sd_pressure_11 = standard_deviation(pressure_11, pressure_11_mean)
sd_pressure_16 = standard_deviation(pressure_16, pressure_16_mean)
sd_pressure_17 = standard_deviation(pressure_17, pressure_17_mean)
sd_pressure_18 = standard_deviation(pressure_18, pressure_18_mean)
sd_pressure_19 = standard_deviation(pressure_19, pressure_19_mean)
sd_pressure_20 = standard_deviation(pressure_20, pressure_20_mean)
sd_pressure_21 = standard_deviation(pressure_21, pressure_21_mean)
sd_pressure_24 = standard_deviation(pressure_24, pressure_24_mean)
sd_pressure_25 = standard_deviation(pressure_25, pressure_25_mean)
sd_pressure_27 = standard_deviation(pressure_27, pressure_27_mean)
sd_pressure_28 = standard_deviation(pressure_28, pressure_28_mean)
sd_pressure_34 = standard_deviation(pressure_34, pressure_34_mean)
sd_pressure_35 = standard_deviation(pressure_35, pressure_35_mean)

# high cmls
sd_pressure_01 = standard_deviation(pressure_01, pressure_01_mean)
sd_pressure_12 = standard_deviation(pressure_12, pressure_12_mean)
sd_pressure_15 = standard_deviation(pressure_15, pressure_15_mean)
sd_pressure_26 = standard_deviation(pressure_26, pressure_26_mean)


pressure_02_median = np.median(pressure_02)
pressure_03_median = np.median(pressure_03)
pressure_04_median = np.median(pressure_04)
pressure_05_median = np.median(pressure_05)
pressure_08_median = np.median(pressure_08)
pressure_09_median = np.median(pressure_09)
pressure_10_median = np.median(pressure_10)
pressure_11_median = np.median(pressure_11)
pressure_16_median = np.median(pressure_16)
pressure_17_median = np.median(pressure_17)
pressure_18_median = np.median(pressure_18)
pressure_19_median = np.median(pressure_19)
pressure_20_median = np.median(pressure_20)
pressure_21_median = np.median(pressure_21)
pressure_24_median = np.median(pressure_24)
pressure_25_median = np.median(pressure_25)
pressure_27_median = np.median(pressure_27)
pressure_28_median = np.median(pressure_28)
pressure_34_median = np.median(pressure_34)
pressure_35_median = np.median(pressure_35)

# high cmls
pressure_01_median = np.median(pressure_01)
pressure_12_median = np.median(pressure_12)
pressure_15_median = np.median(pressure_15)
pressure_26_median = np.median(pressure_26)

pressure_medians = [pressure_01_median, pressure_02_median, pressure_03_median, pressure_04_median, pressure_05_median, pressure_08_median, pressure_09_median, pressure_10_median, pressure_11_median, pressure_12_median, pressure_15_median, pressure_16_median, pressure_17_median, pressure_18_median, pressure_19_median, pressure_20_median, pressure_21_median, pressure_24_median, pressure_25_median, pressure_26_median, pressure_27_median, pressure_28_median, pressure_34_median, pressure_35_median]
pressure_medians = np.array(pressure_medians)


med_pressure_02_err = 1.2533 * sd_pressure_02
med_pressure_03_err = 1.2533 * sd_pressure_03
med_pressure_04_err = 1.2533 * sd_pressure_04
med_pressure_05_err = 1.2533 * sd_pressure_05
med_pressure_08_err = 1.2533 * sd_pressure_08
med_pressure_09_err = 1.2533 * sd_pressure_09
med_pressure_10_err = 1.2533 * sd_pressure_10
med_pressure_11_err = 1.2533 * sd_pressure_11
med_pressure_16_err = 1.2533 * sd_pressure_16
med_pressure_17_err = 1.2533 * sd_pressure_17
med_pressure_18_err = 1.2533 * sd_pressure_18
med_pressure_19_err = 1.2533 * sd_pressure_19
med_pressure_20_err = 1.2533 * sd_pressure_20
med_pressure_21_err = 1.2533 * sd_pressure_21
med_pressure_24_err = 1.2533 * sd_pressure_24
med_pressure_25_err = 1.2533 * sd_pressure_25
med_pressure_27_err = 1.2533 * sd_pressure_27
med_pressure_28_err = 1.2533 * sd_pressure_28
med_pressure_34_err = 1.2533 * sd_pressure_34
med_pressure_35_err = 1.2533 * sd_pressure_35

# high cmls
med_pressure_01_err = 1.2533 * sd_pressure_01
med_pressure_12_err = 1.2533 * sd_pressure_12
med_pressure_15_err = 1.2533 * sd_pressure_15
med_pressure_26_err = 1.2533 * sd_pressure_26

#med_pressure_errs = [med_pressure_02_err, med_pressure_03_err, med_pressure_04_err, med_pressure_05_err, med_pressure_08_err, med_pressure_09_err, med_pressure_10_err, med_pressure_11_err, med_pressure_16_err, med_pressure_17_err, med_pressure_18_err, med_pressure_19_err, med_pressure_20_err, med_pressure_21_err, med_pressure_24_err, med_pressure_25_err, med_pressure_27_err, med_pressure_28_err, med_pressure_34_err, med_pressure_35_err]
med_pressure_errs = [med_pressure_01_err, med_pressure_02_err, med_pressure_03_err, med_pressure_04_err, med_pressure_05_err, med_pressure_08_err, med_pressure_09_err, med_pressure_10_err, med_pressure_11_err, med_pressure_12_err, med_pressure_15_err, med_pressure_16_err, med_pressure_17_err, med_pressure_18_err, med_pressure_19_err, med_pressure_20_err, med_pressure_21_err, med_pressure_24_err, med_pressure_25_err, med_pressure_26_err, med_pressure_27_err, med_pressure_28_err, med_pressure_34_err, med_pressure_35_err]
med_pressure_errs = np.array(med_pressure_errs)



'''
#travel time +/- 10% pressure
'''

pressure_02_mean_p = np.mean(pressure_02_plus)
pressure_03_mean_p = np.mean(pressure_03_plus)
pressure_04_mean_p = np.mean(pressure_04_plus)
pressure_05_mean_p = np.mean(pressure_05_plus)
pressure_08_mean_p = np.mean(pressure_08_plus)
pressure_09_mean_p = np.mean(pressure_09_plus)
pressure_10_mean_p = np.mean(pressure_10_plus)
pressure_11_mean_p = np.mean(pressure_11_plus)
pressure_16_mean_p = np.mean(pressure_16_plus)
pressure_17_mean_p = np.mean(pressure_17_plus)
pressure_18_mean_p = np.mean(pressure_18_plus)
pressure_19_mean_p = np.mean(pressure_19_plus)
pressure_20_mean_p = np.mean(pressure_20_plus)
pressure_21_mean_p = np.mean(pressure_21_plus)
pressure_24_mean_p = np.mean(pressure_24_plus)
pressure_25_mean_p = np.mean(pressure_25_plus)
pressure_27_mean_p = np.mean(pressure_27_plus)
pressure_28_mean_p = np.mean(pressure_28_plus)
pressure_34_mean_p = np.mean(pressure_34_plus)
pressure_35_mean_p = np.mean(pressure_35_plus)

# high cmls
pressure_01_mean_p = np.mean(pressure_01_plus)
pressure_12_mean_p = np.mean(pressure_12_plus)
pressure_15_mean_p = np.mean(pressure_15_plus)
pressure_26_mean_p = np.mean(pressure_26_plus)

#pressure_means_p = np.array([pressure_02_mean_p,pressure_03_mean_p,pressure_04_mean_p,pressure_05_mean_p,pressure_08_mean_p,pressure_09_mean_p,pressure_10_mean_p,pressure_11_mean_p,pressure_16_mean_p,pressure_17_mean_p,pressure_18_mean_p,pressure_19_mean_p,pressure_20_mean_p,pressure_21_mean_p,pressure_24_mean_p,pressure_25_mean_p,pressure_27_mean_p,pressure_28_mean_p,pressure_34_mean_p,pressure_35_mean_p])
pressure_means_p = np.array([pressure_01_mean_p,pressure_02_mean_p,pressure_03_mean_p,pressure_04_mean_p,pressure_05_mean_p,pressure_08_mean_p,pressure_09_mean_p,pressure_10_mean_p,pressure_11_mean_p,pressure_12_mean_p,pressure_15_mean_p,pressure_16_mean_p,pressure_17_mean_p,pressure_18_mean_p,pressure_19_mean_p,pressure_20_mean_p,pressure_21_mean_p,pressure_24_mean_p,pressure_25_mean_p,pressure_26_mean_p,pressure_27_mean_p,pressure_28_mean_p,pressure_34_mean_p,pressure_35_mean_p])


pressure_02_mean_m = np.mean(pressure_02_minus)
pressure_03_mean_m = np.mean(pressure_03_minus)
pressure_04_mean_m = np.mean(pressure_04_minus)
pressure_05_mean_m = np.mean(pressure_05_minus)
pressure_08_mean_m = np.mean(pressure_08_minus)
pressure_09_mean_m = np.mean(pressure_09_minus)
pressure_10_mean_m = np.mean(pressure_10_minus)
pressure_11_mean_m = np.mean(pressure_11_minus)
pressure_16_mean_m = np.mean(pressure_16_minus)
pressure_17_mean_m = np.mean(pressure_17_minus)
pressure_18_mean_m = np.mean(pressure_18_minus)
pressure_19_mean_m = np.mean(pressure_19_minus)
pressure_20_mean_m = np.mean(pressure_20_minus)
pressure_21_mean_m = np.mean(pressure_21_minus)
pressure_24_mean_m = np.mean(pressure_24_minus)
pressure_25_mean_m = np.mean(pressure_25_minus)
pressure_27_mean_m = np.mean(pressure_27_minus)
pressure_28_mean_m = np.mean(pressure_28_minus)
pressure_34_mean_m = np.mean(pressure_34_minus)
pressure_35_mean_m = np.mean(pressure_35_minus)

# high cmls
pressure_01_mean_m = np.mean(pressure_01_minus)
pressure_12_mean_m = np.mean(pressure_12_minus)
pressure_15_mean_m = np.mean(pressure_15_minus)
pressure_26_mean_m = np.mean(pressure_26_minus)

#pressure_means_m = np.array([pressure_02_mean_m,pressure_03_mean_m,pressure_04_mean_m,pressure_05_mean_m,pressure_08_mean_m,pressure_09_mean_m,pressure_10_mean_m,pressure_11_mean_m,pressure_16_mean_m,pressure_17_mean_m,pressure_18_mean_m,pressure_19_mean_m,pressure_20_mean_m,pressure_21_mean_m,pressure_24_mean_m,pressure_25_mean_m,pressure_27_mean_m,pressure_28_mean_m,pressure_34_mean_m,pressure_35_mean_m])
pressure_means_m = np.array([pressure_01_mean_m,pressure_02_mean_m,pressure_03_mean_m,pressure_04_mean_m,pressure_05_mean_m,pressure_08_mean_m,pressure_09_mean_m,pressure_10_mean_m,pressure_11_mean_m,pressure_12_mean_m,pressure_15_mean_m,pressure_16_mean_m,pressure_17_mean_m,pressure_18_mean_m,pressure_19_mean_m,pressure_20_mean_m,pressure_21_mean_m,pressure_24_mean_m,pressure_25_mean_m,pressure_26_mean_m,pressure_27_mean_m,pressure_28_mean_m,pressure_34_mean_m,pressure_35_mean_m])


sd_pressure_02_p = standard_deviation(pressure_02_plus, pressure_02_mean_p)
sd_pressure_03_p = standard_deviation(pressure_03_plus, pressure_03_mean_p)
sd_pressure_04_p = standard_deviation(pressure_04_plus, pressure_04_mean_p)
sd_pressure_05_p = standard_deviation(pressure_05_plus, pressure_05_mean_p)
sd_pressure_08_p = standard_deviation(pressure_08_plus, pressure_08_mean_p)
sd_pressure_09_p = standard_deviation(pressure_09_plus, pressure_09_mean_p)
sd_pressure_10_p = standard_deviation(pressure_10_plus, pressure_10_mean_p)
sd_pressure_11_p = standard_deviation(pressure_11_plus, pressure_11_mean_p)
sd_pressure_16_p = standard_deviation(pressure_16_plus, pressure_16_mean_p)
sd_pressure_17_p = standard_deviation(pressure_17_plus, pressure_17_mean_p)
sd_pressure_18_p = standard_deviation(pressure_18_plus, pressure_18_mean_p)
sd_pressure_19_p = standard_deviation(pressure_19_plus, pressure_19_mean_p)
sd_pressure_20_p = standard_deviation(pressure_20_plus, pressure_20_mean_p)
sd_pressure_21_p = standard_deviation(pressure_21_plus, pressure_21_mean_p)
sd_pressure_24_p = standard_deviation(pressure_24_plus, pressure_24_mean_p)
sd_pressure_25_p = standard_deviation(pressure_25_plus, pressure_25_mean_p)
sd_pressure_27_p = standard_deviation(pressure_27_plus, pressure_27_mean_p)
sd_pressure_28_p = standard_deviation(pressure_28_plus, pressure_28_mean_p)
sd_pressure_34_p = standard_deviation(pressure_34_plus, pressure_34_mean_p)
sd_pressure_35_p = standard_deviation(pressure_35_plus, pressure_35_mean_p)

# high cmls
sd_pressure_01_p = standard_deviation(pressure_11_plus, pressure_01_mean_p)
sd_pressure_12_p = standard_deviation(pressure_12_plus, pressure_12_mean_p)
sd_pressure_15_p = standard_deviation(pressure_15_plus, pressure_15_mean_p)
sd_pressure_26_p = standard_deviation(pressure_26_plus, pressure_26_mean_p)

sd_pressure_02_m = standard_deviation(pressure_02_minus, pressure_02_mean_m)
sd_pressure_03_m = standard_deviation(pressure_03_minus, pressure_03_mean_m)
sd_pressure_04_m = standard_deviation(pressure_04_minus, pressure_04_mean_m)
sd_pressure_05_m = standard_deviation(pressure_05_minus, pressure_05_mean_m)
sd_pressure_08_m = standard_deviation(pressure_08_minus, pressure_08_mean_m)
sd_pressure_09_m = standard_deviation(pressure_09_minus, pressure_09_mean_m)
sd_pressure_10_m = standard_deviation(pressure_10_minus, pressure_10_mean_m)
sd_pressure_11_m = standard_deviation(pressure_11_minus, pressure_11_mean_m)
sd_pressure_16_m = standard_deviation(pressure_16_minus, pressure_16_mean_m)
sd_pressure_17_m = standard_deviation(pressure_17_minus, pressure_17_mean_m)
sd_pressure_18_m = standard_deviation(pressure_18_minus, pressure_18_mean_m)
sd_pressure_19_m = standard_deviation(pressure_19_minus, pressure_19_mean_m)
sd_pressure_20_m = standard_deviation(pressure_20_minus, pressure_20_mean_m)
sd_pressure_21_m = standard_deviation(pressure_21_minus, pressure_21_mean_m)
sd_pressure_24_m = standard_deviation(pressure_24_minus, pressure_24_mean_m)
sd_pressure_25_m = standard_deviation(pressure_25_minus, pressure_25_mean_m)
sd_pressure_27_m = standard_deviation(pressure_27_minus, pressure_27_mean_m)
sd_pressure_28_m = standard_deviation(pressure_28_minus, pressure_28_mean_m)
sd_pressure_34_m = standard_deviation(pressure_34_minus, pressure_34_mean_m)
sd_pressure_35_m = standard_deviation(pressure_35_minus, pressure_35_mean_m)

# high cmls
sd_pressure_01_m = standard_deviation(pressure_01_minus, pressure_01_mean_m)
sd_pressure_12_m = standard_deviation(pressure_12_minus, pressure_12_mean_m)
sd_pressure_15_m = standard_deviation(pressure_15_minus, pressure_15_mean_m)
sd_pressure_26_m = standard_deviation(pressure_26_minus, pressure_26_mean_m)


pressure_02_median_p = np.median(pressure_02_plus)
pressure_03_median_p = np.median(pressure_03_plus)
pressure_04_median_p = np.median(pressure_04_plus)
pressure_05_median_p = np.median(pressure_05_plus)
pressure_08_median_p = np.median(pressure_08_plus)
pressure_09_median_p = np.median(pressure_09_plus)
pressure_10_median_p = np.median(pressure_10_plus)
pressure_11_median_p = np.median(pressure_11_plus)
pressure_16_median_p = np.median(pressure_16_plus)
pressure_17_median_p = np.median(pressure_17_plus)
pressure_18_median_p = np.median(pressure_18_plus)
pressure_19_median_p = np.median(pressure_19_plus)
pressure_20_median_p = np.median(pressure_20_plus)
pressure_21_median_p = np.median(pressure_21_plus)
pressure_24_median_p = np.median(pressure_24_plus)
pressure_25_median_p = np.median(pressure_25_plus)
pressure_27_median_p = np.median(pressure_27_plus)
pressure_28_median_p = np.median(pressure_28_plus)
pressure_34_median_p = np.median(pressure_34_plus)
pressure_35_median_p = np.median(pressure_35_plus)

#high cmls
pressure_01_median_p = np.median(pressure_01_plus)
pressure_12_median_p = np.median(pressure_12_plus)
pressure_15_median_p = np.median(pressure_15_plus)
pressure_26_median_p = np.median(pressure_26_plus)

#pressure_medians_p = np.array([pressure_02_median_p,pressure_03_median_p,pressure_04_median_p,pressure_05_median_p,pressure_08_median_p,pressure_09_median_p,pressure_10_median_p,pressure_11_median_p,pressure_16_median_p,pressure_17_median_p,pressure_18_median_p,pressure_19_median_p,pressure_20_median_p,pressure_21_median_p,pressure_24_median_p,pressure_25_median_p,pressure_27_median_p,pressure_28_median_p,pressure_34_median_p,pressure_35_median_p])
pressure_medians_p = np.array([pressure_01_median_p,pressure_02_median_p,pressure_03_median_p,pressure_04_median_p,pressure_05_median_p,pressure_08_median_p,pressure_09_median_p,pressure_10_median_p,pressure_11_median_p,pressure_12_median_p,pressure_15_median_p,pressure_16_median_p,pressure_17_median_p,pressure_18_median_p,pressure_19_median_p,pressure_20_median_p,pressure_21_median_p,pressure_24_median_p,pressure_25_median_p,pressure_26_median_p,pressure_27_median_p,pressure_28_median_p,pressure_34_median_p,pressure_35_median_p])


pressure_02_median_m = np.median(pressure_02_minus)
pressure_03_median_m = np.median(pressure_03_minus)
pressure_04_median_m = np.median(pressure_04_minus)
pressure_05_median_m = np.median(pressure_05_minus)
pressure_08_median_m = np.median(pressure_08_minus)
pressure_09_median_m = np.median(pressure_09_minus)
pressure_10_median_m = np.median(pressure_10_minus)
pressure_11_median_m = np.median(pressure_11_minus)
pressure_16_median_m = np.median(pressure_16_minus)
pressure_17_median_m = np.median(pressure_17_minus)
pressure_18_median_m = np.median(pressure_18_minus)
pressure_19_median_m = np.median(pressure_19_minus)
pressure_20_median_m = np.median(pressure_20_minus)
pressure_21_median_m = np.median(pressure_21_minus)
pressure_24_median_m = np.median(pressure_24_minus)
pressure_25_median_m = np.median(pressure_25_minus)
pressure_27_median_m = np.median(pressure_27_minus)
pressure_28_median_m = np.median(pressure_28_minus)
pressure_34_median_m = np.median(pressure_34_minus)
pressure_35_median_m = np.median(pressure_35_minus)

# high cmls
pressure_01_median_m = np.median(pressure_01_minus)
pressure_12_median_m = np.median(pressure_12_minus)
pressure_15_median_m = np.median(pressure_15_minus)
pressure_26_median_m = np.median(pressure_26_minus)

#pressure_medians_m = np.array([pressure_02_median_m,pressure_03_median_m,pressure_04_median_m,pressure_05_median_m,pressure_08_median_m,pressure_09_median_m,pressure_10_median_m,pressure_11_median_m,pressure_16_median_m,pressure_17_median_m,pressure_18_median_m,pressure_19_median_m,pressure_20_median_m,pressure_21_median_m,pressure_24_median_m,pressure_25_median_m,pressure_27_median_m,pressure_28_median_m,pressure_34_median_m,pressure_35_median_m])
pressure_medians_m = np.array([pressure_01_median_m,pressure_02_median_m,pressure_03_median_m,pressure_04_median_m,pressure_05_median_m,pressure_08_median_m,pressure_09_median_m,pressure_10_median_m,pressure_11_median_m,pressure_12_median_m,pressure_15_median_m,pressure_16_median_m,pressure_17_median_m,pressure_18_median_m,pressure_19_median_m,pressure_20_median_m,pressure_21_median_m,pressure_24_median_m,pressure_25_median_m,pressure_26_median_m,pressure_27_median_m,pressure_28_median_m,pressure_34_median_m,pressure_35_median_m])


med_pressure_02_err_p = 1.2533 * sd_pressure_02_p
med_pressure_03_err_p = 1.2533 * sd_pressure_03_p
med_pressure_04_err_p = 1.2533 * sd_pressure_04_p
med_pressure_05_err_p = 1.2533 * sd_pressure_05_p
med_pressure_08_err_p = 1.2533 * sd_pressure_08_p
med_pressure_09_err_p = 1.2533 * sd_pressure_09_p
med_pressure_10_err_p = 1.2533 * sd_pressure_10_p
med_pressure_11_err_p = 1.2533 * sd_pressure_11_p
med_pressure_16_err_p = 1.2533 * sd_pressure_16_p
med_pressure_17_err_p = 1.2533 * sd_pressure_17_p
med_pressure_18_err_p = 1.2533 * sd_pressure_18_p
med_pressure_19_err_p = 1.2533 * sd_pressure_19_p
med_pressure_20_err_p = 1.2533 * sd_pressure_20_p
med_pressure_21_err_p = 1.2533 * sd_pressure_21_p
med_pressure_24_err_p = 1.2533 * sd_pressure_24_p
med_pressure_25_err_p = 1.2533 * sd_pressure_25_p
med_pressure_27_err_p = 1.2533 * sd_pressure_27_p
med_pressure_28_err_p = 1.2533 * sd_pressure_28_p
med_pressure_34_err_p = 1.2533 * sd_pressure_34_p
med_pressure_35_err_p = 1.2533 * sd_pressure_35_p

# high cmls
med_pressure_01_err_p = 1.2533 * sd_pressure_01_p
med_pressure_12_err_p = 1.2533 * sd_pressure_12_p
med_pressure_15_err_p = 1.2533 * sd_pressure_15_p
med_pressure_26_err_p = 1.2533 * sd_pressure_26_p

#pressure_med_errs_p = np.array([med_pressure_02_err_p,med_pressure_03_err_p,med_pressure_04_err_p,med_pressure_05_err_p,med_pressure_08_err_p,med_pressure_09_err_p,med_pressure_10_err_p,med_pressure_11_err_p,med_pressure_16_err_p,med_pressure_17_err_p,med_pressure_18_err_p,med_pressure_19_err_p,med_pressure_20_err_p,med_pressure_21_err_p,med_pressure_24_err_p,med_pressure_25_err_p,med_pressure_27_err_p,med_pressure_28_err_p,med_pressure_34_err_p,med_pressure_35_err_p])
pressure_med_errs_p = np.array([med_pressure_01_err,med_pressure_02_err_p,med_pressure_03_err_p,med_pressure_04_err_p,med_pressure_05_err_p,med_pressure_08_err_p,med_pressure_09_err_p,med_pressure_10_err_p,med_pressure_11_err_p,med_pressure_12_err,med_pressure_15_err,med_pressure_16_err_p,med_pressure_17_err_p,med_pressure_18_err_p,med_pressure_19_err_p,med_pressure_20_err_p,med_pressure_21_err_p,med_pressure_24_err_p,med_pressure_25_err_p,med_pressure_26_err,med_pressure_27_err_p,med_pressure_28_err_p,med_pressure_34_err_p,med_pressure_35_err_p])


med_pressure_02_err_m = 1.2533 * sd_pressure_02_m
med_pressure_03_err_m = 1.2533 * sd_pressure_03_m
med_pressure_04_err_m = 1.2533 * sd_pressure_04_m
med_pressure_05_err_m = 1.2533 * sd_pressure_05_m
med_pressure_08_err_m = 1.2533 * sd_pressure_08_m
med_pressure_09_err_m = 1.2533 * sd_pressure_09_m
med_pressure_10_err_m = 1.2533 * sd_pressure_10_m
med_pressure_11_err_m = 1.2533 * sd_pressure_11_m
med_pressure_16_err_m = 1.2533 * sd_pressure_16_m
med_pressure_17_err_m = 1.2533 * sd_pressure_17_m
med_pressure_18_err_m = 1.2533 * sd_pressure_18_m
med_pressure_19_err_m = 1.2533 * sd_pressure_19_m
med_pressure_20_err_m = 1.2533 * sd_pressure_20_m
med_pressure_21_err_m = 1.2533 * sd_pressure_21_m
med_pressure_24_err_m = 1.2533 * sd_pressure_24_m
med_pressure_25_err_m = 1.2533 * sd_pressure_25_m
med_pressure_27_err_m = 1.2533 * sd_pressure_27_m
med_pressure_28_err_m = 1.2533 * sd_pressure_28_m
med_pressure_34_err_m = 1.2533 * sd_pressure_34_m
med_pressure_35_err_m = 1.2533 * sd_pressure_35_m

# high cmls
med_pressure_01_err_m = 1.2533 * sd_pressure_01_m
med_pressure_12_err_m = 1.2533 * sd_pressure_12_m
med_pressure_15_err_m = 1.2533 * sd_pressure_15_m
med_pressure_26_err_m = 1.2533 * sd_pressure_26_m

#pressure_med_errs_m = np.array([med_pressure_02_err_m,med_pressure_03_err_m,med_pressure_04_err_m,med_pressure_05_err_m,med_pressure_08_err_m,med_pressure_09_err_m,med_pressure_10_err_m,med_pressure_11_err_m,med_pressure_16_err_m,med_pressure_17_err_m,med_pressure_18_err_m,med_pressure_19_err_m,med_pressure_20_err_m,med_pressure_21_err_m,med_pressure_24_err_m,med_pressure_25_err_m,med_pressure_27_err_m,med_pressure_28_err_m,med_pressure_34_err_m,med_pressure_35_err_m])
pressure_med_errs_m = np.array([med_pressure_01_err,med_pressure_02_err_m,med_pressure_03_err_m,med_pressure_04_err_m,med_pressure_05_err_m,med_pressure_08_err_m,med_pressure_09_err_m,med_pressure_10_err_m,med_pressure_11_err_m,med_pressure_12_err,med_pressure_15_err,med_pressure_16_err_m,med_pressure_17_err_m,med_pressure_18_err_m,med_pressure_19_err_m,med_pressure_20_err_m,med_pressure_21_err_m,med_pressure_24_err_m,med_pressure_25_err_m,med_pressure_26_err,med_pressure_27_err_m,med_pressure_28_err_m,med_pressure_34_err_m,med_pressure_35_err_m])



'''
#AVERAGE LL RECONNECTION VOLTAGE
'''

LL_02_mean = np.mean(LL_02)
LL_03_mean = np.mean(LL_03)
LL_04_mean = np.mean(LL_04)
LL_05_mean = np.mean(LL_05)
LL_08_mean = np.mean(LL_08)
LL_09_mean = np.mean(LL_09)
LL_10_mean = np.mean(LL_10)
LL_11_mean = np.mean(LL_11)
LL_16_mean = np.mean(LL_16)
LL_17_mean = np.mean(LL_17)
LL_18_mean = np.mean(LL_18)
LL_19_mean = np.mean(LL_19)
LL_20_mean = np.mean(LL_20)
LL_21_mean = np.mean(LL_21)
LL_24_mean = np.mean(LL_24)
LL_25_mean = np.mean(LL_25)
LL_27_mean = np.mean(LL_27)
LL_28_mean = np.mean(LL_28)
LL_34_mean = np.mean(LL_34)
LL_35_mean = np.mean(LL_35)

# high cmls
LL_01_mean = np.mean(LL_10)
LL_12_mean = np.mean(LL_12)
LL_15_mean = np.mean(LL_15)
LL_26_mean = np.mean(LL_26)

#LL_means = np.array([LL_02_mean,LL_03_mean,LL_04_mean,LL_05_mean,LL_08_mean,LL_09_mean,LL_10_mean,LL_11_mean,LL_16_mean,LL_17_mean,LL_18_mean,LL_19_mean,LL_20_mean,LL_21_mean,LL_24_mean,LL_25_mean,LL_27_mean,LL_28_mean,LL_34_mean,LL_35_mean])
LL_means = np.array([LL_01_mean,LL_02_mean,LL_03_mean,LL_04_mean,LL_05_mean,LL_08_mean,LL_09_mean,LL_10_mean,LL_11_mean,LL_12_mean,LL_15_mean,LL_16_mean,LL_17_mean,LL_18_mean,LL_19_mean,LL_20_mean,LL_21_mean,LL_24_mean,LL_25_mean,LL_26_mean,LL_27_mean,LL_28_mean,LL_34_mean,LL_35_mean])


sd_LL_02 = standard_deviation(LL_02, LL_02_mean)
sd_LL_03 = standard_deviation(LL_03, LL_03_mean)
sd_LL_04 = standard_deviation(LL_04, LL_04_mean)
sd_LL_05 = standard_deviation(LL_05, LL_05_mean)
sd_LL_08 = standard_deviation(LL_08, LL_08_mean)
sd_LL_09 = standard_deviation(LL_09, LL_09_mean)
sd_LL_10 = standard_deviation(LL_10, LL_10_mean)
sd_LL_11 = standard_deviation(LL_11, LL_11_mean)
sd_LL_16 = standard_deviation(LL_16, LL_16_mean)
sd_LL_17 = standard_deviation(LL_17, LL_17_mean)
sd_LL_18 = standard_deviation(LL_18, LL_18_mean)
sd_LL_19 = standard_deviation(LL_19, LL_19_mean)
sd_LL_20 = standard_deviation(LL_20, LL_20_mean)
sd_LL_21 = standard_deviation(LL_21, LL_21_mean)
sd_LL_24 = standard_deviation(LL_24, LL_24_mean)
sd_LL_25 = standard_deviation(LL_25, LL_25_mean)
sd_LL_27 = standard_deviation(LL_27, LL_27_mean)
sd_LL_28 = standard_deviation(LL_28, LL_28_mean)
sd_LL_34 = standard_deviation(LL_34, LL_34_mean)
sd_LL_35 = standard_deviation(LL_35, LL_35_mean)

# high cmls
sd_LL_01 = standard_deviation(LL_01, LL_01_mean)
sd_LL_12 = standard_deviation(LL_12, LL_12_mean)
sd_LL_15 = standard_deviation(LL_15, LL_15_mean)
sd_LL_26 = standard_deviation(LL_26, LL_26_mean)


LL_02_median = np.median(LL_02)
LL_03_median = np.median(LL_03)
LL_04_median = np.median(LL_04)
LL_05_median = np.median(LL_05)
LL_08_median = np.median(LL_08)
LL_09_median = np.median(LL_09)
LL_10_median = np.median(LL_10)
LL_11_median = np.median(LL_11)
LL_16_median = np.median(LL_16)
LL_17_median = np.median(LL_17)
LL_18_median = np.median(LL_18)
LL_19_median = np.median(LL_19)
LL_20_median = np.median(LL_20)
LL_21_median = np.median(LL_21)
LL_24_median = np.median(LL_24)
LL_25_median = np.median(LL_25)
LL_27_median = np.median(LL_27)
LL_28_median = np.median(LL_28)
LL_34_median = np.median(LL_34)
LL_35_median = np.median(LL_35)

# high cmls
LL_01_median = np.median(LL_01)
LL_12_median = np.median(LL_12)
LL_15_median = np.median(LL_15)
LL_26_median = np.median(LL_26)

#LL_medians = np.array([LL_02_median,LL_03_median,LL_04_median,LL_05_median,LL_08_median,LL_09_median,LL_10_median,LL_11_median,LL_16_median,LL_17_median,LL_18_median,LL_19_median,LL_20_median,LL_21_median,LL_24_median,LL_25_median,LL_27_median,LL_28_median,LL_34_median,LL_35_median])
LL_medians = np.array([LL_01_median,LL_02_median,LL_03_median,LL_04_median,LL_05_median,LL_08_median,LL_09_median,LL_10_median,LL_11_median,LL_12_median,LL_15_median,LL_16_median,LL_17_median,LL_18_median,LL_19_median,LL_20_median,LL_21_median,LL_24_median,LL_25_median,LL_26_median,LL_27_median,LL_28_median,LL_34_median,LL_35_median])


med_LL_02_err = 1.2533 * sd_LL_02
med_LL_03_err = 1.2533 * sd_LL_03
med_LL_04_err = 1.2533 * sd_LL_04
med_LL_05_err = 1.2533 * sd_LL_05
med_LL_08_err = 1.2533 * sd_LL_08
med_LL_09_err = 1.2533 * sd_LL_09
med_LL_10_err = 1.2533 * sd_LL_10
med_LL_11_err = 1.2533 * sd_LL_11
med_LL_16_err = 1.2533 * sd_LL_16
med_LL_17_err = 1.2533 * sd_LL_17
med_LL_18_err = 1.2533 * sd_LL_18
med_LL_19_err = 1.2533 * sd_LL_19
med_LL_20_err = 1.2533 * sd_LL_20
med_LL_21_err = 1.2533 * sd_LL_21
med_LL_24_err = 1.2533 * sd_LL_24
med_LL_25_err = 1.2533 * sd_LL_25
med_LL_27_err = 1.2533 * sd_LL_27
med_LL_28_err = 1.2533 * sd_LL_28
med_LL_34_err = 1.2533 * sd_LL_34
med_LL_35_err = 1.2533 * sd_LL_35

# high cmls
med_LL_01_err = 1.2533 * sd_LL_01
med_LL_12_err = 1.2533 * sd_LL_12
med_LL_15_err = 1.2533 * sd_LL_15
med_LL_26_err = 1.2533 * sd_LL_26

#med_LL_errs = np.array([med_LL_02_err,med_LL_03_err,med_LL_04_err,med_LL_05_err,med_LL_08_err,med_LL_09_err,med_LL_10_err,med_LL_11_err,med_LL_16_err,med_LL_17_err,med_LL_18_err,med_LL_19_err,med_LL_20_err,med_LL_21_err,med_LL_24_err,med_LL_25_err,med_LL_27_err,med_LL_28_err,med_LL_34_err,med_LL_35_err])
med_LL_errs = np.array([med_LL_01_err,med_LL_02_err,med_LL_03_err,med_LL_04_err,med_LL_05_err,med_LL_08_err,med_LL_09_err,med_LL_10_err,med_LL_11_err,med_LL_12_err,med_LL_15_err,med_LL_16_err,med_LL_17_err,med_LL_18_err,med_LL_19_err,med_LL_20_err,med_LL_21_err,med_LL_24_err,med_LL_25_err,med_LL_26_err,med_LL_27_err,med_LL_28_err,med_LL_34_err,med_LL_35_err])


'''
#travel time +/- 10% LL
'''
LL_02_mean_p = np.mean(LL_02_plus)
LL_03_mean_p = np.mean(LL_03_plus)
LL_04_mean_p = np.mean(LL_04_plus)
LL_05_mean_p = np.mean(LL_05_plus)
LL_08_mean_p = np.mean(LL_08_plus)
LL_09_mean_p = np.mean(LL_09_plus)
LL_10_mean_p = np.mean(LL_10_plus)
LL_11_mean_p = np.mean(LL_11_plus)
LL_16_mean_p = np.mean(LL_16_plus)
LL_17_mean_p = np.mean(LL_17_plus)
LL_18_mean_p = np.mean(LL_18_plus)
LL_19_mean_p = np.mean(LL_19_plus)
LL_20_mean_p = np.mean(LL_20_plus)
LL_21_mean_p = np.mean(LL_21_plus)
LL_24_mean_p = np.mean(LL_24_plus)
LL_25_mean_p = np.mean(LL_25_plus)
LL_27_mean_p = np.mean(LL_27_plus)
LL_28_mean_p = np.mean(LL_28_plus)
LL_34_mean_p = np.mean(LL_34_plus)
LL_35_mean_p = np.mean(LL_35_plus)

# high cmls
LL_01_mean_p = np.mean(LL_01_plus)
LL_12_mean_p = np.mean(LL_12_plus)
LL_15_mean_p = np.mean(LL_15_plus)
LL_26_mean_p = np.mean(LL_26_plus)

#LL_means_p = np.array([LL_02_mean_p,LL_03_mean_p,LL_04_mean_p,LL_05_mean_p,LL_08_mean_p,LL_09_mean_p,LL_10_mean_p,LL_11_mean_p,LL_16_mean_p,LL_17_mean_p,LL_18_mean_p,LL_19_mean_p,LL_20_mean_p,LL_21_mean_p,LL_24_mean_p,LL_25_mean_p,LL_27_mean_p,LL_28_mean_p,LL_34_mean_p,LL_35_mean_p])
LL_means_p = np.array([LL_01_mean_p,LL_02_mean_p,LL_03_mean_p,LL_04_mean_p,LL_05_mean_p,LL_08_mean_p,LL_09_mean_p,LL_10_mean_p,LL_11_mean_p,LL_12_mean_p,LL_15_mean_p,LL_16_mean_p,LL_17_mean_p,LL_18_mean_p,LL_19_mean_p,LL_20_mean_p,LL_21_mean_p,LL_24_mean_p,LL_25_mean_p,LL_26_mean_p,LL_27_mean_p,LL_28_mean_p,LL_34_mean_p,LL_35_mean_p])


sd_LL_02_p = standard_deviation(LL_02_plus, LL_02_mean_p)
sd_LL_03_p = standard_deviation(LL_03_plus, LL_03_mean_p)
sd_LL_04_p = standard_deviation(LL_04_plus, LL_04_mean_p)
sd_LL_05_p = standard_deviation(LL_05_plus, LL_05_mean_p)
sd_LL_08_p = standard_deviation(LL_08_plus, LL_08_mean_p)
sd_LL_09_p = standard_deviation(LL_09_plus, LL_09_mean_p)
sd_LL_10_p = standard_deviation(LL_10_plus, LL_10_mean_p)
sd_LL_11_p = standard_deviation(LL_11_plus, LL_11_mean_p)
sd_LL_16_p = standard_deviation(LL_16_plus, LL_16_mean_p)
sd_LL_17_p = standard_deviation(LL_17_plus, LL_17_mean_p)
sd_LL_18_p = standard_deviation(LL_18_plus, LL_18_mean_p)
sd_LL_19_p = standard_deviation(LL_19_plus, LL_19_mean_p)
sd_LL_20_p = standard_deviation(LL_20_plus, LL_20_mean_p)
sd_LL_21_p = standard_deviation(LL_21_plus, LL_21_mean_p)
sd_LL_24_p = standard_deviation(LL_24_plus, LL_24_mean_p)
sd_LL_25_p = standard_deviation(LL_25_plus, LL_25_mean_p)
sd_LL_27_p = standard_deviation(LL_27_plus, LL_27_mean_p)
sd_LL_28_p = standard_deviation(LL_28_plus, LL_28_mean_p)
sd_LL_34_p = standard_deviation(LL_34_plus, LL_34_mean_p)
sd_LL_35_p = standard_deviation(LL_35_plus, LL_35_mean_p)

# high cmls
sd_LL_01_p = standard_deviation(LL_01_plus, LL_01_mean_p)
sd_LL_12_p = standard_deviation(LL_12_plus, LL_12_mean_p)
sd_LL_15_p = standard_deviation(LL_15_plus, LL_15_mean_p)
sd_LL_26_p = standard_deviation(LL_26_plus, LL_26_mean_p)


LL_02_median_p = np.median(LL_02_plus)
LL_03_median_p = np.median(LL_03_plus)
LL_04_median_p = np.median(LL_04_plus)
LL_05_median_p = np.median(LL_05_plus)
LL_08_median_p = np.median(LL_08_plus)
LL_09_median_p = np.median(LL_09_plus)
LL_10_median_p = np.median(LL_10_plus)
LL_11_median_p = np.median(LL_11_plus)
LL_16_median_p = np.median(LL_16_plus)
LL_17_median_p = np.median(LL_17_plus)
LL_18_median_p = np.median(LL_18_plus)
LL_19_median_p = np.median(LL_19_plus)
LL_20_median_p = np.median(LL_20_plus)
LL_21_median_p = np.median(LL_21_plus)
LL_24_median_p = np.median(LL_24_plus)
LL_25_median_p = np.median(LL_25_plus)
LL_27_median_p = np.median(LL_27_plus)
LL_28_median_p = np.median(LL_28_plus)
LL_34_median_p = np.median(LL_34_plus)
LL_35_median_p = np.median(LL_35_plus)

# high cmls
LL_01_median_p = np.median(LL_01_plus)
LL_12_median_p = np.median(LL_12_plus)
LL_15_median_p = np.median(LL_15_plus)
LL_26_median_p = np.median(LL_26_plus)

#LL_medians_p = np.array([LL_02_median_p,LL_03_median_p,LL_04_median_p,LL_05_median_p,LL_08_median_p,LL_09_median_p,LL_10_median_p,LL_11_median_p,LL_16_median_p,LL_17_median_p,LL_18_median_p,LL_19_median_p,LL_20_median_p,LL_21_median_p,LL_24_median_p,LL_25_median_p,LL_27_median_p,LL_28_median_p,LL_34_median_p,LL_35_median_p])
LL_medians_p = np.array([LL_01_median_p,LL_02_median_p,LL_03_median_p,LL_04_median_p,LL_05_median_p,LL_08_median_p,LL_09_median_p,LL_10_median_p,LL_11_median_p,LL_12_median_p,LL_15_median_p,LL_16_median_p,LL_17_median_p,LL_18_median_p,LL_19_median_p,LL_20_median_p,LL_21_median_p,LL_24_median_p,LL_25_median_p,LL_26_median_p,LL_27_median_p,LL_28_median_p,LL_34_median_p,LL_35_median_p])


med_LL_02_err_p = 1.2533 * sd_LL_02_p
med_LL_03_err_p = 1.2533 * sd_LL_03_p
med_LL_04_err_p = 1.2533 * sd_LL_04_p
med_LL_05_err_p = 1.2533 * sd_LL_05_p
med_LL_08_err_p = 1.2533 * sd_LL_08_p
med_LL_09_err_p = 1.2533 * sd_LL_09_p
med_LL_10_err_p = 1.2533 * sd_LL_10_p
med_LL_11_err_p = 1.2533 * sd_LL_11_p
med_LL_16_err_p = 1.2533 * sd_LL_16_p
med_LL_17_err_p = 1.2533 * sd_LL_17_p
med_LL_18_err_p = 1.2533 * sd_LL_18_p
med_LL_19_err_p = 1.2533 * sd_LL_19_p
med_LL_20_err_p = 1.2533 * sd_LL_20_p
med_LL_21_err_p = 1.2533 * sd_LL_21_p
med_LL_24_err_p = 1.2533 * sd_LL_24_p
med_LL_25_err_p = 1.2533 * sd_LL_25_p
med_LL_27_err_p = 1.2533 * sd_LL_27_p
med_LL_28_err_p = 1.2533 * sd_LL_28_p
med_LL_34_err_p = 1.2533 * sd_LL_34_p
med_LL_35_err_p = 1.2533 * sd_LL_35_p

# high cmls
med_LL_01_err_p = 1.2533 * sd_LL_01_p
med_LL_12_err_p = 1.2533 * sd_LL_12_p
med_LL_15_err_p = 1.2533 * sd_LL_15_p
med_LL_26_err_p = 1.2533 * sd_LL_26_p

#LL_med_errs_p = np.array([med_LL_02_err_p,med_LL_03_err_p,med_LL_04_err_p,med_LL_05_err_p,med_LL_08_err_p,med_LL_09_err_p,med_LL_10_err_p,med_LL_11_err_p,med_LL_16_err_p,med_LL_17_err_p,med_LL_18_err_p,med_LL_19_err_p,med_LL_20_err_p,med_LL_21_err_p,med_LL_24_err_p,med_LL_25_err_p,med_LL_27_err_p,med_LL_28_err_p,med_LL_34_err_p,med_LL_35_err_p])
LL_med_errs_p = np.array([med_LL_01_err_p,med_LL_02_err_p,med_LL_03_err_p,med_LL_04_err_p,med_LL_05_err_p,med_LL_08_err_p,med_LL_09_err_p,med_LL_10_err_p,med_LL_11_err_p,med_LL_12_err,med_LL_15_err_p,med_LL_16_err_p,med_LL_17_err_p,med_LL_18_err_p,med_LL_19_err_p,med_LL_20_err_p,med_LL_21_err_p,med_LL_24_err_p,med_LL_25_err_p,med_LL_26_err_p,med_LL_27_err_p,med_LL_28_err_p,med_LL_34_err_p,med_LL_35_err_p])


LL_02_mean_m = np.mean(LL_02_minus)
LL_03_mean_m = np.mean(LL_03_minus)
LL_04_mean_m = np.mean(LL_04_minus)
LL_05_mean_m = np.mean(LL_05_minus)
LL_08_mean_m = np.mean(LL_08_minus)
LL_09_mean_m = np.mean(LL_09_minus)
LL_10_mean_m = np.mean(LL_10_minus)
LL_11_mean_m = np.mean(LL_11_minus)
LL_16_mean_m = np.mean(LL_16_minus)
LL_17_mean_m = np.mean(LL_17_minus)
LL_18_mean_m = np.mean(LL_18_minus)
LL_19_mean_m = np.mean(LL_19_minus)
LL_20_mean_m = np.mean(LL_20_minus)
LL_21_mean_m = np.mean(LL_21_minus)
LL_24_mean_m = np.mean(LL_24_minus)
LL_25_mean_m = np.mean(LL_25_minus)
LL_27_mean_m = np.mean(LL_27_minus)
LL_28_mean_m = np.mean(LL_28_minus)
LL_34_mean_m = np.mean(LL_34_minus)
LL_35_mean_m = np.mean(LL_35_minus)

# high cmls
LL_01_mean_m = np.mean(LL_11_minus)
LL_12_mean_m = np.mean(LL_12_minus)
LL_15_mean_m = np.mean(LL_15_minus)
LL_26_mean_m = np.mean(LL_26_minus)

#LL_means_m = np.array([LL_02_mean_m,LL_03_mean_m,LL_04_mean_m,LL_05_mean_m,LL_08_mean_m,LL_09_mean_m,LL_10_mean_m,LL_11_mean_m,LL_16_mean_m,LL_17_mean_m,LL_18_mean_m,LL_19_mean_m,LL_20_mean_m,LL_21_mean_m,LL_24_mean_m,LL_25_mean_m,LL_27_mean_m,LL_28_mean_m,LL_34_mean_m,LL_35_mean_m])
LL_means_m = np.array([LL_01_mean_m,LL_02_mean_m,LL_03_mean_m,LL_04_mean_m,LL_05_mean_m,LL_08_mean_m,LL_09_mean_m,LL_10_mean_m,LL_11_mean_m,LL_12_mean_m,LL_15_mean_m,LL_16_mean_m,LL_17_mean_m,LL_18_mean_m,LL_19_mean_m,LL_20_mean_m,LL_21_mean_m,LL_24_mean_m,LL_25_mean_m,LL_26_mean_m,LL_27_mean_m,LL_28_mean_m,LL_34_mean_m,LL_35_mean_m])


sd_LL_02_m = standard_deviation(LL_02_minus, LL_02_mean_m)
sd_LL_03_m = standard_deviation(LL_03_minus, LL_03_mean_m)
sd_LL_04_m = standard_deviation(LL_04_minus, LL_04_mean_m)
sd_LL_05_m = standard_deviation(LL_05_minus, LL_05_mean_m)
sd_LL_08_m = standard_deviation(LL_08_minus, LL_08_mean_m)
sd_LL_09_m = standard_deviation(LL_09_minus, LL_09_mean_m)
sd_LL_10_m = standard_deviation(LL_10_minus, LL_10_mean_m)
sd_LL_11_m = standard_deviation(LL_11_minus, LL_11_mean_m)
sd_LL_16_m = standard_deviation(LL_16_minus, LL_16_mean_m)
sd_LL_17_m = standard_deviation(LL_17_minus, LL_17_mean_m)
sd_LL_18_m = standard_deviation(LL_18_minus, LL_18_mean_m)
sd_LL_19_m = standard_deviation(LL_19_minus, LL_19_mean_m)
sd_LL_20_m = standard_deviation(LL_20_minus, LL_20_mean_m)
sd_LL_21_m = standard_deviation(LL_21_minus, LL_21_mean_m)
sd_LL_24_m = standard_deviation(LL_24_minus, LL_24_mean_m)
sd_LL_25_m = standard_deviation(LL_25_minus, LL_25_mean_m)
sd_LL_27_m = standard_deviation(LL_27_minus, LL_27_mean_m)
sd_LL_28_m = standard_deviation(LL_28_minus, LL_28_mean_m)
sd_LL_34_m = standard_deviation(LL_34_minus, LL_34_mean_m)
sd_LL_35_m = standard_deviation(LL_35_minus, LL_35_mean_m)

# high cmls
sd_LL_01_m = standard_deviation(LL_01_minus, LL_01_mean_m)
sd_LL_12_m = standard_deviation(LL_12_minus, LL_12_mean_m)
sd_LL_15_m = standard_deviation(LL_15_minus, LL_15_mean_m)
sd_LL_26_m = standard_deviation(LL_26_minus, LL_26_mean_m)


LL_02_median_m = np.median(LL_02_minus)
LL_03_median_m = np.median(LL_03_minus)
LL_04_median_m = np.median(LL_04_minus)
LL_05_median_m = np.median(LL_05_minus)
LL_08_median_m = np.median(LL_08_minus)
LL_09_median_m = np.median(LL_09_minus)
LL_10_median_m = np.median(LL_10_minus)
LL_11_median_m = np.median(LL_11_minus)
LL_16_median_m = np.median(LL_16_minus)
LL_17_median_m = np.median(LL_17_minus)
LL_18_median_m = np.median(LL_18_minus)
LL_19_median_m = np.median(LL_19_minus)
LL_20_median_m = np.median(LL_20_minus)
LL_21_median_m = np.median(LL_21_minus)
LL_24_median_m = np.median(LL_24_minus)
LL_25_median_m = np.median(LL_25_minus)
LL_27_median_m = np.median(LL_27_minus)
LL_28_median_m = np.median(LL_28_minus)
LL_34_median_m = np.median(LL_34_minus)
LL_35_median_m = np.median(LL_35_minus)

# high cmls
LL_01_median_m = np.median(LL_01_minus)
LL_12_median_m = np.median(LL_12_minus)
LL_15_median_m = np.median(LL_15_minus)
LL_26_median_m = np.median(LL_26_minus)

#LL_medians_m = np.array([LL_02_median_m,LL_03_median_m,LL_04_median_m,LL_05_median_m,LL_08_median_m,LL_09_median_m,LL_10_median_m,LL_11_median_m,LL_16_median_m,LL_17_median_m,LL_18_median_m,LL_19_median_m,LL_20_median_m,LL_21_median_m,LL_24_median_m,LL_25_median_m,LL_27_median_m,LL_28_median_m,LL_34_median_m,LL_35_median_m])
LL_medians_m = np.array([LL_01_median_m,LL_02_median_m,LL_03_median_m,LL_04_median_m,LL_05_median_m,LL_08_median_m,LL_09_median_m,LL_10_median_m,LL_11_median_m,LL_12_median_m,LL_15_median_m,LL_16_median_m,LL_17_median_m,LL_18_median_m,LL_19_median_m,LL_20_median_m,LL_21_median_m,LL_24_median_m,LL_25_median_m,LL_26_median_m,LL_27_median_m,LL_28_median_m,LL_34_median_m,LL_35_median_m])


med_LL_02_err_m = 1.2533 * sd_LL_02_m
med_LL_03_err_m = 1.2533 * sd_LL_03_m
med_LL_04_err_m = 1.2533 * sd_LL_04_m
med_LL_05_err_m = 1.2533 * sd_LL_05_m
med_LL_08_err_m = 1.2533 * sd_LL_08_m
med_LL_09_err_m = 1.2533 * sd_LL_09_m
med_LL_10_err_m = 1.2533 * sd_LL_10_m
med_LL_11_err_m = 1.2533 * sd_LL_11_m
med_LL_16_err_m = 1.2533 * sd_LL_16_m
med_LL_17_err_m = 1.2533 * sd_LL_17_m
med_LL_18_err_m = 1.2533 * sd_LL_18_m
med_LL_19_err_m = 1.2533 * sd_LL_19_m
med_LL_20_err_m = 1.2533 * sd_LL_20_m
med_LL_21_err_m = 1.2533 * sd_LL_21_m
med_LL_24_err_m = 1.2533 * sd_LL_24_m
med_LL_25_err_m = 1.2533 * sd_LL_25_m
med_LL_27_err_m = 1.2533 * sd_LL_27_m
med_LL_28_err_m = 1.2533 * sd_LL_28_m
med_LL_34_err_m = 1.2533 * sd_LL_34_m
med_LL_35_err_m = 1.2533 * sd_LL_35_m

# high cmls
med_LL_01_err_m = 1.2533 * sd_LL_01_m
med_LL_12_err_m = 1.2533 * sd_LL_12_m
med_LL_15_err_m = 1.2533 * sd_LL_15_m
med_LL_26_err_m = 1.2533 * sd_LL_26_m

#LL_med_errs_m = np.array([med_LL_02_err_m,med_LL_03_err_m,med_LL_04_err_m,med_LL_05_err_m,med_LL_08_err_m,med_LL_09_err_m,med_LL_10_err_m,med_LL_11_err_m,med_LL_16_err_m,med_LL_17_err_m,med_LL_18_err_m,med_LL_19_err_m,med_LL_20_err_m,med_LL_21_err_m,med_LL_24_err_m,med_LL_25_err_m,med_LL_27_err_m,med_LL_28_err_m,med_LL_34_err_m,med_LL_35_err_m])
LL_med_errs_m = np.array([med_LL_01_err_m,med_LL_02_err_m,med_LL_03_err_m,med_LL_04_err_m,med_LL_05_err_m,med_LL_08_err_m,med_LL_09_err_m,med_LL_10_err_m,med_LL_11_err_m,med_LL_12_err_m,med_LL_15_err_m,med_LL_16_err_m,med_LL_17_err_m,med_LL_18_err_m,med_LL_19_err_m,med_LL_20_err_m,med_LL_21_err_m,med_LL_24_err_m,med_LL_25_err_m,med_LL_26_err_m,med_LL_27_err_m,med_LL_28_err_m,med_LL_34_err_m,med_LL_35_err_m])


'''
#AVERAGE HL RECONNECTION VOLTAGE  - +BY
'''

HL_02_mean_pos = np.mean(HL_02_pos)
HL_03_mean_pos = np.mean(HL_03_pos)
HL_04_mean_pos = np.mean(HL_04_pos)
HL_05_mean_pos = np.mean(HL_05_pos)
HL_08_mean_pos = np.mean(HL_08_pos)
HL_09_mean_pos = np.mean(HL_09_pos)
HL_10_mean_pos = np.mean(HL_10_pos)
HL_11_mean_pos = np.mean(HL_11_pos)
HL_16_mean_pos = np.mean(HL_16_pos)
HL_17_mean_pos = np.mean(HL_17_pos)
HL_18_mean_pos = np.mean(HL_18_pos)
HL_19_mean_pos = np.mean(HL_19_pos)
HL_20_mean_pos = np.mean(HL_20_pos)
HL_21_mean_pos = np.mean(HL_21_pos)
HL_24_mean_pos = np.mean(HL_24_pos)
HL_25_mean_pos = np.mean(HL_25_pos)
HL_27_mean_pos = np.mean(HL_27_pos)
HL_28_mean_pos = np.mean(HL_28_pos)
HL_34_mean_pos = np.mean(HL_34_pos)
HL_35_mean_pos = np.mean(HL_35_pos)

# high cmls
HL_01_mean_pos = np.mean(HL_01_pos)
HL_12_mean_pos = np.mean(HL_12_pos)
HL_15_mean_pos = np.mean(HL_15_pos)
HL_26_mean_pos = np.mean(HL_26_pos)

#HL_means_pos = np.array([HL_02_mean_pos,HL_03_mean_pos,HL_04_mean_pos,HL_05_mean_pos,HL_08_mean_pos,HL_09_mean_pos,HL_10_mean_pos,HL_11_mean_pos,HL_16_mean_pos,HL_17_mean_pos,HL_18_mean_pos,HL_19_mean_pos,HL_20_mean_pos,HL_21_mean_pos,HL_24_mean_pos,HL_25_mean_pos,HL_27_mean_pos,HL_28_mean_pos,HL_34_mean_pos,HL_35_mean_pos])
HL_means_pos = np.array([HL_01_mean_pos,HL_02_mean_pos,HL_03_mean_pos,HL_04_mean_pos,HL_05_mean_pos,HL_08_mean_pos,HL_09_mean_pos,HL_10_mean_pos,HL_11_mean_pos,HL_12_mean_pos,HL_15_mean_pos,HL_16_mean_pos,HL_17_mean_pos,HL_18_mean_pos,HL_19_mean_pos,HL_20_mean_pos,HL_21_mean_pos,HL_24_mean_pos,HL_25_mean_pos,HL_26_mean_pos,HL_27_mean_pos,HL_28_mean_pos,HL_34_mean_pos,HL_35_mean_pos])


sd_HL_02_pos = standard_deviation(HL_02_pos, HL_02_mean_pos)
sd_HL_03_pos = standard_deviation(HL_03_pos, HL_03_mean_pos)
sd_HL_04_pos = standard_deviation(HL_04_pos, HL_04_mean_pos)
sd_HL_05_pos = standard_deviation(HL_05_pos, HL_05_mean_pos)
sd_HL_08_pos = standard_deviation(HL_08_pos, HL_08_mean_pos)
sd_HL_09_pos = standard_deviation(HL_09_pos, HL_09_mean_pos)
sd_HL_10_pos = standard_deviation(HL_10_pos, HL_10_mean_pos)
sd_HL_11_pos = standard_deviation(HL_11_pos, HL_11_mean_pos)
sd_HL_16_pos = standard_deviation(HL_16_pos, HL_16_mean_pos)
sd_HL_17_pos = standard_deviation(HL_17_pos, HL_17_mean_pos)
sd_HL_18_pos = standard_deviation(HL_18_pos, HL_18_mean_pos)
sd_HL_19_pos = standard_deviation(HL_19_pos, HL_19_mean_pos)
sd_HL_20_pos = standard_deviation(HL_20_pos, HL_20_mean_pos)
sd_HL_21_pos = standard_deviation(HL_21_pos, HL_21_mean_pos)
sd_HL_24_pos = standard_deviation(HL_24_pos, HL_24_mean_pos)
sd_HL_25_pos = standard_deviation(HL_25_pos, HL_25_mean_pos)
sd_HL_27_pos = standard_deviation(HL_27_pos, HL_27_mean_pos)
sd_HL_28_pos = standard_deviation(HL_28_pos, HL_28_mean_pos)
sd_HL_34_pos = standard_deviation(HL_34_pos, HL_34_mean_pos)
sd_HL_35_pos = standard_deviation(HL_35_pos, HL_35_mean_pos)

# high cmls
sd_HL_01_pos = standard_deviation(HL_01_pos, HL_01_mean_pos)
sd_HL_12_pos = standard_deviation(HL_12_pos, HL_12_mean_pos)
sd_HL_15_pos = standard_deviation(HL_15_pos, HL_15_mean_pos)
sd_HL_26_pos = standard_deviation(HL_26_pos, HL_26_mean_pos)


HL_02_median_pos = np.median(HL_02_pos)
HL_03_median_pos = np.median(HL_03_pos)
HL_04_median_pos = np.median(HL_04_pos)
HL_05_median_pos = np.median(HL_05_pos)
HL_08_median_pos = np.median(HL_08_pos)
HL_09_median_pos = np.median(HL_09_pos)
HL_10_median_pos = np.median(HL_10_pos)
HL_11_median_pos = np.median(HL_11_pos)
HL_16_median_pos = np.median(HL_16_pos)
HL_17_median_pos = np.median(HL_17_pos)
HL_18_median_pos = np.median(HL_18_pos)
HL_19_median_pos = np.median(HL_19_pos)
HL_20_median_pos = np.median(HL_20_pos)
HL_21_median_pos = np.median(HL_21_pos)
HL_24_median_pos = np.median(HL_24_pos)
HL_25_median_pos = np.median(HL_25_pos)
HL_27_median_pos = np.median(HL_27_pos)
HL_28_median_pos = np.median(HL_28_pos)
HL_34_median_pos = np.median(HL_34_pos)
HL_35_median_pos = np.median(HL_35_pos)

# high cmls
HL_01_median_pos = np.median(HL_01_pos)
HL_12_median_pos = np.median(HL_12_pos)
HL_15_median_pos = np.median(HL_15_pos)
HL_26_median_pos = np.median(HL_26_pos)

#HL_medians_pos = np.array([HL_02_median_pos,HL_03_median_pos,HL_04_median_pos,HL_05_median_pos,HL_08_median_pos,HL_09_median_pos,HL_10_median_pos,HL_11_median_pos,HL_16_median_pos,HL_17_median_pos,HL_18_median_pos,HL_19_median_pos,HL_20_median_pos,HL_21_median_pos,HL_24_median_pos,HL_25_median_pos,HL_27_median_pos,HL_28_median_pos,HL_34_median_pos,HL_35_median_pos])
HL_medians_pos = np.array([HL_01_median_pos,HL_02_median_pos,HL_03_median_pos,HL_04_median_pos,HL_05_median_pos,HL_08_median_pos,HL_09_median_pos,HL_10_median_pos,HL_11_median_pos,HL_12_median_pos,HL_15_median_pos,HL_16_median_pos,HL_17_median_pos,HL_18_median_pos,HL_19_median_pos,HL_20_median_pos,HL_21_median_pos,HL_24_median_pos,HL_25_median_pos,HL_26_median_pos,HL_27_median_pos,HL_28_median_pos,HL_34_median_pos,HL_35_median_pos])


med_HL_02_err_pos = 1.2533 * sd_HL_02_pos
med_HL_03_err_pos = 1.2533 * sd_HL_03_pos
med_HL_04_err_pos = 1.2533 * sd_HL_04_pos
med_HL_05_err_pos = 1.2533 * sd_HL_05_pos
med_HL_08_err_pos = 1.2533 * sd_HL_08_pos
med_HL_09_err_pos = 1.2533 * sd_HL_09_pos
med_HL_10_err_pos = 1.2533 * sd_HL_10_pos
med_HL_11_err_pos = 1.2533 * sd_HL_11_pos
med_HL_16_err_pos = 1.2533 * sd_HL_16_pos
med_HL_17_err_pos = 1.2533 * sd_HL_17_pos
med_HL_18_err_pos = 1.2533 * sd_HL_18_pos
med_HL_19_err_pos = 1.2533 * sd_HL_19_pos
med_HL_20_err_pos = 1.2533 * sd_HL_20_pos
med_HL_21_err_pos = 1.2533 * sd_HL_21_pos
med_HL_24_err_pos = 1.2533 * sd_HL_24_pos
med_HL_25_err_pos = 1.2533 * sd_HL_25_pos
med_HL_27_err_pos = 1.2533 * sd_HL_27_pos
med_HL_28_err_pos = 1.2533 * sd_HL_28_pos
med_HL_34_err_pos = 1.2533 * sd_HL_34_pos
med_HL_35_err_pos = 1.2533 * sd_HL_35_pos

# high cmls
med_HL_01_err_pos = 1.2533 * sd_HL_01_pos
med_HL_12_err_pos = 1.2533 * sd_HL_12_pos
med_HL_15_err_pos = 1.2533 * sd_HL_15_pos
med_HL_26_err_pos = 1.2533 * sd_HL_26_pos

#med_HL_errs_pos = np.array([med_HL_02_err_pos,med_HL_03_err_pos,med_HL_04_err_pos,med_HL_05_err_pos,med_HL_08_err_pos,med_HL_09_err_pos,med_HL_10_err_pos,med_HL_11_err_pos,med_HL_16_err_pos,med_HL_17_err_pos,med_HL_18_err_pos,med_HL_19_err_pos,med_HL_20_err_pos,med_HL_21_err_pos,med_HL_24_err_pos,med_HL_25_err_pos,med_HL_27_err_pos,med_HL_28_err_pos,med_HL_34_err_pos,med_HL_35_err_pos])
med_HL_errs_pos = np.array([med_HL_01_err_pos,med_HL_02_err_pos,med_HL_03_err_pos,med_HL_04_err_pos,med_HL_05_err_pos,med_HL_08_err_pos,med_HL_09_err_pos,med_HL_10_err_pos,med_HL_11_err_pos,med_HL_12_err_pos,med_HL_15_err_pos,med_HL_16_err_pos,med_HL_17_err_pos,med_HL_18_err_pos,med_HL_19_err_pos,med_HL_20_err_pos,med_HL_21_err_pos,med_HL_24_err_pos,med_HL_25_err_pos,med_HL_26_err_pos,med_HL_27_err_pos,med_HL_28_err_pos,med_HL_34_err_pos,med_HL_35_err_pos])


'''
#travel time +/- 10% HL
'''

HL_02_mean_p_pos = np.mean(HL_02_plus_pos)
HL_03_mean_p_pos = np.mean(HL_03_plus_pos)
HL_04_mean_p_pos = np.mean(HL_04_plus_pos)
HL_05_mean_p_pos = np.mean(HL_05_plus_pos)
HL_08_mean_p_pos = np.mean(HL_08_plus_pos)
HL_09_mean_p_pos = np.mean(HL_09_plus_pos)
HL_10_mean_p_pos = np.mean(HL_10_plus_pos)
HL_11_mean_p_pos = np.mean(HL_11_plus_pos)
HL_16_mean_p_pos = np.mean(HL_16_plus_pos)
HL_17_mean_p_pos = np.mean(HL_17_plus_pos)
HL_18_mean_p_pos = np.mean(HL_18_plus_pos)
HL_19_mean_p_pos = np.mean(HL_19_plus_pos)
HL_20_mean_p_pos = np.mean(HL_20_plus_pos)
HL_21_mean_p_pos = np.mean(HL_21_plus_pos)
HL_24_mean_p_pos = np.mean(HL_24_plus_pos)
HL_25_mean_p_pos = np.mean(HL_25_plus_pos)
HL_27_mean_p_pos = np.mean(HL_27_plus_pos)
HL_28_mean_p_pos = np.mean(HL_28_plus_pos)
HL_34_mean_p_pos = np.mean(HL_34_plus_pos)
HL_35_mean_p_pos = np.mean(HL_35_plus_pos)

# high cmls
HL_01_mean_p_pos = np.mean(HL_01_plus_pos)
HL_12_mean_p_pos = np.mean(HL_12_plus_pos)
HL_15_mean_p_pos = np.mean(HL_15_plus_pos)
HL_26_mean_p_pos = np.mean(HL_26_plus_pos)

#HL_means_p_pos = np.array([HL_02_mean_p_pos,HL_03_mean_p_pos,HL_04_mean_p_pos,HL_05_mean_p_pos,HL_08_mean_p_pos,HL_09_mean_p_pos,HL_10_mean_p_pos,HL_11_mean_p_pos,HL_16_mean_p_pos,HL_17_mean_p_pos,HL_18_mean_p_pos,HL_19_mean_p_pos,HL_20_mean_p_pos,HL_21_mean_p_pos,HL_24_mean_p_pos,HL_25_mean_p_pos,HL_27_mean_p_pos,HL_28_mean_p_pos,HL_34_mean_p_pos,HL_35_mean_p_pos])
HL_means_p_pos = np.array([HL_01_mean_p_pos,HL_02_mean_p_pos,HL_03_mean_p_pos,HL_04_mean_p_pos,HL_05_mean_p_pos,HL_08_mean_p_pos,HL_09_mean_p_pos,HL_10_mean_p_pos,HL_11_mean_p_pos,HL_12_mean_p_pos,HL_15_mean_p_pos,HL_16_mean_p_pos,HL_17_mean_p_pos,HL_18_mean_p_pos,HL_19_mean_p_pos,HL_20_mean_p_pos,HL_21_mean_p_pos,HL_24_mean_p_pos,HL_25_mean_p_pos,HL_26_mean_p_pos,HL_27_mean_p_pos,HL_28_mean_p_pos,HL_34_mean_p_pos,HL_35_mean_p_pos])


sd_HL_02_p_pos = standard_deviation(HL_02_plus_pos, HL_02_mean_p_pos)
sd_HL_03_p_pos = standard_deviation(HL_03_plus_pos, HL_03_mean_p_pos)
sd_HL_04_p_pos = standard_deviation(HL_04_plus_pos, HL_04_mean_p_pos)
sd_HL_05_p_pos = standard_deviation(HL_05_plus_pos, HL_05_mean_p_pos)
sd_HL_08_p_pos = standard_deviation(HL_08_plus_pos, HL_08_mean_p_pos)
sd_HL_09_p_pos = standard_deviation(HL_09_plus_pos, HL_09_mean_p_pos)
sd_HL_10_p_pos = standard_deviation(HL_10_plus_pos, HL_10_mean_p_pos)
sd_HL_11_p_pos = standard_deviation(HL_11_plus_pos, HL_11_mean_p_pos)
sd_HL_16_p_pos = standard_deviation(HL_16_plus_pos, HL_16_mean_p_pos)
sd_HL_17_p_pos = standard_deviation(HL_17_plus_pos, HL_17_mean_p_pos)
sd_HL_18_p_pos = standard_deviation(HL_18_plus_pos, HL_18_mean_p_pos)
sd_HL_19_p_pos = standard_deviation(HL_19_plus_pos, HL_19_mean_p_pos)
sd_HL_20_p_pos = standard_deviation(HL_20_plus_pos, HL_20_mean_p_pos)
sd_HL_21_p_pos = standard_deviation(HL_21_plus_pos, HL_21_mean_p_pos)
sd_HL_24_p_pos = standard_deviation(HL_24_plus_pos, HL_24_mean_p_pos)
sd_HL_25_p_pos = standard_deviation(HL_25_plus_pos, HL_25_mean_p_pos)
sd_HL_27_p_pos = standard_deviation(HL_27_plus_pos, HL_27_mean_p_pos)
sd_HL_28_p_pos = standard_deviation(HL_28_plus_pos, HL_28_mean_p_pos)
sd_HL_34_p_pos = standard_deviation(HL_34_plus_pos, HL_34_mean_p_pos)
sd_HL_35_p_pos = standard_deviation(HL_35_plus_pos, HL_35_mean_p_pos)

# high cmls
sd_HL_01_p_pos = standard_deviation(HL_01_plus_pos, HL_01_mean_p_pos)
sd_HL_12_p_pos = standard_deviation(HL_12_plus_pos, HL_12_mean_p_pos)
sd_HL_15_p_pos = standard_deviation(HL_15_plus_pos, HL_15_mean_p_pos)
sd_HL_26_p_pos = standard_deviation(HL_26_plus_pos, HL_26_mean_p_pos)


HL_02_median_p_pos = np.median(HL_02_plus_pos)
HL_03_median_p_pos = np.median(HL_03_plus_pos)
HL_04_median_p_pos = np.median(HL_04_plus_pos)
HL_05_median_p_pos = np.median(HL_05_plus_pos)
HL_08_median_p_pos = np.median(HL_08_plus_pos)
HL_09_median_p_pos = np.median(HL_09_plus_pos)
HL_10_median_p_pos = np.median(HL_10_plus_pos)
HL_11_median_p_pos = np.median(HL_11_plus_pos)
HL_16_median_p_pos = np.median(HL_16_plus_pos)
HL_17_median_p_pos = np.median(HL_17_plus_pos)
HL_18_median_p_pos = np.median(HL_18_plus_pos)
HL_19_median_p_pos = np.median(HL_19_plus_pos)
HL_20_median_p_pos = np.median(HL_20_plus_pos)
HL_21_median_p_pos = np.median(HL_21_plus_pos)
HL_24_median_p_pos = np.median(HL_24_plus_pos)
HL_25_median_p_pos = np.median(HL_25_plus_pos)
HL_27_median_p_pos = np.median(HL_27_plus_pos)
HL_28_median_p_pos = np.median(HL_28_plus_pos)
HL_34_median_p_pos = np.median(HL_34_plus_pos)
HL_35_median_p_pos = np.median(HL_35_plus_pos)

# high cmls
HL_01_median_p_pos = np.median(HL_01_plus_pos)
HL_12_median_p_pos = np.median(HL_12_plus_pos)
HL_15_median_p_pos = np.median(HL_15_plus_pos)
HL_26_median_p_pos = np.median(HL_26_plus_pos)

#HL_medians_p_pos = np.array([HL_02_median_p_pos,HL_03_median_p_pos,HL_04_median_p_pos,HL_05_median_p_pos,HL_08_median_p_pos,HL_09_median_p_pos,HL_10_median_p_pos,HL_11_median_p_pos,HL_16_median_p_pos,HL_17_median_p_pos,HL_18_median_p_pos,HL_19_median_p_pos,HL_20_median_p_pos,HL_21_median_p_pos,HL_24_median_p_pos,HL_25_median_p_pos,HL_27_median_p_pos,HL_28_median_p_pos,HL_34_median_p_pos,HL_35_median_p_pos])
HL_medians_p_pos = np.array([HL_01_median_p_pos,HL_02_median_p_pos,HL_03_median_p_pos,HL_04_median_p_pos,HL_05_median_p_pos,HL_08_median_p_pos,HL_09_median_p_pos,HL_10_median_p_pos,HL_11_median_p_pos,HL_12_median_p_pos,HL_15_median_p_pos,HL_16_median_p_pos,HL_17_median_p_pos,HL_18_median_p_pos,HL_19_median_p_pos,HL_20_median_p_pos,HL_21_median_p_pos,HL_24_median_p_pos,HL_25_median_p_pos,HL_26_median_p_pos,HL_27_median_p_pos,HL_28_median_p_pos,HL_34_median_p_pos,HL_35_median_p_pos])


med_HL_02_err_p_pos = 1.2533 * sd_HL_02_p_pos
med_HL_03_err_p_pos = 1.2533 * sd_HL_03_p_pos
med_HL_04_err_p_pos = 1.2533 * sd_HL_04_p_pos
med_HL_05_err_p_pos = 1.2533 * sd_HL_05_p_pos
med_HL_08_err_p_pos = 1.2533 * sd_HL_08_p_pos
med_HL_09_err_p_pos = 1.2533 * sd_HL_09_p_pos
med_HL_10_err_p_pos = 1.2533 * sd_HL_10_p_pos
med_HL_11_err_p_pos = 1.2533 * sd_HL_11_p_pos
med_HL_16_err_p_pos = 1.2533 * sd_HL_16_p_pos
med_HL_17_err_p_pos = 1.2533 * sd_HL_17_p_pos
med_HL_18_err_p_pos = 1.2533 * sd_HL_18_p_pos
med_HL_19_err_p_pos = 1.2533 * sd_HL_19_p_pos
med_HL_20_err_p_pos = 1.2533 * sd_HL_20_p_pos
med_HL_21_err_p_pos = 1.2533 * sd_HL_21_p_pos
med_HL_24_err_p_pos = 1.2533 * sd_HL_24_p_pos
med_HL_25_err_p_pos = 1.2533 * sd_HL_25_p_pos
med_HL_27_err_p_pos = 1.2533 * sd_HL_27_p_pos
med_HL_28_err_p_pos = 1.2533 * sd_HL_28_p_pos
med_HL_34_err_p_pos = 1.2533 * sd_HL_34_p_pos
med_HL_35_err_p_pos = 1.2533 * sd_HL_35_p_pos

# high cmls
med_HL_01_err_p_pos = 1.2533 * sd_HL_01_p_pos
med_HL_12_err_p_pos = 1.2533 * sd_HL_12_p_pos
med_HL_15_err_p_pos = 1.2533 * sd_HL_15_p_pos
med_HL_26_err_p_pos = 1.2533 * sd_HL_26_p_pos

#HL_med_errs_p_pos = np.array([med_HL_02_err_p_pos,med_HL_03_err_p_pos,med_HL_04_err_p_pos,med_HL_05_err_p_pos,med_HL_08_err_p_pos,med_HL_09_err_p_pos,med_HL_10_err_p_pos,med_HL_11_err_p_pos,med_HL_16_err_p_pos,med_HL_17_err_p_pos,med_HL_18_err_p_pos,med_HL_19_err_p_pos,med_HL_20_err_p_pos,med_HL_21_err_p_pos,med_HL_24_err_p_pos,med_HL_25_err_p_pos,med_HL_27_err_p_pos,med_HL_28_err_p_pos,med_HL_34_err_p_pos,med_HL_35_err_p_pos])
HL_med_errs_p_pos = np.array([med_HL_01_err_p_pos,med_HL_02_err_p_pos,med_HL_03_err_p_pos,med_HL_04_err_p_pos,med_HL_05_err_p_pos,med_HL_08_err_p_pos,med_HL_09_err_p_pos,med_HL_10_err_p_pos,med_HL_11_err_p_pos,med_HL_12_err_p_pos,med_HL_15_err_p_pos,med_HL_16_err_p_pos,med_HL_17_err_p_pos,med_HL_18_err_p_pos,med_HL_19_err_p_pos,med_HL_20_err_p_pos,med_HL_21_err_p_pos,med_HL_24_err_p_pos,med_HL_25_err_p_pos,med_HL_26_err_p_pos,med_HL_27_err_p_pos,med_HL_28_err_p_pos,med_HL_34_err_p_pos,med_HL_35_err_p_pos])


HL_02_mean_m_pos = np.mean(HL_02_minus_pos)
HL_03_mean_m_pos = np.mean(HL_03_minus_pos)
HL_04_mean_m_pos = np.mean(HL_04_minus_pos)
HL_05_mean_m_pos = np.mean(HL_05_minus_pos)
HL_08_mean_m_pos = np.mean(HL_08_minus_pos)
HL_09_mean_m_pos = np.mean(HL_09_minus_pos)
HL_10_mean_m_pos = np.mean(HL_10_minus_pos)
HL_11_mean_m_pos = np.mean(HL_11_minus_pos)
HL_16_mean_m_pos = np.mean(HL_16_minus_pos)
HL_17_mean_m_pos = np.mean(HL_17_minus_pos)
HL_18_mean_m_pos = np.mean(HL_18_minus_pos)
HL_19_mean_m_pos = np.mean(HL_19_minus_pos)
HL_20_mean_m_pos = np.mean(HL_20_minus_pos)
HL_21_mean_m_pos = np.mean(HL_21_minus_pos)
HL_24_mean_m_pos = np.mean(HL_24_minus_pos)
HL_25_mean_m_pos = np.mean(HL_25_minus_pos)
HL_27_mean_m_pos = np.mean(HL_27_minus_pos)
HL_28_mean_m_pos = np.mean(HL_28_minus_pos)
HL_34_mean_m_pos = np.mean(HL_34_minus_pos)
HL_35_mean_m_pos = np.mean(HL_35_minus_pos)

# high cmls
HL_01_mean_m_pos = np.mean(HL_01_minus_pos)
HL_12_mean_m_pos = np.mean(HL_12_minus_pos)
HL_15_mean_m_pos = np.mean(HL_15_minus_pos)
HL_26_mean_m_pos = np.mean(HL_26_minus_pos)

#HL_means_m_pos = np.array([HL_02_mean_m_pos,HL_03_mean_m_pos,HL_04_mean_m_pos,HL_05_mean_m_pos,HL_08_mean_m_pos,HL_09_mean_m_pos,HL_10_mean_m_pos,HL_11_mean_m_pos,HL_16_mean_m_pos,HL_17_mean_m_pos,HL_18_mean_m_pos,HL_19_mean_m_pos,HL_20_mean_m_pos,HL_21_mean_m_pos,HL_24_mean_m_pos,HL_25_mean_m_pos,HL_27_mean_m_pos,HL_28_mean_m_pos,HL_34_mean_m_pos,HL_35_mean_m_pos])
HL_means_m_pos = np.array([HL_01_mean_m_pos,HL_02_mean_m_pos,HL_03_mean_m_pos,HL_04_mean_m_pos,HL_05_mean_m_pos,HL_08_mean_m_pos,HL_09_mean_m_pos,HL_10_mean_m_pos,HL_11_mean_m_pos,HL_12_mean_m_pos,HL_15_mean_m_pos,HL_16_mean_m_pos,HL_17_mean_m_pos,HL_18_mean_m_pos,HL_19_mean_m_pos,HL_20_mean_m_pos,HL_21_mean_m_pos,HL_24_mean_m_pos,HL_25_mean_m_pos,HL_26_mean_m_pos,HL_27_mean_m_pos,HL_28_mean_m_pos,HL_34_mean_m_pos,HL_35_mean_m_pos])


sd_HL_02_m_pos = standard_deviation(HL_02_minus_pos, HL_02_mean_m_pos)
sd_HL_03_m_pos = standard_deviation(HL_03_minus_pos, HL_03_mean_m_pos)
sd_HL_04_m_pos = standard_deviation(HL_04_minus_pos, HL_04_mean_m_pos)
sd_HL_05_m_pos = standard_deviation(HL_05_minus_pos, HL_05_mean_m_pos)
sd_HL_08_m_pos = standard_deviation(HL_08_minus_pos, HL_08_mean_m_pos)
sd_HL_09_m_pos = standard_deviation(HL_09_minus_pos, HL_09_mean_m_pos)
sd_HL_10_m_pos = standard_deviation(HL_10_minus_pos, HL_10_mean_m_pos)
sd_HL_11_m_pos = standard_deviation(HL_11_minus_pos, HL_11_mean_m_pos)
sd_HL_16_m_pos = standard_deviation(HL_16_minus_pos, HL_16_mean_m_pos)
sd_HL_17_m_pos = standard_deviation(HL_17_minus_pos, HL_17_mean_m_pos)
sd_HL_18_m_pos = standard_deviation(HL_18_minus_pos, HL_18_mean_m_pos)
sd_HL_19_m_pos = standard_deviation(HL_19_minus_pos, HL_19_mean_m_pos)
sd_HL_20_m_pos = standard_deviation(HL_20_minus_pos, HL_20_mean_m_pos)
sd_HL_21_m_pos = standard_deviation(HL_21_minus_pos, HL_21_mean_m_pos)
sd_HL_24_m_pos = standard_deviation(HL_24_minus_pos, HL_24_mean_m_pos)
sd_HL_25_m_pos = standard_deviation(HL_25_minus_pos, HL_25_mean_m_pos)
sd_HL_27_m_pos = standard_deviation(HL_27_minus_pos, HL_27_mean_m_pos)
sd_HL_28_m_pos = standard_deviation(HL_28_minus_pos, HL_28_mean_m_pos)
sd_HL_34_m_pos = standard_deviation(HL_34_minus_pos, HL_34_mean_m_pos)
sd_HL_35_m_pos = standard_deviation(HL_35_minus_pos, HL_35_mean_m_pos)

# high cmls
sd_HL_01_m_pos = standard_deviation(HL_01_minus_pos, HL_01_mean_m_pos)
sd_HL_12_m_pos = standard_deviation(HL_12_minus_pos, HL_12_mean_m_pos)
sd_HL_15_m_pos = standard_deviation(HL_15_minus_pos, HL_15_mean_m_pos)
sd_HL_26_m_pos = standard_deviation(HL_26_minus_pos, HL_26_mean_m_pos)


HL_02_median_m_pos = np.median(HL_02_minus_pos)
HL_03_median_m_pos = np.median(HL_03_minus_pos)
HL_04_median_m_pos = np.median(HL_04_minus_pos)
HL_05_median_m_pos = np.median(HL_05_minus_pos)
HL_08_median_m_pos = np.median(HL_08_minus_pos)
HL_09_median_m_pos = np.median(HL_09_minus_pos)
HL_10_median_m_pos = np.median(HL_10_minus_pos)
HL_11_median_m_pos = np.median(HL_11_minus_pos)
HL_16_median_m_pos = np.median(HL_16_minus_pos)
HL_17_median_m_pos = np.median(HL_17_minus_pos)
HL_18_median_m_pos = np.median(HL_18_minus_pos)
HL_19_median_m_pos = np.median(HL_19_minus_pos)
HL_20_median_m_pos = np.median(HL_20_minus_pos)
HL_21_median_m_pos = np.median(HL_21_minus_pos)
HL_24_median_m_pos = np.median(HL_24_minus_pos)
HL_25_median_m_pos = np.median(HL_25_minus_pos)
HL_27_median_m_pos = np.median(HL_27_minus_pos)
HL_28_median_m_pos = np.median(HL_28_minus_pos)
HL_34_median_m_pos = np.median(HL_34_minus_pos)
HL_35_median_m_pos = np.median(HL_35_minus_pos)

# # high cmls
HL_01_median_m_pos = np.median(HL_01_minus_pos)
HL_12_median_m_pos = np.median(HL_12_minus_pos)
HL_15_median_m_pos = np.median(HL_15_minus_pos)
HL_26_median_m_pos = np.median(HL_26_minus_pos)

#HL_medians_m_pos = np.array([HL_02_median_m_pos,HL_03_median_m_pos,HL_04_median_m_pos,HL_05_median_m_pos,HL_08_median_m_pos,HL_09_median_m_pos,HL_10_median_m_pos,HL_11_median_m_pos,HL_16_median_m_pos,HL_17_median_m_pos,HL_18_median_m_pos,HL_19_median_m_pos,HL_20_median_m_pos,HL_21_median_m_pos,HL_24_median_m_pos,HL_25_median_m_pos,HL_27_median_m_pos,HL_28_median_m_pos,HL_34_median_m_pos,HL_35_median_m_pos])
HL_medians_m_pos = np.array([HL_01_median_m_pos,HL_02_median_m_pos,HL_03_median_m_pos,HL_04_median_m_pos,HL_05_median_m_pos,HL_08_median_m_pos,HL_09_median_m_pos,HL_10_median_m_pos,HL_11_median_m_pos,HL_12_median_m_pos,HL_15_median_m_pos,HL_16_median_m_pos,HL_17_median_m_pos,HL_18_median_m_pos,HL_19_median_m_pos,HL_20_median_m_pos,HL_21_median_m_pos,HL_24_median_m_pos,HL_25_median_m_pos,HL_26_median_m_pos,HL_27_median_m_pos,HL_28_median_m_pos,HL_34_median_m_pos,HL_35_median_m_pos])


med_HL_02_err_m_pos = 1.2533 * sd_HL_02_m_pos
med_HL_03_err_m_pos = 1.2533 * sd_HL_03_m_pos
med_HL_04_err_m_pos = 1.2533 * sd_HL_04_m_pos
med_HL_05_err_m_pos = 1.2533 * sd_HL_05_m_pos
med_HL_08_err_m_pos = 1.2533 * sd_HL_08_m_pos
med_HL_09_err_m_pos = 1.2533 * sd_HL_09_m_pos
med_HL_10_err_m_pos = 1.2533 * sd_HL_10_m_pos
med_HL_11_err_m_pos = 1.2533 * sd_HL_11_m_pos
med_HL_16_err_m_pos = 1.2533 * sd_HL_16_m_pos
med_HL_17_err_m_pos = 1.2533 * sd_HL_17_m_pos
med_HL_18_err_m_pos = 1.2533 * sd_HL_18_m_pos
med_HL_19_err_m_pos = 1.2533 * sd_HL_19_m_pos
med_HL_20_err_m_pos = 1.2533 * sd_HL_20_m_pos
med_HL_21_err_m_pos = 1.2533 * sd_HL_21_m_pos
med_HL_24_err_m_pos = 1.2533 * sd_HL_24_m_pos
med_HL_25_err_m_pos = 1.2533 * sd_HL_25_m_pos
med_HL_27_err_m_pos = 1.2533 * sd_HL_27_m_pos
med_HL_28_err_m_pos = 1.2533 * sd_HL_28_m_pos
med_HL_34_err_m_pos = 1.2533 * sd_HL_34_m_pos
med_HL_35_err_m_pos = 1.2533 * sd_HL_35_m_pos

# high cmls
med_HL_01_err_m_pos = 1.2533 * sd_HL_01_m_pos
med_HL_12_err_m_pos = 1.2533 * sd_HL_12_m_pos
med_HL_15_err_m_pos = 1.2533 * sd_HL_15_m_pos
med_HL_26_err_m_pos = 1.2533 * sd_HL_26_m_pos

#HL_med_errs_m_pos = np.array([med_HL_02_err_m_pos,med_HL_03_err_m_pos,med_HL_04_err_m_pos,med_HL_05_err_m_pos,med_HL_08_err_m_pos,med_HL_09_err_m_pos,med_HL_10_err_m_pos,med_HL_11_err_m_pos,med_HL_16_err_m_pos,med_HL_17_err_m_pos,med_HL_18_err_m_pos,med_HL_19_err_m_pos,med_HL_20_err_m_pos,med_HL_21_err_m_pos,med_HL_24_err_m_pos,med_HL_25_err_m_pos,med_HL_27_err_m_pos,med_HL_28_err_m_pos,med_HL_34_err_m_pos,med_HL_35_err_m_pos])
HL_med_errs_m_pos = np.array([med_HL_01_err_m_pos,med_HL_02_err_m_pos,med_HL_03_err_m_pos,med_HL_04_err_m_pos,med_HL_05_err_m_pos,med_HL_08_err_m_pos,med_HL_09_err_m_pos,med_HL_10_err_m_pos,med_HL_11_err_m_pos,med_HL_12_err_m_pos,med_HL_15_err_m_pos,med_HL_16_err_m_pos,med_HL_17_err_m_pos,med_HL_18_err_m_pos,med_HL_19_err_m_pos,med_HL_20_err_m_pos,med_HL_21_err_m_pos,med_HL_24_err_m_pos,med_HL_25_err_m_pos,med_HL_26_err_m_pos,med_HL_27_err_m_pos,med_HL_28_err_m_pos,med_HL_34_err_m_pos,med_HL_35_err_m_pos])


'''
#AVERAGE HL RECONNECTION VOLTAGE - - BY
'''

HL_02_mean_neg = np.mean(HL_02_neg)
HL_03_mean_neg = np.mean(HL_03_neg)
HL_04_mean_neg = np.mean(HL_04_neg)
HL_05_mean_neg = np.mean(HL_05_neg)
HL_08_mean_neg = np.mean(HL_08_neg)
HL_09_mean_neg = np.mean(HL_09_neg)
HL_10_mean_neg = np.mean(HL_10_neg)
HL_11_mean_neg = np.mean(HL_11_neg)
HL_16_mean_neg = np.mean(HL_16_neg)
HL_17_mean_neg = np.mean(HL_17_neg)
HL_18_mean_neg = np.mean(HL_18_neg)
HL_19_mean_neg = np.mean(HL_19_neg)
HL_20_mean_neg = np.mean(HL_20_neg)
HL_21_mean_neg = np.mean(HL_21_neg)
HL_24_mean_neg = np.mean(HL_24_neg)
HL_25_mean_neg = np.mean(HL_25_neg)
HL_27_mean_neg = np.mean(HL_27_neg)
HL_28_mean_neg = np.mean(HL_28_neg)
HL_34_mean_neg = np.mean(HL_34_neg)
HL_35_mean_neg = np.mean(HL_35_neg)

# high cmls
HL_01_mean_neg = np.mean(HL_01_neg)
HL_12_mean_neg = np.mean(HL_12_neg)
HL_15_mean_neg = np.mean(HL_15_neg)
HL_26_mean_neg = np.mean(HL_26_neg)

#HL_means_neg = np.array([HL_02_mean_neg,HL_03_mean_neg,HL_04_mean_neg,HL_05_mean_neg,HL_08_mean_neg,HL_09_mean_neg,HL_10_mean_neg,HL_11_mean_neg,HL_16_mean_neg,HL_17_mean_neg,HL_18_mean_neg,HL_19_mean_neg,HL_20_mean_neg,HL_21_mean_neg,HL_24_mean_neg,HL_25_mean_neg,HL_27_mean_neg,HL_28_mean_neg,HL_34_mean_neg,HL_35_mean_neg])
HL_means_neg = np.array([HL_01_mean_neg,HL_02_mean_neg,HL_03_mean_neg,HL_04_mean_neg,HL_05_mean_neg,HL_08_mean_neg,HL_09_mean_neg,HL_10_mean_neg,HL_11_mean_neg,HL_12_mean_neg,HL_15_mean_neg,HL_16_mean_neg,HL_17_mean_neg,HL_18_mean_neg,HL_19_mean_neg,HL_20_mean_neg,HL_21_mean_neg,HL_24_mean_neg,HL_25_mean_neg,HL_26_mean_neg,HL_27_mean_neg,HL_28_mean_neg,HL_34_mean_neg,HL_35_mean_neg])


sd_HL_02_neg = standard_deviation(HL_02_neg, HL_02_mean_neg)
sd_HL_03_neg = standard_deviation(HL_03_neg, HL_03_mean_neg)
sd_HL_04_neg = standard_deviation(HL_04_neg, HL_04_mean_neg)
sd_HL_05_neg = standard_deviation(HL_05_neg, HL_05_mean_neg)
sd_HL_08_neg = standard_deviation(HL_08_neg, HL_08_mean_neg)
sd_HL_09_neg = standard_deviation(HL_09_neg, HL_09_mean_neg)
sd_HL_10_neg = standard_deviation(HL_10_neg, HL_10_mean_neg)
sd_HL_11_neg = standard_deviation(HL_11_neg, HL_11_mean_neg)
sd_HL_16_neg = standard_deviation(HL_16_neg, HL_16_mean_neg)
sd_HL_17_neg = standard_deviation(HL_17_neg, HL_17_mean_neg)
sd_HL_18_neg = standard_deviation(HL_18_neg, HL_18_mean_neg)
sd_HL_19_neg = standard_deviation(HL_19_neg, HL_19_mean_neg)
sd_HL_20_neg = standard_deviation(HL_20_neg, HL_20_mean_neg)
sd_HL_21_neg = standard_deviation(HL_21_neg, HL_21_mean_neg)
sd_HL_24_neg = standard_deviation(HL_24_neg, HL_24_mean_neg)
sd_HL_25_neg = standard_deviation(HL_25_neg, HL_25_mean_neg)
sd_HL_27_neg = standard_deviation(HL_27_neg, HL_27_mean_neg)
sd_HL_28_neg = standard_deviation(HL_28_neg, HL_28_mean_neg)
sd_HL_34_neg = standard_deviation(HL_34_neg, HL_34_mean_neg)
sd_HL_35_neg = standard_deviation(HL_35_neg, HL_35_mean_neg)

# high cmls
sd_HL_01_neg = standard_deviation(HL_11_neg, HL_01_mean_neg)
sd_HL_12_neg = standard_deviation(HL_12_neg, HL_12_mean_neg)
sd_HL_15_neg = standard_deviation(HL_15_neg, HL_15_mean_neg)
sd_HL_26_neg = standard_deviation(HL_26_neg, HL_26_mean_neg)

HL_02_median_neg = np.median(HL_02_neg)
HL_03_median_neg = np.median(HL_03_neg)
HL_04_median_neg = np.median(HL_04_neg)
HL_05_median_neg = np.median(HL_05_neg)
HL_08_median_neg = np.median(HL_08_neg)
HL_09_median_neg = np.median(HL_09_neg)
HL_10_median_neg = np.median(HL_10_neg)
HL_11_median_neg = np.median(HL_11_neg)
HL_16_median_neg = np.median(HL_16_neg)
HL_17_median_neg = np.median(HL_17_neg)
HL_18_median_neg = np.median(HL_18_neg)
HL_19_median_neg = np.median(HL_19_neg)
HL_20_median_neg = np.median(HL_20_neg)
HL_21_median_neg = np.median(HL_21_neg)
HL_24_median_neg = np.median(HL_24_neg)
HL_25_median_neg = np.median(HL_25_neg)
HL_27_median_neg = np.median(HL_27_neg)
HL_28_median_neg = np.median(HL_28_neg)
HL_34_median_neg = np.median(HL_34_neg)
HL_35_median_neg = np.median(HL_35_neg)

# high cmls
HL_01_median_neg = np.median(HL_01_neg)
HL_12_median_neg = np.median(HL_12_neg)
HL_15_median_neg = np.median(HL_15_neg)
HL_26_median_neg = np.median(HL_26_neg)

#HL_medians_neg = np.array([HL_02_median_neg,HL_03_median_neg,HL_04_median_neg,HL_05_median_neg,HL_08_median_neg,HL_09_median_neg,HL_10_median_neg,HL_11_median_neg,HL_16_median_neg,HL_17_median_neg,HL_18_median_neg,HL_19_median_neg,HL_20_median_neg,HL_21_median_neg,HL_24_median_neg,HL_25_median_neg,HL_27_median_neg,HL_28_median_neg,HL_34_median_neg,HL_35_median_neg])
HL_medians_neg = np.array([HL_01_median_neg,HL_02_median_neg,HL_03_median_neg,HL_04_median_neg,HL_05_median_neg,HL_08_median_neg,HL_09_median_neg,HL_10_median_neg,HL_11_median_neg,HL_12_median_neg,HL_15_median_neg,HL_16_median_neg,HL_17_median_neg,HL_18_median_neg,HL_19_median_neg,HL_20_median_neg,HL_21_median_neg,HL_24_median_neg,HL_25_median_neg,HL_26_median_neg,HL_27_median_neg,HL_28_median_neg,HL_34_median_neg,HL_35_median_neg])


med_HL_02_err_neg = 1.2533 * sd_HL_02_neg
med_HL_03_err_neg = 1.2533 * sd_HL_03_neg
med_HL_04_err_neg = 1.2533 * sd_HL_04_neg
med_HL_05_err_neg = 1.2533 * sd_HL_05_neg
med_HL_08_err_neg = 1.2533 * sd_HL_08_neg
med_HL_09_err_neg = 1.2533 * sd_HL_09_neg
med_HL_10_err_neg = 1.2533 * sd_HL_10_neg
med_HL_11_err_neg = 1.2533 * sd_HL_11_neg
med_HL_16_err_neg = 1.2533 * sd_HL_16_neg
med_HL_17_err_neg = 1.2533 * sd_HL_17_neg
med_HL_18_err_neg = 1.2533 * sd_HL_18_neg
med_HL_19_err_neg = 1.2533 * sd_HL_19_neg
med_HL_20_err_neg = 1.2533 * sd_HL_20_neg
med_HL_21_err_neg = 1.2533 * sd_HL_21_neg
med_HL_24_err_neg = 1.2533 * sd_HL_24_neg
med_HL_25_err_neg = 1.2533 * sd_HL_25_neg
med_HL_27_err_neg = 1.2533 * sd_HL_27_neg
med_HL_28_err_neg = 1.2533 * sd_HL_28_neg
med_HL_34_err_neg = 1.2533 * sd_HL_34_neg
med_HL_35_err_neg = 1.2533 * sd_HL_35_neg

# high cmls
med_HL_01_err_neg = 1.2533 * sd_HL_01_neg
med_HL_12_err_neg = 1.2533 * sd_HL_12_neg
med_HL_15_err_neg = 1.2533 * sd_HL_15_neg
med_HL_26_err_neg = 1.2533 * sd_HL_26_neg

#med_HL_errs_neg = np.array([med_HL_02_err_neg,med_HL_03_err_neg,med_HL_04_err_neg,med_HL_05_err_neg,med_HL_08_err_neg,med_HL_09_err_neg,med_HL_10_err_neg,med_HL_11_err_neg,med_HL_16_err_neg,med_HL_17_err_neg,med_HL_18_err_neg,med_HL_19_err_neg,med_HL_20_err_neg,med_HL_21_err_neg,med_HL_24_err_neg,med_HL_25_err_neg,med_HL_27_err_neg,med_HL_28_err_neg,med_HL_34_err_neg,med_HL_35_err_neg])
med_HL_errs_neg = np.array([med_HL_01_err_neg,med_HL_02_err_neg,med_HL_03_err_neg,med_HL_04_err_neg,med_HL_05_err_neg,med_HL_08_err_neg,med_HL_09_err_neg,med_HL_10_err_neg,med_HL_11_err_neg,med_HL_12_err_neg,med_HL_15_err_neg,med_HL_16_err_neg,med_HL_17_err_neg,med_HL_18_err_neg,med_HL_19_err_neg,med_HL_20_err_neg,med_HL_21_err_neg,med_HL_24_err_neg,med_HL_25_err_neg,med_HL_26_err_neg,med_HL_27_err_neg,med_HL_28_err_neg,med_HL_34_err_neg,med_HL_35_err_neg])



'''
#travel time +/- 10% HL
'''

HL_02_mean_p_neg = np.mean(HL_02_plus_neg)
HL_03_mean_p_neg = np.mean(HL_03_plus_neg)
HL_04_mean_p_neg = np.mean(HL_04_plus_neg)
HL_05_mean_p_neg = np.mean(HL_05_plus_neg)
HL_08_mean_p_neg = np.mean(HL_08_plus_neg)
HL_09_mean_p_neg = np.mean(HL_09_plus_neg)
HL_10_mean_p_neg = np.mean(HL_10_plus_neg)
HL_11_mean_p_neg = np.mean(HL_11_plus_neg)
HL_16_mean_p_neg = np.mean(HL_16_plus_neg)
HL_17_mean_p_neg = np.mean(HL_17_plus_neg)
HL_18_mean_p_neg = np.mean(HL_18_plus_neg)
HL_19_mean_p_neg = np.mean(HL_19_plus_neg)
HL_20_mean_p_neg = np.mean(HL_20_plus_neg)
HL_21_mean_p_neg = np.mean(HL_21_plus_neg)
HL_24_mean_p_neg = np.mean(HL_24_plus_neg)
HL_25_mean_p_neg = np.mean(HL_25_plus_neg)
HL_27_mean_p_neg = np.mean(HL_27_plus_neg)
HL_28_mean_p_neg = np.mean(HL_28_plus_neg)
HL_34_mean_p_neg = np.mean(HL_34_plus_neg)
HL_35_mean_p_neg = np.mean(HL_35_plus_neg)

# high cmls
HL_01_mean_p_neg = np.mean(HL_01_plus_neg)
HL_12_mean_p_neg = np.mean(HL_12_plus_neg)
HL_15_mean_p_neg = np.mean(HL_15_plus_neg)
HL_26_mean_p_neg = np.mean(HL_26_plus_neg)

#HL_means_p_neg = np.array([HL_02_mean_p_neg,HL_03_mean_p_neg,HL_04_mean_p_neg,HL_05_mean_p_neg,HL_08_mean_p_neg,HL_09_mean_p_neg,HL_10_mean_p_neg,HL_11_mean_p_neg,HL_16_mean_p_neg,HL_17_mean_p_neg,HL_18_mean_p_neg,HL_19_mean_p_neg,HL_20_mean_p_neg,HL_21_mean_p_neg,HL_24_mean_p_neg,HL_25_mean_p_neg,HL_27_mean_p_neg,HL_28_mean_p_neg,HL_34_mean_p_neg,HL_35_mean_p_neg])
HL_means_p_neg = np.array([HL_01_mean_p_neg,HL_02_mean_p_neg,HL_03_mean_p_neg,HL_04_mean_p_neg,HL_05_mean_p_neg,HL_08_mean_p_neg,HL_09_mean_p_neg,HL_10_mean_p_neg,HL_11_mean_p_neg,HL_12_mean_p_neg,HL_15_mean_p_neg,HL_16_mean_p_neg,HL_17_mean_p_neg,HL_18_mean_p_neg,HL_19_mean_p_neg,HL_20_mean_p_neg,HL_21_mean_p_neg,HL_24_mean_p_neg,HL_25_mean_p_neg,HL_26_mean_p_neg,HL_27_mean_p_neg,HL_28_mean_p_neg,HL_34_mean_p_neg,HL_35_mean_p_neg])


sd_HL_02_p_neg = standard_deviation(HL_02_plus_neg, HL_02_mean_p_neg)
sd_HL_03_p_neg = standard_deviation(HL_03_plus_neg, HL_03_mean_p_neg)
sd_HL_04_p_neg = standard_deviation(HL_04_plus_neg, HL_04_mean_p_neg)
sd_HL_05_p_neg = standard_deviation(HL_05_plus_neg, HL_05_mean_p_neg)
sd_HL_08_p_neg = standard_deviation(HL_08_plus_neg, HL_08_mean_p_neg)
sd_HL_09_p_neg = standard_deviation(HL_09_plus_neg, HL_09_mean_p_neg)
sd_HL_10_p_neg = standard_deviation(HL_10_plus_neg, HL_10_mean_p_neg)
sd_HL_11_p_neg = standard_deviation(HL_11_plus_neg, HL_11_mean_p_neg)
sd_HL_16_p_neg = standard_deviation(HL_16_plus_neg, HL_16_mean_p_neg)
sd_HL_17_p_neg = standard_deviation(HL_17_plus_neg, HL_17_mean_p_neg)
sd_HL_18_p_neg = standard_deviation(HL_18_plus_neg, HL_18_mean_p_neg)
sd_HL_19_p_neg = standard_deviation(HL_19_plus_neg, HL_19_mean_p_neg)
sd_HL_20_p_neg = standard_deviation(HL_20_plus_neg, HL_20_mean_p_neg)
sd_HL_21_p_neg = standard_deviation(HL_21_plus_neg, HL_21_mean_p_neg)
sd_HL_24_p_neg = standard_deviation(HL_24_plus_neg, HL_24_mean_p_neg)
sd_HL_25_p_neg = standard_deviation(HL_25_plus_neg, HL_25_mean_p_neg)
sd_HL_27_p_neg = standard_deviation(HL_27_plus_neg, HL_27_mean_p_neg)
sd_HL_28_p_neg = standard_deviation(HL_28_plus_neg, HL_28_mean_p_neg)
sd_HL_34_p_neg = standard_deviation(HL_34_plus_neg, HL_34_mean_p_neg)
sd_HL_35_p_neg = standard_deviation(HL_35_plus_neg, HL_35_mean_p_neg)

# high cmls
sd_HL_01_p_neg = standard_deviation(HL_01_plus_neg, HL_01_mean_p_neg)
sd_HL_12_p_neg = standard_deviation(HL_12_plus_neg, HL_12_mean_p_neg)
sd_HL_15_p_neg = standard_deviation(HL_15_plus_neg, HL_15_mean_p_neg)
sd_HL_26_p_neg = standard_deviation(HL_26_plus_neg, HL_26_mean_p_neg)


HL_02_median_p_neg = np.median(HL_02_plus_neg)
HL_03_median_p_neg = np.median(HL_03_plus_neg)
HL_04_median_p_neg = np.median(HL_04_plus_neg)
HL_05_median_p_neg = np.median(HL_05_plus_neg)
HL_08_median_p_neg = np.median(HL_08_plus_neg)
HL_09_median_p_neg = np.median(HL_09_plus_neg)
HL_10_median_p_neg = np.median(HL_10_plus_neg)
HL_11_median_p_neg = np.median(HL_11_plus_neg)
HL_16_median_p_neg = np.median(HL_16_plus_neg)
HL_17_median_p_neg = np.median(HL_17_plus_neg)
HL_18_median_p_neg = np.median(HL_18_plus_neg)
HL_19_median_p_neg = np.median(HL_19_plus_neg)
HL_20_median_p_neg = np.median(HL_20_plus_neg)
HL_21_median_p_neg = np.median(HL_21_plus_neg)
HL_24_median_p_neg = np.median(HL_24_plus_neg)
HL_25_median_p_neg = np.median(HL_25_plus_neg)
HL_27_median_p_neg = np.median(HL_27_plus_neg)
HL_28_median_p_neg = np.median(HL_28_plus_neg)
HL_34_median_p_neg = np.median(HL_34_plus_neg)
HL_35_median_p_neg = np.median(HL_35_plus_neg)

# high cmls
HL_01_median_p_neg = np.median(HL_01_plus_neg)
HL_12_median_p_neg = np.median(HL_12_plus_neg)
HL_15_median_p_neg = np.median(HL_15_plus_neg)
HL_26_median_p_neg = np.median(HL_26_plus_neg)

#HL_medians_p_neg = np.array([HL_02_median_p_neg,HL_03_median_p_neg,HL_04_median_p_neg,HL_05_median_p_neg,HL_08_median_p_neg,HL_09_median_p_neg,HL_10_median_p_neg,HL_11_median_p_neg,HL_16_median_p_neg,HL_17_median_p_neg,HL_18_median_p_neg,HL_19_median_p_neg,HL_20_median_p_neg,HL_21_median_p_neg,HL_24_median_p_neg,HL_25_median_p_neg,HL_27_median_p_neg,HL_28_median_p_neg,HL_34_median_p_neg,HL_35_median_p_neg])
HL_medians_p_neg = np.array([HL_01_median_p_neg,HL_02_median_p_neg,HL_03_median_p_neg,HL_04_median_p_neg,HL_05_median_p_neg,HL_08_median_p_neg,HL_09_median_p_neg,HL_10_median_p_neg,HL_11_median_p_neg,HL_12_median_p_neg,HL_15_median_p_neg,HL_16_median_p_neg,HL_17_median_p_neg,HL_18_median_p_neg,HL_19_median_p_neg,HL_20_median_p_neg,HL_21_median_p_neg,HL_24_median_p_neg,HL_25_median_p_neg,HL_26_median_p_neg,HL_27_median_p_neg,HL_28_median_p_neg,HL_34_median_p_neg,HL_35_median_p_neg])


med_HL_02_err_p_neg = 1.2533 * sd_HL_02_p_neg
med_HL_03_err_p_neg = 1.2533 * sd_HL_03_p_neg
med_HL_04_err_p_neg = 1.2533 * sd_HL_04_p_neg
med_HL_05_err_p_neg = 1.2533 * sd_HL_05_p_neg
med_HL_08_err_p_neg = 1.2533 * sd_HL_08_p_neg
med_HL_09_err_p_neg = 1.2533 * sd_HL_09_p_neg
med_HL_10_err_p_neg = 1.2533 * sd_HL_10_p_neg
med_HL_11_err_p_neg = 1.2533 * sd_HL_11_p_neg
med_HL_16_err_p_neg = 1.2533 * sd_HL_16_p_neg
med_HL_17_err_p_neg = 1.2533 * sd_HL_17_p_neg
med_HL_18_err_p_neg = 1.2533 * sd_HL_18_p_neg
med_HL_19_err_p_neg = 1.2533 * sd_HL_19_p_neg
med_HL_20_err_p_neg = 1.2533 * sd_HL_20_p_neg
med_HL_21_err_p_neg = 1.2533 * sd_HL_21_p_neg
med_HL_24_err_p_neg = 1.2533 * sd_HL_24_p_neg
med_HL_25_err_p_neg = 1.2533 * sd_HL_25_p_neg
med_HL_27_err_p_neg = 1.2533 * sd_HL_27_p_neg
med_HL_28_err_p_neg = 1.2533 * sd_HL_28_p_neg
med_HL_34_err_p_neg = 1.2533 * sd_HL_34_p_neg
med_HL_35_err_p_neg = 1.2533 * sd_HL_35_p_neg

# high cmls
med_HL_01_err_p_neg = 1.2533 * sd_HL_01_p_neg
med_HL_12_err_p_neg = 1.2533 * sd_HL_12_p_neg
med_HL_15_err_p_neg = 1.2533 * sd_HL_15_p_neg
med_HL_26_err_p_neg = 1.2533 * sd_HL_26_p_neg

#HL_med_errs_p_neg = np.array([med_HL_02_err_p_neg,med_HL_03_err_p_neg,med_HL_04_err_p_neg,med_HL_05_err_p_neg,med_HL_08_err_p_neg,med_HL_09_err_p_neg,med_HL_10_err_p_neg,med_HL_11_err_p_neg,med_HL_16_err_p_neg,med_HL_17_err_p_neg,med_HL_18_err_p_neg,med_HL_19_err_p_neg,med_HL_20_err_p_neg,med_HL_21_err_p_neg,med_HL_24_err_p_neg,med_HL_25_err_p_neg,med_HL_27_err_p_neg,med_HL_28_err_p_neg,med_HL_34_err_p_neg,med_HL_35_err_p_neg])
HL_med_errs_p_neg = np.array([med_HL_01_err_p_neg,med_HL_02_err_p_neg,med_HL_03_err_p_neg,med_HL_04_err_p_neg,med_HL_05_err_p_neg,med_HL_08_err_p_neg,med_HL_09_err_p_neg,med_HL_10_err_p_neg,med_HL_11_err_p_neg,med_HL_12_err_p_neg,med_HL_15_err_p_neg,med_HL_16_err_p_neg,med_HL_17_err_p_neg,med_HL_18_err_p_neg,med_HL_19_err_p_neg,med_HL_20_err_p_neg,med_HL_21_err_p_neg,med_HL_24_err_p_neg,med_HL_25_err_p_neg,med_HL_26_err_p_neg,med_HL_27_err_p_neg,med_HL_28_err_p_neg,med_HL_34_err_p_neg,med_HL_35_err_p_neg])


HL_02_mean_m_neg = np.mean(HL_02_minus_neg)
HL_03_mean_m_neg = np.mean(HL_03_minus_neg)
HL_04_mean_m_neg = np.mean(HL_04_minus_neg)
HL_05_mean_m_neg = np.mean(HL_05_minus_neg)
HL_08_mean_m_neg = np.mean(HL_08_minus_neg)
HL_09_mean_m_neg = np.mean(HL_09_minus_neg)
HL_10_mean_m_neg = np.mean(HL_10_minus_neg)
HL_11_mean_m_neg = np.mean(HL_11_minus_neg)
HL_16_mean_m_neg = np.mean(HL_16_minus_neg)
HL_17_mean_m_neg = np.mean(HL_17_minus_neg)
HL_18_mean_m_neg = np.mean(HL_18_minus_neg)
HL_19_mean_m_neg = np.mean(HL_19_minus_neg)
HL_20_mean_m_neg = np.mean(HL_20_minus_neg)
HL_21_mean_m_neg = np.mean(HL_21_minus_neg)
HL_24_mean_m_neg = np.mean(HL_24_minus_neg)
HL_25_mean_m_neg = np.mean(HL_25_minus_neg)
HL_27_mean_m_neg = np.mean(HL_27_minus_neg)
HL_28_mean_m_neg = np.mean(HL_28_minus_neg)
HL_34_mean_m_neg = np.mean(HL_34_minus_neg)
HL_35_mean_m_neg = np.mean(HL_35_minus_neg)

# high cmls
HL_01_mean_m_neg = np.mean(HL_01_minus_neg)
HL_12_mean_m_neg = np.mean(HL_12_minus_neg)
HL_15_mean_m_neg = np.mean(HL_15_minus_neg)
HL_26_mean_m_neg = np.mean(HL_26_minus_neg)

#HL_means_m_neg = np.array([HL_02_mean_m_neg,HL_03_mean_m_neg,HL_04_mean_m_neg,HL_05_mean_m_neg,HL_08_mean_m_neg,HL_09_mean_m_neg,HL_10_mean_m_neg,HL_11_mean_m_neg,HL_16_mean_m_neg,HL_17_mean_m_neg,HL_18_mean_m_neg,HL_19_mean_m_neg,HL_20_mean_m_neg,HL_21_mean_m_neg,HL_24_mean_m_neg,HL_25_mean_m_neg,HL_27_mean_m_neg,HL_28_mean_m_neg,HL_34_mean_m_neg,HL_35_mean_m_neg])
HL_means_m_neg = np.array([HL_01_mean_m_neg,HL_02_mean_m_neg,HL_03_mean_m_neg,HL_04_mean_m_neg,HL_05_mean_m_neg,HL_08_mean_m_neg,HL_09_mean_m_neg,HL_10_mean_m_neg,HL_11_mean_m_neg,HL_12_mean_m_neg,HL_15_mean_m_neg,HL_16_mean_m_neg,HL_17_mean_m_neg,HL_18_mean_m_neg,HL_19_mean_m_neg,HL_20_mean_m_neg,HL_21_mean_m_neg,HL_24_mean_m_neg,HL_25_mean_m_neg,HL_26_mean_m_neg,HL_27_mean_m_neg,HL_28_mean_m_neg,HL_34_mean_m_neg,HL_35_mean_m_neg])


sd_HL_02_m_neg = standard_deviation(HL_02_minus_neg, HL_02_mean_m_neg)
sd_HL_03_m_neg = standard_deviation(HL_03_minus_neg, HL_03_mean_m_neg)
sd_HL_04_m_neg = standard_deviation(HL_04_minus_neg, HL_04_mean_m_neg)
sd_HL_05_m_neg = standard_deviation(HL_05_minus_neg, HL_05_mean_m_neg)
sd_HL_08_m_neg = standard_deviation(HL_08_minus_neg, HL_08_mean_m_neg)
sd_HL_09_m_neg = standard_deviation(HL_09_minus_neg, HL_09_mean_m_neg)
sd_HL_10_m_neg = standard_deviation(HL_10_minus_neg, HL_10_mean_m_neg)
sd_HL_11_m_neg = standard_deviation(HL_11_minus_neg, HL_11_mean_m_neg)
sd_HL_16_m_neg = standard_deviation(HL_16_minus_neg, HL_16_mean_m_neg)
sd_HL_17_m_neg = standard_deviation(HL_17_minus_neg, HL_17_mean_m_neg)
sd_HL_18_m_neg = standard_deviation(HL_18_minus_neg, HL_18_mean_m_neg)
sd_HL_19_m_neg = standard_deviation(HL_19_minus_neg, HL_19_mean_m_neg)
sd_HL_20_m_neg = standard_deviation(HL_20_minus_neg, HL_20_mean_m_neg)
sd_HL_21_m_neg = standard_deviation(HL_21_minus_neg, HL_21_mean_m_neg)
sd_HL_24_m_neg = standard_deviation(HL_24_minus_neg, HL_24_mean_m_neg)
sd_HL_25_m_neg = standard_deviation(HL_25_minus_neg, HL_25_mean_m_neg)
sd_HL_27_m_neg = standard_deviation(HL_27_minus_neg, HL_27_mean_m_neg)
sd_HL_28_m_neg = standard_deviation(HL_28_minus_neg, HL_28_mean_m_neg)
sd_HL_34_m_neg = standard_deviation(HL_34_minus_neg, HL_34_mean_m_neg)
sd_HL_35_m_neg = standard_deviation(HL_35_minus_neg, HL_35_mean_m_neg)

# high cmls
sd_HL_01_m_neg = standard_deviation(HL_01_minus_neg, HL_01_mean_m_neg)
sd_HL_12_m_neg = standard_deviation(HL_12_minus_neg, HL_12_mean_m_neg)
sd_HL_15_m_neg = standard_deviation(HL_15_minus_neg, HL_15_mean_m_neg)
sd_HL_26_m_neg = standard_deviation(HL_26_minus_neg, HL_26_mean_m_neg)


HL_02_median_m_neg = np.median(HL_02_minus_neg)
HL_03_median_m_neg = np.median(HL_03_minus_neg)
HL_04_median_m_neg = np.median(HL_04_minus_neg)
HL_05_median_m_neg = np.median(HL_05_minus_neg)
HL_08_median_m_neg = np.median(HL_08_minus_neg)
HL_09_median_m_neg = np.median(HL_09_minus_neg)
HL_10_median_m_neg = np.median(HL_10_minus_neg)
HL_11_median_m_neg = np.median(HL_11_minus_neg)
HL_16_median_m_neg = np.median(HL_16_minus_neg)
HL_17_median_m_neg = np.median(HL_17_minus_neg)
HL_18_median_m_neg = np.median(HL_18_minus_neg)
HL_19_median_m_neg = np.median(HL_19_minus_neg)
HL_20_median_m_neg = np.median(HL_20_minus_neg)
HL_21_median_m_neg = np.median(HL_21_minus_neg)
HL_24_median_m_neg = np.median(HL_24_minus_neg)
HL_25_median_m_neg = np.median(HL_25_minus_neg)
HL_27_median_m_neg = np.median(HL_27_minus_neg)
HL_28_median_m_neg = np.median(HL_28_minus_neg)
HL_34_median_m_neg = np.median(HL_34_minus_neg)
HL_35_median_m_neg = np.median(HL_35_minus_neg)

# high cmls
HL_01_median_m_neg = np.median(HL_01_minus_neg)
HL_12_median_m_neg = np.median(HL_12_minus_neg)
HL_15_median_m_neg = np.median(HL_15_minus_neg)
HL_26_median_m_neg = np.median(HL_26_minus_neg)

#HL_medians_m_neg = np.array([HL_02_median_m_neg,HL_03_median_m_neg,HL_04_median_m_neg,HL_05_median_m_neg,HL_08_median_m_neg,HL_09_median_m_neg,HL_10_median_m_neg,HL_11_median_m_neg,HL_16_median_m_neg,HL_17_median_m_neg,HL_18_median_m_neg,HL_19_median_m_neg,HL_20_median_m_neg,HL_21_median_m_neg,HL_24_median_m_neg,HL_25_median_m_neg,HL_27_median_m_neg,HL_28_median_m_neg,HL_34_median_m_neg,HL_35_median_m_neg])
HL_medians_m_neg = np.array([HL_01_median_m_neg,HL_02_median_m_neg,HL_03_median_m_neg,HL_04_median_m_neg,HL_05_median_m_neg,HL_08_median_m_neg,HL_09_median_m_neg,HL_10_median_m_neg,HL_11_median_m_neg,HL_12_median_m_neg,HL_15_median_m_neg,HL_16_median_m_neg,HL_17_median_m_neg,HL_18_median_m_neg,HL_19_median_m_neg,HL_20_median_m_neg,HL_21_median_m_neg,HL_24_median_m_neg,HL_25_median_m_neg,HL_26_median_m_neg,HL_27_median_m_neg,HL_28_median_m_neg,HL_34_median_m_neg,HL_35_median_m_neg])


med_HL_02_err_m_neg = 1.2533 * sd_HL_02_m_neg
med_HL_03_err_m_neg = 1.2533 * sd_HL_03_m_neg
med_HL_04_err_m_neg = 1.2533 * sd_HL_04_m_neg
med_HL_05_err_m_neg = 1.2533 * sd_HL_05_m_neg
med_HL_08_err_m_neg = 1.2533 * sd_HL_08_m_neg
med_HL_09_err_m_neg = 1.2533 * sd_HL_09_m_neg
med_HL_10_err_m_neg = 1.2533 * sd_HL_10_m_neg
med_HL_11_err_m_neg = 1.2533 * sd_HL_11_m_neg
med_HL_16_err_m_neg = 1.2533 * sd_HL_16_m_neg
med_HL_17_err_m_neg = 1.2533 * sd_HL_17_m_neg
med_HL_18_err_m_neg = 1.2533 * sd_HL_18_m_neg
med_HL_19_err_m_neg = 1.2533 * sd_HL_19_m_neg
med_HL_20_err_m_neg = 1.2533 * sd_HL_20_m_neg
med_HL_21_err_m_neg = 1.2533 * sd_HL_21_m_neg
med_HL_24_err_m_neg = 1.2533 * sd_HL_24_m_neg
med_HL_25_err_m_neg = 1.2533 * sd_HL_25_m_neg
med_HL_27_err_m_neg = 1.2533 * sd_HL_27_m_neg
med_HL_28_err_m_neg = 1.2533 * sd_HL_28_m_neg
med_HL_34_err_m_neg = 1.2533 * sd_HL_34_m_neg
med_HL_35_err_m_neg = 1.2533 * sd_HL_35_m_neg

# high cmls
med_HL_01_err_m_neg = 1.2533 * sd_HL_01_m_neg
med_HL_12_err_m_neg = 1.2533 * sd_HL_12_m_neg
med_HL_15_err_m_neg = 1.2533 * sd_HL_15_m_neg
med_HL_26_err_m_neg = 1.2533 * sd_HL_26_m_neg

#HL_med_errs_m_neg = np.array([med_HL_02_err_m_neg,med_HL_03_err_m_neg,med_HL_04_err_m_neg,med_HL_05_err_m_neg,med_HL_08_err_m_neg,med_HL_09_err_m_neg,med_HL_10_err_m_neg,med_HL_11_err_m_neg,med_HL_16_err_m_neg,med_HL_17_err_m_neg,med_HL_18_err_m_neg,med_HL_19_err_m_neg,med_HL_20_err_m_neg,med_HL_21_err_m_neg,med_HL_24_err_m_neg,med_HL_25_err_m_neg,med_HL_27_err_m_neg,med_HL_28_err_m_neg,med_HL_34_err_m_neg,med_HL_35_err_m_neg])
HL_med_errs_m_neg = np.array([med_HL_01_err_m_neg,med_HL_02_err_m_neg,med_HL_03_err_m_neg,med_HL_04_err_m_neg,med_HL_05_err_m_neg,med_HL_08_err_m_neg,med_HL_09_err_m_neg,med_HL_10_err_m_neg,med_HL_11_err_m_neg,med_HL_12_err_m_neg,med_HL_15_err_m_neg,med_HL_16_err_m_neg,med_HL_17_err_m_neg,med_HL_18_err_m_neg,med_HL_19_err_m_neg,med_HL_20_err_m_neg,med_HL_21_err_m_neg,med_HL_24_err_m_neg,med_HL_25_err_m_neg,med_HL_26_err_m_neg,med_HL_27_err_m_neg,med_HL_28_err_m_neg,med_HL_34_err_m_neg,med_HL_35_err_m_neg])


'''
AVERAGE GERSHMAN RECONNECTION VOLTAGE
'''

gersh_01_mean = np.mean(gersh_01)
gersh_12_mean = np.mean(gersh_12)
gersh_15_mean = np.mean(gersh_15)
gersh_26_mean = np.mean(gersh_26)

gersh_02_mean = np.mean(gersh_02)
gersh_03_mean = np.mean(gersh_03)
gersh_04_mean = np.mean(gersh_04)
gersh_05_mean = np.mean(gersh_05)
gersh_08_mean = np.mean(gersh_08)
gersh_09_mean = np.mean(gersh_09)
gersh_10_mean = np.mean(gersh_10)
gersh_11_mean = np.mean(gersh_11)
gersh_16_mean = np.mean(gersh_16)
gersh_17_mean = np.mean(gersh_17)
gersh_18_mean = np.mean(gersh_18)
gersh_19_mean = np.mean(gersh_19)
gersh_20_mean = np.mean(gersh_20)
gersh_21_mean = np.mean(gersh_21)
gersh_24_mean = np.mean(gersh_24)
gersh_25_mean = np.mean(gersh_25)
gersh_27_mean = np.mean(gersh_27)
gersh_28_mean = np.mean(gersh_28)
gersh_34_mean = np.mean(gersh_34)
gersh_35_mean = np.mean(gersh_35)

gersh_means = np.array([gersh_01_mean,gersh_02_mean,gersh_03_mean,gersh_04_mean,gersh_05_mean,gersh_08_mean,gersh_09_mean,gersh_10_mean,gersh_11_mean,gersh_12_mean,gersh_15_mean,gersh_16_mean,gersh_17_mean,gersh_18_mean,gersh_19_mean,gersh_20_mean,gersh_21_mean,gersh_24_mean,gersh_25_mean,gersh_26_mean,gersh_27_mean,gersh_28_mean,gersh_34_mean,gersh_35_mean])

sd_gersh_01 = standard_deviation(gersh_01,gersh_01_mean)
sd_gersh_12 = standard_deviation(gersh_12,gersh_12_mean)
sd_gersh_15 = standard_deviation(gersh_15,gersh_15_mean)
sd_gersh_26 = standard_deviation(gersh_26,gersh_26_mean)

sd_gersh_02 = standard_deviation(gersh_02,gersh_02_mean)
sd_gersh_03 = standard_deviation(gersh_03,gersh_03_mean)
sd_gersh_04 = standard_deviation(gersh_04,gersh_04_mean)
sd_gersh_05 = standard_deviation(gersh_05,gersh_05_mean)
sd_gersh_08 = standard_deviation(gersh_08,gersh_08_mean)
sd_gersh_09 = standard_deviation(gersh_09,gersh_09_mean)
sd_gersh_10 = standard_deviation(gersh_10,gersh_10_mean)
sd_gersh_11 = standard_deviation(gersh_11,gersh_11_mean)
sd_gersh_16 = standard_deviation(gersh_16,gersh_16_mean)
sd_gersh_17 = standard_deviation(gersh_17,gersh_17_mean)
sd_gersh_18 = standard_deviation(gersh_18,gersh_18_mean)
sd_gersh_19 = standard_deviation(gersh_19,gersh_19_mean)
sd_gersh_20 = standard_deviation(gersh_20,gersh_20_mean)
sd_gersh_21 = standard_deviation(gersh_21,gersh_21_mean)
sd_gersh_24 = standard_deviation(gersh_24,gersh_24_mean)
sd_gersh_25 = standard_deviation(gersh_25,gersh_25_mean)
sd_gersh_27 = standard_deviation(gersh_27,gersh_27_mean)
sd_gersh_28 = standard_deviation(gersh_28,gersh_28_mean)
sd_gersh_34 = standard_deviation(gersh_34,gersh_34_mean)
sd_gersh_35 = standard_deviation(gersh_35,gersh_35_mean)


gersh_01_median = np.median(gersh_01)
gersh_12_median = np.median(gersh_12)
gersh_15_median = np.median(gersh_15)
gersh_26_median = np.median(gersh_26)

gersh_02_median = np.median(gersh_02)
gersh_03_median = np.median(gersh_03)
gersh_04_median = np.median(gersh_04)
gersh_05_median = np.median(gersh_05)
gersh_08_median = np.median(gersh_08)
gersh_09_median = np.median(gersh_09)
gersh_10_median = np.median(gersh_10)
gersh_11_median = np.median(gersh_11)
gersh_16_median = np.median(gersh_16)
gersh_17_median = np.median(gersh_17)
gersh_18_median = np.median(gersh_18)
gersh_19_median = np.median(gersh_19)
gersh_20_median = np.median(gersh_20)
gersh_21_median = np.median(gersh_21)
gersh_24_median = np.median(gersh_24)
gersh_25_median = np.median(gersh_25)
gersh_27_median = np.median(gersh_27)
gersh_28_median = np.median(gersh_28)
gersh_34_median = np.median(gersh_34)
gersh_35_median = np.median(gersh_35)

gersh_medians = np.array([gersh_01_median,gersh_02_median,gersh_03_median,gersh_04_median,gersh_05_median,gersh_08_median,gersh_09_median,gersh_10_median,gersh_11_median,gersh_12_median,gersh_15_median,gersh_16_median,gersh_17_median,gersh_18_median,gersh_19_median,gersh_20_median,gersh_21_median,gersh_24_median,gersh_25_median,gersh_26_median,gersh_27_median,gersh_28_median,gersh_34_median,gersh_35_median])

med_gersh_err_01 = 1.2533 * sd_gersh_01
med_gersh_err_12 = 1.2533 * sd_gersh_12
med_gersh_err_15 = 1.2533 * sd_gersh_15
med_gersh_err_26 = 1.2533 * sd_gersh_26

med_gersh_err_02 = 1.2533 * sd_gersh_02
med_gersh_err_03 = 1.2533 * sd_gersh_03
med_gersh_err_04 = 1.2533 * sd_gersh_04
med_gersh_err_05 = 1.2533 * sd_gersh_05
med_gersh_err_08 = 1.2533 * sd_gersh_08
med_gersh_err_09 = 1.2533 * sd_gersh_09
med_gersh_err_10 = 1.2533 * sd_gersh_10
med_gersh_err_11 = 1.2533 * sd_gersh_11
med_gersh_err_12 = 1.2533 * sd_gersh_12
med_gersh_err_15 = 1.2533 * sd_gersh_15
med_gersh_err_16 = 1.2533 * sd_gersh_16
med_gersh_err_17 = 1.2533 * sd_gersh_17
med_gersh_err_18 = 1.2533 * sd_gersh_18
med_gersh_err_19 = 1.2533 * sd_gersh_19
med_gersh_err_20 = 1.2533 * sd_gersh_20
med_gersh_err_21 = 1.2533 * sd_gersh_21
med_gersh_err_24 = 1.2533 * sd_gersh_24
med_gersh_err_25 = 1.2533 * sd_gersh_25
med_gersh_err_27 = 1.2533 * sd_gersh_27
med_gersh_err_28 = 1.2533 * sd_gersh_28
med_gersh_err_34 = 1.2533 * sd_gersh_34
med_gersh_err_35 = 1.2533 * sd_gersh_35

med_gersh_errs = np.array([med_gersh_err_01,med_gersh_err_02,med_gersh_err_03,med_gersh_err_04,med_gersh_err_05,med_gersh_err_08,med_gersh_err_09,med_gersh_err_10,med_gersh_err_11,med_gersh_err_12,med_gersh_err_15,med_gersh_err_16,med_gersh_err_17,med_gersh_err_18,med_gersh_err_19,med_gersh_err_20,med_gersh_err_21,med_gersh_err_24,med_gersh_err_25,med_gersh_err_26,med_gersh_err_27,med_gersh_err_28,med_gersh_err_34,med_gersh_err_35])

'''
+ %
'''

gersh_01_mean_p = np.mean(gersh_01_plus)
gersh_12_mean_p = np.mean(gersh_12_plus)
gersh_15_mean_p = np.mean(gersh_15_plus)
gersh_26_mean_p = np.mean(gersh_26_plus)

gersh_02_mean_p = np.mean(gersh_02_plus)
gersh_03_mean_p = np.mean(gersh_03_plus)
gersh_04_mean_p = np.mean(gersh_04_plus)
gersh_05_mean_p = np.mean(gersh_05_plus)
gersh_08_mean_p = np.mean(gersh_08_plus)
gersh_09_mean_p = np.mean(gersh_09_plus)
gersh_10_mean_p = np.mean(gersh_10_plus)
gersh_11_mean_p = np.mean(gersh_11_plus)
gersh_16_mean_p = np.mean(gersh_16_plus)
gersh_17_mean_p = np.mean(gersh_17_plus)
gersh_18_mean_p = np.mean(gersh_18_plus)
gersh_19_mean_p = np.mean(gersh_19_plus)
gersh_20_mean_p = np.mean(gersh_20_plus)
gersh_21_mean_p = np.mean(gersh_21_plus)
gersh_24_mean_p = np.mean(gersh_24_plus)
gersh_25_mean_p = np.mean(gersh_25_plus)
gersh_27_mean_p = np.mean(gersh_27_plus)
gersh_28_mean_p = np.mean(gersh_28_plus)
gersh_34_mean_p = np.mean(gersh_34_plus)
gersh_35_mean_p = np.mean(gersh_35_plus)

gersh_means_p = np.array([gersh_01_mean_p,gersh_02_mean_p,gersh_03_mean_p,gersh_04_mean_p,gersh_05_mean_p,gersh_08_mean_p,gersh_09_mean_p,gersh_10_mean_p,gersh_11_mean_p,gersh_12_mean_p,gersh_15_mean_p,gersh_16_mean_p,gersh_17_mean_p,gersh_18_mean_p,gersh_19_mean_p,gersh_20_mean_p,gersh_21_mean_p,gersh_24_mean_p,gersh_25_mean_p,gersh_26_mean_p,gersh_27_mean_p,gersh_28_mean_p,gersh_34_mean_p,gersh_35_mean_p])

sd_gersh_01_p = standard_deviation(gersh_01_plus,gersh_01_mean_p)
sd_gersh_12_p = standard_deviation(gersh_12_plus,gersh_12_mean_p)
sd_gersh_15_p = standard_deviation(gersh_15_plus,gersh_15_mean_p)
sd_gersh_26_p = standard_deviation(gersh_26_plus,gersh_26_mean_p)

sd_gersh_02_p = standard_deviation(gersh_02_plus,gersh_02_mean_p)
sd_gersh_03_p = standard_deviation(gersh_03_plus,gersh_03_mean_p)
sd_gersh_04_p = standard_deviation(gersh_04_plus,gersh_04_mean_p)
sd_gersh_05_p = standard_deviation(gersh_05_plus,gersh_05_mean_p)
sd_gersh_08_p = standard_deviation(gersh_08_plus,gersh_08_mean_p)
sd_gersh_09_p = standard_deviation(gersh_09_plus,gersh_09_mean_p)
sd_gersh_10_p = standard_deviation(gersh_10_plus,gersh_10_mean_p)
sd_gersh_11_p = standard_deviation(gersh_11_plus,gersh_11_mean_p)
sd_gersh_16_p = standard_deviation(gersh_16_plus,gersh_16_mean_p)
sd_gersh_17_p = standard_deviation(gersh_17_plus,gersh_17_mean_p)
sd_gersh_18_p = standard_deviation(gersh_18_plus,gersh_18_mean_p)
sd_gersh_19_p = standard_deviation(gersh_19_plus,gersh_19_mean_p)
sd_gersh_20_p = standard_deviation(gersh_20_plus,gersh_20_mean_p)
sd_gersh_21_p = standard_deviation(gersh_21_plus,gersh_21_mean_p)
sd_gersh_24_p = standard_deviation(gersh_24_plus,gersh_24_mean_p)
sd_gersh_25_p = standard_deviation(gersh_25_plus,gersh_25_mean_p)
sd_gersh_27_p = standard_deviation(gersh_27_plus,gersh_27_mean_p)
sd_gersh_28_p = standard_deviation(gersh_28_plus,gersh_28_mean_p)
sd_gersh_34_p = standard_deviation(gersh_34_plus,gersh_34_mean_p)
sd_gersh_35_p = standard_deviation(gersh_35_plus,gersh_35_mean_p)


gersh_01_median_p = np.median(gersh_01_plus)
gersh_12_median_p = np.median(gersh_12_plus)
gersh_15_median_p = np.median(gersh_15_plus)
gersh_26_median_p = np.median(gersh_26_plus)

gersh_02_median_p = np.median(gersh_02_plus)
gersh_03_median_p = np.median(gersh_03_plus)
gersh_04_median_p = np.median(gersh_04_plus)
gersh_05_median_p = np.median(gersh_05_plus)
gersh_08_median_p = np.median(gersh_08_plus)
gersh_09_median_p = np.median(gersh_09_plus)
gersh_10_median_p = np.median(gersh_10_plus)
gersh_11_median_p = np.median(gersh_11_plus)
gersh_16_median_p = np.median(gersh_16_plus)
gersh_17_median_p = np.median(gersh_17_plus)
gersh_18_median_p = np.median(gersh_18_plus)
gersh_19_median_p = np.median(gersh_19_plus)
gersh_20_median_p = np.median(gersh_20_plus)
gersh_21_median_p = np.median(gersh_21_plus)
gersh_24_median_p = np.median(gersh_24_plus)
gersh_25_median_p = np.median(gersh_25_plus)
gersh_27_median_p = np.median(gersh_27_plus)
gersh_28_median_p = np.median(gersh_28_plus)
gersh_34_median_p = np.median(gersh_34_plus)
gersh_35_median_p = np.median(gersh_35_plus)

gersh_medians_p = np.array([gersh_01_median_p,gersh_02_median_p,gersh_03_median_p,gersh_04_median_p,gersh_05_median_p,gersh_08_median_p,gersh_09_median_p,gersh_10_median_p,gersh_11_median_p,gersh_12_median_p,gersh_15_median_p,gersh_16_median_p,gersh_17_median_p,gersh_18_median_p,gersh_19_median_p,gersh_20_median_p,gersh_21_median_p,gersh_24_median_p,gersh_25_median_p,gersh_26_median_p,gersh_27_median_p,gersh_28_median_p,gersh_34_median_p,gersh_35_median_p])

med_gersh_err_01_p = 1.2533 * sd_gersh_01_p
med_gersh_err_12_p = 1.2533 * sd_gersh_12_p
med_gersh_err_15_p = 1.2533 * sd_gersh_15_p
med_gersh_err_26_p = 1.2533 * sd_gersh_26_p

med_gersh_err_02_p = 1.2533 * sd_gersh_02_p
med_gersh_err_03_p = 1.2533 * sd_gersh_03_p
med_gersh_err_04_p = 1.2533 * sd_gersh_04_p
med_gersh_err_05_p = 1.2533 * sd_gersh_05_p
med_gersh_err_08_p = 1.2533 * sd_gersh_08_p
med_gersh_err_09_p = 1.2533 * sd_gersh_09_p
med_gersh_err_10_p = 1.2533 * sd_gersh_10_p
med_gersh_err_11_p = 1.2533 * sd_gersh_11_p
med_gersh_err_12_p = 1.2533 * sd_gersh_12_p
med_gersh_err_15_p = 1.2533 * sd_gersh_15_p
med_gersh_err_16_p = 1.2533 * sd_gersh_16_p
med_gersh_err_17_p = 1.2533 * sd_gersh_17_p
med_gersh_err_18_p = 1.2533 * sd_gersh_18_p
med_gersh_err_19_p = 1.2533 * sd_gersh_19_p
med_gersh_err_20_p = 1.2533 * sd_gersh_20_p
med_gersh_err_21_p = 1.2533 * sd_gersh_21_p
med_gersh_err_24_p = 1.2533 * sd_gersh_24_p
med_gersh_err_25_p = 1.2533 * sd_gersh_25_p
med_gersh_err_27_p = 1.2533 * sd_gersh_27_p
med_gersh_err_28_p = 1.2533 * sd_gersh_28_p
med_gersh_err_34_p = 1.2533 * sd_gersh_34_p
med_gersh_err_35_p = 1.2533 * sd_gersh_35_p

med_gersh_errs_p = np.array([med_gersh_err_01_p,med_gersh_err_02_p,med_gersh_err_03_p,med_gersh_err_04_p,med_gersh_err_05_p,med_gersh_err_08_p,med_gersh_err_09_p,med_gersh_err_10_p,med_gersh_err_11_p,med_gersh_err_12_p,med_gersh_err_15_p,med_gersh_err_16_p,med_gersh_err_17_p,med_gersh_err_18_p,med_gersh_err_19_p,med_gersh_err_20_p,med_gersh_err_21_p,med_gersh_err_24_p,med_gersh_err_25_p,med_gersh_err_26_p,med_gersh_err_27_p,med_gersh_err_28_p,med_gersh_err_34_p,med_gersh_err_35_p])

'''
- %
'''

gersh_01_mean_m = np.mean(gersh_01_minus)
gersh_12_mean_m = np.mean(gersh_12_minus)
gersh_15_mean_m = np.mean(gersh_15_minus)
gersh_26_mean_m = np.mean(gersh_26_minus)

gersh_02_mean_m = np.mean(gersh_02_minus)
gersh_03_mean_m = np.mean(gersh_03_minus)
gersh_04_mean_m = np.mean(gersh_04_minus)
gersh_05_mean_m = np.mean(gersh_05_minus)
gersh_08_mean_m = np.mean(gersh_08_minus)
gersh_09_mean_m = np.mean(gersh_09_minus)
gersh_10_mean_m = np.mean(gersh_10_minus)
gersh_11_mean_m = np.mean(gersh_11_minus)
gersh_16_mean_m = np.mean(gersh_16_minus)
gersh_17_mean_m = np.mean(gersh_17_minus)
gersh_18_mean_m = np.mean(gersh_18_minus)
gersh_19_mean_m = np.mean(gersh_19_minus)
gersh_20_mean_m = np.mean(gersh_20_minus)
gersh_21_mean_m = np.mean(gersh_21_minus)
gersh_24_mean_m = np.mean(gersh_24_minus)
gersh_25_mean_m = np.mean(gersh_25_minus)
gersh_27_mean_m = np.mean(gersh_27_minus)
gersh_28_mean_m = np.mean(gersh_28_minus)
gersh_34_mean_m = np.mean(gersh_34_minus)
gersh_35_mean_m = np.mean(gersh_35_minus)

gersh_means_m = np.array([gersh_01_mean_m,gersh_02_mean_m,gersh_03_mean_m,gersh_04_mean_m,gersh_05_mean_m,gersh_08_mean_m,gersh_09_mean_m,gersh_10_mean_m,gersh_11_mean_m,gersh_12_mean_m,gersh_15_mean_m,gersh_16_mean_m,gersh_17_mean_m,gersh_18_mean_m,gersh_19_mean_m,gersh_20_mean_m,gersh_21_mean_m,gersh_24_mean_m,gersh_25_mean_m,gersh_26_mean_m,gersh_27_mean_m,gersh_28_mean_m,gersh_34_mean_m,gersh_35_mean_m])

sd_gersh_01_m = standard_deviation(gersh_01_minus,gersh_01_mean_m)
sd_gersh_12_m = standard_deviation(gersh_12_minus,gersh_12_mean_m)
sd_gersh_15_m = standard_deviation(gersh_15_minus,gersh_15_mean_m)
sd_gersh_26_m = standard_deviation(gersh_26_minus,gersh_26_mean_m)

sd_gersh_02_m = standard_deviation(gersh_02_minus,gersh_02_mean_m)
sd_gersh_03_m = standard_deviation(gersh_03_minus,gersh_03_mean_m)
sd_gersh_04_m = standard_deviation(gersh_04_minus,gersh_04_mean_m)
sd_gersh_05_m = standard_deviation(gersh_05_minus,gersh_05_mean_m)
sd_gersh_08_m = standard_deviation(gersh_08_minus,gersh_08_mean_m)
sd_gersh_09_m = standard_deviation(gersh_09_minus,gersh_09_mean_m)
sd_gersh_10_m = standard_deviation(gersh_10_minus,gersh_10_mean_m)
sd_gersh_11_m = standard_deviation(gersh_11_minus,gersh_11_mean_m)
sd_gersh_16_m = standard_deviation(gersh_16_minus,gersh_16_mean_m)
sd_gersh_17_m = standard_deviation(gersh_17_minus,gersh_17_mean_m)
sd_gersh_18_m = standard_deviation(gersh_18_minus,gersh_18_mean_m)
sd_gersh_19_m = standard_deviation(gersh_19_minus,gersh_19_mean_m)
sd_gersh_20_m = standard_deviation(gersh_20_minus,gersh_20_mean_m)
sd_gersh_21_m = standard_deviation(gersh_21_minus,gersh_21_mean_m)
sd_gersh_24_m = standard_deviation(gersh_24_minus,gersh_24_mean_m)
sd_gersh_25_m = standard_deviation(gersh_25_minus,gersh_25_mean_m)
sd_gersh_27_m = standard_deviation(gersh_27_minus,gersh_27_mean_m)
sd_gersh_28_m = standard_deviation(gersh_28_minus,gersh_28_mean_m)
sd_gersh_34_m = standard_deviation(gersh_34_minus,gersh_34_mean_m)
sd_gersh_35_m = standard_deviation(gersh_35_minus,gersh_35_mean_m)


gersh_01_median_m = np.median(gersh_01_minus)
gersh_12_median_m = np.median(gersh_12_minus)
gersh_15_median_m = np.median(gersh_15_minus)
gersh_26_median_m = np.median(gersh_26_minus)

gersh_02_median_m = np.median(gersh_02_minus)
gersh_03_median_m = np.median(gersh_03_minus)
gersh_04_median_m = np.median(gersh_04_minus)
gersh_05_median_m = np.median(gersh_05_minus)
gersh_08_median_m = np.median(gersh_08_minus)
gersh_09_median_m = np.median(gersh_09_minus)
gersh_10_median_m = np.median(gersh_10_minus)
gersh_11_median_m = np.median(gersh_11_minus)
gersh_16_median_m = np.median(gersh_16_minus)
gersh_17_median_m = np.median(gersh_17_minus)
gersh_18_median_m = np.median(gersh_18_minus)
gersh_19_median_m = np.median(gersh_19_minus)
gersh_20_median_m = np.median(gersh_20_minus)
gersh_21_median_m = np.median(gersh_21_minus)
gersh_24_median_m = np.median(gersh_24_minus)
gersh_25_median_m = np.median(gersh_25_minus)
gersh_27_median_m = np.median(gersh_27_minus)
gersh_28_median_m = np.median(gersh_28_minus)
gersh_34_median_m = np.median(gersh_34_minus)
gersh_35_median_m = np.median(gersh_35_minus)

gersh_medians_m = np.array([gersh_01_median_m,gersh_02_median_m,gersh_03_median_m,gersh_04_median_m,gersh_05_median_m,gersh_08_median_m,gersh_09_median_m,gersh_10_median_m,gersh_11_median_m,gersh_12_median_m,gersh_15_median_m,gersh_16_median_m,gersh_17_median_m,gersh_18_median_m,gersh_19_median_m,gersh_20_median_m,gersh_21_median_m,gersh_24_median_m,gersh_25_median_m,gersh_26_median_m,gersh_27_median_m,gersh_28_median_m,gersh_34_median_m,gersh_35_median_m])

med_gersh_err_01_m = 1.2533 * sd_gersh_01_m
med_gersh_err_12_m = 1.2533 * sd_gersh_12_m
med_gersh_err_15_m = 1.2533 * sd_gersh_15_m
med_gersh_err_26_m = 1.2533 * sd_gersh_26_m

med_gersh_err_02_m = 1.2533 * sd_gersh_02_m
med_gersh_err_03_m = 1.2533 * sd_gersh_03_m
med_gersh_err_04_m = 1.2533 * sd_gersh_04_m
med_gersh_err_05_m = 1.2533 * sd_gersh_05_m
med_gersh_err_08_m = 1.2533 * sd_gersh_08_m
med_gersh_err_09_m = 1.2533 * sd_gersh_09_m
med_gersh_err_10_m = 1.2533 * sd_gersh_10_m
med_gersh_err_11_m = 1.2533 * sd_gersh_11_m
med_gersh_err_12_m = 1.2533 * sd_gersh_12_m
med_gersh_err_15_m = 1.2533 * sd_gersh_15_m
med_gersh_err_16_m = 1.2533 * sd_gersh_16_m
med_gersh_err_17_m = 1.2533 * sd_gersh_17_m
med_gersh_err_18_m = 1.2533 * sd_gersh_18_m
med_gersh_err_19_m = 1.2533 * sd_gersh_19_m
med_gersh_err_20_m = 1.2533 * sd_gersh_20_m
med_gersh_err_21_m = 1.2533 * sd_gersh_21_m
med_gersh_err_24_m = 1.2533 * sd_gersh_24_m
med_gersh_err_25_m = 1.2533 * sd_gersh_25_m
med_gersh_err_27_m = 1.2533 * sd_gersh_27_m
med_gersh_err_28_m = 1.2533 * sd_gersh_28_m
med_gersh_err_34_m = 1.2533 * sd_gersh_34_m
med_gersh_err_35_m = 1.2533 * sd_gersh_35_m

med_gersh_errs_m = np.array([med_gersh_err_01_m,med_gersh_err_02_m,med_gersh_err_03_m,med_gersh_err_04_m,med_gersh_err_05_m,med_gersh_err_08_m,med_gersh_err_09_m,med_gersh_err_10_m,med_gersh_err_11_m,med_gersh_err_12_m,med_gersh_err_15_m,med_gersh_err_16_m,med_gersh_err_17_m,med_gersh_err_18_m,med_gersh_err_19_m,med_gersh_err_20_m,med_gersh_err_21_m,med_gersh_err_24_m,med_gersh_err_25_m,med_gersh_err_26_m,med_gersh_err_27_m,med_gersh_err_28_m,med_gersh_err_34_m,med_gersh_err_35_m])


'''
AVERAGE KH
'''
KH_02_mean = np.mean(KH_02)
KH_03_mean = np.mean(KH_03)
KH_04_mean = np.mean(KH_04)
KH_05_mean = np.mean(KH_05)
KH_08_mean = np.mean(KH_08)
KH_09_mean = np.mean(KH_09)
KH_10_mean = np.mean(KH_10)
KH_11_mean = np.mean(KH_11)
KH_16_mean = np.mean(KH_16)
KH_17_mean = np.mean(KH_17)
KH_18_mean = np.mean(KH_18)
KH_19_mean = np.mean(KH_19)
KH_20_mean = np.mean(KH_20)
KH_21_mean = np.mean(KH_21)
KH_24_mean = np.mean(KH_24)
KH_25_mean = np.mean(KH_25)
KH_27_mean = np.mean(KH_27)
KH_28_mean = np.mean(KH_28)
KH_34_mean = np.mean(KH_34)
KH_35_mean = np.mean(KH_35)

KH_01_mean = np.mean(KH_01)
KH_12_mean = np.mean(KH_12)
KH_15_mean = np.mean(KH_15)
KH_26_mean = np.mean(KH_26)

KH_means = np.array([KH_01_mean,KH_02_mean,KH_03_mean,KH_04_mean,KH_05_mean,KH_08_mean,KH_09_mean,KH_10_mean,KH_11_mean,KH_12_mean,KH_15_mean,KH_16_mean,KH_17_mean,KH_18_mean,KH_19_mean,KH_20_mean,KH_21_mean,KH_24_mean,KH_25_mean,KH_26_mean,KH_27_mean,KH_28_mean,KH_34_mean,KH_35_mean])

sd_KH_02 = standard_deviation(KH_02,KH_02_mean)
sd_KH_03 = standard_deviation(KH_03,KH_03_mean)
sd_KH_04 = standard_deviation(KH_04,KH_04_mean)
sd_KH_05 = standard_deviation(KH_05,KH_05_mean)
sd_KH_08 = standard_deviation(KH_08,KH_08_mean)
sd_KH_09 = standard_deviation(KH_09,KH_09_mean)
sd_KH_10 = standard_deviation(KH_10,KH_10_mean)
sd_KH_11 = standard_deviation(KH_11,KH_11_mean)
sd_KH_16 = standard_deviation(KH_16,KH_16_mean)
sd_KH_17 = standard_deviation(KH_17,KH_17_mean)
sd_KH_18 = standard_deviation(KH_18,KH_18_mean)
sd_KH_19 = standard_deviation(KH_19,KH_19_mean)
sd_KH_20 = standard_deviation(KH_20,KH_20_mean)
sd_KH_21 = standard_deviation(KH_21,KH_21_mean)
sd_KH_24 = standard_deviation(KH_24,KH_24_mean)
sd_KH_25 = standard_deviation(KH_25,KH_25_mean)
sd_KH_27 = standard_deviation(KH_27,KH_27_mean)
sd_KH_28 = standard_deviation(KH_28,KH_28_mean)
sd_KH_34 = standard_deviation(KH_34,KH_34_mean)
sd_KH_35 = standard_deviation(KH_35,KH_35_mean)

sd_KH_01 = standard_deviation(KH_01,KH_01_mean)
sd_KH_12 = standard_deviation(KH_12,KH_12_mean)
sd_KH_15 = standard_deviation(KH_15,KH_15_mean)
sd_KH_26 = standard_deviation(KH_26,KH_26_mean)


KH_02_median = np.median(KH_02)
KH_03_median = np.median(KH_03)
KH_04_median = np.median(KH_04)
KH_05_median = np.median(KH_05)
KH_08_median = np.median(KH_08)
KH_09_median = np.median(KH_09)
KH_10_median = np.median(KH_10)
KH_11_median = np.median(KH_11)
KH_16_median = np.median(KH_16)
KH_17_median = np.median(KH_17)
KH_18_median = np.median(KH_18)
KH_19_median = np.median(KH_19)
KH_20_median = np.median(KH_20)
KH_21_median = np.median(KH_21)
KH_24_median = np.median(KH_24)
KH_25_median = np.median(KH_25)
KH_27_median = np.median(KH_27)
KH_28_median = np.median(KH_28)
KH_34_median = np.median(KH_34)
KH_35_median = np.median(KH_35)

KH_01_median = np.median(KH_01)
KH_12_median = np.median(KH_12)
KH_15_median = np.median(KH_15)
KH_26_median = np.median(KH_26)

KH_medians = np.array([KH_01_median,KH_02_median,KH_03_median,KH_04_median,KH_05_median,KH_08_median,KH_09_median,KH_10_median,KH_11_median,KH_12_median,KH_15_median,KH_16_median,KH_17_median,KH_18_median,KH_19_median,KH_20_median,KH_21_median,KH_24_median,KH_25_median,KH_26_median,KH_27_median,KH_28_median,KH_34_median,KH_35_median])

med_KH_02_err = 1.2533 * sd_KH_02
med_KH_03_err = 1.2533 * sd_KH_03
med_KH_04_err = 1.2533 * sd_KH_04
med_KH_05_err = 1.2533 * sd_KH_05
med_KH_08_err = 1.2533 * sd_KH_08
med_KH_09_err = 1.2533 * sd_KH_09
med_KH_10_err = 1.2533 * sd_KH_10
med_KH_11_err = 1.2533 * sd_KH_11
med_KH_16_err = 1.2533 * sd_KH_16
med_KH_17_err = 1.2533 * sd_KH_17
med_KH_18_err = 1.2533 * sd_KH_18
med_KH_19_err = 1.2533 * sd_KH_19
med_KH_20_err = 1.2533 * sd_KH_20
med_KH_21_err = 1.2533 * sd_KH_21
med_KH_24_err = 1.2533 * sd_KH_24
med_KH_25_err = 1.2533 * sd_KH_25
med_KH_27_err = 1.2533 * sd_KH_27
med_KH_28_err = 1.2533 * sd_KH_28
med_KH_34_err = 1.2533 * sd_KH_34
med_KH_35_err = 1.2533 * sd_KH_35

med_KH_01_err = 1.2533 * sd_KH_01
med_KH_12_err = 1.2533 * sd_KH_12
med_KH_15_err = 1.2533 * sd_KH_15
med_KH_26_err = 1.2533 * sd_KH_26

med_KH_errs = np.array([med_KH_01_err,med_KH_02_err,med_KH_03_err,med_KH_04_err,med_KH_05_err,med_KH_08_err,med_KH_09_err,med_KH_10_err,med_KH_11_err,med_KH_12_err,med_KH_15_err,med_KH_16_err,med_KH_17_err,med_KH_18_err,med_KH_19_err,med_KH_20_err,med_KH_21_err,med_KH_24_err,med_KH_25_err,med_KH_26_err,med_KH_27_err,med_KH_28_err,med_KH_34_err,med_KH_35_err])

'''
+
'''

KH_02_mean_plus = np.mean(KH_02_plus)
KH_03_mean_plus = np.mean(KH_03_plus)
KH_04_mean_plus = np.mean(KH_04_plus)
KH_05_mean_plus = np.mean(KH_05_plus)
KH_08_mean_plus = np.mean(KH_08_plus)
KH_09_mean_plus = np.mean(KH_09_plus)
KH_10_mean_plus = np.mean(KH_10_plus)
KH_11_mean_plus = np.mean(KH_11_plus)
KH_16_mean_plus = np.mean(KH_16_plus)
KH_17_mean_plus = np.mean(KH_17_plus)
KH_18_mean_plus = np.mean(KH_18_plus)
KH_19_mean_plus = np.mean(KH_19_plus)
KH_20_mean_plus = np.mean(KH_20_plus)
KH_21_mean_plus = np.mean(KH_21_plus)
KH_24_mean_plus = np.mean(KH_24_plus)
KH_25_mean_plus = np.mean(KH_25_plus)
KH_27_mean_plus = np.mean(KH_27_plus)
KH_28_mean_plus = np.mean(KH_28_plus)
KH_34_mean_plus = np.mean(KH_34_plus)
KH_35_mean_plus = np.mean(KH_35_plus)

KH_01_mean_plus = np.mean(KH_01_plus)
KH_12_mean_plus = np.mean(KH_12_plus)
KH_15_mean_plus = np.mean(KH_15_plus)
KH_26_mean_plus = np.mean(KH_26_plus)

KH_means_plus = np.array([KH_01_mean_plus,KH_02_mean_plus,KH_03_mean_plus,KH_04_mean_plus,KH_05_mean_plus,KH_08_mean_plus,KH_09_mean_plus,KH_10_mean_plus,KH_11_mean,KH_12_mean_plus,KH_15_mean_plus,KH_16_mean_plus,KH_17_mean_plus,KH_18_mean_plus,KH_19_mean_plus,KH_20_mean_plus,KH_21_mean_plus,KH_24_mean_plus,KH_25_mean_plus,KH_26_mean_plus,KH_27_mean_plus,KH_28_mean_plus,KH_34_mean_plus,KH_35_mean_plus])

sd_KH_02_plus = standard_deviation(KH_02_plus,KH_02_mean_plus)
sd_KH_03_plus = standard_deviation(KH_03_plus,KH_03_mean_plus)
sd_KH_04_plus = standard_deviation(KH_04_plus,KH_04_mean_plus)
sd_KH_05_plus = standard_deviation(KH_05_plus,KH_05_mean_plus)
sd_KH_08_plus = standard_deviation(KH_08_plus,KH_08_mean_plus)
sd_KH_09_plus = standard_deviation(KH_09_plus,KH_09_mean_plus)
sd_KH_10_plus = standard_deviation(KH_10_plus,KH_10_mean_plus)
sd_KH_11_plus = standard_deviation(KH_11_plus,KH_11_mean_plus)
sd_KH_16_plus = standard_deviation(KH_16_plus,KH_16_mean_plus)
sd_KH_17_plus = standard_deviation(KH_17_plus,KH_17_mean_plus)
sd_KH_18_plus = standard_deviation(KH_18_plus,KH_18_mean_plus)
sd_KH_19_plus = standard_deviation(KH_19_plus,KH_19_mean_plus)
sd_KH_20_plus = standard_deviation(KH_20_plus,KH_20_mean_plus)
sd_KH_21_plus = standard_deviation(KH_21_plus,KH_21_mean_plus)
sd_KH_24_plus = standard_deviation(KH_24_plus,KH_24_mean_plus)
sd_KH_25_plus = standard_deviation(KH_25_plus,KH_25_mean_plus)
sd_KH_27_plus = standard_deviation(KH_27_plus,KH_27_mean_plus)
sd_KH_28_plus = standard_deviation(KH_28_plus,KH_28_mean_plus)
sd_KH_34_plus = standard_deviation(KH_34_plus,KH_34_mean_plus)
sd_KH_35_plus = standard_deviation(KH_35_plus,KH_35_mean_plus)

sd_KH_01_plus = standard_deviation(KH_01_plus,KH_01_mean_plus)
sd_KH_12_plus = standard_deviation(KH_12_plus,KH_12_mean_plus)
sd_KH_15_plus = standard_deviation(KH_15_plus,KH_15_mean_plus)
sd_KH_26_plus = standard_deviation(KH_26_plus,KH_26_mean_plus)


KH_02_median_plus = np.median(KH_02_plus)
KH_03_median_plus = np.median(KH_03_plus)
KH_04_median_plus = np.median(KH_04_plus)
KH_05_median_plus = np.median(KH_05_plus)
KH_08_median_plus = np.median(KH_08_plus)
KH_09_median_plus = np.median(KH_09_plus)
KH_10_median_plus = np.median(KH_10_plus)
KH_11_median_plus = np.median(KH_11_plus)
KH_16_median_plus = np.median(KH_16_plus)
KH_17_median_plus = np.median(KH_17_plus)
KH_18_median_plus = np.median(KH_18_plus)
KH_19_median_plus = np.median(KH_19_plus)
KH_20_median_plus = np.median(KH_20_plus)
KH_21_median_plus = np.median(KH_21_plus)
KH_24_median_plus = np.median(KH_24_plus)
KH_25_median_plus = np.median(KH_25_plus)
KH_27_median_plus = np.median(KH_27_plus)
KH_28_median_plus = np.median(KH_28_plus)
KH_34_median_plus = np.median(KH_34_plus)
KH_35_median_plus = np.median(KH_35_plus)

KH_01_median_plus = np.median(KH_01_plus)
KH_12_median_plus = np.median(KH_12_plus)
KH_15_median_plus = np.median(KH_15_plus)
KH_26_median_plus = np.median(KH_26_plus)

KH_medians_plus = np.array([KH_01_median_plus,KH_02_median_plus,KH_03_median_plus,KH_04_median_plus,KH_05_median_plus,KH_08_median_plus,KH_09_median_plus,KH_10_median_plus,KH_11_median_plus,KH_12_median_plus,KH_15_median_plus,KH_16_median_plus,KH_17_median_plus,KH_18_median_plus,KH_19_median_plus,KH_20_median_plus,KH_21_median_plus,KH_24_median_plus,KH_25_median_plus,KH_26_median_plus,KH_27_median_plus,KH_28_median_plus,KH_34_median_plus,KH_35_median_plus])

med_KH_02_err_plus = 1.2533 * sd_KH_02_plus
med_KH_03_err_plus = 1.2533 * sd_KH_03_plus
med_KH_04_err_plus = 1.2533 * sd_KH_04_plus
med_KH_05_err_plus = 1.2533 * sd_KH_05_plus
med_KH_08_err_plus = 1.2533 * sd_KH_08_plus
med_KH_09_err_plus = 1.2533 * sd_KH_09_plus
med_KH_10_err_plus = 1.2533 * sd_KH_10_plus
med_KH_11_err_plus = 1.2533 * sd_KH_11_plus
med_KH_16_err_plus = 1.2533 * sd_KH_16_plus
med_KH_17_err_plus = 1.2533 * sd_KH_17_plus
med_KH_18_err_plus = 1.2533 * sd_KH_18_plus
med_KH_19_err_plus = 1.2533 * sd_KH_19_plus
med_KH_20_err_plus = 1.2533 * sd_KH_20_plus
med_KH_21_err_plus = 1.2533 * sd_KH_21_plus
med_KH_24_err_plus = 1.2533 * sd_KH_24_plus
med_KH_25_err_plus = 1.2533 * sd_KH_25_plus
med_KH_27_err_plus = 1.2533 * sd_KH_27_plus
med_KH_28_err_plus = 1.2533 * sd_KH_28_plus
med_KH_34_err_plus = 1.2533 * sd_KH_34_plus
med_KH_35_err_plus = 1.2533 * sd_KH_35_plus

med_KH_01_err_plus = 1.2533 * sd_KH_01_plus
med_KH_12_err_plus = 1.2533 * sd_KH_12_plus
med_KH_15_err_plus = 1.2533 * sd_KH_15_plus
med_KH_26_err_plus = 1.2533 * sd_KH_26_plus

med_KH_errs_plus = np.array([med_KH_01_err_plus,med_KH_02_err_plus,med_KH_03_err_plus,med_KH_04_err_plus,med_KH_05_err_plus,med_KH_08_err_plus,med_KH_09_err_plus,med_KH_10_err_plus,med_KH_11_err_plus,med_KH_12_err_plus,med_KH_15_err_plus,med_KH_16_err_plus,med_KH_17_err_plus,med_KH_18_err_plus,med_KH_19_err_plus,med_KH_20_err_plus,med_KH_21_err_plus,med_KH_24_err_plus,med_KH_25_err_plus,med_KH_26_err_plus,med_KH_27_err_plus,med_KH_28_err_plus,med_KH_34_err_plus,med_KH_35_err_plus])

'''
-
'''
KH_02_mean_minus = np.mean(KH_02_minus)
KH_03_mean_minus = np.mean(KH_03_minus)
KH_04_mean_minus = np.mean(KH_04_minus)
KH_05_mean_minus = np.mean(KH_05_minus)
KH_08_mean_minus = np.mean(KH_08_minus)
KH_09_mean_minus = np.mean(KH_09_minus)
KH_10_mean_minus = np.mean(KH_10_minus)
KH_11_mean_minus = np.mean(KH_11_minus)
KH_16_mean_minus = np.mean(KH_16_minus)
KH_17_mean_minus = np.mean(KH_17_minus)
KH_18_mean_minus = np.mean(KH_18_minus)
KH_19_mean_minus = np.mean(KH_19_minus)
KH_20_mean_minus = np.mean(KH_20_minus)
KH_21_mean_minus = np.mean(KH_21_minus)
KH_24_mean_minus = np.mean(KH_24_minus)
KH_25_mean_minus = np.mean(KH_25_minus)
KH_27_mean_minus = np.mean(KH_27_minus)
KH_28_mean_minus = np.mean(KH_28_minus)
KH_34_mean_minus = np.mean(KH_34_minus)
KH_35_mean_minus = np.mean(KH_35_minus)

KH_01_mean_minus = np.mean(KH_01_minus)
KH_12_mean_minus = np.mean(KH_12_minus)
KH_15_mean_minus = np.mean(KH_15_minus)
KH_26_mean_minus = np.mean(KH_26_minus)

KH_means_minus = np.array([KH_01_mean_minus,KH_02_mean_minus,KH_03_mean_minus,KH_04_mean_minus,KH_05_mean_minus,KH_08_mean_minus,KH_09_mean_minus,KH_10_mean_minus,KH_11_mean_minus,KH_12_mean_minus,KH_15_mean_minus,KH_16_mean_minus,KH_17_mean_minus,KH_18_mean_minus,KH_19_mean_minus,KH_20_mean_minus,KH_21_mean_minus,KH_24_mean_minus,KH_25_mean_minus,KH_26_mean_minus,KH_27_mean_minus,KH_28_mean_minus,KH_34_mean_minus,KH_35_mean_minus])

sd_KH_02_minus = standard_deviation(KH_02_minus,KH_02_mean_minus)
sd_KH_03_minus = standard_deviation(KH_03_minus,KH_03_mean_minus)
sd_KH_04_minus = standard_deviation(KH_04_minus,KH_04_mean_minus)
sd_KH_05_minus = standard_deviation(KH_05_minus,KH_05_mean_minus)
sd_KH_08_minus = standard_deviation(KH_08_minus,KH_08_mean_minus)
sd_KH_09_minus = standard_deviation(KH_09_minus,KH_09_mean_minus)
sd_KH_10_minus = standard_deviation(KH_10_minus,KH_10_mean_minus)
sd_KH_11_minus = standard_deviation(KH_11_minus,KH_11_mean_minus)
sd_KH_16_minus = standard_deviation(KH_16_minus,KH_16_mean_minus)
sd_KH_17_minus = standard_deviation(KH_17_minus,KH_17_mean_minus)
sd_KH_18_minus = standard_deviation(KH_18_minus,KH_18_mean_minus)
sd_KH_19_minus = standard_deviation(KH_19_minus,KH_19_mean_minus)
sd_KH_20_minus = standard_deviation(KH_20_minus,KH_20_mean_minus)
sd_KH_21_minus = standard_deviation(KH_21_minus,KH_21_mean_minus)
sd_KH_24_minus = standard_deviation(KH_24_minus,KH_24_mean_minus)
sd_KH_25_minus = standard_deviation(KH_25_minus,KH_25_mean_minus)
sd_KH_27_minus = standard_deviation(KH_27_minus,KH_27_mean_minus)
sd_KH_28_minus = standard_deviation(KH_28_minus,KH_28_mean_minus)
sd_KH_34_minus = standard_deviation(KH_34_minus,KH_34_mean_minus)
sd_KH_35_minus = standard_deviation(KH_35_minus,KH_35_mean_minus)

sd_KH_01_minus = standard_deviation(KH_01_minus,KH_01_mean_minus)
sd_KH_12_minus = standard_deviation(KH_12_minus,KH_12_mean_minus)
sd_KH_15_minus = standard_deviation(KH_15_minus,KH_15_mean_minus)
sd_KH_26_minus = standard_deviation(KH_26_minus,KH_26_mean_minus)


KH_02_median_minus = np.median(KH_02_minus)
KH_03_median_minus = np.median(KH_03_minus)
KH_04_median_minus = np.median(KH_04_minus)
KH_05_median_minus = np.median(KH_05_minus)
KH_08_median_minus = np.median(KH_08_minus)
KH_09_median_minus = np.median(KH_09_minus)
KH_10_median_minus = np.median(KH_10_minus)
KH_11_median_minus = np.median(KH_11_minus)
KH_16_median_minus = np.median(KH_16_minus)
KH_17_median_minus = np.median(KH_17_minus)
KH_18_median_minus = np.median(KH_18_minus)
KH_19_median_minus = np.median(KH_19_minus)
KH_20_median_minus = np.median(KH_20_minus)
KH_21_median_minus = np.median(KH_21_minus)
KH_24_median_minus = np.median(KH_24_minus)
KH_25_median_minus = np.median(KH_25_minus)
KH_27_median_minus = np.median(KH_27_minus)
KH_28_median_minus = np.median(KH_28_minus)
KH_34_median_minus = np.median(KH_34_minus)
KH_35_median_minus = np.median(KH_35_minus)

KH_01_median_minus = np.median(KH_01_minus)
KH_12_median_minus = np.median(KH_12_minus)
KH_15_median_minus = np.median(KH_15_minus)
KH_26_median_minus = np.median(KH_26_minus)

KH_medians_minus = np.array([KH_01_median_minus,KH_02_median_minus,KH_03_median_minus,KH_04_median_minus,KH_05_median_minus,KH_08_median_minus,KH_09_median_minus,KH_10_median_minus,KH_11_median_minus,KH_12_median_minus,KH_15_median_minus,KH_16_median_minus,KH_17_median_minus,KH_18_median_minus,KH_19_median_minus,KH_20_median_minus,KH_21_median_minus,KH_24_median_minus,KH_25_median_minus,KH_26_median_minus,KH_27_median_minus,KH_28_median_minus,KH_34_median_minus,KH_35_median_minus])

med_KH_02_err_minus = 1.2533 * sd_KH_02_minus
med_KH_03_err_minus = 1.2533 * sd_KH_03_minus
med_KH_04_err_minus = 1.2533 * sd_KH_04_minus
med_KH_05_err_minus = 1.2533 * sd_KH_05_minus
med_KH_08_err_minus = 1.2533 * sd_KH_08_minus
med_KH_09_err_minus = 1.2533 * sd_KH_09_minus
med_KH_10_err_minus = 1.2533 * sd_KH_10_minus
med_KH_11_err_minus = 1.2533 * sd_KH_11_minus
med_KH_16_err_minus = 1.2533 * sd_KH_16_minus
med_KH_17_err_minus = 1.2533 * sd_KH_17_minus
med_KH_18_err_minus = 1.2533 * sd_KH_18_minus
med_KH_19_err_minus = 1.2533 * sd_KH_19_minus
med_KH_20_err_minus = 1.2533 * sd_KH_20_minus
med_KH_21_err_minus = 1.2533 * sd_KH_21_minus
med_KH_24_err_minus = 1.2533 * sd_KH_24_minus
med_KH_25_err_minus = 1.2533 * sd_KH_25_minus
med_KH_27_err_minus = 1.2533 * sd_KH_27_minus
med_KH_28_err_minus = 1.2533 * sd_KH_28_minus
med_KH_34_err_minus = 1.2533 * sd_KH_34_minus
med_KH_35_err_minus = 1.2533 * sd_KH_35_minus

med_KH_01_err_minus = 1.2533 * sd_KH_01_minus
med_KH_12_err_minus = 1.2533 * sd_KH_12_minus
med_KH_15_err_minus = 1.2533 * sd_KH_15_minus
med_KH_26_err_minus = 1.2533 * sd_KH_26_minus

med_KH_errs_minus = np.array([med_KH_01_err_minus,med_KH_02_err_minus,med_KH_03_err_minus,med_KH_04_err_minus,med_KH_05_err_minus,med_KH_08_err_minus,med_KH_09_err_minus,med_KH_10_err_minus,med_KH_11_err_minus,med_KH_12_err_minus,med_KH_15_err_minus,med_KH_16_err_minus,med_KH_17_err_minus,med_KH_18_err_minus,med_KH_19_err_minus,med_KH_20_err_minus,med_KH_21_err_minus,med_KH_24_err_minus,med_KH_25_err_minus,med_KH_26_err_minus,med_KH_27_err_minus,med_KH_28_err_minus,med_KH_34_err_minus,med_KH_35_err_minus])


'''
AVERAGE KH - DAWN
'''
KH_dawn_02_mean = np.mean(KH_dawn_02)
KH_dawn_03_mean = np.mean(KH_dawn_03)
KH_dawn_04_mean = np.mean(KH_dawn_04)
KH_dawn_05_mean = np.mean(KH_dawn_05)
KH_dawn_08_mean = np.mean(KH_dawn_08)
KH_dawn_09_mean = np.mean(KH_dawn_09)
KH_dawn_10_mean = np.mean(KH_dawn_10)
KH_dawn_11_mean = np.mean(KH_dawn_11)
KH_dawn_16_mean = np.mean(KH_dawn_16)
KH_dawn_17_mean = np.mean(KH_dawn_17)
KH_dawn_18_mean = np.mean(KH_dawn_18)
KH_dawn_19_mean = np.mean(KH_dawn_19)
KH_dawn_20_mean = np.mean(KH_dawn_20)
KH_dawn_21_mean = np.mean(KH_dawn_21)
KH_dawn_24_mean = np.mean(KH_dawn_24)
KH_dawn_25_mean = np.mean(KH_dawn_25)
KH_dawn_27_mean = np.mean(KH_dawn_27)
KH_dawn_28_mean = np.mean(KH_dawn_28)
KH_dawn_34_mean = np.mean(KH_dawn_34)
KH_dawn_35_mean = np.mean(KH_dawn_35)

KH_dawn_01_mean = np.mean(KH_dawn_01)
KH_dawn_12_mean = np.mean(KH_dawn_12)
KH_dawn_15_mean = np.mean(KH_dawn_15)
KH_dawn_26_mean = np.mean(KH_dawn_26)

KH_dawn_means = np.array([KH_dawn_01_mean,KH_dawn_02_mean,KH_dawn_03_mean,KH_dawn_04_mean,KH_dawn_05_mean,KH_dawn_08_mean,KH_dawn_09_mean,KH_dawn_10_mean,KH_dawn_11_mean,KH_dawn_12_mean,KH_dawn_15_mean,KH_dawn_16_mean,KH_dawn_17_mean,KH_dawn_18_mean,KH_dawn_19_mean,KH_dawn_20_mean,KH_dawn_21_mean,KH_dawn_24_mean,KH_dawn_25_mean,KH_dawn_26_mean,KH_dawn_27_mean,KH_dawn_28_mean,KH_dawn_34_mean,KH_dawn_35_mean])

sd_KH_dawn_02 = standard_deviation(KH_dawn_02,KH_dawn_02_mean)
sd_KH_dawn_03 = standard_deviation(KH_dawn_03,KH_dawn_03_mean)
sd_KH_dawn_04 = standard_deviation(KH_dawn_04,KH_dawn_04_mean)
sd_KH_dawn_05 = standard_deviation(KH_dawn_05,KH_dawn_05_mean)
sd_KH_dawn_08 = standard_deviation(KH_dawn_08,KH_dawn_08_mean)
sd_KH_dawn_09 = standard_deviation(KH_dawn_09,KH_dawn_09_mean)
sd_KH_dawn_10 = standard_deviation(KH_dawn_10,KH_dawn_10_mean)
sd_KH_dawn_11 = standard_deviation(KH_dawn_11,KH_dawn_11_mean)
sd_KH_dawn_16 = standard_deviation(KH_dawn_16,KH_dawn_16_mean)
sd_KH_dawn_17 = standard_deviation(KH_dawn_17,KH_dawn_17_mean)
sd_KH_dawn_18 = standard_deviation(KH_dawn_18,KH_dawn_18_mean)
sd_KH_dawn_19 = standard_deviation(KH_dawn_19,KH_dawn_19_mean)
sd_KH_dawn_20 = standard_deviation(KH_dawn_20,KH_dawn_20_mean)
sd_KH_dawn_21 = standard_deviation(KH_dawn_21,KH_dawn_21_mean)
sd_KH_dawn_24 = standard_deviation(KH_dawn_24,KH_dawn_24_mean)
sd_KH_dawn_25 = standard_deviation(KH_dawn_25,KH_dawn_25_mean)
sd_KH_dawn_27 = standard_deviation(KH_dawn_27,KH_dawn_27_mean)
sd_KH_dawn_28 = standard_deviation(KH_dawn_28,KH_dawn_28_mean)
sd_KH_dawn_34 = standard_deviation(KH_dawn_34,KH_dawn_34_mean)
sd_KH_dawn_35 = standard_deviation(KH_dawn_35,KH_dawn_35_mean)

sd_KH_dawn_01 = standard_deviation(KH_dawn_01,KH_dawn_01_mean)
sd_KH_dawn_12 = standard_deviation(KH_dawn_12,KH_dawn_12_mean)
sd_KH_dawn_15 = standard_deviation(KH_dawn_15,KH_dawn_15_mean)
sd_KH_dawn_26 = standard_deviation(KH_dawn_26,KH_dawn_26_mean)


KH_dawn_02_median = np.median(KH_dawn_02)
KH_dawn_03_median = np.median(KH_dawn_03)
KH_dawn_04_median = np.median(KH_dawn_04)
KH_dawn_05_median = np.median(KH_dawn_05)
KH_dawn_08_median = np.median(KH_dawn_08)
KH_dawn_09_median = np.median(KH_dawn_09)
KH_dawn_10_median = np.median(KH_dawn_10)
KH_dawn_11_median = np.median(KH_dawn_11)
KH_dawn_16_median = np.median(KH_dawn_16)
KH_dawn_17_median = np.median(KH_dawn_17)
KH_dawn_18_median = np.median(KH_dawn_18)
KH_dawn_19_median = np.median(KH_dawn_19)
KH_dawn_20_median = np.median(KH_dawn_20)
KH_dawn_21_median = np.median(KH_dawn_21)
KH_dawn_24_median = np.median(KH_dawn_24)
KH_dawn_25_median = np.median(KH_dawn_25)
KH_dawn_27_median = np.median(KH_dawn_27)
KH_dawn_28_median = np.median(KH_dawn_28)
KH_dawn_34_median = np.median(KH_dawn_34)
KH_dawn_35_median = np.median(KH_dawn_35)

KH_dawn_01_median = np.median(KH_dawn_01)
KH_dawn_12_median = np.median(KH_dawn_12)
KH_dawn_15_median = np.median(KH_dawn_15)
KH_dawn_26_median = np.median(KH_dawn_26)

KH_dawn_medians = np.array([KH_dawn_01_median,KH_dawn_02_median,KH_dawn_03_median,KH_dawn_04_median,KH_dawn_05_median,KH_dawn_08_median,KH_dawn_09_median,KH_dawn_10_median,KH_dawn_11_median,KH_dawn_12_median,KH_dawn_15_median,KH_dawn_16_median,KH_dawn_17_median,KH_dawn_18_median,KH_dawn_19_median,KH_dawn_20_median,KH_dawn_21_median,KH_dawn_24_median,KH_dawn_25_median,KH_dawn_26_median,KH_dawn_27_median,KH_dawn_28_median,KH_dawn_34_median,KH_dawn_35_median])

med_KH_dawn_02_err = 1.2533 * sd_KH_dawn_02
med_KH_dawn_03_err = 1.2533 * sd_KH_dawn_03
med_KH_dawn_04_err = 1.2533 * sd_KH_dawn_04
med_KH_dawn_05_err = 1.2533 * sd_KH_dawn_05
med_KH_dawn_08_err = 1.2533 * sd_KH_dawn_08
med_KH_dawn_09_err = 1.2533 * sd_KH_dawn_09
med_KH_dawn_10_err = 1.2533 * sd_KH_dawn_10
med_KH_dawn_11_err = 1.2533 * sd_KH_dawn_11
med_KH_dawn_16_err = 1.2533 * sd_KH_dawn_16
med_KH_dawn_17_err = 1.2533 * sd_KH_dawn_17
med_KH_dawn_18_err = 1.2533 * sd_KH_dawn_18
med_KH_dawn_19_err = 1.2533 * sd_KH_dawn_19
med_KH_dawn_20_err = 1.2533 * sd_KH_dawn_20
med_KH_dawn_21_err = 1.2533 * sd_KH_dawn_21
med_KH_dawn_24_err = 1.2533 * sd_KH_dawn_24
med_KH_dawn_25_err = 1.2533 * sd_KH_dawn_25
med_KH_dawn_27_err = 1.2533 * sd_KH_dawn_27
med_KH_dawn_28_err = 1.2533 * sd_KH_dawn_28
med_KH_dawn_34_err = 1.2533 * sd_KH_dawn_34
med_KH_dawn_35_err = 1.2533 * sd_KH_dawn_35

med_KH_dawn_01_err = 1.2533 * sd_KH_dawn_01
med_KH_dawn_12_err = 1.2533 * sd_KH_dawn_12
med_KH_dawn_15_err = 1.2533 * sd_KH_dawn_15
med_KH_dawn_26_err = 1.2533 * sd_KH_dawn_26

med_KH_dawn_errs = np.array([med_KH_dawn_01_err,med_KH_dawn_02_err,med_KH_dawn_03_err,med_KH_dawn_04_err,med_KH_dawn_05_err,med_KH_dawn_08_err,med_KH_dawn_09_err,med_KH_dawn_10_err,med_KH_dawn_11_err,med_KH_dawn_12_err,med_KH_dawn_15_err,med_KH_dawn_16_err,med_KH_dawn_17_err,med_KH_dawn_18_err,med_KH_dawn_19_err,med_KH_dawn_20_err,med_KH_dawn_21_err,med_KH_dawn_24_err,med_KH_dawn_25_err,med_KH_dawn_26_err,med_KH_dawn_27_err,med_KH_dawn_28_err,med_KH_dawn_34_err,med_KH_dawn_35_err])

'''
+
'''

KH_dawn_02_mean_plus = np.mean(KH_dawn_02_plus)
KH_dawn_03_mean_plus = np.mean(KH_dawn_03_plus)
KH_dawn_04_mean_plus = np.mean(KH_dawn_04_plus)
KH_dawn_05_mean_plus = np.mean(KH_dawn_05_plus)
KH_dawn_08_mean_plus = np.mean(KH_dawn_08_plus)
KH_dawn_09_mean_plus = np.mean(KH_dawn_09_plus)
KH_dawn_10_mean_plus = np.mean(KH_dawn_10_plus)
KH_dawn_11_mean_plus = np.mean(KH_dawn_11_plus)
KH_dawn_16_mean_plus = np.mean(KH_dawn_16_plus)
KH_dawn_17_mean_plus = np.mean(KH_dawn_17_plus)
KH_dawn_18_mean_plus = np.mean(KH_dawn_18_plus)
KH_dawn_19_mean_plus = np.mean(KH_dawn_19_plus)
KH_dawn_20_mean_plus = np.mean(KH_dawn_20_plus)
KH_dawn_21_mean_plus = np.mean(KH_dawn_21_plus)
KH_dawn_24_mean_plus = np.mean(KH_dawn_24_plus)
KH_dawn_25_mean_plus = np.mean(KH_dawn_25_plus)
KH_dawn_27_mean_plus = np.mean(KH_dawn_27_plus)
KH_dawn_28_mean_plus = np.mean(KH_dawn_28_plus)
KH_dawn_34_mean_plus = np.mean(KH_dawn_34_plus)
KH_dawn_35_mean_plus = np.mean(KH_dawn_35_plus)

KH_dawn_01_mean_plus = np.mean(KH_dawn_01_plus)
KH_dawn_12_mean_plus = np.mean(KH_dawn_12_plus)
KH_dawn_15_mean_plus = np.mean(KH_dawn_15_plus)
KH_dawn_26_mean_plus = np.mean(KH_dawn_26_plus)

KH_dawn_means_plus = np.array([KH_dawn_01_mean_plus,KH_dawn_02_mean_plus,KH_dawn_03_mean_plus,KH_dawn_04_mean_plus,KH_dawn_05_mean_plus,KH_dawn_08_mean_plus,KH_dawn_09_mean_plus,KH_dawn_10_mean_plus,KH_dawn_11_mean,KH_dawn_12_mean_plus,KH_dawn_15_mean_plus,KH_dawn_16_mean_plus,KH_dawn_17_mean_plus,KH_dawn_18_mean_plus,KH_dawn_19_mean_plus,KH_dawn_20_mean_plus,KH_dawn_21_mean_plus,KH_dawn_24_mean_plus,KH_dawn_25_mean_plus,KH_dawn_26_mean_plus,KH_dawn_27_mean_plus,KH_dawn_28_mean_plus,KH_dawn_34_mean_plus,KH_dawn_35_mean_plus])

sd_KH_dawn_02_plus = standard_deviation(KH_dawn_02_plus,KH_dawn_02_mean_plus)
sd_KH_dawn_03_plus = standard_deviation(KH_dawn_03_plus,KH_dawn_03_mean_plus)
sd_KH_dawn_04_plus = standard_deviation(KH_dawn_04_plus,KH_dawn_04_mean_plus)
sd_KH_dawn_05_plus = standard_deviation(KH_dawn_05_plus,KH_dawn_05_mean_plus)
sd_KH_dawn_08_plus = standard_deviation(KH_dawn_08_plus,KH_dawn_08_mean_plus)
sd_KH_dawn_09_plus = standard_deviation(KH_dawn_09_plus,KH_dawn_09_mean_plus)
sd_KH_dawn_10_plus = standard_deviation(KH_dawn_10_plus,KH_dawn_10_mean_plus)
sd_KH_dawn_11_plus = standard_deviation(KH_dawn_11_plus,KH_dawn_11_mean_plus)
sd_KH_dawn_16_plus = standard_deviation(KH_dawn_16_plus,KH_dawn_16_mean_plus)
sd_KH_dawn_17_plus = standard_deviation(KH_dawn_17_plus,KH_dawn_17_mean_plus)
sd_KH_dawn_18_plus = standard_deviation(KH_dawn_18_plus,KH_dawn_18_mean_plus)
sd_KH_dawn_19_plus = standard_deviation(KH_dawn_19_plus,KH_dawn_19_mean_plus)
sd_KH_dawn_20_plus = standard_deviation(KH_dawn_20_plus,KH_dawn_20_mean_plus)
sd_KH_dawn_21_plus = standard_deviation(KH_dawn_21_plus,KH_dawn_21_mean_plus)
sd_KH_dawn_24_plus = standard_deviation(KH_dawn_24_plus,KH_dawn_24_mean_plus)
sd_KH_dawn_25_plus = standard_deviation(KH_dawn_25_plus,KH_dawn_25_mean_plus)
sd_KH_dawn_27_plus = standard_deviation(KH_dawn_27_plus,KH_dawn_27_mean_plus)
sd_KH_dawn_28_plus = standard_deviation(KH_dawn_28_plus,KH_dawn_28_mean_plus)
sd_KH_dawn_34_plus = standard_deviation(KH_dawn_34_plus,KH_dawn_34_mean_plus)
sd_KH_dawn_35_plus = standard_deviation(KH_dawn_35_plus,KH_dawn_35_mean_plus)

sd_KH_dawn_01_plus = standard_deviation(KH_dawn_01_plus,KH_dawn_01_mean_plus)
sd_KH_dawn_12_plus = standard_deviation(KH_dawn_12_plus,KH_dawn_12_mean_plus)
sd_KH_dawn_15_plus = standard_deviation(KH_dawn_15_plus,KH_dawn_15_mean_plus)
sd_KH_dawn_26_plus = standard_deviation(KH_dawn_26_plus,KH_dawn_26_mean_plus)


KH_dawn_02_median_plus = np.median(KH_dawn_02_plus)
KH_dawn_03_median_plus = np.median(KH_dawn_03_plus)
KH_dawn_04_median_plus = np.median(KH_dawn_04_plus)
KH_dawn_05_median_plus = np.median(KH_dawn_05_plus)
KH_dawn_08_median_plus = np.median(KH_dawn_08_plus)
KH_dawn_09_median_plus = np.median(KH_dawn_09_plus)
KH_dawn_10_median_plus = np.median(KH_dawn_10_plus)
KH_dawn_11_median_plus = np.median(KH_dawn_11_plus)
KH_dawn_16_median_plus = np.median(KH_dawn_16_plus)
KH_dawn_17_median_plus = np.median(KH_dawn_17_plus)
KH_dawn_18_median_plus = np.median(KH_dawn_18_plus)
KH_dawn_19_median_plus = np.median(KH_dawn_19_plus)
KH_dawn_20_median_plus = np.median(KH_dawn_20_plus)
KH_dawn_21_median_plus = np.median(KH_dawn_21_plus)
KH_dawn_24_median_plus = np.median(KH_dawn_24_plus)
KH_dawn_25_median_plus = np.median(KH_dawn_25_plus)
KH_dawn_27_median_plus = np.median(KH_dawn_27_plus)
KH_dawn_28_median_plus = np.median(KH_dawn_28_plus)
KH_dawn_34_median_plus = np.median(KH_dawn_34_plus)
KH_dawn_35_median_plus = np.median(KH_dawn_35_plus)

KH_dawn_01_median_plus = np.median(KH_dawn_01_plus)
KH_dawn_12_median_plus = np.median(KH_dawn_12_plus)
KH_dawn_15_median_plus = np.median(KH_dawn_15_plus)
KH_dawn_26_median_plus = np.median(KH_dawn_26_plus)

KH_dawn_medians_plus = np.array([KH_dawn_01_median_plus,KH_dawn_02_median_plus,KH_dawn_03_median_plus,KH_dawn_04_median_plus,KH_dawn_05_median_plus,KH_dawn_08_median_plus,KH_dawn_09_median_plus,KH_dawn_10_median_plus,KH_dawn_11_median_plus,KH_dawn_12_median_plus,KH_dawn_15_median_plus,KH_dawn_16_median_plus,KH_dawn_17_median_plus,KH_dawn_18_median_plus,KH_dawn_19_median_plus,KH_dawn_20_median_plus,KH_dawn_21_median_plus,KH_dawn_24_median_plus,KH_dawn_25_median_plus,KH_dawn_26_median_plus,KH_dawn_27_median_plus,KH_dawn_28_median_plus,KH_dawn_34_median_plus,KH_dawn_35_median_plus])

med_KH_dawn_02_err_plus = 1.2533 * sd_KH_dawn_02_plus
med_KH_dawn_03_err_plus = 1.2533 * sd_KH_dawn_03_plus
med_KH_dawn_04_err_plus = 1.2533 * sd_KH_dawn_04_plus
med_KH_dawn_05_err_plus = 1.2533 * sd_KH_dawn_05_plus
med_KH_dawn_08_err_plus = 1.2533 * sd_KH_dawn_08_plus
med_KH_dawn_09_err_plus = 1.2533 * sd_KH_dawn_09_plus
med_KH_dawn_10_err_plus = 1.2533 * sd_KH_dawn_10_plus
med_KH_dawn_11_err_plus = 1.2533 * sd_KH_dawn_11_plus
med_KH_dawn_16_err_plus = 1.2533 * sd_KH_dawn_16_plus
med_KH_dawn_17_err_plus = 1.2533 * sd_KH_dawn_17_plus
med_KH_dawn_18_err_plus = 1.2533 * sd_KH_dawn_18_plus
med_KH_dawn_19_err_plus = 1.2533 * sd_KH_dawn_19_plus
med_KH_dawn_20_err_plus = 1.2533 * sd_KH_dawn_20_plus
med_KH_dawn_21_err_plus = 1.2533 * sd_KH_dawn_21_plus
med_KH_dawn_24_err_plus = 1.2533 * sd_KH_dawn_24_plus
med_KH_dawn_25_err_plus = 1.2533 * sd_KH_dawn_25_plus
med_KH_dawn_27_err_plus = 1.2533 * sd_KH_dawn_27_plus
med_KH_dawn_28_err_plus = 1.2533 * sd_KH_dawn_28_plus
med_KH_dawn_34_err_plus = 1.2533 * sd_KH_dawn_34_plus
med_KH_dawn_35_err_plus = 1.2533 * sd_KH_dawn_35_plus

med_KH_dawn_01_err_plus = 1.2533 * sd_KH_dawn_01_plus
med_KH_dawn_12_err_plus = 1.2533 * sd_KH_dawn_12_plus
med_KH_dawn_15_err_plus = 1.2533 * sd_KH_dawn_15_plus
med_KH_dawn_26_err_plus = 1.2533 * sd_KH_dawn_26_plus

med_KH_dawn_errs_plus = np.array([med_KH_dawn_01_err_plus,med_KH_dawn_02_err_plus,med_KH_dawn_03_err_plus,med_KH_dawn_04_err_plus,med_KH_dawn_05_err_plus,med_KH_dawn_08_err_plus,med_KH_dawn_09_err_plus,med_KH_dawn_10_err_plus,med_KH_dawn_11_err_plus,med_KH_dawn_12_err_plus,med_KH_dawn_15_err_plus,med_KH_dawn_16_err_plus,med_KH_dawn_17_err_plus,med_KH_dawn_18_err_plus,med_KH_dawn_19_err_plus,med_KH_dawn_20_err_plus,med_KH_dawn_21_err_plus,med_KH_dawn_24_err_plus,med_KH_dawn_25_err_plus,med_KH_dawn_26_err_plus,med_KH_dawn_27_err_plus,med_KH_dawn_28_err_plus,med_KH_dawn_34_err_plus,med_KH_dawn_35_err_plus])

'''
-
'''
KH_dawn_02_mean_minus = np.mean(KH_dawn_02_minus)
KH_dawn_03_mean_minus = np.mean(KH_dawn_03_minus)
KH_dawn_04_mean_minus = np.mean(KH_dawn_04_minus)
KH_dawn_05_mean_minus = np.mean(KH_dawn_05_minus)
KH_dawn_08_mean_minus = np.mean(KH_dawn_08_minus)
KH_dawn_09_mean_minus = np.mean(KH_dawn_09_minus)
KH_dawn_10_mean_minus = np.mean(KH_dawn_10_minus)
KH_dawn_11_mean_minus = np.mean(KH_dawn_11_minus)
KH_dawn_16_mean_minus = np.mean(KH_dawn_16_minus)
KH_dawn_17_mean_minus = np.mean(KH_dawn_17_minus)
KH_dawn_18_mean_minus = np.mean(KH_dawn_18_minus)
KH_dawn_19_mean_minus = np.mean(KH_dawn_19_minus)
KH_dawn_20_mean_minus = np.mean(KH_dawn_20_minus)
KH_dawn_21_mean_minus = np.mean(KH_dawn_21_minus)
KH_dawn_24_mean_minus = np.mean(KH_dawn_24_minus)
KH_dawn_25_mean_minus = np.mean(KH_dawn_25_minus)
KH_dawn_27_mean_minus = np.mean(KH_dawn_27_minus)
KH_dawn_28_mean_minus = np.mean(KH_dawn_28_minus)
KH_dawn_34_mean_minus = np.mean(KH_dawn_34_minus)
KH_dawn_35_mean_minus = np.mean(KH_dawn_35_minus)

KH_dawn_01_mean_minus = np.mean(KH_dawn_01_minus)
KH_dawn_12_mean_minus = np.mean(KH_dawn_12_minus)
KH_dawn_15_mean_minus = np.mean(KH_dawn_15_minus)
KH_dawn_26_mean_minus = np.mean(KH_dawn_26_minus)

KH_dawn_means_minus = np.array([KH_dawn_01_mean_minus,KH_dawn_02_mean_minus,KH_dawn_03_mean_minus,KH_dawn_04_mean_minus,KH_dawn_05_mean_minus,KH_dawn_08_mean_minus,KH_dawn_09_mean_minus,KH_dawn_10_mean_minus,KH_dawn_11_mean_minus,KH_dawn_12_mean_minus,KH_dawn_15_mean_minus,KH_dawn_16_mean_minus,KH_dawn_17_mean_minus,KH_dawn_18_mean_minus,KH_dawn_19_mean_minus,KH_dawn_20_mean_minus,KH_dawn_21_mean_minus,KH_dawn_24_mean_minus,KH_dawn_25_mean_minus,KH_dawn_26_mean_minus,KH_dawn_27_mean_minus,KH_dawn_28_mean_minus,KH_dawn_34_mean_minus,KH_dawn_35_mean_minus])

sd_KH_dawn_02_minus = standard_deviation(KH_dawn_02_minus,KH_dawn_02_mean_minus)
sd_KH_dawn_03_minus = standard_deviation(KH_dawn_03_minus,KH_dawn_03_mean_minus)
sd_KH_dawn_04_minus = standard_deviation(KH_dawn_04_minus,KH_dawn_04_mean_minus)
sd_KH_dawn_05_minus = standard_deviation(KH_dawn_05_minus,KH_dawn_05_mean_minus)
sd_KH_dawn_08_minus = standard_deviation(KH_dawn_08_minus,KH_dawn_08_mean_minus)
sd_KH_dawn_09_minus = standard_deviation(KH_dawn_09_minus,KH_dawn_09_mean_minus)
sd_KH_dawn_10_minus = standard_deviation(KH_dawn_10_minus,KH_dawn_10_mean_minus)
sd_KH_dawn_11_minus = standard_deviation(KH_dawn_11_minus,KH_dawn_11_mean_minus)
sd_KH_dawn_16_minus = standard_deviation(KH_dawn_16_minus,KH_dawn_16_mean_minus)
sd_KH_dawn_17_minus = standard_deviation(KH_dawn_17_minus,KH_dawn_17_mean_minus)
sd_KH_dawn_18_minus = standard_deviation(KH_dawn_18_minus,KH_dawn_18_mean_minus)
sd_KH_dawn_19_minus = standard_deviation(KH_dawn_19_minus,KH_dawn_19_mean_minus)
sd_KH_dawn_20_minus = standard_deviation(KH_dawn_20_minus,KH_dawn_20_mean_minus)
sd_KH_dawn_21_minus = standard_deviation(KH_dawn_21_minus,KH_dawn_21_mean_minus)
sd_KH_dawn_24_minus = standard_deviation(KH_dawn_24_minus,KH_dawn_24_mean_minus)
sd_KH_dawn_25_minus = standard_deviation(KH_dawn_25_minus,KH_dawn_25_mean_minus)
sd_KH_dawn_27_minus = standard_deviation(KH_dawn_27_minus,KH_dawn_27_mean_minus)
sd_KH_dawn_28_minus = standard_deviation(KH_dawn_28_minus,KH_dawn_28_mean_minus)
sd_KH_dawn_34_minus = standard_deviation(KH_dawn_34_minus,KH_dawn_34_mean_minus)
sd_KH_dawn_35_minus = standard_deviation(KH_dawn_35_minus,KH_dawn_35_mean_minus)

sd_KH_dawn_01_minus = standard_deviation(KH_dawn_01_minus,KH_dawn_01_mean_minus)
sd_KH_dawn_12_minus = standard_deviation(KH_dawn_12_minus,KH_dawn_12_mean_minus)
sd_KH_dawn_15_minus = standard_deviation(KH_dawn_15_minus,KH_dawn_15_mean_minus)
sd_KH_dawn_26_minus = standard_deviation(KH_dawn_26_minus,KH_dawn_26_mean_minus)


KH_dawn_02_median_minus = np.median(KH_dawn_02_minus)
KH_dawn_03_median_minus = np.median(KH_dawn_03_minus)
KH_dawn_04_median_minus = np.median(KH_dawn_04_minus)
KH_dawn_05_median_minus = np.median(KH_dawn_05_minus)
KH_dawn_08_median_minus = np.median(KH_dawn_08_minus)
KH_dawn_09_median_minus = np.median(KH_dawn_09_minus)
KH_dawn_10_median_minus = np.median(KH_dawn_10_minus)
KH_dawn_11_median_minus = np.median(KH_dawn_11_minus)
KH_dawn_16_median_minus = np.median(KH_dawn_16_minus)
KH_dawn_17_median_minus = np.median(KH_dawn_17_minus)
KH_dawn_18_median_minus = np.median(KH_dawn_18_minus)
KH_dawn_19_median_minus = np.median(KH_dawn_19_minus)
KH_dawn_20_median_minus = np.median(KH_dawn_20_minus)
KH_dawn_21_median_minus = np.median(KH_dawn_21_minus)
KH_dawn_24_median_minus = np.median(KH_dawn_24_minus)
KH_dawn_25_median_minus = np.median(KH_dawn_25_minus)
KH_dawn_27_median_minus = np.median(KH_dawn_27_minus)
KH_dawn_28_median_minus = np.median(KH_dawn_28_minus)
KH_dawn_34_median_minus = np.median(KH_dawn_34_minus)
KH_dawn_35_median_minus = np.median(KH_dawn_35_minus)

KH_dawn_01_median_minus = np.median(KH_dawn_01_minus)
KH_dawn_12_median_minus = np.median(KH_dawn_12_minus)
KH_dawn_15_median_minus = np.median(KH_dawn_15_minus)
KH_dawn_26_median_minus = np.median(KH_dawn_26_minus)

KH_dawn_medians_minus = np.array([KH_dawn_01_median_minus,KH_dawn_02_median_minus,KH_dawn_03_median_minus,KH_dawn_04_median_minus,KH_dawn_05_median_minus,KH_dawn_08_median_minus,KH_dawn_09_median_minus,KH_dawn_10_median_minus,KH_dawn_11_median_minus,KH_dawn_12_median_minus,KH_dawn_15_median_minus,KH_dawn_16_median_minus,KH_dawn_17_median_minus,KH_dawn_18_median_minus,KH_dawn_19_median_minus,KH_dawn_20_median_minus,KH_dawn_21_median_minus,KH_dawn_24_median_minus,KH_dawn_25_median_minus,KH_dawn_26_median_minus,KH_dawn_27_median_minus,KH_dawn_28_median_minus,KH_dawn_34_median_minus,KH_dawn_35_median_minus])

med_KH_dawn_02_err_minus = 1.2533 * sd_KH_dawn_02_minus
med_KH_dawn_03_err_minus = 1.2533 * sd_KH_dawn_03_minus
med_KH_dawn_04_err_minus = 1.2533 * sd_KH_dawn_04_minus
med_KH_dawn_05_err_minus = 1.2533 * sd_KH_dawn_05_minus
med_KH_dawn_08_err_minus = 1.2533 * sd_KH_dawn_08_minus
med_KH_dawn_09_err_minus = 1.2533 * sd_KH_dawn_09_minus
med_KH_dawn_10_err_minus = 1.2533 * sd_KH_dawn_10_minus
med_KH_dawn_11_err_minus = 1.2533 * sd_KH_dawn_11_minus
med_KH_dawn_16_err_minus = 1.2533 * sd_KH_dawn_16_minus
med_KH_dawn_17_err_minus = 1.2533 * sd_KH_dawn_17_minus
med_KH_dawn_18_err_minus = 1.2533 * sd_KH_dawn_18_minus
med_KH_dawn_19_err_minus = 1.2533 * sd_KH_dawn_19_minus
med_KH_dawn_20_err_minus = 1.2533 * sd_KH_dawn_20_minus
med_KH_dawn_21_err_minus = 1.2533 * sd_KH_dawn_21_minus
med_KH_dawn_24_err_minus = 1.2533 * sd_KH_dawn_24_minus
med_KH_dawn_25_err_minus = 1.2533 * sd_KH_dawn_25_minus
med_KH_dawn_27_err_minus = 1.2533 * sd_KH_dawn_27_minus
med_KH_dawn_28_err_minus = 1.2533 * sd_KH_dawn_28_minus
med_KH_dawn_34_err_minus = 1.2533 * sd_KH_dawn_34_minus
med_KH_dawn_35_err_minus = 1.2533 * sd_KH_dawn_35_minus

med_KH_dawn_01_err_minus = 1.2533 * sd_KH_dawn_01_minus
med_KH_dawn_12_err_minus = 1.2533 * sd_KH_dawn_12_minus
med_KH_dawn_15_err_minus = 1.2533 * sd_KH_dawn_15_minus
med_KH_dawn_26_err_minus = 1.2533 * sd_KH_dawn_26_minus

med_KH_dawn_errs_minus = np.array([med_KH_dawn_01_err_minus,med_KH_dawn_02_err_minus,med_KH_dawn_03_err_minus,med_KH_dawn_04_err_minus,med_KH_dawn_05_err_minus,med_KH_dawn_08_err_minus,med_KH_dawn_09_err_minus,med_KH_dawn_10_err_minus,med_KH_dawn_11_err_minus,med_KH_dawn_12_err_minus,med_KH_dawn_15_err_minus,med_KH_dawn_16_err_minus,med_KH_dawn_17_err_minus,med_KH_dawn_18_err_minus,med_KH_dawn_19_err_minus,med_KH_dawn_20_err_minus,med_KH_dawn_21_err_minus,med_KH_dawn_24_err_minus,med_KH_dawn_25_err_minus,med_KH_dawn_26_err_minus,med_KH_dawn_27_err_minus,med_KH_dawn_28_err_minus,med_KH_dawn_34_err_minus,med_KH_dawn_35_err_minus])


'''
AVERAGE KH - DUSK
'''
KH_dusk_02_mean = np.mean(KH_dusk_02)
KH_dusk_03_mean = np.mean(KH_dusk_03)
KH_dusk_04_mean = np.mean(KH_dusk_04)
KH_dusk_05_mean = np.mean(KH_dusk_05)
KH_dusk_08_mean = np.mean(KH_dusk_08)
KH_dusk_09_mean = np.mean(KH_dusk_09)
KH_dusk_10_mean = np.mean(KH_dusk_10)
KH_dusk_11_mean = np.mean(KH_dusk_11)
KH_dusk_16_mean = np.mean(KH_dusk_16)
KH_dusk_17_mean = np.mean(KH_dusk_17)
KH_dusk_18_mean = np.mean(KH_dusk_18)
KH_dusk_19_mean = np.mean(KH_dusk_19)
KH_dusk_20_mean = np.mean(KH_dusk_20)
KH_dusk_21_mean = np.mean(KH_dusk_21)
KH_dusk_24_mean = np.mean(KH_dusk_24)
KH_dusk_25_mean = np.mean(KH_dusk_25)
KH_dusk_27_mean = np.mean(KH_dusk_27)
KH_dusk_28_mean = np.mean(KH_dusk_28)
KH_dusk_34_mean = np.mean(KH_dusk_34)
KH_dusk_35_mean = np.mean(KH_dusk_35)

KH_dusk_01_mean = np.mean(KH_dusk_01)
KH_dusk_12_mean = np.mean(KH_dusk_12)
KH_dusk_15_mean = np.mean(KH_dusk_15)
KH_dusk_26_mean = np.mean(KH_dusk_26)

KH_dusk_means = np.array([KH_dusk_01_mean,KH_dusk_02_mean,KH_dusk_03_mean,KH_dusk_04_mean,KH_dusk_05_mean,KH_dusk_08_mean,KH_dusk_09_mean,KH_dusk_10_mean,KH_dusk_11_mean,KH_dusk_12_mean,KH_dusk_15_mean,KH_dusk_16_mean,KH_dusk_17_mean,KH_dusk_18_mean,KH_dusk_19_mean,KH_dusk_20_mean,KH_dusk_21_mean,KH_dusk_24_mean,KH_dusk_25_mean,KH_dusk_26_mean,KH_dusk_27_mean,KH_dusk_28_mean,KH_dusk_34_mean,KH_dusk_35_mean])

sd_KH_dusk_02 = standard_deviation(KH_dusk_02,KH_dusk_02_mean)
sd_KH_dusk_03 = standard_deviation(KH_dusk_03,KH_dusk_03_mean)
sd_KH_dusk_04 = standard_deviation(KH_dusk_04,KH_dusk_04_mean)
sd_KH_dusk_05 = standard_deviation(KH_dusk_05,KH_dusk_05_mean)
sd_KH_dusk_08 = standard_deviation(KH_dusk_08,KH_dusk_08_mean)
sd_KH_dusk_09 = standard_deviation(KH_dusk_09,KH_dusk_09_mean)
sd_KH_dusk_10 = standard_deviation(KH_dusk_10,KH_dusk_10_mean)
sd_KH_dusk_11 = standard_deviation(KH_dusk_11,KH_dusk_11_mean)
sd_KH_dusk_16 = standard_deviation(KH_dusk_16,KH_dusk_16_mean)
sd_KH_dusk_17 = standard_deviation(KH_dusk_17,KH_dusk_17_mean)
sd_KH_dusk_18 = standard_deviation(KH_dusk_18,KH_dusk_18_mean)
sd_KH_dusk_19 = standard_deviation(KH_dusk_19,KH_dusk_19_mean)
sd_KH_dusk_20 = standard_deviation(KH_dusk_20,KH_dusk_20_mean)
sd_KH_dusk_21 = standard_deviation(KH_dusk_21,KH_dusk_21_mean)
sd_KH_dusk_24 = standard_deviation(KH_dusk_24,KH_dusk_24_mean)
sd_KH_dusk_25 = standard_deviation(KH_dusk_25,KH_dusk_25_mean)
sd_KH_dusk_27 = standard_deviation(KH_dusk_27,KH_dusk_27_mean)
sd_KH_dusk_28 = standard_deviation(KH_dusk_28,KH_dusk_28_mean)
sd_KH_dusk_34 = standard_deviation(KH_dusk_34,KH_dusk_34_mean)
sd_KH_dusk_35 = standard_deviation(KH_dusk_35,KH_dusk_35_mean)

sd_KH_dusk_01 = standard_deviation(KH_dusk_01,KH_dusk_01_mean)
sd_KH_dusk_12 = standard_deviation(KH_dusk_12,KH_dusk_12_mean)
sd_KH_dusk_15 = standard_deviation(KH_dusk_15,KH_dusk_15_mean)
sd_KH_dusk_26 = standard_deviation(KH_dusk_26,KH_dusk_26_mean)


KH_dusk_02_median = np.median(KH_dusk_02)
KH_dusk_03_median = np.median(KH_dusk_03)
KH_dusk_04_median = np.median(KH_dusk_04)
KH_dusk_05_median = np.median(KH_dusk_05)
KH_dusk_08_median = np.median(KH_dusk_08)
KH_dusk_09_median = np.median(KH_dusk_09)
KH_dusk_10_median = np.median(KH_dusk_10)
KH_dusk_11_median = np.median(KH_dusk_11)
KH_dusk_16_median = np.median(KH_dusk_16)
KH_dusk_17_median = np.median(KH_dusk_17)
KH_dusk_18_median = np.median(KH_dusk_18)
KH_dusk_19_median = np.median(KH_dusk_19)
KH_dusk_20_median = np.median(KH_dusk_20)
KH_dusk_21_median = np.median(KH_dusk_21)
KH_dusk_24_median = np.median(KH_dusk_24)
KH_dusk_25_median = np.median(KH_dusk_25)
KH_dusk_27_median = np.median(KH_dusk_27)
KH_dusk_28_median = np.median(KH_dusk_28)
KH_dusk_34_median = np.median(KH_dusk_34)
KH_dusk_35_median = np.median(KH_dusk_35)

KH_dusk_01_median = np.median(KH_dusk_01)
KH_dusk_12_median = np.median(KH_dusk_12)
KH_dusk_15_median = np.median(KH_dusk_15)
KH_dusk_26_median = np.median(KH_dusk_26)

KH_dusk_medians = np.array([KH_dusk_01_median,KH_dusk_02_median,KH_dusk_03_median,KH_dusk_04_median,KH_dusk_05_median,KH_dusk_08_median,KH_dusk_09_median,KH_dusk_10_median,KH_dusk_11_median,KH_dusk_12_median,KH_dusk_15_median,KH_dusk_16_median,KH_dusk_17_median,KH_dusk_18_median,KH_dusk_19_median,KH_dusk_20_median,KH_dusk_21_median,KH_dusk_24_median,KH_dusk_25_median,KH_dusk_26_median,KH_dusk_27_median,KH_dusk_28_median,KH_dusk_34_median,KH_dusk_35_median])

med_KH_dusk_02_err = 1.2533 * sd_KH_dusk_02
med_KH_dusk_03_err = 1.2533 * sd_KH_dusk_03
med_KH_dusk_04_err = 1.2533 * sd_KH_dusk_04
med_KH_dusk_05_err = 1.2533 * sd_KH_dusk_05
med_KH_dusk_08_err = 1.2533 * sd_KH_dusk_08
med_KH_dusk_09_err = 1.2533 * sd_KH_dusk_09
med_KH_dusk_10_err = 1.2533 * sd_KH_dusk_10
med_KH_dusk_11_err = 1.2533 * sd_KH_dusk_11
med_KH_dusk_16_err = 1.2533 * sd_KH_dusk_16
med_KH_dusk_17_err = 1.2533 * sd_KH_dusk_17
med_KH_dusk_18_err = 1.2533 * sd_KH_dusk_18
med_KH_dusk_19_err = 1.2533 * sd_KH_dusk_19
med_KH_dusk_20_err = 1.2533 * sd_KH_dusk_20
med_KH_dusk_21_err = 1.2533 * sd_KH_dusk_21
med_KH_dusk_24_err = 1.2533 * sd_KH_dusk_24
med_KH_dusk_25_err = 1.2533 * sd_KH_dusk_25
med_KH_dusk_27_err = 1.2533 * sd_KH_dusk_27
med_KH_dusk_28_err = 1.2533 * sd_KH_dusk_28
med_KH_dusk_34_err = 1.2533 * sd_KH_dusk_34
med_KH_dusk_35_err = 1.2533 * sd_KH_dusk_35

med_KH_dusk_01_err = 1.2533 * sd_KH_dusk_01
med_KH_dusk_12_err = 1.2533 * sd_KH_dusk_12
med_KH_dusk_15_err = 1.2533 * sd_KH_dusk_15
med_KH_dusk_26_err = 1.2533 * sd_KH_dusk_26

med_KH_dusk_errs = np.array([med_KH_dusk_01_err,med_KH_dusk_02_err,med_KH_dusk_03_err,med_KH_dusk_04_err,med_KH_dusk_05_err,med_KH_dusk_08_err,med_KH_dusk_09_err,med_KH_dusk_10_err,med_KH_dusk_11_err,med_KH_dusk_12_err,med_KH_dusk_15_err,med_KH_dusk_16_err,med_KH_dusk_17_err,med_KH_dusk_18_err,med_KH_dusk_19_err,med_KH_dusk_20_err,med_KH_dusk_21_err,med_KH_dusk_24_err,med_KH_dusk_25_err,med_KH_dusk_26_err,med_KH_dusk_27_err,med_KH_dusk_28_err,med_KH_dusk_34_err,med_KH_dusk_35_err])

'''
+
'''

KH_dusk_02_mean_plus = np.mean(KH_dusk_02_plus)
KH_dusk_03_mean_plus = np.mean(KH_dusk_03_plus)
KH_dusk_04_mean_plus = np.mean(KH_dusk_04_plus)
KH_dusk_05_mean_plus = np.mean(KH_dusk_05_plus)
KH_dusk_08_mean_plus = np.mean(KH_dusk_08_plus)
KH_dusk_09_mean_plus = np.mean(KH_dusk_09_plus)
KH_dusk_10_mean_plus = np.mean(KH_dusk_10_plus)
KH_dusk_11_mean_plus = np.mean(KH_dusk_11_plus)
KH_dusk_16_mean_plus = np.mean(KH_dusk_16_plus)
KH_dusk_17_mean_plus = np.mean(KH_dusk_17_plus)
KH_dusk_18_mean_plus = np.mean(KH_dusk_18_plus)
KH_dusk_19_mean_plus = np.mean(KH_dusk_19_plus)
KH_dusk_20_mean_plus = np.mean(KH_dusk_20_plus)
KH_dusk_21_mean_plus = np.mean(KH_dusk_21_plus)
KH_dusk_24_mean_plus = np.mean(KH_dusk_24_plus)
KH_dusk_25_mean_plus = np.mean(KH_dusk_25_plus)
KH_dusk_27_mean_plus = np.mean(KH_dusk_27_plus)
KH_dusk_28_mean_plus = np.mean(KH_dusk_28_plus)
KH_dusk_34_mean_plus = np.mean(KH_dusk_34_plus)
KH_dusk_35_mean_plus = np.mean(KH_dusk_35_plus)

KH_dusk_01_mean_plus = np.mean(KH_dusk_01_plus)
KH_dusk_12_mean_plus = np.mean(KH_dusk_12_plus)
KH_dusk_15_mean_plus = np.mean(KH_dusk_15_plus)
KH_dusk_26_mean_plus = np.mean(KH_dusk_26_plus)

KH_dusk_means_plus = np.array([KH_dusk_01_mean_plus,KH_dusk_02_mean_plus,KH_dusk_03_mean_plus,KH_dusk_04_mean_plus,KH_dusk_05_mean_plus,KH_dusk_08_mean_plus,KH_dusk_09_mean_plus,KH_dusk_10_mean_plus,KH_dusk_11_mean,KH_dusk_12_mean_plus,KH_dusk_15_mean_plus,KH_dusk_16_mean_plus,KH_dusk_17_mean_plus,KH_dusk_18_mean_plus,KH_dusk_19_mean_plus,KH_dusk_20_mean_plus,KH_dusk_21_mean_plus,KH_dusk_24_mean_plus,KH_dusk_25_mean_plus,KH_dusk_26_mean_plus,KH_dusk_27_mean_plus,KH_dusk_28_mean_plus,KH_dusk_34_mean_plus,KH_dusk_35_mean_plus])

sd_KH_dusk_02_plus = standard_deviation(KH_dusk_02_plus,KH_dusk_02_mean_plus)
sd_KH_dusk_03_plus = standard_deviation(KH_dusk_03_plus,KH_dusk_03_mean_plus)
sd_KH_dusk_04_plus = standard_deviation(KH_dusk_04_plus,KH_dusk_04_mean_plus)
sd_KH_dusk_05_plus = standard_deviation(KH_dusk_05_plus,KH_dusk_05_mean_plus)
sd_KH_dusk_08_plus = standard_deviation(KH_dusk_08_plus,KH_dusk_08_mean_plus)
sd_KH_dusk_09_plus = standard_deviation(KH_dusk_09_plus,KH_dusk_09_mean_plus)
sd_KH_dusk_10_plus = standard_deviation(KH_dusk_10_plus,KH_dusk_10_mean_plus)
sd_KH_dusk_11_plus = standard_deviation(KH_dusk_11_plus,KH_dusk_11_mean_plus)
sd_KH_dusk_16_plus = standard_deviation(KH_dusk_16_plus,KH_dusk_16_mean_plus)
sd_KH_dusk_17_plus = standard_deviation(KH_dusk_17_plus,KH_dusk_17_mean_plus)
sd_KH_dusk_18_plus = standard_deviation(KH_dusk_18_plus,KH_dusk_18_mean_plus)
sd_KH_dusk_19_plus = standard_deviation(KH_dusk_19_plus,KH_dusk_19_mean_plus)
sd_KH_dusk_20_plus = standard_deviation(KH_dusk_20_plus,KH_dusk_20_mean_plus)
sd_KH_dusk_21_plus = standard_deviation(KH_dusk_21_plus,KH_dusk_21_mean_plus)
sd_KH_dusk_24_plus = standard_deviation(KH_dusk_24_plus,KH_dusk_24_mean_plus)
sd_KH_dusk_25_plus = standard_deviation(KH_dusk_25_plus,KH_dusk_25_mean_plus)
sd_KH_dusk_27_plus = standard_deviation(KH_dusk_27_plus,KH_dusk_27_mean_plus)
sd_KH_dusk_28_plus = standard_deviation(KH_dusk_28_plus,KH_dusk_28_mean_plus)
sd_KH_dusk_34_plus = standard_deviation(KH_dusk_34_plus,KH_dusk_34_mean_plus)
sd_KH_dusk_35_plus = standard_deviation(KH_dusk_35_plus,KH_dusk_35_mean_plus)

sd_KH_dusk_01_plus = standard_deviation(KH_dusk_01_plus,KH_dusk_01_mean_plus)
sd_KH_dusk_12_plus = standard_deviation(KH_dusk_12_plus,KH_dusk_12_mean_plus)
sd_KH_dusk_15_plus = standard_deviation(KH_dusk_15_plus,KH_dusk_15_mean_plus)
sd_KH_dusk_26_plus = standard_deviation(KH_dusk_26_plus,KH_dusk_26_mean_plus)


KH_dusk_02_median_plus = np.median(KH_dusk_02_plus)
KH_dusk_03_median_plus = np.median(KH_dusk_03_plus)
KH_dusk_04_median_plus = np.median(KH_dusk_04_plus)
KH_dusk_05_median_plus = np.median(KH_dusk_05_plus)
KH_dusk_08_median_plus = np.median(KH_dusk_08_plus)
KH_dusk_09_median_plus = np.median(KH_dusk_09_plus)
KH_dusk_10_median_plus = np.median(KH_dusk_10_plus)
KH_dusk_11_median_plus = np.median(KH_dusk_11_plus)
KH_dusk_16_median_plus = np.median(KH_dusk_16_plus)
KH_dusk_17_median_plus = np.median(KH_dusk_17_plus)
KH_dusk_18_median_plus = np.median(KH_dusk_18_plus)
KH_dusk_19_median_plus = np.median(KH_dusk_19_plus)
KH_dusk_20_median_plus = np.median(KH_dusk_20_plus)
KH_dusk_21_median_plus = np.median(KH_dusk_21_plus)
KH_dusk_24_median_plus = np.median(KH_dusk_24_plus)
KH_dusk_25_median_plus = np.median(KH_dusk_25_plus)
KH_dusk_27_median_plus = np.median(KH_dusk_27_plus)
KH_dusk_28_median_plus = np.median(KH_dusk_28_plus)
KH_dusk_34_median_plus = np.median(KH_dusk_34_plus)
KH_dusk_35_median_plus = np.median(KH_dusk_35_plus)

KH_dusk_01_median_plus = np.median(KH_dusk_01_plus)
KH_dusk_12_median_plus = np.median(KH_dusk_12_plus)
KH_dusk_15_median_plus = np.median(KH_dusk_15_plus)
KH_dusk_26_median_plus = np.median(KH_dusk_26_plus)

KH_dusk_medians_plus = np.array([KH_dusk_01_median_plus,KH_dusk_02_median_plus,KH_dusk_03_median_plus,KH_dusk_04_median_plus,KH_dusk_05_median_plus,KH_dusk_08_median_plus,KH_dusk_09_median_plus,KH_dusk_10_median_plus,KH_dusk_11_median_plus,KH_dusk_12_median_plus,KH_dusk_15_median_plus,KH_dusk_16_median_plus,KH_dusk_17_median_plus,KH_dusk_18_median_plus,KH_dusk_19_median_plus,KH_dusk_20_median_plus,KH_dusk_21_median_plus,KH_dusk_24_median_plus,KH_dusk_25_median_plus,KH_dusk_26_median_plus,KH_dusk_27_median_plus,KH_dusk_28_median_plus,KH_dusk_34_median_plus,KH_dusk_35_median_plus])

med_KH_dusk_02_err_plus = 1.2533 * sd_KH_dusk_02_plus
med_KH_dusk_03_err_plus = 1.2533 * sd_KH_dusk_03_plus
med_KH_dusk_04_err_plus = 1.2533 * sd_KH_dusk_04_plus
med_KH_dusk_05_err_plus = 1.2533 * sd_KH_dusk_05_plus
med_KH_dusk_08_err_plus = 1.2533 * sd_KH_dusk_08_plus
med_KH_dusk_09_err_plus = 1.2533 * sd_KH_dusk_09_plus
med_KH_dusk_10_err_plus = 1.2533 * sd_KH_dusk_10_plus
med_KH_dusk_11_err_plus = 1.2533 * sd_KH_dusk_11_plus
med_KH_dusk_16_err_plus = 1.2533 * sd_KH_dusk_16_plus
med_KH_dusk_17_err_plus = 1.2533 * sd_KH_dusk_17_plus
med_KH_dusk_18_err_plus = 1.2533 * sd_KH_dusk_18_plus
med_KH_dusk_19_err_plus = 1.2533 * sd_KH_dusk_19_plus
med_KH_dusk_20_err_plus = 1.2533 * sd_KH_dusk_20_plus
med_KH_dusk_21_err_plus = 1.2533 * sd_KH_dusk_21_plus
med_KH_dusk_24_err_plus = 1.2533 * sd_KH_dusk_24_plus
med_KH_dusk_25_err_plus = 1.2533 * sd_KH_dusk_25_plus
med_KH_dusk_27_err_plus = 1.2533 * sd_KH_dusk_27_plus
med_KH_dusk_28_err_plus = 1.2533 * sd_KH_dusk_28_plus
med_KH_dusk_34_err_plus = 1.2533 * sd_KH_dusk_34_plus
med_KH_dusk_35_err_plus = 1.2533 * sd_KH_dusk_35_plus

med_KH_dusk_01_err_plus = 1.2533 * sd_KH_dusk_01_plus
med_KH_dusk_12_err_plus = 1.2533 * sd_KH_dusk_12_plus
med_KH_dusk_15_err_plus = 1.2533 * sd_KH_dusk_15_plus
med_KH_dusk_26_err_plus = 1.2533 * sd_KH_dusk_26_plus

med_KH_dusk_errs_plus = np.array([med_KH_dusk_01_err_plus,med_KH_dusk_02_err_plus,med_KH_dusk_03_err_plus,med_KH_dusk_04_err_plus,med_KH_dusk_05_err_plus,med_KH_dusk_08_err_plus,med_KH_dusk_09_err_plus,med_KH_dusk_10_err_plus,med_KH_dusk_11_err_plus,med_KH_dusk_12_err_plus,med_KH_dusk_15_err_plus,med_KH_dusk_16_err_plus,med_KH_dusk_17_err_plus,med_KH_dusk_18_err_plus,med_KH_dusk_19_err_plus,med_KH_dusk_20_err_plus,med_KH_dusk_21_err_plus,med_KH_dusk_24_err_plus,med_KH_dusk_25_err_plus,med_KH_dusk_26_err_plus,med_KH_dusk_27_err_plus,med_KH_dusk_28_err_plus,med_KH_dusk_34_err_plus,med_KH_dusk_35_err_plus])

'''
-
'''
KH_dusk_02_mean_minus = np.mean(KH_dusk_02_minus)
KH_dusk_03_mean_minus = np.mean(KH_dusk_03_minus)
KH_dusk_04_mean_minus = np.mean(KH_dusk_04_minus)
KH_dusk_05_mean_minus = np.mean(KH_dusk_05_minus)
KH_dusk_08_mean_minus = np.mean(KH_dusk_08_minus)
KH_dusk_09_mean_minus = np.mean(KH_dusk_09_minus)
KH_dusk_10_mean_minus = np.mean(KH_dusk_10_minus)
KH_dusk_11_mean_minus = np.mean(KH_dusk_11_minus)
KH_dusk_16_mean_minus = np.mean(KH_dusk_16_minus)
KH_dusk_17_mean_minus = np.mean(KH_dusk_17_minus)
KH_dusk_18_mean_minus = np.mean(KH_dusk_18_minus)
KH_dusk_19_mean_minus = np.mean(KH_dusk_19_minus)
KH_dusk_20_mean_minus = np.mean(KH_dusk_20_minus)
KH_dusk_21_mean_minus = np.mean(KH_dusk_21_minus)
KH_dusk_24_mean_minus = np.mean(KH_dusk_24_minus)
KH_dusk_25_mean_minus = np.mean(KH_dusk_25_minus)
KH_dusk_27_mean_minus = np.mean(KH_dusk_27_minus)
KH_dusk_28_mean_minus = np.mean(KH_dusk_28_minus)
KH_dusk_34_mean_minus = np.mean(KH_dusk_34_minus)
KH_dusk_35_mean_minus = np.mean(KH_dusk_35_minus)

KH_dusk_01_mean_minus = np.mean(KH_dusk_01_minus)
KH_dusk_12_mean_minus = np.mean(KH_dusk_12_minus)
KH_dusk_15_mean_minus = np.mean(KH_dusk_15_minus)
KH_dusk_26_mean_minus = np.mean(KH_dusk_26_minus)

KH_dusk_means_minus = np.array([KH_dusk_01_mean_minus,KH_dusk_02_mean_minus,KH_dusk_03_mean_minus,KH_dusk_04_mean_minus,KH_dusk_05_mean_minus,KH_dusk_08_mean_minus,KH_dusk_09_mean_minus,KH_dusk_10_mean_minus,KH_dusk_11_mean_minus,KH_dusk_12_mean_minus,KH_dusk_15_mean_minus,KH_dusk_16_mean_minus,KH_dusk_17_mean_minus,KH_dusk_18_mean_minus,KH_dusk_19_mean_minus,KH_dusk_20_mean_minus,KH_dusk_21_mean_minus,KH_dusk_24_mean_minus,KH_dusk_25_mean_minus,KH_dusk_26_mean_minus,KH_dusk_27_mean_minus,KH_dusk_28_mean_minus,KH_dusk_34_mean_minus,KH_dusk_35_mean_minus])

sd_KH_dusk_02_minus = standard_deviation(KH_dusk_02_minus,KH_dusk_02_mean_minus)
sd_KH_dusk_03_minus = standard_deviation(KH_dusk_03_minus,KH_dusk_03_mean_minus)
sd_KH_dusk_04_minus = standard_deviation(KH_dusk_04_minus,KH_dusk_04_mean_minus)
sd_KH_dusk_05_minus = standard_deviation(KH_dusk_05_minus,KH_dusk_05_mean_minus)
sd_KH_dusk_08_minus = standard_deviation(KH_dusk_08_minus,KH_dusk_08_mean_minus)
sd_KH_dusk_09_minus = standard_deviation(KH_dusk_09_minus,KH_dusk_09_mean_minus)
sd_KH_dusk_10_minus = standard_deviation(KH_dusk_10_minus,KH_dusk_10_mean_minus)
sd_KH_dusk_11_minus = standard_deviation(KH_dusk_11_minus,KH_dusk_11_mean_minus)
sd_KH_dusk_16_minus = standard_deviation(KH_dusk_16_minus,KH_dusk_16_mean_minus)
sd_KH_dusk_17_minus = standard_deviation(KH_dusk_17_minus,KH_dusk_17_mean_minus)
sd_KH_dusk_18_minus = standard_deviation(KH_dusk_18_minus,KH_dusk_18_mean_minus)
sd_KH_dusk_19_minus = standard_deviation(KH_dusk_19_minus,KH_dusk_19_mean_minus)
sd_KH_dusk_20_minus = standard_deviation(KH_dusk_20_minus,KH_dusk_20_mean_minus)
sd_KH_dusk_21_minus = standard_deviation(KH_dusk_21_minus,KH_dusk_21_mean_minus)
sd_KH_dusk_24_minus = standard_deviation(KH_dusk_24_minus,KH_dusk_24_mean_minus)
sd_KH_dusk_25_minus = standard_deviation(KH_dusk_25_minus,KH_dusk_25_mean_minus)
sd_KH_dusk_27_minus = standard_deviation(KH_dusk_27_minus,KH_dusk_27_mean_minus)
sd_KH_dusk_28_minus = standard_deviation(KH_dusk_28_minus,KH_dusk_28_mean_minus)
sd_KH_dusk_34_minus = standard_deviation(KH_dusk_34_minus,KH_dusk_34_mean_minus)
sd_KH_dusk_35_minus = standard_deviation(KH_dusk_35_minus,KH_dusk_35_mean_minus)

sd_KH_dusk_01_minus = standard_deviation(KH_dusk_01_minus,KH_dusk_01_mean_minus)
sd_KH_dusk_12_minus = standard_deviation(KH_dusk_12_minus,KH_dusk_12_mean_minus)
sd_KH_dusk_15_minus = standard_deviation(KH_dusk_15_minus,KH_dusk_15_mean_minus)
sd_KH_dusk_26_minus = standard_deviation(KH_dusk_26_minus,KH_dusk_26_mean_minus)


KH_dusk_02_median_minus = np.median(KH_dusk_02_minus)
KH_dusk_03_median_minus = np.median(KH_dusk_03_minus)
KH_dusk_04_median_minus = np.median(KH_dusk_04_minus)
KH_dusk_05_median_minus = np.median(KH_dusk_05_minus)
KH_dusk_08_median_minus = np.median(KH_dusk_08_minus)
KH_dusk_09_median_minus = np.median(KH_dusk_09_minus)
KH_dusk_10_median_minus = np.median(KH_dusk_10_minus)
KH_dusk_11_median_minus = np.median(KH_dusk_11_minus)
KH_dusk_16_median_minus = np.median(KH_dusk_16_minus)
KH_dusk_17_median_minus = np.median(KH_dusk_17_minus)
KH_dusk_18_median_minus = np.median(KH_dusk_18_minus)
KH_dusk_19_median_minus = np.median(KH_dusk_19_minus)
KH_dusk_20_median_minus = np.median(KH_dusk_20_minus)
KH_dusk_21_median_minus = np.median(KH_dusk_21_minus)
KH_dusk_24_median_minus = np.median(KH_dusk_24_minus)
KH_dusk_25_median_minus = np.median(KH_dusk_25_minus)
KH_dusk_27_median_minus = np.median(KH_dusk_27_minus)
KH_dusk_28_median_minus = np.median(KH_dusk_28_minus)
KH_dusk_34_median_minus = np.median(KH_dusk_34_minus)
KH_dusk_35_median_minus = np.median(KH_dusk_35_minus)

KH_dusk_01_median_minus = np.median(KH_dusk_01_minus)
KH_dusk_12_median_minus = np.median(KH_dusk_12_minus)
KH_dusk_15_median_minus = np.median(KH_dusk_15_minus)
KH_dusk_26_median_minus = np.median(KH_dusk_26_minus)

KH_dusk_medians_minus = np.array([KH_dusk_01_median_minus,KH_dusk_02_median_minus,KH_dusk_03_median_minus,KH_dusk_04_median_minus,KH_dusk_05_median_minus,KH_dusk_08_median_minus,KH_dusk_09_median_minus,KH_dusk_10_median_minus,KH_dusk_11_median_minus,KH_dusk_12_median_minus,KH_dusk_15_median_minus,KH_dusk_16_median_minus,KH_dusk_17_median_minus,KH_dusk_18_median_minus,KH_dusk_19_median_minus,KH_dusk_20_median_minus,KH_dusk_21_median_minus,KH_dusk_24_median_minus,KH_dusk_25_median_minus,KH_dusk_26_median_minus,KH_dusk_27_median_minus,KH_dusk_28_median_minus,KH_dusk_34_median_minus,KH_dusk_35_median_minus])

med_KH_dusk_02_err_minus = 1.2533 * sd_KH_dusk_02_minus
med_KH_dusk_03_err_minus = 1.2533 * sd_KH_dusk_03_minus
med_KH_dusk_04_err_minus = 1.2533 * sd_KH_dusk_04_minus
med_KH_dusk_05_err_minus = 1.2533 * sd_KH_dusk_05_minus
med_KH_dusk_08_err_minus = 1.2533 * sd_KH_dusk_08_minus
med_KH_dusk_09_err_minus = 1.2533 * sd_KH_dusk_09_minus
med_KH_dusk_10_err_minus = 1.2533 * sd_KH_dusk_10_minus
med_KH_dusk_11_err_minus = 1.2533 * sd_KH_dusk_11_minus
med_KH_dusk_16_err_minus = 1.2533 * sd_KH_dusk_16_minus
med_KH_dusk_17_err_minus = 1.2533 * sd_KH_dusk_17_minus
med_KH_dusk_18_err_minus = 1.2533 * sd_KH_dusk_18_minus
med_KH_dusk_19_err_minus = 1.2533 * sd_KH_dusk_19_minus
med_KH_dusk_20_err_minus = 1.2533 * sd_KH_dusk_20_minus
med_KH_dusk_21_err_minus = 1.2533 * sd_KH_dusk_21_minus
med_KH_dusk_24_err_minus = 1.2533 * sd_KH_dusk_24_minus
med_KH_dusk_25_err_minus = 1.2533 * sd_KH_dusk_25_minus
med_KH_dusk_27_err_minus = 1.2533 * sd_KH_dusk_27_minus
med_KH_dusk_28_err_minus = 1.2533 * sd_KH_dusk_28_minus
med_KH_dusk_34_err_minus = 1.2533 * sd_KH_dusk_34_minus
med_KH_dusk_35_err_minus = 1.2533 * sd_KH_dusk_35_minus

med_KH_dusk_01_err_minus = 1.2533 * sd_KH_dusk_01_minus
med_KH_dusk_12_err_minus = 1.2533 * sd_KH_dusk_12_minus
med_KH_dusk_15_err_minus = 1.2533 * sd_KH_dusk_15_minus
med_KH_dusk_26_err_minus = 1.2533 * sd_KH_dusk_26_minus

med_KH_dusk_errs_minus = np.array([med_KH_dusk_01_err_minus,med_KH_dusk_02_err_minus,med_KH_dusk_03_err_minus,med_KH_dusk_04_err_minus,med_KH_dusk_05_err_minus,med_KH_dusk_08_err_minus,med_KH_dusk_09_err_minus,med_KH_dusk_10_err_minus,med_KH_dusk_11_err_minus,med_KH_dusk_12_err_minus,med_KH_dusk_15_err_minus,med_KH_dusk_16_err_minus,med_KH_dusk_17_err_minus,med_KH_dusk_18_err_minus,med_KH_dusk_19_err_minus,med_KH_dusk_20_err_minus,med_KH_dusk_21_err_minus,med_KH_dusk_24_err_minus,med_KH_dusk_25_err_minus,med_KH_dusk_26_err_minus,med_KH_dusk_27_err_minus,med_KH_dusk_28_err_minus,med_KH_dusk_34_err_minus,med_KH_dusk_35_err_minus])

'''
AVERAGE B PERPENDIDUCLAR
'''
b_perp_02_mean = np.mean(b_perp_02)
b_perp_03_mean = np.mean(b_perp_03)
b_perp_04_mean = np.mean(b_perp_04)
b_perp_05_mean = np.mean(b_perp_05)
b_perp_08_mean = np.mean(b_perp_08)
b_perp_09_mean = np.mean(b_perp_09)
b_perp_10_mean = np.mean(b_perp_10)
b_perp_11_mean = np.mean(b_perp_11)
b_perp_16_mean = np.mean(b_perp_16)
b_perp_17_mean = np.mean(b_perp_17)
b_perp_18_mean = np.mean(b_perp_18)
b_perp_19_mean = np.mean(b_perp_19)
b_perp_20_mean = np.mean(b_perp_20)
b_perp_21_mean = np.mean(b_perp_21)
b_perp_24_mean = np.mean(b_perp_24)
b_perp_25_mean = np.mean(b_perp_25)
b_perp_27_mean = np.mean(b_perp_27)
b_perp_28_mean = np.mean(b_perp_28)
b_perp_34_mean = np.mean(b_perp_34)
b_perp_35_mean = np.mean(b_perp_35)

# high cmls
b_perp_01_mean = np.mean(b_perp_01)
b_perp_12_mean = np.mean(b_perp_12)
b_perp_15_mean = np.mean(b_perp_15)
b_perp_26_mean = np.mean(b_perp_26)

#b_perp_means = np.array([b_perp_02_mean,b_perp_03_mean,b_perp_04_mean,b_perp_05_mean,b_perp_08_mean,b_perp_09_mean,b_perp_10_mean,b_perp_11_mean,b_perp_16_mean,b_perp_17_mean,b_perp_18_mean,b_perp_19_mean,b_perp_20_mean,b_perp_21_mean,b_perp_24_mean,b_perp_25_mean,b_perp_27_mean,b_perp_28_mean,b_perp_34_mean,b_perp_35_mean])
b_perp_means = np.array([b_perp_01_mean,b_perp_02_mean,b_perp_03_mean,b_perp_04_mean,b_perp_05_mean,b_perp_08_mean,b_perp_09_mean,b_perp_10_mean,b_perp_11_mean,b_perp_12_mean,b_perp_15_mean,b_perp_16_mean,b_perp_17_mean,b_perp_18_mean,b_perp_19_mean,b_perp_20_mean,b_perp_21_mean,b_perp_24_mean,b_perp_25_mean,b_perp_26_mean,b_perp_27_mean,b_perp_28_mean,b_perp_34_mean,b_perp_35_mean])


sd_b_perp_02 = standard_deviation(b_perp_02,b_perp_02_mean)
sd_b_perp_03 = standard_deviation(b_perp_03,b_perp_03_mean)
sd_b_perp_04 = standard_deviation(b_perp_04,b_perp_04_mean)
sd_b_perp_05 = standard_deviation(b_perp_05,b_perp_05_mean)
sd_b_perp_08 = standard_deviation(b_perp_08,b_perp_08_mean)
sd_b_perp_09 = standard_deviation(b_perp_09,b_perp_09_mean)
sd_b_perp_10 = standard_deviation(b_perp_10,b_perp_10_mean)
sd_b_perp_11 = standard_deviation(b_perp_11,b_perp_11_mean)
sd_b_perp_16 = standard_deviation(b_perp_16,b_perp_16_mean)
sd_b_perp_17 = standard_deviation(b_perp_17,b_perp_17_mean)
sd_b_perp_18 = standard_deviation(b_perp_18,b_perp_18_mean)
sd_b_perp_19 = standard_deviation(b_perp_19,b_perp_19_mean)
sd_b_perp_20 = standard_deviation(b_perp_20,b_perp_20_mean)
sd_b_perp_21 = standard_deviation(b_perp_21,b_perp_21_mean)
sd_b_perp_24 = standard_deviation(b_perp_24,b_perp_24_mean)
sd_b_perp_25 = standard_deviation(b_perp_25,b_perp_25_mean)
sd_b_perp_27 = standard_deviation(b_perp_27,b_perp_27_mean)
sd_b_perp_28 = standard_deviation(b_perp_28,b_perp_28_mean)
sd_b_perp_34 = standard_deviation(b_perp_34,b_perp_34_mean)
sd_b_perp_35 = standard_deviation(b_perp_35,b_perp_35_mean)

# high cmls
sd_b_perp_01 = standard_deviation(b_perp_01,b_perp_01_mean)
sd_b_perp_12 = standard_deviation(b_perp_12,b_perp_12_mean)
sd_b_perp_15 = standard_deviation(b_perp_15,b_perp_15_mean)
sd_b_perp_26 = standard_deviation(b_perp_26,b_perp_26_mean)


b_perp_02_median = np.median(b_perp_02)
b_perp_03_median = np.median(b_perp_03)
b_perp_04_median = np.median(b_perp_04)
b_perp_05_median = np.median(b_perp_05)
b_perp_08_median = np.median(b_perp_08)
b_perp_09_median = np.median(b_perp_09)
b_perp_10_median = np.median(b_perp_10)
b_perp_11_median = np.median(b_perp_11)
b_perp_16_median = np.median(b_perp_16)
b_perp_17_median = np.median(b_perp_17)
b_perp_18_median = np.median(b_perp_18)
b_perp_19_median = np.median(b_perp_19)
b_perp_20_median = np.median(b_perp_20)
b_perp_21_median = np.median(b_perp_21)
b_perp_24_median = np.median(b_perp_24)
b_perp_25_median = np.median(b_perp_25)
b_perp_27_median = np.median(b_perp_27)
b_perp_28_median = np.median(b_perp_28)
b_perp_34_median = np.median(b_perp_34)
b_perp_35_median = np.median(b_perp_35)

# high cmls
b_perp_01_median = np.median(b_perp_01)
b_perp_12_median = np.median(b_perp_12)
b_perp_15_median = np.median(b_perp_15)
b_perp_26_median = np.median(b_perp_26)

#b_perp_medians = np.array([b_perp_02_median,b_perp_03_median,b_perp_04_median,b_perp_05_median,b_perp_08_median,b_perp_09_median,b_perp_10_median,b_perp_11_median,b_perp_16_median,b_perp_17_median,b_perp_18_median,b_perp_19_median,b_perp_20_median,b_perp_21_median,b_perp_24_median,b_perp_25_median,b_perp_27_median,b_perp_28_median,b_perp_34_median,b_perp_35_median])
b_perp_medians = np.array([b_perp_01_median,b_perp_02_median,b_perp_03_median,b_perp_04_median,b_perp_05_median,b_perp_08_median,b_perp_09_median,b_perp_10_median,b_perp_11_median,b_perp_12_median,b_perp_15_median,b_perp_16_median,b_perp_17_median,b_perp_18_median,b_perp_19_median,b_perp_20_median,b_perp_21_median,b_perp_24_median,b_perp_25_median,b_perp_26_median,b_perp_27_median,b_perp_28_median,b_perp_34_median,b_perp_35_median])


med_b_perp_02_err = 1.2533 * sd_b_perp_02
med_b_perp_03_err = 1.2533 * sd_b_perp_03
med_b_perp_04_err = 1.2533 * sd_b_perp_04
med_b_perp_05_err = 1.2533 * sd_b_perp_05
med_b_perp_08_err = 1.2533 * sd_b_perp_08
med_b_perp_09_err = 1.2533 * sd_b_perp_09
med_b_perp_10_err = 1.2533 * sd_b_perp_10
med_b_perp_11_err = 1.2533 * sd_b_perp_11
med_b_perp_16_err = 1.2533 * sd_b_perp_16
med_b_perp_17_err = 1.2533 * sd_b_perp_17
med_b_perp_18_err = 1.2533 * sd_b_perp_18
med_b_perp_19_err = 1.2533 * sd_b_perp_19
med_b_perp_20_err = 1.2533 * sd_b_perp_20
med_b_perp_21_err = 1.2533 * sd_b_perp_21
med_b_perp_24_err = 1.2533 * sd_b_perp_24
med_b_perp_25_err = 1.2533 * sd_b_perp_25
med_b_perp_27_err = 1.2533 * sd_b_perp_27
med_b_perp_28_err = 1.2533 * sd_b_perp_28
med_b_perp_34_err = 1.2533 * sd_b_perp_34
med_b_perp_35_err = 1.2533 * sd_b_perp_35

# high cmls
med_b_perp_01_err = 1.2533 * sd_b_perp_01
med_b_perp_12_err = 1.2533 * sd_b_perp_12
med_b_perp_15_err = 1.2533 * sd_b_perp_15
med_b_perp_26_err = 1.2533 * sd_b_perp_26

#med_b_perp_errs = np.array([med_b_perp_02_err,med_b_perp_03_err,med_b_perp_04_err,med_b_perp_05_err,med_b_perp_08_err,med_b_perp_09_err,med_b_perp_10_err,med_b_perp_11_err,med_b_perp_16_err,med_b_perp_17_err,med_b_perp_18_err,med_b_perp_19_err,med_b_perp_20_err,med_b_perp_21_err,med_b_perp_24_err,med_b_perp_25_err,med_b_perp_27_err,med_b_perp_28_err,med_b_perp_34_err,med_b_perp_35_err])
med_b_perp_errs = np.array([med_b_perp_01_err,med_b_perp_02_err,med_b_perp_03_err,med_b_perp_04_err,med_b_perp_05_err,med_b_perp_08_err,med_b_perp_09_err,med_b_perp_10_err,med_b_perp_11_err,med_b_perp_12_err,med_b_perp_15_err,med_b_perp_16_err,med_b_perp_17_err,med_b_perp_18_err,med_b_perp_19_err,med_b_perp_20_err,med_b_perp_21_err,med_b_perp_24_err,med_b_perp_25_err,med_b_perp_26_err,med_b_perp_27_err,med_b_perp_28_err,med_b_perp_34_err,med_b_perp_35_err])


'''need +/- bperp stuff here'''

b_perp_02_mean_plus = np.mean(b_perp_02_plus)
b_perp_03_mean_plus = np.mean(b_perp_03_plus)
b_perp_04_mean_plus = np.mean(b_perp_04_plus)
b_perp_05_mean_plus = np.mean(b_perp_05_plus)
b_perp_08_mean_plus = np.mean(b_perp_08_plus)
b_perp_09_mean_plus = np.mean(b_perp_09_plus)
b_perp_10_mean_plus = np.mean(b_perp_10_plus)
b_perp_11_mean_plus = np.mean(b_perp_11_plus)
b_perp_16_mean_plus = np.mean(b_perp_16_plus)
b_perp_17_mean_plus = np.mean(b_perp_17_plus)
b_perp_18_mean_plus = np.mean(b_perp_18_plus)
b_perp_19_mean_plus = np.mean(b_perp_19_plus)
b_perp_20_mean_plus = np.mean(b_perp_20_plus)
b_perp_21_mean_plus = np.mean(b_perp_21_plus)
b_perp_24_mean_plus = np.mean(b_perp_24_plus)
b_perp_25_mean_plus = np.mean(b_perp_25_plus)
b_perp_27_mean_plus = np.mean(b_perp_27_plus)
b_perp_28_mean_plus = np.mean(b_perp_28_plus)
b_perp_34_mean_plus = np.mean(b_perp_34_plus)
b_perp_35_mean_plus = np.mean(b_perp_35_plus)

# high cmls
b_perp_01_mean_plus = np.mean(b_perp_01_plus)
b_perp_12_mean_plus = np.mean(b_perp_12_plus)
b_perp_15_mean_plus = np.mean(b_perp_15_plus)
b_perp_26_mean_plus = np.mean(b_perp_26_plus)

#b_perp_means_plus = np.array([b_perp_02_mean_plus,b_perp_03_mean_plus,b_perp_04_mean_plus,b_perp_05_mean_plus,b_perp_08_mean_plus,b_perp_09_mean_plus,b_perp_10_mean_plus,b_perp_11_mean_plus,b_perp_16_mean_plus,b_perp_17_mean_plus,b_perp_18_mean_plus,b_perp_19_mean_plus,b_perp_20_mean_plus,b_perp_21_mean_plus,b_perp_24_mean_plus,b_perp_25_mean_plus,b_perp_27_mean_plus,b_perp_28_mean_plus,b_perp_34_mean_plus,b_perp_35_mean_plus])
b_perp_means_plus = np.array([b_perp_01_mean_plus,b_perp_02_mean_plus,b_perp_03_mean_plus,b_perp_04_mean_plus,b_perp_05_mean_plus,b_perp_08_mean_plus,b_perp_09_mean_plus,b_perp_10_mean_plus,b_perp_11_mean_plus,b_perp_12_mean_plus,b_perp_15_mean_plus,b_perp_16_mean_plus,b_perp_17_mean_plus,b_perp_18_mean_plus,b_perp_19_mean_plus,b_perp_20_mean_plus,b_perp_21_mean_plus,b_perp_24_mean_plus,b_perp_25_mean_plus,b_perp_26_mean_plus,b_perp_27_mean_plus,b_perp_28_mean_plus,b_perp_34_mean_plus,b_perp_35_mean_plus])


sd_b_perp_02_plus = standard_deviation(b_perp_02_plus,b_perp_02_mean_plus)
sd_b_perp_03_plus = standard_deviation(b_perp_03_plus,b_perp_03_mean_plus)
sd_b_perp_04_plus = standard_deviation(b_perp_04_plus,b_perp_04_mean_plus)
sd_b_perp_05_plus = standard_deviation(b_perp_05_plus,b_perp_05_mean_plus)
sd_b_perp_08_plus = standard_deviation(b_perp_08_plus,b_perp_08_mean_plus)
sd_b_perp_09_plus = standard_deviation(b_perp_09_plus,b_perp_09_mean_plus)
sd_b_perp_10_plus = standard_deviation(b_perp_10_plus,b_perp_10_mean_plus)
sd_b_perp_11_plus = standard_deviation(b_perp_11_plus,b_perp_11_mean_plus)
sd_b_perp_16_plus = standard_deviation(b_perp_16_plus,b_perp_16_mean_plus)
sd_b_perp_17_plus = standard_deviation(b_perp_17_plus,b_perp_17_mean_plus)
sd_b_perp_18_plus = standard_deviation(b_perp_18_plus,b_perp_18_mean_plus)
sd_b_perp_19_plus = standard_deviation(b_perp_19_plus,b_perp_19_mean_plus)
sd_b_perp_20_plus = standard_deviation(b_perp_20_plus,b_perp_20_mean_plus)
sd_b_perp_21_plus = standard_deviation(b_perp_21_plus,b_perp_21_mean_plus)
sd_b_perp_24_plus = standard_deviation(b_perp_24_plus,b_perp_24_mean_plus)
sd_b_perp_25_plus = standard_deviation(b_perp_25_plus,b_perp_25_mean_plus)
sd_b_perp_27_plus = standard_deviation(b_perp_27_plus,b_perp_27_mean_plus)
sd_b_perp_28_plus = standard_deviation(b_perp_28_plus,b_perp_28_mean_plus)
sd_b_perp_34_plus = standard_deviation(b_perp_34_plus,b_perp_34_mean_plus)
sd_b_perp_35_plus = standard_deviation(b_perp_35_plus,b_perp_35_mean_plus)

# high cmls
sd_b_perp_01_plus = standard_deviation(b_perp_01_plus,b_perp_01_mean_plus)
sd_b_perp_12_plus = standard_deviation(b_perp_12_plus,b_perp_12_mean_plus)
sd_b_perp_15_plus = standard_deviation(b_perp_15_plus,b_perp_15_mean_plus)
sd_b_perp_26_plus = standard_deviation(b_perp_26_plus,b_perp_26_mean_plus)


b_perp_02_median_plus = np.median(b_perp_02_plus)
b_perp_03_median_plus = np.median(b_perp_03_plus)
b_perp_04_median_plus = np.median(b_perp_04_plus)
b_perp_05_median_plus = np.median(b_perp_05_plus)
b_perp_08_median_plus = np.median(b_perp_08_plus)
b_perp_09_median_plus = np.median(b_perp_09_plus)
b_perp_10_median_plus = np.median(b_perp_10_plus)
b_perp_11_median_plus = np.median(b_perp_11_plus)
b_perp_16_median_plus = np.median(b_perp_16_plus)
b_perp_17_median_plus = np.median(b_perp_17_plus)
b_perp_18_median_plus = np.median(b_perp_18_plus)
b_perp_19_median_plus = np.median(b_perp_19_plus)
b_perp_20_median_plus = np.median(b_perp_20_plus)
b_perp_21_median_plus = np.median(b_perp_21_plus)
b_perp_24_median_plus = np.median(b_perp_24_plus)
b_perp_25_median_plus = np.median(b_perp_25_plus)
b_perp_27_median_plus = np.median(b_perp_27_plus)
b_perp_28_median_plus = np.median(b_perp_28_plus)
b_perp_34_median_plus = np.median(b_perp_34_plus)
b_perp_35_median_plus = np.median(b_perp_25_plus)

# high cmls
b_perp_01_median_plus = np.median(b_perp_01_plus)
b_perp_12_median_plus = np.median(b_perp_12_plus)
b_perp_15_median_plus = np.median(b_perp_15_plus)
b_perp_26_median_plus = np.median(b_perp_26_plus)

#b_perp_medians_plus = np.array([b_perp_02_median_plus,b_perp_03_median_plus,b_perp_04_median_plus,b_perp_05_median_plus,b_perp_08_median_plus,b_perp_09_median_plus,b_perp_10_median_plus,b_perp_11_median_plus,b_perp_16_median_plus,b_perp_17_median_plus,b_perp_18_median_plus,b_perp_19_median_plus,b_perp_20_median_plus,b_perp_21_median_plus,b_perp_24_median_plus,b_perp_25_median_plus,b_perp_27_median_plus,b_perp_28_median_plus,b_perp_34_median_plus,b_perp_35_median_plus])
b_perp_medians_plus = np.array([b_perp_01_median_plus,b_perp_02_median_plus,b_perp_03_median_plus,b_perp_04_median_plus,b_perp_05_median_plus,b_perp_08_median_plus,b_perp_09_median_plus,b_perp_10_median_plus,b_perp_11_median_plus,b_perp_12_median_plus,b_perp_15_median_plus,b_perp_16_median_plus,b_perp_17_median_plus,b_perp_18_median_plus,b_perp_19_median_plus,b_perp_20_median_plus,b_perp_21_median_plus,b_perp_24_median_plus,b_perp_25_median_plus,b_perp_26_median_plus,b_perp_27_median_plus,b_perp_28_median_plus,b_perp_34_median_plus,b_perp_35_median_plus])


med_b_perp_02_err_plus = 1.2533 * sd_b_perp_02_plus
med_b_perp_03_err_plus = 1.2533 * sd_b_perp_03_plus
med_b_perp_04_err_plus = 1.2533 * sd_b_perp_04_plus
med_b_perp_05_err_plus = 1.2533 * sd_b_perp_05_plus
med_b_perp_08_err_plus = 1.2533 * sd_b_perp_08_plus
med_b_perp_09_err_plus = 1.2533 * sd_b_perp_09_plus
med_b_perp_10_err_plus = 1.2533 * sd_b_perp_10_plus
med_b_perp_11_err_plus = 1.2533 * sd_b_perp_11_plus
med_b_perp_16_err_plus = 1.2533 * sd_b_perp_16_plus
med_b_perp_17_err_plus = 1.2533 * sd_b_perp_17_plus
med_b_perp_18_err_plus = 1.2533 * sd_b_perp_18_plus
med_b_perp_19_err_plus = 1.2533 * sd_b_perp_19_plus
med_b_perp_20_err_plus = 1.2533 * sd_b_perp_20_plus
med_b_perp_21_err_plus = 1.2533 * sd_b_perp_21_plus
med_b_perp_24_err_plus = 1.2533 * sd_b_perp_24_plus
med_b_perp_25_err_plus = 1.2533 * sd_b_perp_25_plus
med_b_perp_27_err_plus = 1.2533 * sd_b_perp_27_plus
med_b_perp_28_err_plus = 1.2533 * sd_b_perp_28_plus
med_b_perp_34_err_plus = 1.2533 * sd_b_perp_34_plus
med_b_perp_35_err_plus = 1.2533 * sd_b_perp_35_plus

# high cmls
med_b_perp_01_err_plus = 1.2533 * sd_b_perp_01_plus
med_b_perp_12_err_plus = 1.2533 * sd_b_perp_12_plus
med_b_perp_15_err_plus = 1.2533 * sd_b_perp_15_plus
med_b_perp_26_err_plus = 1.2533 * sd_b_perp_26_plus

#med_b_perp_errs_plus = np.array([med_b_perp_02_err_plus,med_b_perp_03_err_plus,med_b_perp_04_err_plus,med_b_perp_05_err_plus,med_b_perp_08_err_plus,med_b_perp_09_err_plus,med_b_perp_10_err_plus,med_b_perp_11_err_plus,med_b_perp_16_err_plus,med_b_perp_17_err_plus,med_b_perp_18_err_plus,med_b_perp_19_err_plus,med_b_perp_20_err_plus,med_b_perp_21_err_plus,med_b_perp_24_err_plus,med_b_perp_25_err_plus,med_b_perp_27_err_plus,med_b_perp_28_err_plus,med_b_perp_34_err_plus,med_b_perp_35_err_plus])
med_b_perp_errs_plus = np.array([med_b_perp_01_err_plus,med_b_perp_02_err_plus,med_b_perp_03_err_plus,med_b_perp_04_err_plus,med_b_perp_05_err_plus,med_b_perp_08_err_plus,med_b_perp_09_err_plus,med_b_perp_10_err_plus,med_b_perp_11_err_plus,med_b_perp_12_err_plus,med_b_perp_15_err_plus,med_b_perp_16_err_plus,med_b_perp_17_err_plus,med_b_perp_18_err_plus,med_b_perp_19_err_plus,med_b_perp_20_err_plus,med_b_perp_21_err_plus,med_b_perp_24_err_plus,med_b_perp_25_err_plus,med_b_perp_26_err_plus,med_b_perp_27_err_plus,med_b_perp_28_err_plus,med_b_perp_34_err_plus,med_b_perp_35_err_plus])


b_perp_02_mean_minus = np.mean(b_perp_02_minus)
b_perp_03_mean_minus = np.mean(b_perp_03_minus)
b_perp_04_mean_minus = np.mean(b_perp_04_minus)
b_perp_05_mean_minus = np.mean(b_perp_05_minus)
b_perp_08_mean_minus = np.mean(b_perp_08_minus)
b_perp_09_mean_minus = np.mean(b_perp_09_minus)
b_perp_10_mean_minus = np.mean(b_perp_10_minus)
b_perp_11_mean_minus = np.mean(b_perp_11_minus)
b_perp_16_mean_minus = np.mean(b_perp_16_minus)
b_perp_17_mean_minus = np.mean(b_perp_17_minus)
b_perp_18_mean_minus = np.mean(b_perp_18_minus)
b_perp_19_mean_minus = np.mean(b_perp_19_minus)
b_perp_20_mean_minus = np.mean(b_perp_20_minus)
b_perp_21_mean_minus = np.mean(b_perp_21_minus)
b_perp_24_mean_minus = np.mean(b_perp_24_minus)
b_perp_25_mean_minus = np.mean(b_perp_25_minus)
b_perp_27_mean_minus = np.mean(b_perp_27_minus)
b_perp_28_mean_minus = np.mean(b_perp_28_minus)
b_perp_34_mean_minus = np.mean(b_perp_34_minus)
b_perp_35_mean_minus = np.mean(b_perp_35_minus)

# high cmls
b_perp_01_mean_minus = np.mean(b_perp_01_minus)
b_perp_12_mean_minus = np.mean(b_perp_12_minus)
b_perp_15_mean_minus = np.mean(b_perp_15_minus)
b_perp_26_mean_minus = np.mean(b_perp_26_minus)

#b_perp_means_minus = np.array([b_perp_02_mean_minus,b_perp_03_mean_minus,b_perp_04_mean_minus,b_perp_05_mean_minus,b_perp_08_mean_minus,b_perp_09_mean_minus,b_perp_10_mean_minus,b_perp_11_mean_minus,b_perp_16_mean_minus,b_perp_17_mean_minus,b_perp_18_mean_minus,b_perp_19_mean_minus,b_perp_20_mean_minus,b_perp_21_mean_minus,b_perp_24_mean_minus,b_perp_25_mean_minus,b_perp_27_mean_minus,b_perp_28_mean_minus,b_perp_34_mean_minus,b_perp_35_mean_minus])
b_perp_means_minus = np.array([b_perp_01_mean_minus,b_perp_02_mean_minus,b_perp_03_mean_minus,b_perp_04_mean_minus,b_perp_05_mean_minus,b_perp_08_mean_minus,b_perp_09_mean_minus,b_perp_10_mean_minus,b_perp_11_mean_minus,b_perp_12_mean_minus,b_perp_15_mean_minus,b_perp_16_mean_minus,b_perp_17_mean_minus,b_perp_18_mean_minus,b_perp_19_mean_minus,b_perp_20_mean_minus,b_perp_21_mean_minus,b_perp_24_mean_minus,b_perp_25_mean_minus,b_perp_26_mean_minus,b_perp_27_mean_minus,b_perp_28_mean_minus,b_perp_34_mean_minus,b_perp_35_mean_minus])


sd_b_perp_02_minus = standard_deviation(b_perp_02_minus,b_perp_02_mean_minus)
sd_b_perp_03_minus = standard_deviation(b_perp_03_minus,b_perp_03_mean_minus)
sd_b_perp_04_minus = standard_deviation(b_perp_04_minus,b_perp_04_mean_minus)
sd_b_perp_05_minus = standard_deviation(b_perp_05_minus,b_perp_05_mean_minus)
sd_b_perp_08_minus = standard_deviation(b_perp_08_minus,b_perp_08_mean_minus)
sd_b_perp_09_minus = standard_deviation(b_perp_09_minus,b_perp_09_mean_minus)
sd_b_perp_10_minus = standard_deviation(b_perp_10_minus,b_perp_10_mean_minus)
sd_b_perp_11_minus = standard_deviation(b_perp_11_minus,b_perp_11_mean_minus)
sd_b_perp_16_minus = standard_deviation(b_perp_16_minus,b_perp_16_mean_minus)
sd_b_perp_17_minus = standard_deviation(b_perp_17_minus,b_perp_17_mean_minus)
sd_b_perp_18_minus = standard_deviation(b_perp_18_minus,b_perp_18_mean_minus)
sd_b_perp_19_minus = standard_deviation(b_perp_19_minus,b_perp_19_mean_minus)
sd_b_perp_20_minus = standard_deviation(b_perp_20_minus,b_perp_20_mean_minus)
sd_b_perp_21_minus = standard_deviation(b_perp_21_minus,b_perp_21_mean_minus)
sd_b_perp_24_minus = standard_deviation(b_perp_24_minus,b_perp_24_mean_minus)
sd_b_perp_25_minus = standard_deviation(b_perp_25_minus,b_perp_25_mean_minus)
sd_b_perp_27_minus = standard_deviation(b_perp_27_minus,b_perp_27_mean_minus)
sd_b_perp_28_minus = standard_deviation(b_perp_28_minus,b_perp_28_mean_minus)
sd_b_perp_34_minus = standard_deviation(b_perp_34_minus,b_perp_34_mean_minus)
sd_b_perp_35_minus = standard_deviation(b_perp_35_minus,b_perp_35_mean_minus)

# high cmls
sd_b_perp_01_minus = standard_deviation(b_perp_01_minus,b_perp_01_mean_minus)
sd_b_perp_12_minus = standard_deviation(b_perp_12_minus,b_perp_12_mean_minus)
sd_b_perp_15_minus = standard_deviation(b_perp_15_minus,b_perp_15_mean_minus)
sd_b_perp_26_minus = standard_deviation(b_perp_26_minus,b_perp_26_mean_minus)


b_perp_02_median_minus = np.median(b_perp_02_minus)
b_perp_03_median_minus = np.median(b_perp_03_minus)
b_perp_04_median_minus = np.median(b_perp_04_minus)
b_perp_05_median_minus = np.median(b_perp_05_minus)
b_perp_08_median_minus = np.median(b_perp_08_minus)
b_perp_09_median_minus = np.median(b_perp_09_minus)
b_perp_10_median_minus = np.median(b_perp_10_minus)
b_perp_11_median_minus = np.median(b_perp_11_minus)
b_perp_16_median_minus = np.median(b_perp_16_minus)
b_perp_17_median_minus = np.median(b_perp_17_minus)
b_perp_18_median_minus = np.median(b_perp_18_minus)
b_perp_19_median_minus = np.median(b_perp_19_minus)
b_perp_20_median_minus = np.median(b_perp_20_minus)
b_perp_21_median_minus = np.median(b_perp_21_minus)
b_perp_24_median_minus = np.median(b_perp_24_minus)
b_perp_25_median_minus = np.median(b_perp_25_minus)
b_perp_27_median_minus = np.median(b_perp_27_minus)
b_perp_28_median_minus = np.median(b_perp_28_minus)
b_perp_34_median_minus = np.median(b_perp_34_minus)
b_perp_35_median_minus = np.median(b_perp_35_minus)

# high cmls
b_perp_01_median_minus = np.median(b_perp_01_minus)
b_perp_12_median_minus = np.median(b_perp_12_minus)
b_perp_15_median_minus = np.median(b_perp_15_minus)
b_perp_26_median_minus = np.median(b_perp_26_minus)

#b_perp_medians_minus = np.array([b_perp_02_median_minus,b_perp_03_median_minus,b_perp_04_median_minus,b_perp_05_median_minus,b_perp_08_median_minus,b_perp_09_median_minus,b_perp_10_median_minus,b_perp_11_median_minus,b_perp_16_median_minus,b_perp_17_median_minus,b_perp_18_median_minus,b_perp_19_median_minus,b_perp_20_median_minus,b_perp_21_median_minus,b_perp_24_median_minus,b_perp_25_median_minus,b_perp_27_median_minus,b_perp_28_median_minus,b_perp_34_median_minus,b_perp_35_median_minus])
b_perp_medians_minus = np.array([b_perp_01_median_minus,b_perp_02_median_minus,b_perp_03_median_minus,b_perp_04_median_minus,b_perp_05_median_minus,b_perp_08_median_minus,b_perp_09_median_minus,b_perp_10_median_minus,b_perp_11_median_minus,b_perp_12_median_minus,b_perp_15_median_minus,b_perp_16_median_minus,b_perp_17_median_minus,b_perp_18_median_minus,b_perp_19_median_minus,b_perp_20_median_minus,b_perp_21_median_minus,b_perp_24_median_minus,b_perp_25_median_minus,b_perp_26_median_minus,b_perp_27_median_minus,b_perp_28_median_minus,b_perp_34_median_minus,b_perp_35_median_minus])


med_b_perp_02_err_minus = 1.2533 * sd_b_perp_02_minus
med_b_perp_03_err_minus = 1.2533 * sd_b_perp_03_minus
med_b_perp_04_err_minus = 1.2533 * sd_b_perp_04_minus
med_b_perp_05_err_minus = 1.2533 * sd_b_perp_05_minus
med_b_perp_08_err_minus = 1.2533 * sd_b_perp_08_minus
med_b_perp_09_err_minus = 1.2533 * sd_b_perp_09_minus
med_b_perp_10_err_minus = 1.2533 * sd_b_perp_10_minus
med_b_perp_11_err_minus = 1.2533 * sd_b_perp_11_minus
med_b_perp_16_err_minus = 1.2533 * sd_b_perp_16_minus
med_b_perp_17_err_minus = 1.2533 * sd_b_perp_17_minus
med_b_perp_18_err_minus = 1.2533 * sd_b_perp_18_minus
med_b_perp_19_err_minus = 1.2533 * sd_b_perp_19_minus
med_b_perp_20_err_minus = 1.2533 * sd_b_perp_20_minus
med_b_perp_21_err_minus = 1.2533 * sd_b_perp_21_minus
med_b_perp_24_err_minus = 1.2533 * sd_b_perp_24_minus
med_b_perp_25_err_minus = 1.2533 * sd_b_perp_25_minus
med_b_perp_27_err_minus = 1.2533 * sd_b_perp_27_minus
med_b_perp_28_err_minus = 1.2533 * sd_b_perp_28_minus
med_b_perp_34_err_minus = 1.2533 * sd_b_perp_34_minus
med_b_perp_35_err_minus = 1.2533 * sd_b_perp_35_minus

# high cmls
med_b_perp_01_err_minus = 1.2533 * sd_b_perp_01_minus
med_b_perp_12_err_minus = 1.2533 * sd_b_perp_12_minus
med_b_perp_15_err_minus = 1.2533 * sd_b_perp_15_minus
med_b_perp_26_err_minus = 1.2533 * sd_b_perp_26_minus

#med_b_perp_errs_minus = np.array([med_b_perp_02_err_minus,med_b_perp_03_err_minus,med_b_perp_04_err_minus,med_b_perp_05_err_minus,med_b_perp_08_err_minus,med_b_perp_09_err_minus,med_b_perp_10_err_minus,med_b_perp_11_err_minus,med_b_perp_16_err_minus,med_b_perp_17_err_minus,med_b_perp_18_err_minus,med_b_perp_19_err_minus,med_b_perp_20_err_minus,med_b_perp_21_err_minus,med_b_perp_24_err_minus,med_b_perp_25_err_minus,med_b_perp_27_err_minus,med_b_perp_28_err_minus,med_b_perp_34_err_minus,med_b_perp_35_err_minus])
med_b_perp_errs_minus = np.array([med_b_perp_01_err_minus,med_b_perp_02_err_minus,med_b_perp_03_err_minus,med_b_perp_04_err_minus,med_b_perp_05_err_minus,med_b_perp_08_err_minus,med_b_perp_09_err_minus,med_b_perp_10_err_minus,med_b_perp_11_err_minus,med_b_perp_12_err_minus,med_b_perp_15_err_minus,med_b_perp_16_err_minus,med_b_perp_17_err_minus,med_b_perp_18_err_minus,med_b_perp_19_err_minus,med_b_perp_20_err_minus,med_b_perp_21_err_minus,med_b_perp_24_err_minus,med_b_perp_25_err_minus,med_b_perp_26_err_minus,med_b_perp_27_err_minus,med_b_perp_28_err_minus,med_b_perp_34_err_minus,med_b_perp_35_err_minus])


'''
AVERAGE B TOTAL
'''

btot_02_mean = np.mean(btot_02)
btot_03_mean = np.mean(btot_03)
btot_04_mean = np.mean(btot_04)
btot_05_mean = np.mean(btot_05)
btot_08_mean = np.mean(btot_08)
btot_09_mean = np.mean(btot_09)
btot_10_mean = np.mean(btot_10)
btot_11_mean = np.mean(btot_11)
btot_16_mean = np.mean(btot_16)
btot_17_mean = np.mean(btot_17)
btot_18_mean = np.mean(btot_18)
btot_19_mean = np.mean(btot_19)
btot_20_mean = np.mean(btot_20)
btot_21_mean = np.mean(btot_21)
btot_24_mean = np.mean(btot_24)
btot_25_mean = np.mean(btot_25)
btot_27_mean = np.mean(btot_27)
btot_28_mean = np.mean(btot_28)
btot_34_mean = np.mean(btot_34)
btot_35_mean = np.mean(btot_35)

# high cmls
btot_01_mean = np.mean(btot_01)
btot_12_mean = np.mean(btot_12)
btot_15_mean = np.mean(btot_15)
btot_26_mean = np.mean(btot_26)

#btot_means = np.array([btot_02_mean,btot_03_mean,btot_04_mean,btot_05_mean,btot_08_mean,btot_09_mean,btot_10_mean,btot_11_mean,btot_16_mean,btot_17_mean,btot_18_mean,btot_19_mean,btot_20_mean,btot_21_mean,btot_24_mean,btot_25_mean,btot_27_mean,btot_28_mean,btot_34_mean,btot_35_mean])
btot_means = np.array([btot_01_mean,btot_02_mean,btot_03_mean,btot_04_mean,btot_05_mean,btot_08_mean,btot_09_mean,btot_10_mean,btot_11_mean,btot_12_mean,btot_15_mean,btot_16_mean,btot_17_mean,btot_18_mean,btot_19_mean,btot_20_mean,btot_21_mean,btot_24_mean,btot_25_mean,btot_26_mean,btot_27_mean,btot_28_mean,btot_34_mean,btot_35_mean])


sd_btot_02 = standard_deviation(btot_02,btot_02_mean)
sd_btot_03 = standard_deviation(btot_03,btot_03_mean)
sd_btot_04 = standard_deviation(btot_04,btot_04_mean)
sd_btot_05 = standard_deviation(btot_05,btot_05_mean)
sd_btot_08 = standard_deviation(btot_08,btot_08_mean)
sd_btot_09 = standard_deviation(btot_09,btot_09_mean)
sd_btot_10 = standard_deviation(btot_10,btot_10_mean)
sd_btot_11 = standard_deviation(btot_11,btot_11_mean)
sd_btot_16 = standard_deviation(btot_16,btot_16_mean)
sd_btot_17 = standard_deviation(btot_17,btot_17_mean)
sd_btot_18 = standard_deviation(btot_18,btot_18_mean)
sd_btot_19 = standard_deviation(btot_19,btot_19_mean)
sd_btot_20 = standard_deviation(btot_20,btot_20_mean)
sd_btot_21 = standard_deviation(btot_21,btot_21_mean)
sd_btot_24 = standard_deviation(btot_24,btot_24_mean)
sd_btot_25 = standard_deviation(btot_25,btot_25_mean)
sd_btot_27 = standard_deviation(btot_27,btot_27_mean)
sd_btot_28 = standard_deviation(btot_28,btot_28_mean)
sd_btot_34 = standard_deviation(btot_34,btot_34_mean)
sd_btot_35 = standard_deviation(btot_35,btot_35_mean)

# high cmls
sd_btot_01 = standard_deviation(btot_01,btot_01_mean)
sd_btot_12 = standard_deviation(btot_12,btot_12_mean)
sd_btot_15 = standard_deviation(btot_15,btot_15_mean)
sd_btot_26 = standard_deviation(btot_26,btot_26_mean)


btot_02_median = np.median(btot_02)
btot_03_median = np.median(btot_03)
btot_04_median = np.median(btot_04)
btot_05_median = np.median(btot_05)
btot_08_median = np.median(btot_08)
btot_09_median = np.median(btot_09)
btot_10_median = np.median(btot_10)
btot_11_median = np.median(btot_11)
btot_16_median = np.median(btot_16)
btot_17_median = np.median(btot_17)
btot_18_median = np.median(btot_18)
btot_19_median = np.median(btot_19)
btot_20_median = np.median(btot_20)
btot_21_median = np.median(btot_21)
btot_24_median = np.median(btot_24)
btot_25_median = np.median(btot_25)
btot_27_median = np.median(btot_27)
btot_28_median = np.median(btot_28)
btot_34_median = np.median(btot_34)
btot_35_median = np.median(btot_35)

# high cmls
btot_01_median = np.median(btot_01)
btot_12_median = np.median(btot_12)
btot_15_median = np.median(btot_15)
btot_26_median = np.median(btot_26)

#btot_medians = np.array([btot_02_median,btot_03_median,btot_04_median,btot_05_median,btot_08_median,btot_09_median,btot_10_median,btot_11_median,btot_16_median,btot_17_median,btot_18_median,btot_19_median,btot_20_median,btot_21_median,btot_24_median,btot_25_median,btot_27_median,btot_28_median,btot_34_median,btot_35_median])
btot_medians = np.array([btot_01_median,btot_02_median,btot_03_median,btot_04_median,btot_05_median,btot_08_median,btot_09_median,btot_10_median,btot_11_median,btot_12_median,btot_15_median,btot_16_median,btot_17_median,btot_18_median,btot_19_median,btot_20_median,btot_21_median,btot_24_median,btot_25_median,btot_26_median,btot_27_median,btot_28_median,btot_34_median,btot_35_median])


btot_med_sd_err_02 = 1.2533 * sd_btot_02
btot_med_sd_err_03 = 1.2533 * sd_btot_03
btot_med_sd_err_04 = 1.2533 * sd_btot_04
btot_med_sd_err_05 = 1.2533 * sd_btot_05
btot_med_sd_err_08 = 1.2533 * sd_btot_08
btot_med_sd_err_09 = 1.2533 * sd_btot_09
btot_med_sd_err_10 = 1.2533 * sd_btot_10
btot_med_sd_err_11 = 1.2533 * sd_btot_11
btot_med_sd_err_16 = 1.2533 * sd_btot_16
btot_med_sd_err_17 = 1.2533 * sd_btot_17
btot_med_sd_err_18 = 1.2533 * sd_btot_18
btot_med_sd_err_19 = 1.2533 * sd_btot_19
btot_med_sd_err_20 = 1.2533 * sd_btot_20
btot_med_sd_err_21 = 1.2533 * sd_btot_21
btot_med_sd_err_24 = 1.2533 * sd_btot_24
btot_med_sd_err_25 = 1.2533 * sd_btot_25
btot_med_sd_err_27 = 1.2533 * sd_btot_27
btot_med_sd_err_28 = 1.2533 * sd_btot_28
btot_med_sd_err_34 = 1.2533 * sd_btot_34
btot_med_sd_err_35 = 1.2533 * sd_btot_35

# high cmls
btot_med_sd_err_01 = 1.2533 * sd_btot_01
btot_med_sd_err_12 = 1.2533 * sd_btot_12
btot_med_sd_err_15 = 1.2533 * sd_btot_15
btot_med_sd_err_26 = 1.2533 * sd_btot_26

#btot_med_sd_errs = np.array([btot_med_sd_err_02,btot_med_sd_err_03,btot_med_sd_err_04,btot_med_sd_err_05,btot_med_sd_err_08,btot_med_sd_err_09,btot_med_sd_err_10,btot_med_sd_err_11,btot_med_sd_err_16,btot_med_sd_err_17,btot_med_sd_err_18,btot_med_sd_err_19,btot_med_sd_err_20,btot_med_sd_err_21,btot_med_sd_err_24,btot_med_sd_err_25,btot_med_sd_err_27,btot_med_sd_err_28,btot_med_sd_err_34,btot_med_sd_err_35])
btot_med_sd_errs = np.array([btot_med_sd_err_01,btot_med_sd_err_02,btot_med_sd_err_03,btot_med_sd_err_04,btot_med_sd_err_05,btot_med_sd_err_08,btot_med_sd_err_09,btot_med_sd_err_10,btot_med_sd_err_11,btot_med_sd_err_12,btot_med_sd_err_15,btot_med_sd_err_16,btot_med_sd_err_17,btot_med_sd_err_18,btot_med_sd_err_19,btot_med_sd_err_20,btot_med_sd_err_21,btot_med_sd_err_24,btot_med_sd_err_25,btot_med_sd_err_26,btot_med_sd_err_27,btot_med_sd_err_28,btot_med_sd_err_34,btot_med_sd_err_35])


# ---- + % ------

btot_02_mean_p = np.mean(btot_02_p)
btot_03_mean_p = np.mean(btot_03_p)
btot_04_mean_p = np.mean(btot_04_p)
btot_05_mean_p = np.mean(btot_05_p)
btot_08_mean_p = np.mean(btot_08_p)
btot_09_mean_p = np.mean(btot_09_p)
btot_10_mean_p = np.mean(btot_10_p)
btot_11_mean_p = np.mean(btot_11_p)
btot_16_mean_p = np.mean(btot_16_p)
btot_17_mean_p = np.mean(btot_17_p)
btot_18_mean_p = np.mean(btot_18_p)
btot_19_mean_p = np.mean(btot_19_p)
btot_20_mean_p = np.mean(btot_20_p)
btot_21_mean_p = np.mean(btot_21_p)
btot_24_mean_p = np.mean(btot_24_p)
btot_25_mean_p = np.mean(btot_25_p)
btot_27_mean_p = np.mean(btot_27_p)
btot_28_mean_p = np.mean(btot_28_p)
btot_34_mean_p = np.mean(btot_34_p)
btot_35_mean_p = np.mean(btot_35_p)

# higher cmls
btot_01_mean_p = np.mean(btot_01_p)
btot_12_mean_p = np.mean(btot_12_p)
btot_15_mean_p = np.mean(btot_15_p)
btot_26_mean_p = np.mean(btot_26_p)

#btot_means_p = np.array([btot_02_mean_p,btot_03_mean_p,btot_04_mean_p,btot_05_mean_p,btot_08_mean_p,btot_09_mean_p,btot_10_mean_p,btot_11_mean_p,btot_16_mean,btot_17_mean_p,btot_18_mean_p,btot_19_mean_p,btot_20_mean_p,btot_21_mean_p,btot_24_mean,btot_25_mean_p,btot_27_mean_p,btot_28_mean_p,btot_34_mean_p,btot_35_mean_p])
btot_means_p = np.array([btot_01_mean_p,btot_02_mean_p,btot_03_mean_p,btot_04_mean_p,btot_05_mean_p,btot_08_mean_p,btot_09_mean_p,btot_10_mean_p,btot_11_mean_p,btot_12_mean_p,btot_15_mean_p,btot_16_mean,btot_17_mean_p,btot_18_mean_p,btot_19_mean_p,btot_20_mean_p,btot_21_mean_p,btot_24_mean,btot_25_mean_p,btot_26_mean_p,btot_27_mean_p,btot_28_mean_p,btot_34_mean_p,btot_35_mean_p])


sd_btot_02_p = standard_deviation(btot_02_p,btot_02_mean_p)
sd_btot_03_p = standard_deviation(btot_03_p,btot_03_mean_p)
sd_btot_04_p = standard_deviation(btot_04_p,btot_04_mean_p)
sd_btot_05_p = standard_deviation(btot_05_p,btot_05_mean_p)
sd_btot_08_p = standard_deviation(btot_08_p,btot_08_mean_p)
sd_btot_09_p = standard_deviation(btot_09_p,btot_09_mean_p)
sd_btot_10_p = standard_deviation(btot_10_p,btot_10_mean_p)
sd_btot_11_p = standard_deviation(btot_11_p,btot_11_mean_p)
sd_btot_16_p = standard_deviation(btot_16_p,btot_16_mean_p)
sd_btot_17_p = standard_deviation(btot_17_p,btot_17_mean_p)
sd_btot_18_p = standard_deviation(btot_18_p,btot_18_mean_p)
sd_btot_19_p = standard_deviation(btot_19_p,btot_19_mean_p)
sd_btot_20_p = standard_deviation(btot_20_p,btot_20_mean_p)
sd_btot_21_p = standard_deviation(btot_21_p,btot_21_mean_p)
sd_btot_24_p = standard_deviation(btot_24_p,btot_24_mean_p)
sd_btot_25_p = standard_deviation(btot_25_p,btot_25_mean_p)
sd_btot_27_p = standard_deviation(btot_27_p,btot_27_mean_p)
sd_btot_28_p = standard_deviation(btot_28_p,btot_28_mean_p)
sd_btot_34_p = standard_deviation(btot_34_p,btot_34_mean_p)
sd_btot_35_p = standard_deviation(btot_35_p,btot_35_mean_p)

# higher cmls
sd_btot_01_p = standard_deviation(btot_01_p,btot_01_mean_p)
sd_btot_12_p = standard_deviation(btot_12_p,btot_12_mean_p)
sd_btot_15_p = standard_deviation(btot_15_p,btot_15_mean_p)
sd_btot_26_p = standard_deviation(btot_26_p,btot_26_mean_p)


btot_02_median_p = np.median(btot_02_p)
btot_03_median_p = np.median(btot_03_p)
btot_04_median_p = np.median(btot_04_p)
btot_05_median_p = np.median(btot_05_p)
btot_08_median_p = np.median(btot_08_p)
btot_09_median_p = np.median(btot_09_p)
btot_10_median_p = np.median(btot_10_p)
btot_11_median_p = np.median(btot_11_p)
btot_16_median_p = np.median(btot_16_p)
btot_17_median_p = np.median(btot_17_p)
btot_18_median_p = np.median(btot_18_p)
btot_19_median_p = np.median(btot_19_p)
btot_20_median_p = np.median(btot_20_p)
btot_21_median_p = np.median(btot_21_p)
btot_24_median_p = np.median(btot_24_p)
btot_25_median_p = np.median(btot_25_p)
btot_27_median_p = np.median(btot_27_p)
btot_28_median_p = np.median(btot_28_p)
btot_34_median_p = np.median(btot_34_p)
btot_35_median_p = np.median(btot_35_p)

# higher cmls
btot_01_median_p = np.median(btot_01_p)
btot_12_median_p = np.median(btot_12_p)
btot_15_median_p = np.median(btot_15_p)
btot_26_median_p = np.median(btot_26_p)

#btot_medians_p = np.array([btot_02_median_p,btot_03_median_p,btot_04_median_p,btot_05_median_p,btot_08_median_p,btot_09_median_p,btot_10_median_p,btot_11_median_p,btot_16_median_p,btot_17_median_p,btot_18_median_p,btot_19_median_p,btot_20_median_p,btot_21_median_p,btot_24_median_p,btot_25_median_p,btot_27_median_p,btot_28_median_p,btot_34_median_p,btot_35_median_p])
btot_medians_p = np.array([btot_01_median_p,btot_02_median_p,btot_03_median_p,btot_04_median_p,btot_05_median_p,btot_08_median_p,btot_09_median_p,btot_10_median_p,btot_11_median_p,btot_12_median_p,btot_15_median_p,btot_16_median_p,btot_17_median_p,btot_18_median_p,btot_19_median_p,btot_20_median_p,btot_21_median_p,btot_24_median_p,btot_25_median_p,btot_26_median_p,btot_27_median_p,btot_28_median_p,btot_34_median_p,btot_35_median_p])


btot_med_sd_err_02_p = 1.2533 * sd_btot_02_p
btot_med_sd_err_03_p = 1.2533 * sd_btot_03_p
btot_med_sd_err_04_p = 1.2533 * sd_btot_04_p
btot_med_sd_err_05_p = 1.2533 * sd_btot_05_p
btot_med_sd_err_08_p = 1.2533 * sd_btot_08_p
btot_med_sd_err_09_p = 1.2533 * sd_btot_09_p
btot_med_sd_err_10_p = 1.2533 * sd_btot_10_p
btot_med_sd_err_11_p = 1.2533 * sd_btot_11_p
btot_med_sd_err_16_p = 1.2533 * sd_btot_16_p
btot_med_sd_err_17_p = 1.2533 * sd_btot_17_p
btot_med_sd_err_18_p = 1.2533 * sd_btot_18_p
btot_med_sd_err_19_p = 1.2533 * sd_btot_19_p
btot_med_sd_err_20_p = 1.2533 * sd_btot_20_p
btot_med_sd_err_21_p = 1.2533 * sd_btot_21_p
btot_med_sd_err_24_p = 1.2533 * sd_btot_24_p
btot_med_sd_err_25_p = 1.2533 * sd_btot_25_p
btot_med_sd_err_27_p = 1.2533 * sd_btot_27_p
btot_med_sd_err_28_p = 1.2533 * sd_btot_28_p
btot_med_sd_err_34_p = 1.2533 * sd_btot_34_p
btot_med_sd_err_35_p = 1.2533 * sd_btot_35_p

# higher cmls
btot_med_sd_err_01_p = 1.2533 * sd_btot_01_p
btot_med_sd_err_12_p = 1.2533 * sd_btot_12_p
btot_med_sd_err_15_p = 1.2533 * sd_btot_15_p
btot_med_sd_err_26_p = 1.2533 * sd_btot_26_p

#btot_med_sd_errs_p = np.array([btot_med_sd_err_02_p,btot_med_sd_err_03_p,btot_med_sd_err_04_p,btot_med_sd_err_05_p,btot_med_sd_err_08_p,btot_med_sd_err_09_p,btot_med_sd_err_10_p,btot_med_sd_err_11_p,btot_med_sd_err_16_p,btot_med_sd_err_17_p,btot_med_sd_err_18_p,btot_med_sd_err_19_p,btot_med_sd_err_20_p,btot_med_sd_err_21_p,btot_med_sd_err_24_p,btot_med_sd_err_25_p,btot_med_sd_err_27_p,btot_med_sd_err_28_p,btot_med_sd_err_34_p,btot_med_sd_err_35_p])
btot_med_sd_errs_p = np.array([btot_med_sd_err_01_p,btot_med_sd_err_02_p,btot_med_sd_err_03_p,btot_med_sd_err_04_p,btot_med_sd_err_05_p,btot_med_sd_err_08_p,btot_med_sd_err_09_p,btot_med_sd_err_10_p,btot_med_sd_err_11_p,btot_med_sd_err_12_p,btot_med_sd_err_15_p,btot_med_sd_err_16_p,btot_med_sd_err_17_p,btot_med_sd_err_18_p,btot_med_sd_err_19_p,btot_med_sd_err_20_p,btot_med_sd_err_21_p,btot_med_sd_err_24_p,btot_med_sd_err_25_p,btot_med_sd_err_26_p,btot_med_sd_err_27_p,btot_med_sd_err_28_p,btot_med_sd_err_34_p,btot_med_sd_err_35_p])


# ------- -% ----------

btot_02_mean_m = np.mean(btot_02_m)
btot_03_mean_m = np.mean(btot_03_m)
btot_04_mean_m = np.mean(btot_04_m)
btot_05_mean_m = np.mean(btot_05_m)
btot_08_mean_m = np.mean(btot_08_m)
btot_09_mean_m = np.mean(btot_09_m)
btot_10_mean_m = np.mean(btot_10_m)
btot_11_mean_m = np.mean(btot_11_m)
btot_16_mean_m = np.mean(btot_16_m)
btot_17_mean_m = np.mean(btot_17_m)
btot_18_mean_m = np.mean(btot_18_m)
btot_19_mean_m = np.mean(btot_19_m)
btot_20_mean_m = np.mean(btot_20_m)
btot_21_mean_m = np.mean(btot_21_m)
btot_24_mean_m = np.mean(btot_24_m)
btot_25_mean_m = np.mean(btot_25_m)
btot_27_mean_m = np.mean(btot_27_m)
btot_28_mean_m = np.mean(btot_28_m)
btot_34_mean_m = np.mean(btot_34_m)
btot_35_mean_m = np.mean(btot_35_m)

# higher cmls
btot_01_mean_m = np.mean(btot_01_m)
btot_12_mean_m = np.mean(btot_12_m)
btot_15_mean_m = np.mean(btot_15_m)
btot_26_mean_m = np.mean(btot_26_m)

#btot_means_m = np.array([btot_02_mean_m,btot_03_mean_m,btot_04_mean_m,btot_05_mean_m,btot_08_mean_m,btot_09_mean_m,btot_10_mean_m,btot_11_mean_m,btot_16_mean_m,btot_17_mean_m,btot_18_mean_m,btot_19_mean_m,btot_20_mean_m,btot_21_mean_m,btot_24_mean_m,btot_25_mean_m,btot_27_mean_m,btot_28_mean_m,btot_34_mean_m,btot_35_mean_m])
btot_means_m = np.array([btot_01_mean_m,btot_02_mean_m,btot_03_mean_m,btot_04_mean_m,btot_05_mean_m,btot_08_mean_m,btot_09_mean_m,btot_10_mean_m,btot_11_mean_m,btot_12_mean_m,btot_15_mean_m,btot_16_mean_m,btot_17_mean_m,btot_18_mean_m,btot_19_mean_m,btot_20_mean_m,btot_21_mean_m,btot_24_mean_m,btot_25_mean_m,btot_26_mean_m,btot_27_mean_m,btot_28_mean_m,btot_34_mean_m,btot_35_mean_m])


sd_btot_02_m = standard_deviation(btot_02_m,btot_02_mean_m)
sd_btot_03_m = standard_deviation(btot_03_m,btot_03_mean_m)
sd_btot_04_m = standard_deviation(btot_04_m,btot_04_mean_m)
sd_btot_05_m = standard_deviation(btot_05_m,btot_05_mean_m)
sd_btot_08_m = standard_deviation(btot_08_m,btot_08_mean_m)
sd_btot_09_m = standard_deviation(btot_09_m,btot_09_mean_m)
sd_btot_10_m = standard_deviation(btot_10_m,btot_10_mean_m)
sd_btot_11_m = standard_deviation(btot_11_m,btot_11_mean_m)
sd_btot_16_m = standard_deviation(btot_16_m,btot_16_mean_m)
sd_btot_17_m = standard_deviation(btot_17_m,btot_17_mean_m)
sd_btot_18_m = standard_deviation(btot_18_m,btot_18_mean_m)
sd_btot_19_m = standard_deviation(btot_19_m,btot_19_mean_m)
sd_btot_20_m = standard_deviation(btot_20_m,btot_20_mean_m)
sd_btot_21_m = standard_deviation(btot_21_m,btot_21_mean_m)
sd_btot_24_m = standard_deviation(btot_24_m,btot_24_mean_m)
sd_btot_25_m = standard_deviation(btot_25_m,btot_25_mean_m)
sd_btot_27_m = standard_deviation(btot_27_m,btot_27_mean_m)
sd_btot_28_m = standard_deviation(btot_28_m,btot_28_mean_m)
sd_btot_34_m = standard_deviation(btot_34_m,btot_34_mean_m)
sd_btot_35_m = standard_deviation(btot_35_m,btot_35_mean_m)

# higher cmls
sd_btot_01_m = standard_deviation(btot_01_m,btot_01_mean_m)
sd_btot_12_m = standard_deviation(btot_12_m,btot_12_mean_m)
sd_btot_15_m = standard_deviation(btot_15_m,btot_15_mean_m)
sd_btot_26_m = standard_deviation(btot_26_m,btot_26_mean_m)


btot_02_median_m = np.median(btot_02_m)
btot_03_median_m = np.median(btot_03_m)
btot_04_median_m = np.median(btot_04_m)
btot_05_median_m = np.median(btot_05_m)
btot_08_median_m = np.median(btot_08_m)
btot_09_median_m = np.median(btot_09_m)
btot_10_median_m = np.median(btot_10_m)
btot_11_median_m = np.median(btot_11_m)
btot_16_median_m = np.median(btot_16_m)
btot_17_median_m = np.median(btot_17_m)
btot_18_median_m = np.median(btot_18_m)
btot_19_median_m = np.median(btot_19_m)
btot_20_median_m = np.median(btot_20_m)
btot_21_median_m = np.median(btot_21_m)
btot_24_median_m = np.median(btot_24_m)
btot_25_median_m = np.median(btot_25_m)
btot_27_median_m = np.median(btot_27_m)
btot_28_median_m = np.median(btot_28_m)
btot_34_median_m = np.median(btot_34_m)
btot_35_median_m = np.median(btot_35_m)

# higher cmls
btot_01_median_m = np.median(btot_01_m)
btot_12_median_m = np.median(btot_12_m)
btot_15_median_m = np.median(btot_15_m)
btot_26_median_m = np.median(btot_26_m)

#btot_medians_m = np.array([btot_02_median_m,btot_03_median_m,btot_04_median_m,btot_05_median_m,btot_08_median_m,btot_09_median_m,btot_10_median_m,btot_11_median_m,btot_16_median_m,btot_17_median_m,btot_18_median_m,btot_19_median_m,btot_20_median_m,btot_21_median_m,btot_24_median_m,btot_25_median_m,btot_27_median_m,btot_28_median_m,btot_34_median_m,btot_35_median_m])
btot_medians_m = np.array([btot_01_median_m,btot_02_median_m,btot_03_median_m,btot_04_median_m,btot_05_median_m,btot_08_median_m,btot_09_median_m,btot_10_median_m,btot_11_median_m,btot_12_median_m,btot_15_median_m,btot_16_median_m,btot_17_median_m,btot_18_median_m,btot_19_median_m,btot_20_median_m,btot_21_median_m,btot_24_median_m,btot_25_median_m,btot_26_median_m,btot_27_median_m,btot_28_median_m,btot_34_median_m,btot_35_median_m])


btot_med_sd_err_02_m = 1.2533 * sd_btot_02_m
btot_med_sd_err_03_m = 1.2533 * sd_btot_03_m
btot_med_sd_err_04_m = 1.2533 * sd_btot_04_m
btot_med_sd_err_05_m = 1.2533 * sd_btot_05_m
btot_med_sd_err_08_m = 1.2533 * sd_btot_08_m
btot_med_sd_err_09_m = 1.2533 * sd_btot_09_m
btot_med_sd_err_10_m = 1.2533 * sd_btot_10_m
btot_med_sd_err_11_m = 1.2533 * sd_btot_11_m
btot_med_sd_err_16_m = 1.2533 * sd_btot_16_m
btot_med_sd_err_17_m = 1.2533 * sd_btot_17_m
btot_med_sd_err_18_m = 1.2533 * sd_btot_18_m
btot_med_sd_err_19_m = 1.2533 * sd_btot_19_m
btot_med_sd_err_20_m = 1.2533 * sd_btot_20_m
btot_med_sd_err_21_m = 1.2533 * sd_btot_21_m
btot_med_sd_err_24_m = 1.2533 * sd_btot_24_m
btot_med_sd_err_25_m = 1.2533 * sd_btot_25_m
btot_med_sd_err_27_m = 1.2533 * sd_btot_27_m
btot_med_sd_err_28_m = 1.2533 * sd_btot_28_m
btot_med_sd_err_34_m = 1.2533 * sd_btot_34_m
btot_med_sd_err_35_m = 1.2533 * sd_btot_35_m

# higher cmls
btot_med_sd_err_01_m = 1.2533 * sd_btot_01_m
btot_med_sd_err_12_m = 1.2533 * sd_btot_12_m
btot_med_sd_err_15_m = 1.2533 * sd_btot_15_m
btot_med_sd_err_26_m = 1.2533 * sd_btot_26_m

#btot_med_sd_errs_m = np.array([btot_med_sd_err_02_m,btot_med_sd_err_03_m,btot_med_sd_err_04_m,btot_med_sd_err_05_m,btot_med_sd_err_08_m,btot_med_sd_err_09_m,btot_med_sd_err_10_m,btot_med_sd_err_11_m,btot_med_sd_err_16_m,btot_med_sd_err_17_m,btot_med_sd_err_18_m,btot_med_sd_err_19_m,btot_med_sd_err_20_m,btot_med_sd_err_21_m,btot_med_sd_err_24_m,btot_med_sd_err_25_m,btot_med_sd_err_27_m,btot_med_sd_err_28_m,btot_med_sd_err_34_m,btot_med_sd_err_35_m])
btot_med_sd_errs_m = np.array([btot_med_sd_err_01_m,btot_med_sd_err_02_m,btot_med_sd_err_03_m,btot_med_sd_err_04_m,btot_med_sd_err_05_m,btot_med_sd_err_08_m,btot_med_sd_err_09_m,btot_med_sd_err_10_m,btot_med_sd_err_11_m,btot_med_sd_err_12_m,btot_med_sd_err_15_m,btot_med_sd_err_16_m,btot_med_sd_err_17_m,btot_med_sd_err_18_m,btot_med_sd_err_19_m,btot_med_sd_err_20_m,btot_med_sd_err_21_m,btot_med_sd_err_24_m,btot_med_sd_err_25_m,btot_med_sd_err_26_m,btot_med_sd_err_27_m,btot_med_sd_err_28_m,btot_med_sd_err_34_m,btot_med_sd_err_35_m])


'''
AVERAGE CLOCK ANGLE
'''

clock_02_mean = np.mean(clock_02)
clock_03_mean = np.mean(clock_03)
clock_04_mean = np.mean(clock_04)
clock_05_mean = np.mean(clock_05)
clock_08_mean = np.mean(clock_08)
clock_09_mean = np.mean(clock_09)
clock_10_mean = np.mean(clock_10)
clock_11_mean = np.mean(clock_11)
clock_16_mean = np.mean(clock_16)
clock_17_mean = np.mean(clock_17)
clock_18_mean = np.mean(clock_18)
clock_19_mean = np.mean(clock_19)
clock_20_mean = np.mean(clock_20)
clock_21_mean = np.mean(clock_21)
clock_24_mean = np.mean(clock_24)
clock_25_mean = np.mean(clock_25)
clock_27_mean = np.mean(clock_27)
clock_28_mean = np.mean(clock_28)
clock_34_mean = np.mean(clock_34)
clock_35_mean = np.mean(clock_35)

# high cmls
clock_01_mean = np.mean(clock_01)
clock_12_mean = np.mean(clock_12)
clock_15_mean = np.mean(clock_15)
clock_26_mean = np.mean(clock_26)

#clock_means = [clock_02_mean,clock_03_mean,clock_04_mean, clock_05_mean, clock_08_mean, clock_09_mean, clock_10_mean,clock_11_mean, clock_16_mean, clock_17_mean, clock_18_mean, clock_19_mean, clock_20_mean, clock_21_mean,clock_24_mean,clock_25_mean,clock_27_mean, clock_28_mean,clock_34_mean,clock_25_mean]
clock_means = [clock_01_mean,clock_02_mean,clock_03_mean,clock_04_mean, clock_05_mean, clock_08_mean, clock_09_mean,clock_10_mean,clock_11_mean, clock_12_mean,clock_15_mean,clock_16_mean, clock_17_mean, clock_18_mean, clock_19_mean, clock_20_mean, clock_21_mean,clock_24_mean,clock_25_mean,clock_26_mean,clock_27_mean, clock_28_mean,clock_34_mean,clock_25_mean]
clock_means = np.array(clock_means)


clock_02_median = np.median(clock_02)
clock_03_median = np.median(clock_03)
clock_04_median = np.median(clock_04)
clock_05_median = np.median(clock_05)
clock_08_median = np.median(clock_08)
clock_09_median = np.median(clock_09)
clock_10_median = np.median(clock_10)
clock_11_median = np.median(clock_11)
clock_16_median = np.median(clock_16)
clock_17_median = np.median(clock_17)
clock_18_median = np.median(clock_18)
clock_19_median = np.median(clock_19)
clock_20_median = np.median(clock_20)
clock_21_median = np.median(clock_21)
clock_24_median = np.median(clock_24)
clock_25_median = np.median(clock_25)
clock_27_median = np.median(clock_27)
clock_28_median = np.median(clock_28)
clock_34_median = np.median(clock_34)
clock_35_median = np.median(clock_35)

# hiigh cmls
clock_01_median = np.median(clock_01)
clock_12_median = np.median(clock_12)
clock_15_median = np.median(clock_15)
clock_26_median = np.median(clock_26)

#clock_medians = [clock_02_median,clock_03_median,clock_04_median, clock_05_median, clock_08_median, clock_09_median, clock_10_median,clock_11_median,clock_16_median, clock_17_median, clock_18_median, clock_19_median, clock_20_median, clock_21_median,clock_24_median,clock_25_median,clock_27_median,clock_28_median,clock_34_median,clock_35_median]
clock_medians = [clock_01_median,clock_02_median,clock_03_median,clock_04_median, clock_05_median, clock_08_median, clock_09_median,clock_10_median,clock_11_median,clock_12_median,clock_15_median,clock_16_median, clock_17_median, clock_18_median, clock_19_median, clock_20_median, clock_21_median,clock_24_median,clock_25_median,clock_26_median,clock_27_median,clock_28_median,clock_34_median,clock_35_median]
clock_medians = np.array(clock_medians)


'''
make sure to check which clock angles will need converting - any in which the median/mean causes large/small angles
'''

def clock_format_convert(clock_angle):
    new_clock_angle = []
    for i in range(len(clock_angle)):
        if clock_angle[i] < 0:
            clock = clock_angle[i] + 360
            new_clock_angle.append(clock)
        else:
            new_clock_angle.append(clock_angle[i])
    return new_clock_angle


new_clock_24 = clock_format_convert(clock_24)

new_clock_24_median = np.median(new_clock_24)
# # #new_clock_24_median = new_clock_24_median_1 - 360
clock_medians[17] = new_clock_24_median

new_clock_24_mean = np.mean(new_clock_24)
# # new_clock_24_mean = new_clock_24_mean_1 - 360
clock_means[17] = new_clock_24_mean

sd_clock_02 = standard_deviation(clock_02, clock_02_mean)
sd_clock_03 = standard_deviation(clock_03, clock_03_mean)
sd_clock_04 = standard_deviation(clock_04, clock_04_mean)
sd_clock_05 = standard_deviation(clock_05, clock_05_mean)
sd_clock_08 = standard_deviation(clock_08, clock_08_mean)
sd_clock_09 = standard_deviation(clock_09, clock_09_mean)
sd_clock_10 = standard_deviation(clock_10, clock_10_mean)
sd_clock_11 = standard_deviation(clock_11, clock_11_mean)
sd_clock_16 = standard_deviation(clock_16, clock_16_mean)
sd_clock_17 = standard_deviation(clock_17, clock_17_mean)
sd_clock_18 = standard_deviation(clock_18, clock_18_mean)
sd_clock_19 = standard_deviation(clock_19, clock_19_mean)
sd_clock_20 = standard_deviation(clock_20, clock_20_mean)
sd_clock_21 = standard_deviation(clock_21, clock_21_mean)
sd_clock_24 = standard_deviation(new_clock_24, new_clock_24_mean)
sd_clock_25 = standard_deviation(clock_25, clock_25_mean)
sd_clock_27 = standard_deviation(clock_27, clock_27_mean)
sd_clock_28 = standard_deviation(clock_28, clock_28_mean)
sd_clock_34 = standard_deviation(clock_34, clock_34_mean)
sd_clock_35 = standard_deviation(clock_35, clock_35_mean)

# high cmls
sd_clock_01 = standard_deviation(clock_01, clock_01_mean)
sd_clock_12 = standard_deviation(clock_12, clock_12_mean)
sd_clock_15 = standard_deviation(clock_15, clock_15_mean)
sd_clock_26 = standard_deviation(clock_26, clock_26_mean)

#clock_sds = [sd_clock_02,sd_clock_03,sd_clock_04, sd_clock_05, sd_clock_08, sd_clock_09, sd_clock_10, sd_clock_11,sd_clock_16, sd_clock_17, sd_clock_18, sd_clock_19, sd_clock_20, sd_clock_21,sd_clock_24,sd_clock_25,sd_clock_27,sd_clock_28,sd_clock_34,sd_clock_35]
clock_sds = [sd_clock_01,sd_clock_02,sd_clock_03,sd_clock_04, sd_clock_05, sd_clock_08, sd_clock_09, sd_clock_10, sd_clock_11,sd_clock_12,sd_clock_15,sd_clock_16, sd_clock_17, sd_clock_18, sd_clock_19, sd_clock_20, sd_clock_21,sd_clock_24,sd_clock_25,sd_clock_26,sd_clock_27,sd_clock_28,sd_clock_34,sd_clock_35]
clock_sds = np.array(clock_sds)


med_clock_02_err = 1.2533 * sd_clock_02
med_clock_03_err = 1.2533 * sd_clock_03
med_clock_04_err = 1.2533 * sd_clock_04
med_clock_05_err = 1.2533 * sd_clock_05
med_clock_08_err = 1.2533 * sd_clock_08
med_clock_09_err = 1.2533 * sd_clock_09
med_clock_10_err = 1.2533 * sd_clock_10
med_clock_11_err = 1.2533 * sd_clock_11
med_clock_16_err = 1.2533 * sd_clock_16
med_clock_17_err = 1.2533 * sd_clock_17
med_clock_18_err = 1.2533 * sd_clock_18
med_clock_19_err = 1.2533 * sd_clock_19
med_clock_20_err = 1.2533 * sd_clock_20
med_clock_21_err = 1.2533 * sd_clock_21
med_clock_24_err = 1.2533 * sd_clock_24
med_clock_25_err = 1.2533 * sd_clock_25
med_clock_27_err = 1.2533 * sd_clock_27
med_clock_28_err = 1.2533 * sd_clock_28
med_clock_34_err = 1.2533 * sd_clock_34
med_clock_35_err = 1.2533 * sd_clock_35

# high cmls
med_clock_01_err = 1.2533 * sd_clock_01
med_clock_12_err = 1.2533 * sd_clock_12
med_clock_15_err = 1.2533 * sd_clock_15
med_clock_26_err = 1.2533 * sd_clock_26

#clock_med_errs = [med_clock_02_err,med_clock_03_err,med_clock_04_err, med_clock_05_err, med_clock_08_err, med_clock_09_err, med_clock_10_err, med_clock_11_err,med_clock_16_err, med_clock_16_err, med_clock_18_err, med_clock_19_err, med_clock_20_err, med_clock_21_err,med_clock_24_err,med_clock_25_err,med_clock_27_err,med_clock_28_err,med_clock_34_err,med_clock_35_err]
clock_med_errs = [med_clock_01_err,med_clock_02_err,med_clock_03_err,med_clock_04_err, med_clock_05_err, med_clock_08_err, med_clock_09_err,med_clock_10_err, med_clock_11_err,med_clock_12_err,med_clock_15_err,med_clock_16_err, med_clock_16_err, med_clock_18_err, med_clock_19_err, med_clock_20_err, med_clock_21_err,med_clock_24_err,med_clock_25_err,med_clock_26_err,med_clock_27_err,med_clock_28_err,med_clock_34_err,med_clock_35_err]
clock_med_errs = np.array(clock_med_errs)



'''
travel time +/- 10% clocks - we're up to here! replace 08,09 and 17 - they're not new clock anymore but
24 and 34 will be (check the values but non error ones need to be shifted)
'''


clock_02_mean_p = np.mean(clock_02_plus)
clock_03_mean_p = np.mean(clock_03_plus)
clock_04_mean_p = np.mean(clock_04_plus)
clock_05_mean_p = np.mean(clock_05_plus)
clock_08_mean_p = np.mean(clock_08_plus)
clock_09_mean_p = np.mean(clock_09_plus)
clock_10_mean_p = np.mean(clock_10_plus)
clock_11_mean_p = np.mean(clock_11_plus)
clock_16_mean_p = np.mean(clock_16_plus)
clock_17_mean_p = np.mean(clock_17_plus)
clock_18_mean_p = np.mean(clock_18_plus)
clock_19_mean_p = np.mean(clock_19_plus)
clock_20_mean_p = np.mean(clock_20_plus)
clock_21_mean_p = np.mean(clock_21_plus)
clock_24_mean_p = np.mean(clock_24_plus)
clock_25_mean_p = np.mean(clock_25_plus)
clock_27_mean_p = np.mean(clock_27_plus)
clock_28_mean_p = np.mean(clock_28_plus)
clock_34_mean_p = np.mean(clock_24_plus)
clock_35_mean_p = np.mean(clock_35_plus)

# high cmls
clock_01_mean_p = np.mean(clock_01_plus)
clock_12_mean_p = np.mean(clock_12_plus)
clock_15_mean_p = np.mean(clock_15_plus)
clock_26_mean_p = np.mean(clock_26_plus)

#clock_means_p = [clock_02_mean_p,clock_03_mean_p,clock_04_mean_p,clock_05_mean_p,clock_08_mean_p,clock_09_mean_p,clock_10_mean_p,clock_11_mean_p,clock_16_mean_p,clock_17_mean_p,clock_18_mean_p,clock_19_mean_p,clock_20_mean_p,clock_21_mean_p,clock_24_mean_p,clock_25_mean_p,clock_27_mean_p,clock_28_mean_p,clock_34_mean_p,clock_35_mean_p]
clock_means_p = [clock_01_mean_p,clock_02_mean_p,clock_03_mean_p,clock_04_mean_p,clock_05_mean_p,clock_08_mean_p,clock_09_mean_p,clock_10_mean_p,clock_11_mean_p,clock_12_mean_p,clock_15_mean_p,clock_16_mean_p,clock_17_mean_p,clock_18_mean_p,clock_19_mean_p,clock_20_mean_p,clock_21_mean_p,clock_24_mean_p,clock_25_mean_p,clock_26_mean_p,clock_27_mean_p,clock_28_mean_p,clock_34_mean_p,clock_35_mean_p]
clock_means_p = np.array(clock_means_p)


new_clock_24_plus = clock_format_convert(clock_24_plus)
new_clock_18_plus = clock_format_convert(clock_18_plus)

new_clock_24_mean_p = np.mean(new_clock_24_plus)
new_clock_18_mean_p = np.mean(new_clock_18_plus)
# # new_clock_24_mean = new_clock_24_mean_1 - 360
clock_means_p[17] = new_clock_24_mean_p
clock_means_p[13] = new_clock_18_mean_p


sd_clock_02_p = standard_deviation(clock_02_plus, clock_02_mean_p)
sd_clock_03_p = standard_deviation(clock_03_plus, clock_03_mean_p)
sd_clock_04_p = standard_deviation(clock_04_plus, clock_04_mean_p)
sd_clock_05_p = standard_deviation(clock_05_plus, clock_05_mean_p)
sd_clock_08_p = standard_deviation(clock_08_plus, clock_08_mean_p)
sd_clock_09_p = standard_deviation(clock_09_plus, clock_09_mean_p)
sd_clock_10_p = standard_deviation(clock_10_plus, clock_10_mean_p)
sd_clock_11_p = standard_deviation(clock_11_plus, clock_11_mean_p)
sd_clock_16_p = standard_deviation(clock_16_plus, clock_16_mean_p)
sd_clock_17_p = standard_deviation(clock_17_plus, clock_17_mean_p)
sd_clock_18_p = standard_deviation(new_clock_18_plus, new_clock_18_mean_p)
sd_clock_19_p = standard_deviation(clock_19_plus, clock_19_mean_p)
sd_clock_20_p = standard_deviation(clock_20_plus, clock_20_mean_p)
sd_clock_21_p = standard_deviation(clock_21_plus, clock_21_mean_p)
sd_clock_24_p = standard_deviation(new_clock_24_plus, new_clock_24_mean_p)
sd_clock_25_p = standard_deviation(clock_25_plus, clock_25_mean_p)
sd_clock_27_p = standard_deviation(clock_27_plus, clock_27_mean_p)
sd_clock_28_p = standard_deviation(clock_28_plus, clock_28_mean_p)
sd_clock_34_p = standard_deviation(clock_34_plus, clock_34_mean_p)
sd_clock_35_p = standard_deviation(clock_35_plus, clock_35_mean_p)

# high cmls
sd_clock_01_p = standard_deviation(clock_01_plus, clock_01_mean_p)
sd_clock_12_p = standard_deviation(clock_12_plus, clock_12_mean_p)
sd_clock_15_p = standard_deviation(clock_15_plus, clock_15_mean_p)
sd_clock_26_p = standard_deviation(clock_26_plus, clock_26_mean_p)


clock_02_median_p = np.median(clock_02_plus)
clock_03_median_p = np.median(clock_03_plus)
clock_04_median_p = np.median(clock_04_plus)
clock_05_median_p = np.median(clock_05_plus)
clock_08_median_p = np.median(clock_08_plus)
clock_09_median_p = np.median(clock_09_plus)
clock_10_median_p = np.median(clock_10_plus)
clock_11_median_p = np.median(clock_11_plus)
clock_16_median_p = np.median(clock_16_plus)
clock_17_median_p = np.median(clock_17_plus)
clock_18_median_p = np.median(clock_18_plus)
clock_19_median_p = np.median(clock_19_plus)
clock_20_median_p = np.median(clock_20_plus)
clock_21_median_p = np.median(clock_21_plus)
clock_24_median_p = np.median(clock_24_plus)
clock_25_median_p = np.median(clock_25_plus)
clock_27_median_p = np.median(clock_27_plus)
clock_28_median_p = np.median(clock_28_plus)
clock_34_median_p = np.median(clock_34_plus)
clock_35_median_p = np.median(clock_35_plus)

# high cmls
clock_01_median_p = np.median(clock_01_plus)
clock_12_median_p = np.median(clock_12_plus)
clock_15_median_p = np.median(clock_15_plus)
clock_26_median_p = np.median(clock_26_plus)

#clock_medians_p = [clock_02_median_p,clock_03_median_p,clock_04_median_p,clock_05_median_p,clock_08_median_p,clock_09_median_p,clock_10_median_p,clock_11_median_p,clock_16_median_p,clock_17_median_p,clock_18_median_p,clock_19_median_p,clock_20_median_p,clock_21_median_p,clock_24_median_p,clock_25_median_p,clock_27_median_p,clock_28_median_p,clock_34_median_p,clock_35_median_p]
clock_medians_p = [clock_01_median_p,clock_02_median_p,clock_03_median_p,clock_04_median_p,clock_05_median_p,clock_08_median_p,clock_09_median_p,clock_10_median_p,clock_11_median_p,clock_12_median_p,clock_15_median_p,clock_16_median_p,clock_17_median_p,clock_18_median_p,clock_19_median_p,clock_20_median_p,clock_21_median_p,clock_24_median_p,clock_25_median_p,clock_26_median_p,clock_27_median_p,clock_28_median_p,clock_34_median_p,clock_35_median_p]
clock_medians_p = np.array(clock_medians_p)

new_clock_24_median_p = np.median(new_clock_24_plus)
new_clock_18_median_p = np.median(new_clock_18_plus)
# # #new_clock_24_median = new_clock_24_median_1 - 360
clock_medians_p[17] = new_clock_24_median_p
clock_medians_p[13] = new_clock_18_median_p


med_clock_02_err_p = 1.2533 * sd_clock_02_p
med_clock_03_err_p = 1.2533 * sd_clock_03_p
med_clock_04_err_p = 1.2533 * sd_clock_04_p
med_clock_05_err_p = 1.2533 * sd_clock_05_p
med_clock_08_err_p = 1.2533 * sd_clock_08_p
med_clock_09_err_p = 1.2533 * sd_clock_09_p
med_clock_10_err_p = 1.2533 * sd_clock_10_p
med_clock_11_err_p = 1.2533 * sd_clock_11_p
med_clock_16_err_p = 1.2533 * sd_clock_16_p
med_clock_17_err_p = 1.2533 * sd_clock_17_p
med_clock_18_err_p = 1.2533 * sd_clock_18_p
med_clock_19_err_p = 1.2533 * sd_clock_19_p
med_clock_20_err_p = 1.2533 * sd_clock_20_p
med_clock_21_err_p = 1.2533 * sd_clock_21_p
med_clock_24_err_p = 1.2533 * sd_clock_24_p
med_clock_25_err_p = 1.2533 * sd_clock_25_p
med_clock_27_err_p = 1.2533 * sd_clock_27_p
med_clock_28_err_p = 1.2533 * sd_clock_28_p
med_clock_34_err_p = 1.2533 * sd_clock_34_p
med_clock_35_err_p = 1.2533 * sd_clock_35_p

# high cmls
med_clock_01_err_p = 1.2533 * sd_clock_11_p
med_clock_12_err_p = 1.2533 * sd_clock_12_p
med_clock_15_err_p = 1.2533 * sd_clock_15_p
med_clock_26_err_p = 1.2533 * sd_clock_26_p


#clock_med_errs_p = [med_clock_02_err_p,med_clock_03_err_p,med_clock_04_err_p, med_clock_05_err_p, med_clock_08_err_p, med_clock_09_err_p, med_clock_10_err_p,med_clock_11_err_p, med_clock_16_err_p, med_clock_16_err_p, med_clock_18_err_p, med_clock_19_err_p, med_clock_20_err_p, med_clock_21_err_p,med_clock_24_err_p,med_clock_25_err_p,med_clock_27_err_p,med_clock_28_err_p,med_clock_34_err_p,med_clock_35_err_p]
clock_med_errs_p = [med_clock_01_err_p,med_clock_02_err_p,med_clock_03_err_p,med_clock_04_err_p, med_clock_05_err_p, med_clock_08_err_p, med_clock_09_err_p,med_clock_10_err_p,med_clock_11_err_p, med_clock_12_err_p,med_clock_15_err_p,med_clock_16_err_p, med_clock_16_err_p, med_clock_18_err_p, med_clock_19_err_p, med_clock_20_err_p, med_clock_21_err_p,med_clock_24_err_p,med_clock_25_err_p,med_clock_26_err_p,med_clock_27_err_p,med_clock_28_err_p,med_clock_34_err_p,med_clock_35_err_p]
clock_med_errs_p = np.array(clock_med_errs_p)


new_clock_24_minus = clock_format_convert(clock_24_minus)
new_clock_18_minus = clock_format_convert(clock_18_minus)

clock_02_mean_m = np.mean(clock_02_minus)
clock_03_mean_m = np.mean(clock_03_minus)
clock_04_mean_m = np.mean(clock_04_minus)
clock_05_mean_m = np.mean(clock_05_minus)
clock_08_mean_m = np.mean(clock_08_minus)
clock_09_mean_m = np.mean(clock_09_minus)
clock_10_mean_m = np.mean(clock_10_minus)
clock_11_mean_m = np.mean(clock_11_minus)
clock_16_mean_m = np.mean(clock_16_minus)
clock_17_mean_m = np.mean(clock_17_minus)
clock_18_mean_m = np.mean(new_clock_18_minus)
clock_19_mean_m = np.mean(clock_19_minus)
clock_20_mean_m = np.mean(clock_20_minus)
clock_21_mean_m = np.mean(clock_21_minus)
clock_24_mean_m = np.mean(new_clock_24_minus)
clock_25_mean_m = np.mean(clock_25_minus)
clock_27_mean_m = np.mean(clock_27_minus)
clock_28_mean_m = np.mean(clock_28_minus)
clock_34_mean_m = np.mean(clock_34_minus)
clock_35_mean_m = np.mean(clock_35_minus)

# high cmls
clock_01_mean_m = np.mean(clock_01_minus)
clock_12_mean_m = np.mean(clock_12_minus)
clock_15_mean_m = np.mean(clock_15_minus)
clock_26_mean_m = np.mean(clock_26_minus)

#clock_means_m = [clock_02_mean_m,clock_03_mean_m,clock_04_mean_m,clock_05_mean_m,clock_08_mean_m,clock_09_mean_m,clock_10_mean_m,clock_11_mean_m,clock_16_mean_m,clock_17_mean_m,clock_18_mean_m,clock_19_mean_m,clock_20_mean_m,clock_21_mean_m,clock_24_mean_m,clock_25_mean_m,clock_27_mean_m,clock_28_mean_m,clock_34_mean_m,clock_35_mean_m]
clock_means_m = [clock_01_mean_m,clock_02_mean_m,clock_03_mean_m,clock_04_mean_m,clock_05_mean_m,clock_08_mean_m,clock_09_mean_m,clock_10_mean_m,clock_11_mean_m,clock_12_mean_m,clock_15_mean_m,clock_16_mean_m,clock_17_mean_m,clock_18_mean_m,clock_19_mean_m,clock_20_mean_m,clock_21_mean_m,clock_24_mean_m,clock_25_mean_m,clock_26_mean_m,clock_27_mean_m,clock_28_mean_m,clock_34_mean_m,clock_35_mean_m]
clock_means_m = np.array(clock_means_m)


sd_clock_02_m = standard_deviation(clock_02_minus, clock_02_mean_m)
sd_clock_03_m = standard_deviation(clock_03_minus, clock_03_mean_m)
sd_clock_04_m = standard_deviation(clock_04_minus, clock_04_mean_m)
sd_clock_05_m = standard_deviation(clock_05_minus, clock_05_mean_m)
sd_clock_08_m = standard_deviation(clock_08_minus, clock_08_mean_m)
sd_clock_09_m = standard_deviation(clock_09_minus, clock_09_mean_m)
sd_clock_10_m = standard_deviation(clock_10_minus, clock_10_mean_m)
sd_clock_11_m = standard_deviation(clock_11_minus, clock_11_mean_m)
sd_clock_16_m = standard_deviation(clock_16_minus, clock_16_mean_m)
sd_clock_17_m = standard_deviation(clock_17_minus, clock_17_mean_m)
sd_clock_18_m = standard_deviation(new_clock_18_minus, clock_18_mean_m)
sd_clock_19_m = standard_deviation(clock_19_minus, clock_19_mean_m)
sd_clock_20_m = standard_deviation(clock_20_minus, clock_20_mean_m)
sd_clock_21_m = standard_deviation(clock_21_minus, clock_21_mean_m)
sd_clock_24_m = standard_deviation(new_clock_24_minus, clock_24_mean_m)
sd_clock_25_m = standard_deviation(clock_25_minus, clock_25_mean_m)
sd_clock_27_m = standard_deviation(clock_27_minus, clock_27_mean_m)
sd_clock_28_m = standard_deviation(clock_28_minus, clock_28_mean_m)
sd_clock_34_m = standard_deviation(clock_34_minus, clock_34_mean_m)
sd_clock_35_m = standard_deviation(clock_35_minus, clock_35_mean_m)

# high cmls
sd_clock_01_m = standard_deviation(clock_01_minus, clock_01_mean_m)
sd_clock_12_m = standard_deviation(clock_12_minus, clock_12_mean_m)
sd_clock_15_m = standard_deviation(clock_15_minus, clock_15_mean_m)
sd_clock_26_m = standard_deviation(clock_26_minus, clock_26_mean_m)


clock_02_median_m = np.median(clock_02_minus)
clock_03_median_m = np.median(clock_03_minus)
clock_04_median_m = np.median(clock_04_minus)
clock_05_median_m = np.median(clock_05_minus)
clock_08_median_m = np.median(clock_08_minus)
clock_09_median_m = np.median(clock_09_minus)
clock_10_median_m = np.median(clock_10_minus)
clock_11_median_m = np.median(clock_11_minus)
clock_16_median_m = np.median(clock_16_minus)
clock_17_median_m = np.median(clock_17_minus)
clock_18_median_m = np.median(new_clock_18_minus)
clock_19_median_m = np.median(clock_19_minus)
clock_20_median_m = np.median(clock_20_minus)
clock_21_median_m = np.median(clock_21_minus)
clock_24_median_m = np.median(new_clock_24_minus)
clock_25_median_m = np.median(clock_25_minus)
clock_27_median_m = np.median(clock_27_minus)
clock_28_median_m = np.median(clock_28_minus)
clock_34_median_m = np.median(clock_34_minus)
clock_35_median_m = np.median(clock_35_minus)

# high cmls
clock_01_median_m = np.median(clock_01_minus)
clock_12_median_m = np.median(clock_12_minus)
clock_15_median_m = np.median(clock_15_minus)
clock_26_median_m = np.median(clock_26_minus)

#clock_medians_m = [clock_02_median_m,clock_03_median_m,clock_04_median_m,clock_05_median_m,clock_08_median_m,clock_09_median_m,clock_10_median_m,clock_11_median_m,clock_16_median_m,clock_17_median_m,clock_18_median_m,clock_19_median_m,clock_20_median_m,clock_21_median_m,clock_24_median_m,clock_25_median_m,clock_27_median_m,clock_28_median_m,clock_34_median_m,clock_35_median_m]
clock_medians_m = [clock_01_median_m,clock_02_median_m,clock_03_median_m,clock_04_median_m,clock_05_median_m,clock_08_median_m,clock_09_median_m,clock_10_median_m,clock_11_median_m,clock_12_median_m,clock_15_median_m,clock_16_median_m,clock_17_median_m,clock_18_median_m,clock_19_median_m,clock_20_median_m,clock_21_median_m,clock_24_median_m,clock_25_median_m,clock_26_median_m,clock_27_median_m,clock_28_median_m,clock_34_median_m,clock_35_median_m]
clock_medians_m = np.array(clock_medians_m)

#clock_medians_m[17] = clock_24_median_m - 360

med_clock_02_err_m = 1.2533 * sd_clock_02_m
med_clock_03_err_m = 1.2533 * sd_clock_03_m
med_clock_04_err_m = 1.2533 * sd_clock_04_m
med_clock_05_err_m = 1.2533 * sd_clock_05_m
med_clock_08_err_m = 1.2533 * sd_clock_08_m
med_clock_09_err_m = 1.2533 * sd_clock_09_m
med_clock_10_err_m = 1.2533 * sd_clock_10_m
med_clock_11_err_m = 1.2533 * sd_clock_11_m
med_clock_16_err_m = 1.2533 * sd_clock_16_m
med_clock_17_err_m = 1.2533 * sd_clock_17_m
med_clock_18_err_m = 1.2533 * sd_clock_18_m
med_clock_19_err_m = 1.2533 * sd_clock_19_m
med_clock_20_err_m = 1.2533 * sd_clock_20_m
med_clock_21_err_m = 1.2533 * sd_clock_21_m
med_clock_24_err_m = 1.2533 * sd_clock_24_m
med_clock_25_err_m = 1.2533 * sd_clock_25_m
med_clock_27_err_m = 1.2533 * sd_clock_27_m
med_clock_28_err_m = 1.2533 * sd_clock_28_m
med_clock_34_err_m = 1.2533 * sd_clock_34_m
med_clock_35_err_m = 1.2533 * sd_clock_35_m

# high cmls
med_clock_01_err_m = 1.2533 * sd_clock_01_m
med_clock_12_err_m = 1.2533 * sd_clock_12_m
med_clock_15_err_m = 1.2533 * sd_clock_15_m
med_clock_26_err_m = 1.2533 * sd_clock_26_m

#clock_med_errs_m = [med_clock_02_err_m,med_clock_03_err_m,med_clock_04_err_m, med_clock_05_err_m, med_clock_08_err_m, med_clock_09_err_m, med_clock_10_err_m,med_clock_11_err_m, med_clock_16_err_m, med_clock_16_err_m, med_clock_18_err_m, med_clock_19_err_m, med_clock_20_err_m, med_clock_21_err_m,med_clock_24_err_m,med_clock_25_err_m,med_clock_27_err_m,med_clock_28_err_m,med_clock_34_err_m,med_clock_35_err_m]
clock_med_errs_m = [med_clock_01_err_m,med_clock_02_err_m,med_clock_03_err_m,med_clock_04_err_m, med_clock_05_err_m, med_clock_08_err_m, med_clock_09_err_m,med_clock_10_err_m,med_clock_11_err_m, med_clock_12_err_m,med_clock_15_err_m,med_clock_16_err_m, med_clock_16_err_m, med_clock_18_err_m, med_clock_19_err_m, med_clock_20_err_m, med_clock_21_err_m,med_clock_24_err_m,med_clock_26_err_m,med_clock_25_err_m,med_clock_27_err_m,med_clock_28_err_m,med_clock_34_err_m,med_clock_35_err_m]
clock_med_errs_m = np.array(clock_med_errs_m)




'''
colorcoding by quadrant
'''

Bz_pos = []
Bz_pos_power = []
Bz_pos_swirlpower = []
Bz_pos_duskpower = []
Bz_pos_noonpower = []
Bz_pos_pressure = []
Bz_pos_b_perp = []
Bz_pos_LL = []
Bz_pos_HL_pos = []
Bz_pos_HL_neg = []
Bz_pos_gersh = []
Bz_pos_KH = []
Bz_pos_btot = []
Bz_pos_KH_dawn = []
Bz_pos_KH_dusk = []

Bz_neg = []
Bz_neg_power = []
Bz_neg_swirlpower = []
Bz_neg_duskpower = []
Bz_neg_noonpower = []
Bz_neg_pressure = []
Bz_neg_b_perp = []
Bz_neg_LL = []
Bz_neg_HL_pos = []
Bz_neg_HL_neg = []
Bz_neg_gersh = []
Bz_neg_KH = []
Bz_neg_btot = []
Bz_neg_KH_dawn = []
Bz_neg_KH_dusk = []

By_pos = []
By_pos_power = []
By_pos_swirlpower = []
By_pos_duskpower = []
By_pos_noonpower = []
By_pos_pressure = []
By_pos_b_perp = []
By_pos_LL = []
By_pos_HL_pos = []
By_pos_HL_neg = []
By_pos_gersh = []
By_pos_KH = []
By_pos_btot = []
By_pos_KH_dawn = []
By_pos_KH_dusk = []

By_neg = []
By_neg_power = []
By_neg_swirlpower = []
By_neg_duskpower = []
By_neg_noonpower = []
By_neg_pressure = []
By_neg_b_perp = []
By_neg_LL = []
By_neg_HL_pos = []
By_neg_HL_neg = []
By_neg_gersh = []
By_neg_KH = []
By_neg_btot = []
By_neg_KH_dawn = []
By_neg_KH_dusk = []


for angle in range(len(clock_medians)):
    if clock_medians[angle] > -45 and clock_medians[angle] < 45:
        Bz_pos.append(clock_medians[angle])
        Bz_pos_power.append(polar_medians[angle])
        Bz_pos_pressure.append(pressure_medians[angle])
        Bz_pos_b_perp.append(b_perp_medians[angle])
        Bz_pos_btot.append(btot_medians[angle])
        
        Bz_pos_LL.append(LL_medians[angle])
        Bz_pos_HL_neg.append(HL_medians_neg[angle])
        Bz_pos_HL_pos.append(HL_medians_pos[angle])
        Bz_pos_gersh.append(gersh_medians[angle])
        Bz_pos_KH.append(KH_medians[angle])
        
        Bz_pos_swirlpower.append(medians_swirl[angle])
        Bz_pos_duskpower.append(dusk_medians[angle])
        Bz_pos_noonpower.append(noon_medians[angle])
        
        Bz_pos_KH_dawn.append(KH_dawn_medians[angle])
        Bz_pos_KH_dusk.append(KH_dusk_medians[angle])
        
    elif clock_medians[angle] > 45 and clock_medians[angle] < 135:
        By_pos.append(clock_medians[angle])
        By_pos_power.append(polar_medians[angle])
        
        By_pos_swirlpower.append(medians_swirl[angle])
        By_pos_duskpower.append(dusk_medians[angle])
        By_pos_noonpower.append(noon_medians[angle])
        
        By_pos_pressure.append(pressure_medians[angle])
        By_pos_LL.append(LL_medians[angle])
        By_pos_HL_pos.append(HL_medians_pos[angle])
        By_pos_HL_neg.append(HL_medians_neg[angle])
        By_pos_gersh.append(gersh_medians[angle])
        By_pos_KH.append(KH_medians[angle])
        By_pos_b_perp.append(b_perp_medians[angle])
        By_pos_btot.append(btot_medians[angle])
        
        By_pos_KH_dawn.append(KH_dawn_medians[angle])
        By_pos_KH_dusk.append(KH_dusk_medians[angle])
        
    elif clock_medians[angle] < -45 and clock_medians[angle] > -135:
        By_neg.append(clock_medians[angle])
        By_neg_power.append(polar_medians[angle])
        
        By_neg_swirlpower.append(medians_swirl[angle])
        By_neg_duskpower.append(dusk_medians[angle])
        By_neg_noonpower.append(noon_medians[angle])
        
        By_neg_pressure.append(pressure_medians[angle])
        By_neg_LL.append(LL_medians[angle])
        By_neg_HL_pos.append(HL_medians_pos[angle])
        By_neg_HL_neg.append(HL_medians_neg[angle])
        By_neg_gersh.append(gersh_medians[angle])
        By_neg_KH.append(KH_medians[angle])
        By_neg_b_perp.append(b_perp_medians[angle])
        By_neg_btot.append(btot_medians[angle])
        
        By_neg_KH_dawn.append(KH_dawn_medians[angle])
        By_neg_KH_dusk.append(KH_dusk_medians[angle])
        
        
    else:
        Bz_neg.append(clock_medians[angle])
        Bz_neg_power.append(polar_medians[angle])
        
        Bz_neg_swirlpower.append(medians_swirl[angle])
        Bz_neg_duskpower.append(dusk_medians[angle])
        Bz_neg_noonpower.append(noon_medians[angle])
        
        Bz_neg_pressure.append(pressure_medians[angle])
        Bz_neg_LL.append(LL_medians[angle])
        Bz_neg_HL_pos.append(HL_medians_pos[angle])
        Bz_neg_HL_neg.append(HL_medians_neg[angle])
        Bz_neg_gersh.append(gersh_medians[angle])
        Bz_neg_KH.append(KH_medians[angle])
        Bz_neg_b_perp.append(b_perp_medians[angle])
        Bz_neg_btot.append(btot_medians[angle])
        
        Bz_neg_KH_dawn.append(KH_dawn_medians[angle])
        Bz_neg_KH_dusk.append(KH_dusk_medians[angle])
        
                    
Bz_pos = np.array(Bz_pos)
Bz_pos_power = np.array(Bz_pos_power)
By_pos = np.array(By_pos)
By_pos_power = np.array(By_pos_power)
Bz_neg = np.array(Bz_neg)
Bz_neg_power = np.array(Bz_neg_power)
By_neg = np.array(By_neg)
By_neg_power = np.array(By_neg_power)

Bz_pos_pressure = np.array(Bz_pos_pressure)
Bz_neg_pressure = np.array(Bz_neg_pressure)
By_pos_pressure = np.array(By_pos_pressure)
By_neg_pressure = np.array(By_neg_pressure)

Bz_pos_HL_pos = np.array(Bz_pos_HL_pos)
Bz_neg_HL_pos = np.array(Bz_neg_HL_pos)
By_pos_HL_pos = np.array(By_pos_HL_pos)
By_neg_HL_pos = np.array(By_neg_HL_pos)

Bz_pos_HL_neg = np.array(Bz_pos_HL_neg)
Bz_neg_HL_neg = np.array(Bz_neg_HL_neg)
By_pos_HL_neg = np.array(By_pos_HL_neg)
By_neg_HL_neg = np.array(By_neg_HL_neg)

Bz_pos_LL = np.array(Bz_pos_LL)
Bz_neg_LL = np.array(Bz_neg_LL)
By_pos_LL = np.array(By_pos_LL)
By_neg_LL= np.array(By_neg_LL)

Bz_pos_b_perp = np.array(Bz_pos_b_perp)
Bz_neg_b_perp = np.array(Bz_neg_b_perp)
By_pos_b_perp = np.array(By_pos_b_perp)
By_neg_b_perp = np.array(By_neg_b_perp)

Bz_pos_btot = np.array(Bz_pos_btot)
Bz_neg_btot = np.array(Bz_neg_btot)
By_pos_btot = np.array(By_pos_btot)
By_neg_btot = np.array(By_neg_btot)

Bz_pos_gersh = np.array(Bz_pos_gersh)
Bz_neg_gersh = np.array(Bz_neg_gersh)
By_pos_gersh = np.array(By_pos_gersh)
By_neg_gersh = np.array(By_neg_gersh)

Bz_pos_KH = np.array(Bz_pos_KH)
Bz_neg_KH = np.array(Bz_neg_KH)
By_pos_KH = np.array(By_pos_KH)
By_neg_KH = np.array(By_neg_KH)

Bz_pos_KH_dawn = np.array(Bz_pos_KH_dawn)
Bz_neg_KH_dawn = np.array(Bz_neg_KH_dawn)
By_pos_KH_dawn = np.array(By_pos_KH_dawn)
By_neg_KH_dawn = np.array(By_neg_KH_dawn)

Bz_pos_KH_dusk = np.array(Bz_pos_KH_dusk)
Bz_neg_KH_dusk = np.array(Bz_neg_KH_dusk)
By_pos_KH_dusk = np.array(By_pos_KH_dusk)
By_neg_KH_dusk = np.array(By_neg_KH_dusk)
        
Bz_pos_swirlpower = np.array(Bz_pos_swirlpower)
Bz_pos_duskpower = np.array(Bz_pos_duskpower)
Bz_pos_noonpower = np.array(Bz_pos_noonpower)
Bz_neg_swirlpower = np.array(Bz_neg_swirlpower)
Bz_neg_duskpower = np.array(Bz_neg_duskpower)
Bz_neg_noonpower = np.array(Bz_neg_noonpower)
By_pos_swirlpower = np.array(By_pos_swirlpower)
By_pos_duskpower = np.array(By_pos_duskpower)
By_pos_noonpower = np.array(By_pos_noonpower)
By_neg_swirlpower = np.array(By_neg_swirlpower)
By_neg_duskpower = np.array(By_neg_duskpower)
By_neg_noonpower = np.array(By_neg_noonpower)


if plotting == 'b_tot':
    
    x_new = np.linspace(min(btot_medians),max(btot_medians),len(btot_medians))
    
    a, acov = np.polyfit(btot_medians, medians_swirl, 1, cov=True)
    aa = np.poly1d(a)
    # print("Swirl Fit Gradient:")
    # print(a[0])
    # print("Error of Swirl Fit:")
    # print(np.sqrt(np.diag(acov)))
    y_regress_swirl = aa(btot_medians)

    print(f"Testing '{plotting}'")
    print('Swirl Region Fit')
    ss_res = np.sum((medians_swirl - y_regress_swirl)**2)
    ss_tot = np.sum((medians_swirl - np.mean(medians_swirl))**2)
    r_squared_s = 1 - (ss_res / ss_tot)
    print(f"R squared: {r_squared_s:.4f}")
    
    adj_s = 1-r_squared_s
    adjr_s = 1 - ((adj_s*23)/22)
    print(f"Adjusted R squared: {adjr_s:.4f}")
    
    swirl_grad = round(a[0],2)
    swirl_grad_err = np.sqrt(np.diag(acov))
    swirl_grad_err = round(swirl_grad_err[0], 2)
    
    swirl_perp_fit = (f'{swirl_grad} ± {swirl_grad_err}')
    
    
    a_p,acov_p = np.polyfit(btot_medians_p, medians_swirl,1, cov=True)
    aa_p = np.poly1d(a_p)
    # print(f"Swirl Fit +{error} Gradient:")
    # print(a_p[0])
    # print(f"Error of Swirl Fit +{error}:")
    # print(np.sqrt(np.diag(acov_p)))
    y_regress_swirl_p = aa_p(btot_medians_p)

    print("---------------")
    print('Swirl Region +20% Fit')
    ssp_res = np.sum((medians_swirl - y_regress_swirl_p)**2)
    ssp_tot = np.sum((medians_swirl - np.mean(medians_swirl))**2)
    r_squared_sp = 1 - (ssp_res / ssp_tot)
    print(f"R squared: {r_squared_sp:.4f}")
    
    adj_sp = 1-r_squared_sp
    adjr_sp = 1 - ((adj_sp*23)/22)
    print(f"Adjusted R squared: {adjr_sp:.4f}")
    
    swirl_grad_p = round(a_p[0],2)
    swirl_grad_err_p = np.sqrt(np.diag(acov_p))
    swirl_grad_err_p = round(swirl_grad_err_p[0], 2)
    
    swirl_perp_fit_p = (f'{swirl_grad_p} ± {swirl_grad_err_p}')
    
    
    a_m,acov_m = np.polyfit(btot_medians_m, medians_swirl,1, cov=True)
    aa_m = np.poly1d(a_m)
    # print(f"Swirl Fit -{error} Gradient:")
    # print(a_m[0])
    # print(f"Error of Swirl Fit -{error}:")
    # print(np.sqrt(np.diag(acov_m)))
    y_regress_swirl_m = aa_p(btot_medians_m)

    print("---------------")
    print('Swirl Region -20% Fit')
    ssm_res = np.sum((medians_swirl - y_regress_swirl_m)**2)
    ssm_tot = np.sum((medians_swirl - np.mean(medians_swirl))**2)
    r_squared_sm = 1 - (ssm_res / ssm_tot)
    print(f"R squared: {r_squared_sm:.4f}")
    
    adj_sm = 1-r_squared_sm
    adjr_sm = 1 - ((adj_sm*23)/22)
    print(f"Adjusted R squared: {adjr_sm:.4f}")
    
    swirl_grad_m = round(a_m[0],2)
    swirl_grad_err_m = np.sqrt(np.diag(acov_m))
    swirl_grad_err_m = round(swirl_grad_err_m[0], 2)
    
    swirl_perp_fit_m = (f'{swirl_grad_m} ± {swirl_grad_err_m}')
   
    
    b, bcov = np.polyfit(btot_medians, dusk_medians, 1, cov=True)
    bb = np.poly1d(b)
    # print("Dusk Fit Gradient:")
    # print(b[0])
    # print("Error of Dusk Fit:")
    # print(np.sqrt(np.diag(bcov)))
    y_regress_dusk = bb(btot_medians)

    print("---------------")
    print('Dusk Region Fit')
    sd_res = np.sum((dusk_medians - y_regress_dusk)**2)
    sd_tot = np.sum((dusk_medians - np.mean(dusk_medians))**2)
    r_squared_d = 1 - (sd_res / sd_tot)
    print(f"R squared: {r_squared_d:.4f}")
    
    adj_d = 1-r_squared_d
    adjr_d = 1 - ((adj_d*23)/22)
    print(f"Adjusted R squared: {adjr_d:.4f}")
    
    dusk_grad = round(b[0],2)
    dusk_grad_err = np.sqrt(np.diag(bcov))
    dusk_grad_err = round(dusk_grad_err[0], 2)
    
    dusk_perp_fit = (f'{dusk_grad} ± {dusk_grad_err}')
    
    
    b_p, bcov_p = np.polyfit(btot_medians_p, dusk_medians, 1, cov=True)
    bb_p = np.poly1d(b_p)
    # print(f"Dusk Fit +{error} Gradient:")
    # print(b_p[0])
    # print(f"Error of Dusk Fit +{error}:")
    # print(np.sqrt(np.diag(bcov_p)))
    y_regress_dusk_p = bb_p(btot_medians_p)

    print("---------------")
    print('Dusk Region +20% Fit')
    sdp_res = np.sum((dusk_medians - y_regress_dusk_p)**2)
    sdp_tot = np.sum((dusk_medians - np.mean(dusk_medians))**2)
    r_squared_dp = 1 - (sdp_res / sdp_tot)
    print(f"R squared: {r_squared_dp:.4f}")
    
    adj_dp = 1-r_squared_dp
    adjr_dp = 1 - ((adj_dp*23)/22)
    print(f"Adjusted R squared: {adjr_dp:.4f}")
    
    dusk_grad_p = round(b_p[0],2)
    dusk_grad_err_p = np.sqrt(np.diag(bcov_p))
    dusk_grad_err_p = round(dusk_grad_err_p[0], 2)
    
    dusk_perp_fit_p = (f'{dusk_grad_p} ± {dusk_grad_err_p}')
    
    
    b_m, bcov_m = np.polyfit(btot_medians_m, dusk_medians, 1, cov=True)
    bb_m = np.poly1d(b_m)
    # print(f"Dusk Fit -{error} Gradient:")
    # print(b_m[0])
    # print(f"Error of Dusk Fit -{error}:")
    # print(np.sqrt(np.diag(bcov_m)))
    y_regress_dusk_m = bb_m(btot_medians_m)

    print("---------------")
    print('Dusk Region -20% Fit')
    sdm_res = np.sum((dusk_medians - y_regress_dusk_m)**2)
    sdm_tot = np.sum((dusk_medians - np.mean(dusk_medians))**2)
    r_squared_dm = 1 - (sdm_res / sdm_tot)
    print(f"R squared: {r_squared_dm:.4f}")
    
    adj_dm = 1-r_squared_dm
    adjr_dm = 1 - ((adj_dm*23)/22)
    print(f"Adjusted R squared: {adjr_dm:.4f}")
    
    dusk_grad_m = round(b_m[0],2)
    dusk_grad_err_m = np.sqrt(np.diag(bcov_m))
    dusk_grad_err_m = round(dusk_grad_err_m[0], 2)
    
    dusk_perp_fit_m = (f'{dusk_grad_m} ± {dusk_grad_err_m}')
    
    
    
    c, ccov = np.polyfit(btot_medians, noon_medians, 1, cov=True)
    cc = np.poly1d(c)
    # print("Noon Fit Gradient:")
    # print(c[0])
    # print("Error of Noon Fit:")
    # print(np.sqrt(np.diag(ccov)))
    y_regress_noon = cc(btot_medians)

    print("---------------")
    print('Noon Region Fit')
    sn_res = np.sum((noon_medians - y_regress_noon)**2)
    sn_tot = np.sum((noon_medians - np.mean(noon_medians))**2)
    r_squared_n = 1 - (sn_res / sn_tot)
    print(f"R squared: {r_squared_n:.4f}")
    
    adj_n = 1-r_squared_n
    adjr_n = 1 - ((adj_n*23)/22)
    print(f"Adjusted R squared: {adjr_n:.4f}")
    
    noon_grad = round(c[0],2)
    noon_grad_err = np.sqrt(np.diag(ccov))
    noon_grad_err = round(noon_grad_err[0], 2)
    
    noon_perp_fit = (f'{noon_grad} ± {noon_grad_err}')
    
    
    c_p, ccov_p = np.polyfit(btot_medians_p, noon_medians, 1, cov=True)
    cc_p = np.poly1d(c_p)
    # print(f"Noon Fit +{error} Gradient:")
    # print(c_p[0])
    # print(f"Error of Noon Fit +{error}:")
    # print(np.sqrt(np.diag(ccov_p)))
    y_regress_noon_p = cc_p(btot_medians_p)

    print("---------------")
    print('Noon Region +20% Fit')
    snp_res = np.sum((noon_medians - y_regress_noon_p)**2)
    snp_tot = np.sum((noon_medians - np.mean(noon_medians))**2)
    r_squared_np = 1 - (snp_res / snp_tot)
    print(f"R squared: {r_squared_np:.4f}")
    
    adj_np = 1-r_squared_np
    adjr_np = 1 - ((adj_np*23)/22)
    print(f"Adjusted R squared: {adjr_np:.4f}")

    
    noon_grad_p = round(c_p[0],2)
    noon_grad_err_p = np.sqrt(np.diag(ccov_p))
    noon_grad_err_p = round(noon_grad_err_p[0], 2)
    
    noon_perp_fit_p = (f'{noon_grad_p} ± {noon_grad_err_p}')
    
    
    c_m, ccov_m = np.polyfit(btot_medians_m, noon_medians, 1, cov=True)
    cc_m = np.poly1d(c_m)
    # print(f"Noon Fit -{error} Gradient:")
    # print(c_m[0])
    # print(f"Error of Noon Fit -{error}:")
    # print(np.sqrt(np.diag(ccov_m)))
    y_regress_noon_m = cc_m(btot_medians_m)

    print("---------------")
    print('Noon Region -20% Fit')
    snm_res = np.sum((noon_medians - y_regress_noon_m)**2)
    snm_tot = np.sum((noon_medians - np.mean(noon_medians))**2)
    r_squared_nm = 1 - (snm_res / snm_tot)
    print(f"R squared: {r_squared_nm:.4f}")
    
    adj_nm = 1-r_squared_nm
    adjr_nm = 1 - ((adj_nm*23)/22)
    print(f"Adjusted R squared: {adjr_nm:.4f}")
    
    noon_grad_m = round(c_m[0],2)
    noon_grad_err_m = np.sqrt(np.diag(ccov_m))
    noon_grad_err_m = round(noon_grad_err_m[0], 2)
    
    noon_perp_fit_m = (f'{noon_grad_m} ± {noon_grad_err_m}')
    
    
    y_new_swirl = aa(x_new)
    y_new_dusk = bb(x_new)
    y_new_noon = cc(x_new)
    
    
    y_new_swirl_p = aa_p(x_new)
    y_new_dusk_p = bb_p(x_new)
    y_new_noon_p = cc_p(x_new)
    
    y_new_swirl_m = aa_m(x_new)
    y_new_dusk_m = bb_m(x_new)
    y_new_noon_m = cc_m(x_new)
    
    
    '''
    plotting
    '''
    
    fig = plt.figure(figsize=(28,50))
    ax1 = plt.subplot(4,2,1)
    plt.subplots_adjust(hspace=0.1, wspace=0.15)
    
    # main plot
    ax1.scatter(Bz_pos_btot[0:2], Bz_pos_power[0:2], color='red', s=200, label='+Bz CME')
    ax1.scatter(Bz_pos_btot[2:], Bz_pos_power[2:], color='red', s=200, marker='v', label='+Bz Rarefaction (Deep)')

    ax1.scatter(Bz_neg_btot[0], Bz_neg_power[0], color='orange', s=200, label='-Bz CME') #
    ax1.scatter(Bz_neg_btot[1], Bz_neg_power[1], color='orange', s=200,marker=',', label='-Bz CIR') 

    ax1.scatter(By_pos_btot[0], By_pos_power[0], color='green', s=200, label='+By CME') # compression 1 2-5 incl visit 15
    ax1.scatter(By_pos_btot[1:4], By_pos_power[1:4], color='green', s=200, marker=",", label='+By CIR')
    ax1.scatter(By_pos_btot[4:6], By_pos_power[4:6], color='green', s=200, marker="v", label='+By Rarefaction (Deep)') # rarefraction
    ax1.scatter(By_pos_btot[6], By_pos_power[6], color='green', s=200, marker="^", label ='+By Rarefaction (Shallow)')
    ax1.scatter(By_pos_btot[7:], By_pos_power[7:], color='green', s=200, marker=",")

    ax1.scatter(By_neg_btot[0:4], By_neg_power[0:4], color='blue', s=200, marker="^", label='-By Rarefaction (Shallow)')
    ax1.scatter(By_neg_btot[4:6], By_neg_power[4:6], color='blue', s=200, label='-By CME') # 
    ax1.scatter(By_neg_btot[6], By_neg_power[6], color='blue', s=200, marker="v", label='-By Rarefaction (Deep)')
    ax1.scatter(By_neg_btot[7:], By_neg_power[7:], color='blue', s=200, marker=",", label='-By CIR')

    ax1.errorbar(btot_medians, polar_medians, yerr=polar_med_errs, xerr=btot_med_sd_errs, fmt='.', color='lightgray')
    
    ax1.text(1.97,15,'a',style='italic',fontsize=40)
   
    ax1.set_xlabel('Median B$_{total}$ (nT)',fontsize=22)
    ax1.set_ylabel('Median Total Polar Power (GW)',fontsize=22)
   
    ax1.set_xlim(0,2.1)
    ax1.set_ylim(-20,1030)
    ax1.tick_params(axis='x',which='minor',direction='in')
    ax1.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax1.minorticks_on()
    ax1.tick_params(which='major',direction='in',bottom=True, top=True, left=True, right=True, labelsize=20, length=10)
    ax1.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True, length=5)
    lgnd = plt.legend(framealpha=0.5, labelspacing=0.3, fontsize=18, loc="upper left")
    for handle in lgnd.legend_handles:
        handle.set_sizes([100.0])
    
    
    ax2 = plt.subplot(4,2,2)
    
    #ax2.plot(x_new, y_new_swirl, '--')
    ax2.plot(x_new, y_new_swirl, '--',color='black', linewidth=3, markersize=12, label=(f'Swirl $R^2$: {adjr_s:.2g}'), zorder=10)
    
    
    if error_plot == '20' or error_plot == '10' or error_plot == '50':
        # plus error
        ax2.scatter(btot_medians_p[0:2], medians_swirl[0:2],s=100,color='plum',marker="^")
        ax2.scatter(btot_medians_p[2], medians_swirl[2],s=100,color='plum')
        ax2.scatter(btot_medians_p[3:5], medians_swirl[3:5],s=100,color='plum',marker="^")
        ax2.scatter(btot_medians_p[5:9], medians_swirl[5:9],s=100,color='plum')
        ax2.scatter(btot_medians_p[9], medians_swirl[9],s=100,color='plum',marker="v")
        ax2.scatter(btot_medians_p[10], medians_swirl[10],s=100,color='plum')
        ax2.scatter(btot_medians_p[11], medians_swirl[11],s=100,color='plum',marker="v")
        ax2.scatter(btot_medians_p[12:18], medians_swirl[12:18],s=100,color='plum',marker=",")
        ax2.scatter(btot_medians_p[18:21], medians_swirl[18:21],s=100,color='plum',marker="v")
        ax2.scatter(btot_medians_p[21], medians_swirl[21],s=100,color='plum',marker="^")
        ax2.scatter(btot_medians_p[22:], medians_swirl[22:],s=100,color='plum',marker=",")
    
        # minus error
        ax2.scatter(btot_medians_m[0:2], medians_swirl[0:2],s=100,color='peachpuff',marker="^")
        ax2.scatter(btot_medians_m[2], medians_swirl[2],s=100,color='peachpuff')
        ax2.scatter(btot_medians_m[3:5], medians_swirl[3:5],s=100,color='peachpuff',marker="^")
        ax2.scatter(btot_medians_m[5:9], medians_swirl[5:9],s=100,color='peachpuff')
        ax2.scatter(btot_medians_m[9], medians_swirl[9],s=100,color='peachpuff',marker="v")
        ax2.scatter(btot_medians_m[10], medians_swirl[10],s=100,color='peachpuff')
        ax2.scatter(btot_medians_m[11], medians_swirl[11],s=100,color='peachpuff',marker="v")
        ax2.scatter(btot_medians_m[12:18], medians_swirl[12:18],s=100,color='peachpuff',marker=",")
        ax2.scatter(btot_medians_m[18:21], medians_swirl[18:21],s=100,color='peachpuff',marker="v")
        ax2.scatter(btot_medians_m[21], medians_swirl[21],s=100,color='peachpuff',marker="^")
        ax2.scatter(btot_medians_m[22:], medians_swirl[22:],s=100,color='peachpuff',marker=",")
        
        ax2.plot(x_new, y_new_swirl_p, '--',color='darkviolet', linewidth=3, markersize=12, label=f'+ {error}%: {adjr_sp:.2g}')
        ax2.plot(x_new, y_new_swirl_m, '--',color='sandybrown', linewidth=3, markersize=12, label=f'- {error}%: {adjr_sm:.2g}')
        
        
        ax2.errorbar(btot_medians_p, medians_swirl, yerr=swirl_med_errs, xerr=btot_med_sd_errs_p, fmt='.', color='plum')
        
        ax2.errorbar(btot_medians_m, medians_swirl, yerr=swirl_med_errs, xerr=btot_med_sd_errs_m, fmt='.', color='peachpuff')
    
    
    # main plot
    ax2.scatter(Bz_pos_btot[0:2], Bz_pos_swirlpower[0:2], color='red', s=200,zorder=9)
    ax2.scatter(Bz_pos_btot[2:], Bz_pos_swirlpower[2:], color='red', s=200, marker='v',zorder=9)

    ax2.scatter(Bz_neg_btot[0], Bz_neg_swirlpower[0], color='orange', s=200,zorder=9) #
    ax2.scatter(Bz_neg_btot[1], Bz_neg_swirlpower[1], color='orange', s=200,marker=',',zorder=9) 

    ax2.scatter(By_pos_btot[0], By_pos_swirlpower[0], color='green', s=200,zorder=9) # compression 1 2-5 incl visit 15
    ax2.scatter(By_pos_btot[1:4], By_pos_swirlpower[1:4], color='green', s=200, marker=",",zorder=9)
    ax2.scatter(By_pos_btot[4:6], By_pos_swirlpower[4:6], color='green', s=200, marker="v",zorder=9) # rarefraction
    ax2.scatter(By_pos_btot[6], By_pos_swirlpower[6], color='green', s=200, marker="^",zorder=9)
    ax2.scatter(By_pos_btot[7:], By_pos_swirlpower[7:], color='green', s=200, marker=",",zorder=9)

    ax2.scatter(By_neg_btot[0:4], By_neg_swirlpower[0:4], color='blue', s=200, marker="^",zorder=9)
    ax2.scatter(By_neg_btot[4:6], By_neg_swirlpower[4:6], color='blue', s=200,zorder=9) # 
    ax2.scatter(By_neg_btot[6], By_neg_swirlpower[6], color='blue', s=200, marker="v",zorder=9)
    ax2.scatter(By_neg_btot[7:], By_neg_swirlpower[7:], color='blue', s=200, marker=",",zorder=9)
    
    ax2.errorbar(btot_medians, medians_swirl, yerr=swirl_med_errs, xerr=btot_med_sd_errs, fmt='.', color='lightgray',zorder=9)
    
    ax2.text(1.97,-1.8,'b',style='italic',fontsize=40) #-1.5, -2.1
    
    ax2.set_xlabel('Median B$_{total}$ (nT)',fontsize=22)
    ax2.set_ylabel('Median Swirl Region Power (GW)',fontsize=22)
   
    ax2.set_xlim(0,2.1)
    ax2.set_ylim(-5,90)#85 74
    ax2.tick_params(axis='x',which='minor',direction='in')
    ax2.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax2.minorticks_on()
    ax2.tick_params(which='major',direction='in',bottom=True, top=True, left=True, right=True, labelsize=20, length=10)
    ax2.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True, length=5)
    ax2.legend(framealpha=0.5,fontsize=21, loc ="upper left")
    
    #########
    
    ax3 = plt.subplot(4,2,3)
    
    ax3.plot(x_new, y_new_dusk,'--', color='black', linewidth=3, markersize=12, label=f'Dusk $R^2$: {adjr_d:.2g}', zorder=10)
    
    
    if error_plot == '20' or error_plot == '10' or error_plot == '50':
        ax3.plot(x_new, y_new_dusk_p, '--',color='darkviolet', linewidth=3, markersize=12,label=f'+ {error}%: {adjr_dp:.2g}')
        ax3.plot(x_new, y_new_dusk_m, '--',color='sandybrown', linewidth=3, markersize=12, label=f'- {error}%: {adjr_dm:.2g}')
        
        ax3.errorbar(btot_medians_p, dusk_medians, yerr=dusk_med_errs, xerr=btot_med_sd_errs_p, fmt='.', color='plum')
        ax3.errorbar(btot_medians_m, dusk_medians, yerr=dusk_med_errs, xerr=btot_med_sd_errs_m, fmt='.', color='peachpuff')
    
    
        ax3.scatter(btot_medians_p[0:2], dusk_medians[0:2],s=100,color='plum',marker="^")
        ax3.scatter(btot_medians_p[2], dusk_medians[2],s=100,color='plum')
        ax3.scatter(btot_medians_p[3:5], dusk_medians[3:5],s=100,color='plum',marker="^")
        ax3.scatter(btot_medians_p[5:9], dusk_medians[5:9],s=100,color='plum')
        ax3.scatter(btot_medians_p[9], dusk_medians[9],s=100,color='plum',marker="v")
        ax3.scatter(btot_medians_p[10], dusk_medians[10],s=100,color='plum')
        ax3.scatter(btot_medians_p[11], dusk_medians[11],s=100,color='plum',marker="v")
        ax3.scatter(btot_medians_p[12:18], dusk_medians[12:18],s=100,color='plum',marker=",")
        ax3.scatter(btot_medians_p[18:21], dusk_medians[18:21],s=100,color='plum',marker="v")
        ax3.scatter(btot_medians_p[21], dusk_medians[21],s=100,color='plum',marker="^")
        ax3.scatter(btot_medians_p[22:], dusk_medians[22:],s=100,color='plum',marker=",")
    
        # minus error
        ax3.scatter(btot_medians_m[0:2], dusk_medians[0:2],s=100,color='peachpuff',marker="^")
        ax3.scatter(btot_medians_m[2], dusk_medians[2],s=100,color='peachpuff')
        ax3.scatter(btot_medians_m[3:5], dusk_medians[3:5],s=100,color='peachpuff',marker="^")
        ax3.scatter(btot_medians_m[5:9], dusk_medians[5:9],s=100,color='peachpuff')
        ax3.scatter(btot_medians_m[9], dusk_medians[9],s=100,color='peachpuff',marker="v")
        ax3.scatter(btot_medians_m[10], dusk_medians[10],s=100,color='peachpuff')
        ax3.scatter(btot_medians_m[11], dusk_medians[11],s=100,color='peachpuff',marker="v")
        ax3.scatter(btot_medians_m[12:18], dusk_medians[12:18],s=100,color='peachpuff',marker=",")
        ax3.scatter(btot_medians_m[18:21], dusk_medians[18:21],s=100,color='peachpuff',marker="v")
        ax3.scatter(btot_medians_m[21], dusk_medians[21],s=100,color='peachpuff',marker="^")
        ax3.scatter(btot_medians_m[22:], dusk_medians[22:],s=100,color='peachpuff',marker=",")
    
    
    # main plot
    ax3.scatter(Bz_pos_btot[0:2], Bz_pos_duskpower[0:2], color='red', s=200,zorder=9)
    ax3.scatter(Bz_pos_btot[2:], Bz_pos_duskpower[2:], color='red', s=200, marker='v',zorder=9)

    ax3.scatter(Bz_neg_btot[0], Bz_neg_duskpower[0], color='orange', s=200,zorder=9) #
    ax3.scatter(Bz_neg_btot[1], Bz_neg_duskpower[1], color='orange', s=200,marker=',',zorder=9)

    ax3.scatter(By_pos_btot[0], By_pos_duskpower[0], color='green', s=200,zorder=9) # compression 1 2-5 incl visit 15
    ax3.scatter(By_pos_btot[1:4], By_pos_duskpower[1:4], color='green', s=200, marker=",",zorder=9)
    ax3.scatter(By_pos_btot[4:6], By_pos_duskpower[4:6], color='green', s=200, marker="v",zorder=9) # rarefraction
    ax3.scatter(By_pos_btot[6], By_pos_duskpower[6], color='green', s=200, marker="^",zorder=9)
    ax3.scatter(By_pos_btot[7:], By_pos_duskpower[7:], color='green', s=200, marker=",",zorder=9)

    ax3.scatter(By_neg_btot[0:4], By_neg_duskpower[0:4], color='blue', s=200, marker="^",zorder=9)
    ax3.scatter(By_neg_btot[4:6], By_neg_duskpower[4:6], color='blue', s=200,zorder=9) # 
    ax3.scatter(By_neg_btot[6], By_neg_duskpower[6], color='blue', s=200, marker="v",zorder=9)
    ax3.scatter(By_neg_btot[7:], By_neg_duskpower[7:], color='blue', s=200, marker=",",zorder=9)
    
    ax3.errorbar(btot_medians, dusk_medians, yerr=dusk_med_errs, xerr=btot_med_sd_errs, fmt='.', color='lightgray', zorder=9)
    
    ax3.text(1.97,0,'c',style='italic',fontsize=40)

    ax3.set_xlabel('Median B$_{total}$ (nT)',fontsize=22)
    ax3.set_ylabel('Median Dusk Active Region Power (GW)',fontsize=22)
    
    ax3.set_xlim(0,2.1)
    ax3.set_ylim(-20,510) #500 #440
    ax3.tick_params(axis='x',which='minor',direction='in')
    ax3.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax3.minorticks_on()
    ax3.tick_params(which='major',direction='in',bottom=True, top=True, left=True, right=True, labelsize=20, length=10)
    ax3.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True, length=5)
    ax3.legend(framealpha=0.5,fontsize=21,loc ="upper left")
    
    
    ########
    
    ax4 = plt.subplot(4,2,4)
    
    ax4.plot(x_new, y_new_noon, '--',color='black', linewidth=3, markersize=12, label=f'Noon $R^2$: {adjr_n:.2g}', zorder=10)
    
    if error_plot == '20' or error_plot == '10' or error_plot == '50':
        ax4.plot(x_new, y_new_noon_p, '--',color='darkviolet', linewidth=3, markersize=12,label=f'+ {error}%: {adjr_np:.2g}')
        ax4.plot(x_new, y_new_noon_m, '--',color='sandybrown', linewidth=3, markersize=12,label=f'- {error}%: {adjr_nm:.2g}')
    
        
        ax4.errorbar(btot_medians_p, noon_medians, yerr=noon_med_errs, xerr=btot_med_sd_errs_p, fmt='.', color='plum')
        ax4.errorbar(btot_medians_m, noon_medians, yerr=noon_med_errs, xerr=btot_med_sd_errs_m, fmt='.', color='peachpuff')
        
        # plus error
        ax4.scatter(btot_medians_p[0:2], noon_medians[0:2],s=100,color='plum',marker="^")
        ax4.scatter(btot_medians_p[2], noon_medians[2],s=100,color='plum')
        ax4.scatter(btot_medians_p[3:5], noon_medians[3:5],s=100,color='plum',marker="^")
        ax4.scatter(btot_medians_p[5:9], noon_medians[5:9],s=100,color='plum')
        ax4.scatter(btot_medians_p[9], noon_medians[9],s=100,color='plum',marker="v")
        ax4.scatter(btot_medians_p[10], noon_medians[10],s=100,color='plum')
        ax4.scatter(btot_medians_p[11], noon_medians[11],s=100,color='plum',marker="v")
        ax4.scatter(btot_medians_p[12:18], noon_medians[12:18],s=100,color='plum',marker=",")
        ax4.scatter(btot_medians_p[18:21], noon_medians[18:21],s=100,color='plum',marker="v")
        ax4.scatter(btot_medians_p[21], noon_medians[21],s=100,color='plum',marker="^")
        ax4.scatter(btot_medians_p[22:], noon_medians[22:],s=100,color='plum',marker=",")
    
        # minus error
        ax4.scatter(btot_medians_m[0:2], noon_medians[0:2],s=100,color='peachpuff',marker="^")
        ax4.scatter(btot_medians_m[2], noon_medians[2],s=100,color='peachpuff')
        ax4.scatter(btot_medians_m[3:5], noon_medians[3:5],s=100,color='peachpuff',marker="^")
        ax4.scatter(btot_medians_m[5:9], noon_medians[5:9],s=100,color='peachpuff')
        ax4.scatter(btot_medians_m[9], noon_medians[9],s=100,color='peachpuff',marker="v")
        ax4.scatter(btot_medians_m[10], noon_medians[10],s=100,color='peachpuff')
        ax4.scatter(btot_medians_m[11], noon_medians[11],s=100,color='peachpuff',marker="v")
        ax4.scatter(btot_medians_m[12:18], noon_medians[12:18],s=100,color='peachpuff',marker=",")
        ax4.scatter(btot_medians_m[18:21], noon_medians[18:21],s=100,color='peachpuff',marker="v")
        ax4.scatter(btot_medians_m[21], noon_medians[21],s=100,color='peachpuff',marker="^")
        ax4.scatter(btot_medians_m[22:], noon_medians[22:],s=100,color='peachpuff',marker=",")
    
    # main plot
    ax4.scatter(Bz_pos_btot[0:2], Bz_pos_noonpower[0:2], color='red', s=200,zorder=9)
    ax4.scatter(Bz_pos_btot[2:], Bz_pos_noonpower[2:], color='red', s=200, marker='v',zorder=9)

    ax4.scatter(Bz_neg_btot[0], Bz_neg_noonpower[0], color='orange', s=200,zorder=9) #
    ax4.scatter(Bz_neg_btot[1], Bz_neg_noonpower[1], color='orange', s=200,marker=',',zorder=9) 

    ax4.scatter(By_pos_btot[0], By_pos_noonpower[0], color='green', s=200,zorder=9) # compression 1 2-5 incl visit 15
    ax4.scatter(By_pos_btot[1:4], By_pos_noonpower[1:4], color='green', s=200, marker=",",zorder=9)
    ax4.scatter(By_pos_btot[4:6], By_pos_noonpower[4:6], color='green', s=200, marker="v",zorder=9) # rarefraction
    ax4.scatter(By_pos_btot[6], By_pos_noonpower[6], color='green', s=200, marker="^",zorder=9)
    ax4.scatter(By_pos_btot[7:], By_pos_noonpower[7:], color='green', s=200, marker=",",zorder=9)

    ax4.scatter(By_neg_btot[0:4], By_neg_noonpower[0:4], color='blue', s=200, marker="^",zorder=9)
    ax4.scatter(By_neg_btot[4:6], By_neg_noonpower[4:6], color='blue', s=200,zorder=9)
    ax4.scatter(By_neg_btot[6], By_neg_noonpower[6], color='blue', s=200, marker="v",zorder=9)
    ax4.scatter(By_neg_btot[7:], By_neg_noonpower[7:], color='blue', s=200, marker=",",zorder=9)
    
    
    ax4.errorbar(btot_medians, noon_medians, yerr=noon_med_errs, xerr=btot_med_sd_errs, fmt='.', color='lightgray', zorder=9)
    
    ax4.text(1.97,0,'d',style='italic',fontsize=40) #-3

    ax4.set_ylabel('Median Noon Active Region Power (GW)',fontsize=22)
    ax4.set_xlabel('Median B$_{total}$ (nT)',fontsize=22)
    
    ax4.tick_params(axis='x',which='minor',direction='in')
    ax4.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax4.minorticks_on()
    ax4.tick_params(which='major',direction='in',bottom=True, top=True, left=True, right=True, labelsize=20, length=10)
    ax4.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True, length=5)
    
    ax4.set_xlim(0,2.1)
    ax4.set_ylim(-20,510) #500 440
    
    ax4.legend(framealpha=0.5,fontsize=21,loc ="upper left")
    
    # save plot
    saveloc = (f'{root_saves}median_btot_vs_power_{error_plot}.jpg') 
    #saveloc = (f'/Users/hannah/OneDrive - Lancaster University/aurora/median_btot_vs_power_{error}.jpg') #remember to add f' at front when code {error} back on
    plt.savefig(saveloc,bbox_inches='tight',dpi=400)

if plotting == 'b_perp':
    
    x_new = np.linspace(min(b_perp_medians),max(b_perp_medians),len(b_perp_medians))
    
    a, acov = np.polyfit(b_perp_medians, medians_swirl, 1, cov=True)
    aa = np.poly1d(a)
    # print("Swirl Fit Gradient:")
    # print(a[0])
    # print("Error of Swirl Fit:")
    # print(np.sqrt(np.diag(acov)))
    y_regress_swirl = aa(b_perp_medians)

    print(f"Testing '{plotting}'")
    print('Swirl Region Fit')
    ss_res = np.sum((medians_swirl - y_regress_swirl)**2)
    ss_tot = np.sum((medians_swirl - np.mean(medians_swirl))**2)
    r_squared_s = 1 - (ss_res / ss_tot)
    print(f"R squared: {r_squared_s:.4f}")
    
    adj_s = 1-r_squared_s
    adjr_s = 1 - ((adj_s*23)/22)
    print(f"Adjusted R squared: {adjr_s:.4f}")
    
    swirl_grad = round(a[0],2)
    swirl_grad_err = np.sqrt(np.diag(acov))
    swirl_grad_err = round(swirl_grad_err[0], 2)
    
    swirl_perp_fit = (f'{swirl_grad} ± {swirl_grad_err}')
    
    
    a_p,acov_p = np.polyfit(b_perp_medians_plus, medians_swirl,1, cov=True)
    aa_p = np.poly1d(a_p)
    # print(f"Swirl Fit +{error} Gradient:")
    # print(a_p[0])
    # print(f"Error of Swirl Fit +{error}:")
    # print(np.sqrt(np.diag(acov_p)))
    y_regress_swirl_p = aa_p(b_perp_medians_plus)

    print("---------------")
    print('Swirl Region +20% Fit')
    ssp_res = np.sum((medians_swirl - y_regress_swirl_p)**2)
    ssp_tot = np.sum((medians_swirl - np.mean(medians_swirl))**2)
    r_squared_sp = 1 - (ssp_res / ssp_tot)
    print(f"R squared: {r_squared_sp:.4f}")
    
    adj_sp = 1-r_squared_sp
    adjr_sp = 1 - ((adj_sp*23)/22)
    print(f"Adjusted R squared: {adjr_sp:.4f}")
    
    swirl_grad_p = round(a_p[0],2)
    swirl_grad_err_p = np.sqrt(np.diag(acov_p))
    swirl_grad_err_p = round(swirl_grad_err_p[0], 2)
    
    swirl_perp_fit_p = (f'{swirl_grad_p} ± {swirl_grad_err_p}')
    
    
    a_m,acov_m = np.polyfit(b_perp_medians_minus, medians_swirl,1, cov=True)
    aa_m = np.poly1d(a_m)
    # print(f"Swirl Fit -{error} Gradient:")
    # print(a_m[0])
    # print(f"Error of Swirl Fit -{error}:")
    # print(np.sqrt(np.diag(acov_m)))
    y_regress_swirl_m = aa_p(b_perp_medians_minus)

    print("---------------")
    print('Swirl Region -20% Fit')
    ssm_res = np.sum((medians_swirl - y_regress_swirl_m)**2)
    ssm_tot = np.sum((medians_swirl - np.mean(medians_swirl))**2)
    r_squared_sm = 1 - (ssm_res / ssm_tot)
    print(f"R squared: {r_squared_sm:.4f}")
    
    adj_sm = 1-r_squared_sm
    adjr_sm = 1 - ((adj_sm*23)/22)
    print(f"Adjusted R squared: {adjr_sm:.4f}")
    
    swirl_grad_m = round(a_m[0],2)
    swirl_grad_err_m = np.sqrt(np.diag(acov_m))
    swirl_grad_err_m = round(swirl_grad_err_m[0], 2)
    
    swirl_perp_fit_m = (f'{swirl_grad_m} ± {swirl_grad_err_m}')
   
    
    b, bcov = np.polyfit(b_perp_medians, dusk_medians, 1, cov=True)
    bb = np.poly1d(b)
    # print("Dusk Fit Gradient:")
    # print(b[0])
    # print("Error of Dusk Fit:")
    # print(np.sqrt(np.diag(bcov)))
    y_regress_dusk = bb(b_perp_medians)

    print("---------------")
    print('Dusk Region Fit')
    sd_res = np.sum((dusk_medians - y_regress_dusk)**2)
    sd_tot = np.sum((dusk_medians - np.mean(dusk_medians))**2)
    r_squared_d = 1 - (sd_res / sd_tot)
    print(f"R squared: {r_squared_d:.4f}")
    
    adj_d = 1-r_squared_d
    adjr_d = 1 - ((adj_d*23)/22)
    print(f"Adjusted R squared: {adjr_d:.4f}")
    
    dusk_grad = round(b[0],2)
    dusk_grad_err = np.sqrt(np.diag(bcov))
    dusk_grad_err = round(dusk_grad_err[0], 2)
    
    dusk_perp_fit = (f'{dusk_grad} ± {dusk_grad_err}')
    
    
    b_p, bcov_p = np.polyfit(b_perp_medians_plus, dusk_medians, 1, cov=True)
    bb_p = np.poly1d(b_p)
    # print(f"Dusk Fit +{error} Gradient:")
    # print(b_p[0])
    # print(f"Error of Dusk Fit +{error}:")
    # print(np.sqrt(np.diag(bcov_p)))
    y_regress_dusk_p = bb_p(b_perp_medians_plus)
    
    print("---------------")
    print('Dusk Region +20% Fit')
    sdp_res = np.sum((dusk_medians - y_regress_dusk_p)**2)
    sdp_tot = np.sum((dusk_medians - np.mean(dusk_medians))**2)
    r_squared_dp = 1 - (sdp_res / sdp_tot)
    print(f"R squared: {r_squared_dp:.4f}")
    
    adj_dp = 1-r_squared_dp
    adjr_dp = 1 - ((adj_dp*23)/22)
    print(f"Adjusted R squared: {adjr_dp:.4f}")
    
    dusk_grad_p = round(b_p[0],2)
    dusk_grad_err_p = np.sqrt(np.diag(bcov_p))
    dusk_grad_err_p = round(dusk_grad_err_p[0], 2)
    
    dusk_perp_fit_p = (f'{dusk_grad_p} ± {dusk_grad_err_p}')
    
    
    b_m, bcov_m = np.polyfit(b_perp_medians_minus, dusk_medians, 1, cov=True)
    bb_m = np.poly1d(b_m)
    # print(f"Dusk Fit -{error} Gradient:")
    # print(b_m[0])
    # print(f"Error of Dusk Fit -{error}:")
    # print(np.sqrt(np.diag(bcov_m)))
    y_regress_dusk_m = bb_m(b_perp_medians_minus)

    print("---------------")
    print('Dusk Region -20% Fit')
    sdm_res = np.sum((dusk_medians - y_regress_dusk_m)**2)
    sdm_tot = np.sum((dusk_medians - np.mean(dusk_medians))**2)
    r_squared_dm = 1 - (sdm_res / sdm_tot)
    print(f"R squared: {r_squared_dm:.4f}")
    
    adj_dm = 1-r_squared_dm
    adjr_dm = 1 - ((adj_dm*23)/22)
    print(f"Adjusted R squared: {adjr_dm:.4f}")
    
    dusk_grad_m = round(b_m[0],2)
    dusk_grad_err_m = np.sqrt(np.diag(bcov_m))
    dusk_grad_err_m = round(dusk_grad_err_m[0], 2)
    
    dusk_perp_fit_m = (f'{dusk_grad_m} ± {dusk_grad_err_m}')
    
    
    
    c, ccov = np.polyfit(b_perp_medians, noon_medians, 1, cov=True)
    cc = np.poly1d(c)
    # print("Noon Fit Gradient:")
    # print(c[0])
    # print("Error of Noon Fit:")
    # print(np.sqrt(np.diag(ccov)))
    y_regress_noon = cc(b_perp_medians)

    print("---------------")
    print('Noon Region Fit')
    sn_res = np.sum((noon_medians - y_regress_noon)**2)
    sn_tot = np.sum((noon_medians - np.mean(noon_medians))**2)
    r_squared_n = 1 - (sn_res / sn_tot)
    print(f"R squared: {r_squared_n:.4f}")
    
    adj_n = 1-r_squared_n
    adjr_n = 1 - ((adj_n*23)/22)
    print(f"Adjusted R squared: {adjr_n:.4f}")
    
    noon_grad = round(c[0],2)
    noon_grad_err = np.sqrt(np.diag(ccov))
    noon_grad_err = round(noon_grad_err[0], 2)
    
    noon_perp_fit = (f'{noon_grad} ± {noon_grad_err}')
    
    
    c_p, ccov_p = np.polyfit(b_perp_medians_plus, noon_medians, 1, cov=True)
    cc_p = np.poly1d(c_p)
    # print(f"Noon Fit +{error} Gradient:")
    # print(c_p[0])
    # print(f"Error of Noon Fit +{error}:")
    # print(np.sqrt(np.diag(ccov_p)))
    y_regress_noon_p = cc_p(b_perp_medians_plus)

    print("---------------")
    print('Noon Region +20% Fit')
    snp_res = np.sum((noon_medians - y_regress_noon_p)**2)
    snp_tot = np.sum((noon_medians - np.mean(noon_medians))**2)
    r_squared_np = 1 - (snp_res / snp_tot)
    print(f"R squared: {r_squared_np:.4f}")
    
    adj_np = 1-r_squared_np
    adjr_np = 1 - ((adj_np*23)/22)
    print(f"Adjusted R squared: {adjr_np:.4f}")
    
    noon_grad_p = round(c_p[0],2)
    noon_grad_err_p = np.sqrt(np.diag(ccov_p))
    noon_grad_err_p = round(noon_grad_err_p[0], 2)
    
    noon_perp_fit_p = (f'{noon_grad_p} ± {noon_grad_err_p}')
    
    
    c_m, ccov_m = np.polyfit(b_perp_medians_minus, noon_medians, 1, cov=True)
    cc_m = np.poly1d(c_m)
    # print(f"Noon Fit -{error} Gradient:")
    # print(c_m[0])
    # print(f"Error of Noon Fit -{error}:")
    # print(np.sqrt(np.diag(ccov_m)))
    y_regress_noon_m = cc_m(b_perp_medians_minus)

    print("---------------")
    print('Noon Region -20% Fit')
    snm_res = np.sum((noon_medians - y_regress_noon_m)**2)
    snm_tot = np.sum((noon_medians - np.mean(noon_medians))**2)
    r_squared_nm = 1 - (snm_res / snm_tot)
    print(f"R squared: {r_squared_nm:.4f}")
    
    adj_nm = 1-r_squared_nm
    adjr_nm = 1 - ((adj_nm*23)/22)
    print(f"Adjusted R squared: {adjr_nm:.4f}")
    
    noon_grad_m = round(c_m[0],2)
    noon_grad_err_m = np.sqrt(np.diag(ccov_m))
    noon_grad_err_m = round(noon_grad_err_m[0], 2)
    
    noon_perp_fit_m = (f'{noon_grad_m} ± {noon_grad_err_m}')
    
    
    y_new_swirl = aa(x_new)
    y_new_dusk = bb(x_new)
    y_new_noon = cc(x_new)
    
    
    y_new_swirl_p = aa_p(x_new)
    y_new_dusk_p = bb_p(x_new)
    y_new_noon_p = cc_p(x_new)
    
    y_new_swirl_m = aa_m(x_new)
    y_new_dusk_m = bb_m(x_new)
    y_new_noon_m = cc_m(x_new)
    
    
    '''
    plotting
    '''
    
    fig = plt.figure(figsize=(28,50))
    ax1 = plt.subplot(4,2,1)
    plt.subplots_adjust(hspace=0.1, wspace=0.15)
    
    # main plot
    ax1.scatter(Bz_pos_b_perp[0:2], Bz_pos_power[0:2], color='red', s=200, label='+Bz CME')
    ax1.scatter(Bz_pos_b_perp[2:], Bz_pos_power[2:], color='red', s=200, marker='v', label='+Bz Rarefaction (Deep)')

    ax1.scatter(Bz_neg_b_perp[0], Bz_neg_power[0], color='orange', s=200, label='-Bz CME') #
    ax1.scatter(Bz_neg_b_perp[1], Bz_neg_power[1], color='orange', s=200,marker=',', label='-Bz CIR') 

    ax1.scatter(By_pos_b_perp[0], By_pos_power[0], color='green', s=200, label='+By CME') # compression 1 2-5 incl visit 15
    ax1.scatter(By_pos_b_perp[1:4], By_pos_power[1:4], color='green', s=200, marker=",", label='+By CIR')
    ax1.scatter(By_pos_b_perp[4:6], By_pos_power[4:6], color='green', s=200, marker="v", label='+By Rarefaction (Deep)') # rarefraction
    ax1.scatter(By_pos_b_perp[6], By_pos_power[6], color='green', s=200, marker="^", label ='+By Rarefaction (Shallow)')
    ax1.scatter(By_pos_b_perp[7:], By_pos_power[7:], color='green', s=200, marker=",")

    ax1.scatter(By_neg_b_perp[0:4], By_neg_power[0:4], color='blue', s=200, marker="^", label='-By Rarefaction (Shallow)')
    ax1.scatter(By_neg_b_perp[4:6], By_neg_power[4:6], color='blue', s=200, label='-By CME') # 
    ax1.scatter(By_neg_b_perp[6], By_neg_power[6], color='blue', s=200, marker="v", label='-By Rarefaction (Deep)')
    ax1.scatter(By_neg_b_perp[7:], By_neg_power[7:], color='blue', s=200, marker=",", label='-By CIR')
    
    ax1.errorbar(b_perp_medians, polar_medians, yerr=polar_med_errs, xerr=med_b_perp_errs, fmt='.', color='lightgray')
    
    ax1.text(1.92,3,'a',style='italic',fontsize=40)
    
    ax1.set_xlabel('Median B$_⊥$ (nT)',fontsize=22)
    ax1.set_ylabel('Median Total Polar Power (GW)',fontsize=22)

    ax1.set_xlim(0.0,2.05)
    ax1.set_ylim(-30,1040)
    ax1.tick_params(axis='x',which='minor',direction='in')
    ax1.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax1.minorticks_on()
    ax1.tick_params(which='major',direction='in',bottom=True, top=True, left=True, right=True, labelsize=20, length=10)
    ax1.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True, length=5)
    lgnd = plt.legend(framealpha=0.5, labelspacing=0.3, fontsize=18, loc="upper left")
    for handle in lgnd.legend_handles:
        handle.set_sizes([100.0])
    
    
    ax2 = plt.subplot(4,2,2)
    
    ax2.plot(x_new, y_new_swirl, '--',color='black', linewidth=3, markersize=12, label=(f'Swirl $R^2$: {adjr_s:.2g}'), zorder=10)
    
    if error_plot == '20' or error_plot == '10' or error_plot == '50':
        ax2.errorbar(b_perp_medians_plus, medians_swirl, yerr=swirl_med_errs, xerr=med_b_perp_errs_plus, fmt='.', color='plum')
        
        ax2.errorbar(b_perp_medians_minus, medians_swirl, yerr=swirl_med_errs, xerr=med_b_perp_errs_minus, fmt='.', color='peachpuff')
        
        ax2.plot(x_new, y_new_swirl_p, '--',color='darkviolet', linewidth=3, markersize=12, label=f'+ {error}%: {adjr_sp:.2g}')
        ax2.plot(x_new, y_new_swirl_m, '--',color='sandybrown', linewidth=3, markersize=12, label=f'- {error}%: {adjr_sm:.2g}')
        
        # plus error
        ax2.scatter(b_perp_medians_plus[0:2], medians_swirl[0:2],s=100,color='plum',marker="^")
        ax2.scatter(b_perp_medians_plus[2], medians_swirl[2],s=100,color='plum')
        ax2.scatter(b_perp_medians_plus[3:5], medians_swirl[3:5],s=100,color='plum',marker="^")
        ax2.scatter(b_perp_medians_plus[5:9], medians_swirl[5:9],s=100,color='plum')
        ax2.scatter(b_perp_medians_plus[9], medians_swirl[9],s=100,color='plum',marker="v")
        ax2.scatter(b_perp_medians_plus[10], medians_swirl[10],s=100,color='plum')
        ax2.scatter(b_perp_medians_plus[11], medians_swirl[11],s=100,color='plum',marker="v")
        ax2.scatter(b_perp_medians_plus[12:18], medians_swirl[12:18],s=100,color='plum',marker=",")
        ax2.scatter(b_perp_medians_plus[18:21], medians_swirl[18:21],s=100,color='plum',marker="v")
        ax2.scatter(b_perp_medians_plus[21], medians_swirl[21],s=100,color='plum',marker="^")
        ax2.scatter(b_perp_medians_plus[22:], medians_swirl[22:],s=100,color='plum',marker=",")
    
        # minus error
        ax2.scatter(b_perp_medians_minus[0:2], medians_swirl[0:2],s=100,color='peachpuff',marker="^")
        ax2.scatter(b_perp_medians_minus[2], medians_swirl[2],s=100,color='peachpuff')
        ax2.scatter(b_perp_medians_minus[3:5], medians_swirl[3:5],s=100,color='peachpuff',marker="^")
        ax2.scatter(b_perp_medians_minus[5:9], medians_swirl[5:9],s=100,color='peachpuff')
        ax2.scatter(b_perp_medians_minus[9], medians_swirl[9],s=100,color='peachpuff',marker="v")
        ax2.scatter(b_perp_medians_minus[10], medians_swirl[10],s=100,color='peachpuff')
        ax2.scatter(b_perp_medians_minus[11], medians_swirl[11],s=100,color='peachpuff',marker="v")
        ax2.scatter(b_perp_medians_minus[12:18], medians_swirl[12:18],s=100,color='peachpuff',marker=",")
        ax2.scatter(b_perp_medians_minus[18:21], medians_swirl[18:21],s=100,color='peachpuff',marker="v")
        ax2.scatter(b_perp_medians_minus[21], medians_swirl[21],s=100,color='peachpuff',marker="^")
        ax2.scatter(b_perp_medians_minus[22:], medians_swirl[22:],s=100,color='peachpuff',marker=",")
    
    # main plot
    ax2.scatter(Bz_pos_b_perp[0:2], Bz_pos_swirlpower[0:2], color='red', s=200,zorder=9)
    ax2.scatter(Bz_pos_b_perp[2:], Bz_pos_swirlpower[2:], color='red', s=200, marker='v',zorder=9)

    ax2.scatter(Bz_neg_b_perp[0], Bz_neg_swirlpower[0], color='orange', s=200,zorder=9) #
    ax2.scatter(Bz_neg_b_perp[1], Bz_neg_swirlpower[1], color='orange', s=200,marker=',',zorder=9) 

    ax2.scatter(By_pos_b_perp[0], By_pos_swirlpower[0], color='green', s=200,zorder=9) # compression 1 2-5 incl visit 15
    ax2.scatter(By_pos_b_perp[1:4], By_pos_swirlpower[1:4], color='green', s=200, marker=",",zorder=9)
    ax2.scatter(By_pos_b_perp[4:6], By_pos_swirlpower[4:6], color='green', s=200, marker="v",zorder=9) # rarefraction
    ax2.scatter(By_pos_b_perp[6], By_pos_swirlpower[6], color='green', s=200, marker="^",zorder=9)
    ax2.scatter(By_pos_b_perp[7:], By_pos_swirlpower[7:], color='green', s=200, marker=",",zorder=9)

    ax2.scatter(By_neg_b_perp[0:4], By_neg_swirlpower[0:4], color='blue', s=200, marker="^",zorder=9)
    ax2.scatter(By_neg_b_perp[4:6], By_neg_swirlpower[4:6], color='blue', s=200,zorder=9) # 
    ax2.scatter(By_neg_b_perp[6], By_neg_swirlpower[6], color='blue', s=200, marker="v",zorder=9)
    ax2.scatter(By_neg_b_perp[7:], By_neg_swirlpower[7:], color='blue', s=200, marker=",",zorder=9)
    
    
    ax2.errorbar(b_perp_medians, medians_swirl, yerr=swirl_med_errs, xerr=med_b_perp_errs, fmt='.', color='lightgray',zorder=9)
    
    ax2.text(1.92,0.8,'b',style='italic',fontsize=40) #0.5, 0.75
    
    ax2.set_xlabel('Median B$_⊥$ (nT)',fontsize=22)
    ax2.set_ylabel('Median Swirl Region Power (GW)',fontsize=22)
    
    ax2.set_xlim(0.0,2.05)
    ax2.set_ylim(-2,82) #81 errors, #73
    ax2.tick_params(axis='x',which='minor',direction='in')
    ax2.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax2.minorticks_on()
    ax2.tick_params(which='major',direction='in',bottom=True, top=True, left=True, right=True, labelsize=20, length=10)
    ax2.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True, length=5)
    ax2.legend(framealpha=0.5,fontsize=21, loc ="upper left")
    
    #########
    
    ax3 = plt.subplot(4,2,3)
    
    ax3.plot(x_new, y_new_dusk,'--', color='black', linewidth=3, markersize=12, label=f'Dusk $R^2$: {adjr_d:.2g}', zorder=10)

    if error_plot == '20' or error_plot == '10' or error_plot == '50':
        ax3.errorbar(b_perp_medians_plus, dusk_medians, yerr=dusk_med_errs, xerr=med_b_perp_errs_plus, fmt='.', color='plum')
        
        ax3.errorbar(b_perp_medians_minus, dusk_medians, yerr=dusk_med_errs, xerr=med_b_perp_errs_minus, fmt='.', color='peachpuff')
        
        ax3.plot(x_new, y_new_dusk_p, '--',color='darkviolet', linewidth=3, markersize=12, label=f'+ {error}%: {adjr_dp:.2g}')
        ax3.plot(x_new, y_new_dusk_m, '--',color='sandybrown', linewidth=3, markersize=12, label=f'- {error}%: {adjr_dm:.2g}')
        
        # plus error
        ax3.scatter(b_perp_medians_plus[0:2], dusk_medians[0:2],s=100,color='plum',marker="^")
        ax3.scatter(b_perp_medians_plus[2], dusk_medians[2],s=100,color='plum')
        ax3.scatter(b_perp_medians_plus[3:5], dusk_medians[3:5],s=100,color='plum',marker="^")
        ax3.scatter(b_perp_medians_plus[5:9], dusk_medians[5:9],s=100,color='plum')
        ax3.scatter(b_perp_medians_plus[9], dusk_medians[9],s=100,color='plum',marker="v")
        ax3.scatter(b_perp_medians_plus[10], dusk_medians[10],s=100,color='plum')
        ax3.scatter(b_perp_medians_plus[11], dusk_medians[11],s=100,color='plum',marker="v")
        ax3.scatter(b_perp_medians_plus[12:18], dusk_medians[12:18],s=100,color='plum',marker=",")
        ax3.scatter(b_perp_medians_plus[18:21], dusk_medians[18:21],s=100,color='plum',marker="v")
        ax3.scatter(b_perp_medians_plus[21], dusk_medians[21],s=100,color='plum',marker="^")
        ax3.scatter(b_perp_medians_plus[22:], dusk_medians[22:],s=100,color='plum',marker=",")
        
        # minus error
        ax3.scatter(b_perp_medians_minus[0:2], dusk_medians[0:2],s=100,color='peachpuff',marker="^")
        ax3.scatter(b_perp_medians_minus[2], dusk_medians[2],s=100,color='peachpuff')
        ax3.scatter(b_perp_medians_minus[3:5], dusk_medians[3:5],s=100,color='peachpuff',marker="^")
        ax3.scatter(b_perp_medians_minus[5:9], dusk_medians[5:9],s=100,color='peachpuff')
        ax3.scatter(b_perp_medians_minus[9], dusk_medians[9],s=100,color='peachpuff',marker="v")
        ax3.scatter(b_perp_medians_minus[10], dusk_medians[10],s=100,color='peachpuff')
        ax3.scatter(b_perp_medians_minus[11], dusk_medians[11],s=100,color='peachpuff',marker="v")
        ax3.scatter(b_perp_medians_minus[12:18], dusk_medians[12:18],s=100,color='peachpuff',marker=",")
        ax3.scatter(b_perp_medians_minus[18:21], dusk_medians[18:21],s=100,color='peachpuff',marker="v")
        ax3.scatter(b_perp_medians_minus[21], dusk_medians[21],s=100,color='peachpuff',marker="^")
        ax3.scatter(b_perp_medians_minus[22:], dusk_medians[22:],s=100,color='peachpuff',marker=",")
    
    # main plot
    ax3.scatter(Bz_pos_b_perp[0:2], Bz_pos_duskpower[0:2], color='red', s=200,zorder=9)
    ax3.scatter(Bz_pos_b_perp[2:], Bz_pos_duskpower[2:], color='red', s=200, marker='v',zorder=9)

    ax3.scatter(Bz_neg_b_perp[0], Bz_neg_duskpower[0], color='orange', s=200,zorder=9) #
    ax3.scatter(Bz_neg_b_perp[1], Bz_neg_duskpower[1], color='orange', s=200,marker=',',zorder=9)

    ax3.scatter(By_pos_b_perp[0], By_pos_duskpower[0], color='green', s=200,zorder=9) # compression 1 2-5 incl visit 15
    ax3.scatter(By_pos_b_perp[1:4], By_pos_duskpower[1:4], color='green', s=200, marker=",",zorder=9)
    ax3.scatter(By_pos_b_perp[4:6], By_pos_duskpower[4:6], color='green', s=200, marker="v",zorder=9) # rarefraction
    ax3.scatter(By_pos_b_perp[6], By_pos_duskpower[6], color='green', s=200, marker="^",zorder=9)
    ax3.scatter(By_pos_b_perp[7:], By_pos_duskpower[7:], color='green', s=200, marker=",",zorder=9)

    ax3.scatter(By_neg_b_perp[0:4], By_neg_duskpower[0:4], color='blue', s=200, marker="^",zorder=9)
    ax3.scatter(By_neg_b_perp[4:6], By_neg_duskpower[4:6], color='blue', s=200,zorder=9) # 
    ax3.scatter(By_neg_b_perp[6], By_neg_duskpower[6], color='blue', s=200, marker="v",zorder=9)
    ax3.scatter(By_neg_b_perp[7:], By_neg_duskpower[7:], color='blue', s=200, marker=",",zorder=9)
    
    ax3.errorbar(b_perp_medians, dusk_medians, yerr=dusk_med_errs, xerr=med_b_perp_errs, fmt='.', color='lightgray', zorder=9)
    
    ax3.text(1.92,-7,'c',style='italic',fontsize=40)

    ax3.set_xlabel('Median B$_⊥$ (nT)',fontsize=22)
    ax3.set_ylabel('Median Dusk Active Region Power (GW)',fontsize=22)
 
    ax3.set_xlim(0.0,2.05)
    ax3.set_ylim(-25,508) #495 error, 420
    ax3.tick_params(axis='x',which='minor',direction='in')
    ax3.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax3.minorticks_on()
    ax3.tick_params(which='major',direction='in',bottom=True, top=True, left=True, right=True, labelsize=20, length=10)
    ax3.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True, length=5)
    ax3.legend(framealpha=0.5,fontsize=21,loc ="upper left")
    
    
    ########
    
    ax4 = plt.subplot(4,2,4)
    
    ax4.plot(x_new, y_new_noon, '--',color='black', linewidth=3, markersize=12, label=f'Noon $R^2$: {adjr_n:.2g}', zorder=10)

    if error_plot == '20' or error_plot == '10' or error_plot == '50':
        ax4.errorbar(b_perp_medians_plus, noon_medians, yerr=noon_med_errs, xerr=med_b_perp_errs_plus, fmt='.', color='plum')
        
        ax4.errorbar(b_perp_medians_minus, noon_medians, yerr=noon_med_errs, xerr=med_b_perp_errs_minus, fmt='.', color='peachpuff')
        
        ax4.plot(x_new, y_new_noon_p, '--',color='darkviolet', linewidth=3, markersize=12, label=f'+ {error}%: {adjr_np:.2g}')
        ax4.plot(x_new, y_new_noon_m, '--',color='sandybrown', linewidth=3, markersize=12, label=f'- {error}%: {adjr_nm:.2g}')
        
        # plus error
        ax4.scatter(b_perp_medians_plus[0:2], noon_medians[0:2],s=100,color='plum',marker="^")
        ax4.scatter(b_perp_medians_plus[2], noon_medians[2],s=100,color='plum')
        ax4.scatter(b_perp_medians_plus[3:5], noon_medians[3:5],s=100,color='plum',marker="^")
        ax4.scatter(b_perp_medians_plus[5:9], noon_medians[5:9],s=100,color='plum')
        ax4.scatter(b_perp_medians_plus[9], noon_medians[9],s=100,color='plum',marker="v")
        ax4.scatter(b_perp_medians_plus[10], noon_medians[10],s=100,color='plum')
        ax4.scatter(b_perp_medians_plus[11], noon_medians[11],s=100,color='plum',marker="v")
        ax4.scatter(b_perp_medians_plus[12:18], noon_medians[12:18],s=100,color='plum',marker=",")
        ax4.scatter(b_perp_medians_plus[18:21], noon_medians[18:21],s=100,color='plum',marker="v")
        ax4.scatter(b_perp_medians_plus[21], noon_medians[21],s=100,color='plum',marker="^")
        ax4.scatter(b_perp_medians_plus[22:], noon_medians[22:],s=100,color='plum',marker=",")
    
        # minus error
        ax4.scatter(b_perp_medians_minus[0:2], noon_medians[0:2],s=100,color='peachpuff',marker="^")
        ax4.scatter(b_perp_medians_minus[2], noon_medians[2],s=100,color='peachpuff')
        ax4.scatter(b_perp_medians_minus[3:5], noon_medians[3:5],s=100,color='peachpuff',marker="^")
        ax4.scatter(b_perp_medians_minus[5:9], noon_medians[5:9],s=100,color='peachpuff')
        ax4.scatter(b_perp_medians_minus[9], noon_medians[9],s=100,color='peachpuff',marker="v")
        ax4.scatter(b_perp_medians_minus[10], noon_medians[10],s=100,color='peachpuff')
        ax4.scatter(b_perp_medians_minus[11], noon_medians[11],s=100,color='peachpuff',marker="v")
        ax4.scatter(b_perp_medians_minus[12:18], noon_medians[12:18],s=100,color='peachpuff',marker=",")
        ax4.scatter(b_perp_medians_minus[18:21], noon_medians[18:21],s=100,color='peachpuff',marker="v")
        ax4.scatter(b_perp_medians_minus[21], noon_medians[21],s=100,color='peachpuff',marker="^")
        ax4.scatter(b_perp_medians_minus[22:], noon_medians[22:],s=100,color='peachpuff',marker=",")
    
    # main plot
    ax4.scatter(Bz_pos_b_perp[0:2], Bz_pos_noonpower[0:2], color='red', s=200,zorder=9)
    ax4.scatter(Bz_pos_b_perp[2:], Bz_pos_noonpower[2:], color='red', s=200, marker='v',zorder=9)

    ax4.scatter(Bz_neg_b_perp[0], Bz_neg_noonpower[0], color='orange', s=200,zorder=9) #
    ax4.scatter(Bz_neg_b_perp[1], Bz_neg_noonpower[1], color='orange', s=200,marker=',',zorder=9) 

    ax4.scatter(By_pos_b_perp[0], By_pos_noonpower[0], color='green', s=200,zorder=9) # compression 1 2-5 incl visit 15
    ax4.scatter(By_pos_b_perp[1:4], By_pos_noonpower[1:4], color='green', s=200, marker=",",zorder=9)
    ax4.scatter(By_pos_b_perp[4:6], By_pos_noonpower[4:6], color='green', s=200, marker="v",zorder=9) # rarefraction
    ax4.scatter(By_pos_b_perp[6], By_pos_noonpower[6], color='green', s=200, marker="^",zorder=9)
    ax4.scatter(By_pos_b_perp[7:], By_pos_noonpower[7:], color='green', s=200, marker=",",zorder=9)

    ax4.scatter(By_neg_b_perp[0:4], By_neg_noonpower[0:4], color='blue', s=200, marker="^",zorder=9)
    ax4.scatter(By_neg_b_perp[4:6], By_neg_noonpower[4:6], color='blue', s=200,zorder=9)
    ax4.scatter(By_neg_b_perp[6], By_neg_noonpower[6], color='blue', s=200, marker="v",zorder=9)
    ax4.scatter(By_neg_b_perp[7:], By_neg_noonpower[7:], color='blue', s=200, marker=",",zorder=9)

    
    ax4.errorbar(b_perp_medians, noon_medians, yerr=noon_med_errs, xerr=med_b_perp_errs, fmt='.', color='lightgray', zorder=9)
    
    ax4.text(1.92,-10,'d',style='italic',fontsize=40)

    ax4.set_ylabel('Median Noon Active Region Power (GW)',fontsize=22)
    ax4.set_xlabel('Median B$_⊥$ (nT)',fontsize=22)

    ax4.set_xlim(0.0,2.05)
    ax4.tick_params(axis='x',which='minor',direction='in')
    ax4.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax4.minorticks_on()
    ax4.tick_params(which='major',direction='in',bottom=True, top=True, left=True, right=True, labelsize=20, length=10)
    ax4.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True, length=5)
    ax4.legend(framealpha=0.5,fontsize=21,loc ="upper left")
    
    # save plot
    saveloc = (f'{root_saves}median_b_perp_vs_region_power_{error_plot}.jpg') 
    #saveloc = (f'/Users/hannah/OneDrive - Lancaster University/aurora/median_b_perp_vs_power_{error}.jpg') #remember to add f' at front when code {error} back on
    plt.savefig(saveloc,bbox_inches='tight',dpi=400)
    


if plotting == 'clock_median': 
    
     # plt.scatter(b_perp_medians,pressure_medians, s=20)
     # plt.xlabel('Median B$_⊥$ (nT)')
     # plt.ylabel('Median Pressure (nPa)')
     
     # # save plot
     # saveloc = (f'{root_saves}scatter_pressure_bperp.jpg') 
     # #saveloc = (f'{root_saves}median_clock_vs_region_power_extra3_{error}.jpg')
     # plt.savefig(saveloc,bbox_inches='tight',dpi=400)
     # plt.tick_params(which='major', direction='in', bottom=True, top=True, left=True, right=True)
     
     # plt.scatter(btot_medians,pressure_medians, s=20)
     # plt.xlabel('Median B$_{total}$ (nT)')
     # plt.ylabel('Median Pressure (nPa)')
     # plt.tick_params(which='major', direction='in', bottom=True, top=True, left=True, right=True)
     
     # # save plot
     # saveloc = (f'{root_saves}scatter_pressure_btot.jpg') 
     # #saveloc = (f'{root_saves}median_clock_vs_region_power_extra3_{error}.jpg')
     # plt.savefig(saveloc,bbox_inches='tight',dpi=400)
     
     # plt.scatter(btot_medians,pressure_medians, s=10, label='B$_{total}$')
     # plt.scatter(b_perp_medians,pressure_medians, s=10, label='B$_⊥$')
     # plt.xlabel('Median Magnetic Field Strength (nT)')
     # plt.ylabel('Median Pressure (nPa)')
     # plt.tick_params(which='major', direction='in', bottom=True, top=True, left=True, right=True)
     # plt.legend(loc='upper left')
     
     # # save plot
     # saveloc = (f'{root_saves}scatter_pressure_magboth.jpg') 
     # plt.savefig(saveloc,bbox_inches='tight',dpi=400)
    
     '''
     fitting
     '''
    
     # arrange cml more evenly for fit plotting
     x_new = np.linspace(min(clock_medians), max(clock_medians), len(clock_medians))
    
     # # arrange cml more evenly for fit plotting
     x_new = np.linspace(min(clock_medians), max(clock_medians), len(clock_medians))
     x_new = np.linspace(-180, 180, 360)

     a = np.polyfit(clock_medians, medians_swirl, 3)
     aa = np.poly1d(a)
     
     b = np.polyfit(clock_medians, dusk_medians, 3)
     bb = np.poly1d(b)
     
     c = np.polyfit(clock_medians, noon_medians, 3)
     cc = np.poly1d(c)
     
     y_new_swirl = aa(x_new)
     y_new_noon = cc(x_new)
     y_new_dusk = bb(x_new)

     # def constant_model(x, offset):
     #     # print(x)
     #     # print(offset)
     #     return np.ones_like(x) * offset
     
    #clock_rad = np.radians(clock_medians)
     if key == 'whole':
         def curve_making(x,amp1,amp2,offset):#,offset):
             if np.max(np.abs(x)) > 2 * np.pi:  # likely in degrees
                 x = np.radians(x)
             if key == 'whole' or key == 'whole_nooff':
                 return (((np.sin(x)**2)*amp1) + ((np.cos(x/2)**2)*amp2)) + offset # fit for full system of By and Bz (seperate ampltitudes)
             
     else:
         def curve_making(x,amp1,offset):#,offset):
             if np.max(np.abs(x)) > 2 * np.pi:  # likely in degrees
                 x = np.radians(x)
             if key == 'by' or key == 'by_nooff':
                 return ((np.sin(x)**2) * amp1) +offset# fit for By dependence
             #return ((np.sin((x+np.deg2rad(90))/2)**2) * amplitude) + offset # fit for By dependence
             if key == 'earth' or key == 'earth_nooff':
                 return (np.cos(x/2)**2 * amp1) +offset # fit for Earth-like reconnection
             if key == 'tang' or key == 'tang_nooff':
                 return (np.sin(x/2)**2 * amp1) +offset # for for testing for tangled field lines on nightside
             if key == 'test':
                 return ((np.sin((x+np.deg2rad(90))/2)**2) * amp1) + offset
         
     # def curve_making_extra(x,amp1):
     #     if np.max(np.abs(x)) > 2 * np.pi:  # likely in degrees
     #        x = np.radians(x)
     #     if key == 'by' or key == 'by_nooff':
     #        return ((np.sin(x)**2) * amp1) #+ offset # fit for By dependence
     #     #return ((np.sin((x+np.deg2rad(90))/2)**2) * amplitude) + offset # fit for By dependence
     #     if key == 'earth' or key == 'earth_nooff':
     #        return (np.cos(x/2)**2 * amp1) #+offset # fit for Earth-like reconnection
     #     if key == 'tang' or key == 'tang_nooff':
     #        return (np.sin(x/2)**2 * amp1) #+offset # for for testing for tangled field lines on nightside
     #     if key == 'whole' or key == 'whole_nooff':
     #        return (((np.sin(x)**2)*amp1) + ((np.cos(x/2)**2)*amp2)) #+ offset # fit for full system of By and Bz (seperate ampltitudes)
     
     if key == 'whole_nooff' or key == 'by' or key == 'tang' or key == 'earth' or key == 'test':
         popt_swirl, pcov_swirl = curve_fit(curve_making, clock_medians, medians_swirl)#, p0=p0)#, b#ounds = ((-np.pi/2, -np.pi/2, -np.pi/2, -np.pi/2), (np.pi/2, np.pi/2, np.pi/2, np.pi/2)))#, p0=(0.1, 1, 2, 0.5, 1))#bounds = ((0.1, 0.1, 0.1, 0.1, 0.1), (2, 2, 2, 2, 2)))#p0=(), bounds = ()
         a_opt_swirl,b_opt_swirl = popt_swirl #, 
        
         y_new_swirl = curve_making(x_new, a_opt_swirl, b_opt_swirl)
         y_regress_swirl = curve_making(clock_medians, a_opt_swirl, b_opt_swirl)
    
    
         popt_dusk, pcov_dusk = curve_fit(curve_making, clock_medians, dusk_medians)#,p0=p0)#, p0=p0)#, b#ounds = ((-np.pi/2, -np.pi/2, -np.pi/2, -np.pi/2), (np.pi/2, np.pi/2, np.pi/2, np.pi/2)))#, p0=(0.1, 1, 2, 0.5, 1))#bounds = ((0.1, 0.1, 0.1, 0.1, 0.1), (2, 2, 2, 2, 2)))#p0=(), bounds = ()
         a_opt_dusk, b_opt_dusk = popt_dusk #
         
          
         y_new_dusk = curve_making(x_new, a_opt_dusk, b_opt_dusk)#
         y_regress_dusk = curve_making(clock_medians, a_opt_dusk, b_opt_dusk)
        
        
         popt_noon, pcov_noon = curve_fit(curve_making, clock_medians, noon_medians)#, p0=p0)#, b#ounds = ((-np.pi/2, -np.pi/2, -np.pi/2, -np.pi/2), (np.pi/2, np.pi/2, np.pi/2, np.pi/2)))#, p0=(0.1, 1, 2, 0.5, 1))#bounds = ((0.1, 0.1, 0.1, 0.1, 0.1), (2, 2, 2, 2, 2)))#p0=(), bounds = ()
         a_opt_noon,b_opt_noon = popt_noon #, 
        
         y_new_noon = curve_making(x_new, a_opt_noon, b_opt_noon)
         y_regress_noon = curve_making(clock_medians, a_opt_noon, b_opt_noon)#
         
        
         # +%
         popt_noon_p, pcov_noon_p = curve_fit(curve_making, clock_medians_p, (noon_medians))#, p0=p0)#, b#ounds = ((-np.pi/2, -np.pi/2, -np.pi/2, -np.pi/2), (np.pi/2, np.pi/2, np.pi/2, np.pi/2)))#, p0=(0.1, 1, 2, 0.5, 1))#bounds = ((0.1, 0.1, 0.1, 0.1, 0.1), (2, 2, 2, 2, 2)))#p0=(), bounds = ()
         a_opt_noon_p, b_opt_noon_p = popt_noon_p #
        
         y_new_noon_p = curve_making(x_new, a_opt_noon_p, b_opt_noon_p)
         y_regress_noon_p = curve_making(clock_medians_p, a_opt_noon_p, b_opt_noon_p)
        
         popt_swirl_p, pcov_swirl_p = curve_fit(curve_making, clock_medians_p, (medians_swirl))#, p0=p0)#, b#ounds = ((-np.pi/2, -np.pi/2, -np.pi/2, -np.pi/2), (np.pi/2, np.pi/2, np.pi/2, np.pi/2)))#, p0=(0.1, 1, 2, 0.5, 1))#bounds = ((0.1, 0.1, 0.1, 0.1, 0.1), (2, 2, 2, 2, 2)))#p0=(), bounds = ()
         a_opt_swirl_p, b_opt_swirl_p = popt_swirl_p #
        
         y_new_swirl_p = curve_making(x_new, a_opt_swirl_p, b_opt_swirl_p)
         y_regress_swirl_p = curve_making(clock_medians_p, a_opt_swirl_p, b_opt_swirl_p)
         
         popt_dusk_p, pcov_dusk_p = curve_fit(curve_making, clock_medians_p, (dusk_medians))#, p0=p0)#, b#ounds = ((-np.pi/2, -np.pi/2, -np.pi/2, -np.pi/2), (np.pi/2, np.pi/2, np.pi/2, np.pi/2)))#, p0=(0.1, 1, 2, 0.5, 1))#bounds = ((0.1, 0.1, 0.1, 0.1, 0.1), (2, 2, 2, 2, 2)))#p0=(), bounds = ()
         a_opt_dusk_p, b_opt_dusk_p = popt_dusk_p #
         
         y_new_dusk_p = curve_making(x_new, a_opt_dusk_p, b_opt_dusk_p)
         y_regress_dusk_p = curve_making(clock_medians_p, a_opt_dusk_p, b_opt_dusk_p)
        
        
         # -%
         popt_noon_m, pcov_noon_m = curve_fit(curve_making, clock_medians_m, (noon_medians))#, p0=p0)#, b#ounds = ((-np.pi/2, -np.pi/2, -np.pi/2, -np.pi/2), (np.pi/2, np.pi/2, np.pi/2, np.pi/2)))#, p0=(0.1, 1, 2, 0.5, 1))#bounds = ((0.1, 0.1, 0.1, 0.1, 0.1), (2, 2, 2, 2, 2)))#p0=(), bounds = ()
         a_opt_noon_m, b_opt_noon_m = popt_noon_m # 
        
         y_new_noon_m = curve_making(x_new, a_opt_noon_m, b_opt_noon_m)
         y_regress_noon_m = curve_making(clock_medians_m, a_opt_noon_m, b_opt_noon_m)
        
         popt_swirl_m, pcov_swirl_m = curve_fit(curve_making, clock_medians_m, (medians_swirl))#, p0=p0)#, b#ounds = ((-np.pi/2, -np.pi/2, -np.pi/2, -np.pi/2), (np.pi/2, np.pi/2, np.pi/2, np.pi/2)))#, p0=(0.1, 1, 2, 0.5, 1))#bounds = ((0.1, 0.1, 0.1, 0.1, 0.1), (2, 2, 2, 2, 2)))#p0=(), bounds = ()
         a_opt_swirl_m, b_opt_swirl_m = popt_swirl_m #
        
         y_new_swirl_m = curve_making(x_new, a_opt_swirl_m, b_opt_swirl_m)
         y_regress_swirl_m = curve_making(clock_medians_m, a_opt_swirl_m, b_opt_swirl_m)
         
         popt_dusk_m, pcov_dusk_m = curve_fit(curve_making, clock_medians_m, (dusk_medians))#, p0=p0)#, b#ounds = ((-np.pi/2, -np.pi/2, -np.pi/2, -np.pi/2), (np.pi/2, np.pi/2, np.pi/2, np.pi/2)))#, p0=(0.1, 1, 2, 0.5, 1))#bounds = ((0.1, 0.1, 0.1, 0.1, 0.1), (2, 2, 2, 2, 2)))#p0=(), bounds = ()
         a_opt_dusk_m, b_opt_dusk_m = popt_dusk_m #
        
         
         y_new_dusk_m = curve_making(x_new, a_opt_dusk_m, b_opt_dusk_m)
         y_regress_dusk_m = curve_making(clock_medians_m, a_opt_dusk_m, b_opt_dusk_m)
         
     elif key == 'whole':
         popt_swirl, pcov_swirl = curve_fit(curve_making, clock_medians, medians_swirl)#, p0=p0)#, b#ounds = ((-np.pi/2, -np.pi/2, -np.pi/2, -np.pi/2), (np.pi/2, np.pi/2, np.pi/2, np.pi/2)))#, p0=(0.1, 1, 2, 0.5, 1))#bounds = ((0.1, 0.1, 0.1, 0.1, 0.1), (2, 2, 2, 2, 2)))#p0=(), bounds = ()
         a_opt_swirl,b_opt_swirl, c_opt_swirl = popt_swirl #, 
        
         y_new_swirl = curve_making(x_new, a_opt_swirl, b_opt_swirl, c_opt_swirl)
         y_regress_swirl = curve_making(clock_medians, a_opt_swirl, b_opt_swirl, c_opt_swirl)
    
    
         popt_dusk, pcov_dusk = curve_fit(curve_making, clock_medians, dusk_medians)#,p0=p0)#, p0=p0)#, b#ounds = ((-np.pi/2, -np.pi/2, -np.pi/2, -np.pi/2), (np.pi/2, np.pi/2, np.pi/2, np.pi/2)))#, p0=(0.1, 1, 2, 0.5, 1))#bounds = ((0.1, 0.1, 0.1, 0.1, 0.1), (2, 2, 2, 2, 2)))#p0=(), bounds = ()
         a_opt_dusk, b_opt_dusk, c_opt_dusk = popt_dusk #
          
         y_new_dusk = curve_making(x_new, a_opt_dusk, b_opt_dusk, c_opt_dusk)#
         y_regress_dusk = curve_making(clock_medians, a_opt_dusk, b_opt_dusk, c_opt_dusk)
        
        
         popt_noon, pcov_noon = curve_fit(curve_making, clock_medians, noon_medians)#, p0=p0)#, b#ounds = ((-np.pi/2, -np.pi/2, -np.pi/2, -np.pi/2), (np.pi/2, np.pi/2, np.pi/2, np.pi/2)))#, p0=(0.1, 1, 2, 0.5, 1))#bounds = ((0.1, 0.1, 0.1, 0.1, 0.1), (2, 2, 2, 2, 2)))#p0=(), bounds = ()
         a_opt_noon,b_opt_noon, c_opt_noon = popt_noon #, 
        
         y_new_noon = curve_making(x_new, a_opt_noon, b_opt_noon, c_opt_noon)#
         y_regress_noon = curve_making(clock_medians, a_opt_noon, b_opt_noon, c_opt_noon)
         
        
         # +%
         popt_noon_p, pcov_noon_p = curve_fit(curve_making, clock_medians_p, (noon_medians))#, p0=p0)#, b#ounds = ((-np.pi/2, -np.pi/2, -np.pi/2, -np.pi/2), (np.pi/2, np.pi/2, np.pi/2, np.pi/2)))#, p0=(0.1, 1, 2, 0.5, 1))#bounds = ((0.1, 0.1, 0.1, 0.1, 0.1), (2, 2, 2, 2, 2)))#p0=(), bounds = ()
         a_opt_noon_p, b_opt_noon_p, c_opt_noon_p = popt_noon_p #
        
         y_new_noon_p = curve_making(x_new, a_opt_noon_p, b_opt_noon_p, c_opt_noon_p)
         y_regress_noon_p = curve_making(clock_medians_p, a_opt_noon_p, b_opt_noon_p, c_opt_noon_p)
        
         popt_swirl_p, pcov_swirl_p = curve_fit(curve_making, clock_medians_p, (medians_swirl))#, p0=p0)#, b#ounds = ((-np.pi/2, -np.pi/2, -np.pi/2, -np.pi/2), (np.pi/2, np.pi/2, np.pi/2, np.pi/2)))#, p0=(0.1, 1, 2, 0.5, 1))#bounds = ((0.1, 0.1, 0.1, 0.1, 0.1), (2, 2, 2, 2, 2)))#p0=(), bounds = ()
         a_opt_swirl_p, b_opt_swirl_p, c_opt_swirl_p  = popt_swirl_p #
        
         y_new_swirl_p = curve_making(x_new, a_opt_swirl_p, b_opt_swirl_p, c_opt_swirl_p)
         y_regress_swirl_p = curve_making(clock_medians_p, a_opt_swirl_p, b_opt_swirl_p, c_opt_swirl_p)
         
         popt_dusk_p, pcov_dusk_p = curve_fit(curve_making, clock_medians_p, (dusk_medians))#, p0=p0)#, b#ounds = ((-np.pi/2, -np.pi/2, -np.pi/2, -np.pi/2), (np.pi/2, np.pi/2, np.pi/2, np.pi/2)))#, p0=(0.1, 1, 2, 0.5, 1))#bounds = ((0.1, 0.1, 0.1, 0.1, 0.1), (2, 2, 2, 2, 2)))#p0=(), bounds = ()
         a_opt_dusk_p, b_opt_dusk_p, c_opt_dusk_p = popt_dusk_p #
         
         y_new_dusk_p = curve_making(x_new, a_opt_dusk_p, b_opt_dusk_p, c_opt_dusk_p)
         y_regress_dusk_p = curve_making(clock_medians_p, a_opt_dusk_p, b_opt_dusk_p, c_opt_dusk_p)
        
        
         # -%
         popt_noon_m, pcov_noon_m = curve_fit(curve_making, clock_medians_m, (noon_medians))#, p0=p0)#, b#ounds = ((-np.pi/2, -np.pi/2, -np.pi/2, -np.pi/2), (np.pi/2, np.pi/2, np.pi/2, np.pi/2)))#, p0=(0.1, 1, 2, 0.5, 1))#bounds = ((0.1, 0.1, 0.1, 0.1, 0.1), (2, 2, 2, 2, 2)))#p0=(), bounds = ()
         a_opt_noon_m, b_opt_noon_m, c_opt_noon_m = popt_noon_m # 
        
         y_new_noon_m = curve_making(x_new, a_opt_noon_m, b_opt_noon_m, c_opt_noon_m)
         y_regress_noon_m = curve_making(clock_medians_m, a_opt_noon_m, b_opt_noon_m, c_opt_noon_m)
        
         popt_swirl_m, pcov_swirl_m = curve_fit(curve_making, clock_medians_m, (medians_swirl))#, p0=p0)#, b#ounds = ((-np.pi/2, -np.pi/2, -np.pi/2, -np.pi/2), (np.pi/2, np.pi/2, np.pi/2, np.pi/2)))#, p0=(0.1, 1, 2, 0.5, 1))#bounds = ((0.1, 0.1, 0.1, 0.1, 0.1), (2, 2, 2, 2, 2)))#p0=(), bounds = ()
         a_opt_swirl_m, b_opt_swirl_m , c_opt_swirl_m= popt_swirl_m #
        
         y_new_swirl_m = curve_making(x_new, a_opt_swirl_m, b_opt_swirl_m, c_opt_swirl_m)
         y_regress_swirl_m = curve_making(clock_medians_m, a_opt_swirl_m, b_opt_swirl_m, c_opt_swirl_m)
         
         popt_dusk_m, pcov_dusk_m = curve_fit(curve_making, clock_medians_m, (dusk_medians))#, p0=p0)#, b#ounds = ((-np.pi/2, -np.pi/2, -np.pi/2, -np.pi/2), (np.pi/2, np.pi/2, np.pi/2, np.pi/2)))#, p0=(0.1, 1, 2, 0.5, 1))#bounds = ((0.1, 0.1, 0.1, 0.1, 0.1), (2, 2, 2, 2, 2)))#p0=(), bounds = ()
         a_opt_dusk_m, b_opt_dusk_m, c_opt_dusk_m = popt_dusk_m #
        
         y_new_dusk_m = curve_making(x_new, a_opt_dusk_m, b_opt_dusk_m, c_opt_dusk_m)
         y_regress_dusk_m = curve_making(clock_medians_m, a_opt_dusk_m, b_opt_dusk_m, c_opt_dusk_m)
         
         
     else:
        popt_swirl, pcov_swirl = curve_fit(curve_making, clock_medians, medians_swirl)#, p0=p0)#, b#ounds = ((-np.pi/2, -np.pi/2, -np.pi/2, -np.pi/2), (np.pi/2, np.pi/2, np.pi/2, np.pi/2)))#, p0=(0.1, 1, 2, 0.5, 1))#bounds = ((0.1, 0.1, 0.1, 0.1, 0.1), (2, 2, 2, 2, 2)))#p0=(), bounds = ()
        a_opt_swirl = popt_swirl #, 
       
        y_new_swirl = curve_making(x_new, a_opt_swirl)
        y_regress_swirl = curve_making(clock_medians, a_opt_swirl)
   
   
        popt_dusk, pcov_dusk = curve_fit(curve_making, clock_medians, dusk_medians)#,p0=p0)#, p0=p0)#, b#ounds = ((-np.pi/2, -np.pi/2, -np.pi/2, -np.pi/2), (np.pi/2, np.pi/2, np.pi/2, np.pi/2)))#, p0=(0.1, 1, 2, 0.5, 1))#bounds = ((0.1, 0.1, 0.1, 0.1, 0.1), (2, 2, 2, 2, 2)))#p0=(), bounds = ()
        a_opt_dusk = popt_dusk #
         
        y_new_dusk = curve_making(x_new, a_opt_dusk)#
        y_regress_dusk = curve_making(clock_medians, a_opt_dusk)#
       
       
        popt_noon, pcov_noon = curve_fit(curve_making, clock_medians, noon_medians)#, p0=p0)#, b#ounds = ((-np.pi/2, -np.pi/2, -np.pi/2, -np.pi/2), (np.pi/2, np.pi/2, np.pi/2, np.pi/2)))#, p0=(0.1, 1, 2, 0.5, 1))#bounds = ((0.1, 0.1, 0.1, 0.1, 0.1), (2, 2, 2, 2, 2)))#p0=(), bounds = ()
        a_opt_noon = popt_noon #, 
       
        y_new_noon = curve_making(x_new, a_opt_noon)#
        y_regress_noon = curve_making(clock_medians, a_opt_noon)#
        
       
        # +%
        popt_noon_p, pcov_noon_p = curve_fit(curve_making, clock_medians_p, (noon_medians))#, p0=p0)#, b#ounds = ((-np.pi/2, -np.pi/2, -np.pi/2, -np.pi/2), (np.pi/2, np.pi/2, np.pi/2, np.pi/2)))#, p0=(0.1, 1, 2, 0.5, 1))#bounds = ((0.1, 0.1, 0.1, 0.1, 0.1), (2, 2, 2, 2, 2)))#p0=(), bounds = ()
        a_opt_noon_p = popt_noon_p #
       
        y_new_noon_p = curve_making(x_new, a_opt_noon_p)
        y_regress_noon_p = curve_making(clock_medians_p, a_opt_noon_p)
       
        popt_swirl_p, pcov_swirl_p = curve_fit(curve_making, clock_medians_p, (medians_swirl))#, p0=p0)#, b#ounds = ((-np.pi/2, -np.pi/2, -np.pi/2, -np.pi/2), (np.pi/2, np.pi/2, np.pi/2, np.pi/2)))#, p0=(0.1, 1, 2, 0.5, 1))#bounds = ((0.1, 0.1, 0.1, 0.1, 0.1), (2, 2, 2, 2, 2)))#p0=(), bounds = ()
        a_opt_swirl_p = popt_swirl_p #
       
        y_new_swirl_p = curve_making(x_new, a_opt_swirl_p)
        y_regress_swirl_p = curve_making(clock_medians_p, a_opt_swirl_p)
        
        popt_dusk_p, pcov_dusk_p = curve_fit(curve_making, clock_medians_p, (dusk_medians))#, p0=p0)#, b#ounds = ((-np.pi/2, -np.pi/2, -np.pi/2, -np.pi/2), (np.pi/2, np.pi/2, np.pi/2, np.pi/2)))#, p0=(0.1, 1, 2, 0.5, 1))#bounds = ((0.1, 0.1, 0.1, 0.1, 0.1), (2, 2, 2, 2, 2)))#p0=(), bounds = ()
        a_opt_dusk_p = popt_dusk_p #
        
        y_new_dusk_p = curve_making(x_new, a_opt_dusk_p)
        y_regress_dusk_p = curve_making(clock_medians_p, a_opt_dusk_p)
        
        # -%
        popt_noon_m, pcov_noon_m = curve_fit(curve_making, clock_medians_m, (noon_medians))#, p0=p0)#, b#ounds = ((-np.pi/2, -np.pi/2, -np.pi/2, -np.pi/2), (np.pi/2, np.pi/2, np.pi/2, np.pi/2)))#, p0=(0.1, 1, 2, 0.5, 1))#bounds = ((0.1, 0.1, 0.1, 0.1, 0.1), (2, 2, 2, 2, 2)))#p0=(), bounds = ()
        a_opt_noon_m = popt_noon_m # 
       
        y_new_noon_m = curve_making(x_new, a_opt_noon_m)
        y_regress_noon_m = curve_making(clock_medians_m, a_opt_noon_m)
       
        popt_swirl_m, pcov_swirl_m = curve_fit(curve_making, clock_medians_m, (medians_swirl))#, p0=p0)#, b#ounds = ((-np.pi/2, -np.pi/2, -np.pi/2, -np.pi/2), (np.pi/2, np.pi/2, np.pi/2, np.pi/2)))#, p0=(0.1, 1, 2, 0.5, 1))#bounds = ((0.1, 0.1, 0.1, 0.1, 0.1), (2, 2, 2, 2, 2)))#p0=(), bounds = ()
        a_opt_swirl_m = popt_swirl_m #
       
        y_new_swirl_m = curve_making(x_new, a_opt_swirl_m)
        y_regress_swirl_m = curve_making(clock_medians_m, a_opt_swirl_m)
        
        popt_dusk_m, pcov_dusk_m = curve_fit(curve_making, clock_medians_m, (dusk_medians))#, p0=p0)#, b#ounds = ((-np.pi/2, -np.pi/2, -np.pi/2, -np.pi/2), (np.pi/2, np.pi/2, np.pi/2, np.pi/2)))#, p0=(0.1, 1, 2, 0.5, 1))#bounds = ((0.1, 0.1, 0.1, 0.1, 0.1), (2, 2, 2, 2, 2)))#p0=(), bounds = ()
        a_opt_dusk_m = popt_dusk_m #
       
        y_new_dusk_m = curve_making(x_new, a_opt_dusk_m)
        y_regress_dusk_m = curve_making(clock_medians_m, a_opt_dusk_m)
    
     clock_rad = np.radians(clock_medians)
    
     print(f"Testing '{key}' Coupling Function")
     
     if key == 'whole':
         # whole system + offset
         swirl_amp1 = str(np.round(a_opt_swirl,2))
         swirl_amp2 = str(np.round(b_opt_swirl,2))
         swirl_off = str(np.round(c_opt_swirl,2))
         
         swirl_fit = (
             f"{swirl_amp1}$\sin(\\theta)^2$ + "
             f"{swirl_amp2}$\cos(\\frac{{\\theta}}{{2}})^2$ + "
             f"{swirl_off}"
             )
         
     elif key == 'by':
        # offset, one amplitude
        swirl_amp = str(np.round(a_opt_swirl,2))
        swirl_off = str(np.round(b_opt_swirl,2))

        #  By dependence
        #swirl_fit = (swirl_amp+'$\sin(\\theta)^2$+'+swirl_off+)

     elif key == 'earth':
       # offset, one amplitude
        swirl_amp = str(np.round(a_opt_swirl,2))
        swirl_off = str(np.round(b_opt_swirl,2))
        swirl_fit_err = np.sqrt(np.diag(pcov_swirl))#pcov_swirl#
        swirl_fit_err1 = str(np.round(swirl_fit_err[0], 2))
        swirl_fit_err2 = str(np.round(swirl_fit_err[1], 2))
        
        # Earth-like
        swirl_fit = ('('+swirl_amp+'±' +swirl_fit_err1+')$\cos(\\frac{\\theta}{2})^2$+('+swirl_off+'±'+swirl_fit_err2+')')
        
     elif key == 'tang':
        # offset, one amplitude
        swirl_amp = str(np.round(a_opt_swirl,2))
        swirl_off = str(np.round(b_opt_swirl,2))
        swirl_fit_err = np.sqrt(np.diag(pcov_swirl))#pcov_swirl#
        swirl_fit_err1 = str(np.round(swirl_fit_err[0], 2))
        swirl_fit_err2 = str(np.round(swirl_fit_err[1], 2))
        
        # Tangled Field
        swirl_fit = ('('+swirl_amp+'±' +swirl_fit_err1+')$\sin(\\frac{\\theta}{2})^2$+('+swirl_off+'±'+swirl_fit_err2+')')
     
     print('Swirl Region Fit')
     ss_res = np.sum((medians_swirl - y_regress_swirl)**2)
     ss_tot = np.sum((medians_swirl - np.mean(medians_swirl))**2)
     r_squared_s = 1 - (ss_res / ss_tot)
     print(f"R squared: {r_squared_s:.4f}")
     
     adj_s = 1-r_squared_s
     adjr_s = 1 - ((adj_s*23)/22)
     print(f"Adjusted R squared: {adjr_s:.4f}")
     '''
     pred_s = np.cos(clock_rad/2)**2#np.cos(2*clock_rad)#np.sin(clock_rad)**2
     
     predm_s = np.vstack([np.ones_like(pred_s), pred_s]).T
     # Linear regression: solve for coefficients [intercept, amplitude]
     coeffs, residuals, rank, s = lstsq(predm_s, medians_swirl, rcond=None)
     intercept, amplitude = coeffs
     
     # intercept_new = intercept + (amplitude/2)
     # amplitude_new = -(amplitude/2)

     print(f"Fitted intercept (offset): {intercept}")
     print(f"Fitted amplitude: {amplitude}")
     
     #coeffs_new = (intercept_new, amplitude_new)
     
     print(f"Original A: {a_opt_swirl}")
     print(f"Original c: {b_opt_swirl}")
     
     # A = -2 * amplitude
     # c = intercept + amplitude

     # print(f"Recovered original A: {A}")
     # print(f"Recovered original c: {c}")

     # Predicted values
     y_pred = predm_s @ coeffs

     # Calculate R^2
     ss_res2 = np.sum((medians_swirl - y_pred)**2)
     ss_tot2 = np.sum((medians_swirl - np.mean(medians_swirl))**2)
     r_squared2 = 1 - ss_res2 / ss_tot2
     print(f"R squared (linear): {r_squared2:.4f}")
     '''
     
     # # whole
     # x1 = np.cos(clock_rad)**2
     # x2 = np.cos(clock_rad)
     # X = np.vstack([np.ones_like(clock_rad), x1, x2]).T  # design matrix with intercept
     
     # cos_theta = np.cos(clock_rad)
     
     # coeffs, residuals, rank, s = lstsq(X, medians_swirl, rcond=None)
     # beta0, beta1, beta2 = coeffs
     
     # A1 = -beta1
     # A2 = 2 * beta2
     # c = beta0 - A1 - (A2 / 2)
     
     # print(f"Original A1: {a_opt_swirl}")
     # print(f"Original A2: {b_opt_swirl}")
     # print(f"Original c: {c_opt_swirl}")

     # y_pred = X @ coeffs
     
     # ss_res2 = np.sum((medians_swirl - y_pred) ** 2)
     # ss_tot2 = np.sum((medians_swirl - np.mean(medians_swirl)) ** 2)
     # r_squared2 = 1 - ss_res2 / ss_tot2

     # print(f"R^2 (linear) = {r_squared2:.4f}")
     # print(f"A1 = {A1:.4f}, A2 = {A2:.4f}, c = {c:.4f}")
     
     # z_s = beta1 * x1 + beta2 * x2
     # z_s2 = x1 + x2
     
     # import matplotlib.colors as mcolors

     
     # plt.figure(figsize=(10,6))
     # norm = mcolors.Normalize(vmin=-180, vmax=180)
     # #norm = mcolors.PowerNorm(gamma=0.4, vmin=-180, vmax=180)
     # sc = plt.scatter(z_s, medians_swirl, c=clock_medians, cmap='plasma', norm=norm, s=30)
     # plt.plot(np.sort(z_s2), y_pred[np.argsort(z_s2)], color='black')
     # #plt.plot(np.sort(np.cos(clock_rad)), medians_swirl)
     # for i in range(len(z_s)):
     #     plt.text(z_s[i], medians_swirl[i], f"{np.degrees(clock_rad[i]):.0f}", fontsize=14, ha='center')

     # #plt.plot(z_s, y_pred, color='black')

     # cbar = plt.colorbar(sc)
     # cbar.set_label('Clock Angle ($^o$)')

     # plt.xlabel('Linear Predictor: $A_1\'$cos²(θ) + $A_2\'$cos(θ)') # can always change A1 and A2
     # plt.ylabel('Swirl Power (GW)')
     # #saveloc = (f'{root_saves}linear_example.jpg') 
     # #plt.savefig(saveloc,bbox_inches='tight',dpi=400)
     
     # plt.figure(figsize=(10,6))
     # sc = plt.scatter(cos_theta, medians_swirl, c=clock_medians, cmap='plasma', s=30)
     # cbar = plt.colorbar(sc)
     # cbar.set_label('Clock Angle (°)')

     # # Fit curve: sort x so line is clean
     # sort_idx = np.argsort(cos_theta)
     # plt.plot(cos_theta[sort_idx], y_pred[sort_idx], color='black', label='Fit')

     # plt.xlabel('Cos(θ)')
     # plt.ylabel('Swirl Power (GW)')
     # saveloc = (f'{root_saves}linear_example_poly.jpg') 
     # plt.savefig(saveloc,bbox_inches='tight',dpi=400)
     
     
     

     rmse_s = np.sqrt(ss_res/22)

     print(f"RMSE: {rmse_s:.4f}")
     # res_s = np.sum(medians_swirl - y_regress_swirl)
     # print(f"Residual Sum: {res_s}")
     

     print('----------------------')
     print('Swirl Region +20% Fit')
     ssp_res = np.sum((medians_swirl - y_regress_swirl_p)**2)
     ssp_tot = np.sum((medians_swirl - np.mean(medians_swirl))**2)
     r_squared_sp = 1 - (ssp_res / ssp_tot)
     print(f"R squared: {r_squared_sp:.4f}")
     
     adj_sp = 1-r_squared_sp
     adjr_sp = 1 - ((adj_sp*23)/22)
     print(f"Adjusted R squared: {adjr_sp:.4f}")
     
     rmse_sp = np.sqrt(ssp_res/22)
     sd_sp = standard_deviation(medians_swirl, np.mean(medians_swirl))
     
     print(f"Mean: {np.mean(medians_swirl):.4f}")
     print(f"Standard Deviation {sd_sp:.4f}")
     print(f"RMSE: {rmse_sp:.4f}")

     
     ssm_res = np.sum((medians_swirl - y_regress_swirl_m)**2)
     ssm_tot = np.sum((medians_swirl - np.mean(medians_swirl))**2)
     r_squared_sm = 1 - (ssm_res / ssm_tot)
     print('----------------------')
     print('Swirl Region -20% Fit')
     print(f"R squared: {r_squared_sm:.4f}")
     
     adj_sm = 1-r_squared_sm
     adjr_sm = 1 - ((adj_sm*23)/22)
     print(f"Adjusted R squared: {adjr_sm:.4f}")
     
     rmse_sm = np.sqrt(ssm_res/22)
     sd_sm = standard_deviation(medians_swirl, np.mean(medians_swirl))
     
     print(f"Mean: {np.mean(medians_swirl):.4f}")
     print(f"Standard Deviation {sd_sm:.4f}")
     print(f"RMSE: {rmse_sm:.4f}")
     
    # ------------------
     
     sn_res = np.sum((noon_medians - y_regress_noon)**2)
     sn_tot = np.sum((noon_medians - np.mean(noon_medians))**2)
     r_squared_n = 1 - (sn_res / sn_tot)
     print('----------------------')
     print('Noon Region Fit')
     print(f"R squared: {r_squared_n:.4f}")
     
     adj_n = 1-r_squared_n
     adjr_n = 1 - ((adj_n*23)/22)
     print(f"Adjusted R squared: {adjr_n:.4f}")
     '''
     # whole
     x1_noon = np.cos(clock_rad)**2
     x2_noon = np.cos(clock_rad)
     X_noon = np.vstack([np.ones_like(clock_rad), x1_noon, x2_noon]).T  # design matrix with intercept
     
     coeffs_noon, residuals_noon, rank_noon, s_noon = lstsq(X, noon_medians, rcond=None)
     beta0_n, beta1_n, beta2_n = coeffs_noon
     
     A1_n = -beta1_n
     A2_n = 2 * beta2_n
     c_n = beta0_n - A1_n - (A2_n / 2)
     
     print(f"Original A1: {a_opt_noon}")
     print(f"Original A2: {b_opt_noon}")
     print(f"Original c: {c_opt_noon}")

     y_pred_n = X_noon @ coeffs_noon
     ss_res2_n = np.sum((noon_medians - y_pred_n) ** 2)
     ss_tot2_n = np.sum((noon_medians - np.mean(noon_medians)) ** 2)
     r_squared2n = 1 - ss_res2_n / ss_tot2_n

     print(f"R^2 (linear) = {r_squared2n:.4f}")
     print(f"A1 = {A1_n:.4f}, A2 = {A2_n:.4f}, c = {c_n:.4f}")
     '''
     rmse_n = np.sqrt(sn_res/22)
     sd_n = standard_deviation(noon_medians, np.mean(noon_medians))
     
     print(f"Mean: {np.mean(noon_medians):.4f}")
     print(f"Standard Deviation {sd_n:.4f}")
     print(f"RMSE: {rmse_n:.4f}")
     
     snp_res = np.sum((noon_medians - y_regress_noon_p)**2)
     snp_tot = np.sum((noon_medians - np.mean(noon_medians))**2)
     r_squared_np = 1 - (snp_res / sn_tot)
     print('----------------------')
     print('Noon Region +20% Fit')
     print(f"R squared: {r_squared_np:.4f}")
     
     adj_np = 1-r_squared_np
     adjr_np = 1 - ((adj_np*23)/22)
     print(f"Adjusted R squared: {adjr_np:.4f}")
     
     rmse_np = np.sqrt(snp_res/22)
     sd_np = standard_deviation(noon_medians, np.mean(noon_medians))
     
     print(f"Mean: {np.mean(noon_medians):.4f}")
     print(f"Standard Deviation {sd_np:.4f}")
     print(f"RMSE: {rmse_np:.4f}")
     
     snm_res = np.sum((noon_medians - y_regress_noon_m)**2)
     snm_tot = np.sum((noon_medians - np.mean(noon_medians))**2)
     r_squared_nm = 1 - (snm_res / snm_tot)
     print('----------------------')
     print('Noon -20% Region Fit')
     print(f"R squared: {r_squared_nm:.4f}")
     
     adj_nm = 1-r_squared_nm
     adjr_nm = 1 - ((adj_nm*23)/22)
     print(f"Adjusted R squared: {adjr_nm:.4f}")
     
     rmse_nm = np.sqrt(snm_res/22)
     sd_nm = standard_deviation(noon_medians, np.mean(noon_medians))
     
     print(f"Mean: {np.mean(noon_medians):.4f}")
     print(f"Standard Deviation {sd_nm:.4f}")
     print(f"RMSE: {rmse_nm:.4f}")

     
     sd_res = np.sum((dusk_medians - y_regress_dusk)**2)
     sd_tot = np.sum((dusk_medians - np.mean(dusk_medians))**2)
     r_squared_d = 1 - (sd_res / sd_tot)
     print('----------------------')
     print('Dusk Region Fit')
     print(f"R squared: {r_squared_d:.4f}")
     
     adj_d = 1-r_squared_d
     adjr_d = 1 - ((adj_d*23)/22)
     print(f"Adjusted R squared: {adjr_d:.4f}")
     
     rmse_d = np.sqrt(sd_res/22)
     sd_d = standard_deviation(dusk_medians, np.mean(dusk_medians))
     
     print(f"Mean: {np.mean(dusk_medians):.4f}")
     print(f"Standard Deviation {sd_d:.4f}")
     print(f"RMSE: {rmse_d:.4f}")
    
     
     sdp_res = np.sum((dusk_medians - y_regress_dusk_p)**2)
     sdp_tot = np.sum((dusk_medians - np.mean(dusk_medians))**2)
     r_squared_dp = 1 - (sdp_res / sdp_tot)
     print('----------------------')
     print('Dusk +20% Region Fit')
     print(f"R squared: {r_squared_dp:.4f}")
     
     adj_dp = 1-r_squared_dp
     adjr_dp = 1 - ((adj_dp*23)/22)
     print(f"Adjusted R squared: {adjr_dp:.4f}")
     
     rmse_dp = np.sqrt(sdp_res/22)
     sd_dp = standard_deviation(dusk_medians, np.mean(dusk_medians))
     
     print(f"Mean: {np.mean(dusk_medians):.4f}")
     print(f"Standard Deviation {sd_dp:.4f}")
     print(f"RMSE: {rmse_dp:.4f}")
     
     sdm_res = np.sum((dusk_medians - y_regress_dusk_m)**2)
     sdm_tot = np.sum((dusk_medians - np.mean(dusk_medians))**2)
     r_squared_dm = 1 - (sdm_res / sdm_tot)
     print('----------------------')
     print('Dusk -20% Region Fit')
     print(f"R squared: {r_squared_dm:.4f}")
     
     adj_dm = 1-r_squared_dm
     adjr_dm = 1 - ((adj_dm*23)/22)
     print(f"Adjusted R squared: {adjr_dm:.4f}")
     
     rmse_dm = np.sqrt(sdm_res/22)
     sd_dm = standard_deviation(dusk_medians, np.mean(dusk_medians))
     
     print(f"Mean: {np.mean(dusk_medians):.4f}")
     print(f"Standard Deviation {sd_dm:.4f}")
     print(f"RMSE: {rmse_dm:.4f}")

     '''
     plotting 
     '''
     
     if key == 'whole':
         func_name = 'Multi-Latitude Coupling Function'
         function = r'$A_1 \sin^2(\theta) + A_2 \cos^2\left(\frac{\theta}{2}\right) + C$'
     elif key == 'by':
         func_name = 'High-Latitude Coupling Function'
         function = r'$A \sin^2(\theta) + C$'
     elif key == 'earth':
         func_name = 'Dayside Coupling Function'
         function = r'$A \cos^2\left(\frac{\theta}{2}\right) + C$'
     elif key == 'tang':
         func_name = 'Post-Cusp Coupling Function'
         function = r'$A \sin^2\left(\frac{\theta}{2}\right) + C$'
     elif key == 'test':
         func_name = 'Test Function'
         function = r'$A \sin^2(\theta+90) + C$'
     else:
         print('not valid coupling function')
         
         
     text = f"{func_name}: {function}"


    
     fig = plt.figure(figsize=(28,55))
     ax1 = plt.subplot(4,2,1)
     plt.subplots_adjust(hspace=0.10, wspace=0.15)#hspace=0.20)
     #fig.text(0.5, 0.6914,text,fontsize=30,ha='center',color='dimgray')#, bbox=dict(facecolor='white', edgecolor='darkgray', boxstyle='round,pad=0.5'))
     
     ax1.scatter(Bz_pos[0:2], Bz_pos_power[0:2], color='red', s=200, label='+Bz CME')
     ax1.scatter(Bz_pos[2:], Bz_pos_power[2:], color='red', s=200, marker='v', label='+Bz Rarefaction (Deep)')
     
     ax1.scatter(Bz_neg[0], Bz_neg_power[0], color='orange', s=200, label='-Bz CME') #
     ax1.scatter(Bz_neg[1], Bz_neg_power[1], color='orange', s=200,marker=',', label='-Bz CIR') 
    
     ax1.scatter(By_pos[0], By_pos_power[0], color='green', s=200, label='+By CME') # compression 1 2-5 incl visit 15
     ax1.scatter(By_pos[1:4], By_pos_power[1:4], color='green', s=200, marker=",", label='+By CIR')
     ax1.scatter(By_pos[4:6], By_pos_power[4:6], color='green', s=200, marker="v", label='+By Rarefaction (Deep)') # rarefraction
     ax1.scatter(By_pos[6], By_pos_power[6], color='green', s=200, marker="^", label ='+By Rarefaction (Shallow)')
     ax1.scatter(By_pos[7:], By_pos_power[7:], color='green', s=200, marker=",")
    
     ax1.scatter(By_neg[0:4], By_neg_power[0:4], color='blue', s=200, marker="^", label='-By Rarefaction (Shallow)')
     ax1.scatter(By_neg[4:6], By_neg_power[4:6], color='blue', s=200, label='-By CME') # 
     ax1.scatter(By_neg[6], By_neg_power[6], color='blue', s=200, marker="v", label='-By Rarefaction (Deep)')
     ax1.scatter(By_neg[7:], By_neg_power[7:], color='blue', s=200, marker=",", label='-By CIR')

    
     ax1.text(165,930,'a',style='italic',fontsize=40) 
     #ax1.text(165,-2,'a',style='italic',fontsize=40) # whole
    
     ax1.errorbar(clock_medians, polar_medians, yerr=polar_med_errs, xerr=clock_med_errs, fmt='.', color='lightgray')
     ax1.set_xlabel('Median Clock Angle ($^o$)',fontsize=22)
     ax1.set_ylabel('Median Total Polar Power (GW)',fontsize=22)
     
     ax1.set_xlim(-190, 190)
     ax1.set_ylim(-30,1015)
     ax1.tick_params(axis='x',which='minor',direction='in')
     ax1.yaxis.set_minor_locator(ticker.AutoMinorLocator())
     ax1.minorticks_on()
     ax1.tick_params(which='major',direction='in',bottom=True, top=True, left=True, right=True, labelsize=20, length=10)
     ax1.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True, length=5)
     lgnd = plt.legend(framealpha=0.5, labelspacing=0.3, fontsize=18, loc="upper left")#, ncol=2)
     for handle in lgnd.legend_handles:
         handle.set_sizes([100.0])
    
    
     ax2 = plt.subplot(4,2,2)
     
     # ax2.plot(x_new, y_new_swirl, '--',color='black', linewidth=3, markersize=12, label=(f'Swirl $R^2$: {adjr_s:.2g}'), zorder=10)

     if error_plot == '20' or error_plot == '10' or error_plot == '50':
         # plus error
         ax2.scatter(clock_medians_p[0:2], medians_swirl[0:2],s=100,color='plum',marker="^")
         ax2.scatter(clock_medians_p[2], medians_swirl[2],s=100,color='plum')
         ax2.scatter(clock_medians_p[3:5], medians_swirl[3:5],s=100,color='plum',marker="^")
         ax2.scatter(clock_medians_p[5:9], medians_swirl[5:9],s=100,color='plum')
         ax2.scatter(clock_medians_p[9], medians_swirl[9],s=100,color='plum',marker="v")
         ax2.scatter(clock_medians_p[10], medians_swirl[10],s=100,color='plum')
         ax2.scatter(clock_medians_p[11], medians_swirl[11],s=100,color='plum',marker="v")
         ax2.scatter(clock_medians_p[12:18], medians_swirl[12:18],s=100,color='plum',marker=",")
         ax2.scatter(clock_medians_p[18:21], medians_swirl[18:21],s=100,color='plum',marker="v")
         ax2.scatter(clock_medians_p[21], medians_swirl[21],s=100,color='plum',marker="^")
         ax2.scatter(clock_medians_p[22:], medians_swirl[22:],s=100,color='plum',marker=",")
     
        # minus error
         ax2.scatter(clock_medians_m[0:2], medians_swirl[0:2],s=100,color='peachpuff',marker="^")
         ax2.scatter(clock_medians_m[2], medians_swirl[2],s=100,color='peachpuff')
         ax2.scatter(clock_medians_m[3:5], medians_swirl[3:5],s=100,color='peachpuff',marker="^")
         ax2.scatter(clock_medians_m[5:9], medians_swirl[5:9],s=100,color='peachpuff')
         ax2.scatter(clock_medians_m[9], medians_swirl[9],s=100,color='peachpuff',marker="v")
         ax2.scatter(clock_medians_m[10], medians_swirl[10],s=100,color='peachpuff')
         ax2.scatter(clock_medians_m[11], medians_swirl[11],s=100,color='peachpuff',marker="v")
         ax2.scatter(clock_medians_m[12:18], medians_swirl[12:18],s=100,color='peachpuff',marker=",")
         ax2.scatter(clock_medians_m[18:21], medians_swirl[18:21],s=100,color='peachpuff',marker="v")
         ax2.scatter(clock_medians_m[21], medians_swirl[21],s=100,color='peachpuff',marker="^")
         ax2.scatter(clock_medians_m[22:], medians_swirl[22:],s=100,color='peachpuff',marker=",")
    
         ax2.errorbar(clock_medians_p, medians_swirl, yerr=swirl_med_errs, xerr=clock_med_errs_p, fmt='.', color='plum')

         ax2.plot(x_new, y_new_swirl_p, '--',color='darkviolet', linewidth=3, markersize=12,label=f'+ {error}%: {adjr_sp:.2g}')
     
         ax2.errorbar(clock_medians_m, medians_swirl, yerr=swirl_med_errs, xerr=clock_med_errs_m, fmt='.', color='peachpuff')
 
         ax2.plot(x_new, y_new_swirl_m, '--',color='sandybrown', linewidth=3, markersize=12,label=f'- {error}%: {adjr_sm:.2g}')
         
     
     # main plot
     ax2.scatter(Bz_pos[0:2], Bz_pos_swirlpower[0:2], color='red', s=200,zorder=9)
     ax2.scatter(Bz_pos[2:], Bz_pos_swirlpower[2:], color='red', s=200, marker='v',zorder=9)
     
     ax2.scatter(Bz_neg[0], Bz_neg_swirlpower[0], color='orange', s=200,zorder=9) #
     ax2.scatter(Bz_neg[1], Bz_neg_swirlpower[1], color='orange', s=200,marker=',',zorder=9) 
    
     ax2.scatter(By_pos[0], By_pos_swirlpower[0], color='green', s=200,zorder=9) # compression 1 2-5 incl visit 15
     ax2.scatter(By_pos[1:4], By_pos_swirlpower[1:4], color='green', s=200, marker=",",zorder=9)
     ax2.scatter(By_pos[4:6], By_pos_swirlpower[4:6], color='green', s=200, marker="v",zorder=9) # rarefraction
     ax2.scatter(By_pos[6], By_pos_swirlpower[6], color='green', s=200, marker="^",zorder=9)
     ax2.scatter(By_pos[7:], By_pos_swirlpower[7:], color='green', s=200, marker=",",zorder=9)
    
     ax2.scatter(By_neg[0:4], By_neg_swirlpower[0:4], color='blue', s=200, marker="^",zorder=9)
     ax2.scatter(By_neg[4:6], By_neg_swirlpower[4:6], color='blue', s=200,zorder=9) # 
     ax2.scatter(By_neg[6], By_neg_swirlpower[6], color='blue', s=200, marker="v",zorder=9)
     ax2.scatter(By_neg[7:], By_neg_swirlpower[7:], color='blue', s=200, marker=",",zorder=9)


     #ax2.text(165,-0.2,'b',style='italic',fontsize=40) #-0.2, -0.7 (whole), -0.3, -0.1 (tang nooff)
     ax2.text(165,80.5,'b',style='italic',fontsize=40) #78 (by nooff+off, err), 80.5/81 (earth nooff+off and tangoff w/err, whole nooff), 71.5 (earth&by nooff/off,noerr+whole,noerr,nooff), 86 (tangnoofferr), 77 (tangnooff), 68.8 (no plots)
     x_min_2, x_max_2 = ax2.get_xlim()
     x_mid_2 = (x_min_2 + x_max_2) / 2

     # Place the text at the midpoint
     # ax2.text(x_mid_2,1,f'Swirl Fit: {swirl_fit}',fontsize=20,ha='center')
     # ax2.text(x_mid_2,2,r'Swirl Fit: $A_1 \sin^2(\theta) + A_2 \cos^2\left(\frac{\theta}{2}\right) + C$',fontsize=20,ha='center',bbox=dict(facecolor='white', edgecolor='darkgray', boxstyle='round,pad=0.5'))
    
     ax2.errorbar(clock_medians, medians_swirl, yerr=swirl_med_errs, xerr=clock_med_errs, fmt='.', color='lightgray',zorder=9)
    
     ax2.set_xlabel('Median Clock Angle ($^o$)',fontsize=22)
     ax2.set_ylabel('Median Swirl Region Power (GW)',fontsize=22)
    
     ax2.set_ylim(-3,88) # 88 & 78 (whole&nooff,by&earth-nooff/off&noerr, tangoffnoerr), 94 (tang,nooff), 84 (tang,nooff,noerr)
     ax2.set_xlim(-190, 190)
     ax2.tick_params(axis='x',which='minor',direction='in')
     ax2.yaxis.set_minor_locator(ticker.AutoMinorLocator())
     ax2.minorticks_on()
     ax2.tick_params(which='major',direction='in',bottom=True, top=True, left=True, right=True, labelsize=20, length=10)
     ax2.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True, length=5)
     ax2.legend(framealpha=0.5,fontsize=21, loc ="upper left")
    
     #########
    
     ax3 = plt.subplot(4,2,3)
     
     # ax3.plot(x_new, y_new_dusk,'--', color='black', linewidth=3, markersize=12,label=f'Dusk $R^2$: {adjr_d:.2g}', zorder=10)
      
     if error_plot == '20' or error_plot == '10' or error_plot == '50':
         # plus error
         ax3.scatter(clock_medians_p[0:2], dusk_medians[0:2],s=100,color='plum',marker="^")
         ax3.scatter(clock_medians_p[2], dusk_medians[2],s=100,color='plum')
         ax3.scatter(clock_medians_p[3:5], dusk_medians[3:5],s=100,color='plum',marker="^")
         ax3.scatter(clock_medians_p[5:9], dusk_medians[5:9],s=100,color='plum')
         ax3.scatter(clock_medians_p[9], dusk_medians[9],s=100,color='plum',marker="v")
         ax3.scatter(clock_medians_p[10], dusk_medians[10],s=100,color='plum')
         ax3.scatter(clock_medians_p[11], dusk_medians[11],s=100,color='plum',marker="v")
         ax3.scatter(clock_medians_p[12:18], dusk_medians[12:18],s=100,color='plum',marker=",")
         ax3.scatter(clock_medians_p[18:21], dusk_medians[18:21],s=100,color='plum',marker="v")
         ax3.scatter(clock_medians_p[21], dusk_medians[21],s=100,color='plum',marker="^")
         ax3.scatter(clock_medians_p[22:], dusk_medians[22:],s=100,color='plum',marker=",")
    
         # minus error
         ax3.scatter(clock_medians_m[0:2], dusk_medians[0:2],s=100,color='peachpuff',marker="^")
         ax3.scatter(clock_medians_m[2], dusk_medians[2],s=100,color='peachpuff')
         ax3.scatter(clock_medians_m[3:5], dusk_medians[3:5],s=100,color='peachpuff',marker="^")
         ax3.scatter(clock_medians_m[5:9], dusk_medians[5:9],s=100,color='peachpuff')
         ax3.scatter(clock_medians_m[9], dusk_medians[9],s=100,color='peachpuff',marker="v")
         ax3.scatter(clock_medians_m[10], dusk_medians[10],s=100,color='peachpuff')
         ax3.scatter(clock_medians_m[11], dusk_medians[11],s=100,color='peachpuff',marker="v")
         ax3.scatter(clock_medians_m[12:18], dusk_medians[12:18],s=100,color='peachpuff',marker=",")
         ax3.scatter(clock_medians_m[18:21], dusk_medians[18:21],s=100,color='peachpuff',marker="v")
         ax3.scatter(clock_medians_m[21], dusk_medians[21],s=100,color='peachpuff',marker="^")
         ax3.scatter(clock_medians_m[22:], dusk_medians[22:],s=100,color='peachpuff',marker=",")
    
         ax3.errorbar(clock_medians_p, dusk_medians, yerr=dusk_med_errs, xerr=clock_med_errs_p, fmt='.', color='plum')
    
         ax3.plot(x_new, y_new_dusk_p, '--',color='darkviolet', linewidth=3, markersize=12,label=f'+ {error}%: {adjr_dp:.2g}')
         
         ax3.errorbar(clock_medians_m, dusk_medians, yerr=dusk_med_errs, xerr=clock_med_errs_m, fmt='.', color='peachpuff')
    
         ax3.plot(x_new, y_new_dusk_m, '--',color='sandybrown', linewidth=3, markersize=12,label=f'- {error}%: {adjr_dm:.2g}')
     
     
     # main plot
     ax3.scatter(Bz_pos[0:2], Bz_pos_duskpower[0:2], color='red', s=200,zorder=9)
     ax3.scatter(Bz_pos[2:], Bz_pos_duskpower[2:], color='red', s=200, marker='v',zorder=9)
     
     ax3.scatter(Bz_neg[0], Bz_neg_duskpower[0], color='orange', s=200,zorder=9) #
     ax3.scatter(Bz_neg[1], Bz_neg_duskpower[1], color='orange', s=200,marker=',',zorder=9)
    
     ax3.scatter(By_pos[0], By_pos_duskpower[0], color='green', s=200,zorder=9) # compression 1 2-5 incl visit 15
     ax3.scatter(By_pos[1:4], By_pos_duskpower[1:4], color='green', s=200, marker=",",zorder=9)
     ax3.scatter(By_pos[4:6], By_pos_duskpower[4:6], color='green', s=200, marker="v",zorder=9) # rarefraction
     ax3.scatter(By_pos[6], By_pos_duskpower[6], color='green', s=200, marker="^",zorder=9)
     ax3.scatter(By_pos[7:], By_pos_duskpower[7:], color='green', s=200, marker=",",zorder=9)
    
     ax3.scatter(By_neg[0:4], By_neg_duskpower[0:4], color='blue', s=200, marker="^",zorder=9)
     ax3.scatter(By_neg[4:6], By_neg_duskpower[4:6], color='blue', s=200,zorder=9) # 
     ax3.scatter(By_neg[6], By_neg_duskpower[6], color='blue', s=200, marker="v",zorder=9)
     ax3.scatter(By_neg[7:], By_neg_duskpower[7:], color='blue', s=200, marker=",",zorder=9)
     
     
     ax3.errorbar(clock_medians, dusk_medians, yerr=dusk_med_errs, xerr=clock_med_errs, fmt='.', color='lightgray', zorder=9)
     
     #ax3.text(165,-4,'c',style='italic',fontsize=40) # -4 /-2 (whole), -15 ()
     ax3.text(165,397,'c',style='italic',fontsize=40) #405 (by&earth nooff+off), 396/5 (earth no err/err, nooff), 397/5 (by no err/err,nooff), 480 (whole,nooff), 422 (whole,nooff,noerr), 395 (by/earth/tang off)

     ax3.set_xlabel('Median Clock Angle ($^o$)', fontsize=22)
     ax3.set_ylabel('Median Dusk Active Region Power (GW)',fontsize=22)
     
     ax3.set_xlim(-190, 190)
     #ax3.set_ylim(-20,520) #520 & 460 (whole), none (tang,by&earth nooff +off)
     ax3.tick_params(axis='x',which='minor',direction='in')
     ax3.yaxis.set_minor_locator(ticker.AutoMinorLocator())
     ax3.minorticks_on()
     ax3.tick_params(which='major',direction='in',bottom=True, top=True, left=True, right=True, labelsize=20, length=10)
     ax3.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True, length=5)
     ax3.legend(framealpha=0.5,fontsize=21,loc ="upper left")
    
    
     ########
    
     ax4 = plt.subplot(4,2,4)
     
     # ax4.plot(x_new, y_new_noon, '--',color='black', linewidth=3, markersize=12,label=f'Noon $R^2$: {adjr_n:.2g}', zorder=10)
    
     if error_plot == '20' or error_plot == '10' or error_plot == '50':
         # plus error
         ax4.scatter(clock_medians_p[0:2], noon_medians[0:2],s=100,color='plum',marker="^")
         ax4.scatter(clock_medians_p[2], noon_medians[2],s=100,color='plum')
         ax4.scatter(clock_medians_p[3:5], noon_medians[3:5],s=100,color='plum',marker="^")
         ax4.scatter(clock_medians_p[5:9], noon_medians[5:9],s=100,color='plum')
         ax4.scatter(clock_medians_p[9], noon_medians[9],s=100,color='plum',marker="v")
         ax4.scatter(clock_medians_p[10], noon_medians[10],s=100,color='plum')
         ax4.scatter(clock_medians_p[11], noon_medians[11],s=100,color='plum',marker="v")
         ax4.scatter(clock_medians_p[12:18], noon_medians[12:18],s=100,color='plum',marker=",")
         ax4.scatter(clock_medians_p[18:21], noon_medians[18:21],s=100,color='plum',marker="v")
         ax4.scatter(clock_medians_p[21], noon_medians[21],s=100,color='plum',marker="^")
         ax4.scatter(clock_medians_p[22:], noon_medians[22:],s=100,color='plum',marker=",")

         # minus error
         ax4.scatter(clock_medians_m[0:2], noon_medians[0:2],s=100,color='peachpuff',marker="^")
         ax4.scatter(clock_medians_m[2], noon_medians[2],s=100,color='peachpuff')
         ax4.scatter(clock_medians_m[3:5], noon_medians[3:5],s=100,color='peachpuff',marker="^")
         ax4.scatter(clock_medians_m[5:9], noon_medians[5:9],s=100,color='peachpuff')
         ax4.scatter(clock_medians_m[9], noon_medians[9],s=100,color='peachpuff',marker="v")
         ax4.scatter(clock_medians_m[10], noon_medians[10],s=100,color='peachpuff')
         ax4.scatter(clock_medians_m[11], noon_medians[11],s=100,color='peachpuff',marker="v")
         ax4.scatter(clock_medians_m[12:18], noon_medians[12:18],s=100,color='peachpuff',marker=",")
         ax4.scatter(clock_medians_m[18:21], noon_medians[18:21],s=100,color='peachpuff',marker="v")
         ax4.scatter(clock_medians_m[21], noon_medians[21],s=100,color='peachpuff',marker="^")
         ax4.scatter(clock_medians_m[22:], noon_medians[22:],s=100,color='peachpuff',marker=",")
    
         ax4.errorbar(clock_medians_p, noon_medians, yerr=noon_med_errs, xerr=clock_med_errs_p, fmt='.', color='plum')
         ax4.plot(x_new, y_new_noon_p, '--',color='darkviolet', linewidth=3, markersize=12,label=f'+ {error}%: {adjr_np:.2g}')
     
         ax4.errorbar(clock_medians_m, noon_medians, yerr=noon_med_errs, xerr=clock_med_errs_m, fmt='.', color='peachpuff')
         ax4.plot(x_new, y_new_noon_m, '--',color='sandybrown', linewidth=3, markersize=12,label=f'- {error}%: {adjr_nm:.2g}')
     
     
     # main plot
     ax4.scatter(Bz_pos[0:2], Bz_pos_noonpower[0:2], color='red', s=200,zorder=9)
     ax4.scatter(Bz_pos[2:], Bz_pos_noonpower[2:], color='red', s=200, marker='v',zorder=9)
     
     ax4.scatter(Bz_neg[0], Bz_neg_noonpower[0], color='orange', s=200,zorder=9) #
     ax4.scatter(Bz_neg[1], Bz_neg_noonpower[1], color='orange', s=200,marker=',',zorder=9) 
    
     ax4.scatter(By_pos[0], By_pos_noonpower[0], color='green', s=200,zorder=9) # compression 1 2-5 incl visit 15
     ax4.scatter(By_pos[1:4], By_pos_noonpower[1:4], color='green', s=200, marker=",",zorder=9)
     ax4.scatter(By_pos[4:6], By_pos_noonpower[4:6], color='green', s=200, marker="v",zorder=9) # rarefraction
     ax4.scatter(By_pos[6], By_pos_noonpower[6], color='green', s=200, marker="^",zorder=9)
     ax4.scatter(By_pos[7:], By_pos_noonpower[7:], color='green', s=200, marker=",",zorder=9)
    
     ax4.scatter(By_neg[0:4], By_neg_noonpower[0:4], color='blue', s=200, marker="^",zorder=9)
     ax4.scatter(By_neg[4:6], By_neg_noonpower[4:6], color='blue', s=200,zorder=9)
     ax4.scatter(By_neg[6], By_neg_noonpower[6], color='blue', s=200, marker="v",zorder=9)
     ax4.scatter(By_neg[7:], By_neg_noonpower[7:], color='blue', s=200, marker=",",zorder=9)


     #ax4.text(165,1,'d',style='italic',fontsize=40) #-2, 0.5 (whole), -12 (tang,nooff)
     ax4.text(165,400,'d',style='italic',fontsize=40) #400 (by&earth nooff+off and tang off), 480 (whole, nooff), 422 (whole,noerr,nooff)
     
     ax4.errorbar(clock_medians, noon_medians, yerr=noon_med_errs, xerr=clock_med_errs, fmt='.', color='lightgray', zorder=9)

     ax4.set_ylabel('Median Noon Active Region Power (GW)',fontsize=22)
     ax4.set_xlabel('Median Clock Angle ($^o$)',fontsize=22)
     
     ax4.set_xlim(-190, 190)
     #ax4.set_ylim(-15,520)# 520, 460 (whole), none (tang,by&earth for off and nooff)
     ax4.tick_params(axis='x',which='minor',direction='in')
     ax4.yaxis.set_minor_locator(ticker.AutoMinorLocator())
     ax4.minorticks_on()
     ax4.tick_params(which='major',direction='in',bottom=True, top=True, left=True, right=True, labelsize=20, length=10)
     ax4.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True, length=5)
     ax4.legend(framealpha=0.5,fontsize=21,loc ="upper left")
    
    
     # save plot
     saveloc = (f'{root_saves}median_clock_vs_region_power_{key}_{error_plot}.jpg') 
     #saveloc = (f'{root_saves}median_clock_vs_region_power_extra3_{error}.jpg')
     plt.savefig(saveloc,bbox_inches='tight',dpi=400)
    

if plotting == 'pressure':
    
    # pressure_medians = np.log10(pressure_medians)
    # pressure_medians_m = np.log10(pressure_medians_m)
    # pressure_medians_p = np.log10(pressure_medians_p)
    
    # med_pressure_errs = np.log10(med_pressure_errs)
    # pressure_med_errs_m = np.log10(pressure_med_errs_m)
    # pressure_med_errs_p = np.log10(pressure_med_errs_p)

    
    a, acov = np.polyfit(np.log10(pressure_medians), medians_swirl, 1, cov=True)
    aa = np.poly1d(a)
    # print("Swirl Fit Gradient:")
    # print(a[0])
    # print("Error of Swirl Fit:")
    # print(np.sqrt(np.diag(acov)))
    y_regress_swirl = a[0]*np.log10(pressure_medians)+a[1]
    
    print(f"Testing '{plotting}'")
    print('Swirl Region Fit')
    ss_res = np.sum((medians_swirl - y_regress_swirl)**2)
    ss_tot = np.sum((medians_swirl - np.mean(medians_swirl))**2)
    r_squared_s = 1 - (ss_res / ss_tot)
    print(f"R squared: {r_squared_s:.4f}")
    
    adj_s = 1-r_squared_s
    adjr_s = 1 - ((adj_s*23)/22)
    print(f"Adjusted R squared: {adjr_s:.4f}")
    
    swirl_grad = round(a[0],2)
    swirl_grad_err = np.sqrt(np.diag(acov))
    swirl_grad_err = round(swirl_grad_err[0], 2)
    
    swirl_pres_fit = (f'{swirl_grad} ± {swirl_grad_err}')
    
    
    a_p,acov_p = np.polyfit(np.log10(pressure_medians_p), medians_swirl,1, cov=True)
    aa_p = np.poly1d(a_p)
    # print(f"Swirl Fit +{error} Gradient:")
    # print(a_p[0])
    # print(f"Error of Swirl Fit +{error}:")
    # print(np.sqrt(np.diag(acov_p)))
    y_regress_swirl_p = a_p[0]*np.log10(pressure_medians_p)+a_p[1]
    
    print("---------------")
    print('Swirl Region +20% Fit')
    ssp_res = np.sum((medians_swirl - y_regress_swirl_p)**2)
    ssp_tot = np.sum((medians_swirl - np.mean(medians_swirl))**2)
    r_squared_sp = 1 - (ssp_res / ssp_tot)
    print(f"R squared: {r_squared_sp:.4f}")
    
    adj_sp = 1-r_squared_sp
    adjr_sp = 1 - ((adj_sp*23)/22)
    print(f"Adjusted R squared: {adjr_sp:.4f}")
    
    swirl_grad_p = round(a_p[0],2)
    swirl_grad_err_p = np.sqrt(np.diag(acov_p))
    swirl_grad_err_p = round(swirl_grad_err_p[0], 2)
    
    swirl_pres_fit_p = (f'{swirl_grad_p} ± {swirl_grad_err_p}')
    
    
    a_m,acov_m = np.polyfit(np.log10(pressure_medians_m), medians_swirl,1, cov=True)
    aa_m = np.poly1d(a_m)
    # print(f"Swirl Fit -{error} Gradient:")
    # print(a_m[0])
    # print(f"Error of Swirl Fit -{error}:")
    # print(np.sqrt(np.diag(acov_m)))
    y_regress_swirl_m = a_m[0]*np.log10(pressure_medians_m)+a_m[1]
    
    print("---------------")
    print('Swirl Region -20% Fit')
    ssm_res = np.sum((medians_swirl - y_regress_swirl_m)**2)
    ssm_tot = np.sum((medians_swirl - np.mean(medians_swirl))**2)
    r_squared_sm = 1 - (ssm_res / ssm_tot)
    print(f"R squared: {r_squared_sm:.4f}")
    
    adj_sm = 1-r_squared_sm
    adjr_sm = 1 - ((adj_sm*23)/22)
    print(f"Adjusted R squared: {adjr_sm:.4f}")
    
    swirl_grad_m = round(a_m[0],2)
    swirl_grad_err_m = np.sqrt(np.diag(acov_m))
    swirl_grad_err_m = round(swirl_grad_err_m[0], 2)
    
    swirl_pres_fit_m = (f'{swirl_grad_m} ± {swirl_grad_err_m}')
   
    
    b, bcov = np.polyfit(np.log10(pressure_medians), dusk_medians, 1, cov=True)
    bb = np.poly1d(b)
    # print("Dusk Fit Gradient:")
    # print(b[0])
    # print("Error of Dusk Fit:")
    # print(np.sqrt(np.diag(bcov)))
    y_regress_dusk = b[0]*np.log10(pressure_medians)+b[1]
    
    print("---------------")
    print('Dusk Region Fit')
    sd_res = np.sum((dusk_medians - y_regress_dusk)**2)
    sd_tot = np.sum((dusk_medians - np.mean(dusk_medians))**2)
    r_squared_d = 1 - (sd_res / sd_tot)
    print(f"R squared: {r_squared_d:.4f}")
    
    adj_d = 1-r_squared_d
    adjr_d = 1 - ((adj_d*23)/22)
    print(f"Adjusted R squared: {adjr_d:.4f}")
    
    dusk_grad = round(b[0],2)
    dusk_grad_err = np.sqrt(np.diag(bcov))
    dusk_grad_err = round(dusk_grad_err[0], 2)
    
    dusk_pres_fit = (f'{dusk_grad} ± {dusk_grad_err}')
    
    
    b_p, bcov_p = np.polyfit(np.log10(pressure_medians_p), dusk_medians, 1, cov=True)
    bb_p = np.poly1d(b_p)
    # print(f"Dusk Fit +{error} Gradient:")
    # print(b_p[0])
    # print(f"Error of Dusk Fit +{error}:")
    # print(np.sqrt(np.diag(bcov_p)))
    y_regress_dusk_p = b_p[0]*np.log10(pressure_medians_p)+b_p[1]
    
    print("---------------")
    print('Dusk Region +20% Fit')
    sdp_res = np.sum((dusk_medians - y_regress_dusk_p)**2)
    sdp_tot = np.sum((dusk_medians - np.mean(dusk_medians))**2)
    r_squared_dp = 1 - (sdp_res / sdp_tot)
    print(f"R squared: {r_squared_dp:.4f}")
    
    adj_dp = 1-r_squared_dp
    adjr_dp = 1 - ((adj_dp*23)/22)
    print(f"Adjusted R squared: {adjr_dp:.4f}")
    
    dusk_grad_p = round(b_p[0],2)
    dusk_grad_err_p = np.sqrt(np.diag(bcov_p))
    dusk_grad_err_p = round(dusk_grad_err_p[0], 2)
    
    dusk_pres_fit_p = (f'{dusk_grad_p} ± {dusk_grad_err_p}')
    
    
    b_m, bcov_m = np.polyfit(np.log10(pressure_medians_m), dusk_medians, 1, cov=True)
    bb_m = np.poly1d(b_m)
    # print(f"Dusk Fit -{error} Gradient:")
    # print(b_m[0])
    # print(f"Error of Dusk Fit -{error}:")
    # print(np.sqrt(np.diag(bcov_m)))
    y_regress_dusk_m = b_m[0]*np.log10(pressure_medians_m)+b_m[1]
    
    print("---------------")
    print('Dusk Region -20% Fit')
    sdm_res = np.sum((dusk_medians - y_regress_dusk_m)**2)
    sdm_tot = np.sum((dusk_medians - np.mean(dusk_medians))**2)
    r_squared_dm = 1 - (sdm_res / sdm_tot)
    print(f"R squared: {r_squared_dm:.4f}")
    
    adj_dm = 1-r_squared_dm
    adjr_dm = 1 - ((adj_dm*23)/22)
    print(f"Adjusted R squared: {adjr_dm:.4f}")
    
    dusk_grad_m = round(b_m[0],2)
    dusk_grad_err_m = np.sqrt(np.diag(bcov_m))
    dusk_grad_err_m = round(dusk_grad_err_m[0], 2)
    
    dusk_pres_fit_m = (f'{dusk_grad_m} ± {dusk_grad_err_m}')
    
    
    
    c, ccov = np.polyfit(np.log10(pressure_medians), noon_medians, 1, cov=True)
    cc = np.poly1d(c)
    # print("Noon Fit Gradient:")
    # print(c[0])
    # print("Error of Noon Fit:")
    # print(np.sqrt(np.diag(ccov)))
    y_regress_noon = c[0]*np.log10(pressure_medians)+c[1]
    
    y_regress_noon = cc(btot_medians)

    print("---------------")
    print('Noon Region Fit')
    sn_res = np.sum((noon_medians - y_regress_noon)**2)
    sn_tot = np.sum((noon_medians - np.mean(noon_medians))**2)
    r_squared_n = 1 - (sn_res / sn_tot)
    print(f"R squared: {r_squared_n:.4f}")
    
    adj_n = 1-r_squared_n
    adjr_n = 1 - ((adj_n*23)/22)
    print(f"Adjusted R squared: {adjr_n:.4f}")
    
    noon_grad = round(c[0],2)
    noon_grad_err = np.sqrt(np.diag(ccov))
    noon_grad_err = round(noon_grad_err[0], 2)
    
    noon_pres_fit = (f'{noon_grad} ± {noon_grad_err}')
    
    
    c_p, ccov_p = np.polyfit(np.log10(pressure_medians_p), noon_medians, 1, cov=True)
    cc_p = np.poly1d(c_p)
    # print(f"Noon Fit +{error} Gradient:")
    # print(c_p[0])
    # print(f"Error of Noon Fit +{error}:")
    # print(np.sqrt(np.diag(ccov_p)))
    y_regress_noon_p = c_p[0]*np.log10(pressure_medians_p)+c_p[1]

    print("---------------")
    print('Noon Region +20% Fit')
    snp_res = np.sum((noon_medians - y_regress_noon_p)**2)
    snp_tot = np.sum((noon_medians - np.mean(noon_medians))**2)
    r_squared_np = 1 - (snp_res / snp_tot)
    print(f"R squared: {r_squared_np:.4f}")
    
    adj_np = 1-r_squared_np
    adjr_np = 1 - ((adj_np*23)/22)
    print(f"Adjusted R squared: {adjr_np:.4f}")
    
    noon_grad_p = round(c_p[0],2)
    noon_grad_err_p = np.sqrt(np.diag(ccov_p))
    noon_grad_err_p = round(noon_grad_err_p[0], 2)
    
    noon_pres_fit_p = (f'{noon_grad_p} ± {noon_grad_err_p}')
    
    
    c_m, ccov_m = np.polyfit(np.log10(pressure_medians_m), noon_medians, 1, cov=True)
    cc_m = np.poly1d(c_m)
    # print(f"Noon Fit -{error} Gradient:")
    # print(c_m[0])
    # print(f"Error of Noon Fit -{error}:")
    # print(np.sqrt(np.diag(ccov_m)))
    y_regress_noon_m = c_m[0]*np.log10(pressure_medians_m)+c_m[1]
    
    print("---------------")
    print('Noon Region -20% Fit')
    snm_res = np.sum((noon_medians - y_regress_noon_m)**2)
    snm_tot = np.sum((noon_medians - np.mean(noon_medians))**2)
    r_squared_nm = 1 - (snm_res / snm_tot)
    print(f"R squared: {r_squared_nm:.4f}")
    
    adj_nm = 1-r_squared_nm
    adjr_nm = 1 - ((adj_nm*23)/22)
    print(f"Adjusted R squared: {adjr_nm:.4f}")
    
    noon_grad_m = round(c_m[0],2)
    noon_grad_err_m = np.sqrt(np.diag(ccov_m))
    noon_grad_err_m = round(noon_grad_err_m[0], 2)
    
    noon_pres_fit_m = (f'{noon_grad_m} ± {noon_grad_err_m}')
    
    
    # arrange cml more evenly for fit plotting
    x_new = np.linspace(min(pressure_medians_m), max(pressure_medians_m), len(pressure_medians))

    
    y_new_swirl = a[0]*np.log10(x_new)+a[1]
    y_new_dusk = b[0]*np.log10(x_new)+b[1]
    y_new_noon = c[0]*np.log10(x_new)+c[1]
    
    y_new_swirl_p = a_p[0]*np.log10(x_new)+a_p[1]
    y_new_dusk_p = b_p[0]*np.log10(x_new)+b_p[1]
    y_new_noon_p = c_p[0]*np.log10(x_new)+c_p[1]
    
    y_new_swirl_m = a_m[0]*np.log10(x_new)+a_m[1]
    y_new_dusk_m = b_m[0]*np.log10(x_new)+b_m[1]
    y_new_noon_m = c_m[0]*np.log10(x_new)+c_m[1]
    

    fig = plt.figure(figsize=(28,50))
    ax1 = plt.subplot(4,2,1)
    plt.subplots_adjust(hspace=0.1, wspace=0.15)

    # main plot
    #ax1.semilogx(pressure_medians, polar_medians/(1e12))   
    
    ax1.scatter(Bz_pos_pressure[0:2], Bz_pos_power[0:2], color='red', s=200, label='+Bz CME')
    ax1.scatter(Bz_pos_pressure[2:], Bz_pos_power[2:], color='red', s=200, marker='v', label='+Bz Rarefaction (Deep)')

    ax1.scatter(Bz_neg_pressure[0], Bz_neg_power[0], color='orange', s=200, label='-Bz CME') #
    ax1.scatter(Bz_neg_pressure[1], Bz_neg_power[1], color='orange', s=200,marker=',', label='-Bz CIR') 

    ax1.scatter(By_pos_pressure[0], By_pos_power[0], color='green', s=200, label='+By CME') # compression 1 2-5 incl visit 15
    ax1.scatter(By_pos_pressure[1:4], By_pos_power[1:4], color='green', s=200, marker=",", label='+By CIR')
    ax1.scatter(By_pos_pressure[4:6], By_pos_power[4:6], color='green', s=200, marker="v", label='+By Rarefaction (Deep)') # rarefraction
    ax1.scatter(By_pos_pressure[6], By_pos_power[6], color='green', s=200, marker="^", label ='+By Rarefaction (Shallow)')
    ax1.scatter(By_pos_pressure[7:], By_pos_power[7:], color='green', s=200, marker=",")

    ax1.scatter(By_neg_pressure[0:4], By_neg_power[0:4], color='blue', s=200, marker="^", label='-By Rarefaction (Shallow)')
    ax1.scatter(By_neg_pressure[4:6], By_neg_power[4:6], color='blue', s=200, label='-By CME') # 
    ax1.scatter(By_neg_pressure[6], By_neg_power[6], color='blue', s=200, marker="v", label='-By Rarefaction (Deep)')
    ax1.scatter(By_neg_pressure[7:], By_neg_power[7:], color='blue', s=200, marker=",", label='-By CIR')

    
    ax1.errorbar(pressure_medians, polar_medians, yerr=polar_med_errs, xerr=med_pressure_errs, fmt='.', color='lightgray')
    
    ax1.text(3.3e-1,-11,'a',style='italic',fontsize=40) #0.32
    
    ax1.set_xlabel('Median Pressure (nPa)', fontsize=22)
    ax1.set_ylabel('Median Total Polar Power (GW)', fontsize=22)
   
    ax1.set_xlim(3.4e-3, 4.5e-1)
    #ax1.set_xlim(0.005, 0.36)
    #ax1.set_ylim(-20,990)
    ax1.set_xscale('log')
    locmaj = matplotlib.ticker.LogLocator(base=10,numticks=30) 
    ax1.xaxis.set_major_locator(locmaj)
    ax1.tick_params(axis='x',which='minor',direction='in')
    ax1.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax1.minorticks_on()
    ax1.tick_params(which='major',direction='in',bottom=True, top=True, left=True, right=True, labelsize=20, length=10)
    ax1.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True, length=5)
    lgnd = plt.legend(framealpha=0.5, labelspacing=0.3, fontsize=18, loc="upper left")
    for handle in lgnd.legend_handles:
        handle.set_sizes([100.0])


    ax2 = plt.subplot(4,2,2)
    
    #ax2.semilogx(pressure_medians, medians_swirl/(1e12)) 
    ax2.plot(x_new, y_new_swirl, '--',color='black',linewidth=3, markersize=12, label=(f'$R^2$: {adjr_s:.2g}'), zorder=10)
    
    ax2.errorbar(pressure_medians, medians_swirl, yerr=swirl_med_errs, xerr=med_pressure_errs, fmt='.', color='lightgray', zorder=10)
    
    if error_plot == '20' or error_plot == '10' or error_plot == '50':
    # plus error
        ax2.scatter(pressure_medians_p[0:2], medians_swirl[0:2],s=100,color='plum',marker="^")
        ax2.scatter(pressure_medians_p[2], medians_swirl[2],s=100,color='plum')
        ax2.scatter(pressure_medians_p[3:5], medians_swirl[3:5],s=100,color='plum',marker="^")
        ax2.scatter(pressure_medians_p[5:9], medians_swirl[5:9],s=100,color='plum')
        ax2.scatter(pressure_medians_p[9], medians_swirl[9],s=100,color='plum',marker="v")
        ax2.scatter(pressure_medians_p[10], medians_swirl[10],s=100,color='plum')
        ax2.scatter(pressure_medians_p[11], medians_swirl[11],s=100,color='plum',marker="v")
        ax2.scatter(pressure_medians_p[12:18], medians_swirl[12:18],s=100,color='plum',marker=",")
        ax2.scatter(pressure_medians_p[18:21], medians_swirl[18:21],s=100,color='plum',marker="v")
        ax2.scatter(pressure_medians_p[21], medians_swirl[21],s=100,color='plum',marker="^")
        ax2.scatter(pressure_medians_p[22:], medians_swirl[22:],s=100,color='plum',marker=",")
    
        # minus error
        ax2.scatter(pressure_medians_m[0:2], medians_swirl[0:2],s=100,color='peachpuff',marker="^")
        ax2.scatter(pressure_medians_m[2], medians_swirl[2],s=100,color='peachpuff')
        ax2.scatter(pressure_medians_m[3:5], medians_swirl[3:5],s=100,color='peachpuff',marker="^")
        ax2.scatter(pressure_medians_m[5:9], medians_swirl[5:9],s=100,color='peachpuff')
        ax2.scatter(pressure_medians_m[9], medians_swirl[9],s=100,color='peachpuff',marker="v")
        ax2.scatter(pressure_medians_m[10], medians_swirl[10],s=100,color='peachpuff')
        ax2.scatter(pressure_medians_m[11], medians_swirl[11],s=100,color='peachpuff',marker="v")
        ax2.scatter(pressure_medians_m[12:18], medians_swirl[12:18],s=100,color='peachpuff',marker=",")
        ax2.scatter(pressure_medians_m[18:21], medians_swirl[18:21],s=100,color='peachpuff',marker="v")
        ax2.scatter(pressure_medians_m[21], medians_swirl[21],s=100,color='peachpuff',marker="^")
        ax2.scatter(pressure_medians_m[22:], medians_swirl[22:],s=100,color='peachpuff',marker=",")
        
        ax2.errorbar(pressure_medians_p, medians_swirl, yerr=swirl_med_errs, xerr=pressure_med_errs_p, fmt='.', color='plum')
        
        ax2.plot(x_new, y_new_swirl_p, '--',color='darkviolet', linewidth=3, markersize=12,label=f'+ {error}%: {adjr_sp:.2g}')
    
        ax2.errorbar(pressure_medians_m, medians_swirl, yerr=swirl_med_errs, xerr=pressure_med_errs_m, fmt='.', color='peachpuff')
        
        ax2.plot(x_new, y_new_swirl_m, '--',color='sandybrown', linewidth=3, markersize=12,label=f'- {error}%: {adjr_sm:.2g}')
    
    
    # main plot
    ax2.scatter(Bz_pos_pressure[0:2], Bz_pos_swirlpower[0:2], color='red', s=200,zorder=9)
    ax2.scatter(Bz_pos_pressure[2:], Bz_pos_swirlpower[2:], color='red', s=200, marker='v',zorder=9)

    ax2.scatter(Bz_neg_pressure[0], Bz_neg_swirlpower[0], color='orange', s=200,zorder=9) #
    ax2.scatter(Bz_neg_pressure[1], Bz_neg_swirlpower[1], color='orange', s=200,marker=',',zorder=9) 

    ax2.scatter(By_pos_pressure[0], By_pos_swirlpower[0], color='green', s=200,zorder=9) # compression 1 2-5 incl visit 15
    ax2.scatter(By_pos_pressure[1:4], By_pos_swirlpower[1:4], color='green', s=200, marker=",",zorder=9)
    ax2.scatter(By_pos_pressure[4:6], By_pos_swirlpower[4:6], color='green', s=200, marker="v",zorder=9) # rarefraction
    ax2.scatter(By_pos_pressure[6], By_pos_swirlpower[6], color='green', s=200, marker="^",zorder=9)
    ax2.scatter(By_pos_pressure[7:], By_pos_swirlpower[7:], color='green', s=200, marker=",",zorder=9)

    ax2.scatter(By_neg_pressure[0:4], By_neg_swirlpower[0:4], color='blue', s=200, marker="^",zorder=9)
    ax2.scatter(By_neg_pressure[4:6], By_neg_swirlpower[4:6], color='blue', s=200,zorder=9) # 
    ax2.scatter(By_neg_pressure[6], By_neg_swirlpower[6], color='blue', s=200, marker="v",zorder=9)
    ax2.scatter(By_neg_pressure[7:], By_neg_swirlpower[7:], color='blue', s=200, marker=",",zorder=9)
    
    ax2.text(3.3e-1,6.48,'b',style='italic',fontsize=40) #0.33

    ax2.set_xlabel('Median Pressure (nPa)', fontsize=22)
    ax2.set_ylabel('Median Swirl Region Power (GW)', fontsize=22)
    
    ax2.set_xlim(3.4e-3, 4.5e-1)
    #ax2.set_xlim(0.0, 0.36)
    ax2.set_xscale('log')
    locmaj = matplotlib.ticker.LogLocator(base=10,numticks=30) 
    ax2.xaxis.set_major_locator(locmaj)
    ax2.tick_params(axis='x',which='minor',direction='in')
    ax2.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax2.minorticks_on()
    ax2.tick_params(which='major',direction='in',bottom=True, top=True, left=True, right=True, labelsize=20, length=10)
    ax2.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True, length=5)
    ax2.legend(framealpha=0.5,fontsize=21,loc ="upper left")
    
    
    ax3 = plt.subplot(4,2,3)
    
    #ax3.semilogx(pressure_medians, dusk_medians/(1e12))     
    ax3.plot(x_new, y_new_dusk, '--',color='black',linewidth=3, markersize=12, label=(f'$R^2$: {adjr_d:.2g}'), zorder=10)
    ax3.errorbar(pressure_medians, dusk_medians, yerr=dusk_med_errs, xerr=med_pressure_errs, fmt='.', color='lightgray', zorder=10)
    
    if error_plot == '20' or error_plot == '10' or error_plot == '50':
        # plus error
        ax3.scatter(pressure_medians_p[0:2], dusk_medians[0:2],s=100,color='plum',marker="^")
        ax3.scatter(pressure_medians_p[2], dusk_medians[2],s=100,color='plum')
        ax3.scatter(pressure_medians_p[3:5], dusk_medians[3:5],s=100,color='plum',marker="^")
        ax3.scatter(pressure_medians_p[5:9], dusk_medians[5:9],s=100,color='plum')
        ax3.scatter(pressure_medians_p[9], dusk_medians[9],s=100,color='plum',marker="v")
        ax3.scatter(pressure_medians_p[10], dusk_medians[10],s=100,color='plum')
        ax3.scatter(pressure_medians_p[11], dusk_medians[11],s=100,color='plum',marker="v")
        ax3.scatter(pressure_medians_p[12:18], dusk_medians[12:18],s=100,color='plum',marker=",")
        ax3.scatter(pressure_medians_p[18:21], dusk_medians[18:21],s=100,color='plum',marker="v")
        ax3.scatter(pressure_medians_p[21], dusk_medians[21],s=100,color='plum',marker="^")
        ax3.scatter(pressure_medians_p[22:], dusk_medians[22:],s=100,color='plum',marker=",")
    
        # minus error
        ax3.scatter(pressure_medians_m[0:2], dusk_medians[0:2],s=100,color='peachpuff',marker="^")
        ax3.scatter(pressure_medians_m[2], dusk_medians[2],s=100,color='peachpuff')
        ax3.scatter(pressure_medians_m[3:5], dusk_medians[3:5],s=100,color='peachpuff',marker="^")
        ax3.scatter(pressure_medians_m[5:9], dusk_medians[5:9],s=100,color='peachpuff')
        ax3.scatter(pressure_medians_m[9], dusk_medians[9],s=100,color='peachpuff',marker="v")
        ax3.scatter(pressure_medians_m[10], dusk_medians[10],s=100,color='peachpuff')
        ax3.scatter(pressure_medians_m[11], dusk_medians[11],s=100,color='peachpuff',marker="v")
        ax3.scatter(pressure_medians_m[12:18], dusk_medians[12:18],s=100,color='peachpuff',marker=",")
        ax3.scatter(pressure_medians_m[18:21], dusk_medians[18:21],s=100,color='peachpuff',marker="v")
        ax3.scatter(pressure_medians_m[21], dusk_medians[21],s=100,color='peachpuff',marker="^")
        ax3.scatter(pressure_medians_m[22:], dusk_medians[22:],s=100,color='peachpuff',marker=",")
        
        ax3.errorbar(pressure_medians_p, dusk_medians, yerr=dusk_med_errs, xerr=pressure_med_errs_p, fmt='.', color='plum')
        
        ax3.plot(x_new, y_new_dusk_p, '--',color='darkviolet', linewidth=3, markersize=12,label=f'+ {error}%: {adjr_dp:.2g}')
        
        ax3.errorbar(pressure_medians_m, dusk_medians, yerr=dusk_med_errs, xerr=pressure_med_errs_m, fmt='.', color='peachpuff')
        
        ax3.plot(x_new, y_new_dusk_m, '--',color='sandybrown', linewidth=3, markersize=12,label=f'- {error}%: {adjr_dm:.2g}')
    
    # main plot
    ax3.scatter(Bz_pos_pressure[0:2], Bz_pos_duskpower[0:2], color='red', s=200,zorder=9)
    ax3.scatter(Bz_pos_pressure[2:], Bz_pos_duskpower[2:], color='red', s=200, marker='v',zorder=9)

    ax3.scatter(Bz_neg_pressure[0], Bz_neg_duskpower[0], color='orange', s=200,zorder=9) #
    ax3.scatter(Bz_neg_pressure[1], Bz_neg_duskpower[1], color='orange', s=200,marker=',',zorder=9)

    ax3.scatter(By_pos_pressure[0], By_pos_duskpower[0], color='green', s=200,zorder=9) # compression 1 2-5 incl visit 15
    ax3.scatter(By_pos_pressure[1:4], By_pos_duskpower[1:4], color='green', s=200, marker=",",zorder=9)
    ax3.scatter(By_pos_pressure[4:6], By_pos_duskpower[4:6], color='green', s=200, marker="v",zorder=9) # rarefraction
    ax3.scatter(By_pos_pressure[6], By_pos_duskpower[6], color='green', s=200, marker="^",zorder=9)
    ax3.scatter(By_pos_pressure[7:], By_pos_duskpower[7:], color='green', s=200, marker=",",zorder=9)

    ax3.scatter(By_neg_pressure[0:4], By_neg_duskpower[0:4], color='blue', s=200, marker="^",zorder=9)
    ax3.scatter(By_neg_pressure[4:6], By_neg_duskpower[4:6], color='blue', s=200,zorder=9) # 
    ax3.scatter(By_neg_pressure[6], By_neg_duskpower[6], color='blue', s=200, marker="v",zorder=9)
    ax3.scatter(By_neg_pressure[7:], By_neg_duskpower[7:], color='blue', s=200, marker=",",zorder=9)


    ax3.text(3.3e-1,-13,'c',style='italic',fontsize=40)
    
    ax3.set_xlabel('Median Pressure (nPa)', fontsize=22)
    ax3.set_ylabel('Median Dusk Active Region Power (GW)', fontsize=22)

    ax3.set_xlim(3.4e-3, 4.5e-1)
    #ax3.set_xlim(0.0, 0.36)
    ax3.set_xscale('log')
    locmaj = matplotlib.ticker.LogLocator(base=10,numticks=30) 
    ax3.xaxis.set_major_locator(locmaj)
    ax3.tick_params(axis='x',which='minor',direction='in')
    ax3.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax3.minorticks_on()
    ax3.tick_params(which='major',direction='in',bottom=True, top=True, left=True, right=True, labelsize=20, length=10)
    ax3.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True, length=5)
    ax3.legend(framealpha=0.5,fontsize=21,loc ="upper left")
    
    
    ax4 = plt.subplot(4,2,4)

    #ax4.semilogx(pressure_medians, noon_medians/(1e12)) 
    ax4.plot(x_new, y_new_noon, '--',color='black', linewidth=3, markersize=12,label=(f'$R^2$: {adjr_n:.2g}'), zorder=10)
    
    if error_plot == '20' or error_plot == '10' or error_plot == '50':
        # plus error
        ax4.scatter(pressure_medians_p[0:2], noon_medians[0:2],s=100,color='plum',marker="^")
        ax4.scatter(pressure_medians_p[2], noon_medians[2],s=100,color='plum')
        ax4.scatter(pressure_medians_p[3:5], noon_medians[3:5],s=100,color='plum',marker="^")
        ax4.scatter(pressure_medians_p[5:9], noon_medians[5:9],s=100,color='plum')
        ax4.scatter(pressure_medians_p[9], noon_medians[9],s=100,color='plum',marker="v")
        ax4.scatter(pressure_medians_p[10], noon_medians[10],s=100,color='plum')
        ax4.scatter(pressure_medians_p[11], noon_medians[11],s=100,color='plum',marker="v")
        ax4.scatter(pressure_medians_p[12:18], noon_medians[12:18],s=100,color='plum',marker=",")
        ax4.scatter(pressure_medians_p[18:21], noon_medians[18:21],s=100,color='plum',marker="v")
        ax4.scatter(pressure_medians_p[21], noon_medians[21],s=100,color='plum',marker="^")
        ax4.scatter(pressure_medians_p[22:], noon_medians[22:],s=100,color='plum',marker=",")
    
        # minus error
        ax4.scatter(pressure_medians_m[0:2], noon_medians[0:2],s=100,color='peachpuff',marker="^")
        ax4.scatter(pressure_medians_m[2], noon_medians[2],s=100,color='peachpuff')
        ax4.scatter(pressure_medians_m[3:5], noon_medians[3:5],s=100,color='peachpuff',marker="^")
        ax4.scatter(pressure_medians_m[5:9], noon_medians[5:9],s=100,color='peachpuff')
        ax4.scatter(pressure_medians_m[9], noon_medians[9],s=100,color='peachpuff',marker="v")
        ax4.scatter(pressure_medians_m[10], noon_medians[10],s=100,color='peachpuff')
        ax4.scatter(pressure_medians_m[11], noon_medians[11],s=100,color='peachpuff',marker="v")
        ax4.scatter(pressure_medians_m[12:18], noon_medians[12:18],s=100,color='peachpuff',marker=",")
        ax4.scatter(pressure_medians_m[18:21], noon_medians[18:21],s=100,color='peachpuff',marker="v")
        ax4.scatter(pressure_medians_m[21], noon_medians[21],s=100,color='peachpuff',marker="^")
        ax4.scatter(pressure_medians_m[22:], noon_medians[22:],s=100,color='peachpuff',marker=",")
        
        
        ax4.errorbar(pressure_medians_p, noon_medians, yerr=noon_med_errs, xerr=pressure_med_errs_p, fmt='.', color='plum')
        
        ax4.plot(x_new, y_new_noon_p, '--',color='darkviolet',linewidth=3, markersize=12, label=f'+ {error}%: {adjr_np:.2g}')
        
        ax4.errorbar(pressure_medians_m, noon_medians, yerr=noon_med_errs, xerr=pressure_med_errs_m, fmt='.', color='peachpuff')
        
        ax4.plot(x_new, y_new_noon_m, '--',color='sandybrown', linewidth=3, markersize=12,label=f'- {error}%: {adjr_nm:.2g}')
    
    
    # main plot
    ax4.scatter(Bz_pos_pressure[0:2], Bz_pos_noonpower[0:2], color='red', s=200,zorder=9)
    ax4.scatter(Bz_pos_pressure[2:], Bz_pos_noonpower[2:], color='red', s=200, marker='v',zorder=9)

    ax4.scatter(Bz_neg_pressure[0], Bz_neg_noonpower[0], color='orange', s=200,zorder=9) #
    ax4.scatter(Bz_neg_pressure[1], Bz_neg_noonpower[1], color='orange', s=200,marker=',',zorder=9) 

    ax4.scatter(By_pos_pressure[0], By_pos_noonpower[0], color='green', s=200,zorder=9) # compression 1 2-5 incl visit 15
    ax4.scatter(By_pos_pressure[1:4], By_pos_noonpower[1:4], color='green', s=200, marker=",",zorder=9)
    ax4.scatter(By_pos_pressure[4:6], By_pos_noonpower[4:6], color='green', s=200, marker="v",zorder=9) # rarefraction
    ax4.scatter(By_pos_pressure[6], By_pos_noonpower[6], color='green', s=200, marker="^",zorder=9)
    ax4.scatter(By_pos_pressure[7:], By_pos_noonpower[7:], color='green', s=200, marker=",",zorder=9)

    ax4.scatter(By_neg_pressure[0:4], By_neg_noonpower[0:4], color='blue', s=200, marker="^",zorder=9)
    ax4.scatter(By_neg_pressure[4:6], By_neg_noonpower[4:6], color='blue', s=200,zorder=9)
    ax4.scatter(By_neg_pressure[6], By_neg_noonpower[6], color='blue', s=200, marker="v",zorder=9)
    ax4.scatter(By_neg_pressure[7:], By_neg_noonpower[7:], color='blue', s=200, marker=",",zorder=9)
    
    ax4.errorbar(pressure_medians, noon_medians, yerr=noon_med_errs, xerr=med_pressure_errs, fmt='.', color='lightgray', zorder=9)#
    ax4.text(3.3e-1,-12,'d',style='italic',fontsize=40)

    ax4.set_ylabel('Median Noon Active Region Power (GW)', fontsize=22)
    ax4.set_xlabel('Median Pressure (nPa)', fontsize=22)

    ax4.set_xlim(3.4e-3, 4.5e-1)
    ax4.set_xscale('log')
    locmaj = matplotlib.ticker.LogLocator(base=10,numticks=30) 
    ax4.xaxis.set_major_locator(locmaj)
    ax4.tick_params(axis='x',which='minor',direction='in')
    ax4.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax4.minorticks_on()
    ax4.tick_params(which='major',direction='in',bottom=True, top=True, left=True, right=True, labelsize=20, length=10)
    ax4.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True, length=5)
    ax4.legend(framealpha=0.5,fontsize=21,loc ="upper left")
    
    
    # save plot
    saveloc = (f'{root_saves}median_pressure_vs_region_power_{error_plot}.jpg')
    #saveloc = (f'/Users/hannah/OneDrive - Lancaster University/aurora/median_pressure_vs_region_power_log_{error}.jpg')
    plt.savefig(saveloc,bbox_inches='tight',dpi=400)

  

if plotting == 'LL':  
    
    x_new = np.linspace(min(LL_medians_m/(1e6)), (max(LL_medians/1e6)), len(LL_medians)) #_p
    
    a, acov = np.polyfit(np.log10(LL_medians/(1e6)), medians_swirl, 1, cov=True)
    aa = np.poly1d(a)
    # print("Swirl Fit Gradient:")
    # print(a[0])
    # print("Error of Swirl Fit:")
    # print(np.sqrt(np.diag(acov)))
    y_regress_swirl = a[0]*np.log10(LL_medians/1e6)+a[1]
    
    print(f"Testing '{plotting}'")
    print('Swirl Region Fit')
    ss_res = np.sum((medians_swirl - y_regress_swirl)**2)
    ss_tot = np.sum((medians_swirl - np.mean(medians_swirl))**2)
    r_squared_s = 1 - (ss_res / ss_tot)
    print(f"R squared: {r_squared_s:.4f}")

    adj_s = 1-r_squared_s
    adjr_s = 1 - ((adj_s*23)/22)
    print(f"Adjusted R squared: {adjr_s:.4f}")
    
    
    swirl_grad = round(a[0],2)
    swirl_grad_err = np.sqrt(np.diag(acov))
    swirl_grad_err = round(swirl_grad_err[0], 2)
    
    swirl_LL_fit = (f'{swirl_grad} ± {swirl_grad_err}')
    
    
    a_p,acov_p = np.polyfit(np.log10(LL_medians_p/(1e6)), medians_swirl,1, cov=True)
    aa_p = np.poly1d(a_p)
    # print(f"Swirl Fit +{error} Gradient:")
    # print(a_p[0])
    # print(f"Error of Swirl Fit +{error}:")
    # print(np.sqrt(np.diag(acov_p)))
    y_regress_swirl_p = a_p[0]*np.log10(LL_medians_p/1e6)+a_p[1]
    
    print("---------------")
    print('Swirl Region +20% Fit')
    ssp_res = np.sum((medians_swirl - y_regress_swirl_p)**2)
    ssp_tot = np.sum((medians_swirl - np.mean(medians_swirl))**2)
    r_squared_sp = 1 - (ssp_res / ssp_tot)
    print(f"R squared: {r_squared_sp:.4f}")

    adj_sp = 1-r_squared_sp
    adjr_sp = 1 - ((adj_sp*23)/22)
    print(f"Adjusted R squared: {adjr_sp:.4f}")
    
    
    swirl_grad_p = round(a_p[0],2)
    swirl_grad_err_p = np.sqrt(np.diag(acov_p))
    swirl_grad_err_p = round(swirl_grad_err_p[0], 2)
    
    swirl_LL_fit_p = (f'{swirl_grad_p} ± {swirl_grad_err_p}')
    
    
    a_m,acov_m = np.polyfit(np.log10(LL_medians_m/(1e6)), medians_swirl,1, cov=True)
    aa_m = np.poly1d(a_m)
    # print(f"Swirl Fit -{error} Gradient:")
    # print(a_m[0])
    # print(f"Error of Swirl Fit -{error}:")
    # print(np.sqrt(np.diag(acov_m)))
    y_regress_swirl_m = a_m[0]*np.log10(LL_medians_m/1e6)+a_m[1]
    
    print("---------------")
    print('Swirl Region -20% Fit')
    ssm_res = np.sum((medians_swirl - y_regress_swirl_m)**2)
    ssm_tot = np.sum((medians_swirl - np.mean(medians_swirl))**2)
    r_squared_sm = 1 - (ssm_res / ssm_tot)
    print(f"R squared: {r_squared_sm:.4f}")

    adj_sm = 1-r_squared_sm
    adjr_sm = 1 - ((adj_sm*23)/22)
    print(f"Adjusted R squared: {adjr_sm:.4f}")
    
    
    swirl_grad_m = round(a_m[0],2)
    swirl_grad_err_m = np.sqrt(np.diag(acov_m))
    swirl_grad_err_m = round(swirl_grad_err_m[0], 2)
    
    swirl_LL_fit_m = (f'{swirl_grad_m} ± {swirl_grad_err_m}')
   
    
    
    b, bcov = np.polyfit(np.log10(LL_medians/(1e6)), dusk_medians, 1, cov=True)
    bb = np.poly1d(b)
    # print("Dusk Fit Gradient:")
    # print(b[0])
    # print("Error of Dusk Fit:")
    # print(np.sqrt(np.diag(bcov)))
    y_regress_dusk = b[0]*np.log10(LL_medians/1e6)+b[1]
    
    print("---------------")
    print('Dusk Region Fit')
    sd_res = np.sum((dusk_medians - y_regress_dusk)**2)
    sd_tot = np.sum((dusk_medians - np.mean(dusk_medians))**2)
    r_squared_d = 1 - (sd_res / sd_tot)
    print(f"R squared: {r_squared_d:.4f}")

    adj_d = 1-r_squared_d
    adjr_d = 1 - ((adj_d*23)/22)
    print(f"Adjusted R squared: {adjr_d:.4f}")
    
    dusk_grad = round(b[0],2)
    dusk_grad_err = np.sqrt(np.diag(bcov))
    dusk_grad_err = round(dusk_grad_err[0], 2)
    
    dusk_LL_fit = (f'{dusk_grad} ± {dusk_grad_err}')
    
    
    b_p, bcov_p = np.polyfit(np.log10(LL_medians_p/(1e6)), dusk_medians, 1, cov=True)
    bb_p = np.poly1d(b_p)
    # print(f"Dusk Fit +{error} Gradient:")
    # print(b_p[0])
    # print(f"Error of Dusk Fit +{error}:")
    # print(np.sqrt(np.diag(bcov_p)))
    y_regress_dusk_p = b_p[0]*np.log10(LL_medians_p/1e6)+b_p[1]
    
    print("---------------")
    print('Dusk Region +20% Fit')
    sdp_res = np.sum((dusk_medians - y_regress_dusk_p)**2)
    sdp_tot = np.sum((dusk_medians - np.mean(dusk_medians))**2)
    r_squared_dp = 1 - (sdp_res / sdp_tot)
    print(f"R squared: {r_squared_dp:.4f}")

    adj_dp = 1-r_squared_dp
    adjr_dp = 1 - ((adj_dp*23)/22)
    print(f"Adjusted R squared: {adjr_dp:.4f}")
    
    dusk_grad_p = round(b_p[0],2)
    dusk_grad_err_p = np.sqrt(np.diag(bcov_p))
    dusk_grad_err_p = round(dusk_grad_err_p[0], 2)
    
    dusk_LL_fit_p = (f'{dusk_grad_p} ± {dusk_grad_err_p}')
    
    
    b_m, bcov_m = np.polyfit(np.log10(LL_medians_m/(1e6)), dusk_medians, 1, cov=True)
    bb_m = np.poly1d(b_m)
    # print(f"Dusk Fit -{error} Gradient:")
    # print(b_m[0])
    # print(f"Error of Dusk Fit -{error}:")
    # print(np.sqrt(np.diag(bcov_m)))
    y_regress_dusk_m = b_m[0]*np.log10(LL_medians_m/1e6)+b_m[1]
    
    print("---------------")
    print('Dusk Region -20% Fit')
    sdm_res = np.sum((dusk_medians - y_regress_dusk_m)**2)
    sdm_tot = np.sum((dusk_medians - np.mean(dusk_medians))**2)
    r_squared_dm = 1 - (sdm_res / sdm_tot)
    print(f"R squared: {r_squared_dm:.4f}")

    adj_dm = 1-r_squared_dm
    adjr_dm = 1 - ((adj_dm*23)/22)
    print(f"Adjusted R squared: {adjr_dm:.4f}")
    
    dusk_grad_m = round(b_m[0],2)
    dusk_grad_err_m = np.sqrt(np.diag(bcov_m))
    dusk_grad_err_m = round(dusk_grad_err_m[0], 2)
    
    dusk_LL_fit_m = (f'{dusk_grad_m} ± {dusk_grad_err_m}')
    
    
    
    c, ccov = np.polyfit(np.log10(LL_medians/(1e6)), noon_medians, 1, cov=True)
    cc = np.poly1d(c)
    # print("Noon Fit Gradient:")
    # print(c[0])
    # print("Error of Noon Fit:")
    # print(np.sqrt(np.diag(ccov)))
    y_regress_noon = c[0]*np.log10(LL_medians/1e6)+c[1]
    
    print("---------------")
    print('Noon Region Fit')
    sn_res = np.sum((noon_medians - y_regress_noon)**2)
    sn_tot = np.sum((noon_medians - np.mean(noon_medians))**2)
    r_squared_n = 1 - (sn_res / sn_tot)
    print(f"R squared: {r_squared_n:.4f}")

    adj_n = 1-r_squared_n
    adjr_n = 1 - ((adj_n*23)/22)
    print(f"Adjusted R squared: {adjr_n:.4f}")
    
    
    noon_grad = round(c[0],2)
    noon_grad_err = np.sqrt(np.diag(ccov))
    noon_grad_err = round(noon_grad_err[0], 2)
    
    noon_LL_fit = (f'{noon_grad} ± {noon_grad_err}')
    
    
    c_p, ccov_p = np.polyfit(np.log10(LL_medians_p/(1e6)), noon_medians, 1, cov=True)
    cc_p = np.poly1d(c_p)
    # print(f"Noon Fit +{error} Gradient:")
    # print(c_p[0])
    # print(f"Error of Noon Fit +{error}:")
    # print(np.sqrt(np.diag(ccov_p)))
    y_regress_noon_p = c_p[0]*np.log10(LL_medians_p/1e6)+c_p[1]
    
    print("---------------")
    print('Noon Region +20% Fit')
    snp_res = np.sum((noon_medians - y_regress_noon_p)**2)
    snp_tot = np.sum((noon_medians - np.mean(noon_medians))**2)
    r_squared_np = 1 - (snp_res / snp_tot)
    print(f"R squared: {r_squared_np:.4f}")

    adj_np = 1-r_squared_np
    adjr_np = 1 - ((adj_np*23)/22)
    print(f"Adjusted R squared: {adjr_np:.4f}")
    
    
    noon_grad_p = round(c_p[0],2)
    noon_grad_err_p = np.sqrt(np.diag(ccov_p))
    noon_grad_err_p = round(noon_grad_err_p[0], 2)
    
    noon_LL_fit_p = (f'{noon_grad_p} ± {noon_grad_err_p}')
    
    
    c_m, ccov_m = np.polyfit(np.log10(LL_medians_m/(1e6)), noon_medians, 1, cov=True)
    cc_m = np.poly1d(c_m)
    # print(f"Noon Fit -{error} Gradient:")
    # print(c_m[0])
    # print(f"Error of Noon Fit -{error}:")
    # print(np.sqrt(np.diag(ccov_m)))
    y_regress_noon_m = c_p[0]*np.log10(LL_medians_m/1e6)+c_m[1]
    
    noon_grad_m = round(c_m[0],2)
    noon_grad_err_m = np.sqrt(np.diag(ccov_m))
    noon_grad_err_m = round(noon_grad_err_m[0], 2)
    
    noon_LL_fit_m = (f'{noon_grad_m} ± {noon_grad_err_m}')
    
    print("---------------")
    print('Noon Region -20% Fit')
    snm_res = np.sum((noon_medians - y_regress_noon_m)**2)
    snm_tot = np.sum((noon_medians - np.mean(noon_medians))**2)
    r_squared_nm = 1 - (snm_res / snm_tot)
    print(f"R squared: {r_squared_nm:.4f}")

    adj_nm = 1-r_squared_nm
    adjr_nm = 1 - ((adj_nm*23)/22)
    print(f"Adjusted R squared: {adjr_nm:.4f}")
    
    
    y_new_swirl = a[0]*np.log10(x_new)+a[1]
    y_new_dusk = b[0]*np.log10(x_new)+b[1]
    y_new_noon = c[0]*np.log10(x_new)+c[1]
    
    y_new_swirl_p = a_p[0]*np.log10(x_new)+a_p[1]
    y_new_dusk_p = b_p[0]*np.log10(x_new)+b_p[1]
    y_new_noon_p = c_p[0]*np.log10(x_new)+c_p[1]
    
    y_new_swirl_m = a_m[0]*np.log10(x_new)+a_m[1]
    y_new_dusk_m = b_m[0]*np.log10(x_new)+b_m[1]
    y_new_noon_m = c_m[0]*np.log10(x_new)+c_m[1]
    
    
      
    fig = plt.figure(figsize=(28,50))
    ax1 = plt.subplot(4,2,1)
    plt.subplots_adjust(hspace=0.1, wspace=0.15)
    
    # main plot
    ax1.scatter(Bz_pos_LL[0:2]/(1e6), Bz_pos_power[0:2], color='red', s=200, label='+Bz CME')
    ax1.scatter(Bz_pos_LL[2:]/(1e6), Bz_pos_power[2:], color='red', s=200, marker='v', label='+Bz Rarefaction (Deep)')

    ax1.scatter(Bz_neg_LL[0]/(1e6), Bz_neg_power[0], color='orange', s=200, label='-Bz CME') #
    ax1.scatter(Bz_neg_LL[1]/(1e6), Bz_neg_power[1], color='orange', s=200,marker=',', label='-Bz CIR') 

    ax1.scatter(By_pos_LL[0]/(1e6), By_pos_power[0], color='green', s=200, label='+By CME') # compression 1 2-5 incl visit 15
    ax1.scatter(By_pos_LL[1:4]/(1e6), By_pos_power[1:4], color='green', s=200, marker=",", label='+By CIR')
    ax1.scatter(By_pos_LL[4:6]/(1e6), By_pos_power[4:6], color='green', s=200, marker="v", label='+By Rarefaction (Deep)') # rarefraction
    ax1.scatter(By_pos_LL[6]/(1e6), By_pos_power[6], color='green', s=200, marker="^", label ='+By Rarefaction (Shallow)')
    ax1.scatter(By_pos_LL[7:]/(1e6), By_pos_power[7:], color='green', s=200, marker=",")

    ax1.scatter(By_neg_LL[0:4]/(1e6), By_neg_power[0:4], color='blue', s=200, marker="^", label='-By Rarefaction (Shallow)')
    ax1.scatter(By_neg_LL[4:6]/(1e6), By_neg_power[4:6], color='blue', s=200, label='-By CME') # 
    ax1.scatter(By_neg_LL[6]/(1e6), By_neg_power[6], color='blue', s=200, marker="v", label='-By Rarefaction (Deep)')
    ax1.scatter(By_neg_LL[7:]/(1e6), By_neg_power[7:], color='blue', s=200, marker=",", label='-By CIR')
    
    ax1.errorbar(LL_medians/(1e6), polar_medians, yerr=polar_med_errs, xerr=med_LL_errs/(1e6), fmt='.', color='lightgray')

    # leg1 = ax1.legend(labels=[r'$ϕ_{LL} = V_{sw}B_⊥(\frac{MP_{loc}}{2})cos^4(\frac{θ_c}{2})$           '],fontsize=22,loc='lower center')
    # for handle in leg1.legend_handles:
    #     handle.set_color('white')
    # ax1.add_artist(leg1)
    
    #ax1.text(0.05,1,'a',style='italic',fontsize=30) #1.4 RHS, 0.05 LHS (1 is top)]
    ax1.text(2.2,830,'a',style='italic',fontsize=40) #1.4 RHS, 0.05 LHS (1 is top)
    
    ax1.set_xlabel('Low-Latitude Reconnection Voltage (MV)', fontsize=22)
    ax1.set_ylabel('Median Total Polar Power (GW)', fontsize=22)
    ax1.set_xscale('log')
    locmaj = matplotlib.ticker.LogLocator(base=10,numticks=30) 
    ax1.xaxis.set_major_locator(locmaj)
    ax1.tick_params(axis='x',which='minor',direction='in')
    ax1.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax1.minorticks_on()
    ax1.set_ylim(-5,900)
    ax1.tick_params(which='major',direction='in',bottom=True, top=True, left=True, right=True, labelsize=20, length=10)
    ax1.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True, length=5)
    lgnd = plt.legend(framealpha=0.5, labelspacing=0.3, fontsize=18, loc="upper left") #if using equation legend
    for handle in lgnd.legend_handles:
        handle.set_sizes([100.0])
    
    
    
    ax2 = plt.subplot(4,2,2)
    
    ax2.plot(x_new, y_new_swirl, '--',color='black',linewidth=3, markersize=12, label=(f'Swirl $R^2$: {adjr_s:.2g}'), zorder=10)
    

    if error_plot == '20' or error_plot == '10' or error_plot == '50':
        
        ax2.plot(x_new, y_new_swirl_p, '--',color='darkviolet',linewidth=3, markersize=12, label=f'+ {error}%: {adjr_sp:.2g}')
        ax2.plot(x_new, y_new_swirl_m, '--',color='sandybrown',linewidth=3, markersize=12, label=f'- {error}%: {adjr_sm:.2g}')
        # plus error
        ax2.scatter(LL_medians_p[0:2]/(1e6), medians_swirl[0:2],s=100,color='plum',marker="^")
        ax2.scatter(LL_medians_p[2]/(1e6), medians_swirl[2],s=100,color='plum')
        ax2.scatter(LL_medians_p[3:5]/(1e6), medians_swirl[3:5],s=100,color='plum',marker="^")
        ax2.scatter(LL_medians_p[5:9]/(1e6), medians_swirl[5:9],s=100,color='plum')
        ax2.scatter(LL_medians_p[9]/(1e6), medians_swirl[9],s=100,color='plum',marker="v")
        ax2.scatter(LL_medians_p[10]/(1e6), medians_swirl[10],s=100,color='plum')
        ax2.scatter(LL_medians_p[11]/(1e6), medians_swirl[11],s=100,color='plum')
        ax2.scatter(LL_medians_p[12:18]/(1e6), medians_swirl[12:18],s=100,color='plum',marker=",")
        ax2.scatter(LL_medians_p[18:21]/(1e6), medians_swirl[18:21],s=100,color='plum',marker="v")
        ax2.scatter(LL_medians_p[21]/(1e6), medians_swirl[21],s=100,color='plum',marker="^")
        ax2.scatter(LL_medians_p[22:]/(1e6), medians_swirl[22:],s=100,color='plum',marker=",")
    
        # minus error
        ax2.scatter(LL_medians_m[0:2]/(1e6), medians_swirl[0:2],s=100,color='peachpuff',marker="^")
        ax2.scatter(LL_medians_m[2]/(1e6), medians_swirl[2],s=100,color='peachpuff')
        ax2.scatter(LL_medians_m[3:5]/(1e6), medians_swirl[3:5],s=100,color='peachpuff',marker="^")
        ax2.scatter(LL_medians_m[5:9]/(1e6), medians_swirl[5:9],s=100,color='peachpuff')
        ax2.scatter(LL_medians_m[9]/(1e6), medians_swirl[9],s=100,color='peachpuff',marker="v")
        ax2.scatter(LL_medians_m[10]/(1e6), medians_swirl[10],s=100,color='peachpuff')
        ax2.scatter(LL_medians_m[11]/(1e6), medians_swirl[11],s=100,color='peachpuff')
        ax2.scatter(LL_medians_m[12:18]/(1e6), medians_swirl[12:18],s=100,color='peachpuff',marker=",")
        ax2.scatter(LL_medians_m[18:21]/(1e6), medians_swirl[18:21],s=100,color='peachpuff',marker="v")
        ax2.scatter(LL_medians_m[21]/(1e6), medians_swirl[21],s=100,color='peachpuff',marker="^")
        ax2.scatter(LL_medians_m[22:]/(1e6), medians_swirl[22:],s=100,color='peachpuff',marker=",")
        
        ax2.errorbar(LL_medians_p/(1e6), medians_swirl, yerr=swirl_med_errs, xerr=LL_med_errs_p/(1e6), fmt='.', color='plum')
        
        ax2.errorbar(LL_medians_m/(1e6), medians_swirl, yerr=swirl_med_errs, xerr=LL_med_errs_m/(1e6), fmt='.', color='peachpuff')
    
    # main plot
    ax2.scatter(Bz_pos_LL[0:2]/(1e6), Bz_pos_swirlpower[0:2], color='red', s=200,zorder=9)
    ax2.scatter(Bz_pos_LL[2:]/(1e6), Bz_pos_swirlpower[2:], color='red', s=200, marker='v',zorder=9)

    ax2.scatter(Bz_neg_LL[0]/(1e6), Bz_neg_swirlpower[0], color='orange', s=200,zorder=9) #
    ax2.scatter(Bz_neg_LL[1]/(1e6), Bz_neg_swirlpower[1], color='orange', s=200,marker=',',zorder=9) 

    ax2.scatter(By_pos_LL[0]/(1e6), By_pos_swirlpower[0], color='green', s=200,zorder=9) # compression 1 2-5 incl visit 15
    ax2.scatter(By_pos_LL[1:4]/(1e6), By_pos_swirlpower[1:4], color='green', s=200, marker=",",zorder=9)
    ax2.scatter(By_pos_LL[4:6]/(1e6), By_pos_swirlpower[4:6], color='green', s=200, marker="v",zorder=9) # rarefraction
    ax2.scatter(By_pos_LL[6]/(1e6), By_pos_swirlpower[6], color='green', s=200, marker="^",zorder=9)
    ax2.scatter(By_pos_LL[7:]/(1e6), By_pos_swirlpower[7:], color='green', s=200, marker=",",zorder=9)

    ax2.scatter(By_neg_LL[0:4]/(1e6), By_neg_swirlpower[0:4], color='blue', s=200, marker="^",zorder=9)
    ax2.scatter(By_neg_LL[4:6]/(1e6), By_neg_swirlpower[4:6], color='blue', s=200,zorder=9) # 
    ax2.scatter(By_neg_LL[6]/(1e6), By_neg_swirlpower[6], color='blue', s=200, marker="v",zorder=9)
    ax2.scatter(By_neg_LL[7:]/(1e6), By_neg_swirlpower[7:], color='blue', s=200, marker=",",zorder=9)
    
    ax2.errorbar(LL_medians/(1e6), medians_swirl, yerr=swirl_med_errs, xerr=med_LL_errs/(1e6), fmt='.', color='lightgray', zorder=9)
    
    ax2.text(2.2,76,'b',style='italic',fontsize=40) #69.5
    
    ax2.set_xlabel('Low-Latitude Reconnection Voltage (MV)', fontsize=22)
    ax2.set_ylabel('Median Swirl Region Power (GW)', fontsize=22)
    ax2.set_xscale('log')
    ax2.set_ylim(-5,83)
    locmaj = matplotlib.ticker.LogLocator(base=10,numticks=30) 
    ax2.xaxis.set_major_locator(locmaj)
    ax2.tick_params(axis='x',which='minor',direction='in')
    ax2.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax2.minorticks_on()
    ax2.tick_params(which='major',direction='in',bottom=True, top=True, left=True, right=True, labelsize=20, length=10)
    ax2.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True, length=5)
    ax2.legend(framealpha=0.5,fontsize=21,loc ="upper left")
    
    
    ax3 = plt.subplot(4,2,3)
    
    ax3.plot(x_new, y_new_dusk,'--', color='black',linewidth=3, markersize=12, label=f'Dusk $R^2$: {adjr_d:.2g}', zorder=10)

    if error_plot == '20' or error_plot == '10' or error_plot == '50':
        ax3.plot(x_new, y_new_dusk_p, '--',color='darkviolet',linewidth=3, markersize=12, label=f'+ {error}%: {adjr_dp:.2g}')
        ax3.plot(x_new, y_new_dusk_m, '--',color='sandybrown',linewidth=3, markersize=12, label=f'- {error}%: {adjr_dm:.2g}')
        
        # plus error
        ax3.scatter(LL_medians_p[0:2]/(1e6), dusk_medians[0:2],s=100,color='plum',marker="^")
        ax3.scatter(LL_medians_p[2]/(1e6), dusk_medians[2],s=100,color='plum')
        ax3.scatter(LL_medians_p[3:5]/(1e6), dusk_medians[3:5],s=100,color='plum',marker="^")
        ax3.scatter(LL_medians_p[5:9]/(1e6), dusk_medians[5:9],s=100,color='plum')
        ax3.scatter(LL_medians_p[9]/(1e6), dusk_medians[9],s=100,color='plum',marker="v")
        ax3.scatter(LL_medians_p[10]/(1e6), dusk_medians[10],s=100,color='plum')
        ax3.scatter(LL_medians_p[11]/(1e6), dusk_medians[11],s=100,color='plum')
        ax3.scatter(LL_medians_p[12:18]/(1e6), dusk_medians[12:18],s=100,color='plum',marker=",")
        ax3.scatter(LL_medians_p[18:21]/(1e6), dusk_medians[18:21],s=100,color='plum',marker="v")
        ax3.scatter(LL_medians_p[21]/(1e6), dusk_medians[21],s=100,color='plum',marker="^")
        ax3.scatter(LL_medians_p[22:]/(1e6), dusk_medians[22:],s=100,color='plum',marker=",")
    
        # minus error
        ax3.scatter(LL_medians_m[0:2]/(1e6), dusk_medians[0:2],s=100,color='peachpuff',marker="^")
        ax3.scatter(LL_medians_m[2]/(1e6), dusk_medians[2],s=100,color='peachpuff')
        ax3.scatter(LL_medians_m[3:5]/(1e6), dusk_medians[3:5],s=100,color='peachpuff',marker="^")
        ax3.scatter(LL_medians_m[5:9]/(1e6), dusk_medians[5:9],s=100,color='peachpuff')
        ax3.scatter(LL_medians_m[9]/(1e6), dusk_medians[9],s=100,color='peachpuff',marker="v")
        ax3.scatter(LL_medians_m[10]/(1e6), dusk_medians[10],s=100,color='peachpuff')
        ax3.scatter(LL_medians_m[11]/(1e6), dusk_medians[11],s=100,color='peachpuff')
        ax3.scatter(LL_medians_m[12:18]/(1e6), dusk_medians[12:18],s=100,color='peachpuff',marker=",")
        ax3.scatter(LL_medians_m[18:21]/(1e6), dusk_medians[18:21],s=100,color='peachpuff',marker="v")
        ax3.scatter(LL_medians_m[21]/(1e6), dusk_medians[21],s=100,color='peachpuff',marker="^")
        ax3.scatter(LL_medians_m[22:]/(1e6), dusk_medians[22:],s=100,color='peachpuff',marker=",")
        
        ax3.errorbar(LL_medians_p/(1e6), dusk_medians, yerr=dusk_med_errs, xerr=LL_med_errs_p/(1e6), fmt='.', color='plum')
        
        ax3.errorbar(LL_medians_m/(1e6), dusk_medians, yerr=dusk_med_errs, xerr=LL_med_errs_m/(1e6), fmt='.', color='peachpuff')
    
    
    # main plot
    ax3.scatter(Bz_pos_LL[0:2]/(1e6), Bz_pos_duskpower[0:2], color='red', s=200,zorder=9)
    ax3.scatter(Bz_pos_LL[2:]/(1e6), Bz_pos_duskpower[2:], color='red', s=200, marker='v',zorder=9)

    ax3.scatter(Bz_neg_LL[0]/(1e6), Bz_neg_duskpower[0], color='orange', s=200,zorder=9) #
    ax3.scatter(Bz_neg_LL[1]/(1e6), Bz_neg_duskpower[1], color='orange', s=200,marker=',',zorder=9)

    ax3.scatter(By_pos_LL[0]/(1e6), By_pos_duskpower[0], color='green', s=200,zorder=9) # compression 1 2-5 incl visit 15
    ax3.scatter(By_pos_LL[1:4]/(1e6), By_pos_duskpower[1:4], color='green', s=200, marker=",",zorder=9)
    ax3.scatter(By_pos_LL[4:6]/(1e6), By_pos_duskpower[4:6], color='green', s=200, marker="v",zorder=9) # rarefraction
    ax3.scatter(By_pos_LL[6]/(1e6), By_pos_duskpower[6], color='green', s=200, marker="^",zorder=9)
    ax3.scatter(By_pos_LL[7:]/(1e6), By_pos_duskpower[7:], color='green', s=200, marker=",",zorder=9)

    ax3.scatter(By_neg_LL[0:4]/(1e6), By_neg_duskpower[0:4], color='blue', s=200, marker="^",zorder=9)
    ax3.scatter(By_neg_LL[4:6]/(1e6), By_neg_duskpower[4:6], color='blue', s=200,zorder=9) # 
    ax3.scatter(By_neg_LL[6]/(1e6), By_neg_duskpower[6], color='blue', s=200, marker="v",zorder=9)
    ax3.scatter(By_neg_LL[7:]/(1e6), By_neg_duskpower[7:], color='blue', s=200, marker=",",zorder=9)
    
    ax3.errorbar(LL_medians/(1e6), dusk_medians, yerr=dusk_med_errs, xerr=med_LL_errs/(1e6), fmt='.', color='lightgray', zorder=9)
    
    ax3.text(2.2,399,'c',style='italic',fontsize=40) #1.4 RHS, 0.05 LHS (1 is top)
    
    ax3.set_xlabel('Low-Latitude Reconnection Voltage (MV)', fontsize=22)
    ax3.set_ylabel('Median Dusk Region Power (GW)', fontsize=22)
    ax3.set_xscale('log')
    locmaj = matplotlib.ticker.LogLocator(base=10,numticks=30) 
    ax3.xaxis.set_major_locator(locmaj)
    ax3.tick_params(axis='x',which='minor',direction='in')
    ax3.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax3.minorticks_on()
    ax3.tick_params(which='major',direction='in',bottom=True, top=True, left=True, right=True, labelsize=20, length=10)
    ax3.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True, length=5)
    ax3.legend(framealpha=0.5,fontsize=21,loc ="upper left")
    
    
    ax4 = plt.subplot(4,2,4)
    
    ax4.plot(x_new, y_new_noon, '--',color='black', linewidth=3, markersize=12, label=f'Noon $R^2$: {adjr_n:.2g}', zorder=10)
    
    ax4.plot(x_new, y_new_noon_p, '--',color='darkviolet', linewidth=3, markersize=12,label=f'+ {error}%: {adjr_np:.2g}')
    ax4.plot(x_new, y_new_noon_m, '--',color='sandybrown', linewidth=3, markersize=12, label=f'- {error}%: {adjr_nm:.2g}')
    
    # plus error
    ax4.scatter(LL_medians_p[0:2]/(1e6), noon_medians[0:2],s=100,color='plum',marker="^")
    ax4.scatter(LL_medians_p[2]/(1e6), noon_medians[2],s=100,color='plum')
    ax4.scatter(LL_medians_p[3:5]/(1e6), noon_medians[3:5],s=100,color='plum',marker="^")
    ax4.scatter(LL_medians_p[5:9]/(1e6), noon_medians[5:9],s=100,color='plum')
    ax4.scatter(LL_medians_p[9]/(1e6), noon_medians[9],s=100,color='plum',marker="v")
    ax4.scatter(LL_medians_p[10]/(1e6), noon_medians[10],s=100,color='plum')
    ax4.scatter(LL_medians_p[11]/(1e6), noon_medians[11],s=100,color='plum')
    ax4.scatter(LL_medians_p[12:18]/(1e6), noon_medians[12:18],s=100,color='plum',marker=",")
    ax4.scatter(LL_medians_p[18:21]/(1e6), noon_medians[18:21],s=100,color='plum',marker="v")
    ax4.scatter(LL_medians_p[21]/(1e6), noon_medians[21],s=100,color='plum',marker="^")
    ax4.scatter(LL_medians_p[22:]/(1e6), noon_medians[22:],s=100,color='plum',marker=",")

    # minus error
    ax4.scatter(LL_medians_m[0:2]/(1e6), noon_medians[0:2],s=100,color='peachpuff',marker="^")
    ax4.scatter(LL_medians_m[2]/(1e6), noon_medians[2],s=100,color='peachpuff')
    ax4.scatter(LL_medians_m[3:5]/(1e6), noon_medians[3:5],s=100,color='peachpuff',marker="^")
    ax4.scatter(LL_medians_m[5:9]/(1e6), noon_medians[5:9],s=100,color='peachpuff')
    ax4.scatter(LL_medians_m[9]/(1e6), noon_medians[9],s=100,color='peachpuff',marker="v")
    ax4.scatter(LL_medians_m[10]/(1e6), noon_medians[10],s=100,color='peachpuff')
    ax4.scatter(LL_medians_m[11]/(1e6), noon_medians[11],s=100,color='peachpuff')
    ax4.scatter(LL_medians_m[12:18]/(1e6), noon_medians[12:18],s=100,color='peachpuff',marker=",")
    ax4.scatter(LL_medians_m[18:21]/(1e6), noon_medians[18:21],s=100,color='peachpuff',marker="v")
    ax4.scatter(LL_medians_m[21]/(1e6), noon_medians[21],s=100,color='peachpuff',marker="^")
    ax4.scatter(LL_medians_m[22:]/(1e6), noon_medians[22:],s=100,color='peachpuff',marker=",")
    
    ax4.errorbar(LL_medians_p/(1e6), noon_medians, yerr=noon_med_errs, xerr=LL_med_errs_p/(1e6), fmt='.', color='plum')
    
    ax4.errorbar(LL_medians_m/(1e6), noon_medians, yerr=noon_med_errs, xerr=LL_med_errs_m/(1e6), fmt='.', color='peachpuff')
    
    
    # main plot
    ax4.scatter(Bz_pos_LL[0:2]/(1e6), Bz_pos_noonpower[0:2], color='red', s=200,zorder=9)
    ax4.scatter(Bz_pos_LL[2:]/(1e6), Bz_pos_noonpower[2:], color='red', s=200, marker='v',zorder=9)

    ax4.scatter(Bz_neg_LL[0]/(1e6), Bz_neg_noonpower[0], color='orange', s=200,zorder=9) #
    ax4.scatter(Bz_neg_LL[1]/(1e6), Bz_neg_noonpower[1], color='orange', s=200,marker=',',zorder=9) 

    ax4.scatter(By_pos_LL[0]/(1e6), By_pos_noonpower[0], color='green', s=200,zorder=9) # compression 1 2-5 incl visit 15
    ax4.scatter(By_pos_LL[1:4]/(1e6), By_pos_noonpower[1:4], color='green', s=200, marker=",",zorder=9)
    ax4.scatter(By_pos_LL[4:6]/(1e6), By_pos_noonpower[4:6], color='green', s=200, marker="v",zorder=9) # rarefraction
    ax4.scatter(By_pos_LL[6]/(1e6), By_pos_noonpower[6], color='green', s=200, marker="^",zorder=9)
    ax4.scatter(By_pos_LL[7:]/(1e6), By_pos_noonpower[7:], color='green', s=200, marker=",",zorder=9)

    ax4.scatter(By_neg_LL[0:4]/(1e6), By_neg_noonpower[0:4], color='blue', s=200, marker="^",zorder=9)
    ax4.scatter(By_neg_LL[4:6]/(1e6), By_neg_noonpower[4:6], color='blue', s=200,zorder=9)
    ax4.scatter(By_neg_LL[6]/(1e6), By_neg_noonpower[6], color='blue', s=200, marker="v",zorder=9)
    ax4.scatter(By_neg_LL[7:]/(1e6), By_neg_noonpower[7:], color='blue', s=200, marker=",",zorder=9)
    
    ax4.errorbar(LL_medians/(1e6), noon_medians, yerr=noon_med_errs, xerr=med_LL_errs/(1e6), fmt='.', color='lightgray', zorder=9)
    
    ax4.text(2.2,403,'d',style='italic',fontsize=40) #405
    
    ax4.set_xlabel('Low-Latitude Reconnection Voltage (MV)', fontsize=22)
    ax4.set_ylabel('Median Noon Region Power (GW)', fontsize=22)
    ax4.set_xscale('log')
    locmaj = matplotlib.ticker.LogLocator(base=10,numticks=30) 
    ax4.xaxis.set_major_locator(locmaj)
    ax4.tick_params(axis='x',which='minor',direction='in')
    ax4.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax4.minorticks_on()
    ax4.tick_params(which='major',direction='in',bottom=True, top=True, left=True, right=True, labelsize=20, length=10)
    ax4.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True, length=5)
    ax4.legend(framealpha=0.5,fontsize=21,loc ="upper left")
    
    
    # save plot
    saveloc = (f'{root_saves}median_LL_rec_V_vs_region_power_{error_plot}.jpg')
    #saveloc = (f'/Users/hannah/OneDrive - Lancaster University/aurora/median_LL_rec_V_vs_region_power_{error}.jpg')
    plt.savefig(saveloc,bbox_inches='tight',dpi=400)
    
if plotting == 'HL_neg':
    
     
    x_new = np.linspace(min(HL_medians_neg/(1e6)), max(HL_medians_m_neg/(1e6)), len(HL_medians_neg))
    
    a, acov = np.polyfit(np.log10(HL_medians_neg/(1e6)), medians_swirl, 1, cov=True)
    aa = np.poly1d(a)
    # print("Swirl Fit Gradient:")
    # print(a[0])
    # print("Error of Swirl Fit:")
    # print(np.sqrt(np.diag(acov)))
    y_regress_swirl = a[0]*np.log10(HL_medians_neg/1e6)+a[1]
    
    print(f"Testing '{plotting}'")
    print('Swirl Region Fit')
    ss_res = np.sum((medians_swirl - y_regress_swirl)**2)
    ss_tot = np.sum((medians_swirl - np.mean(medians_swirl))**2)
    r_squared_s = 1 - (ss_res / ss_tot)
    print(f"R squared: {r_squared_s:.4f}")

    adj_s = 1-r_squared_s
    adjr_s = 1 - ((adj_s*23)/22)
    print(f"Adjusted R squared: {adjr_s:.4f}")
    
    swirl_grad = round(a[0],2)
    swirl_grad_err = np.sqrt(np.diag(acov))
    swirl_grad_err = round(swirl_grad_err[0], 2)
    
    swirl_HL_fit = (f'{swirl_grad} ± {swirl_grad_err}')
    
    
    a_p,acov_p = np.polyfit(np.log10(HL_medians_p_neg/(1e6)), medians_swirl,1, cov=True)
    aa_p = np.poly1d(a_p)
    # print(f"Swirl Fit +{error} Gradient:")
    # print(a_p[0])
    # print(f"Error of Swirl Fit +{error}:")
    # print(np.sqrt(np.diag(acov_p)))
    y_regress_swirl_p = a_p[0]*np.log10(HL_medians_p_neg/1e6)+a_p[1]
    
    print("---------------")
    print('Swirl Region +20% Fit')
    ssp_res = np.sum((medians_swirl - y_regress_swirl_p)**2)
    ssp_tot = np.sum((medians_swirl - np.mean(medians_swirl))**2)
    r_squared_sp = 1 - (ssp_res / ssp_tot)
    print(f"R squared: {r_squared_sp:.4f}")

    adj_sp = 1-r_squared_sp
    adjr_sp = 1 - ((adj_sp*23)/22)
    print(f"Adjusted R squared: {adjr_sp:.4f}")
    
    swirl_grad_p = round(a_p[0],2)
    swirl_grad_err_p = np.sqrt(np.diag(acov_p))
    swirl_grad_err_p = round(swirl_grad_err_p[0], 2)
    
    swirl_HL_fit_p = (f'{swirl_grad_p} ± {swirl_grad_err_p}')
    
    
    a_m,acov_m = np.polyfit(np.log10(HL_medians_m_neg/(1e6)), medians_swirl,1, cov=True)
    # aa_m = np.poly1d(a_m)
    # print(f"Swirl Fit -{error} Gradient:")
    # print(a_m[0])
    # print(f"Error of Swirl Fit -{error}:")
    # print(np.sqrt(np.diag(acov_m)))
    y_regress_swirl_m = a_m[0]*np.log10(HL_medians_m_neg/1e6)+a_m[1]
    
    print("---------------")
    print('Swirl Region -20% Fit')
    ssm_res = np.sum((medians_swirl - y_regress_swirl_m)**2)
    ssm_tot = np.sum((medians_swirl - np.mean(medians_swirl))**2)
    r_squared_sm = 1 - (ssm_res / ssm_tot)
    print(f"R squared: {r_squared_sm:.4f}")

    adj_sm = 1-r_squared_sm
    adjr_sm = 1 - ((adj_sm*23)/22)
    print(f"Adjusted R squared: {adjr_sm:.4f}")
    
    
    swirl_grad_m = round(a_m[0],2)
    swirl_grad_err_m = np.sqrt(np.diag(acov_m))
    swirl_grad_err_m = round(swirl_grad_err_m[0], 2)
    
    swirl_HL_fit_m = (f'{swirl_grad_m} ± {swirl_grad_err_m}')
   
    
    b, bcov = np.polyfit(np.log10(HL_medians_neg/(1e6)), dusk_medians, 1, cov=True)
    bb = np.poly1d(b)
    # print("Dusk Fit Gradient:")
    # print(b[0])
    # print("Error of Dusk Fit:")
    # print(np.sqrt(np.diag(bcov)))
    y_regress_dusk = b[0]*np.log10(HL_medians_neg/1e6)+b[1]
    
    print("---------------")
    print('Dusk Region Fit')
    sd_res = np.sum((dusk_medians - y_regress_dusk)**2)
    sd_tot = np.sum((dusk_medians - np.mean(dusk_medians))**2)
    r_squared_d = 1 - (sd_res / sd_tot)
    print(f"R squared: {r_squared_d:.4f}")

    adj_d = 1-r_squared_d
    adjr_d = 1 - ((adj_d*23)/22)
    print(f"Adjusted R squared: {adjr_d:.4f}")
    
    
    dusk_grad = round(b[0],2)
    dusk_grad_err = np.sqrt(np.diag(bcov))
    dusk_grad_err = round(dusk_grad_err[0], 2)
    
    dusk_HL_fit = (f'{dusk_grad} ± {dusk_grad_err}')
    
    
    b_p, bcov_p = np.polyfit(np.log10(HL_medians_p_neg/(1e6)), dusk_medians, 1, cov=True)
    bb_p = np.poly1d(b_p)
    # print(f"Dusk Fit +{error} Gradient:")
    # print(b_p[0])
    # print(f"Error of Dusk Fit +{error}:")
    # print(np.sqrt(np.diag(bcov_p)))
    y_regress_dusk_p = b_p[0]*np.log10(HL_medians_p_neg/1e6)+b_p[1]
    
    print("---------------")
    print('Dusk Region +20% Fit')
    sdp_res = np.sum((dusk_medians - y_regress_dusk_p)**2)
    sdp_tot = np.sum((dusk_medians - np.mean(dusk_medians))**2)
    r_squared_dp = 1 - (sdp_res / sdp_tot)
    print(f"R squared: {r_squared_dp:.4f}")

    adj_dp = 1-r_squared_dp
    adjr_dp = 1 - ((adj_dp*23)/22)
    print(f"Adjusted R squared: {adjr_dp:.4f}")
    
    
    dusk_grad_p = round(b_p[0],2)
    dusk_grad_err_p = np.sqrt(np.diag(bcov_p))
    dusk_grad_err_p = round(dusk_grad_err_p[0], 2)
    
    dusk_HL_fit_p = (f'{dusk_grad_p} ± {dusk_grad_err_p}')
    
    
    b_m, bcov_m = np.polyfit(np.log10(HL_medians_m_neg/(1e6)), dusk_medians, 1, cov=True)
    bb_m = np.poly1d(b_m)
    # print(f"Dusk Fit -{error} Gradient:")
    # print(b_m[0])
    # print(f"Error of Dusk Fit -{error}:")
    # print(np.sqrt(np.diag(bcov_m)))
    y_regress_dusk_m = b_m[0]*np.log10(HL_medians_m_neg/1e6)+b_m[1]
    
    print("---------------")
    print('Dusk Region -20% Fit')
    sdm_res = np.sum((dusk_medians - y_regress_dusk_m)**2)
    sdm_tot = np.sum((dusk_medians - np.mean(dusk_medians))**2)
    r_squared_dm = 1 - (sdm_res / sdm_tot)
    print(f"R squared: {r_squared_dm:.4f}")

    adj_dm = 1-r_squared_dm
    adjr_dm = 1 - ((adj_dm*23)/22)
    print(f"Adjusted R squared: {adjr_dm:.4f}")
    
    dusk_grad_m = round(b_m[0],2)
    dusk_grad_err_m = np.sqrt(np.diag(bcov_m))
    dusk_grad_err_m = round(dusk_grad_err_m[0], 2)
    
    dusk_HL_fit_m = (f'{dusk_grad_m} ± {dusk_grad_err_m}')
    
    
    
    c, ccov = np.polyfit(np.log10(HL_medians_neg/(1e6)), noon_medians, 1, cov=True)
    cc = np.poly1d(c)
    # print("Noon Fit Gradient:")
    # print(c[0])
    # print("Error of Noon Fit:")
    # print(np.sqrt(np.diag(ccov)))
    y_regress_noon = c[0]*np.log10(HL_medians_neg/1e6)+c[1]
    
    print("---------------")
    print('Noon Region Fit')
    sn_res = np.sum((noon_medians - y_regress_noon)**2)
    sn_tot = np.sum((noon_medians - np.mean(noon_medians))**2)
    r_squared_n = 1 - (sn_res / sn_tot)
    print(f"R squared: {r_squared_n:.4f}")

    adj_n = 1-r_squared_n
    adjr_n = 1 - ((adj_n*23)/22)
    print(f"Adjusted R squared: {adjr_n:.4f}")
    
    noon_grad = round(c[0],2)
    noon_grad_err = np.sqrt(np.diag(ccov))
    noon_grad_err = round(noon_grad_err[0], 2)
    
    noon_HL_fit = (f'{noon_grad} ± {noon_grad_err}')
    
    
    c_p, ccov_p = np.polyfit(np.log10(HL_medians_p_neg/(1e6)), noon_medians, 1, cov=True)
    cc_p = np.poly1d(c_p)
    # print(f"Noon Fit +{error} Gradient:")
    # print(c_p[0])
    # print(f"Error of Noon Fit +{error}:")
    # print(np.sqrt(np.diag(ccov_p)))
    y_regress_noon_p = c_p[0]*np.log10(HL_medians_p_neg/1e6)+c_p[1]
    
    print("---------------")
    print('Noon Region +20% Fit')
    snp_res = np.sum((noon_medians - y_regress_noon_p)**2)
    snp_tot = np.sum((noon_medians - np.mean(noon_medians))**2)
    r_squared_np = 1 - (snp_res / snp_tot)
    print(f"R squared: {r_squared_np:.4f}")

    adj_np = 1-r_squared_np
    adjr_np = 1 - ((adj_np*23)/22)
    print(f"Adjusted R squared: {adjr_np:.4f}")
    
    noon_grad_p = round(c_p[0],2)
    noon_grad_err_p = np.sqrt(np.diag(ccov_p))
    noon_grad_err_p = round(noon_grad_err_p[0], 2)
    
    noon_HL_fit_p = (f'{noon_grad_p} ± {noon_grad_err_p}')
    
    
    c_m, ccov_m = np.polyfit(np.log10(HL_medians_m_neg/(1e6)), noon_medians, 1, cov=True)
    cc_m = np.poly1d(c_m)
    # print(f"Noon Fit -{error} Gradient:")
    # print(c_m[0])
    # print(f"Error of Noon Fit -{error}:")
    # print(np.sqrt(np.diag(ccov_m)))
    y_regress_noon_m = c_m[0]*np.log10(HL_medians_m_neg/1e6)+c_m[1]
    
    print("---------------")
    print('Noon Region -20% Fit')
    snm_res = np.sum((noon_medians - y_regress_noon_m)**2)
    snm_tot = np.sum((noon_medians - np.mean(noon_medians))**2)
    r_squared_nm = 1 - (snm_res / snm_tot)
    print(f"R squared: {r_squared_nm:.4f}")

    adj_nm = 1-r_squared_nm
    adjr_nm = 1 - ((adj_nm*23)/22)
    print(f"Adjusted R squared: {adjr_nm:.4f}")
    
    
    noon_grad_m = round(c_m[0],2)
    noon_grad_err_m = np.sqrt(np.diag(ccov_m))
    noon_grad_err_m = round(noon_grad_err_m[0], 2)
    
    noon_HL_fit_m = (f'{noon_grad_m} ± {noon_grad_err_m}')

    
    y_new_swirl = a[0]*np.log10(x_new)+a[1]
    y_new_dusk = b[0]*np.log10(x_new)+b[1]
    y_new_noon = c[0]*np.log10(x_new)+c[1]
    
    y_new_swirl_p = a_p[0]*np.log10(x_new)+a_p[1]
    y_new_dusk_p = b_p[0]*np.log10(x_new)+b_p[1]
    y_new_noon_p = c_p[0]*np.log10(x_new)+c_p[1]
    
    y_new_swirl_m = a_m[0]*np.log10(x_new)+a_m[1]
    y_new_dusk_m = b_m[0]*np.log10(x_new)+b_m[1]
    y_new_noon_m = c_m[0]*np.log10(x_new)+c_m[1]
    
    
    fig = plt.figure(figsize=(28,50))
    ax1 = plt.subplot(4,2,1)
    plt.subplots_adjust(hspace=0.1, wspace=0.15)
    
    
    # main plot
    ax1.scatter(Bz_pos_HL_neg[0:2]/(1e6), Bz_pos_power[0:2], color='red', s=200, label='+Bz CME')
    ax1.scatter(Bz_pos_HL_neg[2:]/(1e6), Bz_pos_power[2:], color='red', s=200, marker='v', label='+Bz Rarefaction (Deep)')

    ax1.scatter(Bz_neg_HL_neg[0]/(1e6), Bz_neg_power[0], color='orange', s=200, label='-Bz CME') #
    ax1.scatter(Bz_neg_HL_neg[1]/(1e6), Bz_neg_power[1], color='orange', s=200,marker=',', label='-Bz CIR') 

    ax1.scatter(By_pos_HL_neg[0]/(1e6), By_pos_power[0], color='green', s=200, label='+By CME') # compression 1 2-5 incl visit 15
    ax1.scatter(By_pos_HL_neg[1:4]/(1e6), By_pos_power[1:4], color='green', s=200, marker=",", label='+By CIR')
    ax1.scatter(By_pos_HL_neg[4:6]/(1e6), By_pos_power[4:6], color='green', s=200, marker="v", label='+By Rarefaction (Deep)') # rarefraction
    ax1.scatter(By_pos_HL_neg[6]/(1e6), By_pos_power[6], color='green', s=200, marker="^", label ='+By Rarefaction (Shallow)')
    ax1.scatter(By_pos_HL_neg[7:]/(1e6), By_pos_power[7:], color='green', s=200, marker=",")

    ax1.scatter(By_neg_HL_neg[0:4]/(1e6), By_neg_power[0:4], color='blue', s=200, marker="^", label='-By Rarefaction (Shallow)')
    ax1.scatter(By_neg_HL_neg[4:6]/(1e6), By_neg_power[4:6], color='blue', s=200, label='-By CME') # 
    ax1.scatter(By_neg_HL_neg[6]/(1e6), By_neg_power[6], color='blue', s=200, marker="v", label='-By Rarefaction (Deep)')
    ax1.scatter(By_neg_HL_neg[7:]/(1e6), By_neg_power[7:], color='blue', s=200, marker=",", label='-By CIR')
    
    ax1.errorbar(HL_medians_neg/(1e6), polar_medians, yerr=polar_med_errs, xerr=med_HL_errs_neg/(1e6), fmt='.', color='lightgray')
    
    # leg1 = ax1.legend(labels=[r'$ϕ_{HL,-B_y} = \frac{1}{2}V_{sw}B_⊥(\frac{MP_{loc}}{2})cos^4(\frac{θ_c + 90}{2})$'],fontsize=22,loc='lower center')
    # for handle in leg1.legend_handles:
    #     handle.set_color('white')
    # #ax1.add_artist(leg1)
    
    #ax1.text(0.6,0.1, r'$ϕ_{HL,-B_y} = \frac{1}{2}V_{sw}B_⊥(\frac{MP_{loc}}{2})cos^4(\frac{θ_c + 90}{2})$',fontsize=16)
    # try this instead 
    ax1.text(3.7e-1,15,'a',style='italic',fontsize=40) #1.4 RHS, 0.05 LHS (1 is top)
    
    ax1.set_xlabel('High-Latitude Dayside $(-B_y)$ Reconnection Voltage', fontsize=22)
    ax1.set_ylabel('Median Total Polar Power (GW)', fontsize=22)
    
    ax1.set_xlim(8e-7, 0.85e0)
    ax1.set_ylim(-25,1140)
    
    ax1.set_xscale('log')
    locmaj = matplotlib.ticker.LogLocator(base=10,numticks=30) 
    ax1.xaxis.set_major_locator(locmaj)
    ax1.tick_params(axis='x',which='minor',direction='in')
    ax1.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax1.minorticks_on()
    ax1.tick_params(which='major',direction='in',bottom=True, top=True, left=True, right=True, labelsize=20, length=10)
    ax1.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True, length=5)
    
    lgnd = plt.legend(framealpha=0.5, labelspacing=0.3, ncols=2, fontsize=18, loc="upper right")
    for handle in lgnd.legend_handles:
        handle.set_sizes([100.0])
    
    
    ax2 = plt.subplot(4,2,2)
    
    ax2.plot(x_new, y_new_swirl, '--',color='black', linewidth=3, markersize=12,label=(f'Swirl $R^2$: {adjr_s:.2g}'), zorder=10)
    
    if error_plot == '20' or error_plot == '10' or error_plot == '50':
        ax2.plot(x_new, y_new_swirl_p, '--',color='darkviolet', linewidth=3, markersize=12,label=f'+ {error}%: {adjr_sp:.2g}')
        ax2.plot(x_new, y_new_swirl_m, '--',color='sandybrown', linewidth=3, markersize=12,label=f'- {error}%: {adjr_sm:.2g}')
        
        # plus error
        ax2.scatter(HL_medians_p_neg[0:2]/(1e6), medians_swirl[0:2],s=100,color='plum',marker="^")
        ax2.scatter(HL_medians_p_neg[2]/(1e6), medians_swirl[2],s=100,color='plum')
        ax2.scatter(HL_medians_p_neg[3:5]/(1e6), medians_swirl[3:5],s=100,color='plum',marker="^")
        ax2.scatter(HL_medians_p_neg[5:9]/(1e6), medians_swirl[5:9],s=100,color='plum')
        ax2.scatter(HL_medians_p_neg[9]/(1e6), medians_swirl[9],s=100,color='plum',marker="v")
        ax2.scatter(HL_medians_p_neg[10]/(1e6), medians_swirl[10],s=100,color='plum')
        ax2.scatter(HL_medians_p_neg[11]/(1e6), medians_swirl[11],s=100,color='plum',marker="v")
        ax2.scatter(HL_medians_p_neg[12:18]/(1e6), medians_swirl[12:18],s=100,color='plum',marker=",")
        ax2.scatter(HL_medians_p_neg[18:21]/(1e6), medians_swirl[18:21],s=100,color='plum',marker="v")
        ax2.scatter(HL_medians_p_neg[21]/(1e6), medians_swirl[21],s=100,color='plum',marker="^")
        ax2.scatter(HL_medians_p_neg[22:]/(1e6), medians_swirl[22:],s=100,color='plum',marker=",")
    
        # minus error
        ax2.scatter(HL_medians_m_neg[0:2]/(1e6), medians_swirl[0:2],s=100,color='peachpuff',marker="^")
        ax2.scatter(HL_medians_m_neg[2]/(1e6), medians_swirl[2],s=100,color='peachpuff')
        ax2.scatter(HL_medians_m_neg[3:5]/(1e6), medians_swirl[3:5],s=100,color='peachpuff',marker="^")
        ax2.scatter(HL_medians_m_neg[5:9]/(1e6), medians_swirl[5:9],s=100,color='peachpuff')
        ax2.scatter(HL_medians_m_neg[9]/(1e6), medians_swirl[9],s=100,color='peachpuff',marker="v")
        ax2.scatter(HL_medians_m_neg[10]/(1e6), medians_swirl[10],s=100,color='peachpuff')
        ax2.scatter(HL_medians_m_neg[11]/(1e6), medians_swirl[11],s=100,color='peachpuff',marker="v")
        ax2.scatter(HL_medians_m_neg[12:18]/(1e6), medians_swirl[12:18],s=100,color='peachpuff',marker=",")
        ax2.scatter(HL_medians_m_neg[18:21]/(1e6), medians_swirl[18:21],s=100,color='peachpuff',marker="v")
        ax2.scatter(HL_medians_m_neg[21]/(1e6), medians_swirl[21],s=100,color='peachpuff',marker="^")
        ax2.scatter(HL_medians_m_neg[22:]/(1e6), medians_swirl[22:],s=100,color='peachpuff',marker=",")
        
        ax2.errorbar(HL_medians_p_neg/(1e6), medians_swirl, yerr=swirl_med_errs, xerr=HL_med_errs_p_neg/(1e6), fmt='.', color='plum')
        
        ax2.errorbar(HL_medians_m_neg/(1e6), medians_swirl, yerr=swirl_med_errs, xerr=HL_med_errs_m_neg/(1e6), fmt='.', color='peachpuff')
    
    
    # main plot
    ax2.scatter(Bz_pos_HL_neg[0:2]/(1e6), Bz_pos_swirlpower[0:2], color='red', s=200,zorder=9)
    ax2.scatter(Bz_pos_HL_neg[2:]/(1e6), Bz_pos_swirlpower[2:], color='red', s=200, marker='v',zorder=9)

    ax2.scatter(Bz_neg_HL_neg[0]/(1e6), Bz_neg_swirlpower[0], color='orange', s=200,zorder=9) #
    ax2.scatter(Bz_neg_HL_neg[1]/(1e6), Bz_neg_swirlpower[1], color='orange', s=200,marker=',',zorder=9) 

    ax2.scatter(By_pos_HL_neg[0]/(1e6), By_pos_swirlpower[0], color='green', s=200,zorder=9) # compression 1 2-5 incl visit 15
    ax2.scatter(By_pos_HL_neg[1:4]/(1e6), By_pos_swirlpower[1:4], color='green', s=200, marker=",",zorder=9)
    ax2.scatter(By_pos_HL_neg[4:6]/(1e6), By_pos_swirlpower[4:6], color='green', s=200, marker="v",zorder=9) # rarefraction
    ax2.scatter(By_pos_HL_neg[6]/(1e6), By_pos_swirlpower[6], color='green', s=200, marker="^",zorder=9)
    ax2.scatter(By_pos_HL_neg[7:]/(1e6), By_pos_swirlpower[7:], color='green', s=200, marker=",",zorder=9)

    ax2.scatter(By_neg_HL_neg[0:4]/(1e6), By_neg_swirlpower[0:4], color='blue', s=200, marker="^",zorder=9)
    ax2.scatter(By_neg_HL_neg[4:6]/(1e6), By_neg_swirlpower[4:6], color='blue', s=200,zorder=9) # 
    ax2.scatter(By_neg_HL_neg[6]/(1e6), By_neg_swirlpower[6], color='blue', s=200, marker="v",zorder=9)
    ax2.scatter(By_neg_HL_neg[7:]/(1e6), By_neg_swirlpower[7:], color='blue', s=200, marker=",",zorder=9)
    
    ax2.errorbar(HL_medians_neg/(1e6), medians_swirl, yerr=swirl_med_errs, xerr=med_HL_errs_neg/(1e6), fmt='.', color='lightgray', zorder=9)
    
    ax2.text(3.4e-1,0.9,'b',style='italic',fontsize=40) #3.5, 3.8 / 0.4
    
    ax2.set_xlabel('High-Latitude Dayside $(-B_y)$ Reconnection Voltage (MV)', fontsize=22)
    ax2.set_ylabel('Median Swirl Region Power (GW)', fontsize=22)
    
    ax2.set_xlim(4e-7, 0.85e0)
    ax2.set_ylim(-2,87) # or no
    ax2.set_xscale('log')
    locmaj = matplotlib.ticker.LogLocator(base=10,numticks=30) 
    ax2.xaxis.set_major_locator(locmaj)
    ax2.tick_params(axis='x',which='minor',direction='in')
    ax2.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax2.minorticks_on()
    ax2.tick_params(which='major',direction='in',bottom=True, top=True, left=True, right=True, labelsize=20, length=10)
    ax2.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True, length=5)

    ax2.legend(framealpha=0.5,fontsize=21,loc ="upper right")
    
    
    ax3 = plt.subplot(4,2,3)
    
    ax3.plot(x_new, y_new_dusk,'--', color='black', linewidth=3, markersize=12,label=f'Dusk $R^2$: {adjr_d:.2g}', zorder=10)
   
    if error_plot == '20' or error_plot == '10' or error_plot == '50':
        ax3.plot(x_new, y_new_dusk_p, '--',color='darkviolet',linewidth=3, markersize=12, label=f'+ {error}%: {adjr_dp:.2g}')
        ax3.plot(x_new, y_new_dusk_m, '--',color='sandybrown', linewidth=3, markersize=12,label=f'- {error}%: {adjr_dm:.2g}')
        
        # plus error
        ax3.scatter(HL_medians_p_neg[0:2]/(1e6), dusk_medians[0:2],s=100,color='plum',marker="^")
        ax3.scatter(HL_medians_p_neg[2]/(1e6), dusk_medians[2],s=100,color='plum')
        ax3.scatter(HL_medians_p_neg[3:5]/(1e6), dusk_medians[3:5],s=100,color='plum',marker="^")
        ax3.scatter(HL_medians_p_neg[5:9]/(1e6), dusk_medians[5:9],s=100,color='plum')
        ax3.scatter(HL_medians_p_neg[9]/(1e6), dusk_medians[9],s=100,color='plum',marker="v")
        ax3.scatter(HL_medians_p_neg[10]/(1e6), dusk_medians[10],s=100,color='plum')
        ax3.scatter(HL_medians_p_neg[11]/(1e6), dusk_medians[11],s=100,color='plum',marker="v")
        ax3.scatter(HL_medians_p_neg[12:18]/(1e6), dusk_medians[12:18],s=100,color='plum',marker=",")
        ax3.scatter(HL_medians_p_neg[18:21]/(1e6), dusk_medians[18:21],s=100,color='plum',marker="v")
        ax3.scatter(HL_medians_p_neg[21]/(1e6), dusk_medians[21],s=100,color='plum',marker="^")
        ax3.scatter(HL_medians_p_neg[22:]/(1e6), dusk_medians[22:],s=100,color='plum',marker=",")
    
        # minus error
        ax3.scatter(HL_medians_m_neg[0:2]/(1e6), dusk_medians[0:2],s=100,color='peachpuff',marker="^")
        ax3.scatter(HL_medians_m_neg[2]/(1e6), dusk_medians[2],s=100,color='peachpuff')
        ax3.scatter(HL_medians_m_neg[3:5]/(1e6), dusk_medians[3:5],s=100,color='peachpuff',marker="^")
        ax3.scatter(HL_medians_m_neg[5:9]/(1e6), dusk_medians[5:9],s=100,color='peachpuff')
        ax3.scatter(HL_medians_m_neg[9]/(1e6), dusk_medians[9],s=100,color='peachpuff',marker="v")
        ax3.scatter(HL_medians_m_neg[10]/(1e6), dusk_medians[10],s=100,color='peachpuff')
        ax3.scatter(HL_medians_m_neg[11]/(1e6), dusk_medians[11],s=100,color='peachpuff',marker="v")
        ax3.scatter(HL_medians_m_neg[12:18]/(1e6), dusk_medians[12:18],s=100,color='peachpuff',marker=",")
        ax3.scatter(HL_medians_m_neg[18:21]/(1e6), dusk_medians[18:21],s=100,color='peachpuff',marker="v")
        ax3.scatter(HL_medians_m_neg[21]/(1e6), dusk_medians[21],s=100,color='peachpuff',marker="^")
        ax3.scatter(HL_medians_m_neg[22:]/(1e6), dusk_medians[22:],s=100,color='peachpuff',marker=",")
        
        ax3.errorbar(HL_medians_p_neg/(1e6), dusk_medians, yerr=dusk_med_errs, xerr=HL_med_errs_p_neg/(1e6), fmt='.', color='plum')
        
        ax3.errorbar(HL_medians_m_neg/(1e6), dusk_medians, yerr=dusk_med_errs, xerr=HL_med_errs_m_neg/(1e6), fmt='.', color='peachpuff')
    
    
    # main plot
    ax3.scatter(Bz_pos_HL_neg[0:2]/(1e6), Bz_pos_duskpower[0:2], color='red', s=200,zorder=9)
    ax3.scatter(Bz_pos_HL_neg[2:]/(1e6), Bz_pos_duskpower[2:], color='red', s=200, marker='v',zorder=9)

    ax3.scatter(Bz_neg_HL_neg[0]/(1e6), Bz_neg_duskpower[0], color='orange', s=200,zorder=9) #
    ax3.scatter(Bz_neg_HL_neg[1]/(1e6), Bz_neg_duskpower[1], color='orange', s=200,marker=',',zorder=9)

    ax3.scatter(By_pos_HL_neg[0]/(1e6), By_pos_duskpower[0], color='green', s=200,zorder=9) # compression 1 2-5 incl visit 15
    ax3.scatter(By_pos_HL_neg[1:4]/(1e6), By_pos_duskpower[1:4], color='green', s=200, marker=",",zorder=9)
    ax3.scatter(By_pos_HL_neg[4:6]/(1e6), By_pos_duskpower[4:6], color='green', s=200, marker="v",zorder=9) # rarefraction
    ax3.scatter(By_pos_HL_neg[6]/(1e6), By_pos_duskpower[6], color='green', s=200, marker="^",zorder=9)
    ax3.scatter(By_pos_HL_neg[7:]/(1e6), By_pos_duskpower[7:], color='green', s=200, marker=",",zorder=9)

    ax3.scatter(By_neg_HL_neg[0:4]/(1e6), By_neg_duskpower[0:4], color='blue', s=200, marker="^",zorder=9)
    ax3.scatter(By_neg_HL_neg[4:6]/(1e6), By_neg_duskpower[4:6], color='blue', s=200,zorder=9) # 
    ax3.scatter(By_neg_HL_neg[6]/(1e6), By_neg_duskpower[6], color='blue', s=200, marker="v",zorder=9)
    ax3.scatter(By_neg_HL_neg[7:]/(1e6), By_neg_duskpower[7:], color='blue', s=200, marker=",",zorder=9)
    
    ax3.errorbar(HL_medians_neg/(1e6), dusk_medians, yerr=dusk_med_errs, xerr=med_HL_errs_neg/(1e6), fmt='.', color='lightgray', zorder=9)
    
    ax3.text(3.4e-1,-14.5,'c',style='italic',fontsize=40) #1.4 RHS, 0.05 LHS (1 is top)
    
    ax3.set_xlabel('High-Latitude Dayside $(-B_y)$ Reconnection Voltage (MV)', fontsize=22)
    ax3.set_ylabel('Median Dusk Region Power (GW)', fontsize=22)

    ax3.set_xlim(4e-7, 0.85e0)
    #ax3.set_ylim(-20,480)
    
    ax3.set_xscale('log')
    locmaj = matplotlib.ticker.LogLocator(base=10,numticks=30) 
    ax3.xaxis.set_major_locator(locmaj)
    ax3.tick_params(axis='x',which='minor',direction='in')
    ax3.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax3.minorticks_on()
    ax3.tick_params(which='major',direction='in',bottom=True, top=True, left=True, right=True, labelsize=20, length=10)
    ax3.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True, length=5)
    ax3.legend(framealpha=0.5,fontsize=21,loc ="upper right")
    
    
    
    ax4 = plt.subplot(4,2,4)
    
    ax4.plot(x_new, y_new_noon, '--',color='black',linewidth=3, markersize=12, label=f'Noon $R^2$: {adjr_n:.2g}', zorder=10)
    
    if error_plot == '20' or error_plot == '10' or error_plot == '50':
        ax4.plot(x_new, y_new_noon_p, '--',color='darkviolet',linewidth=3, markersize=12, label=f'+ {error}%: {adjr_np:.2g}')
        ax4.plot(x_new, y_new_noon_m, '--',color='sandybrown',linewidth=3, markersize=12, label=f'- {error}%: {adjr_nm:.2g}')
        
        # plus error
        ax4.scatter(HL_medians_p_neg[0:2]/(1e6), noon_medians[0:2],s=100,color='plum',marker="^")
        ax4.scatter(HL_medians_p_neg[2]/(1e6), noon_medians[2],s=100,color='plum')
        ax4.scatter(HL_medians_p_neg[3:5]/(1e6), noon_medians[3:5],s=100,color='plum',marker="^")
        ax4.scatter(HL_medians_p_neg[5:9]/(1e6), noon_medians[5:9],s=100,color='plum')
        ax4.scatter(HL_medians_p_neg[9]/(1e6), noon_medians[9],s=100,color='plum',marker="v")
        ax4.scatter(HL_medians_p_neg[10]/(1e6), noon_medians[10],s=100,color='plum')
        ax4.scatter(HL_medians_p_neg[11]/(1e6), noon_medians[11],s=100,color='plum',marker="v")
        ax4.scatter(HL_medians_p_neg[12:18]/(1e6), noon_medians[12:18],s=100,color='plum',marker=",")
        ax4.scatter(HL_medians_p_neg[18:21]/(1e6), noon_medians[18:21],s=100,color='plum',marker="v")
        ax4.scatter(HL_medians_p_neg[21]/(1e6), noon_medians[21],s=100,color='plum',marker="^")
        ax4.scatter(HL_medians_p_neg[22:]/(1e6), noon_medians[22:],s=100,color='plum',marker=",")
    
        # minus error
        ax4.scatter(HL_medians_m_neg[0:2]/(1e6), noon_medians[0:2],s=100,color='peachpuff',marker="^")
        ax4.scatter(HL_medians_m_neg[2]/(1e6), noon_medians[2],s=100,color='peachpuff')
        ax4.scatter(HL_medians_m_neg[3:5]/(1e6), noon_medians[3:5],s=100,color='peachpuff',marker="^")
        ax4.scatter(HL_medians_m_neg[5:9]/(1e6), noon_medians[5:9],s=100,color='peachpuff')
        ax4.scatter(HL_medians_m_neg[9]/(1e6), noon_medians[9],s=100,color='peachpuff',marker="v")
        ax4.scatter(HL_medians_m_neg[10]/(1e6), noon_medians[10],s=100,color='peachpuff')
        ax4.scatter(HL_medians_m_neg[11]/(1e6), noon_medians[11],s=100,color='peachpuff',marker="v")
        ax4.scatter(HL_medians_m_neg[12:18]/(1e6), noon_medians[12:18],s=100,color='peachpuff',marker=",")
        ax4.scatter(HL_medians_m_neg[18:21]/(1e6), noon_medians[18:21],s=100,color='peachpuff',marker="v")
        ax4.scatter(HL_medians_m_neg[21]/(1e6), noon_medians[21],s=100,color='peachpuff',marker="^")
        ax4.scatter(HL_medians_m_neg[22:]/(1e6), noon_medians[22:],s=100,color='peachpuff',marker=",")
        
        ax4.errorbar(HL_medians_p_neg/(1e6), noon_medians, yerr=noon_med_errs, xerr=HL_med_errs_p_neg/(1e6), fmt='.', color='plum')
        
        ax4.errorbar(HL_medians_m_neg/(1e6), noon_medians, yerr=noon_med_errs, xerr=HL_med_errs_m_neg/(1e6), fmt='.', color='peachpuff')
    
    # main plot
    ax4.scatter(Bz_pos_HL_neg[0:2]/(1e6), Bz_pos_noonpower[0:2], color='red', s=200,zorder=9)
    ax4.scatter(Bz_pos_HL_neg[2:]/(1e6), Bz_pos_noonpower[2:], color='red', s=200, marker='v',zorder=9)

    ax4.scatter(Bz_neg_HL_neg[0]/(1e6), Bz_neg_noonpower[0], color='orange', s=200,zorder=9) #
    ax4.scatter(Bz_neg_HL_neg[1]/(1e6), Bz_neg_noonpower[1], color='orange', s=200,marker=',',zorder=9) 

    ax4.scatter(By_pos_HL_neg[0]/(1e6), By_pos_noonpower[0], color='green', s=200,zorder=9) # compression 1 2-5 incl visit 15
    ax4.scatter(By_pos_HL_neg[1:4]/(1e6), By_pos_noonpower[1:4], color='green', s=200, marker=",",zorder=9)
    ax4.scatter(By_pos_HL_neg[4:6]/(1e6), By_pos_noonpower[4:6], color='green', s=200, marker="v",zorder=9) # rarefraction
    ax4.scatter(By_pos_HL_neg[6]/(1e6), By_pos_noonpower[6], color='green', s=200, marker="^",zorder=9)
    ax4.scatter(By_pos_HL_neg[7:]/(1e6), By_pos_noonpower[7:], color='green', s=200, marker=",",zorder=9)

    ax4.scatter(By_neg_HL_neg[0:4]/(1e6), By_neg_noonpower[0:4], color='blue', s=200, marker="^",zorder=9)
    ax4.scatter(By_neg_HL_neg[4:6]/(1e6), By_neg_noonpower[4:6], color='blue', s=200,zorder=9)
    ax4.scatter(By_neg_HL_neg[6]/(1e6), By_neg_noonpower[6], color='blue', s=200, marker="v",zorder=9)
    ax4.scatter(By_neg_HL_neg[7:]/(1e6), By_neg_noonpower[7:], color='blue', s=200, marker=",",zorder=9)
    
    ax4.errorbar(HL_medians_neg/(1e6), noon_medians, yerr=noon_med_errs, xerr=med_HL_errs_neg/(1e6), fmt='.', color='lightgray', zorder=9)
    
    ax4.set_xlabel('High-Latitude Dayside $(-B_y)$ Reconnection Voltage (MV)', fontsize=22)
    ax4.set_ylabel('Median Noon Region Power (GW)', fontsize=22)
    
    ax4.set_xscale('log')
    locmaj = matplotlib.ticker.LogLocator(base=10,numticks=30) 
    ax4.xaxis.set_major_locator(locmaj)
    ax4.tick_params(axis='x',which='minor',direction='in')
    ax4.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax4.minorticks_on()
    ax4.tick_params(which='major',direction='in',bottom=True, top=True, left=True, right=True, labelsize=20, length=10)
    ax4.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True, length=5)
    ax4.legend(framealpha=0.5,fontsize=21,loc ="upper right")
    
    ax4.text(3.4e-1,-13,'d',style='italic',fontsize=40) #1.4 RHS, 0.05 LHS (1 is top)
    ax4.set_xlim(4e-7, 0.85e0) #9e-8 e7
    #ax4.set_ylim(-20,510)
    # save plot
    saveloc = (f'{root_saves}median_HL_neg_rec_V_vs_region_power_{error_plot}.jpg')
    #saveloc = (f'/Users/hannah/OneDrive - Lancaster University/aurora/median_HL_neg_rec_V_vs_region_power_{error}.jpg')
    plt.savefig(saveloc,bbox_inches='tight',dpi=400)
    

if plotting == 'HL_pos':
    
     
    x_new = np.linspace(min(HL_medians_p_pos/(1e6)), max(HL_medians_p_pos/(1e6)), len(HL_medians_pos)) 
    
    a, acov = np.polyfit(np.log10(HL_medians_pos/(1e6)), medians_swirl, 1, cov=True)
    aa = np.poly1d(a)
    # print("Swirl Fit Gradient:")
    # print(a[0])
    # print("Error of Swirl Fit:")
    # print(np.sqrt(np.diag(acov)))
    y_regress_swirl = a[0]*np.log10(HL_medians_pos/1e6)+a[1]
    
    print(f"Testing '{plotting}'")
    print('Swirl Region Fit')
    ss_res = np.sum((medians_swirl - y_regress_swirl)**2)
    ss_tot = np.sum((medians_swirl - np.mean(medians_swirl))**2)
    r_squared_s = 1 - (ss_res / ss_tot)
    print(f"R squared: {r_squared_s:.4f}")

    adj_s = 1-r_squared_s
    adjr_s = 1 - ((adj_s*23)/22)
    print(f"Adjusted R squared: {adjr_s:.4f}")
    
    swirl_grad = round(a[0],2)
    swirl_grad_err = np.sqrt(np.diag(acov))
    swirl_grad_err = round(swirl_grad_err[0], 2)
    
    swirl_HL_fit = (f'{swirl_grad} ± {swirl_grad_err}')
    
    
    a_p,acov_p = np.polyfit(np.log10(HL_medians_p_pos/(1e6)), medians_swirl,1, cov=True)
    # aa_p = np.poly1d(a_p)
    # print(f"Swirl Fit +{error} Gradient:")
    # print(a_p[0])
    # print(f"Error of Swirl Fit +{error}:")
    # print(np.sqrt(np.diag(acov_p)))
    y_regress_swirl_p = a_p[0]*np.log10(HL_medians_p_pos/1e6)+a_p[1]
    
    print("---------------")
    print('Swirl Region +20% Fit')
    ssp_res = np.sum((medians_swirl - y_regress_swirl_p)**2)
    ssp_tot = np.sum((medians_swirl - np.mean(medians_swirl))**2)
    r_squared_sp = 1 - (ssp_res / ssp_tot)
    print(f"R squared: {r_squared_sp:.4f}")

    adj_sp = 1-r_squared_sp
    adjr_sp = 1 - ((adj_sp*23)/22)
    print(f"Adjusted R squared: {adjr_sp:.4f}")
    
    swirl_grad_p = round(a_p[0],2)
    swirl_grad_err_p = np.sqrt(np.diag(acov_p))
    swirl_grad_err_p = round(swirl_grad_err_p[0], 2)
    
    swirl_HL_fit_p = (f'{swirl_grad_p} ± {swirl_grad_err_p}')
    
    
    a_m,acov_m = np.polyfit(np.log10(HL_medians_m_pos/(1e6)), medians_swirl,1, cov=True)
    aa_m = np.poly1d(a_m)
    # print(f"Swirl Fit -{error} Gradient:")
    # print(a_m[0])
    # print(f"Error of Swirl Fit -{error}:")
    # print(np.sqrt(np.diag(acov_m)))
    y_regress_swirl_m = a_m[0]*np.log10(HL_medians_m_pos/1e6)+a_m[1]
    
    print("---------------")
    print('Swirl Region -20% Fit')
    ssm_res = np.sum((medians_swirl - y_regress_swirl_m)**2)
    ssm_tot = np.sum((medians_swirl - np.mean(medians_swirl))**2)
    r_squared_sm = 1 - (ssm_res / ssm_tot)
    print(f"R squared: {r_squared_sm:.4f}")

    adj_sm = 1-r_squared_sm
    adjr_sm = 1 - ((adj_sm*23)/22)
    print(f"Adjusted R squared: {adjr_sm:.4f}")
    
    swirl_grad_m = round(a_m[0],2)
    swirl_grad_err_m = np.sqrt(np.diag(acov_m))
    swirl_grad_err_m = round(swirl_grad_err_m[0], 2)
    
    swirl_HL_fit_m = (f'{swirl_grad_m} ± {swirl_grad_err_m}')
   
    
    b, bcov = np.polyfit(np.log10(HL_medians_pos/(1e6)), dusk_medians, 1, cov=True)
    bb = np.poly1d(b)
    # print("Dusk Fit Gradient:")
    # print(b[0])
    # print("Error of Dusk Fit:")
    # print(np.sqrt(np.diag(bcov)))
    y_regress_dusk = b[0]*np.log10(HL_medians_pos/1e6)+b[1]
    
    print("---------------")
    print('Dusk Region Fit')
    sd_res = np.sum((dusk_medians - y_regress_dusk)**2)
    sd_tot = np.sum((dusk_medians - np.mean(dusk_medians))**2)
    r_squared_d = 1 - (sd_res / sd_tot)
    print(f"R squared: {r_squared_d:.4f}")

    adj_d = 1-r_squared_d
    adjr_d = 1 - ((adj_d*23)/22)
    print(f"Adjusted R squared: {adjr_d:.4f}")
    
    dusk_grad = round(b[0],2)
    dusk_grad_err = np.sqrt(np.diag(bcov))
    dusk_grad_err = round(dusk_grad_err[0], 2)
    
    dusk_HL_fit = (f'{dusk_grad} ± {dusk_grad_err}')
    
    
    b_p, bcov_p = np.polyfit(np.log10(HL_medians_p_pos/(1e6)), dusk_medians, 1, cov=True)
    bb_p = np.poly1d(b_p)
    # print(f"Dusk Fit +{error} Gradient:")
    # print(b_p[0])
    # print(f"Error of Dusk Fit +{error}:")
    # print(np.sqrt(np.diag(bcov_p)))
    y_regress_dusk_p = b_p[0]*np.log10(HL_medians_p_pos/1e6)+b_p[1]
    
    print("---------------")
    print('Dusk Region +20% Fit')
    sdp_res = np.sum((dusk_medians - y_regress_dusk_p)**2)
    sdp_tot = np.sum((dusk_medians - np.mean(dusk_medians))**2)
    r_squared_dp = 1 - (sdp_res / sdp_tot)
    print(f"R squared: {r_squared_dp:.4f}")

    adj_dp = 1-r_squared_dp
    adjr_dp = 1 - ((adj_dp*23)/22)
    print(f"Adjusted R squared: {adjr_dp:.4f}")
    
    dusk_grad_p = round(b_p[0],2)
    dusk_grad_err_p = np.sqrt(np.diag(bcov_p))
    dusk_grad_err_p = round(dusk_grad_err_p[0], 2)
    
    dusk_HL_fit_p = (f'{dusk_grad_p} ± {dusk_grad_err_p}')
    
    
    b_m, bcov_m = np.polyfit(np.log10(HL_medians_m_pos/(1e6)), dusk_medians, 1, cov=True)
    bb_m = np.poly1d(b_m)
    # print(f"Dusk Fit -{error} Gradient:")
    # print(b_m[0])
    # print(f"Error of Dusk Fit -{error}:")
    # print(np.sqrt(np.diag(bcov_m)))
    y_regress_dusk_m = b_m[0]*np.log10(HL_medians_m_pos/1e6)+b_m[1]
    
    print("---------------")
    print('Dusk Region -20% Fit')
    sdm_res = np.sum((dusk_medians - y_regress_dusk_m)**2)
    sdm_tot = np.sum((dusk_medians - np.mean(dusk_medians))**2)
    r_squared_dm = 1 - (sdm_res / sdm_tot)
    print(f"R squared: {r_squared_dm:.4f}")

    adj_dm = 1-r_squared_dm
    adjr_dm = 1 - ((adj_dm*23)/22)
    print(f"Adjusted R squared: {adjr_dm:.4f}")
    
    dusk_grad_m = round(b_m[0],2)
    dusk_grad_err_m = np.sqrt(np.diag(bcov_m))
    dusk_grad_err_m = round(dusk_grad_err_m[0], 2)
    
    dusk_HL_fit_m = (f'{dusk_grad_m} ± {dusk_grad_err_m}')
    
    
    
    c, ccov = np.polyfit(np.log10(HL_medians_pos/(1e6)), noon_medians, 1, cov=True)
    cc = np.poly1d(c)
    # print("Noon Fit Gradient:")
    # print(c[0])
    # print("Error of Noon Fit:")
    # print(np.sqrt(np.diag(ccov)))
    y_regress_noon = c[0]*np.log10(HL_medians_pos/1e6)+c[1]
    
    print("---------------")
    print('Noon Region Fit')
    sn_res = np.sum((noon_medians - y_regress_noon)**2)
    sn_tot = np.sum((noon_medians - np.mean(noon_medians))**2)
    r_squared_n = 1 - (sn_res / sn_tot)
    print(f"R squared: {r_squared_n:.4f}")

    adj_n = 1-r_squared_n
    adjr_n = 1 - ((adj_n*23)/22)
    print(f"Adjusted R squared: {adjr_n:.4f}")
    
    noon_grad = round(c[0],2)
    noon_grad_err = np.sqrt(np.diag(ccov))
    noon_grad_err = round(noon_grad_err[0], 2)
    
    noon_HL_fit = (f'{noon_grad} ± {noon_grad_err}')
    
    
    c_p, ccov_p = np.polyfit(np.log10(HL_medians_p_pos/(1e6)), noon_medians, 1, cov=True)
    cc_p = np.poly1d(c_p)
    # print(f"Noon Fit +{error} Gradient:")
    # print(c_p[0])
    # print(f"Error of Noon Fit +{error}:")
    # print(np.sqrt(np.diag(ccov_p)))
    y_regress_noon_p = c_p[0]*np.log10(HL_medians_p_pos/1e6)+c_p[1]
    
    print("---------------")
    print('Noon Region +20% Fit')
    snp_res = np.sum((noon_medians - y_regress_noon_p)**2)
    snp_tot = np.sum((noon_medians - np.mean(noon_medians))**2)
    r_squared_np = 1 - (snp_res / snp_tot)
    print(f"R squared: {r_squared_np:.4f}")

    adj_np = 1-r_squared_np
    adjr_np = 1 - ((adj_np*23)/22)
    print(f"Adjusted R squared: {adjr_np:.4f}")
    
    noon_grad_p = round(c_p[0],2)
    noon_grad_err_p = np.sqrt(np.diag(ccov_p))
    noon_grad_err_p = round(noon_grad_err_p[0], 2)
    
    noon_HL_fit_p = (f'{noon_grad_p} ± {noon_grad_err_p}')
    
    
    c_m, ccov_m = np.polyfit(np.log10(HL_medians_m_pos/(1e6)), noon_medians, 1, cov=True)
    cc_m = np.poly1d(c_m)
    # print(f"Noon Fit -{error} Gradient:")
    # print(c_m[0])
    # print(f"Error of Noon Fit -{error}:")
    # print(np.sqrt(np.diag(ccov_m)))
    y_regress_noon_m = c_m[0]*np.log10(HL_medians_m_pos/1e6)+c_m[1]
    
    print("---------------")
    print('Noon Region -20% Fit')
    snm_res = np.sum((noon_medians - y_regress_noon_m)**2)
    snm_tot = np.sum((noon_medians - np.mean(noon_medians))**2)
    r_squared_nm = 1 - (snm_res / snm_tot)
    print(f"R squared: {r_squared_nm:.4f}")

    adj_nm = 1-r_squared_nm
    adjr_nm = 1 - ((adj_nm*23)/22)
    print(f"Adjusted R squared: {adjr_nm:.4f}")
    
    noon_grad_m = round(c_m[0],2)
    noon_grad_err_m = np.sqrt(np.diag(ccov_m))
    noon_grad_err_m = round(noon_grad_err_m[0], 2)
    
    noon_HL_fit_m = (f'{noon_grad_m} ± {noon_grad_err_m}')

    
    y_new_swirl = a[0]*np.log10(x_new)+a[1]
    y_new_dusk = b[0]*np.log10(x_new)+b[1]
    y_new_noon = c[0]*np.log10(x_new)+c[1]
    
    y_new_swirl_p = a_p[0]*np.log10(x_new)+a_p[1]
    y_new_dusk_p = b_p[0]*np.log10(x_new)+b_p[1]
    y_new_noon_p = c_p[0]*np.log10(x_new)+c_p[1]
    
    y_new_swirl_m = a_m[0]*np.log10(x_new)+a_m[1]
    y_new_dusk_m = b_m[0]*np.log10(x_new)+b_m[1]
    y_new_noon_m = c_m[0]*np.log10(x_new)+c_m[1]
    
    
    fig = plt.figure(figsize=(28,50))
    ax1 = plt.subplot(4,2,1)
    plt.subplots_adjust(hspace=0.1, wspace=0.15)
    
    # main plot
    ax1.scatter(Bz_pos_HL_pos[0:2]/1e6, Bz_pos_power[0:2], color='red', s=200, label='+Bz CME')
    ax1.scatter(Bz_pos_HL_pos[2:]/1e6, Bz_pos_power[2:], color='red', s=200, marker='v', label='+Bz Rarefaction (Deep)')

    ax1.scatter(Bz_neg_HL_pos[0]/1e6, Bz_neg_power[0], color='orange', s=200, label='-Bz CME') #
    ax1.scatter(Bz_neg_HL_pos[1]/1e6, Bz_neg_power[1], color='orange', s=200,marker=',', label='-Bz CIR') 

    ax1.scatter(By_pos_HL_pos[0]/1e6, By_pos_power[0], color='green', s=200, label='+By CME') # compression 1 2-5 incl visit 15
    ax1.scatter(By_pos_HL_pos[1:4]/1e6, By_pos_power[1:4], color='green', s=200, marker=",", label='+By CIR')
    ax1.scatter(By_pos_HL_pos[4:6]/1e6, By_pos_power[4:6], color='green', s=200, marker="v", label='+By Rarefaction (Deep)') # rarefraction
    ax1.scatter(By_pos_HL_pos[6]/1e6, By_pos_power[6], color='green', s=200, marker="^", label ='+By Rarefaction (Shallow)')
    ax1.scatter(By_pos_HL_pos[7:]/1e6, By_pos_power[7:], color='green', s=200, marker=",")

    ax1.scatter(By_neg_HL_pos[0:4]/1e6, By_neg_power[0:4], color='blue', s=200, marker="^", label='-By Rarefaction (Shallow)')
    ax1.scatter(By_neg_HL_pos[4:6]/1e6, By_neg_power[4:6], color='blue', s=200, label='-By CME') # 
    ax1.scatter(By_neg_HL_pos[6]/1e6, By_neg_power[6], color='blue', s=200, marker="v", label='-By Rarefaction (Deep)')
    ax1.scatter(By_neg_HL_pos[7:]/1e6, By_neg_power[7:], color='blue', s=200, marker=",", label='-By CIR')
    
    ax1.errorbar(HL_medians_pos/(1e6), polar_medians, yerr=polar_med_errs, xerr=med_HL_errs_pos/(1e6), fmt='.', color='lightgray')
    
    ax1.text(0.97e0,10,'a',style='italic',fontsize=40) #1.4 RHS, 0.05 LHS (1 is top)
    
    ax1.set_xlabel('High-Latitude Dayside $(+B_y)$ Reconnection Voltage', fontsize=22)
    ax1.set_ylabel('Median Total Polar Power (GW)', fontsize=22)

    ax1.set_xlim(2e-7, 2e0)
    ax1.set_ylim(-20,970)
    
    ax1.set_xscale('log')
    locmaj = matplotlib.ticker.LogLocator(base=10,numticks=30) 
    ax1.xaxis.set_major_locator(locmaj)
    ax1.tick_params(axis='x',which='minor',direction='in')
    ax1.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax1.minorticks_on()
    ax1.tick_params(which='major',direction='in',bottom=True, top=True, left=True, right=True, labelsize=20, length=10)
    ax1.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True, length=5)
    lgnd = plt.legend(framealpha=0.5, labelspacing=0.3, ncols=2, fontsize=18, loc="upper left")
    for handle in lgnd.legend_handles:
        handle.set_sizes([100.0])
    
    
    ax2 = plt.subplot(4,2,2)
    
    ax2.plot(x_new, y_new_swirl, '--',color='black', label=(f'Swirl $R^2$: {adjr_s:.2g}'), zorder=10)
    
    if error_plot == '20' or error_plot == '10' or error_plot == '50':
        
        ax2.plot(x_new, y_new_swirl_p, '--',color='darkviolet',linewidth=3, markersize=12, label=f'+ {error}%: {adjr_sp:.2g}')
        ax2.plot(x_new, y_new_swirl_m, '--',color='sandybrown', linewidth=3, markersize=12, label=f'- {error}%: {adjr_sm:.2g}')
        
        # plus error
        ax2.scatter(HL_medians_p_pos[0:2]/(1e6), medians_swirl[0:2],s=100,color='plum',marker="^")
        ax2.scatter(HL_medians_p_pos[2]/(1e6), medians_swirl[2],s=100,color='plum')
        ax2.scatter(HL_medians_p_pos[3:5]/(1e6), medians_swirl[3:5],s=100,color='plum',marker="^")
        ax2.scatter(HL_medians_p_pos[5:9]/(1e6), medians_swirl[5:9],s=100,color='plum')
        ax2.scatter(HL_medians_p_pos[9]/(1e6), medians_swirl[9],s=100,color='plum',marker="v")
        ax2.scatter(HL_medians_p_pos[10]/(1e6), medians_swirl[10],s=100,color='plum')
        ax2.scatter(HL_medians_p_pos[11]/(1e6), medians_swirl[11],s=100,color='plum',marker="v")
        ax2.scatter(HL_medians_p_pos[12:18]/(1e6), medians_swirl[12:18],s=100,color='plum',marker=",")
        ax2.scatter(HL_medians_p_pos[18:21]/(1e6), medians_swirl[18:21],s=100,color='plum',marker="v")
        ax2.scatter(HL_medians_p_pos[21]/(1e6), medians_swirl[21],s=100,color='plum',marker="^")
        ax2.scatter(HL_medians_p_pos[22:]/(1e6), medians_swirl[22:],s=100,color='plum',marker=",")
    
        # minus error
        ax2.scatter(HL_medians_m_pos[0:2]/(1e6), medians_swirl[0:2],s=100,color='peachpuff',marker="^")
        ax2.scatter(HL_medians_m_pos[2]/(1e6), medians_swirl[2],s=100,color='peachpuff')
        ax2.scatter(HL_medians_m_pos[3:5]/(1e6), medians_swirl[3:5],s=100,color='peachpuff',marker="^")
        ax2.scatter(HL_medians_m_pos[5:9]/(1e6), medians_swirl[5:9],s=100,color='peachpuff')
        ax2.scatter(HL_medians_m_pos[9]/(1e6), medians_swirl[9],s=100,color='peachpuff',marker="v")
        ax2.scatter(HL_medians_m_pos[10]/(1e6), medians_swirl[10],s=100,color='peachpuff')
        ax2.scatter(HL_medians_m_pos[11]/(1e6), medians_swirl[11],s=100,color='peachpuff',marker="v")
        ax2.scatter(HL_medians_m_pos[12:18]/(1e6), medians_swirl[12:18],s=100,color='peachpuff',marker=",")
        ax2.scatter(HL_medians_m_pos[18:21]/(1e6), medians_swirl[18:21],s=100,color='peachpuff',marker="v")
        ax2.scatter(HL_medians_m_pos[21]/(1e6), medians_swirl[21],s=100,color='peachpuff',marker="^")
        ax2.scatter(HL_medians_m_pos[22:]/(1e6), medians_swirl[22:],s=100,color='peachpuff',marker=",")
        
        ax2.errorbar(HL_medians_p_pos/(1e6), medians_swirl, yerr=swirl_med_errs, xerr=HL_med_errs_p_pos/(1e6), fmt='.', color='plum')
        
        ax2.errorbar(HL_medians_m_pos/(1e6), medians_swirl, yerr=swirl_med_errs, xerr=HL_med_errs_m_pos/(1e6), fmt='.', color='peachpuff')
    
    # main plot
    ax2.scatter(Bz_pos_HL_pos[0:2]/(1e6), Bz_pos_swirlpower[0:2], color='red', s=200,zorder=9)
    ax2.scatter(Bz_pos_HL_pos[2:]/(1e6), Bz_pos_swirlpower[2:], color='red', s=200, marker='v',zorder=9)

    ax2.scatter(Bz_neg_HL_pos[0]/(1e6), Bz_neg_swirlpower[0], color='orange', s=200,zorder=9) #
    ax2.scatter(Bz_neg_HL_pos[1]/(1e6), Bz_neg_swirlpower[1], color='orange', s=200,marker=',',zorder=9) 

    ax2.scatter(By_pos_HL_pos[0]/(1e6), By_pos_swirlpower[0], color='green', s=200,zorder=9) # compression 1 2-5 incl visit 15
    ax2.scatter(By_pos_HL_pos[1:4]/(1e6), By_pos_swirlpower[1:4], color='green', s=200, marker=",",zorder=9)
    ax2.scatter(By_pos_HL_pos[4:6]/(1e6), By_pos_swirlpower[4:6], color='green', s=200, marker="v",zorder=9) # rarefraction
    ax2.scatter(By_pos_HL_pos[6]/(1e6), By_pos_swirlpower[6], color='green', s=200, marker="^",zorder=9)
    ax2.scatter(By_pos_HL_pos[7:]/(1e6), By_pos_swirlpower[7:], color='green', s=200, marker=",",zorder=9)

    ax2.scatter(By_neg_HL_pos[0:4]/(1e6), By_neg_swirlpower[0:4], color='blue', s=200, marker="^",zorder=9)
    ax2.scatter(By_neg_HL_pos[4:6]/(1e6), By_neg_swirlpower[4:6], color='blue', s=200,zorder=9) # 
    ax2.scatter(By_neg_HL_pos[6]/(1e6), By_neg_swirlpower[6], color='blue', s=200, marker="v",zorder=9)
    ax2.scatter(By_neg_HL_pos[7:]/(1e6), By_neg_swirlpower[7:], color='blue', s=200, marker=",",zorder=9)
    
    ax2.errorbar(HL_medians_pos/(1e6), medians_swirl, yerr=swirl_med_errs, xerr=med_HL_errs_pos/(1e6), fmt='.', color='lightgray', zorder=9)
    
    ax2.text(0.6e0,1.4,'b',style='italic',fontsize=40) #0.97
    
    ax2.set_xlabel('High-Latitude Dayside $(+B_y)$ Reconnection Voltage (MV)', fontsize=22)
    ax2.set_ylabel('Median Swirl Region Power (GW)', fontsize=22)

    ax2.set_xlim(2e-7, 2e0) #2e-8, 2e-5
    ax2.set_ylim(-1,75)
    ax2.set_xscale('log')
    locmaj = matplotlib.ticker.LogLocator(base=10,numticks=30) 
    ax2.xaxis.set_major_locator(locmaj)
    ax2.tick_params(axis='x',which='minor',direction='in')
    ax2.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax2.minorticks_on()
    ax2.tick_params(which='major',direction='in',bottom=True, top=True, left=True, right=True, labelsize=20, length=10)
    ax2.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True, length=5)
    
    ax2.legend(framealpha=0.5,fontsize=21,loc ="upper left")
    
    
    ax3 = plt.subplot(4,2,3)
    
    ax3.plot(x_new, y_new_dusk,'--', color='black',linewidth=3, markersize=12, label=f'Dusk $R^2$: {adjr_d:.2g}', zorder=10)
    
    if error_plot == '20' or error_plot == '10' or error_plot == '50':
        
        ax3.plot(x_new, y_new_dusk_p, '--',color='darkviolet',linewidth=3, markersize=12, label=f'+ {error}%: {adjr_dp:.2g}')
        ax3.plot(x_new, y_new_dusk_m, '--',color='sandybrown', linewidth=3, markersize=12,label=f'- {error}%: {adjr_dm:.2g}')
        
        # plus error
        ax3.scatter(HL_medians_p_pos[0:2]/(1e6), dusk_medians[0:2],s=100,color='plum',marker="^")
        ax3.scatter(HL_medians_p_pos[2]/(1e6), dusk_medians[2],s=100,color='plum')
        ax3.scatter(HL_medians_p_pos[3:5]/(1e6), dusk_medians[3:5],s=100,color='plum',marker="^")
        ax3.scatter(HL_medians_p_pos[5:9]/(1e6), dusk_medians[5:9],s=100,color='plum')
        ax3.scatter(HL_medians_p_pos[9]/(1e6), dusk_medians[9],s=100,color='plum',marker="v")
        ax3.scatter(HL_medians_p_pos[10]/(1e6), dusk_medians[10],s=100,color='plum')
        ax3.scatter(HL_medians_p_pos[11]/(1e6), dusk_medians[11],s=100,color='plum',marker="v")
        ax3.scatter(HL_medians_p_pos[12:18]/(1e6), dusk_medians[12:18],s=100,color='plum',marker=",")
        ax3.scatter(HL_medians_p_pos[18:21]/(1e6), dusk_medians[18:21],s=100,color='plum',marker="v")
        ax3.scatter(HL_medians_p_pos[21]/(1e6), dusk_medians[21],s=100,color='plum',marker="^")
        ax3.scatter(HL_medians_p_pos[22:]/(1e6), dusk_medians[22:],s=100,color='plum',marker=",")
    
        # minus error
        ax3.scatter(HL_medians_m_pos[0:2]/(1e6), dusk_medians[0:2],s=100,color='peachpuff',marker="^")
        ax3.scatter(HL_medians_m_pos[2]/(1e6), dusk_medians[2],s=100,color='peachpuff')
        ax3.scatter(HL_medians_m_pos[3:5]/(1e6), dusk_medians[3:5],s=100,color='peachpuff',marker="^")
        ax3.scatter(HL_medians_m_pos[5:9]/(1e6), dusk_medians[5:9],s=100,color='peachpuff')
        ax3.scatter(HL_medians_m_pos[9]/(1e6), dusk_medians[9],s=100,color='peachpuff',marker="v")
        ax3.scatter(HL_medians_m_pos[10]/(1e6), dusk_medians[10],s=100,color='peachpuff')
        ax3.scatter(HL_medians_m_pos[11]/(1e6), dusk_medians[11],s=100,color='peachpuff',marker="v")
        ax3.scatter(HL_medians_m_pos[12:18]/(1e6), dusk_medians[12:18],s=100,color='peachpuff',marker=",")
        ax3.scatter(HL_medians_m_pos[18:21]/(1e6), dusk_medians[18:21],s=100,color='peachpuff',marker="v")
        ax3.scatter(HL_medians_m_pos[21]/(1e6), dusk_medians[21],s=100,color='peachpuff',marker="^")
        ax3.scatter(HL_medians_m_pos[22:]/(1e6), dusk_medians[22:],s=100,color='peachpuff',marker=",")
        
        ax3.errorbar(HL_medians_p_pos/(1e6), dusk_medians, yerr=dusk_med_errs, xerr=HL_med_errs_p_pos/(1e6), fmt='.', color='plum')
        
        ax3.errorbar(HL_medians_m_pos/(1e6), dusk_medians, yerr=dusk_med_errs, xerr=HL_med_errs_m_pos/(1e6), fmt='.', color='peachpuff')
    
    
    # main plot
    ax3.scatter(Bz_pos_HL_pos[0:2]/(1e6), Bz_pos_duskpower[0:2], color='red', s=200,zorder=9)
    ax3.scatter(Bz_pos_HL_pos[2:]/(1e6), Bz_pos_duskpower[2:], color='red', s=200, marker='v',zorder=9)

    ax3.scatter(Bz_neg_HL_pos[0]/(1e6), Bz_neg_duskpower[0], color='orange', s=200,zorder=9) #
    ax3.scatter(Bz_neg_HL_pos[1]/(1e6), Bz_neg_duskpower[1], color='orange', s=200,marker=',',zorder=9)

    ax3.scatter(By_pos_HL_pos[0]/(1e6), By_pos_duskpower[0], color='green', s=200,zorder=9) # compression 1 2-5 incl visit 15
    ax3.scatter(By_pos_HL_pos[1:4]/(1e6), By_pos_duskpower[1:4], color='green', s=200, marker=",",zorder=9)
    ax3.scatter(By_pos_HL_pos[4:6]/(1e6), By_pos_duskpower[4:6], color='green', s=200, marker="v",zorder=9) # rarefraction
    ax3.scatter(By_pos_HL_pos[6]/(1e6), By_pos_duskpower[6], color='green', s=200, marker="^",zorder=9)
    ax3.scatter(By_pos_HL_pos[7:]/(1e6), By_pos_duskpower[7:], color='green', s=200, marker=",",zorder=9)

    ax3.scatter(By_neg_HL_pos[0:4]/(1e6), By_neg_duskpower[0:4], color='blue', s=200, marker="^",zorder=9)
    ax3.scatter(By_neg_HL_pos[4:6]/(1e6), By_neg_duskpower[4:6], color='blue', s=200,zorder=9) # 
    ax3.scatter(By_neg_HL_pos[6]/(1e6), By_neg_duskpower[6], color='blue', s=200, marker="v",zorder=9)
    ax3.scatter(By_neg_HL_pos[7:]/(1e6), By_neg_duskpower[7:], color='blue', s=200, marker=",",zorder=9)
    
    ax3.errorbar(HL_medians_pos/(1e6), dusk_medians, yerr=dusk_med_errs, xerr=med_HL_errs_pos/(1e6), fmt='.', color='lightgray', zorder=9)
    
    ax3.text(0.6e0,-13,'c',style='italic',fontsize=40) # 0.7, 0.97
    
    ax3.set_xlabel('High-Latitude Dayside $(+B_y)$ Reconnection Voltage (MV)', fontsize=22)
    ax3.set_ylabel('Median Dusk Region Power (GW)', fontsize=22)
    
    ax3.set_xlim(2e-7, 2e0)
    #ax3.set_ylim(-20,520)
    ax3.set_xscale('log')
    locmaj = matplotlib.ticker.LogLocator(base=10,numticks=30) 
    ax3.xaxis.set_major_locator(locmaj)
    ax3.tick_params(axis='x',which='minor',direction='in')
    ax3.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax3.minorticks_on()
    ax3.tick_params(which='major',direction='in',bottom=True, top=True, left=True, right=True, labelsize=20, length=10)
    ax3.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True, length=5)
    ax3.legend(framealpha=0.5,fontsize=21,loc ="upper left")
    
    
    ax4 = plt.subplot(4,2,4)
    
    ax4.plot(x_new, y_new_noon, '--',color='black', linewidth=3, markersize=12,label=f'Noon $R^2$: {adjr_n:.2g}', zorder=10)
    
    if error_plot == '20' or error_plot == '10' or error_plot == '50':
        ax4.plot(x_new, y_new_noon_p, '--',color='darkviolet',linewidth=3, markersize=12, label=f'+ {error}%: {adjr_np:.2g}')
        ax4.plot(x_new, y_new_noon_m, '--',color='sandybrown',linewidth=3, markersize=12, label=f'- {error}%: {adjr_nm:.2g}')
        
        # plus error
        ax4.scatter(HL_medians_p_pos[0:2]/(1e6), noon_medians[0:2],s=100,color='plum',marker="^")
        ax4.scatter(HL_medians_p_pos[2]/(1e6), noon_medians[2],s=100,color='plum')
        ax4.scatter(HL_medians_p_pos[3:5]/(1e6), noon_medians[3:5],s=100,color='plum',marker="^")
        ax4.scatter(HL_medians_p_pos[5:9]/(1e6), noon_medians[5:9],s=100,color='plum')
        ax4.scatter(HL_medians_p_pos[9]/(1e6), noon_medians[9],s=100,color='plum',marker="v")
        ax4.scatter(HL_medians_p_pos[10]/(1e6), noon_medians[10],s=100,color='plum')
        ax4.scatter(HL_medians_p_pos[11]/(1e6), noon_medians[11],s=100,color='plum',marker="v")
        ax4.scatter(HL_medians_p_pos[12:18]/(1e6), noon_medians[12:18],s=100,color='plum',marker=",")
        ax4.scatter(HL_medians_p_pos[18:21]/(1e6), noon_medians[18:21],s=100,color='plum',marker="v")
        ax4.scatter(HL_medians_p_pos[21]/(1e6), noon_medians[21],s=100,color='plum',marker="^")
        ax4.scatter(HL_medians_p_pos[22:]/(1e6), noon_medians[22:],s=100,color='plum',marker=",")
    
        # minus error
        ax4.scatter(HL_medians_m_pos[0:2]/(1e6), noon_medians[0:2],s=100,color='peachpuff',marker="^")
        ax4.scatter(HL_medians_m_pos[2]/(1e6), noon_medians[2],s=100,color='peachpuff')
        ax4.scatter(HL_medians_m_pos[3:5]/(1e6), noon_medians[3:5],s=100,color='peachpuff',marker="^")
        ax4.scatter(HL_medians_m_pos[5:9]/(1e6), noon_medians[5:9],s=100,color='peachpuff')
        ax4.scatter(HL_medians_m_pos[9]/(1e6), noon_medians[9],s=100,color='peachpuff',marker="v")
        ax4.scatter(HL_medians_m_pos[10]/(1e6), noon_medians[10],s=100,color='peachpuff')
        ax4.scatter(HL_medians_m_pos[11]/(1e6), noon_medians[11],s=100,color='peachpuff',marker="v")
        ax4.scatter(HL_medians_m_pos[12:18]/(1e6), noon_medians[12:18],s=100,color='peachpuff',marker=",")
        ax4.scatter(HL_medians_m_pos[18:21]/(1e6), noon_medians[18:21],s=100,color='peachpuff',marker="v")
        ax4.scatter(HL_medians_m_pos[21]/(1e6), noon_medians[21],s=100,color='peachpuff',marker="^")
        ax4.scatter(HL_medians_m_pos[22:]/(1e6), noon_medians[22:],s=100,color='peachpuff',marker=",")
        
        ax4.errorbar(HL_medians_p_pos/(1e6), noon_medians, yerr=noon_med_errs, xerr=HL_med_errs_p_pos/(1e6), fmt='.', color='plum')
        
        ax4.errorbar(HL_medians_m_pos/(1e6), noon_medians, yerr=noon_med_errs, xerr=HL_med_errs_m_pos/(1e6), fmt='.', color='peachpuff')
    
    
    # main plot
    ax4.scatter(Bz_pos_HL_pos[0:2]/(1e6), Bz_pos_noonpower[0:2], color='red', s=200,zorder=9)
    ax4.scatter(Bz_pos_HL_pos[2:]/(1e6), Bz_pos_noonpower[2:], color='red', s=200, marker='v',zorder=9)

    ax4.scatter(Bz_neg_HL_pos[0]/(1e6), Bz_neg_noonpower[0], color='orange', s=200,zorder=9) #
    ax4.scatter(Bz_neg_HL_pos[1]/(1e6), Bz_neg_noonpower[1], color='orange', s=200,marker=',',zorder=9) 

    ax4.scatter(By_pos_HL_pos[0]/(1e6), By_pos_noonpower[0], color='green', s=200,zorder=9) # compression 1 2-5 incl visit 15
    ax4.scatter(By_pos_HL_pos[1:4]/(1e6), By_pos_noonpower[1:4], color='green', s=200, marker=",",zorder=9)
    ax4.scatter(By_pos_HL_pos[4:6]/(1e6), By_pos_noonpower[4:6], color='green', s=200, marker="v",zorder=9) # rarefraction
    ax4.scatter(By_pos_HL_pos[6]/(1e6), By_pos_noonpower[6], color='green', s=200, marker="^",zorder=9)
    ax4.scatter(By_pos_HL_pos[7:]/(1e6), By_pos_noonpower[7:], color='green', s=200, marker=",",zorder=9)

    ax4.scatter(By_neg_HL_pos[0:4]/(1e6), By_neg_noonpower[0:4], color='blue', s=200, marker="^",zorder=9)
    ax4.scatter(By_neg_HL_pos[4:6]/(1e6), By_neg_noonpower[4:6], color='blue', s=200,zorder=9)
    ax4.scatter(By_neg_HL_pos[6]/(1e6), By_neg_noonpower[6], color='blue', s=200, marker="v",zorder=9)
    ax4.scatter(By_neg_HL_pos[7:]/(1e6), By_neg_noonpower[7:], color='blue', s=200, marker=",",zorder=9)
    
    ax4.errorbar(HL_medians_pos/(1e6), noon_medians, yerr=noon_med_errs, xerr=med_HL_errs_pos/(1e6), fmt='.', color='lightgray', zorder=9)
    
    ax4.set_xlabel('High-Latitude Dayside $(+B_y)$ Reconnection Voltage (MV)', fontsize=22)
    ax4.set_ylabel('Median Noon Region Power (GW)', fontsize=22)

    ax4.set_xlim(2e-7, 2e0) #6e-8, 3e-5
    #ax4.set_ylim(-20,520)
    ax4.set_xscale('log')
    locmaj = matplotlib.ticker.LogLocator(base=10,numticks=30) 
    ax4.xaxis.set_major_locator(locmaj)
    ax4.tick_params(axis='x',which='minor',direction='in')
    ax4.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax4.minorticks_on()
    ax4.tick_params(which='major',direction='in',bottom=True, top=True, left=True, right=True, labelsize=20, length=10)
    ax4.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True, length=5)
    ax4.legend(framealpha=0.5,fontsize=21,loc ="upper left")
    
    ax4.text(0.6e0,-12,'d',style='italic',fontsize=40) #0.7, 0.97
    
    # save plot
    saveloc = (f'{root_saves}median_HL_rec_pos_V_vs_region_power_{error_plot}.jpg') 
    #saveloc = (f'/Users/hannah/OneDrive - Lancaster University/aurora/median_HL_rec_pos_V_vs_region_power_{error}.jpg')
    plt.savefig(saveloc,bbox_inches='tight',dpi=400)
    
if plotting == 'Gershman':  # out of date
    
    x_new = np.linspace(min(gersh_medians_p/1e9),max(gersh_medians/1e9), len(gersh_medians)) # p for min and m for max
    
    a, acov = np.polyfit(np.log(gersh_medians/1e9), medians_swirl, 1, cov=True)
    aa = np.poly1d(a)
    print("Swirl Fit Gradient:")
    print(a[0])
    print("Error of Swirl Fit:")
    print(np.sqrt(np.diag(acov)))
    
    swirl_grad = round(a[0],2)
    swirl_grad_err = np.sqrt(np.diag(acov))
    swirl_grad_err = round(swirl_grad_err[0], 2)
    
    swirl_g_fit = (f'{swirl_grad} ± {swirl_grad_err}')
    
    
    a_p,acov_p = np.polyfit(np.log10(gersh_medians_p/1e9), medians_swirl,1, cov=True)
    aa_p = np.poly1d(a_p)
    print(f"Swirl Fit +{error} Gradient:")
    print(a_p[0])
    print(f"Error of Swirl Fit +{error}:")
    print(np.sqrt(np.diag(acov_p)))
    
    swirl_grad_p = round(a_p[0],2)
    swirl_grad_err_p = np.sqrt(np.diag(acov_p))
    swirl_grad_err_p = round(swirl_grad_err_p[0], 2)
    
    swirl_g_fit_p = (f'{swirl_grad_p} ± {swirl_grad_err_p}')
    
    
    a_m,acov_m = np.polyfit(np.log10(gersh_medians_m/1e9), medians_swirl,1, cov=True)
    aa_m = np.poly1d(a_m)
    print(f"Swirl Fit -{error} Gradient:")
    print(a_m[0])
    print(f"Error of Swirl Fit -{error}:")
    print(np.sqrt(np.diag(acov_m)))
    
    swirl_grad_m = round(a_m[0],2)
    swirl_grad_err_m = np.sqrt(np.diag(acov_m))
    swirl_grad_err_m = round(swirl_grad_err_m[0], 2)
    
    swirl_g_fit_m = (f'{swirl_grad_m} ± {swirl_grad_err_m}')
   
    
    
    b, bcov = np.polyfit(np.log10(gersh_medians/1e9), dusk_medians, 1, cov=True)
    bb = np.poly1d(b)
    print("Dusk Fit Gradient:")
    print(b[0])
    print("Error of Dusk Fit:")
    print(np.sqrt(np.diag(bcov)))
    
    dusk_grad = round(b[0],2)
    dusk_grad_err = np.sqrt(np.diag(bcov))
    dusk_grad_err = round(dusk_grad_err[0], 2)
    
    dusk_g_fit = (f'{dusk_grad} ± {dusk_grad_err}')
    
    
    b_p, bcov_p = np.polyfit(np.log10(gersh_medians_p/1e9), dusk_medians, 1, cov=True)
    bb_p = np.poly1d(b_p)
    print(f"Dusk Fit +{error} Gradient:")
    print(b_p[0])
    print(f"Error of Dusk Fit +{error}:")
    print(np.sqrt(np.diag(bcov_p)))
    
    dusk_grad_p = round(b_p[0],2)
    dusk_grad_err_p = np.sqrt(np.diag(bcov_p))
    dusk_grad_err_p = round(dusk_grad_err_p[0], 2)
    
    dusk_g_fit_p = (f'{dusk_grad_p} ± {dusk_grad_err_p}')
    
    
    b_m, bcov_m = np.polyfit(np.log10(gersh_medians_m/1e9), dusk_medians, 1, cov=True)
    bb_m = np.poly1d(b_m)
    print(f"Dusk Fit -{error} Gradient:")
    print(b_m[0])
    print(f"Error of Dusk Fit -{error}:")
    print(np.sqrt(np.diag(bcov_m)))
    
    dusk_grad_m = round(b_m[0],2)
    dusk_grad_err_m = np.sqrt(np.diag(bcov_m))
    dusk_grad_err_m = round(dusk_grad_err_m[0], 2)
    
    dusk_g_fit_m = (f'{dusk_grad_m} ± {dusk_grad_err_m}')
    
    
    
    c, ccov = np.polyfit(np.log10(gersh_medians/1e9), noon_medians, 1, cov=True)
    cc = np.poly1d(c)
    print("Noon Fit Gradient:")
    print(c[0])
    print("Error of Noon Fit:")
    print(np.sqrt(np.diag(ccov)))
    
    noon_grad = round(c[0],2)
    noon_grad_err = np.sqrt(np.diag(ccov))
    noon_grad_err = round(noon_grad_err[0], 2)
    
    noon_g_fit = (f'{noon_grad} ± {noon_grad_err}')
    
    
    c_p, ccov_p = np.polyfit(np.log10(gersh_medians_p/1e9), noon_medians, 1, cov=True)
    cc_p = np.poly1d(c_p)
    print(f"Noon Fit +{error} Gradient:")
    print(c_p[0])
    print(f"Error of Noon Fit +{error}:")
    print(np.sqrt(np.diag(ccov_p)))
    
    noon_grad_p = round(c_p[0],2)
    noon_grad_err_p = np.sqrt(np.diag(ccov_p))
    noon_grad_err_p = round(noon_grad_err_p[0], 2)
    
    noon_g_fit_p = (f'{noon_grad_p} ± {noon_grad_err_p}')
    
    
    c_m, ccov_m = np.polyfit(np.log10(gersh_medians_m/1e9), noon_medians, 1, cov=True)
    cc_m = np.poly1d(c_m)
    print(f"Noon Fit -{error} Gradient:")
    print(c_m[0])
    print(f"Error of Noon Fit -{error}:")
    print(np.sqrt(np.diag(ccov_m)))
    
    noon_grad_m = round(c_m[0],2)
    noon_grad_err_m = np.sqrt(np.diag(ccov_m))
    noon_grad_err_m = round(noon_grad_err_m[0], 2)
    
    noon_g_fit_m = (f'{noon_grad_m} ± {noon_grad_err_m}')
    
    
    y_new_swirl = a[0]*np.log10(x_new)+a[1]
    y_new_dusk = b[0]*np.log10(x_new)+b[1]
    y_new_noon = c[0]*np.log10(x_new)+c[1]
    
    y_new_swirl_p = a_p[0]*np.log10(x_new)+a_p[1]
    y_new_dusk_p = b_p[0]*np.log10(x_new)+b_p[1]
    y_new_noon_p = c_p[0]*np.log10(x_new)+c_p[1]
    
    y_new_swirl_m = a_m[0]*np.log10(x_new)+a_m[1]
    y_new_dusk_m = b_m[0]*np.log10(x_new)+b_m[1]
    y_new_noon_m = c_m[0]*np.log10(x_new)+c_m[1]
    
    
      
    fig = plt.figure(figsize=(28,50))
    ax1 = plt.subplot(4,2,1)
    plt.subplots_adjust(hspace=0.1, wspace=0.15)
    
    # main plot
    ax1.scatter(Bz_pos_gersh/1e9, Bz_pos_power, color='red', s=200, label='+Bz CME')

    ax1.scatter(Bz_neg_gersh/1e9, Bz_neg_power, color='orange', s=200, label='-Bz CME') # all compression 1
       
    ax1.scatter(By_pos_gersh[0:2]/1e9, By_pos_power[0:2], color='green', s=200, label='+By CME') # compression 1 2-5 incl visit 15
    ax1.scatter(By_pos_gersh[2:6]/1e9, By_pos_power[2:6], color='green', s=200, marker=",", label='+By CIR')
    ax1.scatter(By_pos_gersh[6:8]/1e9, By_pos_power[6:8], color='green', s=200, marker="v", label='+By Rarefraction (Deep)') # rarefraction
    ax1.scatter(By_pos_gersh[8]/1e9, By_pos_power[8], color='green', s=200, marker="^", label ='+By Rarefaction (Shallow)')
    ax1.scatter(By_pos_gersh[9:]/1e9, By_pos_power[9:], color='green', s=200, marker=",")
       
    ax1.scatter(By_neg_gersh[0:4]/1e9, By_neg_power[0:4], color='blue', s=200, marker="^", label='-By Rarefaction (Shallow)')
    ax1.scatter(By_neg_gersh[4]/1e9, By_neg_power[4], color='blue', s=200, label='-By CME') # 
    ax1.scatter(By_neg_gersh[5]/1e9, By_neg_power[5], color='blue', s=200, marker="v", label='-By Rarefaction (Deep)')
    ax1.scatter(By_neg_gersh[6]/1e9, By_neg_power[6], color='blue', s=200)
    ax1.scatter(By_neg_gersh[7]/1e9, By_neg_power[7], color='blue', s=200, marker="v")
    ax1.scatter(By_neg_gersh[8:10]/1e9, By_neg_power[8:10], color='blue', s=200, marker=",", label='-By CIR')
    ax1.scatter(By_neg_gersh[10]/1e9, By_neg_power[10], color='blue', s=200, marker="v")#, label='-By CIR (2nd)') 
    
    ax1.errorbar(gersh_medians/1e9, polar_medians, yerr=polar_med_errs, xerr=med_gersh_errs/1e9, fmt='.', color='lightgray')
    
    ax1.text(1.5e4,10,'a',style='italic',fontsize=40) #1.4 RHS, 0.05 LHS (1 is top)
    
    ax1.set_xlabel('Gershman Reconnection Power (GW)', fontsize=22)
    ax1.set_ylabel('Median Total Polar Power (GW)', fontsize=22)
    ax1.set_xscale('log')
    ax1.set_ylim(-20, 960)
    locmaj = matplotlib.ticker.LogLocator(base=10,numticks=30) 
    ax1.xaxis.set_major_locator(locmaj)
    ax1.tick_params(axis='x',which='minor',direction='in')
    ax1.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax1.minorticks_on()
    ax1.tick_params(which='major',direction='in',bottom=True, top=True, left=True, right=True, labelsize=20, length=10)
    ax1.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True, length=5)
    lgnd = plt.legend(framealpha=0.5, labelspacing=0.3, fontsize=18, loc="upper right") #if using equation legend
    for handle in lgnd.legend_handles:
        handle.set_sizes([100.0])
    
    
    
    ax2 = plt.subplot(4,2,2)
    
    ax2.plot(x_new, y_new_swirl, '--', color='black', label=f'Swirl: {swirl_g_fit}', zorder=10)
    
    # plus error
    ax2.scatter(gersh_medians_p[0:2]/1e9, medians_swirl[0:2],s=100,color='thistle',marker="^")
    ax2.scatter(gersh_medians_p[2]/1e9, medians_swirl[2],s=100,color='thistle')
    ax2.scatter(gersh_medians_p[3:5]/1e9, medians_swirl[3:5],s=100,color='thistle',marker="^")
    ax2.scatter(gersh_medians_p[5:9]/1e9, medians_swirl[5:9],s=100,color='thistle')
    ax2.scatter(gersh_medians_p[9]/1e9, medians_swirl[9],s=100,color='thistle',marker="v")
    ax2.scatter(gersh_medians_p[10]/1e9, medians_swirl[10],s=100,color='thistle')
    ax2.scatter(gersh_medians_p[11]/1e9, medians_swirl[11],s=100,color='thistle',marker="v")
    ax2.scatter(gersh_medians_p[12:18]/1e9, medians_swirl[12:18],s=100,color='thistle',marker=",")
    ax2.scatter(gersh_medians_p[18:21]/1e9, medians_swirl[18:21],s=100,color='thistle',marker="v")
    ax2.scatter(gersh_medians_p[21]/1e9, medians_swirl[21],s=100,color='thistle',marker="^")
    ax2.scatter(gersh_medians_p[22:]/1e9, medians_swirl[22:],s=100,color='thistle',marker=",")

    # minus error
    ax2.scatter(gersh_medians_m[0:2]/1e9, medians_swirl[0:2],s=100,color='peachpuff',marker="^")
    ax2.scatter(gersh_medians_m[2]/1e9, medians_swirl[2],s=100,color='peachpuff')
    ax2.scatter(gersh_medians_m[3:5]/1e9, medians_swirl[3:5],s=100,color='peachpuff',marker="^")
    ax2.scatter(gersh_medians_m[5:9]/1e9, medians_swirl[5:9],s=100,color='peachpuff')
    ax2.scatter(gersh_medians_m[9]/1e9, medians_swirl[9],s=100,color='peachpuff',marker="v")
    ax2.scatter(gersh_medians_m[10]/1e9, medians_swirl[10],s=100,color='peachpuff')
    ax2.scatter(gersh_medians_m[11]/1e9, medians_swirl[11],s=100,color='peachpuff',marker="v")
    ax2.scatter(gersh_medians_m[12:18]/1e9, medians_swirl[12:18],s=100,color='peachpuff',marker=",")
    ax2.scatter(gersh_medians_m[18:21]/1e9, medians_swirl[18:21],s=100,color='peachpuff',marker="v")
    ax2.scatter(gersh_medians_m[21]/1e9, medians_swirl[21],s=100,color='peachpuff',marker="^")
    ax2.scatter(gersh_medians_m[22:]/1e9, medians_swirl[22:],s=100,color='peachpuff',marker=",")
    
    ax2.errorbar(gersh_medians_p/1e9, medians_swirl, yerr=swirl_med_errs, xerr=med_gersh_errs_p/1e9, fmt='.', color='thistle')
    
    ax2.plot(x_new, y_new_swirl_p, '--',color='thistle', label=f'+ {error}: {swirl_g_fit_p}')
    
    ax2.errorbar(gersh_medians_m/1e9, medians_swirl, yerr=swirl_med_errs, xerr=med_gersh_errs_m/1e9, fmt='.', color='peachpuff')
    
    ax2.plot(x_new, y_new_swirl_m, '--',color='peachpuff', label=f'- {error}%: {swirl_g_fit_m}')
    

    # main plot
    ax2.scatter(Bz_pos_gersh/1e9, Bz_pos_swirlpower, color='red', s=200)#, label='+Bz CME')

    ax2.scatter(Bz_neg_gersh/1e9, Bz_neg_swirlpower, color='orange', s=200)#, label='-Bz CME') # all compression 1

    ax2.scatter(By_pos_gersh[0:2]/1e9, By_pos_swirlpower[0:2], color='green', s=200)#, label='+By CME') # compression 1 2-5 incl visit 15
    ax2.scatter(By_pos_gersh[2:6]/1e9, By_pos_swirlpower[2:6], color='green', s=200, marker=",")#, label='+By CIR')
    ax2.scatter(By_pos_gersh[6:8]/1e9, By_pos_swirlpower[6:8], color='green', s=200, marker="v")#, label='+By Rarefraction (Deep)') # rarefraction
    ax2.scatter(By_pos_gersh[8]/1e9, By_pos_swirlpower[8], color='green', s=200, marker="^")#, label ='+By Rarefaction (Shallow)')
    ax2.scatter(By_pos_gersh[9:]/1e9, By_pos_swirlpower[9:], color='green', s=200, marker=",")

    ax2.scatter(By_neg_gersh[0:4]/1e9, By_neg_swirlpower[0:4], color='blue', s=200, marker="^")#, label='-By Rarefaction (Shallow)')
    ax2.scatter(By_neg_gersh[4]/1e9, By_neg_swirlpower[4], color='blue', s=200)#, label='-By CME') # 
    ax2.scatter(By_neg_gersh[5]/1e9, By_neg_swirlpower[5], color='blue', s=200, marker="v")#, label='-By Rarefaction (Deep)')
    ax2.scatter(By_neg_gersh[6]/1e9, By_neg_swirlpower[6], color='blue', s=200)
    ax2.scatter(By_neg_gersh[7]/1e9, By_neg_swirlpower[7], color='blue', s=200, marker="v")
    ax2.scatter(By_neg_gersh[8:10]/1e9, By_neg_swirlpower[8:10], color='blue', s=200, marker=",")#, label='-By CIR')
    ax2.scatter(By_neg_gersh[10]/1e9, By_neg_swirlpower[10], color='blue', s=200, marker="v")#, label='-By CIR (2nd)')
    
    ax2.errorbar(gersh_medians/1e9, medians_swirl, yerr=swirl_med_errs, xerr=med_gersh_errs/1e9, fmt='.', color='lightgray', zorder=9)
    
    ax2.text(1.5e4,-2.4,'b',style='italic',fontsize=40) #6.6
    
    ax2.set_xlabel('Gershman Reconnection Power (GW)', fontsize=22)
    ax2.set_ylabel('Median Swirl Region Power (GW)', fontsize=22)
    
    ax2.set_ylim(-5,85)
    
    ax2.set_xscale('log')
    locmaj = matplotlib.ticker.LogLocator(base=10,numticks=30) 
    ax2.xaxis.set_major_locator(locmaj)
    ax2.tick_params(axis='x',which='minor',direction='in')
    ax2.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax2.minorticks_on()
    ax2.tick_params(which='major',direction='in',bottom=True, top=True, left=True, right=True, labelsize=20, length=10)
    ax2.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True, length=5)
    
    ax2.legend(framealpha=0.5,fontsize=18,loc ="upper right")
    
    
    ax3 = plt.subplot(4,2,3)
    
    ax3.plot(x_new, y_new_dusk, '--', color='black', label=f'Dusk: {dusk_g_fit}', zorder=10)
    
    # plus error
    ax3.scatter(gersh_medians_p[0:2]/1e9, dusk_medians[0:2],s=100,color='thistle',marker="^")
    ax3.scatter(gersh_medians_p[2]/1e9, dusk_medians[2],s=100,color='thistle')
    ax3.scatter(gersh_medians_p[3:5]/1e9, dusk_medians[3:5],s=100,color='thistle',marker="^")
    ax3.scatter(gersh_medians_p[5:9]/1e9, dusk_medians[5:9],s=100,color='thistle')
    ax3.scatter(gersh_medians_p[9]/1e9, dusk_medians[9],s=100,color='thistle',marker="v")
    ax3.scatter(gersh_medians_p[10]/1e9, dusk_medians[10],s=100,color='thistle')
    ax3.scatter(gersh_medians_p[11]/1e9, dusk_medians[11],s=100,color='thistle',marker="v")
    ax3.scatter(gersh_medians_p[12:18]/1e9, dusk_medians[12:18],s=100,color='thistle',marker=",")
    ax3.scatter(gersh_medians_p[18:21]/1e9, dusk_medians[18:21],s=100,color='thistle',marker="v")
    ax3.scatter(gersh_medians_p[21]/1e9, dusk_medians[21],s=100,color='thistle',marker="^")
    ax3.scatter(gersh_medians_p[22:]/1e9, dusk_medians[22:],s=100,color='thistle',marker=",")

    # minus error
    ax3.scatter(gersh_medians_m[0:2]/1e9, dusk_medians[0:2],s=100,color='peachpuff',marker="^")
    ax3.scatter(gersh_medians_m[2]/1e9, dusk_medians[2],s=100,color='peachpuff')
    ax3.scatter(gersh_medians_m[3:5]/1e9, dusk_medians[3:5],s=100,color='peachpuff',marker="^")
    ax3.scatter(gersh_medians_m[5:9]/1e9, dusk_medians[5:9],s=100,color='peachpuff')
    ax3.scatter(gersh_medians_m[9]/1e9, dusk_medians[9],s=100,color='peachpuff',marker="v")
    ax3.scatter(gersh_medians_m[10]/1e9, dusk_medians[10],s=100,color='peachpuff')
    ax3.scatter(gersh_medians_m[11]/1e9, dusk_medians[11],s=100,color='peachpuff',marker="v")
    ax3.scatter(gersh_medians_m[12:18]/1e9, dusk_medians[12:18],s=100,color='peachpuff',marker=",")
    ax3.scatter(gersh_medians_m[18:21]/1e9, dusk_medians[18:21],s=100,color='peachpuff',marker="v")
    ax3.scatter(gersh_medians_m[21]/1e9, dusk_medians[21],s=100,color='peachpuff',marker="^")
    ax3.scatter(gersh_medians_m[22:]/1e9, dusk_medians[22:],s=100,color='peachpuff',marker=",")
    
    ax3.errorbar(gersh_medians_p/1e9, dusk_medians, yerr=dusk_med_errs, xerr=med_gersh_errs_p/1e9, fmt='.', color='thistle')
    
    ax3.plot(x_new, y_new_dusk_p, '--',color='thistle', label=f'+ {error}: {dusk_g_fit_p}')
    
    ax3.errorbar(gersh_medians_m/1e9, dusk_medians, yerr=dusk_med_errs, xerr=med_gersh_errs_m/1e9, fmt='.', color='peachpuff')
    
    ax3.plot(x_new, y_new_dusk_m, '--',color='peachpuff', label=f'- {error}%: {dusk_g_fit_m}')
    
    # main plot
    ax3.scatter(Bz_pos_gersh/1e9, Bz_pos_duskpower, color='red', s=200)#, label='+Bz CME')

    ax3.scatter(Bz_neg_gersh/1e9, Bz_neg_duskpower, color='orange', s=200)#, label='-Bz CME') # all compression 1

    ax3.scatter(By_pos_gersh[0:2]/1e9, By_pos_duskpower[0:2], color='green', s=200)#, label='+By CME') # compression 1 2-5 incl visit 15
    ax3.scatter(By_pos_gersh[2:6]/1e9, By_pos_duskpower[2:6], color='green', s=200, marker=",")#, label='+By CIR')
    ax3.scatter(By_pos_gersh[6:8]/1e9, By_pos_duskpower[6:8], color='green', s=200, marker="v")#, label='+By Rarefraction (Deep)') # rarefraction
    ax3.scatter(By_pos_gersh[8]/1e9, By_pos_duskpower[8], color='green', s=200, marker="^")#, label ='+By Rarefaction (Shallow)')
    ax3.scatter(By_pos_gersh[9:]/1e9, By_pos_duskpower[9:], color='green', s=200, marker=",")

    ax3.scatter(By_neg_gersh[0:4]/1e9, By_neg_duskpower[0:4], color='blue', s=200, marker="^")#, label='-By Rarefaction (Shallow)')
    ax3.scatter(By_neg_gersh[4]/1e9, By_neg_duskpower[4], color='blue', s=200)#, label='-By CME') # 
    ax3.scatter(By_neg_gersh[5]/1e9, By_neg_duskpower[5], color='blue', s=200, marker="v")#, label='-By Rarefaction (Deep)')
    ax3.scatter(By_neg_gersh[6]/1e9, By_neg_duskpower[6], color='blue', s=200)
    ax3.scatter(By_neg_gersh[7]/1e9, By_neg_duskpower[7], color='blue', s=200, marker="v")
    ax3.scatter(By_neg_gersh[8:10]/1e9, By_neg_duskpower[8:10], color='blue', s=200, marker=",")#, label='-By CIR')
    ax3.scatter(By_neg_gersh[10]/1e9, By_neg_duskpower[10], color='blue', s=200, marker="v")#, label='-By CIR (2nd)')
    
    ax3.errorbar(gersh_medians/1e9, dusk_medians, yerr=dusk_med_errs, xerr=med_gersh_errs/1e9, fmt='.', color='lightgray', zorder=9)
    
    ax3.text(1.5e4,-14,'c',style='italic',fontsize=40) #9.5e3 for all but 1.3e4 for 50%
    
    ax3.set_xlabel('Gershman Reconnection Power (GW)', fontsize=22)
    ax3.set_ylabel('Median Dusk Region Power (GW)', fontsize=22)
    
    #ax3.set_xlim(0.00002, 350)
    #ax3.set_ylim(-15,500)
    
    ax3.set_xscale('log')
    locmaj = matplotlib.ticker.LogLocator(base=10,numticks=30) 
    ax3.xaxis.set_major_locator(locmaj)
    ax3.tick_params(axis='x',which='minor',direction='in')
    ax3.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax3.minorticks_on()
    ax3.tick_params(which='major',direction='in',bottom=True, top=True, left=True, right=True, labelsize=20, length=10)
    ax3.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True, length=5)
    ax3.legend(framealpha=0.5,fontsize=18,loc ="upper right")
    
    
    ax4 = plt.subplot(4,2,4)
    
    ax4.plot(x_new, y_new_noon, '--', color='black', label=f'Noon: {noon_g_fit}', zorder=10)
    
    # plus error
    ax4.scatter(gersh_medians_p[0:2]/1e9, noon_medians[0:2],s=100,color='thistle',marker="^")
    ax4.scatter(gersh_medians_p[2]/1e9, noon_medians[2],s=100,color='thistle')
    ax4.scatter(gersh_medians_p[3:5]/1e9, noon_medians[3:5],s=100,color='thistle',marker="^")
    ax4.scatter(gersh_medians_p[5:9]/1e9, noon_medians[5:9],s=100,color='thistle')
    ax4.scatter(gersh_medians_p[9]/1e9, noon_medians[9],s=100,color='thistle',marker="v")
    ax4.scatter(gersh_medians_p[10]/1e9, noon_medians[10],s=100,color='thistle')
    ax4.scatter(gersh_medians_p[11]/1e9, noon_medians[11],s=100,color='thistle',marker="v")
    ax4.scatter(gersh_medians_p[12:18]/1e9, noon_medians[12:18],s=100,color='thistle',marker=",")
    ax4.scatter(gersh_medians_p[18:21]/1e9, noon_medians[18:21],s=100,color='thistle',marker="v")
    ax4.scatter(gersh_medians_p[21]/1e9, noon_medians[21],s=100,color='thistle',marker="^")
    ax4.scatter(gersh_medians_p[22:]/1e9, noon_medians[22:],s=100,color='thistle',marker=",")

    # minus error
    ax4.scatter(gersh_medians_m[0:2]/1e9, noon_medians[0:2],s=100,color='peachpuff',marker="^")
    ax4.scatter(gersh_medians_m[2]/1e9, noon_medians[2],s=100,color='peachpuff')
    ax4.scatter(gersh_medians_m[3:5]/1e9, noon_medians[3:5],s=100,color='peachpuff',marker="^")
    ax4.scatter(gersh_medians_m[5:9]/1e9, noon_medians[5:9],s=100,color='peachpuff')
    ax4.scatter(gersh_medians_m[9]/1e9, noon_medians[9],s=100,color='peachpuff',marker="v")
    ax4.scatter(gersh_medians_m[10]/1e9, noon_medians[10],s=100,color='peachpuff')
    ax4.scatter(gersh_medians_m[11]/1e9, noon_medians[11],s=100,color='peachpuff',marker="v")
    ax4.scatter(gersh_medians_m[12:18]/1e9, noon_medians[12:18],s=100,color='peachpuff',marker=",")
    ax4.scatter(gersh_medians_m[18:21]/1e9, noon_medians[18:21],s=100,color='peachpuff',marker="v")
    ax4.scatter(gersh_medians_m[21]/1e9, noon_medians[21],s=100,color='peachpuff',marker="^")
    ax4.scatter(gersh_medians_m[22:]/1e9, noon_medians[22:],s=100,color='peachpuff',marker=",")
    
    ax4.errorbar(gersh_medians_p/1e9, noon_medians, yerr=noon_med_errs, xerr=med_gersh_errs_p/1e9, fmt='.', color='thistle')
    
    ax4.plot(x_new, y_new_noon_p, '--',color='thistle', label=f'+ {error}: {noon_g_fit_p}')
    
    ax4.errorbar(gersh_medians_m/1e9, noon_medians, yerr=noon_med_errs, xerr=med_gersh_errs_m/1e9, fmt='.', color='peachpuff')
    
    ax4.plot(x_new, y_new_noon_m, '--',color='peachpuff', label=f'- {error}%: {noon_g_fit_m}')
    
    # main plot
    ax4.scatter(Bz_pos_gersh/1e9, Bz_pos_noonpower, color='red', s=200)#, label='+Bz CME')

    ax4.scatter(Bz_neg_gersh/1e9, Bz_neg_noonpower, color='orange', s=200)#, label='-Bz CME') # all compression 1

    ax4.scatter(By_pos_gersh[0:2]/1e9, By_pos_noonpower[0:2], color='green', s=200)#, label='+By CME') # compression 1 2-5 incl visit 15
    ax4.scatter(By_pos_gersh[2:6]/1e9, By_pos_noonpower[2:6], color='green', s=200, marker=",")#, label='+By CIR')
    ax4.scatter(By_pos_gersh[6:8]/1e9, By_pos_noonpower[6:8], color='green', s=200, marker="v")#, label='+By Rarefraction (Deep)') # rarefraction
    ax4.scatter(By_pos_gersh[8]/1e9, By_pos_noonpower[8], color='green', s=200, marker="^")#, label ='+By Rarefaction (Shallow)')
    ax4.scatter(By_pos_gersh[9:]/1e9, By_pos_noonpower[9:], color='green', s=200, marker=",")

    ax4.scatter(By_neg_gersh[0:4]/1e9, By_neg_noonpower[0:4], color='blue', s=200, marker="^")#, label='-By Rarefaction (Shallow)')
    ax4.scatter(By_neg_gersh[4]/1e9, By_neg_noonpower[4], color='blue', s=200)#, label='-By CME') # 
    ax4.scatter(By_neg_gersh[5]/1e9, By_neg_noonpower[5], color='blue', s=200, marker="v")#, label='-By Rarefaction (Deep)')
    ax4.scatter(By_neg_gersh[6]/1e9, By_neg_noonpower[6], color='blue', s=200)
    ax4.scatter(By_neg_gersh[7]/1e9, By_neg_noonpower[7], color='blue', s=200, marker="v")
    ax4.scatter(By_neg_gersh[8:10]/1e9, By_neg_noonpower[8:10], color='blue', s=200, marker=",")#, label='-By CIR')
    ax4.scatter(By_neg_gersh[10]/1e9, By_neg_noonpower[10], color='blue', s=200, marker="v")#, label='-By CIR (2nd)'
    
    ax4.errorbar(gersh_medians/1e9, noon_medians, yerr=noon_med_errs, xerr=med_gersh_errs/1e9, fmt='.', color='lightgray', zorder=9)
    
    ax4.text(1.5e4,-6,'d',style='italic',fontsize=40) #9.5e3 for all but 1.3e4 for 50%
    
    ax4.set_xlabel('Gershman Reconnection Power (GW)', fontsize=22)
    ax4.set_ylabel('Median Noon Region Power (GW)', fontsize=22)
    
    ax4.set_ylim(-20,440)
    #ax4.set_xlim(0.00002, 350)
    ax4.set_xscale('log')
    locmaj = matplotlib.ticker.LogLocator(base=10,numticks=30) 
    ax4.xaxis.set_major_locator(locmaj)
    ax4.tick_params(axis='x',which='minor',direction='in')
    ax4.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax4.minorticks_on()
    ax4.tick_params(which='major',direction='in',bottom=True, top=True, left=True, right=True, labelsize=20, length=10)
    ax4.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True, length=5)
    
    ax4.legend(framealpha=0.5,fontsize=18,loc ="upper right")
    
    
    # save plot
    saveloc = ('/Users/hannah/OneDrive - Lancaster University/aurora/median_gershman_V_vs_region_power_20.jpg')
    #saveloc = (f'/Users/hannah/OneDrive - Lancaster University/aurora/median_gershman_V_vs_region_power_{error}.jpg')
    plt.savefig(saveloc,bbox_inches='tight',dpi=400)
    
    
if plotting == 'KH_dawn':  
    
    x_new = np.linspace(min(KH_dawn_medians_minus/1e9),max(KH_dawn_medians/1e9), len(KH_dawn_medians))
    
    a, acov = np.polyfit(np.log10(KH_dawn_medians/1e9), medians_swirl, 1, cov=True)
    aa = np.poly1d(a)
    # print("Swirl Fit Gradient:")
    # print(a[0])
    # print("Error of Swirl Fit:")
    # print(np.sqrt(np.diag(acov)))
    y_regress_swirl = a[0]*np.log10(KH_dawn_medians/1e9)+a[1]
    
    print(f"Testing '{plotting}'")
    print('Swirl Region Fit')
    ss_res = np.sum((medians_swirl - y_regress_swirl)**2)
    ss_tot = np.sum((medians_swirl - np.mean(medians_swirl))**2)
    r_squared_s = 1 - (ss_res / ss_tot)
    print(f"R squared: {r_squared_s:.4f}")
    
    adj_s = 1-r_squared_s
    adjr_s = 1 - ((adj_s*23)/22)
    print(f"Adjusted R squared: {adjr_s:.4f}")
    
    swirl_grad = round(a[0],2)
    swirl_grad_err = np.sqrt(np.diag(acov))
    swirl_grad_err = round(swirl_grad_err[0], 2)
    
    swirl_g_fit = (f'{swirl_grad} ± {swirl_grad_err}')


    a_p,acov_p = np.polyfit(np.log10(KH_dawn_medians_plus/1e9), medians_swirl,1, cov=True)
    aa_p = np.poly1d(a_p)
    # print(f"Swirl Fit +{error} Gradient:")
    # print(a_p[0])
    # print(f"Error of Swirl Fit +{error}:")
    # print(np.sqrt(np.diag(acov_p)))
    y_regress_swirl_p = a_p[0]*np.log10(KH_dawn_medians_plus/1e9)+a_p[1]
    
    print("---------------")
    print('Swirl Region +20% Fit')
    ssp_res = np.sum((medians_swirl - y_regress_swirl_p)**2)
    ssp_tot = np.sum((medians_swirl - np.mean(medians_swirl))**2)
    r_squared_sp = 1 - (ssp_res / ssp_tot)
    print(f"R squared: {r_squared_sp:.4f}")
    
    adj_sp = 1-r_squared_sp
    adjr_sp = 1 - ((adj_sp*23)/22)
    print(f"Adjusted R squared: {adjr_sp:.4f}")
    
    swirl_grad_p = round(a_p[0],2)
    swirl_grad_err_p = np.sqrt(np.diag(acov_p))
    swirl_grad_err_p = round(swirl_grad_err_p[0], 2)
    
    swirl_g_fit_p = (f'{swirl_grad_p} ± {swirl_grad_err_p}')
    
    
    a_m,acov_m = np.polyfit(np.log10(KH_dawn_medians_minus/1e9), medians_swirl,1, cov=True)
    aa_m = np.poly1d(a_m)
    # print(f"Swirl Fit -{error} Gradient:")
    # print(a_m[0])
    # print(f"Error of Swirl Fit -{error}:")
    # print(np.sqrt(np.diag(acov_m)))
    y_regress_swirl_m = a_m[0]*np.log10(KH_dawn_medians_minus/1e9)+a_m[1]
    
    print("---------------")
    print('Swirl Region -20% Fit')
    ssm_res = np.sum((medians_swirl - y_regress_swirl_m)**2)
    ssm_tot = np.sum((medians_swirl - np.mean(medians_swirl))**2)
    r_squared_sm = 1 - (ssm_res / ssm_tot)
    print(f"R squared: {r_squared_sm:.4f}")
    
    adj_sm = 1-r_squared_sm
    adjr_sm = 1 - ((adj_sm*23)/22)
    print(f"Adjusted R squared: {adjr_sm:.4f}")
    
    swirl_grad_m = round(a_m[0],2)
    swirl_grad_err_m = np.sqrt(np.diag(acov_m))
    swirl_grad_err_m = round(swirl_grad_err_m[0], 2)
    
    swirl_g_fit_m = (f'{swirl_grad_m} ± {swirl_grad_err_m}')
   
    
    
    b, bcov = np.polyfit(np.log10(KH_dawn_medians/1e9), dusk_medians, 1, cov=True)
    bb = np.poly1d(b)
    # print("Dusk Fit Gradient:")
    # print(b[0])
    # print("Error of Dusk Fit:")
    # print(np.sqrt(np.diag(bcov)))
    y_regress_dusk = b[0]*np.log10(KH_dawn_medians/1e9)+b[1]
    
    print("---------------")
    print('Dusk Region Fit')
    sd_res = np.sum((dusk_medians - y_regress_dusk)**2)
    sd_tot = np.sum((dusk_medians - np.mean(dusk_medians))**2)
    r_squared_d = 1 - (sd_res / sd_tot)
    print(f"R squared: {r_squared_d:.4f}")
    
    adj_d = 1-r_squared_d
    adjr_d = 1 - ((adj_d*23)/22)
    print(f"Adjusted R squared: {adjr_d:.4f}")
    
    dusk_grad = round(b[0],2)
    dusk_grad_err = np.sqrt(np.diag(bcov))
    dusk_grad_err = round(dusk_grad_err[0], 2)
    
    dusk_g_fit = (f'{dusk_grad} ± {dusk_grad_err}')
    
    
    b_p, bcov_p = np.polyfit(np.log10(KH_dawn_medians_plus/1e9), dusk_medians, 1, cov=True)
    bb_p = np.poly1d(b_p)
    # print(f"Dusk Fit +{error} Gradient:")
    # print(b_p[0])
    # print(f"Error of Dusk Fit +{error}:")
    # print(np.sqrt(np.diag(bcov_p)))
    y_regress_dusk_p = b_p[0]*np.log10(KH_dawn_medians_plus/1e9)+b_p[1]
    
    print("---------------")
    print('Dusk Region +20% Fit')
    sdp_res = np.sum((dusk_medians - y_regress_dusk_p)**2)
    sdp_tot = np.sum((dusk_medians - np.mean(dusk_medians))**2)
    r_squared_dp = 1 - (sdp_res / sdp_tot)
    print(f"R squared: {r_squared_dp:.4f}")
    
    adj_dp = 1-r_squared_dp
    adjr_dp = 1 - ((adj_dp*23)/22)
    print(f"Adjusted R squared: {adjr_dp:.4f}")
    
    dusk_grad_p = round(b_p[0],2)
    dusk_grad_err_p = np.sqrt(np.diag(bcov_p))
    dusk_grad_err_p = round(dusk_grad_err_p[0], 2)
    
    dusk_g_fit_p = (f'{dusk_grad_p} ± {dusk_grad_err_p}')
    
    
    b_m, bcov_m = np.polyfit(np.log10(KH_dawn_medians_minus/1e9), dusk_medians, 1, cov=True)
    bb_m = np.poly1d(b_m)
    # print(f"Dusk Fit -{error} Gradient:")
    # print(b_m[0])
    # print(f"Error of Dusk Fit -{error}:")
    # print(np.sqrt(np.diag(bcov_m)))
    y_regress_dusk_m = b_m[0]*np.log10(KH_dawn_medians_minus/1e9)+b_m[1]
    
    print("---------------")
    print('Dusk Region -20% Fit')
    sdm_res = np.sum((dusk_medians - y_regress_dusk_m)**2)
    sdm_tot = np.sum((dusk_medians - np.mean(dusk_medians))**2)
    r_squared_dm = 1 - (sdm_res / sdm_tot)
    print(f"R squared: {r_squared_dm:.4f}")
    
    adj_dm = 1-r_squared_dm
    adjr_dm = 1 - ((adj_dm*23)/22)
    print(f"Adjusted R squared: {adjr_dm:.4f}")
    
    dusk_grad_m = round(b_m[0],2)
    dusk_grad_err_m = np.sqrt(np.diag(bcov_m))
    dusk_grad_err_m = round(dusk_grad_err_m[0], 2)
    
    dusk_g_fit_m = (f'{dusk_grad_m} ± {dusk_grad_err_m}')
    
    
    
    c, ccov = np.polyfit(np.log10(KH_dawn_medians/1e9), noon_medians, 1, cov=True)
    cc = np.poly1d(c)
    # print("Noon Fit Gradient:")
    # print(c[0])
    # print("Error of Noon Fit:")
    # print(np.sqrt(np.diag(ccov)))
    y_regress_noon = c[0]*np.log10(KH_dawn_medians/1e9)+c[1]
    
    print("---------------")
    print('Noon Region Fit')
    sn_res = np.sum((noon_medians - y_regress_noon)**2)
    sn_tot = np.sum((noon_medians - np.mean(noon_medians))**2)
    r_squared_n = 1 - (sn_res / sn_tot)
    print(f"R squared: {r_squared_n:.4f}")
    
    adj_n = 1-r_squared_n
    adjr_n = 1 - ((adj_n*23)/22)
    print(f"Adjusted R squared: {adjr_n:.4f}")
    
    noon_grad = round(c[0],2)
    noon_grad_err = np.sqrt(np.diag(ccov))
    noon_grad_err = round(noon_grad_err[0], 2)
    
    noon_g_fit = (f'{noon_grad} ± {noon_grad_err}')
    
    
    c_p, ccov_p = np.polyfit(np.log10(KH_dawn_medians_plus/1e9), noon_medians, 1, cov=True)
    cc_p = np.poly1d(c_p)
    # print(f"Noon Fit +{error} Gradient:")
    # print(c_p[0])
    # print(f"Error of Noon Fit +{error}:")
    # print(np.sqrt(np.diag(ccov_p)))
    y_regress_noon_p = c_p[0]*np.log10(KH_dawn_medians_plus/1e9)+c_p[1]
    
    print("---------------")
    print('Noon Region +20% Fit')
    snp_res = np.sum((noon_medians - y_regress_noon_p)**2)
    snp_tot = np.sum((noon_medians - np.mean(noon_medians))**2)
    r_squared_np = 1 - (snp_res / snp_tot)
    print(f"R squared: {r_squared_np:.4f}")
    
    adj_np = 1-r_squared_np
    adjr_np = 1 - ((adj_np*23)/22)
    print(f"Adjusted R squared: {adjr_np:.4f}")
    
    noon_grad_p = round(c_p[0],2)
    noon_grad_err_p = np.sqrt(np.diag(ccov_p))
    noon_grad_err_p = round(noon_grad_err_p[0], 2)
    
    noon_g_fit_p = (f'{noon_grad_p} ± {noon_grad_err_p}')
    
    
    c_m, ccov_m = np.polyfit(np.log10(KH_dawn_medians_minus/1e9), noon_medians, 1, cov=True)
    cc_m = np.poly1d(c_m)
    # print(f"Noon Fit -{error} Gradient:")
    # print(c_m[0])
    # print(f"Error of Noon Fit -{error}:")
    # print(np.sqrt(np.diag(ccov_m)))
    y_regress_noon_m = c_m[0]*np.log10(KH_dawn_medians_minus/1e9)+c_m[1]
    
    print("---------------")
    print('Noon Region -20% Fit')
    snm_res = np.sum((noon_medians - y_regress_noon_m)**2)
    snm_tot = np.sum((noon_medians - np.mean(noon_medians))**2)
    r_squared_nm = 1 - (snm_res / snm_tot)
    print(f"R squared: {r_squared_nm:.4f}")
    
    adj_nm = 1-r_squared_nm
    adjr_nm = 1 - ((adj_nm*23)/22)
    print(f"Adjusted R squared: {adjr_nm:.4f}")
    
    noon_grad_m = round(c_m[0],2)
    noon_grad_err_m = np.sqrt(np.diag(ccov_m))
    noon_grad_err_m = round(noon_grad_err_m[0], 2)
    
    noon_g_fit_m = (f'{noon_grad_m} ± {noon_grad_err_m}')
    
    
    y_new_swirl = a[0]*np.log10(x_new)+a[1]
    y_new_dusk = b[0]*np.log10(x_new)+b[1]
    y_new_noon = c[0]*np.log10(x_new)+c[1]
    
    y_new_swirl_p = a_p[0]*np.log10(x_new)+a_p[1]
    y_new_dusk_p = b_p[0]*np.log10(x_new)+b_p[1]
    y_new_noon_p = c_p[0]*np.log10(x_new)+c_p[1]
    
    y_new_swirl_m = a_m[0]*np.log10(x_new)+a_m[1]
    y_new_dusk_m = b_m[0]*np.log10(x_new)+b_m[1]
    y_new_noon_m = c_m[0]*np.log10(x_new)+c_m[1]
    
    
      
    fig = plt.figure(figsize=(28,50))
    ax1 = plt.subplot(4,2,1)
    plt.subplots_adjust(hspace=0.1, wspace=0.15)
    
    # main 
    ax1.scatter(Bz_pos_KH_dawn[0:2]/1e9, Bz_pos_power[0:2], color='red', s=200, label='+Bz CME')
    ax1.scatter(Bz_pos_KH_dawn[2:]/1e9, Bz_pos_power[2:], color='red', s=200, marker='v', label='+Bz Rarefaction (Deep)')

    ax1.scatter(Bz_neg_KH_dawn[0]/1e9, Bz_neg_power[0], color='orange', s=200, label='-Bz CME') #
    ax1.scatter(Bz_neg_KH_dawn[1]/1e9, Bz_neg_power[1], color='orange', s=200,marker=',', label='-Bz CIR') 

    ax1.scatter(By_pos_KH_dawn[0]/1e9, By_pos_power[0], color='green', s=200, label='+By CME') # compression 1 2-5 incl visit 15
    ax1.scatter(By_pos_KH_dawn[1:4]/1e9, By_pos_power[1:4], color='green', s=200, marker=",", label='+By CIR')
    ax1.scatter(By_pos_KH_dawn[4:6]/1e9, By_pos_power[4:6], color='green', s=200, marker="v", label='+By Rarefaction (Deep)') # rarefraction
    ax1.scatter(By_pos_KH_dawn[6]/1e9, By_pos_power[6], color='green', s=200, marker="^", label ='+By Rarefaction (Shallow)')
    ax1.scatter(By_pos_KH_dawn[7:]/1e9, By_pos_power[7:], color='green', s=200, marker=",")

    ax1.scatter(By_neg_KH_dawn[0:4]/1e9, By_neg_power[0:4], color='blue', s=200, marker="^", label='-By Rarefaction (Shallow)')
    ax1.scatter(By_neg_KH_dawn[4:6]/1e9, By_neg_power[4:6], color='blue', s=200, label='-By CME') # 
    ax1.scatter(By_neg_KH_dawn[6]/1e9, By_neg_power[6], color='blue', s=200, marker="v", label='-By Rarefaction (Deep)')
    ax1.scatter(By_neg_KH_dawn[7:]/1e9, By_neg_power[7:], color='blue', s=200, marker=",", label='-By CIR')
    
    ax1.errorbar(KH_dawn_medians/1e9, polar_medians, yerr=polar_med_errs, xerr=med_KH_dawn_errs/1e9, fmt='.', color='lightgray')
    
    ax1.text(4.8e4,6.5,'a',style='italic',fontsize=40) #1.4 RHS, 0.05 LHS (1 is top)
    
    ax1.set_xlabel('Kelvin Helmholtz Reconnection Power (GW)', fontsize=22)
    ax1.set_ylabel('Median Total Polar Power (GW)', fontsize=22)
    
    ax1.set_xlim(1.5e3, 6e4)

    ax1.set_xscale('log')
    ax1.xaxis.set_major_locator(LogLocator(base=10.0))
    ax1.xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10), numticks=100))

    ax1.minorticks_on()
    ax1.set_ylim(-25,950)
    #ax1.set_xscale('log')
    ax1.tick_params(which='major',direction='in',bottom=True, top=True, left=True, right=True, labelsize=20, length=10)
    ax1.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True, labelsize=10, length=5)
    lgnd = plt.legend(framealpha=0.5, labelspacing=0.3, fontsize=18, loc="upper left") #if using equation legend
    for handle in lgnd.legend_handles:
        handle.set_sizes([100.0])
    
    
    
    ax2 = plt.subplot(4,2,2)
    
    ax2.plot(x_new, y_new_swirl, '--',color='black', linewidth=3, markersize=12,label=(f'Swirl $R^2$: {adjr_s:.2g}'), zorder=10)
    
    if error_plot == '20' or error_plot == '10' or error_plot == '50':
        ax2.plot(x_new, y_new_swirl_p, '--',color='darkviolet',linewidth=3, markersize=12, label=f'+ {error}%: {adjr_sp:.2g}')
        ax2.plot(x_new, y_new_swirl_m, '--',color='sandybrown', linewidth=3, markersize=12,label=f'- {error}%: {adjr_sm:.2g}')
        
        # plus error
        ax2.scatter(KH_dawn_medians_plus[0:2]/1e9, medians_swirl[0:2],s=100,color='plum',marker="^")
        ax2.scatter(KH_dawn_medians_plus[2]/1e9, medians_swirl[2],s=100,color='plum')
        ax2.scatter(KH_dawn_medians_plus[3:5]/1e9, medians_swirl[3:5],s=100,color='plum',marker="^")
        ax2.scatter(KH_dawn_medians_plus[5:9]/1e9, medians_swirl[5:9],s=100,color='plum')
        ax2.scatter(KH_dawn_medians_plus[9]/1e9, medians_swirl[9],s=100,color='plum',marker="v")
        ax2.scatter(KH_dawn_medians_plus[10]/1e9, medians_swirl[10],s=100,color='plum')
        ax2.scatter(KH_dawn_medians_plus[11]/1e9, medians_swirl[11],s=100,color='plum',marker="v")
        ax2.scatter(KH_dawn_medians_plus[12:18]/1e9, medians_swirl[12:18],s=100,color='plum',marker=",")
        ax2.scatter(KH_dawn_medians_plus[18:21]/1e9, medians_swirl[18:21],s=100,color='plum',marker="v")
        ax2.scatter(KH_dawn_medians_plus[21]/1e9, medians_swirl[21],s=100,color='plum',marker="^")
        ax2.scatter(KH_dawn_medians_plus[22:]/1e9, medians_swirl[22:],s=100,color='plum',marker=",")
    
        # minus error
        ax2.scatter(KH_dawn_medians_minus[0:2]/1e9, medians_swirl[0:2],s=100,color='peachpuff',marker="^")
        ax2.scatter(KH_dawn_medians_minus[2]/1e9, medians_swirl[2],s=100,color='peachpuff')
        ax2.scatter(KH_dawn_medians_minus[3:5]/1e9, medians_swirl[3:5],s=100,color='peachpuff',marker="^")
        ax2.scatter(KH_dawn_medians_minus[5:9]/1e9, medians_swirl[5:9],s=100,color='peachpuff')
        ax2.scatter(KH_dawn_medians_minus[9]/1e9, medians_swirl[9],s=100,color='peachpuff',marker="v")
        ax2.scatter(KH_dawn_medians_minus[10]/1e9, medians_swirl[10],s=100,color='peachpuff')
        ax2.scatter(KH_dawn_medians_minus[11]/1e9, medians_swirl[11],s=100,color='peachpuff',marker="v")
        ax2.scatter(KH_dawn_medians_minus[12:18]/1e9, medians_swirl[12:18],s=100,color='peachpuff',marker=",")
        ax2.scatter(KH_dawn_medians_minus[18:21]/1e9, medians_swirl[18:21],s=100,color='peachpuff',marker="v")
        ax2.scatter(KH_dawn_medians_minus[21]/1e9, medians_swirl[21],s=100,color='peachpuff',marker="^")
        ax2.scatter(KH_dawn_medians_minus[22:]/1e9, medians_swirl[22:],s=100,color='peachpuff',marker=",")
        
        ax2.errorbar(KH_dawn_medians_plus/1e9, medians_swirl, yerr=swirl_med_errs, xerr=med_KH_dawn_errs_plus/1e9, fmt='.', color='plum')
        
        ax2.errorbar(KH_dawn_medians_minus/1e9, medians_swirl, yerr=swirl_med_errs, xerr=med_KH_dawn_errs_minus/1e9, fmt='.', color='peachpuff')
    
    # main plot
    ax2.scatter(Bz_pos_KH_dawn[0:2]/1e9, Bz_pos_swirlpower[0:2], color='red', s=200,zorder=9)
    ax2.scatter(Bz_pos_KH_dawn[2:]/1e9, Bz_pos_swirlpower[2:], color='red', s=200, marker='v',zorder=9)

    ax2.scatter(Bz_neg_KH_dawn[0]/1e9, Bz_neg_swirlpower[0], color='orange', s=200,zorder=9) #
    ax2.scatter(Bz_neg_KH_dawn[1]/1e9, Bz_neg_swirlpower[1], color='orange', s=200,marker=',',zorder=9) 

    ax2.scatter(By_pos_KH_dawn[0]/1e9, By_pos_swirlpower[0], color='green', s=200,zorder=9) # compression 1 2-5 incl visit 15
    ax2.scatter(By_pos_KH_dawn[1:4]/1e9, By_pos_swirlpower[1:4], color='green', s=200, marker=",",zorder=9)
    ax2.scatter(By_pos_KH_dawn[4:6]/1e9, By_pos_swirlpower[4:6], color='green', s=200, marker="v",zorder=9) # rarefraction
    ax2.scatter(By_pos_KH_dawn[6]/1e9, By_pos_swirlpower[6], color='green', s=200, marker="^",zorder=9)
    ax2.scatter(By_pos_KH_dawn[7:]/1e9, By_pos_swirlpower[7:], color='green', s=200, marker=",",zorder=9)

    ax2.scatter(By_neg_KH_dawn[0:4]/1e9, By_neg_swirlpower[0:4], color='blue', s=200, marker="^",zorder=9)
    ax2.scatter(By_neg_KH_dawn[4:6]/1e9, By_neg_swirlpower[4:6], color='blue', s=200,zorder=9) # 
    ax2.scatter(By_neg_KH_dawn[6]/1e9, By_neg_swirlpower[6], color='blue', s=200, marker="v",zorder=9)
    ax2.scatter(By_neg_KH_dawn[7:]/1e9, By_neg_swirlpower[7:], color='blue', s=200, marker=",",zorder=9)
    
    ax2.errorbar(KH_dawn_medians/1e9, medians_swirl, yerr=swirl_med_errs, xerr=med_KH_dawn_errs/1e9, fmt='.', color='lightgray', zorder=9)
    
    ax2.text(4.8e4,6.5,'b',style='italic',fontsize=40) #1070 no error, 1100 error
    
    ax2.set_xlabel('Kelvin Helmholtz Reconnection Power (GW)', fontsize=22)
    ax2.set_ylabel('Median Swirl Region Power (GW)', fontsize=22)

    ax2.set_xlim(1.5e3, 6e4)
    #ax2.set_ylim(-5,120)
    
    ax2.set_xscale('log')
    locmaj = matplotlib.ticker.LogLocator(base=10,numticks=30) 
    ax2.xaxis.set_major_locator(locmaj)
    ax2.tick_params(axis='x',which='minor',direction='in')
    ax2.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax2.minorticks_on()
    ax2.tick_params(which='major',direction='in',bottom=True, top=True, left=True, right=True, labelsize=20, length=10)
    ax2.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True, length=5)
    ax2.legend(framealpha=0.5,fontsize=21,loc ="upper left")
    
    
    ax3 = plt.subplot(4,2,3)
    
    ax3.plot(x_new, y_new_dusk,'--', color='black',linewidth=3, markersize=12, label=f'Dusk $R^2$: {adjr_d:.2g}', zorder=10)
    
    if error_plot == '20' or error_plot == '10' or error_plot == '50':
        ax3.plot(x_new, y_new_dusk_p, '--',color='darkviolet',linewidth=3, markersize=12, label=f'+ {error}%: {adjr_dp:.2g}')
        ax3.plot(x_new, y_new_dusk_m, '--',color='sandybrown', linewidth=3, markersize=12,label=f'- {error}%: {adjr_dm:.2g}')
        
        # plus error
        ax3.scatter(KH_dawn_medians_plus[0:2]/1e9, dusk_medians[0:2],s=100,color='plum',marker="^")
        ax3.scatter(KH_dawn_medians_plus[2]/1e9, dusk_medians[2],s=100,color='plum')
        ax3.scatter(KH_dawn_medians_plus[3:5]/1e9, dusk_medians[3:5],s=100,color='plum',marker="^")
        ax3.scatter(KH_dawn_medians_plus[5:9]/1e9, dusk_medians[5:9],s=100,color='plum')
        ax3.scatter(KH_dawn_medians_plus[9]/1e9, dusk_medians[9],s=100,color='plum',marker="v")
        ax3.scatter(KH_dawn_medians_plus[10]/1e9, dusk_medians[10],s=100,color='plum')
        ax3.scatter(KH_dawn_medians_plus[11]/1e9, dusk_medians[11],s=100,color='plum',marker="v")
        ax3.scatter(KH_dawn_medians_plus[12:18]/1e9, dusk_medians[12:18],s=100,color='plum',marker=",")
        ax3.scatter(KH_dawn_medians_plus[18:21]/1e9, dusk_medians[18:21],s=100,color='plum',marker="v")
        ax3.scatter(KH_dawn_medians_plus[21]/1e9, dusk_medians[21],s=100,color='plum',marker="^")
        ax3.scatter(KH_dawn_medians_plus[22:]/1e9, dusk_medians[22:],s=100,color='plum',marker=",")
    
        # minus error
        ax3.scatter(KH_dawn_medians_minus[0:2]/1e9, dusk_medians[0:2],s=100,color='peachpuff',marker="^")
        ax3.scatter(KH_dawn_medians_minus[2]/1e9, dusk_medians[2],s=100,color='peachpuff')
        ax3.scatter(KH_dawn_medians_minus[3:5]/1e9, dusk_medians[3:5],s=100,color='peachpuff',marker="^")
        ax3.scatter(KH_dawn_medians_minus[5:9]/1e9, dusk_medians[5:9],s=100,color='peachpuff')
        ax3.scatter(KH_dawn_medians_minus[9]/1e9, dusk_medians[9],s=100,color='peachpuff',marker="v")
        ax3.scatter(KH_dawn_medians_minus[10]/1e9, dusk_medians[10],s=100,color='peachpuff')
        ax3.scatter(KH_dawn_medians_minus[11]/1e9, dusk_medians[11],s=100,color='peachpuff',marker="v")
        ax3.scatter(KH_dawn_medians_minus[12:18]/1e9, dusk_medians[12:18],s=100,color='peachpuff',marker=",")
        ax3.scatter(KH_dawn_medians_minus[18:21]/1e9, dusk_medians[18:21],s=100,color='peachpuff',marker="v")
        ax3.scatter(KH_dawn_medians_minus[21]/1e9, dusk_medians[21],s=100,color='peachpuff',marker="^")
        ax3.scatter(KH_dawn_medians_minus[22:]/1e9, dusk_medians[22:],s=100,color='peachpuff',marker=",")
        
        ax3.errorbar(KH_dawn_medians_plus/1e9, dusk_medians, yerr=dusk_med_errs, xerr=med_KH_dawn_errs_plus/1e9, fmt='.', color='plum')
        
        ax3.errorbar(KH_dawn_medians_minus/1e9, dusk_medians, yerr=dusk_med_errs, xerr=med_KH_dawn_errs_minus/1e9, fmt='.', color='peachpuff')
    
    # main plot
    ax3.scatter(Bz_pos_KH_dawn[0:2]/1e9, Bz_pos_duskpower[0:2], color='red', s=200,zorder=9)
    ax3.scatter(Bz_pos_KH_dawn[2:]/1e9, Bz_pos_duskpower[2:], color='red', s=200, marker='v',zorder=9)

    ax3.scatter(Bz_neg_KH_dawn[0]/1e9, Bz_neg_duskpower[0], color='orange', s=200,zorder=9) #
    ax3.scatter(Bz_neg_KH_dawn[1]/1e9, Bz_neg_duskpower[1], color='orange', s=200,marker=',',zorder=9)

    ax3.scatter(By_pos_KH_dawn[0]/1e9, By_pos_duskpower[0], color='green', s=200,zorder=9) # compression 1 2-5 incl visit 15
    ax3.scatter(By_pos_KH_dawn[1:4]/1e9, By_pos_duskpower[1:4], color='green', s=200, marker=",",zorder=9)
    ax3.scatter(By_pos_KH_dawn[4:6]/1e9, By_pos_duskpower[4:6], color='green', s=200, marker="v",zorder=9) # rarefraction
    ax3.scatter(By_pos_KH_dawn[6]/1e9, By_pos_duskpower[6], color='green', s=200, marker="^",zorder=9)
    ax3.scatter(By_pos_KH_dawn[7:]/1e9, By_pos_duskpower[7:], color='green', s=200, marker=",",zorder=9)

    ax3.scatter(By_neg_KH_dawn[0:4]/1e9, By_neg_duskpower[0:4], color='blue', s=200, marker="^",zorder=9)
    ax3.scatter(By_neg_KH_dawn[4:6]/1e9, By_neg_duskpower[4:6], color='blue', s=200,zorder=9) # 
    ax3.scatter(By_neg_KH_dawn[6]/1e9, By_neg_duskpower[6], color='blue', s=200, marker="v",zorder=9)
    ax3.scatter(By_neg_KH_dawn[7:]/1e9, By_neg_duskpower[7:], color='blue', s=200, marker=",",zorder=9)
    
    ax3.errorbar(KH_dawn_medians/1e9, dusk_medians, yerr=dusk_med_errs, xerr=med_KH_dawn_errs/1e9, fmt='.', color='lightgray', zorder=9)
    
    ax3.text(4.8e4,-13,'c',style='italic',fontsize=40)#1070 no error, 1100 error
    
    ax3.set_xlabel('Kelvin Helmholtz Reconnection Power (GW)', fontsize=22)
    ax3.set_ylabel('Median Dusk Region Power (GW)', fontsize=22)
    
    ax3.set_xlim(1.5e3, 6e4)
    #ax3.set_xlim(0.00002, 350)
    #ax3.set_ylim(-15,500)
    ax3.set_xscale('log')
    locmaj = matplotlib.ticker.LogLocator(base=10,numticks=30) 
    ax3.xaxis.set_major_locator(locmaj)
    ax3.tick_params(axis='x',which='minor',direction='in')
    ax3.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax3.minorticks_on()
    ax3.tick_params(which='major',direction='in',bottom=True, top=True, left=True, right=True, labelsize=20, length=10)
    ax3.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True, length=5)

    ax3.legend(framealpha=0.5,fontsize=21,loc ="upper left")
    
    
    ax4 = plt.subplot(4,2,4)
    
    ax4.plot(x_new, y_new_noon, '--',color='black',linewidth=3, markersize=12, label=f'Noon $R^2$: {adjr_n:.2g}', zorder=10)
    
    if error_plot == '20' or error_plot == '10' or error_plot == '50':
        ax4.plot(x_new, y_new_noon_p, '--',color='darkviolet',linewidth=3, markersize=12, label=f'+ {error}%: {adjr_np:.2g}')
        ax4.plot(x_new, y_new_noon_m, '--',color='sandybrown',linewidth=3, markersize=12, label=f'- {error}%: {adjr_nm:.2g}')
        
        # plus error
        ax4.scatter(KH_dawn_medians_plus[0:2]/1e9, noon_medians[0:2],s=100,color='plum',marker="^")
        ax4.scatter(KH_dawn_medians_plus[2]/1e9, noon_medians[2],s=100,color='plum')
        ax4.scatter(KH_dawn_medians_plus[3:5]/1e9, noon_medians[3:5],s=100,color='plum',marker="^")
        ax4.scatter(KH_dawn_medians_plus[5:9]/1e9, noon_medians[5:9],s=100,color='plum')
        ax4.scatter(KH_dawn_medians_plus[9]/1e9, noon_medians[9],s=100,color='plum',marker="v")
        ax4.scatter(KH_dawn_medians_plus[10]/1e9, noon_medians[10],s=100,color='plum')
        ax4.scatter(KH_dawn_medians_plus[11]/1e9, noon_medians[11],s=100,color='plum',marker="v")
        ax4.scatter(KH_dawn_medians_plus[12:18]/1e9, noon_medians[12:18],s=100,color='plum',marker=",")
        ax4.scatter(KH_dawn_medians_plus[18:21]/1e9, noon_medians[18:21],s=100,color='plum',marker="v")
        ax4.scatter(KH_dawn_medians_plus[21]/1e9, noon_medians[21],s=100,color='plum',marker="^")
        ax4.scatter(KH_dawn_medians_plus[22:]/1e9, noon_medians[22:],s=100,color='plum',marker=",")
    
        # minus error
        ax4.scatter(KH_dawn_medians_minus[0:2]/1e9, noon_medians[0:2],s=100,color='peachpuff',marker="^")
        ax4.scatter(KH_dawn_medians_minus[2]/1e9, noon_medians[2],s=100,color='peachpuff')
        ax4.scatter(KH_dawn_medians_minus[3:5]/1e9, noon_medians[3:5],s=100,color='peachpuff',marker="^")
        ax4.scatter(KH_dawn_medians_minus[5:9]/1e9, noon_medians[5:9],s=100,color='peachpuff')
        ax4.scatter(KH_dawn_medians_minus[9]/1e9, noon_medians[9],s=100,color='peachpuff',marker="v")
        ax4.scatter(KH_dawn_medians_minus[10]/1e9, noon_medians[10],s=100,color='peachpuff')
        ax4.scatter(KH_dawn_medians_minus[11]/1e9, noon_medians[11],s=100,color='peachpuff',marker="v")
        ax4.scatter(KH_dawn_medians_minus[12:18]/1e9, noon_medians[12:18],s=100,color='peachpuff',marker=",")
        ax4.scatter(KH_dawn_medians_minus[18:21]/1e9, noon_medians[18:21],s=100,color='peachpuff',marker="v")
        ax4.scatter(KH_dawn_medians_minus[21]/1e9, noon_medians[21],s=100,color='peachpuff',marker="^")
        ax4.scatter(KH_dawn_medians_minus[22:]/1e9, noon_medians[22:],s=100,color='peachpuff',marker=",")
        
        ax4.errorbar(KH_dawn_medians_plus/1e9, noon_medians, yerr=noon_med_errs, xerr=med_KH_dawn_errs_plus/1e9, fmt='.', color='plum')
        
        ax4.errorbar(KH_dawn_medians_minus/1e9, noon_medians, yerr=noon_med_errs, xerr=med_KH_dawn_errs_minus/1e9, fmt='.', color='peachpuff')
    
    
    # main plot
    ax4.scatter(Bz_pos_KH_dawn[0:2]/1e9, Bz_pos_noonpower[0:2], color='red', s=200,zorder=9)
    ax4.scatter(Bz_pos_KH_dawn[2:]/1e9, Bz_pos_noonpower[2:], color='red', s=200, marker='v',zorder=9)

    ax4.scatter(Bz_neg_KH_dawn[0]/1e9, Bz_neg_noonpower[0], color='orange', s=200,zorder=9) #
    ax4.scatter(Bz_neg_KH_dawn[1]/1e9, Bz_neg_noonpower[1], color='orange', s=200,marker=',',zorder=9) 

    ax4.scatter(By_pos_KH_dawn[0]/1e9, By_pos_noonpower[0], color='green', s=200,zorder=9) # compression 1 2-5 incl visit 15
    ax4.scatter(By_pos_KH_dawn[1:4]/1e9, By_pos_noonpower[1:4], color='green', s=200, marker=",",zorder=9)
    ax4.scatter(By_pos_KH_dawn[4:6]/1e9, By_pos_noonpower[4:6], color='green', s=200, marker="v",zorder=9) # rarefraction
    ax4.scatter(By_pos_KH_dawn[6]/1e9, By_pos_noonpower[6], color='green', s=200, marker="^",zorder=9)
    ax4.scatter(By_pos_KH_dawn[7:]/1e9, By_pos_noonpower[7:], color='green', s=200, marker=",",zorder=9)

    ax4.scatter(By_neg_KH_dawn[0:4]/1e9, By_neg_noonpower[0:4], color='blue', s=200, marker="^",zorder=9)
    ax4.scatter(By_neg_KH_dawn[4:6]/1e9, By_neg_noonpower[4:6], color='blue', s=200,zorder=9)
    ax4.scatter(By_neg_KH_dawn[6]/1e9, By_neg_noonpower[6], color='blue', s=200, marker="v",zorder=9)
    ax4.scatter(By_neg_KH_dawn[7:]/1e9, By_neg_noonpower[7:], color='blue', s=200, marker=",",zorder=9)
    
    ax4.errorbar(KH_dawn_medians/1e9, noon_medians, yerr=noon_med_errs, xerr=med_KH_dawn_errs/1e9, fmt='.', color='lightgray', zorder=9)
    
    ax4.text(4.8e4,-11,'d',style='italic',fontsize=40) #1070 no error, 1100 error
    
    ax4.set_xlabel('Kelvin Helmholtz Reconnection Power (GW)', fontsize=22)
    ax4.set_ylabel('Median Noon Region Power (GW)', fontsize=22)
    
    ax4.set_xlim(1.5e3, 6e4)
    #ax4.set_ylim(-15,420)
    #ax4.set_xlim(0.00002, 350)
    ax4.set_xscale('log')
    locmaj = matplotlib.ticker.LogLocator(base=10,numticks=30) 
    ax4.xaxis.set_major_locator(locmaj)
    ax4.tick_params(axis='x',which='minor',direction='in')
    ax4.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax4.minorticks_on()
    ax4.tick_params(which='major',direction='in',bottom=True, top=True, left=True, right=True, labelsize=20, length=10)
    ax4.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True, length=5)
    ax4.legend(framealpha=0.5,fontsize=21,loc ="upper left")
    
    
    # save plot
    saveloc = (f'{root_saves}median_KH_dawn_V_vs_region_power_{error_plot}.jpg') 
    #saveloc = (f'/Users/hannah/OneDrive - Lancaster University/aurora/median_KH_dawn_V_vs_region_power_{error}.jpg')
    plt.savefig(saveloc,bbox_inches='tight',dpi=400)
    
    
    
if plotting == 'KH_dusk':  
    
    x_new = np.linspace(min(KH_dusk_medians_plus/1e9),max(KH_dusk_medians/1e9), len(KH_dusk_medians))
    
    a, acov = np.polyfit(np.log10(KH_dusk_medians/1e9), medians_swirl, 1, cov=True)
    aa = np.poly1d(a)
    print("Swirl Fit Gradient:")
    # print(a[0])
    # print("Error of Swirl Fit:")
    # print(np.sqrt(np.diag(acov)))
    y_regress_swirl = a[0]*np.log10(KH_dusk_medians/1e9)+a[1]
    
    print(f"Testing '{plotting}'")
    print('Swirl Region Fit')
    ss_res = np.sum((medians_swirl - y_regress_swirl)**2)
    ss_tot = np.sum((medians_swirl - np.mean(medians_swirl))**2)
    r_squared_s = 1 - (ss_res / ss_tot)
    print(f"R squared: {r_squared_s:.4f}")

    adj_s = 1-r_squared_s
    adjr_s = 1 - ((adj_s*23)/22)
    print(f"Adjusted R squared: {adjr_s:.4f}")
    
    swirl_grad = round(a[0],2)
    swirl_grad_err = np.sqrt(np.diag(acov))
    swirl_grad_err = round(swirl_grad_err[0], 2)
    
    swirl_g_fit = (f'{swirl_grad} ± {swirl_grad_err}')
    
    
    a_p,acov_p = np.polyfit(np.log10(KH_dusk_medians_plus/1e9), medians_swirl,1, cov=True)
    aa_p = np.poly1d(a_p)
    # print(f"Swirl Fit +{error} Gradient:")
    # print(a_p[0])
    # print(f"Error of Swirl Fit +{error}:")
    # print(np.sqrt(np.diag(acov_p)))
    y_regress_swirl_p = a_p[0]*np.log10(KH_dusk_medians_plus/1e9)+a_p[1]
    
    print("---------------")
    print('Swirl Region +20% Fit')
    ssp_res = np.sum((medians_swirl - y_regress_swirl_p)**2)
    ssp_tot = np.sum((medians_swirl - np.mean(medians_swirl))**2)
    r_squared_sp = 1 - (ssp_res / ssp_tot)
    print(f"R squared: {r_squared_sp:.4f}")

    adj_sp = 1-r_squared_sp
    adjr_sp = 1 - ((adj_sp*23)/22)
    print(f"Adjusted R squared: {adjr_sp:.4f}")
    
    swirl_grad_p = round(a_p[0],2)
    swirl_grad_err_p = np.sqrt(np.diag(acov_p))
    swirl_grad_err_p = round(swirl_grad_err_p[0], 2)
    
    swirl_g_fit_p = (f'{swirl_grad_p} ± {swirl_grad_err_p}')
    
    
    a_m,acov_m = np.polyfit(np.log10(KH_dusk_medians_minus/1e9), medians_swirl,1, cov=True)
    aa_m = np.poly1d(a_m)
    # print(f"Swirl Fit -{error} Gradient:")
    # print(a_m[0])
    # print(f"Error of Swirl Fit -{error}:")
    # print(np.sqrt(np.diag(acov_m)))
    y_regress_swirl_m = a_m[0]*np.log10(KH_dusk_medians_minus/1e9)+a_m[1]
    
    print("---------------")
    print('Swirl Region -20% Fit')
    ssm_res = np.sum((medians_swirl - y_regress_swirl_m)**2)
    ssm_tot = np.sum((medians_swirl - np.mean(medians_swirl))**2)
    r_squared_sm = 1 - (ssm_res / ssm_tot)
    print(f"R squared: {r_squared_sm:.4f}")

    adj_sm = 1-r_squared_sm
    adjr_sm = 1 - ((adj_sm*23)/22)
    print(f"Adjusted R squared: {adjr_sm:.4f}")
    
    swirl_grad_m = round(a_m[0],2)
    swirl_grad_err_m = np.sqrt(np.diag(acov_m))
    swirl_grad_err_m = round(swirl_grad_err_m[0], 2)
    
    swirl_g_fit_m = (f'{swirl_grad_m} ± {swirl_grad_err_m}')
   
    
    
    b, bcov = np.polyfit(np.log10(KH_dusk_medians/1e9), dusk_medians, 1, cov=True)
    bb = np.poly1d(b)
    # print("Dusk Fit Gradient:")
    # print(b[0])
    # print("Error of Dusk Fit:")
    # print(np.sqrt(np.diag(bcov)))
    y_regress_dusk = b[0]*np.log10(KH_dusk_medians/1e9)+b[1]
    
    print("---------------")
    print('Dusk Region Fit')
    sd_res = np.sum((dusk_medians - y_regress_dusk)**2)
    sd_tot = np.sum((dusk_medians - np.mean(dusk_medians))**2)
    r_squared_d = 1 - (sd_res / sd_tot)
    print(f"R squared: {r_squared_d:.4f}")

    adj_d = 1-r_squared_d
    adjr_d = 1 - ((adj_d*23)/22)
    print(f"Adjusted R squared: {adjr_d:.4f}")
    
    dusk_grad = round(b[0],2)
    dusk_grad_err = np.sqrt(np.diag(bcov))
    dusk_grad_err = round(dusk_grad_err[0], 2)
    
    dusk_g_fit = (f'{dusk_grad} ± {dusk_grad_err}')
    
    
    b_p, bcov_p = np.polyfit(np.log10(KH_dusk_medians_plus/1e9), dusk_medians, 1, cov=True)
    bb_p = np.poly1d(b_p)
    # print(f"Dusk Fit +{error} Gradient:")
    # print(b_p[0])
    # print(f"Error of Dusk Fit +{error}:")
    # print(np.sqrt(np.diag(bcov_p)))
    y_regress_dusk_p = b_p[0]*np.log10(KH_dusk_medians_plus/1e9)+b_p[1]
    
    print("---------------")
    print('Dusk Region +20% Fit')
    sdp_res = np.sum((dusk_medians - y_regress_dusk_p)**2)
    sdp_tot = np.sum((dusk_medians - np.mean(dusk_medians))**2)
    r_squared_dp = 1 - (sdp_res / sdp_tot)
    print(f"R squared: {r_squared_dp:.4f}")

    adj_dp = 1-r_squared_dp
    adjr_dp = 1 - ((adj_dp*23)/22)
    print(f"Adjusted R squared: {adjr_dp:.4f}")
    
    dusk_grad_p = round(b_p[0],2)
    dusk_grad_err_p = np.sqrt(np.diag(bcov_p))
    dusk_grad_err_p = round(dusk_grad_err_p[0], 2)
    
    dusk_g_fit_p = (f'{dusk_grad_p} ± {dusk_grad_err_p}')
    
    
    b_m, bcov_m = np.polyfit(np.log10(KH_dusk_medians_minus/1e9), dusk_medians, 1, cov=True)
    bb_m = np.poly1d(b_m)
    # print(f"Dusk Fit -{error} Gradient:")
    # print(b_m[0])
    # print(f"Error of Dusk Fit -{error}:")
    # print(np.sqrt(np.diag(bcov_m)))
    y_regress_dusk_m = b_m[0]*np.log10(KH_dusk_medians_minus/1e9)+b_m[1]
    
    print("---------------")
    print('Dusk Region -20% Fit')
    sdm_res = np.sum((dusk_medians - y_regress_dusk_m)**2)
    sdm_tot = np.sum((dusk_medians - np.mean(dusk_medians))**2)
    r_squared_dm = 1 - (sdm_res / sdm_tot)
    print(f"R squared: {r_squared_dm:.4f}")

    adj_dm = 1-r_squared_dm
    adjr_dm = 1 - ((adj_dm*23)/22)
    print(f"Adjusted R squared: {adjr_dm:.4f}")
    
    dusk_grad_m = round(b_m[0],2)
    dusk_grad_err_m = np.sqrt(np.diag(bcov_m))
    dusk_grad_err_m = round(dusk_grad_err_m[0], 2)
    
    dusk_g_fit_m = (f'{dusk_grad_m} ± {dusk_grad_err_m}')
    
    
    
    c, ccov = np.polyfit(np.log10(KH_dusk_medians/1e9), noon_medians, 1, cov=True)
    cc = np.poly1d(c)
    # print("Noon Fit Gradient:")
    # print(c[0])
    # print("Error of Noon Fit:")
    # print(np.sqrt(np.diag(ccov)))
    y_regress_noon = c[0]*np.log10(KH_dusk_medians/1e9)+c[1]
    
    print("---------------")
    print('Noon Region Fit')
    sn_res = np.sum((noon_medians - y_regress_noon)**2)
    sn_tot = np.sum((noon_medians - np.mean(noon_medians))**2)
    r_squared_n = 1 - (sn_res / sn_tot)
    print(f"R squared: {r_squared_n:.4f}")

    adj_n = 1-r_squared_n
    adjr_n = 1 - ((adj_n*23)/22)
    print(f"Adjusted R squared: {adjr_n:.4f}")
    
    noon_grad = round(c[0],2)
    noon_grad_err = np.sqrt(np.diag(ccov))
    noon_grad_err = round(noon_grad_err[0], 2)
    
    noon_g_fit = (f'{noon_grad} ± {noon_grad_err}')
    
    
    c_p, ccov_p = np.polyfit(np.log10(KH_dusk_medians_plus/1e9), noon_medians, 1, cov=True)
    cc_p = np.poly1d(c_p)
    # print(f"Noon Fit +{error} Gradient:")
    # print(c_p[0])
    # print(f"Error of Noon Fit +{error}:")
    # print(np.sqrt(np.diag(ccov_p)))
    y_regress_noon_p = c_p[0]*np.log10(KH_dusk_medians_plus/1e9)+c_p[1]
    
    print("---------------")
    print('Noon Region +20% Fit')
    snp_res = np.sum((noon_medians - y_regress_noon_p)**2)
    snp_tot = np.sum((noon_medians - np.mean(noon_medians))**2)
    r_squared_np = 1 - (snp_res / snp_tot)
    print(f"R squared: {r_squared_np:.4f}")

    adj_np = 1-r_squared_np
    adjr_np = 1 - ((adj_np*23)/22)
    print(f"Adjusted R squared: {adjr_np:.4f}")
    
    noon_grad_p = round(c_p[0],2)
    noon_grad_err_p = np.sqrt(np.diag(ccov_p))
    noon_grad_err_p = round(noon_grad_err_p[0], 2)
    
    noon_g_fit_p = (f'{noon_grad_p} ± {noon_grad_err_p}')
    
    
    c_m, ccov_m = np.polyfit(np.log10(KH_dusk_medians_minus/1e9), noon_medians, 1, cov=True)
    cc_m = np.poly1d(c_m)
    # print(f"Noon Fit -{error} Gradient:")
    # print(c_m[0])
    # print(f"Error of Noon Fit -{error}:")
    # print(np.sqrt(np.diag(ccov_m)))
    y_regress_noon_m = c_m[0]*np.log10(KH_dusk_medians_minus/1e9)+c_m[1]
    
    print("---------------")
    print('Noon Region -20% Fit')
    snm_res = np.sum((noon_medians - y_regress_noon_m)**2)
    snm_tot = np.sum((noon_medians - np.mean(noon_medians))**2)
    r_squared_nm = 1 - (snm_res / snm_tot)
    print(f"R squared: {r_squared_nm:.4f}")

    adj_nm = 1-r_squared_nm
    adjr_nm = 1 - ((adj_nm*23)/22)
    print(f"Adjusted R squared: {adjr_nm:.4f}")
    
    noon_grad_m = round(c_m[0],2)
    noon_grad_err_m = np.sqrt(np.diag(ccov_m))
    noon_grad_err_m = round(noon_grad_err_m[0], 2)
    
    noon_g_fit_m = (f'{noon_grad_m} ± {noon_grad_err_m}')
    
    
    y_new_swirl = a[0]*np.log10(x_new)+a[1]
    y_new_dusk = b[0]*np.log10(x_new)+b[1]
    y_new_noon = c[0]*np.log10(x_new)+c[1]
    
    y_new_swirl_p = a_p[0]*np.log10(x_new)+a_p[1]
    y_new_dusk_p = b_p[0]*np.log10(x_new)+b_p[1]
    y_new_noon_p = c_p[0]*np.log10(x_new)+c_p[1]
    
    y_new_swirl_m = a_m[0]*np.log10(x_new)+a_m[1]
    y_new_dusk_m = b_m[0]*np.log10(x_new)+b_m[1]
    y_new_noon_m = c_m[0]*np.log10(x_new)+c_m[1]
    
    
      
    
    fig = plt.figure(figsize=(28,50))
    ax1 = plt.subplot(4,2,1)
    plt.subplots_adjust(hspace=0.1, wspace=0.15)

    # main plot
    ax1.scatter(Bz_pos_KH_dusk[0:2]/1e9, Bz_pos_power[0:2], color='red', s=200, label='+Bz CME')
    ax1.scatter(Bz_pos_KH_dusk[2:]/1e9, Bz_pos_power[2:], color='red', s=200, marker='v', label='+Bz Rarefaction (Deep)')

    ax1.scatter(Bz_neg_KH_dusk[0]/1e9, Bz_neg_power[0], color='orange', s=200, label='-Bz CME') #
    ax1.scatter(Bz_neg_KH_dusk[1]/1e9, Bz_neg_power[1], color='orange', s=200,marker=',', label='-Bz CIR') 

    ax1.scatter(By_pos_KH_dusk[0]/1e9, By_pos_power[0], color='green', s=200, label='+By CME') # compression 1 2-5 incl visit 15
    ax1.scatter(By_pos_KH_dusk[1:4]/1e9, By_pos_power[1:4], color='green', s=200, marker=",", label='+By CIR')
    ax1.scatter(By_pos_KH_dusk[4:6]/1e9, By_pos_power[4:6], color='green', s=200, marker="v", label='+By Rarefaction (Deep)') # rarefraction
    ax1.scatter(By_pos_KH_dusk[6]/1e9, By_pos_power[6], color='green', s=200, marker="^", label ='+By Rarefaction (Shallow)')
    ax1.scatter(By_pos_KH_dusk[7:]/1e9, By_pos_power[7:], color='green', s=200, marker=",")

    ax1.scatter(By_neg_KH_dusk[0:4]/1e9, By_neg_power[0:4], color='blue', s=200, marker="^", label='-By Rarefaction (Shallow)')
    ax1.scatter(By_neg_KH_dusk[4:6]/1e9, By_neg_power[4:6], color='blue', s=200, label='-By CME') # 
    ax1.scatter(By_neg_KH_dusk[6]/1e9, By_neg_power[6], color='blue', s=200, marker="v", label='-By Rarefaction (Deep)')
    ax1.scatter(By_neg_KH_dusk[7:]/1e9, By_neg_power[7:], color='blue', s=200, marker=",", label='-By CIR')
    
    ax1.errorbar(KH_dusk_medians/1e9, polar_medians, yerr=polar_med_errs, xerr=med_KH_dusk_errs/1e9, fmt='.', color='lightgray')
    
    ax1.text(5.5e2,10,'a',style='italic',fontsize=40) #1.4 RHS, 0.05 LHS (1 is top)
    
    ax1.set_xlabel('Kelvin Helmholtz Reconnection Power (GW)', fontsize=22)
    ax1.set_ylabel('Median Total Polar Power (GW)', fontsize=22)
    
    ax1.set_xlim(0.3e1,8e2)
    # ax1.set_ylim(-20, 970)
    ax1.set_xscale('log')
    ax1.set_ylim(-25,950)
    locmaj = matplotlib.ticker.LogLocator(base=10,numticks=30) 
    ax1.xaxis.set_major_locator(locmaj)
    ax1.tick_params(axis='x',which='minor',direction='in')
    ax1.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax1.minorticks_on()
    ax1.tick_params(which='major',direction='in',bottom=True, top=True, left=True, right=True, labelsize=20, length=10)
    ax1.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True, length=5)
    lgnd = plt.legend(framealpha=0.5, labelspacing=0.3, fontsize=18, loc="upper left")
    for handle in lgnd.legend_handles:
        handle.set_sizes([100.0])
    
    
    
    ax2 = plt.subplot(4,2,2)
    
    ax2.plot(x_new, y_new_swirl, '--',color='black', linewidth=3, markersize=12,label=(f'Swirl $R^2$: {adjr_s:.2g}'), zorder=10)
    
    if error_plot == '20' or error_plot == '10' or error_plot == '50':
        
        ax2.plot(x_new, y_new_swirl_p, '--',color='darkviolet',linewidth=3, markersize=12, label=f'+ {error}%: {adjr_sp:.2g}')
        ax2.plot(x_new, y_new_swirl_m, '--',color='sandybrown', linewidth=3, markersize=12,label=f'- {error}%: {adjr_sm:.2g}')
        
        # plus error
        ax2.scatter(KH_dusk_medians_plus[0:2]/1e9, medians_swirl[0:2],s=100,color='plum',marker="^")
        ax2.scatter(KH_dusk_medians_plus[2]/1e9, medians_swirl[2],s=100,color='plum')
        ax2.scatter(KH_dusk_medians_plus[3:5]/1e9, medians_swirl[3:5],s=100,color='plum',marker="^")
        ax2.scatter(KH_dusk_medians_plus[5:9]/1e9, medians_swirl[5:9],s=100,color='plum')
        ax2.scatter(KH_dusk_medians_plus[9]/1e9, medians_swirl[9],s=100,color='plum',marker="v")
        ax2.scatter(KH_dusk_medians_plus[10]/1e9, medians_swirl[10],s=100,color='plum')
        ax2.scatter(KH_dusk_medians_plus[11]/1e9, medians_swirl[11],s=100,color='plum',marker="v")
        ax2.scatter(KH_dusk_medians_plus[12:18]/1e9, medians_swirl[12:18],s=100,color='plum',marker=",")
        ax2.scatter(KH_dusk_medians_plus[18:21]/1e9, medians_swirl[18:21],s=100,color='plum',marker="v")
        ax2.scatter(KH_dusk_medians_plus[21]/1e9, medians_swirl[21],s=100,color='plum',marker="^")
        ax2.scatter(KH_dusk_medians_plus[22:]/1e9, medians_swirl[22:],s=100,color='plum',marker=",")
    
        # minus error
        ax2.scatter(KH_dusk_medians_minus[0:2]/1e9, medians_swirl[0:2],s=100,color='peachpuff',marker="^")
        ax2.scatter(KH_dusk_medians_minus[2]/1e9, medians_swirl[2],s=100,color='peachpuff')
        ax2.scatter(KH_dusk_medians_minus[3:5]/1e9, medians_swirl[3:5],s=100,color='peachpuff',marker="^")
        ax2.scatter(KH_dusk_medians_minus[5:9]/1e9, medians_swirl[5:9],s=100,color='peachpuff')
        ax2.scatter(KH_dusk_medians_minus[9]/1e9, medians_swirl[9],s=100,color='peachpuff',marker="v")
        ax2.scatter(KH_dusk_medians_minus[10]/1e9, medians_swirl[10],s=100,color='peachpuff')
        ax2.scatter(KH_dusk_medians_minus[11]/1e9, medians_swirl[11],s=100,color='peachpuff',marker="v")
        ax2.scatter(KH_dusk_medians_minus[12:18]/1e9, medians_swirl[12:18],s=100,color='peachpuff',marker=",")
        ax2.scatter(KH_dusk_medians_minus[18:21]/1e9, medians_swirl[18:21],s=100,color='peachpuff',marker="v")
        ax2.scatter(KH_dusk_medians_minus[21]/1e9, medians_swirl[21],s=100,color='peachpuff',marker="^")
        ax2.scatter(KH_dusk_medians_minus[22:]/1e9, medians_swirl[22:],s=100,color='peachpuff',marker=",")
        
        ax2.errorbar(KH_dusk_medians_plus/1e9, medians_swirl, yerr=swirl_med_errs, xerr=med_KH_dusk_errs_plus/1e9, fmt='.', color='plum')
        
        ax2.errorbar(KH_dusk_medians_minus/1e9, medians_swirl, yerr=swirl_med_errs, xerr=med_KH_dusk_errs_minus/1e9, fmt='.', color='peachpuff')
    
    
    # main plot
    ax2.scatter(Bz_pos_KH_dusk[0:2]/1e9, Bz_pos_swirlpower[0:2], color='red', s=200,zorder=9)
    ax2.scatter(Bz_pos_KH_dusk[2:]/1e9, Bz_pos_swirlpower[2:], color='red', s=200, marker='v',zorder=9)

    ax2.scatter(Bz_neg_KH_dusk[0]/1e9, Bz_neg_swirlpower[0], color='orange', s=200,zorder=9) #
    ax2.scatter(Bz_neg_KH_dusk[1]/1e9, Bz_neg_swirlpower[1], color='orange', s=200,marker=',',zorder=9) 

    ax2.scatter(By_pos_KH_dusk[0]/1e9, By_pos_swirlpower[0], color='green', s=200,zorder=9) # compression 1 2-5 incl visit 15
    ax2.scatter(By_pos_KH_dusk[1:4]/1e9, By_pos_swirlpower[1:4], color='green', s=200, marker=",",zorder=9)
    ax2.scatter(By_pos_KH_dusk[4:6]/1e9, By_pos_swirlpower[4:6], color='green', s=200, marker="v",zorder=9) # rarefraction
    ax2.scatter(By_pos_KH_dusk[6]/1e9, By_pos_swirlpower[6], color='green', s=200, marker="^",zorder=9)
    ax2.scatter(By_pos_KH_dusk[7:]/1e9, By_pos_swirlpower[7:], color='green', s=200, marker=",",zorder=9)

    ax2.scatter(By_neg_KH_dusk[0:4]/1e9, By_neg_swirlpower[0:4], color='blue', s=200, marker="^",zorder=9)
    ax2.scatter(By_neg_KH_dusk[4:6]/1e9, By_neg_swirlpower[4:6], color='blue', s=200,zorder=9) # 
    ax2.scatter(By_neg_KH_dusk[6]/1e9, By_neg_swirlpower[6], color='blue', s=200, marker="v",zorder=9)
    ax2.scatter(By_neg_KH_dusk[7:]/1e9, By_neg_swirlpower[7:], color='blue', s=200, marker=",",zorder=9)
    
    ax2.errorbar(KH_dusk_medians/1e9, medians_swirl, yerr=swirl_med_errs, xerr=med_KH_dusk_errs/1e9, fmt='.', color='lightgray', zorder=9)
    
    ax2.text(5.5e2,6.8,'b',style='italic',fontsize=40) #155 for no plus/minus, 160 otherwise
    
    ax2.set_xlabel('Kelvin Helmholtz Reconnection Power (GW)', fontsize=22)
    ax2.set_ylabel('Median Swirl Region Power (GW)', fontsize=22)
    
    ax2.set_xlim(0.3e1,8e2)
    #ax2.set_xlim(0.00002, 350)
    #ax2.set_ylim(-5,120)
    ax2.set_xscale('log')
    locmaj = matplotlib.ticker.LogLocator(base=10,numticks=30) 
    ax2.xaxis.set_major_locator(locmaj)
    ax2.tick_params(axis='x',which='minor',direction='in')
    ax2.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax2.minorticks_on()
    ax2.tick_params(which='major',direction='in',bottom=True, top=True, left=True, right=True, labelsize=20, length=10)
    ax2.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True, length=5)
    ax2.legend(framealpha=0.5,fontsize=21,loc='upper left')#, title="Swirl Region", title_fontsize=14) #bbox_to_anchor=(-0.67, 0.2)
    
    
    ax3 = plt.subplot(4,2,3)
    
    ax3.plot(x_new, y_new_dusk,'--', color='black',linewidth=3, markersize=12, label=f'Dusk $R^2$: {adjr_d:.2g}', zorder=10)
    
    if error_plot == '20' or error_plot == '10' or error_plot == '50':
    
        ax3.plot(x_new, y_new_dusk_p, '--',color='darkviolet',linewidth=3, markersize=12, label=f'+ {error}%: {adjr_dp:.2g}')
        ax3.plot(x_new, y_new_dusk_m, '--',color='sandybrown',linewidth=3, markersize=12, label=f'- {error}%: {adjr_dm:.2g}')
        
        # plus error
        ax3.scatter(KH_dusk_medians_plus[0:2]/1e9, dusk_medians[0:2],s=100,color='plum',marker="^")
        ax3.scatter(KH_dusk_medians_plus[2]/1e9, dusk_medians[2],s=100,color='plum')
        ax3.scatter(KH_dusk_medians_plus[3:5]/1e9, dusk_medians[3:5],s=100,color='plum',marker="^")
        ax3.scatter(KH_dusk_medians_plus[5:9]/1e9, dusk_medians[5:9],s=100,color='plum')
        ax3.scatter(KH_dusk_medians_plus[9]/1e9, dusk_medians[9],s=100,color='plum',marker="v")
        ax3.scatter(KH_dusk_medians_plus[10]/1e9, dusk_medians[10],s=100,color='plum')
        ax3.scatter(KH_dusk_medians_plus[11]/1e9, dusk_medians[11],s=100,color='plum',marker="v")
        ax3.scatter(KH_dusk_medians_plus[12:18]/1e9, dusk_medians[12:18],s=100,color='plum',marker=",")
        ax3.scatter(KH_dusk_medians_plus[18:21]/1e9, dusk_medians[18:21],s=100,color='plum',marker="v")
        ax3.scatter(KH_dusk_medians_plus[21]/1e9, dusk_medians[21],s=100,color='plum',marker="^")
        ax3.scatter(KH_dusk_medians_plus[22:]/1e9, dusk_medians[22:],s=100,color='plum',marker=",")
    
        # minus error
        ax3.scatter(KH_dusk_medians_minus[0:2]/1e9, dusk_medians[0:2],s=100,color='peachpuff',marker="^")
        ax3.scatter(KH_dusk_medians_minus[2]/1e9, dusk_medians[2],s=100,color='peachpuff')
        ax3.scatter(KH_dusk_medians_minus[3:5]/1e9, dusk_medians[3:5],s=100,color='peachpuff',marker="^")
        ax3.scatter(KH_dusk_medians_minus[5:9]/1e9, dusk_medians[5:9],s=100,color='peachpuff')
        ax3.scatter(KH_dusk_medians_minus[9]/1e9, dusk_medians[9],s=100,color='peachpuff',marker="v")
        ax3.scatter(KH_dusk_medians_minus[10]/1e9, dusk_medians[10],s=100,color='peachpuff')
        ax3.scatter(KH_dusk_medians_minus[11]/1e9, dusk_medians[11],s=100,color='peachpuff',marker="v")
        ax3.scatter(KH_dusk_medians_minus[12:18]/1e9, dusk_medians[12:18],s=100,color='peachpuff',marker=",")
        ax3.scatter(KH_dusk_medians_minus[18:21]/1e9, dusk_medians[18:21],s=100,color='peachpuff',marker="v")
        ax3.scatter(KH_dusk_medians_minus[21]/1e9, dusk_medians[21],s=100,color='peachpuff',marker="^")
        ax3.scatter(KH_dusk_medians_minus[22:]/1e9, dusk_medians[22:],s=100,color='peachpuff',marker=",")
        
        ax3.errorbar(KH_dusk_medians_plus/1e9, dusk_medians, yerr=dusk_med_errs, xerr=med_KH_dusk_errs_plus/1e9, fmt='.', color='plum')
    
        ax3.errorbar(KH_dusk_medians_minus/1e9, dusk_medians, yerr=dusk_med_errs, xerr=med_KH_dusk_errs_minus/1e9, fmt='.', color='peachpuff')
    
    
    # main plot
    ax3.scatter(Bz_pos_KH_dusk[0:2]/1e9, Bz_pos_duskpower[0:2], color='red', s=200,zorder=9)
    ax3.scatter(Bz_pos_KH_dusk[2:]/1e9, Bz_pos_duskpower[2:], color='red', s=200, marker='v',zorder=9)

    ax3.scatter(Bz_neg_KH_dusk[0]/1e9, Bz_neg_duskpower[0], color='orange', s=200,zorder=9) #
    ax3.scatter(Bz_neg_KH_dusk[1]/1e9, Bz_neg_duskpower[1], color='orange', s=200,marker=',',zorder=9)

    ax3.scatter(By_pos_KH_dusk[0]/1e9, By_pos_duskpower[0], color='green', s=200,zorder=9) # compression 1 2-5 incl visit 15
    ax3.scatter(By_pos_KH_dusk[1:4]/1e9, By_pos_duskpower[1:4], color='green', s=200, marker=",",zorder=9)
    ax3.scatter(By_pos_KH_dusk[4:6]/1e9, By_pos_duskpower[4:6], color='green', s=200, marker="v",zorder=9) # rarefraction
    ax3.scatter(By_pos_KH_dusk[6]/1e9, By_pos_duskpower[6], color='green', s=200, marker="^",zorder=9)
    ax3.scatter(By_pos_KH_dusk[7:]/1e9, By_pos_duskpower[7:], color='green', s=200, marker=",",zorder=9)

    ax3.scatter(By_neg_KH_dusk[0:4]/1e9, By_neg_duskpower[0:4], color='blue', s=200, marker="^",zorder=9)
    ax3.scatter(By_neg_KH_dusk[4:6]/1e9, By_neg_duskpower[4:6], color='blue', s=200,zorder=9) # 
    ax3.scatter(By_neg_KH_dusk[6]/1e9, By_neg_duskpower[6], color='blue', s=200, marker="v",zorder=9)
    ax3.scatter(By_neg_KH_dusk[7:]/1e9, By_neg_duskpower[7:], color='blue', s=200, marker=",",zorder=9)
    
    ax3.errorbar(KH_dusk_medians/1e9, dusk_medians, yerr=dusk_med_errs, xerr=med_KH_dusk_errs/1e9, fmt='.', color='lightgray', zorder=9)
    
    ax3.text(5.5e2,-13,'c',style='italic',fontsize=40) #155 for no plus/minus, 160 otherwise
    
    ax3.set_xlabel('Kelvin Helmholtz Reconnection Power (GW)', fontsize=22)
    ax3.set_ylabel('Median Dusk Region Power (GW)', fontsize=22)
    
    ax3.set_xlim(0.3e1,8e2)
    #ax3.set_xlim(0.00002, 350)
    #ax3.set_ylim(-15,500)
    ax3.set_xscale('log')
    locmaj = matplotlib.ticker.LogLocator(base=10,numticks=30) 
    ax3.xaxis.set_major_locator(locmaj)
    ax3.tick_params(axis='x',which='minor',direction='in')
    ax3.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax3.minorticks_on()
    ax3.tick_params(which='major',direction='in',bottom=True, top=True, left=True, right=True, labelsize=20, length=10)
    ax3.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True, length=5)
    ax3.legend(framealpha=0.5,fontsize=21,loc='upper left'), #title="Dusk Active Region", title_fontsize=14)#legend(framealpha=0.5,fontsize=14,loc ="upper left") box_to_anchor=(1.05, 0.8),
    
    
    ax4 = plt.subplot(4,2,4)
    
    ax4.plot(x_new, y_new_noon, '--',color='black', linewidth=3, markersize=12,label=f'Noon $R^2$: {adjr_n:.2g}', zorder=10)
    
    if error_plot == '20' or error_plot == '10' or error_plot == '50':
    
        ax4.plot(x_new, y_new_noon_p, '--',color='darkviolet',linewidth=3, markersize=12, label=f'+ {error}%: {adjr_np:.2g}')
        ax4.plot(x_new, y_new_noon_m, '--',color='sandybrown',linewidth=3, markersize=12, label=f'- {error}%: {adjr_nm:.2g}')
        
        # plus error
        ax4.scatter(KH_dusk_medians_plus[0:2]/1e9, noon_medians[0:2],s=100,color='plum',marker="^")
        ax4.scatter(KH_dusk_medians_plus[2]/1e9, noon_medians[2],s=100,color='plum')
        ax4.scatter(KH_dusk_medians_plus[3:5]/1e9, noon_medians[3:5],s=100,color='plum',marker="^")
        ax4.scatter(KH_dusk_medians_plus[5:9]/1e9, noon_medians[5:9],s=100,color='plum')
        ax4.scatter(KH_dusk_medians_plus[9]/1e9, noon_medians[9],s=100,color='plum',marker="v")
        ax4.scatter(KH_dusk_medians_plus[10]/1e9, noon_medians[10],s=100,color='plum')
        ax4.scatter(KH_dusk_medians_plus[11]/1e9, noon_medians[11],s=100,color='plum',marker="v")
        ax4.scatter(KH_dusk_medians_plus[12:18]/1e9, noon_medians[12:18],s=100,color='plum',marker=",")
        ax4.scatter(KH_dusk_medians_plus[18:21]/1e9, noon_medians[18:21],s=100,color='plum',marker="v")
        ax4.scatter(KH_dusk_medians_plus[21]/1e9, noon_medians[21],s=100,color='plum',marker="^")
        ax4.scatter(KH_dusk_medians_plus[22:]/1e9, noon_medians[22:],s=100,color='plum',marker=",")
    
        # minus error
        ax4.scatter(KH_dusk_medians_minus[0:2]/1e9, noon_medians[0:2],s=100,color='peachpuff',marker="^")
        ax4.scatter(KH_dusk_medians_minus[2]/1e9, noon_medians[2],s=100,color='peachpuff')
        ax4.scatter(KH_dusk_medians_minus[3:5]/1e9, noon_medians[3:5],s=100,color='peachpuff',marker="^")
        ax4.scatter(KH_dusk_medians_minus[5:9]/1e9, noon_medians[5:9],s=100,color='peachpuff')
        ax4.scatter(KH_dusk_medians_minus[9]/1e9, noon_medians[9],s=100,color='peachpuff',marker="v")
        ax4.scatter(KH_dusk_medians_minus[10]/1e9, noon_medians[10],s=100,color='peachpuff')
        ax4.scatter(KH_dusk_medians_minus[11]/1e9, noon_medians[11],s=100,color='peachpuff',marker="v")
        ax4.scatter(KH_dusk_medians_minus[12:18]/1e9, noon_medians[12:18],s=100,color='peachpuff',marker=",")
        ax4.scatter(KH_dusk_medians_minus[18:21]/1e9, noon_medians[18:21],s=100,color='peachpuff',marker="v")
        ax4.scatter(KH_dusk_medians_minus[21]/1e9, noon_medians[21],s=100,color='peachpuff',marker="^")
        ax4.scatter(KH_dusk_medians_minus[22:]/1e9, noon_medians[22:],s=100,color='peachpuff',marker=",")
        
        ax4.errorbar(KH_dusk_medians_plus/1e9, noon_medians, yerr=noon_med_errs, xerr=med_KH_dusk_errs_plus/1e9, fmt='.', color='plum')
    
        ax4.errorbar(KH_dusk_medians_minus/1e9, noon_medians, yerr=noon_med_errs, xerr=med_KH_dusk_errs_minus/1e9, fmt='.', color='peachpuff')
    
    # main plot
    ax4.scatter(Bz_pos_KH_dusk[0:2]/1e9, Bz_pos_noonpower[0:2], color='red', s=200,zorder=9)
    ax4.scatter(Bz_pos_KH_dusk[2:]/1e9, Bz_pos_noonpower[2:], color='red', s=200, marker='v',zorder=9)

    ax4.scatter(Bz_neg_KH_dusk[0]/1e9, Bz_neg_noonpower[0], color='orange', s=200,zorder=9) #
    ax4.scatter(Bz_neg_KH_dusk[1]/1e9, Bz_neg_noonpower[1], color='orange', s=200,marker=',',zorder=9) 

    ax4.scatter(By_pos_KH_dusk[0]/1e9, By_pos_noonpower[0], color='green', s=200,zorder=9) # compression 1 2-5 incl visit 15
    ax4.scatter(By_pos_KH_dusk[1:4]/1e9, By_pos_noonpower[1:4], color='green', s=200, marker=",",zorder=9)
    ax4.scatter(By_pos_KH_dusk[4:6]/1e9, By_pos_noonpower[4:6], color='green', s=200, marker="v",zorder=9) # rarefraction
    ax4.scatter(By_pos_KH_dusk[6]/1e9, By_pos_noonpower[6], color='green', s=200, marker="^",zorder=9)
    ax4.scatter(By_pos_KH_dusk[7:]/1e9, By_pos_noonpower[7:], color='green', s=200, marker=",",zorder=9)

    ax4.scatter(By_neg_KH_dusk[0:4]/1e9, By_neg_noonpower[0:4], color='blue', s=200, marker="^",zorder=9)
    ax4.scatter(By_neg_KH_dusk[4:6]/1e9, By_neg_noonpower[4:6], color='blue', s=200,zorder=9)
    ax4.scatter(By_neg_KH_dusk[6]/1e9, By_neg_noonpower[6], color='blue', s=200, marker="v",zorder=9)
    ax4.scatter(By_neg_KH_dusk[7:]/1e9, By_neg_noonpower[7:], color='blue', s=200, marker=",",zorder=9)
    
    ax4.errorbar(KH_dusk_medians/1e9, noon_medians, yerr=noon_med_errs, xerr=med_KH_dusk_errs/1e9, fmt='.', color='lightgray', zorder=9)
    
    ax4.text(5.5e2,-12,'d',style='italic',fontsize=40) #155 for no plus/minus, 160 otherwise
    
    ax4.set_xlabel('Kelvin Helmholtz Reconnection Power (GW)', fontsize=22)
    ax4.set_ylabel('Median Noon Region Power (GW)', fontsize=22)
   
    #ax4.set_ylim(1e1,1e3)
    ax4.set_xlim(0.3e1,8e2)
    ax4.set_xscale('log')
    locmaj = matplotlib.ticker.LogLocator(base=10,numticks=30) 
    ax4.xaxis.set_major_locator(locmaj)
    ax4.tick_params(axis='x',which='minor',direction='in')
    ax4.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax4.minorticks_on()
    ax4.tick_params(which='major',direction='in',bottom=True, top=True, left=True, right=True, labelsize=20, length=10)
    ax4.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True, length=5)
    ax4.legend(framealpha=0.5,fontsize=21,loc='upper left')#, title="Noon Active Region", title_fontsize=14)#legend(framealpha=0.5,fontsize=14,loc ="upper left") #bbox_to_anchor=(-0.685, 0.3)
    
    
    # save plot
    saveloc = (f'{root_saves}median_KH_V_dusk_vs_region_power_{error_plot}.jpg') 
    #saveloc = (f'/Users/hannah/OneDrive - Lancaster University/aurora/median_KH_dusk_V_vs_region_power_{error}.jpg')
    plt.savefig(saveloc,bbox_inches='tight',dpi=400)




# space incase python is being a scroll dic
