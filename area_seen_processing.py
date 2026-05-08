#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 11:42:13 2025

@author: hannah
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

root_folder = '/Users/hannah/OneDrive - Lancaster University/aurora/python_scripts/dataframes/'

# high CML visits
visit_01 = pd.read_csv(root_folder+'new_total_powers_regions_01_CR_adj.csv',delimiter=',')
visit_12 = pd.read_csv(root_folder+'new_total_powers_regions_12_CR_adj.csv',delimiter=',')
visit_15 = pd.read_csv(root_folder+'new_total_powers_regions_15_CR_adj.csv',delimiter=',')
visit_26 = pd.read_csv(root_folder+'new_total_powers_regions_26_CR_adj.csv',delimiter=',')

#visit_11 = pd.read_csv(root_folder+'new_total_powers_regions_11_march.csv',delimiter=',')

cml_01 = visit_01['CML'].to_numpy()
cml_12 = visit_12['CML'].to_numpy()
cml_15 = visit_15['CML'].to_numpy()
cml_26 = visit_26['CML'].to_numpy()

#cml_11 = visit_11['CML'].to_numpy()


cmls_high = [*cml_01, *cml_12, *cml_15]#, *cml_26]


swirl_tot_01 = visit_01['Total_Region_Swirl'].to_numpy()
swirl_px_01 = visit_01['Total_Pixels_Swirl'].to_numpy()
swirl_tot_12 = visit_12['Total_Region_Swirl'].to_numpy()
swirl_px_12 = visit_12['Total_Pixels_Swirl'].to_numpy()
swirl_tot_15 = visit_15['Total_Region_Swirl'].to_numpy()
swirl_px_15 = visit_15['Total_Pixels_Swirl'].to_numpy()
swirl_tot_26 = visit_26['Total_Region_Swirl'].to_numpy()
swirl_px_26 = visit_26['Total_Pixels_Swirl'].to_numpy()

# swirl_tot_11 = visit_11['Total_Region_Swirl'].to_numpy()
# swirl_px_11 = visit_11['Total_Pixels_Swirl'].to_numpy()

swirl_tot_high = [*swirl_tot_01, *swirl_tot_12, *swirl_tot_15]#, *swirl_tot_26]
swirl_px_high = [*swirl_px_01, *swirl_px_12, *swirl_px_15]#, *swirl_px_26]


dusk_tot_01 = visit_01['Total_Region_Dusk'].to_numpy()
dusk_px_01 = visit_01['Total_Pixels_Dusk'].to_numpy()
dusk_tot_12 = visit_12['Total_Region_Dusk'].to_numpy()
dusk_px_12 = visit_12['Total_Pixels_Dusk'].to_numpy()
dusk_tot_15 = visit_15['Total_Region_Dusk'].to_numpy()
dusk_px_15 = visit_15['Total_Pixels_Dusk'].to_numpy()
dusk_tot_26 = visit_26['Total_Region_Dusk'].to_numpy()
dusk_px_26 = visit_26['Total_Pixels_Dusk'].to_numpy()

# dusk_tot_11 = visit_11['Total_Region_Dusk'].to_numpy()
# dusk_px_11 = visit_11['Total_Pixels_Dusk'].to_numpy()

dusk_tot_high = [*dusk_tot_01, *dusk_tot_12, *dusk_tot_15]#, *dusk_tot_26]
dusk_px_high = [*dusk_px_01, *dusk_px_12, *dusk_px_15]#, *dusk_px_26]


noon_tot_01 = visit_01['Total_Region_Noon'].to_numpy()
noon_px_01 = visit_01['Total_Pixels_Noon'].to_numpy()
noon_tot_12 = visit_12['Total_Region_Noon'].to_numpy()
noon_px_12 = visit_12['Total_Pixels_Noon'].to_numpy()
noon_tot_15 = visit_15['Total_Region_Noon'].to_numpy()
noon_px_15 = visit_15['Total_Pixels_Noon'].to_numpy()
noon_tot_26 = visit_26['Total_Region_Noon'].to_numpy()
noon_px_26 = visit_26['Total_Pixels_Noon'].to_numpy()

# noon_tot_11 = visit_11['Total_Region_Noon'].to_numpy()
# noon_px_11 = visit_11['Total_Pixels_Noon'].to_numpy()

noon_tot_high = [*noon_tot_01, *noon_tot_12, *noon_tot_15]#, *noon_tot_26]
noon_px_high = [*noon_px_01, *noon_px_12, *noon_px_15]#, *noon_px_26]


# percentage calculation for pixels seen
percent_swirl_high = [int(d) / int(x) for d,x in zip(swirl_px_high,swirl_tot_high)]
percentage_swirl_high = [s * 100 for s in percent_swirl_high]
percent_noon_high = [int(e) / int (y) for e,y in zip(noon_px_high,noon_tot_high)]
percentage_noon_high = [t * 100 for t in percent_noon_high]
percent_dusk_high = [int(f) / int(z) for f,z in zip(dusk_px_high,dusk_tot_high)]
percentage_dusk_high = [u * 100 for u in percent_dusk_high]

# make new dataframe
high_df = pd.DataFrame()
high_df = high_df.assign(CML=cmls_high)
high_df = high_df.assign(Swirl_Percentage=percentage_swirl_high)
high_df = high_df.assign(Dusk_Percentage=percentage_dusk_high)
high_df = high_df.assign(Noon_Percentage=percentage_noon_high)

high_df.sort_values(by='CML', inplace=True)
high_df.to_csv(root_folder+'percentage_seen_high_CR_adj.csv',index=False) #_march


# standard CML visits ['02','03','04','05','08','09','10','11','16','17','18','19','20','21','24','25','27','28','34','35']
visit_02 = pd.read_csv(root_folder+'new_total_powers_regions_02_CR_adj.csv',delimiter=',')
visit_03 = pd.read_csv(root_folder+'new_total_powers_regions_03_CR_adj.csv',delimiter=',')
visit_04 = pd.read_csv(root_folder+'new_total_powers_regions_04_CR_adj.csv',delimiter=',')
visit_05 = pd.read_csv(root_folder+'new_total_powers_regions_05_CR_adj.csv',delimiter=',')
visit_08 = pd.read_csv(root_folder+'new_total_powers_regions_08_CR_adj.csv',delimiter=',')
visit_09 = pd.read_csv(root_folder+'new_total_powers_regions_09_CR_adj.csv',delimiter=',')
visit_10 = pd.read_csv(root_folder+'new_total_powers_regions_10_CR_adj.csv',delimiter=',')
visit_11 = pd.read_csv(root_folder+'new_total_powers_regions_11_CR_adj.csv',delimiter=',') #march2
visit_16 = pd.read_csv(root_folder+'new_total_powers_regions_16_CR_adj.csv',delimiter=',')
visit_17 = pd.read_csv(root_folder+'new_total_powers_regions_17_CR_adj.csv',delimiter=',')
visit_18 = pd.read_csv(root_folder+'new_total_powers_regions_18_CR_adj.csv',delimiter=',')
visit_19 = pd.read_csv(root_folder+'new_total_powers_regions_19_CR_adj.csv',delimiter=',')
visit_20 = pd.read_csv(root_folder+'new_total_powers_regions_20_CR_adj.csv',delimiter=',')
visit_21 = pd.read_csv(root_folder+'new_total_powers_regions_21_CR_adj.csv',delimiter=',')
visit_24 = pd.read_csv(root_folder+'new_total_powers_regions_24_CR_adj.csv',delimiter=',')
visit_25 = pd.read_csv(root_folder+'new_total_powers_regions_25_CR_adj.csv',delimiter=',')
visit_27 = pd.read_csv(root_folder+'new_total_powers_regions_27_CR_adj.csv',delimiter=',')
visit_28 = pd.read_csv(root_folder+'new_total_powers_regions_28_CR_adj.csv',delimiter=',')
visit_34 = pd.read_csv(root_folder+'new_total_powers_regions_34_CR_adj.csv',delimiter=',')
visit_35 = pd.read_csv(root_folder+'new_total_powers_regions_35_CR_adj.csv',delimiter=',')

cml_02 = visit_02['CML'].to_numpy()
cml_03 = visit_03['CML'].to_numpy()
cml_04 = visit_04['CML'].to_numpy()
cml_05 = visit_05['CML'].to_numpy()
cml_08 = visit_08['CML'].to_numpy()
cml_09 = visit_09['CML'].to_numpy()
cml_10 = visit_10['CML'].to_numpy()
cml_11 = visit_11['CML'].to_numpy()
cml_16 = visit_16['CML'].to_numpy()
cml_17 = visit_17['CML'].to_numpy()
cml_18 = visit_18['CML'].to_numpy()
cml_19 = visit_19['CML'].to_numpy()
cml_20 = visit_20['CML'].to_numpy()
cml_21 = visit_21['CML'].to_numpy()
cml_24 = visit_24['CML'].to_numpy()
cml_25 = visit_25['CML'].to_numpy()
cml_27 = visit_27['CML'].to_numpy()
cml_28 = visit_28['CML'].to_numpy()
cml_34 = visit_34['CML'].to_numpy()
cml_35 = visit_35['CML'].to_numpy()

cmls_standard = [*cml_02, *cml_03, *cml_04, *cml_05, *cml_08, *cml_09, *cml_10, *cml_11, *cml_16, *cml_17, *cml_18, *cml_19, *cml_20, *cml_21, *cml_24, *cml_25, *cml_27, *cml_28,*cml_34, *cml_35]# *cml_11,


swirl_tot_02 = visit_02['Total_Region_Swirl'].to_numpy()
swirl_px_02 = visit_02['Total_Pixels_Swirl'].to_numpy()
swirl_tot_03 = visit_03['Total_Region_Swirl'].to_numpy()
swirl_px_03 = visit_03['Total_Pixels_Swirl'].to_numpy()
swirl_tot_04 = visit_04['Total_Region_Swirl'].to_numpy()
swirl_px_04 = visit_04['Total_Pixels_Swirl'].to_numpy()
swirl_tot_05 = visit_05['Total_Region_Swirl'].to_numpy()
swirl_px_05 = visit_05['Total_Pixels_Swirl'].to_numpy()
swirl_tot_08 = visit_08['Total_Region_Swirl'].to_numpy()
swirl_px_08 = visit_08['Total_Pixels_Swirl'].to_numpy()
swirl_tot_09 = visit_09['Total_Region_Swirl'].to_numpy()
swirl_px_09 = visit_09['Total_Pixels_Swirl'].to_numpy()
swirl_tot_10 = visit_10['Total_Region_Swirl'].to_numpy()
swirl_px_10 = visit_10['Total_Pixels_Swirl'].to_numpy()
swirl_tot_11 = visit_11['Total_Region_Swirl'].to_numpy()
swirl_px_11 = visit_11['Total_Pixels_Swirl'].to_numpy()
swirl_tot_16 = visit_16['Total_Region_Swirl'].to_numpy()
swirl_px_16 = visit_16['Total_Pixels_Swirl'].to_numpy()
swirl_tot_17 = visit_17['Total_Region_Swirl'].to_numpy()
swirl_px_17 = visit_17['Total_Pixels_Swirl'].to_numpy()
swirl_tot_18 = visit_18['Total_Region_Swirl'].to_numpy()
swirl_px_18 = visit_18['Total_Pixels_Swirl'].to_numpy()
swirl_tot_19 = visit_19['Total_Region_Swirl'].to_numpy()
swirl_px_19 = visit_19['Total_Pixels_Swirl'].to_numpy()
swirl_tot_20 = visit_20['Total_Region_Swirl'].to_numpy()
swirl_px_20 = visit_20['Total_Pixels_Swirl'].to_numpy()
swirl_tot_21 = visit_21['Total_Region_Swirl'].to_numpy()
swirl_px_21 = visit_21['Total_Pixels_Swirl'].to_numpy()
swirl_tot_24 = visit_24['Total_Region_Swirl'].to_numpy()
swirl_px_24 = visit_24['Total_Pixels_Swirl'].to_numpy()
swirl_tot_25 = visit_25['Total_Region_Swirl'].to_numpy()
swirl_px_25 = visit_25['Total_Pixels_Swirl'].to_numpy()
swirl_tot_27 = visit_27['Total_Region_Swirl'].to_numpy()
swirl_px_27 = visit_27['Total_Pixels_Swirl'].to_numpy()
swirl_tot_28 = visit_28['Total_Region_Swirl'].to_numpy()
swirl_px_28 = visit_28['Total_Pixels_Swirl'].to_numpy()
swirl_tot_34 = visit_34['Total_Region_Swirl'].to_numpy()
swirl_px_34 = visit_34['Total_Pixels_Swirl'].to_numpy()
swirl_tot_35 = visit_35['Total_Region_Swirl'].to_numpy()
swirl_px_35 = visit_35['Total_Pixels_Swirl'].to_numpy()

swirl_tot_standard = [*swirl_tot_02, *swirl_tot_03, *swirl_tot_04, *swirl_tot_05, *swirl_tot_08, *swirl_tot_09, *swirl_tot_10, *swirl_tot_11, *swirl_tot_16, *swirl_tot_17, *swirl_tot_18, *swirl_tot_19, *swirl_tot_20, *swirl_tot_21, *swirl_tot_24, *swirl_tot_25, *swirl_tot_27,*swirl_tot_28, *swirl_tot_34, *swirl_tot_35]#*swirl_tot_11,
swirl_px_standard = [*swirl_px_02, *swirl_px_03, *swirl_px_04, *swirl_px_05, *swirl_px_08, *swirl_px_09, *swirl_px_10, *swirl_px_11,*swirl_px_16, *swirl_px_17, *swirl_px_18, *swirl_px_19, *swirl_px_20, *swirl_px_21, *swirl_px_24, *swirl_px_25, *swirl_px_27, *swirl_px_28,*swirl_px_34, *swirl_px_35]#*swirl_px_11,



dusk_tot_02 = visit_02['Total_Region_Dusk'].to_numpy()
dusk_px_02 = visit_02['Total_Pixels_Dusk'].to_numpy()
dusk_tot_03 = visit_03['Total_Region_Dusk'].to_numpy()
dusk_px_03 = visit_03['Total_Pixels_Dusk'].to_numpy()
dusk_tot_04 = visit_04['Total_Region_Dusk'].to_numpy()
dusk_px_04 = visit_04['Total_Pixels_Dusk'].to_numpy()
dusk_tot_05 = visit_05['Total_Region_Dusk'].to_numpy()
dusk_px_05 = visit_05['Total_Pixels_Dusk'].to_numpy()
dusk_tot_08 = visit_08['Total_Region_Dusk'].to_numpy()
dusk_px_08 = visit_08['Total_Pixels_Dusk'].to_numpy()
dusk_tot_09 = visit_09['Total_Region_Dusk'].to_numpy()
dusk_px_09 = visit_09['Total_Pixels_Dusk'].to_numpy()
dusk_tot_10 = visit_10['Total_Region_Dusk'].to_numpy()
dusk_px_10 = visit_10['Total_Pixels_Dusk'].to_numpy()
dusk_tot_11 = visit_11['Total_Region_Dusk'].to_numpy()
dusk_px_11 = visit_11['Total_Pixels_Dusk'].to_numpy()
dusk_tot_16 = visit_16['Total_Region_Dusk'].to_numpy()
dusk_px_16 = visit_16['Total_Pixels_Dusk'].to_numpy()
dusk_tot_17 = visit_17['Total_Region_Dusk'].to_numpy()
dusk_px_17 = visit_17['Total_Pixels_Dusk'].to_numpy()
dusk_tot_18 = visit_18['Total_Region_Dusk'].to_numpy()
dusk_px_18 = visit_18['Total_Pixels_Dusk'].to_numpy()
dusk_tot_19 = visit_19['Total_Region_Dusk'].to_numpy()
dusk_px_19 = visit_19['Total_Pixels_Dusk'].to_numpy()
dusk_tot_20 = visit_20['Total_Region_Dusk'].to_numpy()
dusk_px_20 = visit_20['Total_Pixels_Dusk'].to_numpy()
dusk_tot_21 = visit_21['Total_Region_Dusk'].to_numpy()
dusk_px_21 = visit_21['Total_Pixels_Dusk'].to_numpy()
dusk_tot_24 = visit_24['Total_Region_Dusk'].to_numpy()
dusk_px_24 = visit_24['Total_Pixels_Dusk'].to_numpy()
dusk_tot_25 = visit_25['Total_Region_Dusk'].to_numpy()
dusk_px_25 = visit_25['Total_Pixels_Dusk'].to_numpy()
dusk_tot_27 = visit_27['Total_Region_Dusk'].to_numpy()
dusk_px_27 = visit_27['Total_Pixels_Dusk'].to_numpy()
dusk_tot_28 = visit_28['Total_Region_Dusk'].to_numpy()
dusk_px_28 = visit_28['Total_Pixels_Dusk'].to_numpy()
dusk_tot_34 = visit_34['Total_Region_Dusk'].to_numpy()
dusk_px_34 = visit_34['Total_Pixels_Dusk'].to_numpy()
dusk_tot_35 = visit_35['Total_Region_Dusk'].to_numpy()
dusk_px_35 = visit_35['Total_Pixels_Dusk'].to_numpy()

dusk_tot_standard = [*dusk_tot_02, *dusk_tot_03, *dusk_tot_04, *dusk_tot_05, *dusk_tot_08, *dusk_tot_09, *dusk_tot_10, *dusk_tot_11, *dusk_tot_16, *dusk_tot_17, *dusk_tot_18, *dusk_tot_19, *dusk_tot_20, *dusk_tot_21, *dusk_tot_24, *dusk_tot_25, *dusk_tot_27, *dusk_tot_28,*dusk_tot_34, *dusk_tot_35]#*dusk_tot_11,
dusk_px_standard = [*dusk_px_02, *dusk_px_03, *dusk_px_04, *dusk_px_05, *dusk_px_08, *dusk_px_09, *dusk_px_10, *dusk_px_11, *dusk_px_16, *dusk_px_17, *dusk_px_18, *dusk_px_19, *dusk_px_20, *dusk_px_21, *dusk_px_24, *dusk_px_25, *dusk_px_27,*dusk_px_28, *dusk_px_34, *dusk_px_35]#*dusk_px_11,


noon_tot_02 = visit_02['Total_Region_Noon'].to_numpy()
noon_px_02 = visit_02['Total_Pixels_Noon'].to_numpy()
noon_tot_03 = visit_03['Total_Region_Noon'].to_numpy()
noon_px_03 = visit_03['Total_Pixels_Noon'].to_numpy()
noon_tot_04 = visit_04['Total_Region_Noon'].to_numpy()
noon_px_04 = visit_04['Total_Pixels_Noon'].to_numpy()
noon_tot_05 = visit_05['Total_Region_Noon'].to_numpy()
noon_px_05 = visit_05['Total_Pixels_Noon'].to_numpy()
noon_tot_08 = visit_08['Total_Region_Noon'].to_numpy()
noon_px_08 = visit_08['Total_Pixels_Noon'].to_numpy()
noon_tot_09 = visit_09['Total_Region_Noon'].to_numpy()
noon_px_09 = visit_09['Total_Pixels_Noon'].to_numpy()
noon_tot_10 = visit_10['Total_Region_Noon'].to_numpy()
noon_px_10 = visit_10['Total_Pixels_Noon'].to_numpy()
noon_tot_11 = visit_11['Total_Region_Noon'].to_numpy()
noon_px_11 = visit_11['Total_Pixels_Noon'].to_numpy()
noon_tot_16 = visit_16['Total_Region_Noon'].to_numpy()
noon_px_16 = visit_16['Total_Pixels_Noon'].to_numpy()
noon_tot_17 = visit_17['Total_Region_Noon'].to_numpy()
noon_px_17 = visit_17['Total_Pixels_Noon'].to_numpy()
noon_tot_18 = visit_18['Total_Region_Noon'].to_numpy()
noon_px_18 = visit_18['Total_Pixels_Noon'].to_numpy()
noon_tot_19 = visit_19['Total_Region_Noon'].to_numpy()
noon_px_19 = visit_19['Total_Pixels_Noon'].to_numpy()
noon_tot_20 = visit_20['Total_Region_Noon'].to_numpy()
noon_px_20 = visit_20['Total_Pixels_Noon'].to_numpy()
noon_tot_21 = visit_21['Total_Region_Noon'].to_numpy()
noon_px_21 = visit_21['Total_Pixels_Noon'].to_numpy()
noon_tot_24 = visit_24['Total_Region_Noon'].to_numpy()
noon_px_24 = visit_24['Total_Pixels_Noon'].to_numpy()
noon_tot_25 = visit_25['Total_Region_Noon'].to_numpy()
noon_px_25 = visit_25['Total_Pixels_Noon'].to_numpy()
noon_tot_27 = visit_27['Total_Region_Noon'].to_numpy()
noon_px_27 = visit_27['Total_Pixels_Noon'].to_numpy()
noon_tot_28 = visit_28['Total_Region_Noon'].to_numpy()
noon_px_28 = visit_28['Total_Pixels_Noon'].to_numpy()
noon_tot_34 = visit_34['Total_Region_Noon'].to_numpy()
noon_px_34 = visit_34['Total_Pixels_Noon'].to_numpy()
noon_tot_35 = visit_35['Total_Region_Noon'].to_numpy()
noon_px_35 = visit_35['Total_Pixels_Noon'].to_numpy()

noon_tot_standard = [*noon_tot_02, *noon_tot_03, *noon_tot_04, *noon_tot_05, *noon_tot_08, *noon_tot_09, *noon_tot_10, *noon_tot_11, *noon_tot_16, *noon_tot_17, *noon_tot_18, *noon_tot_19, *noon_tot_20, *noon_tot_21, *noon_tot_24, *noon_tot_25, *noon_tot_27,*noon_tot_28, *noon_tot_34, *noon_tot_35]# *noon_tot_11,
noon_px_standard = [*noon_px_02, *noon_px_03, *noon_px_04, *noon_px_05, *noon_px_08, *noon_px_09, *noon_px_10, *noon_px_11, *noon_px_16, *noon_px_17, *noon_px_18, *noon_px_19, *noon_px_20, *noon_px_21, *noon_px_24, *noon_px_25, *noon_px_27, *noon_px_28,*noon_px_34, *noon_px_35]# *noon_px_11,



# percentage calculation for pixels seen
percent_swirl_standard = [int(b) / int(m) for b,m in zip(swirl_px_standard,swirl_tot_standard)]
percentage_swirl_standard = [p * 100 for p in percent_swirl_standard]
percent_noon_standard = [int(a) / int (n) for a,n in zip(noon_px_standard,noon_tot_standard)]
percentage_noon_standard = [q * 100 for q in percent_noon_standard]
percent_dusk_standard = [int(c) / int(o) for c,o in zip(dusk_px_standard,dusk_tot_standard)]
percentage_dusk_standard = [r * 100 for r in percent_dusk_standard]

# make new dataframe
standard_df = pd.DataFrame()
standard_df = standard_df.assign(CML=cmls_standard)
standard_df = standard_df.assign(Swirl_Percentage=percentage_swirl_standard)
standard_df = standard_df.assign(Dusk_Percentage=percentage_dusk_standard)
standard_df = standard_df.assign(Noon_Percentage=percentage_noon_standard)
standard_df.sort_values(by='CML', inplace=True)

standard_df.to_csv(root_folder+'percentage_seen_standard_CR_adj.csv',index=False)
