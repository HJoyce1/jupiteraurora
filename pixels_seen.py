#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 11:42:13 2025

@author: hannah
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy.polynomial.polynomial as poly

regions = 'high'

root_folder = '/Users/hannah/OneDrive - Lancaster University/aurora/python_scripts/dataframes/'

df_high = pd.read_csv(root_folder+'percentage_seen_high_CR_adj.csv',delimiter=',')
df_standard = pd.read_csv(root_folder+'percentage_seen_standard_CR_adj.csv',delimiter=',')

cml_high = df_high['CML'].to_numpy()
cml_standard = df_standard['CML'].to_numpy()

swirl_percentage_high = df_high['Swirl_Percentage'].to_numpy()
noon_percentage_high = df_high['Noon_Percentage'].to_numpy()
dusk_percentage_high = df_high['Dusk_Percentage'].to_numpy()

swirl_percentage_standard = df_standard['Swirl_Percentage'].to_numpy()
noon_percentage_standard = df_standard['Noon_Percentage'].to_numpy()
dusk_percentage_standard = df_standard['Dusk_Percentage'].to_numpy()


if regions == 'high':
    x_new = np.linspace(min(cml_high), max(cml_high), len(cml_high))
    cml = cml_high
    dusk_percentage = dusk_percentage_high
    swirl_percentage = swirl_percentage_high
    noon_percentage = noon_percentage_high
    
    new_dusk = np.where((dusk_percentage < 100) & (cml < 205), 100, dusk_percentage)
    
    a = np.polyfit(cml_high, new_dusk, 13)
    aa = np.poly1d(a)

    b = np.polyfit(cml_high, swirl_percentage_high, 5)
    bb = np.poly1d(b)

    c = np.polyfit(cml_high, noon_percentage_high, 1)
    cc = np.poly1d(c)
    
    y_new_dusk = aa(x_new)
    y_new_noon = cc(x_new)
    y_new_swirl = bb(x_new)
    
else:
    x_new = np.linspace(min(cml_standard), max(cml_standard), len(cml_standard))
    cml = cml_standard
    dusk_percentage = dusk_percentage_standard
    swirl_percentage = swirl_percentage_standard
    noon_percentage = noon_percentage_standard

    a = np.polyfit(cml_standard, dusk_percentage_standard, 5)
    aa = np.poly1d(a)

    b = np.polyfit(cml_standard, swirl_percentage_standard, 5)
    bb = np.poly1d(b)

    c = np.polyfit(cml_standard, noon_percentage_standard, 1)
    cc = np.poly1d(c)
    
    y_new_dusk = aa(x_new)
    y_new_noon = cc(x_new)
    y_new_swirl = bb(x_new)
    

fig = plt.figure(figsize=(15,8))
plt.rcParams['font.size'] = '14'
ax = plt.subplot(1,1,1)

plt.scatter(cml,swirl_percentage,color='orange',s=25,label='Swirl Region')
plt.plot(x_new, y_new_swirl, color='red',linewidth=2.8)

plt.scatter(cml,noon_percentage,color='blue',s=10,label='Noon Active Region')
plt.plot(x_new, y_new_noon, color='darkblue')

plt.scatter(cml,dusk_percentage,color='limegreen',s=10,label='Dusk Active Region')
plt.plot(x_new, y_new_dusk, color='darkgreen')

lgnd = plt.legend(framealpha=0.5, labelspacing=0.5, fontsize=14, loc="lower left")
#ax.legend(framealpha=0.5)
for handle in lgnd.legend_handles:
    handle.set_sizes([100.0])
ax.set_xlabel('CML (\N{DEGREE SIGN})')
ax.set_ylabel('Percentage Seen (%)')
ax.tick_params(direction='in',bottom=True, top=True, left=True, right=True)

saveloc = (f'/Users/hannah/OneDrive - Lancaster University/aurora/percentage_{regions}_region_CR_adj.jpg')
plt.savefig(saveloc,bbox_inches='tight',dpi=400)

# data = {'Dusk_Active_Fit': a, 'Noon_Active_Fit_Lower': d, 'Noon_Active_Fit_Upper': e, 'Swirl_Fit_Bottom': f, 'Swirl_Fit_Lower': b, 'Swirl_Fit_Upper': c}
data = {'Dusk_Active_Fit': a, 'Noon_Active_Fit': c, 'Swirl_Fit': b}, #'Swirl_Fit_Upper': c} 

dictionary = dict(Dusk_Active_Fit=a, Noon_Active_Fit=c, Swirl_Fit=b)#, Swirl_Fit_Upper=c)
polynomials_df = pd.DataFrame({k : pd.Series(v) for k, v in dictionary.items()})

polynomials_df.to_csv(f'{root_folder}polynomials_{regions}_CR_adj.csv',index=False)