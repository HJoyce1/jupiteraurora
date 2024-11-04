#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 14:34:01 2024

@author: hannah

file to plot out cml vs percentage seen - also saves out dataframe of coefficients for scaling power values
so that powers can be scaled up appropriately to account for HST viewing angle
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

root_folder = '/Users/hannah/OneDrive - Lancaster University/aurora/python_scripts/dataframes/'

df = pd.read_csv(root_folder+'main_regions_pixels_seen.csv',delimiter=',')
df.sort_values(by='CML', inplace=True)

cml = df['CML'].to_numpy()

#cml_new = df['CML'].to_numpy()

swirl_percentage = df['Percentage_Seen_Swirl_Region'].to_numpy()
noon_percentage = df['Percentage_Seen_Noon_Active_Region'].to_numpy()
dusk_percentage = df['Percentage_Seen_Dusk_Active_Region'].to_numpy()


'''
NOON
'''

lower_noon = []
lower_cml_noon = []
upper_noon = []
upper_cml_noon = [] 

# set all points above certain cml to 100%
new_noon = np.where((noon_percentage < 100) & (cml > 114), 100, noon_percentage)

for j in new_noon:
        if j < 100:
            lower_noon.append(j)
        else:
            upper_noon.append(j)
for i in cml:
        if i < 114:
            lower_cml_noon.append(i)
        else:
            upper_cml_noon.append(i)
        
# convert to arrays
lower_cml_noon = np.array(lower_cml_noon)
lower_noon = np.array(lower_noon)
upper_cml_noon = np.array(upper_cml_noon)
upper_noon = np.array(upper_noon)

#######################

'''
SWIRL
'''
      
# set all points above certain cml to 100%
new_swirl = np.where((swirl_percentage < 100) & (cml > 200), 100, swirl_percentage)


'''
FITTING
'''

a = np.polyfit(cml, dusk_percentage, 3)
aa = np.poly1d(a)

b = np.polyfit(cml, new_swirl, 5)
bb = np.poly1d(b)

e = np.polyfit(upper_cml_noon, new_noon, 1)
ee = np.poly1d(e)


# arrange cml more evenly for fit plotting
x_new = np.linspace(min(cml), max(cml), len(cml))
x_noon_u =  np.linspace(min(upper_cml_noon), max(upper_cml_noon), len(upper_cml_noon))

y_new_dusk = aa(x_new)
y_new_swirl = bb(x_new)
y_noon_u = ee(x_noon_u)


fig = plt.figure(figsize=(15,8))
plt.rcParams['font.size'] = '14'
ax = plt.subplot(1,1,1)

plt.scatter(cml,swirl_percentage,color='orange',s=10,label='Swirl Region')
plt.plot(x_new, y_new_swirl, color='red')

plt.scatter(cml,noon_percentage,color='blue',s=10,label='Noon Active Region')
plt.plot(x_noon_u, y_noon_u, color='darkblue')

plt.scatter(cml,dusk_percentage,color='limegreen',s=10,label='Dusk Active Region')
plt.plot(x_new, y_new_dusk, color='darkgreen')

plt.rcParams["legend.markerscale"] = 4
ax.legend(framealpha=0.5)
ax.set_xlabel('CML (\N{DEGREE SIGN})')
ax.set_ylabel('Percentage Seen (%)')
ax.tick_params(direction='in',bottom=True, top=True, left=True, right=True)

saveloc = ('/Users/hannah/OneDrive - Lancaster University/aurora/percentage_rmain_region.jpg')
plt.savefig(saveloc,bbox_inches='tight',dpi=400)

data = {'Dusk_Active_Fit': a, 'Noon_Active_Fit_Upper': e, 'Swirl_Fit': b}, #'Swirl_Fit_Upper': c} 

dictionary = dict(Dusk_Active_Fit=a, Noon_Active_Fit_Upper=e, Swirl_Fit=b)#, Swirl_Fit_Upper=c)
polynomials_df = pd.DataFrame({k : pd.Series(v) for k, v in dictionary.items()})

polynomials_df.to_csv(root_folder+'polynomials_main_regions.csv',index=False)
