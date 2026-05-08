"""
Created on Thu Apr 25 11:22:56 2024

@author: hannah

this script scales the power values up of a single visit (can input visit wanted at top) 
to appropriate values using coefficients calculated to compensate for the viewing angle of HST
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

visit='26'
regions = 'high' 
# standard - ['02','03','04','05','08','09','10','11','16','17','18','19','20','21','24','25','27','28','34','35']
# high - ['01','12','15','26']

root_folder = '/Users/hannah/OneDrive - Lancaster University/aurora/python_scripts/dataframes/'

powers_df = pd.read_csv(f'{root_folder}new_total_powers_regions_{visit}_CR_adj_26update.csv',delimiter=',')#{'root_folder}new_total_powers_regions_{visit}_march.csv', delimiter=',')#
polynomials_df = pd.read_csv(f'{root_folder}polynomials_{regions}_CR_adj.csv', delimiter=',')

cml = powers_df['CML'].to_numpy()

swirl_powers = powers_df['Total_Power_Swirl'].to_numpy()
dusk_active_powers = powers_df['Total_Power_Dusk_Active'].to_numpy()
noon_active_powers = powers_df['Total_Power_Noon_Active'].to_numpy()

dusk_coeff = polynomials_df['Dusk_Active_Fit'].to_numpy()

swirl_coeff = polynomials_df['Swirl_Fit'].to_numpy()
noon_coeff = polynomials_df['Noon_Active_Fit'].to_numpy()

# --------


noon_scale = []
for x in range(len(cml)):
    noon_eq = (noon_coeff[0]*(cml[x]))+ noon_coeff[1]
    noon_scaling = noon_eq/100
    noon_scale.append(noon_scaling)
    
    
if regions == 'high':
    dusk_scale = []
    for x in range(len(cml)):
        dusk_eq = (dusk_coeff[0]*(cml[x]**13)) + (dusk_coeff[1]*(cml[x]**12)) + (dusk_coeff[2]*(cml[x]**11))+(dusk_coeff[3]*(cml[x]**10))+ (dusk_coeff[4]*(cml[x]**9))+(dusk_coeff[5]*(cml[x]**8))+ (dusk_coeff[6]*(cml[x]**7))+(dusk_coeff[7]*(cml[x]**6))+(dusk_coeff[8]*(cml[x]**5))+(dusk_coeff[9]*(cml[x]**4))+(dusk_coeff[10]*(cml[x]**3))+(dusk_coeff[11]*(cml[x]**2))+(dusk_coeff[12]*(cml[x]))+(dusk_coeff[13])
        dusk_scaling = dusk_eq/100
        dusk_scale.append(dusk_scaling)
        
    swirl_scale = []
    for x in range(len(cml)):
        swirl_eq = (swirl_coeff[0]*(cml[x]**5)) + (swirl_coeff[1]*(cml[x]**4)) + (swirl_coeff[2]*(cml[x]**3))+(swirl_coeff[3]*(cml[x]**2))+(swirl_coeff[4]*(cml[x]))+(swirl_coeff[5])
        swirl_scaling = swirl_eq/100
        swirl_scale.append(swirl_scaling)
    
else:
    dusk_scale = []
    for x in range(len(cml)):
        #dusk_eq = (dusk_coeff[0]*(cml[x]**3)) + (dusk_coeff[1]*(cml[x]**2)) + (dusk_coeff[2]*(cml[x]))+(dusk_coeff[3])
        dusk_eq = (dusk_coeff[0]*(cml[x]**5)) + (dusk_coeff[1]*(cml[x]**4)) + (dusk_coeff[2]*(cml[x]**3)) + (dusk_coeff[3]*(cml[x]**2)) + (dusk_coeff[4]*(cml[x])) + (dusk_coeff[5])#*(cml[x])) #+ dusk_coeff[6]
        dusk_scaling = dusk_eq/100
        dusk_scale.append(dusk_scaling)
        
    swirl_scale = []
    for x in range(len(cml)):
        swirl_eq = (swirl_coeff[0]*(cml[x]**5)) + (swirl_coeff[1]*(cml[x]**4)) + (swirl_coeff[2]*(cml[x]**3))+(swirl_coeff[3]*(cml[x]**2))+(swirl_coeff[4]*(cml[x]))+(swirl_coeff[5])
        swirl_scaling = swirl_eq/100
        swirl_scale.append(swirl_scaling)
    # swirl_scale = []
    # for x in range(len(cml)):
    #     swirl_eq = (swirl_coeff[0]*(cml[x]**4)) + (swirl_coeff[1]*(cml[x]**3)) + (swirl_coeff[2]*(cml[x]**2))+(swirl_coeff[3]*(cml[x]))+(swirl_coeff[4])
    #     swirl_scaling = swirl_eq/100
    #     swirl_scale.append(swirl_scaling)
        
        

swirl_scale = np.array(swirl_scale)
dusk_scale = np.array(dusk_scale)
noon_scale = np.array(noon_scale)
        
        
# --------

new_swirl_powers = swirl_powers/swirl_scale
new_dusk_powers = dusk_active_powers/dusk_scale
new_noon_powers = noon_active_powers/noon_scale

swirl_mean = np.mean(new_swirl_powers)
dusk_mean = np.mean(new_dusk_powers)
noon_mean = np.mean(new_noon_powers)

var_swirl = []
for g in range(len(new_swirl_powers)):
    var = ((new_swirl_powers[g] - swirl_mean)**2 / len(new_swirl_powers))
    var_swirl.append(var)

#k = np.sum(var_swirl)
variance_swirl = np.sum(var_swirl)
sd_swirl = np.sqrt(variance_swirl)

# fix standard deviartion for all versions


var_noon = []
for w in range(len(new_noon_powers)):
    var = (new_noon_powers[w] - noon_mean)**2 
    var_noon.append(var)

variance_noon = np.mean(var_noon)
sd_noon = np.sqrt(variance_noon)

var_dusk = []
for u in range(len(new_dusk_powers)):
    var = (new_dusk_powers[u] - dusk_mean)**2 
    var_dusk.append(var)

variance_dusk = np.mean(var_dusk)
sd_dusk = np.sqrt(variance_dusk)



# dusk_err = new_dusk_powers*0.01
# dusk_err_2 = new_dusk_powers*(1-dusk_scale)
# dusk_error = np.sqrt((dusk_err*2)+(dusk_err_2**2))

# #noon_err = new_noon_powers*0.03
# #noon_error = noon_err

# swirl_err = new_swirl_powers*0.05
# swirl_err_2 = new_swirl_powers*(1-swirl_scale)
# swirl_error = np.sqrt((swirl_err**2)+(swirl_err_2**2))

if regions == 'high':
    swirl_error = (new_swirl_powers*0.05) / swirl_scale
    dusk_error = (new_dusk_powers*0.10) / dusk_scale
else:
    swirl_error = (new_swirl_powers*0.03) / swirl_scale
    dusk_error = (new_dusk_powers*0.01) / dusk_scale

# plotting 
fig = plt.figure(figsize=(16,12))
ax = plt.subplot(3,1,1)
#ax.plot(dates,medians,'.', markersize=5)
ax.plot(cml, swirl_powers, color='red', label='Swirl Powers')
ax.errorbar(cml, new_swirl_powers, yerr=swirl_error,fmt='.', markersize=6,color='lightgrey', label='Error')
ax.plot(cml, new_swirl_powers, color='darkorange', label='Corrected Swirl Powers')
#ax.set_ylim(6*10**10,10*10**11)
ax1 = ax.twinx()
ax.tick_params(direction='in',bottom=True, top=True, left=True, right=False)
ax1.plot(cml, swirl_scale, '--', color='lightgray', label='Scaling Factor')
ax.xaxis.set_major_locator(plt.MaxNLocator(20))
ax.legend()
ax.set_xlabel('CML (\N{DEGREE SIGN})')
ax.set_ylabel('Power (GW)')
ax1.set_ylabel('Scaling Factor')
ax1.set_ylim(0,1.2)
ax1.legend(bbox_to_anchor = [0.175 ,0.7])
ax1.tick_params(direction='in',bottom=True, top=True, left=False, right=True)

ax2 = plt.subplot(3,1,2)
#ax.plot(dates,medians,'.', markersize=5)
ax2.plot(cml, dusk_active_powers, color='darkolivegreen', label='Dusk Active Powers')
ax2.errorbar(cml, new_dusk_powers, yerr=dusk_error,fmt='.', markersize=6,color='lightgrey', label='Error')
ax2.plot(cml, new_dusk_powers, color='limegreen', label='Corrected Dusk Active Powers')
#ax2.set_ylim(6*10**10,7*10**12)
ax3 = ax2.twinx()
ax2.tick_params(direction='in',bottom=True, top=True, left=True, right=False)
ax2.xaxis.set_major_locator(plt.MaxNLocator(20))
ax2.legend()
ax3.plot(cml, dusk_scale, '--', color='lightgray', label='Scaling Factor')
ax2.set_xlabel('CML (\N{DEGREE SIGN})')
ax2.set_ylabel('Power (GW)')
ax3.set_ylabel('Scaling Factor')
ax3.set_ylim(0,1.2)
ax3.legend(bbox_to_anchor = [0.175 ,0.7])
ax3.tick_params(direction='in',bottom=True, top=True, left=False, right=True)

ax4 = plt.subplot(3,1,3)
#ax.plot(dates,medians,'.', markersize=5)
ax4.plot(cml,noon_active_powers, color='blue', label='Noon Active Powers')
#ax4.errorbar(cml, new_noon_powers, yerr=noon_error,fmt='.', markersize=6,color='lightgrey', label='Error')
ax4.plot(cml, new_noon_powers, color='indigo', label='Corrected Noon Active Powers')
ax4.tick_params(direction='in',bottom=True, top=True, left=True, right=False)
#ax2.set_ylim(6*10**10,7*10**12)
ax5 = ax4.twinx()
ax5.plot(cml, noon_scale, '--', color='lightgray', label='Scaling Factor')
ax4.xaxis.set_major_locator(plt.MaxNLocator(20))
ax4.legend()
ax4.set_xlabel('CML (\N{DEGREE SIGN})')
ax4.set_ylabel('Power (GW)')
ax5.set_ylabel('Scaling Factor')
ax5.set_ylim(0,1.2)
ax5.legend(bbox_to_anchor = [0.175 ,0.7])
ax5.tick_params(direction='in',bottom=True, top=True, left=False, right=True)


# save plot
saveloc = (f'/Users/hannah/OneDrive - Lancaster University/aurora/scaling_graph_CR_adj_{visit}.jpg')
plt.savefig(saveloc,bbox_inches='tight',dpi=400)

# save power data to dataframe
new_powers_df = pd.DataFrame()

new_powers_df = new_powers_df.assign(CML=cml)
new_powers_df = new_powers_df.assign(Corrected_Swirl_Powers=new_swirl_powers)
new_powers_df = new_powers_df.assign(Swirl_Error=swirl_error)
new_powers_df = new_powers_df.assign(Corrected_Dusk_Active_Powers=new_dusk_powers)
new_powers_df = new_powers_df.assign(Dusk_Active_Error=dusk_error)
new_powers_df = new_powers_df.assign(Corrected_Noon_Active_Powers=new_noon_powers)
#new_powers_df = new_powers_df.assign(Noon_Active_Error=noon_error)

new_powers_df.to_csv(f'{root_folder}new_normalised_powers_{visit}_CR_adj_26update.csv',index=False)
    