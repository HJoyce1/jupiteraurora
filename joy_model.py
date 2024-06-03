"""
Created on Tue Oct  3 11:02:06 2023

@author: hannah
"""
# this code is designed to determine the magnetopause location  and local time
# using the joy model (joy et al 2002)

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

'''
need to write a bit of code to extract pressure from data ie:
    
# read in solar wind paramters file
df = pd.read_csv('/Users/hannah/OneDrive - Lancaster University/aurora/solar_wind_data.csv',delimiter=',')

# need to convert columns to arrays to work with averages
ram_pressure_array = df[['RAM_PRESSURE_PROTONS_NPA']].to_numpy() # pressure 

ram_pressure_array[15000] # <-- need a way to get this location?

don't need this anymore provided function carries over bs_loc and mp_loc
'''

pressure = 0.05

# magnetopause location
def magnetopause_location(pressure,plotting):
    # coefficients, joy et al 2002
    a_mag = -0.134 + (0.488*(pressure**(-0.25)))
    b_mag = -0.581 - 0.225*((pressure**(-0.25)))
    c_mag = -0.186 - 0.016*((pressure**(-0.25)))
    d_mag = -0.014 + (0.096*pressure)
    e_mag = -0.814 - (0.811*pressure)
    f_mag = -0.050 + (0.168*pressure)
    
    # set values for x axis - in RJ
    x_plot = np.arange(-250,200,.01, dtype=float) # .01 for resolution purposes
    x_plot = np.divide(x_plot,120)
    
    # make values for quadratic to make curves
    b_plot = d_mag + (f_mag*x_plot)
    a_plot = e_mag
    c_plot = a_mag + (b_mag*x_plot) + (c_mag*(x_plot**2))

    # dawn side - negative values but +ve quadratic
    y_plot_plus = ((-1*b_plot) + np.sqrt((b_plot**2)-4*a_plot*c_plot))/(2*a_plot)
    # dusk side - postive values but -ve quadratic
    y_plot_minus = ((-1.*b_plot) - np.sqrt((b_plot**2)-4*a_plot*c_plot))/(2*a_plot)
    # rescale x and y to jovian radii as cacls assume R/120
    y_plot_plus = 120*y_plot_plus
    y_plot_minus = 120*y_plot_minus
    
    # radial distance / local time ??
    x_plot = x_plot*(120)
    
    # stick lines together to make parabola
    x = np.append(x_plot,np.flip(x_plot))
    y = np.append(y_plot_minus,np.flip(y_plot_plus))
   
    if plotting == 'yes':
        fig = plt.figure()
        ax1 = fig.add_subplot(1,1,1)
        ax1.plot(x,y)
        ax1.tick_params(which='both',direction='in',bottom=True, top=True, left=True, right=True)
        ax1.set_xlim(-250,150)
    
    # calculate location of y pointa where cuurves 'meet'
    # both of these should be the same value
    max_plus_y = np.nanargmax(y_plot_plus)
    min_minus_y = np.nanargmin(y_plot_minus)
    # if not the same, raise value error
    if max_plus_y != min_minus_y:
        raise ValueError
    
    ## can check values at position identified if wanted
    #print(y[min_minus_y])
    #print(y[max_plus_y])

    # identify x value at tip of curve (where y plots meet)
    # gives location of magnetpopause in RJ
    mp_loc = x[min_minus_y]
    # print(mp_loc)
    return mp_loc


def multi_nose(pressure):
    nose_locs_mp = []
    nose_locs_bs = []
    for k in pressure:
        pp =  magnetopause_location(k,'no')
        nose_locs_mp.append(pp)
        qq = bow_shock_location(k, 'no')
        nose_locs_bs.append(qq)
    return(nose_locs_mp), (nose_locs_bs)
        
    '''
    this one appends all mp_loc
    '''
    

# bow shock location
def bow_shock_location(pressure,plotting):
    # coefficients, joy et al 2002
    a_mag = -1.107 + (1.591*(pressure**(-0.25)))
    b_mag = -0.566 - 0.812*((pressure**(-0.25)))
    c_mag = 0.048 - 0.059*((pressure**(-0.25)))
    d_mag = 0.077 - (0.38*pressure)
    e_mag = -0.874 - (0.299*pressure)
    f_mag = -0.055 + (0.059*pressure)
    
    # set values for x axis - in RJ
    x_plot = np.arange(-250,200,.01, dtype=float) # .01 for resolution purposes
    x_plot = np.divide(x_plot,120)
    
    # make values for quadratic to make curves
    b_plot = d_mag + (f_mag*x_plot)
    a_plot = e_mag
    c_plot = a_mag + (b_mag*x_plot) + (c_mag*(x_plot**2))

    # dawn side - negative values but +ve quadratic
    y_plot_plus = ((-1*b_plot) + np.sqrt((b_plot**2)-4*a_plot*c_plot))/(2*a_plot)
    # dusk side - positive values but -ve quadratic 
    y_plot_minus = ((-1.*b_plot) - np.sqrt((b_plot**2)-4*a_plot*c_plot))/(2*a_plot)
    # rescale x and y to jovian radii as cacls assume R/120
    y_plot_plus = 120*y_plot_plus
    y_plot_minus = 120*y_plot_minus
    
    # radial distance / local time ??
    x_plot = x_plot*(120)
    
    # stick lines together to make parabola
    x = np.append(x_plot,np.flip(x_plot))
    y = np.append(y_plot_minus,np.flip(y_plot_plus))
   
    if plotting =='yes':
        fig = plt.figure()
        ax1 = fig.add_subplot(1,1,1)
        ax1.plot(x,y)
        ax1.tick_params(which='both',direction='in',bottom=True, top=True, left=True, right=True)
        ax1.set_xlim(0,150)

    # calculate location of y points where curves 'meet'
    # both of these should be the same value
    max_plus_y = np.nanargmax(y_plot_plus)
    min_minus_y = np.nanargmin(y_plot_minus)
    # if not the same, raise value error
    if max_plus_y != min_minus_y:
        raise ValueError
        
        '''
        look at equation in paper and do quadratic formula because y and z = 0
        '''
    
    ## can check values at position identified if wanted
    #print(y[min_minus_y])
    #print(y[max_plus_y])

    # identify x value at tip of curve (where y plots meet)
    # gives location of bow shock in RJ
    bs_loc = x[min_minus_y]
    #print(bs_loc)
    #breakpoint()
    return bs_loc

# to run function if ehecking from this file
# magnetopause_location(pressure,'yes')
# bow_shock_location(pressure,'yes')
