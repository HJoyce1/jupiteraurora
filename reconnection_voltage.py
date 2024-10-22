"""
Created on Fri Jan  5 11:19:52 2024

@author: hannah

module to calculate low latitude reconnection & high latitude reconnection
"""
    
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

def reconnection_voltages(clock_angle,B_perp,sw_speed,mp_loc_km,times):
    RJ = 71492 # km

    shifted_clock = []
    for p in range(len(clock_angle)):
        shift = clock_angle[p] + 90
        shifted_clock.append(shift)
    
    shift_angle_rad = np.deg2rad(shifted_clock)
    clock_angle_rad = np.deg2rad(clock_angle)
    
    LL_rec = []
    for q in range(len(clock_angle)):
        LL = (sw_speed[q]*10**3)*(B_perp[q]*10**(-9))*((mp_loc_km[q]*10**3)/2)*(np.cos(clock_angle_rad[q]/2)**4)
        LL_rec.append(LL)
        
    LL_rec = np.array(LL_rec)
    
    # # old - for cusp
    # HL_rec = []
    # for q in range(len(clock_angle)):
    #     HL = 0.5*(sw_speed[q]*10**3)*(B_perp[q]*10**(-9))*((mp_loc_km[q]*10**3)/2)*(np.sin(clock_angle_rad[q]/2)**4)
    #     HL_rec.append(HL)
    # HL_rec = np.array(HL_rec)
    
      
    HL_rec_pos = []
    for q in range(len(clock_angle)):
        HL = 0.5*(sw_speed[q]*10**3)*(B_perp[q]*10**(-9))*((mp_loc_km[q]*10**3)/2)*((np.sin(shift_angle_rad[q]/2))**4)
        HL_rec_pos.append(HL)
        
    HL_rec_pos = np.array(HL_rec_pos)
    
    HL_rec_neg = []
    for q in range(len(clock_angle)):
        HL = 0.5*(sw_speed[q]*10**3)*(B_perp[q]*10**(-9))*((mp_loc_km[q]*10**3)/2)*((np.cos(shift_angle_rad[q]/2))**4)
        HL_rec_neg.append(HL)
        
    HL_rec_neg = np.array(HL_rec_neg)
    
    HL_rec_abs = []
    for q in range(len(clock_angle)):
        HL = 0.5*(sw_speed[q]*10**3)*(B_perp[q]*10**(-9))*((mp_loc_km[q]*10**3)/2)*(abs((np.sin(clock_angle_rad[q])/2)**4))
        HL_rec_abs.append(HL)
        
    HL_rec_abs = np.array(HL_rec_abs)
    
    
    #plotting currently causing python to run out of memory, will need to look at this
    fig = plt.figure(figsize=(30,20))
    ax7 = plt.subplot(4,1,1)
    ax7.plot(times,(LL_rec/1e6), '.', markersize=0.2) #times[:end_first]
    #ax7.plot(times[mag_second:],LL_rec_b, '.', markersize=0.2)
    ax7.set_xlabel('Time (DOY)')
    ax7.set_ylabel('Low-Latitude Reconnection Voltage (MV)')
    ax7.xaxis.set_major_formatter(mdates.DateFormatter('%j'))
    ax7.xaxis.set_major_locator(plt.MaxNLocator(15))
    
    ax8 = plt.subplot(4,1,2)
    ax8.plot(times,(HL_rec_pos/1e6), '.', markersize=0.2) #times[:end_first]
    #ax7.plot(times[mag_second:],LL_rec_b, '.', markersize=0.2)
    ax8.set_xlabel('Time (DOY)')
    ax8.set_ylabel('High-Latitude (+By) Reconnection Voltage (MV)')
    ax8.xaxis.set_major_formatter(mdates.DateFormatter('%j'))
    ax8.xaxis.set_major_locator(plt.MaxNLocator(15))
    
    ax1 = plt.subplot(4,1,3)
    ax1.plot(times,(HL_rec_neg/1e6), '.', markersize=0.2) #times[:end_first]
    #ax7.plot(times[mag_second:],LL_rec_b, '.', markersize=0.2)
    ax1.set_xlabel('Time (DOY)')
    ax1.set_ylabel('High-Latitude (-By) Reconnection Voltage (MV)')
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%j'))
    ax1.xaxis.set_major_locator(plt.MaxNLocator(15))
    
    ax2 = plt.subplot(4,1,4)
    ax2.plot(times,(HL_rec_abs/1e6), '.', markersize=0.2) #times[:end_first]
    #ax7.plot(times[mag_second:],LL_rec_b, '.', markersize=0.2)
    ax2.set_xlabel('Time (DOY)')
    ax2.set_ylabel('High-Latitude Reconnection Voltage (MV)')
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%j'))
    ax2.xaxis.set_major_locator(plt.MaxNLocator(15))
    
    saveloc = ('/Users/hannah/OneDrive - Lancaster University/aurora/high_lat_volt_by_abssin.jpg')
    plt.savefig(saveloc,bbox_inches='tight',dpi=400)

    return LL_rec, HL_rec_pos, HL_rec_neg
