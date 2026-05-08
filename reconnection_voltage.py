#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  5 11:19:52 2024

@author: hannah

module to calculate low latitude reconnection & high latitude reconnection
"""
    
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

def reconnection_voltages(clock_angle,B_perp,sw_speed,mp_loc_km,times):

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
    
    # HL_rec_abs = []
    # for q in range(len(clock_angle)):
    #     HL = 0.5*(sw_speed[q]*10**3)*(B_perp[q]*10**(-9))*((mp_loc_km[q]*10**3)/2)*((np.sin(clock_angle_rad[q]))**4)
    #     HL_rec_abs.append(HL)
        
    # HL_rec_abs = np.array(HL_rec_abs)

    return LL_rec, HL_rec_pos, HL_rec_neg


def reconnection_gershman(Bx_sw,By_sw,Bz_sw,velocity_sw,pressure_sw,mp_loc_km,B_perp,clock_angle):
    
    # ------- functions -----
    def clock_format_convert(clock_angle):
        new_clock_angle = []
        for i in range(len(clock_angle)):
            if clock_angle[i] < 0:
                clock = clock_angle[i] + 360
                new_clock_angle.append(clock)
            else:
                new_clock_angle.append(clock_angle[i])
        return new_clock_angle

    def total_mag(Bx,By,Bz):
        Btot = np.sqrt((Bx**2)+(By**2)+(Bz**2))
        return Btot

    # ------ constants ------

    # constant parameters
    RJ = 71492 # jupiter radius in km
    RJ_m = RJ * 1e3# km
    mp = 1.67262192 * 1e-27 # mass of proton
    q = 1.602 * 1e-19 # charge of proton
    dF = 0.85 # factor associated with plama depletion layer formation at subsolar magnetopause (Gershman 2024, Masters 2018)
    eta = 4 # shock compression factor
    #kappa = 0.881 # derived from hydrodynamic modelling, Spreiter and Alkasne 1970
    mu_0 = 4*np.pi*1e-7 # permeability of free space


    # ----- unit conversions ------

    velocity = velocity_sw * 1e3 # kms-1 > ms-1
    pressure = pressure_sw * 1e-9 # nPa -> Pa
    Bx = Bx_sw * 1e-9
    By = By_sw * 1e-9
    Bz = Bz_sw * 1e-9

    mp_loc_m = mp_loc_km * 1e3

    mp_loc_RJ = []
    for a in range(len(mp_loc_km)):
        mp_loc = mp_loc_m[a] / RJ_m # km cancel out
        mp_loc_RJ.append(mp_loc)
    mp_loc_RJ = np.array(mp_loc_RJ)

    # ------ shear angle -------

    # empty grid for clock angle
    clock_angle = []
    clock_angle_deg = []
    # calculate clock angle for each moving avg
    for j in range(len(Bz)):
    # use arctan2 electric boogaloo
    # 0 is Bz+, +/-180 is Bz-, +90 is By+, -90 is By-
        theta_c = np.arctan2(By[j],Bz[j])
        theta_c_deg = np.degrees(theta_c)
        clock_angle.append(theta_c)
        clock_angle_deg.append(theta_c_deg)
    clock_angle = np.array(clock_angle)
    clock_angle_deg = np.array(clock_angle_deg)

    shear_angle = []
    clock = clock_format_convert(clock_angle_deg)
    clock_to_shear = [x+180 for x in clock]
    shear_ang = np.array(clock_to_shear)
    
    for g in range(len(Bz_sw)):
        shear_clock = shear_ang[g]
        # if loc == 'low' or loc == 'neg_By':
        if shear_clock >= 360:
            shear = shear_clock - 360
        else:
            shear = shear_clock
        shear_angle.append(shear)
        
    shear_angle = np.array(shear_angle)


    # --------- from upstream parameters -------

    Btot = []
    for b in range(len(Bz)):
        Btots = total_mag(Bx[b],By[b],Bz[b])
        Btot.append(Btots)
    Btot = np.array(Btot)

    # get upstream mass density from velocity and pressure values
    # pressure = 0.5 * mass density * velocity**2
    mass_density_u = []
    for i in range(len(pressure)):
        density = (pressure[i])/(velocity[i]**2)
        mass_density_u.append(density)
    mass_density_u = np.array(mass_density_u)

    B_2 = 1e-9
    B_1 = Btot#10e-9
    B_asym = []
    for f in range(len(Btot)):
        B = (2*B_1[f]*B_2) / (B_1[f]+B_2)
        B_asym.append(B)
        
    #B_asym = (2*B_1*B_2) / (B_1+B_2)

    # # B/√mu*rho
    # V_A = []
    # for j in range(len(Btot)):
    #     V = Btot[j]/(np.sqrt(mu_0*mass_density_u[j]))
    #     V_A.append(V)
    # V_A = np.array(V_A)

    # M_A = []
    # for k in range(len(V_A)):
    #     M = velocity[k] / V_A[k]
    #     M_A.append(M)
    # M_A = np.array(M_A)

    mass_density_1 = []
    for o in range(len(mass_density_u)):
        mass_den = dF*eta*mass_density_u[o]
        mass_density_1.append(mass_den)
    mass_density_1 = np.array(mass_density_1)


    # ------ number density -> mass density -------
    
    mass_ratio = 25 * mp
    number_density_2_cm = []
    for qq in range(len(mp_loc_RJ)):
        number_d = (1987*((mp_loc_RJ[qq]/6)**(-8.2))) + (14*((mp_loc_RJ[qq]/6)**(-3.2))) + (0.05*((mp_loc_RJ[qq]/6)**(-0.65)))
        number_density_2_cm.append(number_d)
    number_density_2_cm = np.array(number_density_2_cm)
    number_density_2 = number_density_2_cm * 1e6 # convert from cm^-3 to m^-3

    mass_density_2 = number_density_2 * mass_ratio


    # ------ mass density asym -----
    
    mass_density_asym = []
    for p in range(len(mass_density_1)):
        asym = ((mass_density_1[p]*B_2) + (mass_density_2[p]*B_1[p]))
        mass_asym = asym / (B_1[p] + B_2)
        mass_density_asym.append(mass_asym)
    mass_density_asym = np.array(mass_density_asym)

    di = []
    mq = (32*mp)/(2*q)
    for r in range(len(mass_density_asym)):
        d = mq*(1/(np.sqrt(mu_0*mass_density_asym[r])))
        di.append(d)
    di = np.array(di)


    # ----------- Alfven Velocity asym  -----------

    V_A_asym = []
    for s in range(len(mass_density_asym)):
        Va = B_asym[s] / (np.sqrt(mu_0*mass_density_asym[s]))
        V_A_asym.append(Va)
    V_A_asym = np.array(V_A_asym)


    # --------- final calculation for reconnection power, in w -------
        
    shear_angle_rad = np.deg2rad(shear_angle)
    rec_volt = []
    for w in range(len(shear_angle)):
        voltage =  3.6*1e8 * V_A_asym[w] * B_asym[w]**2 * di[w] * mp_loc_m[w] * np.sin(shear_angle_rad[w]/2)**3
        rec_volt.append(voltage)
    rec_volt = np.array(rec_volt)
    
    return rec_volt

def kelvin_helmholtz(velocity_sw,pressure_sw,mp_loc_km):
    
    # -------------- constant paramters -------------
    
    RJ = 71492 # jupiter radius in km
    RJ_m = RJ *10**3
    Ly = 4*RJ_m  # region helvin helmholtz can occur

    # -------------- unit conversions -------------
    
    velocity = velocity_sw * 1e3 # kms-1 > ms-1
    pressure = pressure_sw * 1e-9 # nPa -> Pa
    mp_loc_m = mp_loc_km * 1e3

    mp_loc_RJ = []
    for a in range(len(mp_loc_km)):
        mp_loc = mp_loc_m[a] / RJ_m # km cancel out
        mp_loc_RJ.append(mp_loc)
    mp_loc_RJ = np.array(mp_loc_RJ)
    
    # -------------- scale height calculation ------------

    # from Bagnal and Delamere 2010
    scale_h = []
    # h = a1 + a2r + a3r**2 + a4r**3 + a5r**4
    # r = log10(R/6)
    for l in range(len(mp_loc_RJ)):
        scale = -0.116 + (2.14*(np.log10(mp_loc_RJ[l]/6))) -(2.05*(np.log10(mp_loc_RJ[l]/6))**2) + (0.491*(np.log10(mp_loc_RJ[l]/6))**3) + (0.126 * (np.log10(mp_loc_RJ[l]/6))**4)
        scale_h.append(scale)
    scale_h = np.array(scale_h)
    
    Lx = (10**scale_h) * RJ_m *2  # finalise scale height as variable in equation
    

    # ------------ mass density calculation ------------

    # pressure = 0.5 * mass density * velocity**2
    mass_density_u = []
    for i in range(len(pressure)):
        density = ((pressure[i])/(velocity[i]**2))
        mass_density_u.append(density)
    mass_density_u = np.array(mass_density_u)

    vel_shear = velocity/2 # velocity shear

    
    # ---------- kelvin helmholtz calculation ------------

    kelv_helm = []
    for k in range(len(mass_density_u)):
        kh = 0.022*mass_density_u[k]*(vel_shear[k]**3)*Ly*Lx[k]
        kelv_helm.append(kh)
    kelv_helm = np.array(kelv_helm)
    
    return kelv_helm

def kelvin_helmholtz_dd(velocity_sw,pressure_sw,mp_loc_km):
    
    # -------------- constant paramters -------------
    
    RJ = 71492 # jupiter radius in km
    RJ_m = RJ *10**3
  
    # -------------- unit conversions -------------
    
    velocity = velocity_sw * 1e3 # kms-1 > ms-1
    pressure = pressure_sw * 1e-9 # nPa -> Pa
    mp_loc_m = mp_loc_km * 1e3
    Ly = 2*mp_loc_m # region helvin helmholtz can occursmall


    mp_loc_RJ = mp_loc_m / RJ_m # km cancel out
    
    # -------------- scale height calculation ------------

    # from Bagnal and Delamere 2010
    scale_h = []
    # h = a1 + a2r + a3r**2 + a4r**3 + a5r**4
    # r = log10(R/6)
    #log_term = np.log10(mp_loc_RJ / 6)
    for l in range(len(mp_loc_RJ)):
        scale = -0.116 + (2.14*(np.log10(mp_loc_RJ[l]/6))) -(2.05*(np.log10(mp_loc_RJ[l]/6))**2) + (0.491*(np.log10(mp_loc_RJ[l]/6))**3) + (0.126 * (np.log10(mp_loc_RJ[l]/6))**4)
        scale_h.append(scale)
    scale_h = np.array(scale_h)
    
    Lx = (10**scale_h) * RJ_m *2  # finalise scale height as variable in equation
    

    # ------------ mass density calculation ------------

    # pressure = mass density * velocity**2
    # mass_density_u = []
    # for i in range(len(pressure)):
    #     density = ((pressure[i])/(velocity[i]**2))
    #     mass_density_u.append(density)
    # mass_density_u = np.array(mass_density_u)
    
    mass_density_u = pressure / velocity**2
    mass_den = (3/2)*mass_density_u


    # ----- velocity shear calculation ----------------
    
    vel_shear_dawn = []
    vel_shear_dusk = []# *1e3 #(velocity +(150*1e3))#velocity/2 #250 * 10e3
    for j in range(len(velocity)):
        dawn = (350*1e3) - (0.5 * -velocity[j]) # -- = +
        dusk = -(350*1e3) - (0.5 * -velocity[j]) 
        vel_shear_dawn.append(dawn)
        vel_shear_dusk.append(dusk)
    vel_shear_dawn = np.array(vel_shear_dawn)
    vel_shear_dusk = np.array(vel_shear_dusk)

    
    # ---------- kelvin helmholtz calculation ------------

    kelv_helm_dawn = []
    for k in range(len(mass_density_u)):
        kh1= 0.022*mass_den[k]*(vel_shear_dawn[k]**3)*Ly[k]*Lx[k]
        kelv_helm_dawn.append(kh1)
    kelv_helm_dawn = np.array(kelv_helm_dawn)
    
    kelv_helm_dusk = []
    for m in range(len(mass_density_u)):
        kh2 = 0.022*mass_den[m]*(abs(vel_shear_dusk[m]**3))*Ly[m]*Lx[m]
        kelv_helm_dusk.append(kh2)
    kelv_helm_dusk = np.array(kelv_helm_dusk)
        
    return kelv_helm_dawn, kelv_helm_dusk