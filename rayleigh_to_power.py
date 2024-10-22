"""
Created on Fri Mar  1 14:34:46 2024

@author: hannah

this script ccontains a module that converts intensities (in kR) to power in (W) as well as returning the 
area values of the apparent area (seen by HST) and absolute area of planet
"""

def rayleigh_to_power(intensity, image, visit):
    import numpy as np
    
    # if intensity has not been converted...
    # need to flip on horizontal axis to convert to other-handed co-ords
    intensity = np.flip(intensity,1)
    intensity = np.where(intensity <= 0, np.nan, intensity)
    
    
    radius_e = (71942 + 240) #*10**3
    radius_p = (66854 + 240) #*10**3
    
# =============================================================================
#   first get absolute area of each pixel - dependent on colatitude only
# =============================================================================

    # colatitudes
    colat = np.linspace(0.25, 180, num=720)
    colat = colat/180*np.pi # radians conversion
    
    # need to get radius at pixel location
    ellipticity = ((radius_e/radius_p)**2) - 1
    
    theta = np.radians(0.25)
    
    abs_area = []
    for clat in range(len(colat)):
        r_to_surface = radius_e/np.sqrt(1+((ellipticity)*(np.cos(colat[clat]))**2))
    
        area = (r_to_surface**2) * np.sin(colat[clat]) * (theta*theta)
        abs_area.append(area)
        
# =============================================================================
#   next get apparent area as seen by detector
# =============================================================================

    loaded_array = np.load(f'/Users/hannah/OneDrive - Lancaster University/aurora/python_scripts/elev_angles/elevation_angles_visit_{visit}.npz',allow_pickle = True)
    angles_all = loaded_array['arr_0']
    angles = angles_all[:,:,image]
        
    # longitudes
    lons = np.linspace(0, 360, num=1441)
    lons = lons[:-1]
    lons = lons/180*np.pi
    
    area_lon_mesh = np.meshgrid(abs_area, lons)
    area_mesh = area_lon_mesh[0]
    
    # latitudes
    lats = np.linspace(-90, 90, num=721)
    lats = lats[:-1]
    lats = lats/180*np.pi
            
    app_area = np.multiply(area_mesh, angles)
        
    app_area = np.flip(app_area,axis=0) # vertical flip bc intensities got flipped

# =============================================================================
#   now can convert intensities into auroral power
# =============================================================================

    radiance = (10**13/(4*np.pi))*intensity
    flux = radiance*4*np.pi*(app_area.T * 1*10**6)
    power = flux*(1.6*10**-18)
    
    power = np.flip(power, axis=1) # flip back to stay in same handed system as gif 
    
    return power, app_area, abs_area


def rayleigh_to_power_noarea(intensity, area):
    # if don't need to do area calculation - not sure what I made this bit for?
    import numpy as np
    
    intensity_nan = np.where(intensity <= 0, np.nan, intensity)
    
    radiance = (10**13/(4*np.pi))*intensity_nan
    flux = radiance*4*np.pi*(area.T * 1*10**6)
    power = flux*(1.6*10**-18)
    
    return power
