import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cmocean.cm as cmo
from scipy.optimize import fmin

from .tools import plog
rmsd = lambda x,y : np.sqrt(np.nanmean( (x-y)**2 ))

def getGeoMagStrength(ADCP):
    lat = np.nanmedian(ADCP.Latitude.values)
    lon = np.nanmedian(ADCP.Longitude.values)
    date = pd.to_datetime(np.nanmedian(ADCP.time.values.astype(float)))
    year = date.year
    month = date.month
    day = date.day

    url = str('https://geomag.bgs.ac.uk/web_service/GMModels/igrf/14/?'+
          'latitude='+str(lat)+'&longitude='+str(lon)+
          '&date='+str(year)+'-'+str(month)+'-'+str(day)+
          '&resultFormat=xml')
    import urllib
    import xml.etree.ElementTree as ET
    with urllib.request.urlopen(url) as resp:
        xml_bytes = resp.read()
    root = ET.fromstring(xml_bytes)

    total_intensity = float(root.find('field-value/total-intensity').text) * 1e-9 * 10000 * 1000 # To tesla, then to gauss then to millgauss
    declination = float(root.find('field-value/declination').text)
    plog("Collecting IGRF14 Model data:")
    plog(f"   Total intensity (nT): {total_intensity} milligauss")
    plog(f"   Declination (deg east): {declination}")

    ADCP['Magnetic_total_intensity_nT'] = total_intensity
    ADCP['Magnetic_declination'] = declination
    
    return ADCP

def add_declination(ADCP):
    if 'Magnetic_declination' not in ADCP:
        ADCP = getGeoMagStrength(ADCP)
    
    ADCP['Heading'] = ADCP['Heading'] + ADCP['Magnetic_declination']
    return ADCP

def correct_heading(ADCP, options):    
    if "Heading_uncorrected" in ADCP:
        ADCP["Heading"] = ("time", ADCP["Heading_uncorrected"].values)
        plog('Reverting to "Heading_uncorrected" variable before recalibration.')

    plog("Beginning compass calibration.")
    ## STEP 1 : collect magnetic target
    # Technically not that big a deal unless covering multiple different areas, in which case we should change the fielf over time, which we don't. Basically, any old constant will do.
    # Get local geomagnetic target strength:
    ADCP = getGeoMagStrength(ADCP)
    target = ADCP['Magnetic_total_intensity_nT'].values
    
    ## STEP 2 : Extract ADCP magnetometer data
    if options['ADCP_mounting_direction'] == 'top':
        sign = -1
    elif options['ADCP_mounting_direction'] == 'bottom':
        sign = 1
    
    MagX = ADCP['MagnetometerX'].values
    MagY = ADCP['MagnetometerY'].values
    MagZ = ADCP['MagnetometerZ'].values
    
    AccX = ADCP['AccelerometerX'].values
    AccY = ADCP['AccelerometerY'].values
    AccZ = ADCP['AccelerometerZ'].values
        
    roll = ADCP['Roll'].values
    pitch = ADCP['Pitch'].values
    
    
    ## STEP 3 : define lambda functions
    rmsd    = lambda x,y,z : np.sqrt( np.mean( ( norm(x[idx],y[idx],z[idx]) - target)**2 ) ) # Score function for Hard iron # Notice the [idx] here which removes bad points from the score function!
    rvar    = lambda x,y,z : np.nanvar(norm(x[idx],y[idx],z[idx])) # Score function for soft iron # Notice the [idx] here which removes bad points from the score function!
    norm    = lambda x,y,z : np.sqrt(x**2 + y**2 + z**2) 

    def circ(x):
        x[x < -180] = x[x < -180] + 360
        x[x >  180] = x[x >  180] - 360
        return x
    
    cosd    = lambda x : np.cos(np.deg2rad(x))
    sind    = lambda x : np.sin(np.deg2rad(x))
    atan2d  = lambda x,y : np.rad2deg(np.arctan2(x,y))
    rot_x   = lambda x,y,z : x*cosd(pitch) + y*sind(roll)*sind(pitch) + z*cosd(roll)*sind(pitch)
    rot_y   = lambda x,y,z : y*cosd(roll) - z*sind(roll)
    wrap    = lambda x : (x+360)%360
    heading = lambda x,y,z : wrap( atan2d(rot_x(x,sign * y,sign * z), rot_y(x,sign * y,sign * z)) - 90 )


    # Offset only calibration
    def calibrate_offset(x,y,z,coeffs):
        return x-coeffs[0], y-coeffs[1], z-coeffs[2]
    def minimisation_offset(coeffs):
        x,y,z = calibrate_offset(MagX,MagY,MagZ,coeffs)
        return rvar(x,y,z)

    # Symmetrical calibration functions
    def calibrate_circularise_symmetrical(x,y,z,offset,sphere):
        A = np.reshape(symsphere(sphere),(3,3))
        B = offset
        out = A @ np.array([x-B[0], y-B[1], z-B[2]])
        return out[0,:],out[1,:],out[2,:]
    def minimisation_circularise_symmetrical(coeffs, *args):
        if len(args) == 1:
            offset = args[0]
        else:
            offset = [0,0,0]
        x,y,z = calibrate_circularise_symmetrical(MagX,MagY,MagZ,offset,coeffs)
        return rmsd(x,y,z)
    symsphere = lambda coeffs : [coeffs[0], coeffs[3], coeffs[4], 
                                  coeffs[3], coeffs[1], coeffs[5],
                                  coeffs[4], coeffs[5], coeffs[2]]
    
    # Fully spherical calibration functions
    def calibrate_circularise_full(x,y,z,coeffs):
        A = np.reshape(fullsphere(coeffs[:9]),(3,3))
        B = coeffs[-3:]
        out = A @ np.array([x-B[0], y-B[1], z-B[2]])
        return out[0,:],out[1,:],out[2,:]
        
    def minimisation_circularise_full(coeffs):
        x,y,z = calibrate_circularise_full(MagX,MagY,MagZ,coeffs)
        return rmsd(x,y,z)
        
    fullsphere = lambda coeffs : [coeffs[0], coeffs[3], coeffs[4], 
                                  coeffs[6], coeffs[1], coeffs[5],
                                  coeffs[7], coeffs[8], coeffs[2]]
    
    
    idx = np.full(len(MagX),True)
    idx = (norm(MagX,MagY,MagZ) < 9999) & ...
        (np.abs(norm(AccX,AccY,AccZ) - 1) < 0.001) & ...
        (np.abs(norm(MagX,MagY,MagZ) - np.nanmedian(norm(MagX,MagY,MagZ))) < 3*np.nanstd(norm(MagX,MagY,MagZ)))
    
    iterations = 1
    split = False
    
    for iteration in range(iterations):        
        if split:
            if iteration > 0:
                MagX,MagY,MagZ = calibrate_circularise_symmetrical(MagX, MagY, MagZ, offset, sphere)
                
            ## Calculate hard iron coefficients
            offset = fmin(minimisation_offset,np.array([0,0,0]), xtol=0.000001, ftol=0.000001)
            plog(f'Offsets : {offset}')

            ## Calculate ellipsoid coefficients
            sphere = fmin(minimisation_circularise_symmetrical,np.array([1,1,1,0,0,0]), args=(offset,), xtol=0.000001, ftol=0.000001)
            plog(f'Spheroid adjustment : {sphere}')
            
            magx,magy,magz = calibrate_circularise_symmetrical(MagX, MagY, MagZ, offset, sphere)
        else:        
            if iteration > 0:
                MagX,MagY,MagZ = calibrate_circularise_full(MagX, MagY, MagZ, coeffs)
            coeffs = fmin(minimisation_circularise_full,np.array([1,1,1,0,0,0,0,0,0,0,0,0]), xtol=0.000001, ftol=0.000001)
            magx,magy,magz = calibrate_circularise_full(MagX, MagY, MagZ, coeffs)
            plog(f'Offsets : {coeffs[-3:]}')
            plog(f'Spheroid adjustment : {np.reshape(fullsphere(coeffs[:9]),(3,3))}')
        
        

        cal_heading = heading(magx, magy, magz)
        mag_bias = norm(magx,magy,magz)-target


        ## Figures
        fig = plt.figure(figsize=(15,5))

        plt.subplot(141)
        _ = plt.hist(circ(cal_heading - ADCP.Heading.values), bins=np.arange(-20,20,1))
        plt.axvline(0, color='k')
        plt.xlabel('Heading correction')

        plt.subplot(142)
        plt.axvline(0,color='k')
        plt.axhline(0,color='k')
        idz = ~np.isnan(MagZ)
        h1 = plt.scatter(np.sqrt(MagY[idz]**2 + MagZ[idz]**2), MagX[idz], 1, 'r')
        h2 = plt.scatter(np.sqrt(MagY[idx&idz]**2 + MagZ[idx&idz]**2), MagX[idx&idz], 1, 'k')
        h3 = plt.scatter(np.sqrt(magy[idx&idz]**2 + magz[idx&idz]**2),magx[idx&idz], 1, 'c')
        plt.xlabel('norm(MagY, MagZ)')
        plt.ylabel('MagX')
        plt.axis('square')

        plt.subplot(143)
        plt.axvline(0,color='k')
        plt.axhline(0,color='k')
        idz = ~np.isnan(MagZ)
        h1 = plt.scatter(np.sqrt(MagX[idz]**2 + MagZ[idz]**2), MagY[idz], 1, 'r')
        h2 = plt.scatter(np.sqrt(MagX[idx&idz]**2 + MagZ[idx&idz]**2), MagY[idx&idz], 1, 'k')
        h3 = plt.scatter(np.sqrt(magx[idx&idz]**2 + magz[idx&idz]**2),magy[idx&idz], 1, 'c')
        plt.xlabel('norm(MagX, MagZ)')
        plt.ylabel('MagY')
        plt.axis('square')

        plt.subplot(144)
        plt.axvline(0,color='k')
        plt.axhline(0,color='k')
        idz = ~np.isnan(MagZ)
        h1 = plt.scatter(np.sqrt(MagX[idz]**2 + MagY[idz]**2), MagZ[idz], 1, 'r')
        h2 = plt.scatter(np.sqrt(MagX[idx&idz]**2 + MagY[idx&idz]**2), MagZ[idx&idz], 1, 'k')
        h3 = plt.scatter(np.sqrt(magx[idx&idz]**2 + magy[idx&idz]**2),magz[idx&idz], 1, 'c')
        plt.xlabel('norm(MagX, MagY)')
        plt.ylabel('MagZ')
        plt.axis('square')
        fig.legend(handles=[h1,h2,h3],labels=['Discarded values','Raw values','Corrected values'], loc='outside upper center')

    ADCP["Heading_uncorrected"] = ("time", ADCP["Heading"].values)
    ADCP["Heading"] = ("time", cal_heading)
    
    plog("Preserving original heading as \"Heading_uncorrected\".")
    plog("Compass calibration complete.")
    return ADCP