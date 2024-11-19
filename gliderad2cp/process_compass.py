####################################################
# CORRECT ADCP HEADING WITHOUT GNSS
####################################################
from glob import glob
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import cmocean.cm as cmo

from scipy.optimize import fmin 
from tqdm import tqdm

from datetime import timedelta
from scipy.interpolate import interp1d, interp2d


def interpt(x,y,xi):
    _gg = np.isfinite(x.astype('float')+y)
    return interp1d(x[_gg].astype('float'), y[_gg], bounds_error=False, fill_value=np.NaN)(xi.astype('float'))
def align_time(y,coeffs):
    time = tax.astype('float')/1e9
    aligned = interpt(time,y,time+coeffs[0])    
    return aligned
rmsd = lambda x,y : np.sqrt(np.nanmean( (x-y)**2 ))

def process(ADCP):    
    
    ## STEP 1 : collect magnetic target
    # Technically not that big a deal unless covering multiple different areas, in which case we should change the fielf over time, which we don't. Basically, any old constant will do.
    # Get local geomagnetic target strength:
    def getGeoMagStrength():    
        lat = 57.682374
        lon = 11.888402
        date = pd.to_datetime(tax[0])
        year = date.year
        month = date.month
        day = date.day
        
        url = str('http://geomag.bgs.ac.uk/web_service/GMModels/igrf/13/?'+
              'latitude='+str(lat)+'&longitude='+str(lon)+
              '&date='+str(year)+'-'+str(month)+'-'+str(day)+
              '&resultFormat=csv')
        import urllib
        magdata = urllib.request.urlopen(url)
        
        string = 'empty'
        while not not string:
            out = magdata.readline().decode("utf-8")
            if 'total-intensity units="nT"' in out:
                string = out
                print(out)
                break
        target = float(string.split('>')[1].split('<')[0])
        nT2milligauss = 10**-9 * 10000 * 1000 # To tesla, then to gauss then to millgauss
        print('Target = '+str(target*nT2milligauss))
        return target*nT2milligauss
    
    target = getGeoMagStrength()
    
    
    ## STEP 2 : Extract ADCP magnetometer data
    if top_mounted:
        sign = -1
    else:
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
    

    
    def calibrate_offset(x,y,z,coeffs):
        return x-coeffs[0], y-coeffs[1], z-coeffs[2]
    def minimisation_offset(coeffs):
        x,y,z = calibrate_offset(MagX,MagY,MagZ,coeffs)
        return rvar(x,y,z)
    
    
    
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
    idx = (norm(MagX,MagY,MagZ) < 9999) & (np.abs(norm(AccX,AccY,AccZ) - 1) < 0.001) 
    
    iterations = 1
    split = True
    
    for iteration in range(iterations):
        
        
        if split:
            if iteration > 0:
                MagX,MagY,MagZ = calibrate_circularise_symmetrical(MagX, MagY, MagZ, offset, sphere)
                
            ## Calculate hard iron coefficients
            offset = fmin(minimisation_offset,np.array([0,0,0]), xtol=0.000001, ftol=0.000001)
            print('Offsets : ',offset)

            ## Calculate ellipsoid coefficients
            sphere = fmin(minimisation_circularise_symmetrical,np.array([1,1,1,0,0,0]), args=(offset,), xtol=0.000001, ftol=0.000001)
            print('Spheroid adjustment : ',sphere)
            
            magx,magy,magz = calibrate_circularise_symmetrical(MagX, MagY, MagZ, offset, sphere)
        else:        
            if iteration > 0:
                MagX,MagY,MagZ = calibrate_circularise_full(MagX, MagY, MagZ, coeffs)
            coeffs = fmin(minimisation_circularise_full,np.array([1,1,1,0,0,0,0,0,0]), xtol=0.000001, ftol=0.000001)
            magx,magy,magz = calibrate_circularise_full(MagX, MagY, MagZ, coeffs)
        
        

        cal_heading = heading(magx, magy, magz)
        
        mag_bias = norm(magx,magy,magz)-target


        ## Figures
        plt.figure(figsize=(15,15))

        plt.subplot(421)
        plt.plot(ADCP['time'],circ(cal_heading - ADCP.Heading.values),'-k')
        plt.plot(ADCP['time'][idx],circ(cal_heading[idx] - ADCP.Heading.values[idx]),'-c')
        plt.plot(ADCP['time'],circ(ADCP.Heading.values - heading(MagX, MagY, MagZ)),'-r.',alpha=0.3)
        plt.ylim([-20,20])

        plt.subplot(422)
        plt.plot(ADCP['time'],norm(AccX,AccY,AccZ),'-k')
        plt.plot(ADCP['time'][idx],norm(AccX,AccY,AccZ)[idx],'-b')
        plt.ylim(1 - np.array([-1,1])*0.05)

        plt.subplot(434)
        plt.plot(ADCP['time'],norm(MagX, MagY, MagZ),'-r')
        plt.plot(ADCP['time'][idx],norm(MagX, MagY, MagZ)[idx],'-k')
        plt.plot(ADCP['time'][idx],norm(magx,magy,magz)[idx],'-c')
        plt.axhline(target)

        plt.subplot(435)
        # bins = np.linspace(500,550,100)
        _ = plt.hist(norm(MagX, MagY, MagZ),100, color='r', alpha=0.5)
        _ = plt.hist(norm(MagX, MagY, MagZ)[idx],100, color='k', alpha=0.5)
        _ = plt.hist(norm(magx,magy,magz)[idx],100, color='b', alpha=0.5)


        plt.subplot(436)
        _ = plt.hist(circ(cal_heading - ADCP.Heading.values), np.linspace(-20,20,100))


        def setSquare(gca):
            gca.set_box_aspect(1)

        plt.subplot(437)
        plt.axvline(0,color='k')
        plt.axhline(0,color='k')
        idz = MagZ > 0
        plt.scatter(MagX[idz], MagY[idz], 1, 'r')
        plt.scatter(MagX[idx&idz], MagY[idx&idz], 1, 'k')
        plt.scatter(magx[idx&idz],magy[idx&idz], 1, 'c')
        plt.title('MagZ > 0')
        plt.xlabel('MagX')
        plt.ylabel('MagY')
        plt.axis('square')

        plt.subplot(438)
        plt.axvline(0,color='k')
        plt.axhline(0,color='k')
        idz = MagZ < 0
        plt.scatter(MagX[idz], MagY[idz], 1, 'r')
        plt.scatter(MagX[idx&idz], MagY[idx&idz], 1, 'k')
        plt.scatter(magx[idx&idz],magy[idx&idz], 1, 'c')
        plt.title('MagZ < 0')
        plt.xlabel('MagX')
        plt.ylabel('MagY')
        plt.axis('square')

        plt.subplot(439)
        plt.axvline(0,color='k')
        plt.axhline(0,color='k')
        idz = ~np.isnan(MagZ)
        plt.scatter(np.sqrt(MagX[idz]**2 + MagY[idz]**2), MagZ[idz], 1, 'r')
        plt.scatter(np.sqrt(MagX[idx&idz]**2 + MagY[idx&idz]**2), MagZ[idx&idz], 1, 'k')
        plt.scatter(np.sqrt(magx[idx&idz]**2 + magy[idx&idz]**2),magz[idx&idz], 1, 'c')
        plt.xlabel('norm(MagX, MagY)')
        plt.ylabel('MagZ')
        plt.axis('square')
    
        plt.figure(figsize=(15,7))
        plt.scatter(cal_heading,norm(magx,magy,magz),2,'k')
        plt.scatter(cal_heading[idx],norm(magx,magy,magz)[idx],2,'c')
    
        plot_validation(cal_heading)
    
    return cal_heading



def correct_heading(ADCP, glider, options): # TODO Move to compass script
    if options["correctADCPHeading"]:
        if "Heading_old" in ADCP:
            ADCP["Heading"] = ("time", ADCP["Heading_old"].values)
            _log.info("Resetting to original heading")

        ADCP["Heading_old"] = ("time", ADCP["Heading"].values)
        ADCP["Heading"] = (
            _heading_correction(ADCP, glider, options) + ADCP["declination"]
        )
        plog("Corrected heading and accounted for declination")
    else:
        plog("Uncorrected heading and declination NOT added.")
    return ADCP




def get_declination(data, key): # TODO: move to compass file
    """
    Function retrieves declination data from NOAA for the average lon, lat and datetime of glider data.
    Requires an API key. Register at https://www.ngdc.noaa.gov/geomag/CalcSurvey.shtml

    :param data: pd.DataFrame of glider data including time, longitude and latitude
    :param key: API key for the NOAA geomag service
    :return: dataframe with added declination column
    """
    if "declination" in list(data):
        _log.info("declination data already present")
        return data
    time = data.time.mean()
    year = time.year
    month = time.month
    day = time.day
    lat = data.latitude.mean()
    lon = data.longitude.mean()
    url = (
        f"https://www.ngdc.noaa.gov/geomag-web/calculators/calculateDeclination?startYear={year}&startMonth={month}"
        f"&startDay={day}&lat1={lat}&lon1={lon}&key={key}&resultFormat=json"
    )
    result = json.load(request.urlopen(url))
    declination = result["result"][0]["declination"]
    _log.info(f"declination of {declination} added to data")
    data["declination"] = declination