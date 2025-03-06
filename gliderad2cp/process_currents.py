"""
Calculate depth-resolved currents from shear data calculated by the gliderAD2CP.process_shear.process() function.
This functionality can handle loiter dives and irregular profiles but it requires that your ADCP is 
always on so that DAC is representative of the baroclinic profile collected by the ADCP.


gliderad2cp.process_currents
------------------------------
process
    The main function which converts shear data from gliderAD2CP.process_shear.process() to referenced velocities.
get_DAC
    Calculate dive-averaged currents using the ADCP as a DVL.
_grid_shear
    Grid shear according to specifications.
_grid_velocity
    Assemble shear measurements to reconstruct velocity profiles.
_reference_velocity
    Reference vertical velocity profiles to dive-averaged currents, paying attenting to the time spent at each depth.
    
Notes
-------
.process() runs the following functions in this order

1. get_DAC - Calculate dive-averaged currents using the ADCP as a DVL.
2. _grid_shear - Grid shear according to specifications.
3. _grid_velocity - Assemble shear measurements to reconstruct velocity profiles.
4. _reference_velocity - Reference vertical velocity profiles to dive-averaged currents, paying attenting to the time spent at each depth.

(Optional steps:) - Not coded up yet

5. assess_shear_bias - Not coded up yet
6. correct_shear_bias - Not coded up yet.

"""

import numpy as np
import pandas as pd
import warnings
import xarray as xr
import gsw
from scipy.interpolate import interp1d
from .tools import plog, grid2d, interp, get_options

warnings.filterwarnings(action="ignore", message="Mean of empty slice")
warnings.filterwarnings(action="ignore", message="invalid value encountered in divide")
warnings.filterwarnings(action="ignore", message="invalid value encountered in true_divide")
warnings.filterwarnings(action="ignore", message="Degrees of freedom <= 0 for slice.")
warnings.filterwarnings(action="ignore", message="All-NaN slice encountered")
warnings.simplefilter(action='ignore', category=FutureWarning)


"""
Processing functions
"""
def get_DAC(ADCP, gps_predive, gps_postdive):
    """
    Calculate dive-averaged currents using the ADCP as a DVL.


    Inputs
    ----------
    ADCP : xarray.Dataframe
        Output from the gliderAD2CP.process_shear() function.
    gps_predive : numpy.array
        Mx3 numpy array where M in the number of dives performed by the glider.
        Columns are np.datetime64[ns], longitude and latitude.
    gps_postdive : numpy.array
        Mx3 numpy array where M in the number of dives performed by the glider. 
        Columns are np.datetime64[ns], longitude and latitude.
         
    Outputs
    -------
    DAC : dictionary containing : 
        dac
            Mx2 np.array containing eastward DAC (column 0) and northward DAC (column 1), for M dives identified.
        dive_time
            Mx1 np.array containing mean time for the dive in np.datetime64[ns] format.
        dive_duration
            Mx1 np.array containing mean duration for the dive in seconds.
        dxy_gps, dxy_dvl
            Mx2 np.arrays containing eastward (col 0) and northward (col 1) displacement in m,
            calculated from GPS and from ADCP as a DVL respectively.
    
    """
    plog('Calculating dive-averaged currents from dive and surfacing GPS coordinates.')
    ## Process time data
    time = ADCP.time.values.astype("float") / 1e9
    deltat = np.append(np.diff(time), np.nan)
    f = np.nanmedian(deltat)
    deltat[deltat > f*10] = np.nan
    deltat = np.nan_to_num(deltat)
    
    ## Create interpolannt to query DVL travel at time t.
    travel_e = interp1d(
            time,
            np.cumsum(np.nan_to_num(-(ADCP["E"].mean(dim="gridded_bin") * deltat).values)),
            bounds_error=False, fill_value=np.nan
            )
    travel_n = interp1d(
            time,
            np.cumsum(np.nan_to_num(-(ADCP["N"].mean(dim="gridded_bin") * deltat).values)),
            bounds_error=False, fill_value=np.nan
            )
    
    ## Preallocate arrays
    r,c = np.shape(gps_predive)
    
    dive_time = np.full(r, np.nan) # Duration in seconds
    dive_duration = np.full(r, np.nan) # Duration in seconds
    
    dac = np.full([r,2], np.nan) # dive-averaged currents in m.s-1
    dxy_gps =  np.full([r,2], np.nan) # x/y displacement in m, over earth, for the dive calculated from GPS coordinates
    dxy_dvl =  np.full([r,2], np.nan) # x/y displacement in m, through water, for the dive calculated from ADCP as DVL
    
    ## Useful tidbits
    pre_t = gps_predive[:, 0].astype('datetime64[s]').astype(np.float64) # s
    post_t = gps_postdive[:, 0].astype('datetime64[s]').astype(np.float64) # s

    lon2m = lambda x, y : gsw.distance([x-0.0005, x + 0.0005], [y, y])*1000
    lat2m = lambda x, y : gsw.distance([x, x], [y-0.0005, y+0.0005])*1000
    
    ## Calculate displacements
    for idx in range(len(gps_predive)):
        dt = post_t - pre_t[idx] # s
        dt[dt < 0] = np.nan
        
        if all(np.isnan(dt)):
            plog(f'    Could not find a surfacing after dive starting: {gps_predive[idx,0]}.')
            continue
        
        idy = np.nanargmin(dt)
        dive_duration[idx] = dt[idy] # s
        dive_time[idx] = (post_t[idy]+pre_t[idx])/2
        
        dxy_dvl[idx,0] = travel_e(post_t[idy]) - travel_e(pre_t[idx])
        dxy_dvl[idx,1] = travel_n(post_t[idy]) - travel_n(pre_t[idx])
        
        dxy_gps[idx,0] = (gps_postdive[idy,1] - gps_predive[idx,1]) * lon2m(gps_predive[idx,1],gps_predive[idx,2])
        dxy_gps[idx,1] = (gps_postdive[idy,2] - gps_predive[idx,2]) * lat2m(gps_predive[idx,1],gps_predive[idx,2])
    
    ## Calculate dive-averaged currents
    dac = (dxy_gps - dxy_dvl) / np.vstack([dive_duration,dive_duration]).T
    dive_time = pd.to_datetime(dive_time*1e9).values
    
    ## Create output dictionary
    DAC = {'dac': dac, 'dive_time': dive_time, 'dive_duration': dive_duration, 'xy_displacement_over_land': dxy_gps ,'xy_displacement_through_water': dxy_dvl}
    
    return DAC

def _grid_shear(ADCP, options, xi, yi):
    """
    Grid shear according to specifications.


    Inputs
    ----------
    ADCP : xarray.Dataframe
        Output from the gliderAD2CP.process_shear() function.
    options : dict
        Set of options for gliderAD2CP, created by the gliderad2cp.tools.get_options() function.
    xi : numpy.array
        Array of bins in x-axis.
    yi : numpy.array
        Array of bins in y-axis.
         
    Outputs
    -------
    currents : xr.Dataset
        Dataset containing gridded shear and various statistical metrics, as well as time spent per bin by the glider.    
        
    """
    
    plog(f'Gridding shear data.')
    plog(f'    Creating output structure and time-based metrics.')
    
    # Create structure
    currents = xr.Dataset(coords={"depth": yi, "profile_index": xi})
    
    # Get coords for gridding function below in the right dimensions
    x, y = xr.broadcast(ADCP.profile_number, ADCP.bin_depth)
    t, _ = xr.broadcast(ADCP.time, ADCP.bin_depth)
    x = x.values.flatten()
    y = y.values.flatten()
    t = t.values.astype('float').flatten()

    # Create time coordinates
    time = np.nanmedian(grid2d(x, y, t, xi=xi, yi=yi, fn=np.nanmedian)[0], axis=0)
    currents = currents.assign_coords(time=("profile_index", pd.to_datetime(time)))
    
    # Calculate time in bin
    # This metric is important when referencing to DAC later.
    # We weight the baroclinic velocity by time spent in each bin when referencing to DAC.
    # This is relevant on loiter dives for example, where we may spend 24 hours loitering 
    # and as a result the DAC is heavily biased to the velocity at the loiter depth.
    # Calculating time in each bin is very easy (differentiate sampling frequency and sum in bins)
    # ... this works great to calculate loiters etc.... 
    # ... BUT we have very long periods at the surface which are NOT included in the DAC.
    # Here we set durations to 0 when dP/dt is under a threshold near the surface.
    # It's not perfect, but it's fairly cross-compatible assumption.    
    duration = np.gradient(ADCP.time.values.astype('float'))/1e9 # Time spent at that depth
    # Set surface values to 0. Imperfect because pressure sensor drift but captures a good amout. We lose some good data but very minimal.
    duration[ADCP.Depth.values < 1] = 0 
    # Same with a looser depth threshold but only if speeds are less than 4cm/s.
    duration[
        (ADCP.Depth.values < 3) & (np.abs(np.gradient(ADCP.Depth.values, ADCP.time.values.astype('float')/1e9)) < 0.04)
         ] = 0
    currents['time_in_bin'] = (('depth', 'profile_index'),
                               grid2d(ADCP.profile_number.values, ADCP.Depth.values, 
                                 duration,
                                 xi=xi, yi=yi, fn=np.nansum)[0]
                              )

    # Gridding function for shear
    plog(f'    Gridding flight metrics.')
    grid = lambda var_name, fn_name : grid2d(x, y, ADCP[var_name].values.flatten(), xi=xi, yi=yi, fn=fn_name)[0]
    count_valid = lambda x : np.count_nonzero(np.isfinite(x))
    cos_mean = lambda head : np.nanmean(np.cos(np.deg2rad(head))) # N
    sin_mean = lambda head : np.nanmean(np.sin(np.deg2rad(head))) # E
    
    currents['heading_N'] = (('depth', 'profile_index'),
                       grid2d(ADCP.profile_number.values, ADCP.Depth.values, 
                         ADCP.Heading.values + options['heading_offset'],
                         xi=xi, yi=yi, fn=cos_mean)[0]
                      )
    currents['heading_E'] = (('depth', 'profile_index'),
                           grid2d(ADCP.profile_number.values, ADCP.Depth.values, 
                             ADCP.Heading.values + options['heading_offset'],
                             xi=xi, yi=yi, fn=sin_mean)[0]
                          )
    currents['speed_through_water'] = (('depth', 'profile_index'),
                           grid2d(ADCP.profile_number.values, ADCP.Depth.values, 
                             np.sqrt(ADCP.X.mean(dim='gridded_bin')**2 +ADCP.Y.mean(dim='gridded_bin')**2 +ADCP.Z.mean(dim='gridded_bin')**2),
                             xi=xi, yi=yi, fn=np.nanmean)[0]
                          )
    
    # Grid shear data and metrics
    plog(f'    Gridding eastward shear metrics.')
    currents['shear_E_mean'] = (('depth', 'profile_index'),   grid('Sh_E', np.nanmean))
    currents['shear_E_median'] = (('depth', 'profile_index'), grid('Sh_E', np.nanmedian))
    currents['shear_E_stddev'] = (('depth', 'profile_index'), grid('Sh_E', np.nanstd))
    currents['shear_E_count'] = (('depth', 'profile_index'),  grid('Sh_E', count_valid))

    plog(f'    Gridding northward shear metrics.')
    currents['shear_N_mean'] = (('depth', 'profile_index'),   grid('Sh_N', np.nanmean))
    currents['shear_N_median'] = (('depth', 'profile_index'), grid('Sh_N', np.nanmedian))
    currents['shear_N_stddev'] = (('depth', 'profile_index'), grid('Sh_N', np.nanstd))
    currents['shear_N_count'] = (('depth', 'profile_index'),  grid('Sh_N', count_valid))

    plog(f'    Calculating standard error as percentage of shear.')
    currents['shear_E_pct_error'] = np.abs(
        100 * 
        currents.shear_E_stddev / 
        np.sqrt(currents.shear_E_count)
        / currents.shear_E_mean
    )

    currents['shear_N_pct_error'] = np.abs(
        100 * 
        currents.shear_N_stddev / 
        np.sqrt(currents.shear_N_count)
        / currents.shear_N_mean
    )
    
    return currents

def _grid_velocity(currents, method='integrate'):
    """
    Assemble shear measurements to reconstruct velocity profiles.
    Currently only able to integrate in the vertical, but I have aspirations to add the least-squared method too.
    Untested and likely not compatible with irregular grids.


    Inputs
    ----------
    currents : xr.Dataset
        Dataset containing gridded shear produced by the process_currents.grid_shear() function.
    method : str, default = 'integrate'
        Which method to use to reconstruct velocity profiles. Currently only "integrate" is available.
        "lsq" coming one day. maybe.
        
         
    Outputs
    -------
    currents : xr.Dataset
        Dataset containing gridded shear, statistical metrics, time spent per bin by the glider, and unreferenced velocity profiles.
        
    """
    def __integrate():
        def __integrate_calc(Sh):
            # Integrate shear vertically
            _bd = ~np.isfinite(Sh)  # Preserve what are originally NaN values to recover later as need conversion to 0 for cumsum-
            Sh = np.nan_to_num(Sh)  # Replace NaNs with 0 for cumsum
            y_res = np.gradient(currents.depth.values)
            y_res = np.broadcast_to(y_res[:,np.newaxis],Sh.shape)
            V = np.cumsum(Sh * y_res, axis=0)  # Cumulative sum of shear to recover velocity profile
            V[_bd] = np.nan  # Return NaNs to their rightful place.
            V = V - np.tile(
                np.nanmean(V, axis=0), (np.shape(V)[0], 1)
            )  # Make mean of baroclinic profiles equal to 0
            return V        
        currents['velocity_E_no_reference'] = (('depth', 'profile_index'), __integrate_calc(currents.shear_E_median.values))
        currents['velocity_N_no_reference'] = (('depth', 'profile_index'), __integrate_calc(currents.shear_N_median.values))
    def __lsq():
        # Space for Martin Visbeck's least squared method
        return False
    
    plog(f'Calculating unreferenced velocity profiles from shear.')
    
    if method == 'integrate':
        fn = __integrate
    elif method == 'lsq':
        plog('WARNING: least squared method not yet coded up. Oops. Any volunteers?')
        fn = __lsq
        
    fn()
    
    return currents

def _reference_velocity(currents, DAC):
    """
    Reference vertical velocity profiles to dive-averaged currents, paying attention to the time spent at each depth.
    This ensures that dives with irregular vertical speeds or periods of loitering at depth are still reproduced correctly.
    As a result, one cannot expect the *depth*-averaged velocity of a dive to match the *dive*-averaged currents.


    Inputs
    ----------
    currents : xr.Dataset
        Dataset containing gridded shear produced by the process_currents.grid_shear() function.
    DAC : dictionary containing fields :
        dac : 
            Mx2 np.array containing eastward DAC (column 0) and northward DAC (column 1), for M dives identified.
        dive_time : 
            Mx1 np.array containing mean time for the dive in np.datetime64[ns] format.
        
         
    Outputs
    -------
    currents : xr.Dataset
        Dataset containing gridded shear, statistical metrics, time spent per bin by the glider, and unreferenced velocity profiles.
        
    """
    
    plog('Referencing velocity profiles to dive-averaged currents using baroclinic velocity weighting.')
    # DAC = \int_{0}^{z} \frac{velocity * time spent in bin}{time spent in profile} + referencing
    # DAC - \int_{0}^{z} \frac{velocity * time spent in bin}{time spent in profile} = referencing
    
    time_in_bin, depth = xr.broadcast(currents.time_in_bin, currents.depth)
    time_in_bin = time_in_bin.values
    depth = depth.values
    time_in_bin_mean = np.nanmean(time_in_bin, axis=0)
    
    
    for DAC_idx, direction in enumerate(['E','N']):
        currents[f'DAC_{direction}'] = ('profile_index', interp(DAC['dive_time'].astype('float'), DAC['dac'][:,DAC_idx], currents.time.values.astype('float')))
        currents[f'displacement_though_water_{direction}'] = ('profile_index', interp(DAC['dive_time'].astype('float'), DAC['xy_displacement_through_water'][:,DAC_idx], currents.time.values.astype('float')))

        
        velocity = currents[f'velocity_{direction}_no_reference'].values * time_in_bin / np.broadcast_to(time_in_bin_mean[np.newaxis,:],time_in_bin.shape)
        reference = currents[f'DAC_{direction}'].values - np.nanmean(velocity, axis=0)

        currents[f'reference_{direction}_from_DAC'] = ('profile_index', reference)

        currents[f'velocity_{direction}_DAC_reference'] = currents[f'velocity_{direction}_no_reference'] + currents[f'reference_{direction}_from_DAC']

    return currents


"""
Main
"""
def process(ADCP, gps_predive, gps_postdive, options=None):
    """
    Calculate depth-resolved currents from shear data calculated by the gliderAD2CP.process_shear.process() function.
    This functionality can handle loiter dives and irregular profiles but it requires that your ADCP is 
    always on so that DAC is representative of the baroclinic profile collected by the ADCP.
    

    Inputs
    ----------
    ADCP : xarray.Dataframe
        Output from the gliderAD2CP.process_shear() function.
    gps_predive : numpy.array
        Mx3 numpy array where M in the number of dives performed by the glider.
        Columns are np.datetime64[ns], longitude and latitude.
    gps_postdive : numpy.array
        Mx3 numpy array where M in the number of dives performed by the glider. 
        Columns are np.datetime64[ns], longitude and latitude.
    options : dict
        Set of options for gliderAD2CP, created by the gliderad2cp.tools.get_options() function.
    
        
    Outputs
    -------
    currents : xr.Dataset
        Dataset containing gridded shear, statistical metrics, time spent per bin by the glider, unreferenced velocity profiles and DAC-referenced velocities.

    Notes
    ------
    .process() runs the following functions in this order

    1. get_DAC - Calculate dive-averaged currents using the ADCP as a DVL.
    2. _grid_shear - Grid shear according to specifications.
    3. _grid_velocity - Assemble shear measurements to reconstruct velocity profiles.
    4. _reference_velocity - Reference vertical velocity profiles to dive-averaged currents, paying attention to the time spent at each depth.
    
    """
    # Load default options if not present.
    if not options:
        options = get_options(verbose=False)
        plog('Using default set of options. See gliderad2cp.tools.get_options() for settings.')
    
    xi = options['xaxis']
    yi = options['yaxis']
    
    if yi is None:
        yi = ADCP.attrs['avg_cellSize'].astype('float')
        plog(f'Setting y-axis resolution to {yi}, based on ADCP sampling cell size.')
        
    if np.size(xi) == 1:
        xi = np.arange(np.nanmin(ADCP.profile_number.values), np.nanmax(ADCP.profile_number.values) + xi, xi)
    if np.size(yi) == 1:
        yi = np.arange(0, np.nanmax(np.nanmax(ADCP.bin_depth)) + yi, yi)
        
    DAC = get_DAC(ADCP, gps_predive, gps_postdive)
    
    currents = _grid_shear(ADCP, options, xi, yi)
    currents = _grid_velocity(currents, method=options['shear_to_velocity_method'])
    currents = _reference_velocity(currents,DAC)
    
    return currents, DAC