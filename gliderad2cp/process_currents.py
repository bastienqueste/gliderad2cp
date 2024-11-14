import logging
import warnings
import xarray as xr
import numpy as np
import gsw
import pandas as pd
from scipy.interpolate import interp1d

warnings.filterwarnings(action="ignore", message="Mean of empty slice")
warnings.filterwarnings(action="ignore", message="invalid value encountered in divide")
warnings.filterwarnings(action="ignore", message="invalid value encountered in true_divide")
warnings.filterwarnings(action="ignore", message="Degrees of freedom <= 0 for slice.")

# def verify_depth_bias(out, yaxis, options, E="ADCP_E", N="ADCP_N"):
#     north = np.gradient(out["latitude"], axis=1) > 0
#     south = np.gradient(out["latitude"], axis=1) < 0
#     depths = np.linspace(0, np.max(yaxis) - 5, 20)
#     drange = np.mean(np.diff(depths)) / 2
#     bins = np.linspace(-1, 1, 100) * 0.5

#     variables = [E, N]
#     if options["debug_plots"]:
#         plt.figure(figsize=(20, 7))

#     SF = drange / 3

#     for idx, var in enumerate(variables):
#         if options["debug_plots"]:
#             plt.subplot(1, 3, idx + 1)
#             plt.axvline(0, color="k")

#         for d in depths:
#             depth = np.abs(out["Pressure"] - d) < drange

#             Nvals = out[var][(north & depth)]
#             Svals = out[var][(south & depth)]
#             N, _ = np.histogram(Nvals, bins=bins, density=True)
#             S, _ = np.histogram(Svals, bins=bins, density=True)

#             Nm = np.nanmean(Nvals)
#             Sm = np.nanmean(Svals)
#             Ns = np.nanstd(Nvals)
#             Ss = np.nanstd(Svals)
#             Nn = np.count_nonzero(np.isfinite(Nvals))
#             Sn = np.count_nonzero(np.isfinite(Svals))
#             Nse = Ns / np.sqrt(Nn)
#             Sse = Ss / np.sqrt(Sn)

#             N[N == 0] = np.nan
#             S[S == 0] = np.nan
#             if options["debug_plots"]:
#                 plt.fill_between(
#                     bins[1:], SF * N - float(d), -float(d), color="r", alpha=0.5
#                 )
#                 plt.fill_between(
#                     bins[1:], SF * S - float(d), -float(d), color="b", alpha=0.5
#                 )
#                 plt.plot(bins[1:], SF * N - float(d), "-r")
#                 plt.plot(bins[1:], SF * S - float(d), "-b")

#                 plt.plot(
#                     [Nm - 2 * Nse, Nm + 2 * Nse],
#                     np.array([1, 1]) * -float(d),
#                     "-k",
#                     linewidth=3,
#                     alpha=1,
#                 )
#                 plt.plot(
#                     [Sm - 2 * Sse, Sm + 2 * Sse],
#                     np.array([1, 1]) * -float(d),
#                     "-k",
#                     linewidth=3,
#                     alpha=1,
#                 )
#                 # plt.plot([Nm,Sm],np.array([1,1])*-float(d),'k',marker='.',linestyle='none')
#         if options["debug_plots"]:
#             plt.ylabel("Depth (m) / " + str(SF) + "*PDF")
#             plt.xlabel("Velocity")
#             plt.legend(("Zero", "Northward travel", "Southward travel"))
#             plt.title(var)
#     if options["debug_plots"]:
#         plt.subplot(133)
#         plt.axvline(0, color="k")

#     for d in depths:
#         depth = np.abs(out["Pressure"] - d) < drange

#         Nvals = np.sqrt(out["ADCP_E"] ** 2 + out["ADCP_N"] ** 2)[(north & depth)]
#         Svals = np.sqrt(out["ADCP_E"] ** 2 + out["ADCP_N"] ** 2)[(south & depth)]
#         N, _ = np.histogram(Nvals, bins=bins, density=True)
#         S, _ = np.histogram(Svals, bins=bins, density=True)

#         Nm = np.nanmean(Nvals)
#         Sm = np.nanmean(Svals)
#         Ns = np.nanstd(Nvals)
#         Ss = np.nanstd(Svals)
#         Nn = np.count_nonzero(np.isfinite(Nvals))
#         Sn = np.count_nonzero(np.isfinite(Svals))
#         Nse = Ns / np.sqrt(Nn)
#         Sse = Ss / np.sqrt(Sn)

#         N[N == 0] = np.nan
#         S[S == 0] = np.nan
#         if options["debug_plots"]:
#             plt.fill_between(
#                 bins[1:], SF * N - float(d), -float(d), color="r", alpha=0.5
#             )
#             plt.fill_between(
#                 bins[1:], SF * S - float(d), -float(d), color="b", alpha=0.5
#             )
#             plt.plot(bins[1:], SF * N - float(d), "-r")
#             plt.plot(bins[1:], SF * S - float(d), "-b")

#             plt.plot(
#                 [Nm - 2 * Nse, Nm + 2 * Nse],
#                 np.array([1, 1]) * -float(d),
#                 "-k",
#                 linewidth=3,
#                 alpha=1,
#             )
#             plt.plot(
#                 [Sm - 2 * Sse, Sm + 2 * Sse],
#                 np.array([1, 1]) * -float(d),
#                 "-k",
#                 linewidth=3,
#                 alpha=1,
#             )
#     if options["debug_plots"]:
#         plt.ylabel("Depth (m) / " + str(SF) + "*PDF")
#         plt.xlabel("Velocity")
#         plt.legend(("Zero", "Northward travel", "Southward travel"))
#         plt.title("MAG")
#         if options["plots_directory"]:
#             save_plot(options["plots_directory"], "verify_depth_bias")


# def calc_bias(out, yaxis, taxis, days, options):
#     """
#     Corrects gridded horizontal velocities for vertical shear bias.

#     :param out: xr.DataSet of gridded horizontal velocities
#     :param yaxis: depth bins of gridded velocity data
#     :param taxis: time of gridded velocity data
#     :param days: days to plot
#     :param options: options dict for processing
#     :return: xr.DataSet of gridded horizontal velocities corrected for shear bias
#     """

#     def get_bias(glider_speed, coeff):
#         r, c = np.shape(glider_speed)
#         bias = np.nancumsum(glider_speed, axis=0)
#         bias[~np.isfinite(glider_speed)] = np.nan
#         bias = bias - np.tile(np.nanmean(bias, axis=0), [r, 1])
#         return bias * coeff

#     def score(E, N):
#         def rmsd_h(x):
#             return np.sqrt(np.nanmean(x**2, axis=1))

#         def rmsd(x):
#             return np.sqrt(np.nanmean(x**2))

#         def y_weighting(x):
#             return x * 0 + 1

#         return rmsd((rmsd_h(E) + rmsd_h(N)) * y_weighting(yaxis)) * 1e6

#     def fn(coeff):
#         return score(
#             out["ADCP_E"] + get_bias(out["speed_e"], coeff),
#             out["ADCP_N"] + get_bias(out["speed_n"], coeff),
#         )

#     for _i in range(100):
#         R = fmin(
#             fn,
#             1,
#             maxiter=100,
#             ftol=0.00001,
#             disp=False,
#             retall=False,
#         )
#     _log.info(R)
#     coeff = R[0]

#     ADCP_E_old = out["ADCP_E"].copy()
#     ADCP_N_old = out["ADCP_N"].copy()

#     out["ADCP_E"] = out["ADCP_E"] + get_bias(out["speed_e"], coeff)
#     out["ADCP_N"] = out["ADCP_N"] + get_bias(out["speed_n"], coeff)
#     if options["debug_plots"]:
#         plt.figure(figsize=(15, 5))
#         plt.subplot(131)
#         plt.plot(np.nanvar(ADCP_E_old, axis=1), yaxis, "-r")
#         plt.plot(np.nanvar(out["ADCP_E"], axis=1), yaxis, "-g")
#         plt.gca().invert_yaxis()
#         plt.subplot(132)
#         plt.plot(np.nanvar(ADCP_N_old, axis=1), yaxis, "-r")
#         plt.plot(np.nanvar(out["ADCP_N"], axis=1), yaxis, "-g")
#         plt.gca().invert_yaxis()
#         plt.subplot(133)
#         plt.plot(
#             np.nanvar(np.sqrt(ADCP_E_old**2 + ADCP_N_old**2), axis=1), yaxis, "-r"
#         )
#         plt.plot(
#             np.nanvar(np.sqrt(out["ADCP_E"] ** 2 + out["ADCP_N"] ** 2), axis=1),
#             yaxis,
#             "-g",
#         )
#         plt.gca().invert_yaxis()
#         if options["plots_directory"]:
#             save_plot(options["plots_directory"], "calc_bias")
#     if options["debug_plots"]:
#         verify_depth_bias(out, yaxis, options, E="ADCP_E", N="ADCP_N")

#         plt.figure(figsize=(20, 20))
#         ## PLOT 1
#         plt.subplot(6, 1, 1)
#         plt.pcolormesh(taxis, yaxis, out["ADCP_E"], cmap="RdBu", shading="auto")
#         plt.clim(np.array([-1, 1]) * 0.5)
#         plt.colorbar()
#         [plt.axvline(x, color="k", alpha=0.3) for x in days]
#         plt.contour(
#             taxis,
#             yaxis,
#             out["salinity"],
#             np.linspace(35.5, 38.5, 6),
#             colors="k",
#             alpha=0.3,
#         )
#         plt.gca().invert_yaxis()
#         plt.xlabel("Yo number")
#         plt.xlabel("Depth")
#         plt.title("Eastward velocity (m.s-1)")

#         plt.subplot(6, 1, 2)
#         plt.pcolormesh(taxis, yaxis, out["ADCP_N"], cmap="RdBu", shading="auto")
#         plt.clim(np.array([-1, 1]) * 0.5)
#         plt.colorbar()
#         [plt.axvline(x, color="k", alpha=0.3) for x in days]
#         plt.contour(
#             taxis,
#             yaxis,
#             out["salinity"],
#             np.linspace(35.5, 38.5, 6),
#             colors="k",
#             alpha=0.3,
#         )
#         plt.gca().invert_yaxis()
#         plt.xlabel("Yo number")
#         plt.xlabel("Depth")
#         plt.title("Northward velocity (m.s-1)")
#         if options["plots_directory"]:
#             save_plot(options["plots_directory"], "calc_bias_2")
#     return out
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################

_log = logging.getLogger(__name__)

'''
Quality of life functions
'''
def plog(msg):
    '''
    Output information to the logger with timestamp.
    '''
    _log.info(str(dt.now().replace(microsecond=0)) + " : " + msg)
    return None

def grid2d(x, y, v, xi=1, yi=1, fn="median"):
    '''
    Quick data binning function relying on pandas.
    x,y,v are flat np.arrays of x, y coordinates and values at these points.
    xi and yi are either np.arrays of bins to be binned into, or the spacing used between min and max of x or y respectively.
    fn defines the function applied to all points from v that fall into the same bin.
    '''
    if np.size(xi) == 1:
        xi = np.arange(np.nanmin(x), np.nanmax(x) + xi, xi)
    if np.size(yi) == 1:
        yi = np.arange(np.nanmin(y), np.nanmax(y) + yi, yi)

    raw = pd.DataFrame({"x": x, "y": y, "v": v}).dropna()

    grid = np.full([np.size(yi), np.size(xi)], np.nan)

    raw["xbins"], xbin_iter = pd.cut(raw.x, xi, retbins=True, labels=False, right=False)
    raw["ybins"], ybin_iter = pd.cut(raw.y, yi, retbins=True, labels=False, right=False)

    _tmp = raw.groupby(["xbins", "ybins"])["v"].agg(fn)
    grid[
        _tmp.index.get_level_values(1).astype(int),
        _tmp.index.get_level_values(0).astype(int),
    ] = _tmp.values

    XI, YI = np.meshgrid(xi, yi, indexing="ij")
    return grid, XI.T, YI.T

def interp(x, y, xi):
    '''
    Interpolation shortcut that removes NaN data and ignores boundary errors.
    '''
    _gg = np.isfinite(x + y)
    return interp1d(x[_gg], y[_gg], bounds_error=False, fill_value=np.nan)(xi)

'''
Processing functions
'''
def get_DAC(ADCP, gps_predive, gps_postdive):
    """
    Calculate dive-averaged currents using the ADCP as a DVL.

    ...

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
    ## Process time data
    time = ADCP.time.values.astype("float") / 1e9
    deltat = np.append(np.diff(time), np.nan)
    f = np.nanmedian(deltat)
    deltat[deltat > f*10] = np.NaN
    deltat = np.nan_to_num(deltat)
    
    ## Create interpolannt to query DVL travel at time t.
    
    #Using the 3-beam configuration:
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
    
    # # Using the 4-beam configuration:
    # spd_thr_water = np.sqrt(ADCP["X4"] ** 2 + ADCP["Y4"] ** 2 + ADCP["ZZ4"] ** 2).mean("bin")
    # hspd = np.cos(np.deg2rad(np.abs(ADCP["Pitch"]) + 5)) * spd_thr_water
    # angle = np.deg2rad(90 - ADCP["Heading"])
    # e = hspd * np.cos(angle)
    # n = hspd * np.sin(angle)
    # e = np.nan_to_num(e) * deltat
    # n = np.nan_to_num(n) * deltat
    # travel_e = interp1d(
    #         time,
    #         np.cumsum(np.nan_to_num(e) * deltat),
    #         bounds_error=False, fill_value=np.nan
    #         )
    # travel_n = interp1d(
    #         time,
    #         np.cumsum(np.nan_to_num(e) * deltat),
    #         bounds_error=False, fill_value=np.nan
    #         )
    
    r,c = np.shape(gps_predive)
    
    dive_time = np.full(r, np.NaN) # Duration in seconds
    dive_duration = np.full(r, np.NaN) # Duration in seconds
    
    dac = np.full([r,2], np.NaN) # dive-averaged currents in m.s-1
    dxy_gps =  np.full([r,2], np.NaN) # x/y displacement in m, over earth, for the dive calculated from GPS coordinates
    dxy_dvl =  np.full([r,2], np.NaN) # x/y displacement in m, through water, for the dive calculated from ADCP as DVL
    
    pre_t  = gps_predive[:,0].astype('float') / 1e9 # s
    post_t = gps_postdive[:,0].astype('float') / 1e9 # s
    
    def lon2m(x, y):
        return gsw.distance([x-0.0005, x + 0.0005], [y, y])*1000
    def lat2m(x, y):
        return gsw.distance([x, x], [y-0.0005, y+0.0005])*1000
    
    for idx in range(len(gps_predive)):
        dt = post_t - pre_t[idx] # s
        dt[dt < 0] = np.NaN
        
        if all(np.isnan(dt)):
            print(f'Could not find a surfacing after dive starting: {gps_predive[idx,0]}.')
            continue
        
        idy = np.nanargmin(dt)
        dive_duration[idx] = dt[idy] # s
        dive_time[idx] = (post_t[idy]+pre_t[idx])/2
        
        dxy_dvl[idx,0] = travel_e(post_t[idy]) - travel_e(pre_t[idx])
        dxy_dvl[idx,1] = travel_n(post_t[idy]) - travel_n(pre_t[idx])
        
        dxy_gps[idx,0] = (gps_postdive[idy,1] - gps_predive[idx,1]) * lon2m(gps_predive[idx,1],gps_predive[idx,2])
        dxy_gps[idx,1] = (gps_postdive[idy,2] - gps_predive[idx,2]) * lat2m(gps_predive[idx,1],gps_predive[idx,2])
    
    dac = (dxy_gps - dxy_dvl) / np.vstack([dive_duration,dive_duration]).T
    dive_time = pd.to_datetime(dive_time*1e9).values
    
    return dac, dive_time, dive_duration, dxy_gps, dxy_dvl

def grid_shear(ADCP, xi, yi):
    """
    Grid shear according to specifications.

    ...

    Inputs
    ----------
    ADCP : xarray.Dataframe
        Output from the gliderAD2CP.process_shear() function.
    xi : numpy.array
        Array indicating bins in x-axis.
    yi : numpy.array
        Array indicating bins in y-axis.
         
    Outputs
    -------
    currents : xr.Dataset
        Dataset containing gridded shear and various statistical metrics, as well as time spent per bin by the glider.    
        
    """
    
    print(f'Creating output structure and time-based metrics...')
    
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
    grid = lambda var_name, fn_name : grid2d(x, y, ADCP[var_name].values.flatten(), xi=xi, yi=yi, fn=fn_name)[0]
    count_valid = lambda x : np.count_nonzero(np.isfinite(x))
    # Grid shear data and metrics
    print(f'Gridding eastward shear metrics...')
    currents['shear_E_mean'] = (('depth', 'profile_index'),   grid('Sh_E', np.nanmean))
    currents['shear_E_median'] = (('depth', 'profile_index'), grid('Sh_E', np.nanmedian))
    currents['shear_E_stddev'] = (('depth', 'profile_index'), grid('Sh_E', np.nanstd))
    currents['shear_E_count'] = (('depth', 'profile_index'),  grid('Sh_E', count_valid))

    print(f'Gridding northward shear metrics...')
    currents['shear_N_mean'] = (('depth', 'profile_index'),   grid('Sh_N', np.nanmean))
    currents['shear_N_median'] = (('depth', 'profile_index'), grid('Sh_N', np.nanmedian))
    currents['shear_N_stddev'] = (('depth', 'profile_index'), grid('Sh_N', np.nanstd))
    currents['shear_N_count'] = (('depth', 'profile_index'),  grid('Sh_N', count_valid))

    print(f'Calculating standard error as percentage of shear...')
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

def grid_velocity(currents, method='integrate'):
    """
    Assemble shear measurements to reconstruct velocity profiles.
    Currently only able to integrate in the vertical, but I have aspirations to add the least-squared method too.
    Untested and likely not compatible with irregular grids.

    ...

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
    def integrate():
        def integrate_calc(Sh):
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
        currents['velocity_E_no_reference'] = (('depth', 'profile_index'), integrate_calc(currents.shear_E_median.values))
        currents['velocity_N_no_reference'] = (('depth', 'profile_index'), integrate_calc(currents.shear_N_median.values))
     
    def lsq():
        # Space for Martin Visbeck's least squared method
        return False
    
    if method == 'integrate':
        fn = integrate
    elif method == 'lsq':
        print('WARNING: least squared method not yet coded up. Oops. Any volunteers?')
        fn = lsq
        
    fn()
    
    return currents

def reference_velocity(currents, DAC, dive_time):
    """
    Reference vertical velocity profiles to dive-averaged currents, paying attenting to the time spent at each depth.
    This ensures that dives with irregular vertical speeds or periods of loitering at depth are still reproduced correctly.
    As a result, one cannot expect the *depth*-averaged velocity of a dive to match the *dive*-averaged currents.

    ...

    Inputs
    ----------
    currents : xr.Dataset
        Dataset containing gridded shear produced by the process_currents.grid_shear() function.
    DAC
        Mx2 np.array containing eastward DAC (column 0) and northward DAC (column 1), for M dives identified.
    dive_time
        Mx1 np.array containing mean time for the dive in np.datetime64[ns] format.
        
         
    Outputs
    -------
    currents : xr.Dataset
        Dataset containing gridded shear, statistical metrics, time spent per bin by the glider, and unreferenced velocity profiles.
        
    """
    
    # DAC = \int_{0}^{z} \frac{velocity * time spent in bin}{time spent in profile} + referencing
    # DAC - \int_{0}^{z} \frac{velocity * time spent in bin}{time spent in profile} = referencing
    time_in_bin, depth = xr.broadcast(currents.time_in_bin, currents.depth)
    time_in_bin = time_in_bin.values
    depth = depth.values
    time_in_bin_mean = np.nanmean(time_in_bin, axis=0)
    
    
    for DAC_idx, direction in enumerate(['E','N']):
        currents[f'DAC_{direction}'] = ('profile_index', interp(dive_time.astype('float'), DAC[:,DAC_idx], currents.time.values.astype('float')))
        
        velocity = currents[f'velocity_{direction}_no_reference'].values * time_in_bin / np.broadcast_to(time_in_bin_mean[np.newaxis,:],time_in_bin.shape)
                
        reference = currents[f'DAC_{direction}'].values - np.nanmean(velocity, axis=0)
        
        currents[f'reference_{direction}_from_DAC'] = ('profile_index', reference)
        
        currents[f'velocity_{direction}_DAC_reference'] = currents[f'velocity_{direction}_no_reference'] + currents[f'reference_{direction}_from_DAC']
        
    return currents

'''
Main
'''
def process(ADCP, gps_predive, gps_postdive, xi=1, yi=None, correct_shear_bias=True):
    """
    Calculate depth-resolved currents from shear data calculated by the gliderAD2CP.process_shear() function.
    This functionality can handle loiter dives and irregular profiles but it requires that your ADCP is 
    always on so that DAC is representative of the baroclinic profile collected by the ADCP.
    ...

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
    xi : int or numpy.array
        Defines the x-axis of the output product. 
        If int, the x-axis will span from the first to last profile number in xi-increments. Use xi=1 per profile or xi=2 per dive.
        If numpy.array, this defines the x-axis the data is gridded on.
    yi : int or numpy.array or None
        Defines the y-axis of the output product. 
        If int, the y-axis will span from the minimum to maximum depth in yi-increments.
        If numpy.array, this defines the y-axis the data is gridded on.
        If None, the ADCP sampling cell size is provided as yi.
    correct_shear_bias : Bool, default=True
        Correct for shear bias using an elaboration on Todd et al. (2017)'s approach.
        
        
    Outputs
    -------
    currents : xr.Dataset
        Dataset containing gridded shear, statistical metrics, time spent per bin by the glider, unreferenced velocity profiles and DAC-referenced velocities.


    .process() runs the following methods in this order
    -------
    1. get_DAC
    2. grid_shear
    3. grid_velocity
    4. reference_velocity
    
    (Optional steps:) - Not coded up yet
    (5. assess_shear_bias) - Not coded up yet
    (6. correct_shear_bias) - Not coded up yet.
    
    
    Methods
    -------
    grid2d
        Quick data binning function relying on pandas.
    interp
        Interpolation shortcut that removes NaN data and ignores boundary errors.
    get_DAC
        Calculate dive-averaged currents using the ADCP as a DVL.
    grid_shear
        Grid shear according to specifications.
    grid_velocity
        Assemble shear measurements to reconstruct velocity profiles.
    reference_velocity
        Reference vertical velocity profiles to dive-averaged currents, paying attenting to the time spent at each depth.
    process
        This function.    
    """
    
    if yi is None:
        yi = ADCP.attrs['avg_cellSize'].astype('float')
        plog(f'Setting yi to {yi}, based on ADCP sampling cell size.')
    if np.size(xi) == 1:
        xi = np.arange(np.nanmin(ADCP.profile_number.values), np.nanmax(ADCP.profile_number.values) + xi, xi)
    if np.size(yi) == 1:
        yi = np.arange(0, np.nanmax(np.nanmax(ADCP.bin_depth)) + yi, yi)
        
    DAC, dive_time, _, _, _ = get_DAC(ADCP, gps_predive, gps_postdive)
    
    currents = grid_shear(ADCP, xi, yi)
    currents = grid_velocity(currents)
    currents = reference_velocity(currents,DAC,dive_time)
    
    return currents
