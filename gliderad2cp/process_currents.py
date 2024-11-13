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


# def reference_shear(ADCP, glider, dE, dN, dT, xaxis, yaxis, taxis, options):
#     """
#     Reference the estimates of vertical shear of horizontal velocity using a per-profile average velocity

#     :param ADCP: xr.DataSet of ADCP data
#     :param options: options dictionary
#     :param dE: eastward surface drift velocity
#     :param dN: northward surface drift velocity
#     :param dT: timestamps surface drift velocity
#     :param xaxis: profile numbers of gridded velocity data
#     :param yaxis: depth bins of gridded velocity data
#     :param taxis: time of gridded velocity data
#     :param options: options dict for processing
#     :return:  xr.DataSet of referenced gridded N and E velocities
#     """
#     out = {}

#     var = ["E", "N"]
#     if options["debug_plots"]:
#         plt.figure(figsize=(20, 20))

#     days = np.unique(glider.time.round("D"))
#     for pstep in range(len(var)):
#         letter = var[pstep]
#         # Grid shear to average out sensor + zooplankton noise
#         Sh, XI, YI = grid2d(
#             np.tile(ADCP.profile_number.values, (len(ADCP.gridded_bin), 1)).T.flatten(),
#             ADCP.bin_depth.values.flatten(),
#             ADCP["Sh_" + letter].values.flatten(),
#             xi=xaxis,
#             yi=yaxis,
#             fn="mean",
#         )

#         # Integrate shear vertically
#         _bd = ~np.isfinite(
#             Sh
#         )  # Preserve what are originally NaN values to recover later as need conversion to 0 for cumsum-
#         Sh = np.nan_to_num(Sh)  # Replace NaNs with 0 for cumsum
#         V = (
#             np.cumsum(Sh, axis=0) * y_res
#         )  # Cumulative sum of shear to recover velocity profile
#         V[_bd] = np.nan  # Return NaNs to their rightful place.
#         V = V - np.tile(
#             np.nanmean(V, axis=0), (np.shape(V)[0], 1)
#         )  # Make mean of baroclinic profiles equal to 0

#         # Grid DAC
#         DAC, XI, YI = grid2d(
#             glider.profile_number.values,
#             glider.pressure.values,
#             glider["DAC_" + letter].values,
#             xi=xaxis,
#             yi=yaxis,
#             fn="mean",
#         )

#         # Grid vertical speed
#         dPdz, XI, YI = grid2d(
#             glider.profile_number.values,
#             glider.pressure.values,
#             glider["speed_vert"].values,
#             xi=xaxis,
#             yi=yaxis,
#             fn="mean",
#         )

#         # Grid salinity
#         SA, XI, YI = grid2d(
#             glider.profile_number.values,
#             glider.pressure.values,
#             glider.salinity.values,
#             xi=xaxis,
#             yi=yaxis,
#             fn="median",
#         )

#         # Seconds spent in each depth bin, to weight referencing
#         SpB = y_res / dPdz
#         SpB[np.isinf(SpB)] = 0
#         strictness = 1
#         SpB_std = np.nanstd(SpB.flatten())
#         SpB[np.abs(SpB) > (strictness * SpB_std)] = strictness * SpB_std

#         # Baroclinic velocity, weighted by depth residence time, should be equal to DAC
#         # So the reference to add to a baroclinic profile of mean = 0 is the DAC - the weighted baroclinic velocity.
#         Ref = np.nanmean(DAC, axis=0) - np.nansum(V * SpB, axis=0) / np.nansum(
#             SpB, axis=0
#         )

#         # Now we reference the velocity
#         V = V + np.tile(Ref, (np.shape(V)[0], 1))
   



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
    
    ## Process time data
    time = ADCP.time.values.astype("float") * 1e-9
    deltat = np.append(np.diff(ADCP["time"].values.astype("float") / 1e9), np.nan)
    f = np.nanmedian(deltat)
    deltat[deltat > f*10] = np.NaN
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
    
    r,c = np.shape(gps_predive)
    
    dive_time = np.full(r, np.NaN) # Duration in seconds
    dive_duration = np.full(r, np.NaN) # Duration in seconds
    
    dac = np.full([r,2], np.NaN) # dive-averaged currents in m.s-1
    dxy_gps =  np.full([r,2], np.NaN) # x/y displacement in m, over earth, for the dive calculated from GPS coordinates
    dxy_dvl =  np.full([r,2], np.NaN) # x/y displacement in m, through water, for the dive calculated from ADCP as DVL
    
    pre_t  = gps_predive[:,0].astype('float') * 1e-9
    post_t = gps_postdive[:,0].astype('float') * 1e-9
    
    
    def lon2m(x, y):
        return gsw.distance([x-0.0005, x + 0.0005], [y, y])*1000
    def lat2m(x, y):
        return gsw.distance([x, x], [y-0.0005, y+0.0005])*1000
    
    for idx in range(len(gps_predive)):
        dt = post_t - pre_t[idx]
        dt[dt < 0] = np.NaN
        
        if all(np.isnan(dt)):
            print(f'Could not find a surfacing after dive starting: {gps_predive[idx,0]}.')
            continue
        
        idy = np.nanargmin(dt)
        dive_duration[idx] = dt[idy] / 1e9
        dive_time[idx] = (post_t[idy]+pre_t[idx])/2
        
        dxy_dvl[idx,0] = travel_e(post_t[idy]) - travel_e(pre_t[idx])
        dxy_dvl[idx,1] = travel_n(post_t[idy]) - travel_n(pre_t[idx])
        
        dxy_gps[idx,0] = (gps_postdive[idy,1] - gps_predive[idx,1]) * lon2m(gps_predive[idx,1],gps_predive[idx,2])
        dxy_gps[idx,1] = (gps_postdive[idy,2] - gps_predive[idx,2]) * lat2m(gps_predive[idx,1],gps_predive[idx,2])
    
    dac = (dxy_gps - dxy_dvl) / np.vstack([dive_duration,dive_duration]).T
    dive_time = pd.to_datetime(dive_time*1e9)
    
    return dac, dive_time, dive_duration, dxy_gps, dxy_dvl

def grid_shear(ADCP, xi, yi):
    currents = xr.Dataset(coords={"profile_index": xi, "depth": yi})
    
    x, y = xr.broadcast(ADCP.profile_number, ADCP.bin_depth)
    t, _ = xr.broadcast(ADCP.time, ADCP.bin_depth)
    x = x.values.flatten()
    y = y.values.flatten()
    t = t.values.astype('float').flatten()
    
    grid = lambda var_name, fn_name : grid2d(x, y, ADCP[var_name].values.flatten(), xi=xi, yi=yi, fn=fn_name)[0].T
    count_valid = lambda x : np.count_nonzero(np.isfinite(x))
    
    time = np.nanmean(grid2d(x, y, t, xi=xi, yi=yi, fn=np.nanmean)[0], axis=0)
    currents = currents.assign_coords(time=("profile_index", pd.to_datetime(time)))
    
    
    abs_vert_spd, _, YI = grid2d(ADCP.profile_number.values, ADCP.Depth.values, 
                                 np.abs(np.gradient(ADCP['Depth'].values,ADCP.time.values.astype('float')/1e9)), xi=xi, yi=yi, fn=np.nanmean)
    
    currents['time_in_bin'] = (('profile_index', 'depth'),  np.gradient(YI, axis=0).T / abs_vert_spd.T)
    
    print(f'Gridding eastward shear metrics...')
    currents['shear_E_mean'] = (('profile_index', 'depth'),   grid('Sh_E', np.nanmean))
    currents['shear_E_median'] = (('profile_index', 'depth'), grid('Sh_E', np.nanmedian))
    currents['shear_E_stddev'] = (('profile_index', 'depth'), grid('Sh_E', np.nanstd))
    currents['shear_E_count'] = (('profile_index', 'depth'),  grid('Sh_E', count_valid))
    
    print(f'Gridding northward shear metrics...')
    currents['shear_N_mean'] = (('profile_index', 'depth'),   grid('Sh_N', np.nanmean))
    currents['shear_N_median'] = (('profile_index', 'depth'), grid('Sh_N', np.nanmedian))
    currents['shear_N_stddev'] = (('profile_index', 'depth'), grid('Sh_N', np.nanstd))
    currents['shear_N_count'] = (('profile_index', 'depth'),  grid('Sh_N', count_valid))
    
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
    def integrate():
        def integrate_calc(Sh):
            # Integrate shear vertically
            _bd = ~np.isfinite(Sh)  # Preserve what are originally NaN values to recover later as need conversion to 0 for cumsum-
            Sh = np.nan_to_num(Sh)  # Replace NaNs with 0 for cumsum
            V = ( np.cumsum(Sh, axis=1) * np.gradient(currents.depth.values) )  # Cumulative sum of shear to recover velocity profile
            V[_bd] = np.nan  # Return NaNs to their rightful place.
            V = V - np.tile(
                np.nanmean(V, axis=1), (np.shape(V)[1], 1)
            ).T  # Make mean of baroclinic profiles equal to 0
            return V        
        currents['velocity_E_no_reference'] = (("profile_index", "depth"), integrate_calc(currents.shear_E_median.values))
        currents['velocity_N_no_reference'] = (("profile_index", "depth"), integrate_calc(currents.shear_N_median.values))
     
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

'''
Main
'''
def process(ADCP, gps_predive, gps_postdive, xi=1, yi=None, correct_shear_bias=True):
    """
    Calculate depth-resolved currents from shear data calculated by the gliderAD2CP.process_shear() function.

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

    __main__ runs the following methods in this order
    -------
    1. get_DAC
    2. grid_shear
    3. grid_velocity
    4. reference_velocity
    
    Optional steps:
    5. assess_shear_bias
    6. correct_shear_bias
    
    
    Methods
    -------
    grid2d
        Prints the person's name and age.
    plot_subsurface_movement
        Prints the person's name and age.
    get_DAC
        Prints the person's name and age.
    grid_shear
        Prints the person's name and age.
    integrate_shear
    reference_shear
    assess_shear_bias
    
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
    
#     currents = reference_shear(currents,DAC)
    
    
    return currents, DAC, dive_time
