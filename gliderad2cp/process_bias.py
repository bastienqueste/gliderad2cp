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

# def assess_bias(currents,DAC):
    
    
    
    

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