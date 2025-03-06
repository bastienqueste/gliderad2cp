"""
Correct velocity profiles obtained by ADCP gliders for shear-bias (as described in Todd et al., 2017; Sec.3.b.2 https://doi.org/10.1175/JTECH-D-16-0156.1)
Applies a different approach to Todd et al. in that we only utilise variance at depth as validation.
Instead, we aim to reduce the slope of the linear regression :

profile_mean(shear - shear_bias_correction) = profile_mean_naturally_occuring_shear + slope * displacement_through_water.

We determine two separate shear_bias_correction coefficients, one in the along-glider direction and one in the across glider direction.
This is akin to subtracting a portion of the glider's travel through water.

A second (potentially better) correction is also coded up which subtracts a shear_bias_correction which scales with glider speed. 
Preliminary work indicates that shear bias is instrument-velocity dependent, but until this is proven/published, this correction remains the
non-default correction. NB: the author of this code uses this one all the time.


gliderad2cp.process_bias
--------------------------
process
    The main function which will run through the whole shear bias correction automatically for you.
visualise
    Function outputting two plots which can be used to assess shear bias and its correction.
regress_bias
    Uses a minimisation function to identify optimal correction parameters
correct_bias
    Outputs DAC-referenced velocities that have been corrected for shear bias.
__linear_regression
    A weighted, linear least squares regression function with optional plotting and indexing capability.

Notes
-------
.process() runs the following functions in this order

1. visualise
2. regress_bias
3. correct_bias
4. visualise

"""

import warnings
import numpy as np
from scipy.stats import t
from scipy.optimize import fmin 
import matplotlib.pyplot as plt
from .tools import plog

warnings.filterwarnings(action="ignore", message="Mean of empty slice")
warnings.filterwarnings(action="ignore", message="invalid value encountered in divide")
warnings.filterwarnings(action="ignore", message="invalid value encountered in true_divide")
warnings.filterwarnings(action="ignore", message="Degrees of freedom <= 0 for slice.")

def _linear_regression(x,y,options,c=None,mask=True,plot=False):
    """
    Perform an optionally-weighted linear least squares regression on x and y, providing 95% confidence intervals and optionally plotting the data.

    Inputs
    ----------
    x,y,c : np.array
        Arrays containing x and y values for the linear regression, and an optional color argument for the scatter plot.
    mask : bool array
        Boolean array of same size as x, indicating which values to keep in the regression.
    plot : bool
        If true, produce a plot of the linear regression.
        
    Outputs
    -------
    slope : 
        Returns the slope of the linear regression
    
    """
    _gd = np.isfinite(x+y) & mask & (np.abs(y) < 0.01) & (np.abs(x) > 100)
    x = np.reshape(x[_gd],[-1,1])
    y = np.reshape(y[_gd],[-1,1])
    
    if options['weight_shear_bias_regression']:
        weights = x.flatten() # Weight by distance travelled, helps with weird shear spikes present on short/shallow dives
    else:
        weights = np.ones(len(x)) # Define equal weights
    
    # Add an intercept term (column of 1s) to X
    X = np.hstack( [np.ones((x.shape[0], 1)), x] )

    # Define weights (example: higher weight for larger x values)
    # weights = np.ones(len(x)) # Define weights for each observation
    W = np.diag(weights) # Create a diagonal matrix from weights

    # Calculate coefficients using the weighted normal equation
    # β = (X^T W X)^(-1) X^T W y
    X_T = X.T # Transpose of X
    beta = np.linalg.inv(X_T @ W @ X) @ (X_T @ W @ y) # Solve for β
    
    # Calculate residuals and residual variance
    y_pred = X @ beta  # Predicted values
    residuals = y - y_pred # Residuals
    n, p = X.shape # Number of observations, parameters
    sigma_squared = (residuals.T @ residuals) / (n - p) # Residual variance

    # Variance-Covariance Matrix
    var_cov_matrix = sigma_squared * np.linalg.inv(X_T @ X)

    # Extract standard error for the slope (2nd coefficient)
    slope_se = np.sqrt(var_cov_matrix[1, 1]) # Standard error of the slope
    intercept_se = np.sqrt(var_cov_matrix[0, 0]) # Standard error of the slope

    # Calculate 95% confidence interval
    t_critical = t.ppf(0.975, df=n-p) # Two-tailed t critical value
    
    intercept = beta[0][0].astype('float')
    slope = beta[1][0].astype('float')
    
    slope_ci = (slope - t_critical * slope_se, slope + t_critical * slope_se)
    intercept_ci = (intercept - t_critical * intercept_se, intercept + t_critical * intercept_se)
    
    
    if plot:
        if c is None:
            plt.plot(x, y, '.', color='k', markersize=1, label='Regressed data')
        else:
            plt.scatter(x, y, 3, c, label='Regressed data')
            plt.colorbar()
        plt.plot(x, intercept_ci[0] + slope_ci[0]*x, 'r', alpha=0.5, linestyle=':')
        plt.plot(x, intercept_ci[0] + slope_ci[1]*x, 'r', alpha=0.5, linestyle=':')
        plt.plot(x, intercept_ci[1] + slope_ci[0]*x, 'r', alpha=0.5, linestyle=':')
        plt.plot(x, intercept_ci[1] + slope_ci[1]*x, 'r', alpha=0.5, linestyle=':')
        
        plt.plot(x, intercept + slope*x, 'r', label='Regression')
        plt.legend()
        
        color = 'lightgrey'
        if np.abs(3*t_critical*slope_se) < np.abs(slope):
            color = 'orange'
        if np.abs(5*t_critical*slope_se) < np.abs(slope):
            color = 'red'
        plt.title(f"slope (95%): {slope:.2E} +/- {t_critical*slope_se:.2E}", color=color)

    return slope


def visualise(currents, options, plot_all=True):
    """
    Produce collection of plots used to estimate shear bias.


    Inputs
    ----------
    currents : xr.Dataset
        Dataset containing gridded shear and velocity data produced by the gliderad2cp.process_currents functions.
    options : dict
        Set of options for gliderAD2CP, created by the gliderad2cp.tools.get_options() function.
        
        
    Outputs
    -------
    None
    
    """
    corr_present = False
    if 'velocity_E_DAC_reference_sb_corrected' in currents: # Assume N also present
        corr_present = True
    
    upper,lower = options['shear_bias_regression_depth_slice']
    # Plot linear regression of shear versus displacement
    plt.figure(figsize=(10,20))
    plt.tight_layout()
    directions = ['N','E']
    idx = 0
    for l1 in directions:
        for l2 in directions:
            idx += 1
            plt.subplot(4,2,idx)
            if corr_present:
                _linear_regression(
                    currents[f'displacement_though_water_{l1}'].values, 
                    (currents[f'velocity_{l2}_DAC_reference_sb_corrected']).sel(depth=slice(upper,lower)).differentiate('depth').mean(dim='depth', skipna=True).values,
                    options,
                    plot=True
                )
            else:
                _linear_regression(
                    currents[f'displacement_though_water_{l1}'].values, 
                    (currents[f'velocity_{l2}_no_reference']).sel(depth=slice(upper,lower)).differentiate('depth').mean(dim='depth', skipna=True).values,
                    options,
                    plot=True
                )
            plt.xlabel(f'Displacement {l1} per profile')
            plt.ylabel(f'Mean shear {l2} per profile')
    
    if plot_all:
        # Plots of referenced velocity variance with depth
        plt.subplot(2,2,4)
        plt.plot(currents.velocity_E_DAC_reference.var(dim='profile_index'), currents.depth, ':r', label='E$_{orig}$')
        plt.plot(currents.velocity_N_DAC_reference.var(dim='profile_index'), currents.depth, ':b', label='N$_{orig}$')
        if corr_present:
            plt.plot(currents.velocity_E_DAC_reference_sb_corrected.var(dim='profile_index'), currents.depth, '-r', label='E$_{corr}$')
            plt.plot(currents.velocity_N_DAC_reference_sb_corrected.var(dim='profile_index'), currents.depth, '-b', label='N$_{corr}$')
        plt.legend()
        plt.xlabel('Velocity Variance')
        plt.ylabel('Depth')
        plt.gca().invert_yaxis()
        YL = plt.ylim()
        
        cl = np.max(np.abs([currents.velocity_E_DAC_reference.max(),currents.velocity_E_DAC_reference.min(),currents.velocity_N_DAC_reference.max(),currents.velocity_N_DAC_reference.min()]))
        cl = np.array([-1,1])*cl
        
        # Velocity sections
        plt.subplot(8,2,9)
        plt.set_cmap('bwr')
        plt.pcolormesh(currents.profile_index, currents.depth, currents.velocity_E_DAC_reference)
        plt.colorbar(label='E$_{orig}$')
        plt.clim(cl)
        XL = plt.xlim()
        plt.ylim(YL)

        if corr_present:
            plt.subplot(8,2,11)
            plt.set_cmap('bwr')
            plt.pcolormesh(currents.profile_index, currents.depth, currents.velocity_E_DAC_reference_sb_corrected)
            plt.colorbar(label='E$_{corr}$')
            plt.clim(cl)
            plt.xlim(XL)
            plt.ylim(YL)

        plt.subplot(8,2,13)
        plt.set_cmap('bwr')
        plt.pcolormesh(currents.profile_index, currents.depth, currents.velocity_N_DAC_reference)
        plt.colorbar(label='N$_{orig}$')
        plt.clim(cl)
        plt.xlim(XL)
        plt.ylim(YL)

        if corr_present:
            plt.subplot(8,2,15)
            plt.set_cmap('bwr')
            plt.pcolormesh(currents.profile_index, currents.depth, currents.velocity_N_DAC_reference_sb_corrected)
            plt.colorbar(label='N$_{corr}$')
            plt.xlim(XL)
            plt.ylim(YL)
            plt.clim(cl)
    
    plt.tight_layout()
    
    return None

def regress_bias(currents,options):
    """
    Determines shear bias correction coefficients by empirically minimimising the slope of various linear regressions.
    Linear regressions are performed on displacement through water and mean shear of the water column in combination of directions.
    

    Inputs
    ----------
    currents : xr.Dataset
        Dataset containing gridded shear and velocity data produced by the gliderad2cp.process_currents functions.
    options : dict
        Set of options for gliderAD2CP, created by the gliderad2cp.tools.get_options() function.
        
        
    Outputs
    -------
    bias_along_glider, bias_across_glider : float
        Floats which define the scaling factor applied to the shear bias correction in the along and across glider directions.
        
    """

    def __score(coeffs, *args):
        currs = correct_bias(args[0], args[1], coeffs[0], coeffs[1])
        directions = ['N','E']
        score = 0.
        idx = 0
        
        upper,lower = options['shear_bias_regression_depth_slice']
        for l1 in directions:
            for l2 in directions:
                idx += 1
                score = score + np.abs(_linear_regression(
                        currs[f'displacement_though_water_{l1}'].values, 
                        currs[f'velocity_{l2}_DAC_reference_sb_corrected'].sel(depth=slice(upper,lower)).differentiate('depth').mean(dim='depth', skipna=True).values,
                        options,
                    ))*1e9
        return score
    
    results = fmin(__score, [0,0], args=(currents,options), maxiter=300, disp=True, xtol=1e-9, ftol=1e-9)
    
    plog(f'Final results of shear bias regression: ')
    plog(f'    Along-glider bias :  {results[0]:.2E}')
    plog(f'    Across-glider bias : {results[1]:.2E}')
    
    return results[0], results[1]

    
def correct_bias(currents, options, bias_along_glider, bias_across_glider):
    """
    Calculates the artificial velocity profile created by the shear bias and outputs two new variables:
        "velocity_{direction}_DAC_reference_sb_corrected" - DAC referenced, shear-bias corrected velocity profiles.
        "shear_bias_velocity_{direction}" - the velocity profile error induced by the shear bias.
    

    Inputs
    ----------
    currents : xr.Dataset
        Dataset containing gridded shear and velocity data produced by the gliderad2cp.process_currents functions.
    options : dict
        Set of options for gliderAD2CP, created by the gliderad2cp.tools.get_options() function.
    bias_along_glider, bias_across_glider : float
        Floats which define the scaling factor applied to the shear bias correction in the along and across glider directions.
        
        
    Outputs
    -------
    currents : xr.Dataset
        Dataset containing gridded shear, statistical metrics, time spent per bin by the glider, and unreferenced velocity profiles.
        
    """
    mask = np.isfinite(currents[f'velocity_E_DAC_reference'].values)
    mask[mask == False] = np.nan
    
    if options['velocity_dependent_shear_bias_correction']:
        currents['shear_bias_velocity_E'] =  (currents.speed_through_water * (currents.heading_E*bias_along_glider - currents.heading_N*bias_across_glider) * currents.depth.differentiate('depth')).cumsum('depth') * mask
        currents['shear_bias_velocity_N'] =  (currents.speed_through_water * (currents.heading_N*bias_along_glider - currents.heading_E*bias_across_glider) * currents.depth.differentiate('depth')).cumsum('depth') * mask
    else:
        currents['shear_bias_velocity_E'] =  ((currents.heading_E*bias_along_glider - currents.heading_N*bias_across_glider) * currents.depth.differentiate('depth')).cumsum('depth') * mask
        currents['shear_bias_velocity_N'] =  ((currents.heading_N*bias_along_glider - currents.heading_E*bias_across_glider) * currents.depth.differentiate('depth')).cumsum('depth') * mask
    
    directions = ['N','E']
    for l in directions:
        currents[f'shear_bias_velocity_{l}'] = currents[f'shear_bias_velocity_{l}'] - currents[f'shear_bias_velocity_{l}'].mean(dim='depth')        
        currents[f'velocity_{l}_DAC_reference_sb_corrected'] = currents[f'velocity_{l}_DAC_reference'] - currents[f'shear_bias_velocity_{l}']
    
    return currents

def process(currents, options):
    """
    Correct velocity profiles obtained by ADCP gliders for shear-bias (as described in Todd et al., 2017; Sec.3.b.2 https://doi.org/10.1175/JTECH-D-16-0156.1)
    Applies a different approach to Todd et al. in that we only utilise variance at depth as validation.
    Instead, we aim to reduce the slope of the linear regression :

    profile_mean(shear - shear_bias_correction) = profile_mean_naturally_occuring_shear + slope * displacement_through_water.

    We determine two separate shear_bias_correction coefficients, one in the along-glider direction and one in the across glider direction.
    This is akin to subtracting a portion of the glider's travel through water.

    A second (potentially better) correction is also coded up which subtracts a shear_bias_correction which scales with glider speed. 
    Preliminary work indicates that shear bias is instrument-velocity dependent, but until this is proven/published, this correction remains the
    non-default correction. NB: the author of this code uses this one all the time.
    

    Inputs
    ----------
    currents : xr.Dataset
        Dataset containing gridded shear and velocity data produced by the gliderad2cp.process_currents functions.
    options : dict
        Set of options for gliderAD2CP, created by the gliderad2cp.tools.get_options() function.
        
        
    Outputs
    -------
    currents : xr.Dataset
        Dataset containing gridded shear, statistical metrics, time spent per bin by the glider, and unreferenced velocity profiles.
        

    Notes
    -------
    .process() runs the following functions in this order
    
    1. visualise -  Produce collection of plots used to visalise estimated shear bias
    2. regress_bias - Determine shear bias correction coefficients by empirically minimimising the slope of various linear regressions
    3. correct_bias -  Calculate the artificial velocity profile created by the shear bias and correct it in a new variable
    4. visualise - Produce collection of plots used to visalise estimated shear bias after correction

    """
    
    visualise(currents, options, plot_all=False)
    
    bias_along_glider, bias_across_glider = regress_bias(currents, options)
    
    currents = correct_bias(currents, options, bias_along_glider, bias_across_glider)
    
    visualise(currents, options)
    
    return currents