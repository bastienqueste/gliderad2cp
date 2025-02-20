"""
gliderad2cp.tools
--------------------------
Helper funcitons for gliderad2cp
"""
import pandas as pd
import numpy as np
import logging
from datetime import datetime as dt
from scipy.interpolate import interp1d


def get_options(verbose=True, **kwargs):
    """
    Returns a dictionary containing options compatible with gliderAD2CP.
    Run with no kwargs to output a default dictionary, or with kwargs to return a full dictionary with the chose options
    """
    options = {
        'correct_compass_calibration' : [False, 'compass correction algorithm is awaiting publication and will be added upon acceptance. Contact Bastien Queste if you require.'],
        'shear_to_velocity_method' : ['integrate',],
        'ADCP_mounting_direction' : ['auto', 'top', 'bottom'],
        'QC_correlation_threshold' : [80, 'minimum acceptable along-beam correlation value.'],
        'QC_amplitude_threshold' : [80, 'maximum acceptable along-beam amplitude.'],
        'QC_velocity_threshold' : [0.8, 'maximum acceptable along-beam velocity in m.s-1.'],
        'QC_SNR_threshold' : [3, 'minimum acceptable dB above the noise floor.'],
        'velocity_regridding_distance_from_glider' : ['auto', 'array of depth-offsets from the glider, in m, at which to interpolate beam velocities onto isobars to avoid shear-smearing. Negative for bottom-mounted ADCPs.'],
        'xaxis' : [1, 'x-axis resolution in number of profiles of the final gridded products.'],
        'yaxis' : [None, 'If None: ADCP cell size. If int: y-axis resolution in metres of the final gridded products.'],
        'weight_shear_bias_regression' : [False, True, 'Give greater weight to dives with greater travel distance which can increase signal to noise.'],
        'velocity_dependent_shear_bias_correction' : [False, True, 'Determine velocity dependent shear-bias correction coefficients rather than constant coefficients.'],
        'shear_bias_regression_depth_slice' : [(0, 1000), 'A tuple containing the upper and lower depth limits over which to determine shear bias. Helpful to avoid increased noise due to surface variability. For deep diving gliders (500,1000) is good.'],
        'pitch_offset' : [0, 'value to be added to pitch to correct for transducer-compass misalignment'],
        'roll_offset' : [0, 'value to be added to roll to correct for transducer-compass misalignment'],
        'heading_offset' : [0, 'value to be added to heading to correct for transducer-compass misalignment'],
        }
        
    default = dict()
    
    if verbose:
        print('Available options are: ')
        for i, (k,v) in enumerate(options.items()):
            print('    %s : %s' % (k, v))
        print('Default setting is the first value.')
        print('Options prefaced with a minus are not yet implemented.')
        
    for i, (k,v) in enumerate(options.items()):
        default[k] = v[0]
    
    for i, (k,v) in enumerate(kwargs.items()):
        default[k] = v    
    
    return default

"""
Quality of life functions
"""
_log = logging.getLogger(__name__)
def plog(msg):
    """
    Output information to the logger with timestamp.
    """
    print(msg)
    _log.info(str(dt.now().replace(microsecond=0)) + ' : ' + str(msg))
    return None

def grid2d(x, y, v, xi=1, yi=1, fn='median'):
    """
    Quick data binning function relying on pandas.
    x,y,v are flat np.arrays of x, y coordinates and values at these points.
    xi and yi are either np.arrays of bins to be binned into, or the spacing used between min and max of x or y respectively.
    fn defines the function applied to all points from v that fall into the same bin.
    """
    if np.size(xi) == 1:
        xi = np.arange(np.nanmin(x), np.nanmax(x) + xi, xi)
    if np.size(yi) == 1:
        yi = np.arange(np.nanmin(y), np.nanmax(y) + yi, yi)

    raw = pd.DataFrame({'x': x, 'y': y, 'v': v}).dropna()

    grid = np.full([np.size(yi), np.size(xi)], np.nan)

    raw['xbins'], xbin_iter = pd.cut(raw.x, xi, retbins=True, labels=False, right=False)
    raw['ybins'], ybin_iter = pd.cut(raw.y, yi, retbins=True, labels=False, right=False)

    _tmp = raw.groupby(['xbins', 'ybins'])['v'].agg(fn)
    grid[
        _tmp.index.get_level_values(1).astype(int),
        _tmp.index.get_level_values(0).astype(int),
    ] = _tmp.values

    XI, YI = np.meshgrid(xi, yi, indexing='ij')
    return grid, XI.T, YI.T

def interp(x, y, xi):
    """
    Interpolation shortcut that removes NaN data and ignores boundary errors.
    """
    _gg = np.isfinite(x + y)
    return interp1d(x[_gg], y[_gg], bounds_error=False, fill_value=np.nan)(xi)