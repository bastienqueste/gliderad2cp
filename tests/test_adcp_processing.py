from gliderad2cp import process_currents, process_shear, tools, process_bias
import xarray as xr
import numpy as np


def test_process():
    options = tools.get_options(xaxis=2, verbose=True)

    data = xr.open_dataset('/home/callum/Documents/hack/data_scratch/SEA055_M82.nc')
    data = data[
        ["dead_reckoning", "time", "temperature", "salinity", "latitude", "pressure", "longitude", "profile_index",
         "declination", "dive_num"]]
    data['profile_number'] = data['profile_index']
    adcp_path = f"/home/callum/Documents/hack/data_scratch/sea055_M82.ad2cp.00000.nc"

    ds_adcp = process_shear.process(adcp_path, data, options)


    gps_predive = []
    gps_postdive = []

    dives = np.round(np.unique(data.dive_num))

    _idx = np.arange(len(data.dead_reckoning.values))
    dr = np.sign(np.gradient(data.dead_reckoning.values))

    for dn in dives:
        _gd = data.dive_num.values == dn
        if all(np.unique(dr[_gd]) == 0):
            continue

        _post = -dr.copy()
        _post[_post != 1] = np.nan
        _post[~_gd] = np.nan

        _pre = dr.copy()
        _pre[_pre != 1] = np.nan
        _pre[~_gd] = np.nan

        if any(np.isfinite(_post)):
            # The last -1 value is when deadreckoning is set to 0, ie. GPS fix. This is post-dive.
            last  = int(np.nanmax(_idx * _post))
            gps_postdive.append(np.array([data.time[last].values, data.longitude[last].values, data.latitude[last].values]))

        if any(np.isfinite(_pre)):
            # The first +1 value is when deadreckoning is set to 1, the index before that is the last GPS fix. This is pre-dive.
            first = int(np.nanmin(_idx * _pre))-1 # Note the -1 here.
            gps_predive.append(np.array([data.time[first].values, data.longitude[first].values, data.latitude[first].values]))

    gps_predive = np.vstack(gps_predive)
    gps_postdive = np.vstack(gps_postdive)



    currents, DAC = process_currents.process(ds_adcp, gps_predive, gps_postdive, options)
    bias_along_glider, bias_across_glider = process_bias.regress_bias(currents, options)

    currents = process_bias.correct_bias(currents, options, bias_along_glider, bias_across_glider)
