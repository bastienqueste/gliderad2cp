import pandas as pd

from gliderad2cp import process_currents, process_shear, tools, process_bias
from gliderad2cp.download_example_data import data_source
import numpy as np

options = tools.get_options(xaxis=2, verbose=True)
profile_range = "160_to_210"
glider_pqt_path = data_source.get_url(f"glider_profiles_{profile_range}.pqt")

data = pd.read_parquet(glider_pqt_path)
data = data.rename({'dive_number': "dive_num"}, axis=1)
data = data[
    ["dead_reckoning", "time", "temperature", "salinity", "latitude", "pressure", "longitude", "profile_number",
     "declination", "dive_num"]]

data_source.fetch(f"adcp_profiles_{profile_range}.nc")
adcp_path = str(data_source.path / f"adcp_profiles_{profile_range}.nc")


def test_process_shear():
    adcp = process_shear.process(adcp_path, glider_pqt_path, options=options)


def test_process_currents():
    adcp = process_shear.process(adcp_path, glider_pqt_path, options=options)

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
            gps_postdive.append(np.array([data.time[last].to_datetime64(), data.longitude[last], data.latitude[last]]))

        if any(np.isfinite(_pre)):
            # The first +1 value is when deadreckoning is set to 1, the index before that is the last GPS fix. This is pre-dive.
            first = int(np.nanmin(_idx * _pre))-1 # Note the -1 here.
            gps_predive.append(np.array([data.time[first].to_datetime64(), data.longitude[first], data.latitude[first]]))

    gps_predive = np.vstack(gps_predive)
    gps_postdive = np.vstack(gps_postdive)
    ds_adcp = process_shear.process(adcp_path, data, options)
    currents, DAC = process_currents.process(ds_adcp, gps_predive, gps_postdive, options)
    bias_along_glider, bias_across_glider = process_bias.regress_bias(currents, options)

    currents = process_bias.correct_bias(currents, options, bias_along_glider, bias_across_glider)

