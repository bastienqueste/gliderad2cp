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
    dr = np.ones(len(data.dead_reckoning.values))
    for dn in dives:
        _gd = data.dive_num.values == dn
        _post = dr.copy()
        _post[~_gd] = np.nan

        _pre = dr.copy()
        _pre[~_gd] = np.nan

        last  = int(np.nanmax(_idx * _post)) - 10
        gps_postdive.append(np.array([data.time[last].to_datetime64(), data.longitude[last], data.latitude[last]]))

        first = int(np.nanmin(_idx * _pre)) + 10
        gps_predive.append(np.array([data.time[first].to_datetime64(), data.longitude[first], data.latitude[first]]))

    gps_predive = np.vstack(gps_predive)
    gps_postdive = np.vstack(gps_postdive)
    ds_adcp = process_shear.process(adcp_path, data, options)
    currents, DAC = process_currents.process(ds_adcp, gps_predive, gps_postdive, options)
    bias_along_glider, bias_across_glider = process_bias.regress_bias(currents, options)

    currents = process_bias.correct_bias(currents, options, bias_along_glider, bias_across_glider)

