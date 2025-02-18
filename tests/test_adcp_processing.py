import pandas as pd

from gliderad2cp import process_currents, process_shear, tools, process_bias
from gliderad2cp.download_example_data import load_sample_dataset
import numpy as np

options = tools.get_options(xaxis=2, verbose=True)
glider_pqt_path = load_sample_dataset(f"glider_profiles_160_to_210.pqt")

data = pd.read_parquet(glider_pqt_path)
data = data.rename({'dive_number': "dive_num"}, axis=1)
data = data[
    ["dead_reckoning", "time", "temperature", "salinity", "latitude", "pressure", "longitude", "profile_number",
     "declination", "dive_num"]]

adcp_path = load_sample_dataset(f"adcp_profiles_160_to_210.nc")


def test_process_shear():
    process_shear.process(adcp_path, glider_pqt_path, options=options)


def test_process_currents():
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
    currents, dac = process_currents.process(ds_adcp, gps_predive, gps_postdive, options)
    __ = process_bias.process(currents,options)
