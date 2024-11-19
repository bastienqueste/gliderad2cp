import pandas as pd
import numpy as np
import gliderad2cp.process_adcp
from gliderad2cp.download_example_data import data_source
from gliderad2cp import process_adcp, process_currents


def test_load_csv():
    df_pqt = gliderad2cp.process_adcp.load(
        data_source.get_url("glider_profiles_160_to_210.pqt")
    )
    csv_file = "test.csv"
    df_pqt.to_csv(csv_file)
    gliderad2cp.process_adcp.load(csv_file)


def test_processing():
    profile_range = "160_to_210"
    glider_pqt_path = data_source.get_url(f"glider_profiles_{profile_range}.pqt")
    data_source.fetch(f"adcp_profiles_{profile_range}.nc")
    adcp_path = str(data_source.path / f"adcp_profiles_{profile_range}.nc")
    options = {
        "correctADCPHeading": True,
        "ADCP_discardFirstBins": 0,
        "ADCP_correlationThreshold": 70,
        "ADCP_amplitudeThreshold": 75,
        "ADCP_velocityThreshold": 0.8,
        "correctXshear": False,
        "correctYshear": False,
        "correctZshear": False,
        "correctZZshear": False,
        "ADCP_regrid_correlation_threshold": 20,
    }
    df = pd.read_parquet(glider_pqt_path).set_index('time')
    data = df.to_xarray()
    data['time'] = pd.DatetimeIndex(data['time'].values)
    data['dive_num'] = data['dive_number']
    data['dead_reckoning'].values = data['dead_reckoning'].astype(int).values

    ds_adcp, df_glider = process_adcp.shear_from_adcp(adcp_path, glider_pqt_path, options=options)

    gps_predive = []
    gps_postdive = []

    dives = np.round(np.unique(data.dive_num))

    _idx = np.arange(len(data.dead_reckoning.values))
    for dn in dives:
        _dn = data.dive_num.values == dn
        _dr = data.dead_reckoning.values == 0
        _gd = (_dn & _dr).astype('float')
        _gd[_gd < 1] = np.nan

        if all(np.isnan(_gd)):
            continue

        first = int(np.nanmin(_idx * _gd))
        last = int(np.nanmax(_idx * _gd))

        gps_predive.append(np.array([data.time[last].values, data.longitude[last].values, data.latitude[last].values]))
        gps_postdive.append(
            np.array([data.time[first].values, data.longitude[first].values, data.latitude[first].values]))
    gps_predive = np.vstack(gps_predive)
    gps_postdive = np.vstack(gps_postdive)

    currents = process_currents.process(ds_adcp, gps_predive, gps_postdive, xi=1, yi=6)



