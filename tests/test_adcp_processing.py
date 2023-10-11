import xarray as xr
import numpy as np
import pandas as pd
from pathlib import Path
import sys
import pooch
module_dir = Path(__file__).parent.parent.absolute()
sys.path.append(str(module_dir))
from gliderad2cp import process_adcp


def test_processing():
    profile_range = "160_to_210"
    data_source = pooch.create(
        path=pooch.os_cache("gliderad2cp"),
        base_url="https://zenodo.org/record/8431329/files/",
        registry={
            f"adcp_profiles_{profile_range}.nc": "sha256:323ff3cc6402b6c7034a57369ee637c1398af38c2d5f876c0456dbbf9928ab6f",
            f"glider_profiles_{profile_range}.pqt": "sha256:ee83f2d0f3bac1c937da4115c5904b7429a3531406654747dd64845a3aeeb7b5",
            f"processed_velocity_{profile_range}.nc": "sha256:cb6f0ccd580db111ad6b54da9f8db831632f461740eb78b26141964b6abe97b6",
        },
    )
    glider_pqt_path = data_source.get_url(f"glider_profiles_{profile_range}.pqt")
    data_source.fetch(f"adcp_profiles_{profile_range}.nc")
    adcp_path = str(data_source.path / f"adcp_profiles_{profile_range}.nc")
    options = {
        'debug_plots': False,
        'correctADCPHeading': True,
        'ADCP_discardFirstBins': 0,
        'ADCP_correlationThreshold': 70,
        'ADCP_amplitudeThreshold': 75,
        'ADCP_velocityThreshold': 0.8,
        'correctXshear': False,
        'correctYshear': False,
        'correctZshear': False,
        'correctZZshear': False,
        'ADCP_regrid_correlation_threshold': 20,
    }
    ADCP, data, options = process_adcp.load_adcp_glider_data(adcp_path, glider_pqt_path, options)
    ADCP = process_adcp.remapADCPdepth(ADCP, options)
    ADCP = process_adcp.correct_heading(ADCP, data, options)
    ADCP = process_adcp.soundspeed_correction(ADCP)
    ADCP = process_adcp.remove_outliers(ADCP, options)
    ADCP = process_adcp.correct_shear(ADCP, options)
    ADCP = process_adcp.correct_backscatter(ADCP, data, options)
    ADCP = process_adcp.regridADCPdata(ADCP, options)
    ADCP = process_adcp.calcXYZfrom3beam(ADCP, options)
    ADCP = process_adcp.calcENUfromXYZ(ADCP, data, options)

    # get your gridded shear here
    xaxis, yaxis, taxis, days = process_adcp.grid_shear_data(ADCP, data, options)
    out = process_adcp.grid_data(ADCP, data, {}, xaxis, yaxis)

    ds = process_adcp.make_dataset(out)
    ds_min = ds[['Sh_E', 'Sh_N', 'Sh_U']]
    data_source.fetch(f"processed_velocity_{profile_range}.nc")
    ds_min_test = xr.open_dataset(str(data_source.path / f"processed_velocity_{profile_range}.nc"))
    for var in list(ds_min):
        assert np.allclose(ds_min[var], ds_min_test[var], equal_nan=True, atol=1e-7, rtol=1e-3)
    # integrate the gridded shear from here

    extra_data = pd.read_parquet(glider_pqt_path)
    extra_data.index = data.index
    data["speed_vert"] = extra_data["speed_vert"]
    data["speed_horz"] = extra_data["speed_horz"]
    data["DeadReckoning"] = extra_data["DeadReckoning"]
    data["NAV_RESOURCE"] = extra_data["NAV_RESOURCE"]
    data["diveNum"] = extra_data["diveNum"]
    data = process_adcp.get_DAC(ADCP, data, options)
    dE, dN, dT = process_adcp.getSurfaceDrift(data, options)
    ADCP = process_adcp.bottom_track(ADCP, adcp_path, options)
    out = process_adcp.reference_shear(ADCP, data, dE, dN, dT, xaxis, yaxis, taxis, options)
    out = process_adcp.grid_data(ADCP, data, out, xaxis, yaxis)
    out = process_adcp.calc_bias(out, yaxis, taxis, days, options)

    ds = process_adcp.make_dataset(out)
    ds_min = ds[['ADCP_E', 'ADCP_N']]
    for var in list(ds_min):
        assert np.allclose(ds_min[var], ds_min_test[var], equal_nan=True, atol=1e-7, rtol=1e-3)

