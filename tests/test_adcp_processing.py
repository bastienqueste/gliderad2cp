import pandas as pd
from gliderad2cp import process_currents, process_shear, process_bias, tools

def bypass():
    return True

#def test_load_csv():
#    df_pqt = gliderad2cp.process_adcp.load(
#        data_source.get_url("glider_profiles_160_to_210.pqt")
#    )
#    csv_file = "test.csv"
#    df_pqt.to_csv(csv_file)
#    gliderad2cp.process_adcp.load(csv_file)


# def test_processing():
#     profile_range = "160_to_210"
#     glider_pqt_path = data_source.get_url(f"glider_profiles_{profile_range}.pqt")
#     data_source.fetch(f"adcp_profiles_{profile_range}.nc")
#     adcp_path = str(data_source.path / f"adcp_profiles_{profile_range}.nc")
#     options = {
#         "correctADCPHeading": True,
#         "ADCP_discardFirstBins": 0,
#         "ADCP_correlationThreshold": 70,
#         "ADCP_amplitudeThreshold": 75,
#         "ADCP_velocityThreshold": 0.8,
#         "correctXshear": False,
#         "correctYshear": False,
#         "correctZshear": False,
#         "correctZZshear": False,
#         "ADCP_regrid_correlation_threshold": 20,
#     }
#     ADCP, data, options = process_adcp.load_adcp_glider_data(
#         adcp_path, glider_pqt_path, options
#     )
#     ADCP = process_adcp.remapADCPdepth(ADCP, options)
#     ADCP = process_adcp.correct_heading(ADCP, data, options)
#     ADCP = process_adcp.soundspeed_correction(ADCP)
#     ADCP = process_adcp.remove_outliers(ADCP, options)
#     #ADCP = process_adcp.correct_shear(ADCP, options)
#     #ADCP = process_adcp.correct_backscatter(ADCP, data)
#     ADCP = process_adcp.regridADCPdata(ADCP, options)
#     ADCP = process_adcp.calcXYZfrom3beam(ADCP, options)
#     ADCP = process_adcp.calcENUfromXYZ(ADCP, options)

#     # get your gridded shear here
#     xaxis, yaxis, taxis, days = process_adcp.grid_shear_data(data)
#     out = process_adcp.grid_data(ADCP, data, {}, xaxis, yaxis)

#     ds = process_adcp.make_dataset(out)

#     # integrate the gridded shear from here

#     extra_data = pd.read_parquet(glider_pqt_path)
#     extra_data.index = data.index
#     data["speed_vert"] = extra_data["speed_vert"]
#     data["speed_horz"] = extra_data["speed_horz"]
#     data["dead_reckoning"] = extra_data["dead_reckoning"]
#     data["nav_resource"] = extra_data["nav_resource"]
#     data["dive_number"] = extra_data["dive_number"]
#     data = process_adcp.get_DAC(ADCP, data)
#     dE, dN, dT = process_adcp.getSurfaceDrift(data)
#     ADCP = process_adcp.bottom_track(ADCP, adcp_path, options)
#     out = process_adcp.reference_shear(ADCP, data, xaxis, yaxis)
#     out = process_adcp.grid_data(ADCP, data, out, xaxis, yaxis)
#     out = process_adcp.calc_bias(out, yaxis)
#     ds = process_adcp.make_dataset(out)


# def test_wrapup_functions():
#     profile_range = "160_to_210"
#     glider_pqt_path = data_source.get_url(f"glider_profiles_{profile_range}.pqt")
#     data_source.fetch(f"adcp_profiles_{profile_range}.nc")
#     adcp_path = str(data_source.path / f"adcp_profiles_{profile_range}.nc")
#     options = {
#         "correctADCPHeading": True,
#         "ADCP_discardFirstBins": 0,
#         "ADCP_correlationThreshold": 70,
#         "ADCP_amplitudeThreshold": 75,
#         "ADCP_velocityThreshold": 0.8,
#         "correctXshear": False,
#         "correctYshear": False,
#         "correctZshear": False,
#         "correctZZshear": False,
#         "ADCP_regrid_correlation_threshold": 20,
#     }
#     ADCP, data, options = process_adcp.load_adcp_glider_data(
#         adcp_path, glider_pqt_path, options
#     )
#     ADCP, data = process_adcp.shear_from_adcp(adcp_path, glider_pqt_path, options=options)
#     ds = process_adcp.grid_shear(ADCP, data)
#     process_adcp.velocity_from_shear(adcp_path, glider_pqt_path, ADCP, data, options=options)
