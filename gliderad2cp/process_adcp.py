import datetime
import json
import logging
import warnings
from datetime import datetime as dt
from glob import glob
from urllib import request

import gsw
import numpy as np
import pandas as pd
import xarray as xr
from scipy.interpolate import interp1d
from scipy.optimize import fmin

warnings.filterwarnings(action="ignore", message="Mean of empty slice")
warnings.filterwarnings(action="ignore", message="invalid value encountered in divide")
warnings.filterwarnings(
    action="ignore", message="invalid value encountered in true_divide"
)
warnings.filterwarnings(action="ignore", message="Degrees of freedom <= 0 for slice.")

y_res = 1  # TODO: move to options
_log = logging.getLogger(__name__)
default_options = {
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
    "top_mounted": False
}


def get_declination(data, key):
    """
    Function retrieves declination data from NOAA for the average lon, lat and datetime of glider data.
    Requires an API key. Register at https://www.ngdc.noaa.gov/geomag/CalcSurvey.shtml

    :param data: pd.DataFrame of glider data including time, longitude and latitude
    :param key: API key for the NOAA geomag service
    :return: dataframe with added declination column
    """
    if "declination" in list(data):
        _log.info("declination data already present")
        return data
    time = data.time.mean()
    year = time.year
    month = time.month
    day = time.day
    lat = data.latitude.mean()
    lon = data.longitude.mean()
    url = (
        f"https://www.ngdc.noaa.gov/geomag-web/calculators/calculateDeclination?startYear={year}&startMonth={month}"
        f"&startDay={day}&lat1={lat}&lon1={lon}&key={key}&resultFormat=json"
    )
    result = json.load(request.urlopen(url))
    declination = result["result"][0]["declination"]
    _log.info(f"declination of {declination} added to data")
    data["declination"] = declination


def load(glider_file):

    """
    Load glider data for processing by gliderad2cp.

    :param glider_file: str or Path. Relative or absolute path to a csv or pqt file of glider data. Or an xarray dataset
    :returns: a pandas dataframe of glider data
    """
    # Read in pyglider nercdf. Note does not check for exact pyglider format, simply assumes
    if type(glider_file) is xr.core.dataset.Dataset:
        _log.info("Input recognised as xarray dataset")
        data = glider_file.to_dataframe()
        data["time"] = data.index
        data["date_float"] = data["time"].values.astype("float")
        if "profile_index" in data:
            data["profile_number"] = data["profile_index"]
        p = data["pressure"]
        SA = gsw.conversions.SA_from_SP(data.salinity, p, data.longitude, data.latitude)
        CT = gsw.CT_from_t(SA, data["temperature"], p)
        data["soundspeed"] = gsw.sound_speed(SA, CT, p)
        data.index.name = None
        return data       

    # Read in csv file or pandas file
    if str(glider_file)[-4:] == ".csv":
        data = pd.read_csv(glider_file, parse_dates=["time"])
    else:
        data = pd.read_parquet(glider_file)
    _log.info(f"Loaded {glider_file}")
    sel_cols = [
        "time",
        "temperature",
        "salinity",
        "latitude",
        "longitude",
        "profile_number",
        "declination",
        "pressure",
        "dive_number"
    ]
    data = data[sel_cols]
    time_ms = data.time.values
    if time_ms.dtype != "<M8[ns]":
        divisor = 1e3
        if time_ms.dtype == "<M8[us]":
            divisor = 1e6
        time_float = time_ms.astype("float") / divisor
        base = datetime.datetime(1970, 1, 1)
        time_ms = []
        for seconds in time_float:
            time_ms.append(base + datetime.timedelta(seconds=seconds + 0.0000001))
        time_nanoseconds = pd.to_datetime(time_ms)
        data["time"] = time_nanoseconds
    data["date_float"] = data["time"].values.astype("float")
    p = data["pressure"]
    SA = gsw.conversions.SA_from_SP(data.salinity, p, data.longitude, data.latitude)
    CT = gsw.CT_from_t(SA, data["temperature"], p)
    data["soundspeed"] = gsw.sound_speed(SA, CT, p)
    data.index = data.time
    data.index.name = None
    return data


def grid2d(x, y, v, xi=1, yi=1, fn="median"):
    if np.size(xi) == 1:
        xi = np.arange(np.nanmin(x), np.nanmax(x) + xi, xi)
    if np.size(yi) == 1:
        yi = np.arange(np.nanmin(y), np.nanmax(y) + yi, yi)

    raw = pd.DataFrame({"x": x, "y": y, "v": v}).dropna()

    grid = np.full([np.size(yi), np.size(xi)], np.nan)

    raw["xbins"], xbin_iter = pd.cut(raw.x, xi, retbins=True, labels=False, right=False)
    raw["ybins"], ybin_iter = pd.cut(raw.y, yi, retbins=True, labels=False, right=False)

    _tmp = raw.groupby(["xbins", "ybins"])["v"].agg(fn)
    grid[
        _tmp.index.get_level_values(1).astype(int),
        _tmp.index.get_level_values(0).astype(int),
    ] = _tmp.values

    XI, YI = np.meshgrid(xi, yi, indexing="ij")
    return grid, XI.T, YI.T


def RunningMean(x, N):
    # is np.convolve better/faster?

    grid = np.ones((len(x) + 2 * N, 1 + 2 * N)) * np.nan
    for istep in range(np.shape(grid)[1]):
        grid[istep: len(x) + istep, istep] = x
    return np.nanmean(grid, axis=1)[N:-N]


def interp(x, y, xi):
    _gg = np.isfinite(x + y)
    return interp1d(x[_gg], y[_gg], bounds_error=False, fill_value=np.nan)(xi)


def rmsd(x):
    return np.sqrt(np.nanmean(x**2))


def plog(msg):
    _log.info(str(dt.now().replace(microsecond=0)) + " : " + msg)
    return None


def load_adcp_glider_data(adcp_file_path, glider_file_path, options):
    glider = load(glider_file_path)

    ADCP = xr.open_mfdataset(adcp_file_path, group="Data/Average")
    ADCP_settings = xr.open_mfdataset(glob(adcp_file_path)[0], group="Config")
    ADCP.attrs = ADCP_settings.attrs
    adcp_time_float = ADCP.time.values.astype("float")
    ADCP = ADCP.drop_vars(['MatlabTimeStamp'])
    # if ADCP.time.dtype == '<M8[ns]':
    #    adcp_time_float = adcp_time_float / 1e6
    plog("Finished loading ADCP data")
    # Coordinates
    ADCP = ADCP.assign_coords(
        Latitude=(
            "time",
            interp(glider["date_float"], glider["latitude"], adcp_time_float),
        )
    )
    ADCP = ADCP.assign_coords(
        Longitude=(
            "time",
            interp(glider["date_float"], glider["longitude"], adcp_time_float),
        )
    )

    # Profile and depth
    ADCP = ADCP.assign_coords(
        profile_number=(
            "time",
            np.round(
                interp(glider["date_float"], glider["profile_number"], adcp_time_float)
            ),
        )
    )
    ADCP = ADCP.assign_coords(
        Depth=("time", -gsw.z_from_p(ADCP["Pressure"].values, ADCP["Latitude"].values))
    )

    # overwrite temperature with glider temperature?
    ADCP["salinity"] = (
        "time",
        interp(
            glider["time"].values.astype("float"),
            glider["salinity"],
            ADCP.time.values.astype("float"),
        ),
    )
    ADCP["declination"] = (
        "time",
        interp(
            glider["time"].values.astype("float"),
            glider["declination"],
            ADCP.time.values.astype("float"),
        ),
    )
    ADCP["glider_soundspeed"] = (
        "time",
        interp(
            glider["date_float"].values,
            glider["soundspeed"],
            ADCP.time.values.astype("float"),
        ),
    )

    # Get rid of pointless dimensions and make them coordinates instead
    ADCP = ADCP.assign_coords(
        bin=("Velocity Range", np.arange(len(ADCP["Velocity Range"].values)))
    )
    ADCP = ADCP.swap_dims({"Velocity Range": "bin"})

    ADCP = ADCP.assign_coords(
        bin=("Correlation Range", np.arange(len(ADCP["Correlation Range"].values)))
    )
    ADCP = ADCP.swap_dims({"Correlation Range": "bin"})

    ADCP = ADCP.assign_coords(
        bin=("Amplitude Range", np.arange(len(ADCP["Amplitude Range"].values)))
    )
    ADCP = ADCP.swap_dims({"Amplitude Range": "bin"})

    plog("Added glider variables")

    # Detect if ADCP is top mounted using the magnetometer
    if ADCP.AccelerometerZ.values.mean() < 0:
        options["top_mounted"] = False
    else:
        options["top_mounted"] = True
    plog(f'top mounted: {options["top_mounted"]}')
    return ADCP, glider, options


def remapADCPdepth(ADCP, options):
    """
    Remaps velocity estimates to depth bins using the attitude of the glider

    :param ADCP: xr.DataSet of ADCP data
    :param options: options dictionary
    :return: ADCP with velocities remapped to depth bins in the variables "Dx"where x is the beam number
    """
    # All *_Range coordinates are distance along beam. Verified with data.
    if options["top_mounted"]:
        direction = 1
        theta_rad_1 = np.arccos(
            np.cos(np.deg2rad(47.5 - ADCP["Pitch"])) * np.cos(np.deg2rad(ADCP["Roll"]))
        )
        theta_rad_2 = np.arccos(
            np.cos(np.deg2rad(25 - ADCP["Roll"])) * np.cos(np.deg2rad(ADCP["Pitch"]))
        )
        theta_rad_3 = np.arccos(
            np.cos(np.deg2rad(47.5 + ADCP["Pitch"])) * np.cos(np.deg2rad(ADCP["Roll"]))
        )
        theta_rad_4 = np.arccos(
            np.cos(np.deg2rad(25 + ADCP["Roll"])) * np.cos(np.deg2rad(ADCP["Pitch"]))
        )
    else:
        direction = -1
        theta_rad_1 = np.arccos(
            np.cos(np.deg2rad(47.5 + ADCP["Pitch"])) * np.cos(np.deg2rad(ADCP["Roll"]))
        )
        theta_rad_2 = np.arccos(
            np.cos(np.deg2rad(25 + ADCP["Roll"])) * np.cos(np.deg2rad(ADCP["Pitch"]))
        )
        theta_rad_3 = np.arccos(
            np.cos(np.deg2rad(47.5 - ADCP["Pitch"])) * np.cos(np.deg2rad(ADCP["Roll"]))
        )
        theta_rad_4 = np.arccos(
            np.cos(np.deg2rad(25 - ADCP["Roll"])) * np.cos(np.deg2rad(ADCP["Pitch"]))
        )
    # Upward facing ADCP, so beam 1 ~= 30 deg on the way up, beam 3 on the way down, when flying at 17.4 degree pitch.
    # Returns angles of each beam from the UP direction

    # z_bin_distance = ADCP["Velocity Range"].values
    
    # THIS IS THE  BEAM REMAPPING CHANGE HERE :: the addition of the cos factor
    z_bin_distance = ADCP['Velocity Range'].values/np.cos(np.deg2rad(30.1))

    ADCP["D1"] = (
        ["time", "bin"],
        np.tile(ADCP["Depth"], (len(ADCP.bin), 1)).T
        - direction
        * np.tile(z_bin_distance, (len(ADCP.time), 1))
        * np.tile(np.cos(theta_rad_1), (len(ADCP.bin), 1)).T,
    )
    ADCP["D2"] = (
        ["time", "bin"],
        np.tile(ADCP["Depth"], (len(ADCP.bin), 1)).T
        - direction
        * np.tile(z_bin_distance, (len(ADCP.time), 1))
        * np.tile(np.cos(theta_rad_2), (len(ADCP.bin), 1)).T,
    )
    ADCP["D3"] = (
        ["time", "bin"],
        np.tile(ADCP["Depth"], (len(ADCP.bin), 1)).T
        - direction
        * np.tile(z_bin_distance, (len(ADCP.time), 1))
        * np.tile(np.cos(theta_rad_3), (len(ADCP.bin), 1)).T,
    )
    ADCP["D4"] = (
        ["time", "bin"],
        np.tile(ADCP["Depth"], (len(ADCP.bin), 1)).T
        - direction
        * np.tile(z_bin_distance, (len(ADCP.time), 1))
        * np.tile(np.cos(theta_rad_4), (len(ADCP.bin), 1)).T,
    )
    plog(
        "Please verify the location of ADCP velocity bins relative to pitch and depth of the sensor to make sure ADCP direction has been properly identified."
    )
    return ADCP


def _heading_correction(ADCP, glider, options):  # TODO: replace with external function
    # # Get local geomagnetic target strength:
    def getGeoMagStrength():
        lat = np.nanmedian(glider.latitude)
        lon = np.nanmedian(glider.longitude)
        date = pd.to_datetime(np.nanmean(glider.time.values.astype("float")))
        year = date.year
        month = date.month
        day = date.day

        url = str(
            "http://geomag.bgs.ac.uk/web_service/GMModels/igrf/13/?"
            + "latitude="
            + str(lat)
            + "&longitude="
            + str(lon)
            + "&date="
            + str(year)
            + "-"
            + str(month)
            + "-"
            + str(day)
            + "&resultFormat=csv"
        )
        magdata = request.urlopen(url)

        string = "empty"
        while not not string:
            out = magdata.readline().decode("utf-8")
            if 'total-intensity units="nT"' in out:
                string = out
                _log.info(out)
                break
        target = float(string.split(">")[1].split("<")[0])
        nT2milligauss = (
            10**-9 * 10000 * 1000
        )  # To tesla, then to gauss then to millgauss
        _log.info("Target = " + str(target * nT2milligauss))
        return target * nT2milligauss

    target = getGeoMagStrength()

    if options["top_mounted"]:
        sign = -1
    else:
        sign = 1

    MagX = ADCP["MagnetometerX"]
    MagY = ADCP["MagnetometerY"]
    MagZ = ADCP["MagnetometerZ"]

    simple = False
    verysimple = False
    softonly = False

    roll = ADCP["Roll"]
    pitch = ADCP["Pitch"]

    def norm(x, y, z):
        return np.sqrt(x**2 + y**2 + z**2)

    def rmsd(x, y, z):
        return np.sqrt(np.mean((norm(x, y, z) - target) ** 2))

    def cosd(x):
        return np.cos(np.deg2rad(x))

    def sind(x):
        return np.sin(np.deg2rad(x))

    def atan2d(x, y):
        return np.rad2deg(np.arctan2(x, y))

    def rot_x(x, y, z):
        return (
            x * cosd(pitch)
            + y * sind(roll) * sind(pitch)
            + z * cosd(roll) * sind(pitch)
        )

    def rot_y(y, z):
        return y * cosd(roll) - z * sind(roll)

    def wrap(x):
        return (x + 360) % 360

    def heading(x, y, z):
        return wrap(
            atan2d(rot_x(x, sign * y, sign * z), rot_y(sign * y, sign * z)) - 90
        )

    def calibrate(x, y, z, coeffs):
        if simple:
            coeffs[[1, 2, 3, 5, 6, 7]] = 0
        if verysimple:
            coeffs[:9] = 0
            coeffs[[0, 4, 8]] = 1
        if softonly:
            coeffs[-3:] = 0

        A = np.reshape(coeffs[:9], (3, 3))
        B = coeffs[-3:]
        out = A @ np.array([x - B[0], y - B[1], z - B[2]])
        return out[0, :], out[1, :], out[2, :]

    def minimisation(coeffs):
        x, y, z = calibrate(MagX, MagY, MagZ, coeffs)
        return rmsd(x, y, z)

    coeffs = fmin(
        minimisation,
        np.array([1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0]),
        disp=False,
        retall=False,
    )
    _log.info(np.reshape(coeffs[:9], (3, 3)))
    _log.info(coeffs[-3:])

    magx, magy, magz = calibrate(MagX.values, MagY.values, MagZ.values, coeffs)
    cal_heading = heading(magx, magy, magz)
    return cal_heading


def correct_heading(ADCP, glider, options):
    if options["correctADCPHeading"]:
        if "Heading_old" in ADCP:
            ADCP["Heading"] = ("time", ADCP["Heading_old"].values)
            _log.info("Resetting to original heading")

        ADCP["Heading_old"] = ("time", ADCP["Heading"].values)
        ADCP["Heading"] = (
            _heading_correction(ADCP, glider, options) + ADCP["declination"]
        )
        plog("Corrected heading and accounted for declination")
    else:
        plog("Uncorrected heading and declination NOT added.")
    return ADCP


def soundspeed_correction(ADCP):
    if "NoSal_SpeedOfSound" not in ADCP:
        ADCP = ADCP.rename({"SpeedOfSound": "NoSal_SpeedOfSound"})
        ADCP = ADCP.rename({"glider_soundspeed": "SpeedOfSound"})
        for beam in ["1", "2", "3", "4"]:
            # V_new = V_old * (c_new/c_old)
            ADCP["VelocityBeam" + beam] = ADCP["VelocityBeam" + beam] * (
                ADCP["SpeedOfSound"] / ADCP["NoSal_SpeedOfSound"]
            )
            plog("Corrected beam " + beam + " velocity for sound speed.")
    else:
        plog("Speed of sound correction already applied")
    return ADCP


def remove_outliers(ADCP, options):
    """
    Function discards velocity estimates outside of user specified bounds in `options`.
    "ADCP_correlationThreshold": minimum % correlation within each ensemble
    "ADCP_amplitudeThreshold": minimum return amplitude
    "ADCP_velocityThreshold": maximum relative velocity

    :param ADCP: xr.DataSet of ADCP data
    :param options: options dictionary
    :return: ADCP with outliers removed
    """

    # From Tanaka:
    # the velocity was 0.5 m s-1 or less,
    # the echo intensity (or amplitude) was 75 dB or less,
    # the percent-good was 80 or greater,
    # the signal-to-noise ratio (SNR) was 20 or greater.
    def prct(ind):
        out = np.count_nonzero(ind) / n * 100
        out = np.round(out * 10) / 10
        return str(out)

    for beam in ["1", "2", "3", "4"]:
        C = ADCP["CorrelationBeam" + beam].values.copy()
        n = len(C.flatten())
        ind = C < options["ADCP_correlationThreshold"]
        C[ind] = np.nan
        C[np.isfinite(C)] = 1
        plog("Beam " + beam + " correlation: " + prct(ind) + "% removed")

        A = ADCP["AmplitudeBeam" + beam].values.copy()
        ind = A > options["ADCP_amplitudeThreshold"]
        A[ind] = np.nan
        # A[A < 40] = np.nan
        A[np.isfinite(A)] = 1
        plog("Beam " + beam + " amplitude: " + prct(ind) + "% removed")

        V = ADCP["VelocityBeam" + beam].values.copy()
        ind = np.abs(V) > options["ADCP_velocityThreshold"]
        V[ind] = np.nan
        V[np.isfinite(V)] = 1
        plog("Beam " + beam + " velocity: " + prct(ind) + "% removed")

        ADCP["VelocityBeam" + beam] = ADCP["VelocityBeam" + beam] * C * A * V

    return ADCP


# Defining coordinate transform functions for the 4 beam ADCP configuration
def quad_beam2xyzz_mat():
    M11 = 0.6782
    M12 = 0.0000
    M13 = -0.6782
    M14 = 0.0000
    M21 = 0.0000
    M22 = -1.1831
    M23 = 0.0000
    M24 = 1.1831
    M31 = 0.7400
    M32 = 0.0000
    M33 = 0.7400
    M34 = 0.0000
    M41 = 0.0000
    M42 = 0.5518
    M43 = 0.0000
    M44 = 0.5518
    T = np.array(
        [
            [M11, M12, M13, M14],
            [M21, M22, M23, M24],
            [M31, M32, M33, M34],
            [M41, M42, M43, M44],
        ]
    )
    # if not top_mounted:
    #     T[1:,:] = -T[1:,:]
    return T


def quad_beam2xyzz(B1, B2, B3, B4):
    T = quad_beam2xyzz_mat()
    try:
        r, c = np.shape(B1.values)
    except ValueError:
        c = 1
        r = len(B1.values)
    V = np.array(
        [
            B1.values.flatten(),
            B2.values.flatten(),
            B3.values.flatten(),
            B4.values.flatten(),
        ]
    )
    XYZZ = V * 0
    for col in np.arange(V.shape[1]):
        XYZZ[:, col] = T @ V[:, col]
    return (
        np.reshape(XYZZ[0, :], [r, c]),
        np.reshape(XYZZ[1, :], [r, c]),
        np.reshape(XYZZ[2, :], [r, c]),
        np.reshape(XYZZ[3, :], [r, c]),
    )


def quad_xyzz2beam(X, Y, Z, ZZ):
    T = np.linalg.inv(quad_beam2xyzz_mat())
    r, c = np.shape(X.values)
    XYZZ = np.array(
        [
            X.values.flatten(),
            Y.values.flatten(),
            Z.values.flatten(),
            ZZ.values.flatten(),
        ]
    )
    V = XYZZ * 0
    for col in np.arange(XYZZ.shape[1]):
        V[:, col] = T @ XYZZ[:, col]
    return (
        np.reshape(V[0, :], [r, c]),
        np.reshape(V[1, :], [r, c]),
        np.reshape(V[2, :], [r, c]),
        np.reshape(V[3, :], [r, c]),
    )


def do_xyzz2beam(ADCP):
    V = quad_xyzz2beam(ADCP["X4"], ADCP["Y4"], ADCP["Z4"], ADCP["ZZ4"])
    ADCP["VelocityBeam1"] = (["time", "bin"], V[0])
    ADCP["VelocityBeam2"] = (["time", "bin"], V[1])
    ADCP["VelocityBeam3"] = (["time", "bin"], V[2])
    ADCP["VelocityBeam4"] = (["time", "bin"], V[3])
    return ADCP


def do_beam2xyzz(ADCP):
    XYZZ = quad_beam2xyzz(
        ADCP["VelocityBeam1"],
        ADCP["VelocityBeam2"],
        ADCP["VelocityBeam3"],
        ADCP["VelocityBeam4"],
    )
    ADCP["X4"] = (["time", "bin"], XYZZ[0])
    ADCP["Y4"] = (["time", "bin"], XYZZ[1])
    ADCP["Z4"] = (["time", "bin"], XYZZ[2])
    ADCP["ZZ4"] = (["time", "bin"], XYZZ[3])
    return ADCP


# def _shear_correction(ADCP, var, correct=True):
#     def get_correction_array(row):
#         ### SEPARATION CRITERIA
#         spd_thr_water = np.sqrt(
#             ADCP["X4"] ** 2 + ADCP["Y4"] ** 2 + ((ADCP["Z4"] + ADCP["ZZ4"]) / 2) ** 2
#         )

#         spd = spd_thr_water.values[:, 0]  # .values[:,5] #.mean('bin')

#         range_interval = 0.01
#         full_range = np.arange(
#             0.2, 0.3, range_interval
#         )  # OMAN BEST: np.arange(0.2,0.3,range_interval)

#         variable = ADCP[var]
#         ref = np.cumsum(
#             np.insert(
#                 variable.isel(time=(spd > full_range[0]) & (spd < full_range[-1]))
#                 .diff("bin")
#                 .mean("time")
#                 .values,
#                 0,
#                 0,
#             )
#         )
#         ref = ref - np.nanmean(ref)

#         return ref

#     if correct:
#         ADCP[var] = ADCP[var] - get_correction_array(0)
#         _ = get_correction_array(3)
#         plog("Corrected " + var)
#     else:
#         _ = get_correction_array(0)

#     return None


# def correct_shear(ADCP, options):
#     ADCP = do_beam2xyzz(ADCP)
#     if options["correctZZshear"]:
#         ADCP = do_beam2xyzz(ADCP)
#         _shear_correction(ADCP, "ZZ4")
#         ADCP = do_xyzz2beam(ADCP)
#     if options["correctZshear"]:
#         ADCP = do_beam2xyzz(ADCP)
#         _shear_correction(ADCP, "Z4")
#         ADCP = do_xyzz2beam(ADCP)
#     if options["correctYshear"]:
#         ADCP = do_beam2xyzz(ADCP)
#         _shear_correction(ADCP, "Y4")
#         ADCP = do_xyzz2beam(ADCP)
#     if options["correctXshear"]:
#         ADCP = do_beam2xyzz(ADCP)
#         _shear_correction(ADCP, "X4")
#         ADCP = do_xyzz2beam(ADCP)

#     return ADCP


# def correct_backscatter(ADCP, glider):
#     # https://www.sciencedirect.com/science/article/pii/S0278434304002171#fig2
#     # 2.1. Intensity correction

#     # For a single acoustic target, intensity falls proportional to distance squared.
#     # For volume scattering, on the other hand, approximately half of the loss is
#     # cancelled by the increasing number of (zooplankton) targets ensonified by the
#     # spreading beam. Defining  and  as the corrected and measured ADCP beam intensity,
#     # respectively, the correction for volume scattering is given by

#     # I(c) = I(m) + 2( 10 * log(r) + alpha * r )

#     # with r the distance from the transducer to each of the depth bins, alpha the
#     # absorption, or attenuation, coefficient, and where the extra factor of 2 accounts for
#     # the return trip the sound pulse must make. For single target scattering, the correction
#     # recommended by RDI, the 10 in Eq. (1) is replaced with 20; this correction has been used,
#     # for example, by Thomson and Allen (2000), who examined macrozooplankton migration and
#     # fish abundance. The rate of absorption in seawater is a function of pressure, temperature,
#     # salinity, pH, and sound frequency; using the empirical equation derived by Francois and
#     # Garrison (1982) for a depth of 100 m, , and  for the ADCP frequency of 307 kHz.

#     def francoisgarrison(freq=None, S=None, T=None, pH=8.1, z=None):
#         # Underwater acoustics: calculation of the absorption coefficient
#         #
#         # SYNTAX: a = francoisgarrison(freq,T,S,pH,z)
#         # Units: [freq] = kHz, [T] = Celsius, [S] = ppt, [z] = m, [a] = dB/m
#         #
#         # Reference:
#         # Sound absorption based on ocean measurements. Part II: Boric acid contribution
#         # and equation for total absorption
#         # R.E. Francois and G. R. Garrison
#         # J. Acoust. Soc. Am. 72(6), December 1982

#         # *******************************************************************************
#         # Faro, Ter Set  7 23:46:58 WEST 2021
#         # Written by Orlando Camargo Rodriguez
#         # *******************************************************************************
#         c = 1412.0 + 3.21 * T + 1.19 * S + 0.0167 * z
#         theta = 273.0 + T
#         fxf = freq**2
#         f1 = 2.8 * np.sqrt(S / 35.0) * 10 ** (4.0 - 1245.0 / theta)
#         f2 = 8.17 * 10 ** (8.0 - 1990.0 / theta) / (1.0 + 0.0018 * (S - 35.0))
#         A1 = 8.86 / c * 10 ** (0.78 * pH - 5)
#         A2 = 21.44 * S / c * (1.0 + 0.025 * T)
#         A3 = 3.964e-4 - 1.146e-5 * T + 1.45e-7 * T**2 - 6.5e-10 * T**3
#         ind = T <= 20
#         A3[ind] = (
#             4.937e-4 - 2.59e-5 * T[ind] + 9.11e-7 * T[ind] ** 2 - 1.50e-8 * T[ind] ** 3
#         )
#         P1 = 1.0
#         P2 = 1 - 1.35e-4 * z + 6.2e-9 * z**2
#         P3 = 1 - 3.83e-5 * z + 4.9e-10 * z**2
#         alpha = (
#             A1 * P1 * f1 * fxf / (f1**2 + fxf)
#             + A2 * P2 * f2 * fxf / (f2**2 + fxf)
#             + A3 * P3 * fxf
#         )
#         return alpha / 1000

#     a = np.cos(np.deg2rad(47.4))  # Beam 1 and 3 angle from Z 47.5
#     b = np.cos(np.deg2rad(25))  # Beam 2 and 4 angle from Z 25

#     ADCP["BeamRange1"] = ("bin", ADCP["Velocity Range"].values / a)
#     ADCP["BeamRange2"] = ("bin", ADCP["Velocity Range"].values / b)
#     ADCP["BeamRange3"] = ("bin", ADCP["Velocity Range"].values / a)
#     ADCP["BeamRange4"] = ("bin", ADCP["Velocity Range"].values / b)

#     ADCP["AcousticAttenuation"] = (
#         "time",
#         francoisgarrison(
#             freq=1000,
#             S=interp(
#                 glider.time.values.astype("float"),
#                 glider.salinity.interpolate("index").bfill().values,
#                 ADCP.time.values.astype("float"),
#             ),
#             T=interp(
#                 glider.time.values.astype("float"),
#                 glider.temperature.interpolate("index").bfill().values,
#                 ADCP.time.values.astype("float"),
#             ),
#             pH=8.1,
#             z=ADCP["Depth"].values,
#         ),
#     )

#     for beam in ["1", "2", "3", "4"]:
#         ADCP["AmplitudeNew" + beam] = (
#             ["time", "bin"],
#             (
#                 ADCP["AmplitudeBeam" + beam]
#                 + 2 * ADCP["AcousticAttenuation"] * ADCP["BeamRange" + beam]
#             ).values,
#         )  # + 20 int(particule attenuation,range)
#     return ADCP


def regridADCPdata(ADCP, options, depth_offsets=None):
    ## This is to avoid shear smearing because of tilted ADCP
    bin_size = ADCP.attrs["avg_cellSize"]
    blanking_distance = ADCP.attrs["avg_blankingDistance"]

    ## FN: calculate isobar offsets relative to glider depth for each ping
    def calc_ideal_depth_offsets(bin_size, blanking_distance):
        if options["top_mounted"]:
            direction = 1
        else:
            direction = -1
        threshold = options["ADCP_regrid_correlation_threshold"]
        means = [
            np.nanmean(ADCP["CorrelationBeam" + str(x)], axis=0) for x in [1, 2, 3, 4]
        ]

        max_bin = np.argmin(
            abs(np.nanmean([means[x] for x in range(4)], axis=0) - threshold)
        )
        max_distance = blanking_distance + max_bin * bin_size + 0.5 * bin_size
        return np.arange(0, max_distance + bin_size, bin_size / 2) * direction

    ## Calculate desired interpolation depth offsets
    if depth_offsets is None:
        depth_offsets = calc_ideal_depth_offsets(bin_size, blanking_distance)

    _log.info("Using the following depth offsets:")
    _log.info(depth_offsets)
    _log.info(" ")
    _log.info("Running gridding on all 4 beams:")

    ## Extract to np array for speed

    for beam in ["1", "2", "3", "4"]:
        plog("Calculating beam " + beam)

        def interp1d_np(x, y):
            _gd = np.isfinite(y)
            if np.count_nonzero(_gd) > 1:
                xi = interp1d(x[_gd], y[_gd], bounds_error=False, fill_value=np.nan)(
                    depth_offsets
                )
            else:
                xi = depth_offsets * np.nan
            return xi

        ADCP.load()
        ADCP["V" + beam] = xr.apply_ufunc(
            interp1d_np,
            ADCP["Depth"] - ADCP["D" + beam],
            ADCP["VelocityBeam" + beam],
            input_core_dims=[["bin"], ["bin"]],
            output_core_dims=[["gridded_bin"]],
            exclude_dims={"bin"},
            vectorize=True,
            output_sizes={"gridded_bin": len(depth_offsets)},
        )

    ADCP = ADCP.assign_coords({"depth_offset": (["gridded_bin"], depth_offsets)})
    ADCP = ADCP.assign_coords(
        {
            "bin_depth": (
                ["time", "gridded_bin"],
                np.tile(
                    ADCP["Depth"].values.astype("float"), (len(ADCP.gridded_bin), 1)
                ).T
                - np.tile(depth_offsets, (len(ADCP.time), 1)),
            )
        }
    )

    return ADCP


def calcXYZfrom3beam(ADCP, options):
    """
    Coordinate transform that converts velocity estimates along beams (V1-4) to glider relative velocities X, Y , Z

    :param ADCP: xr.DataSet of ADCP data
    :param options: options dictionary
    :return: ADCP with X, Y, Z velocities
    """

    def sin(x):
        return np.sin(np.deg2rad(x))

    def cos(x):
        return np.cos(np.deg2rad(x))

    tf = 47.5  # theta front - Beam 1 and 3 angle from Z
    ts = 25  # theta side - Beam 2 and 4 angle from Z

    c = 1  # for convex transducer head

    def a(t):
        return 1 / (2 * sin(t))

    def b(t):
        return 1 / (4 * cos(t))

    V1 = ADCP["V1"].values
    V2 = ADCP["V2"].values
    V3 = ADCP["V3"].values
    V4 = ADCP["V4"].values

    # Because of duff beam part of the time, we must use the 3 beam solution that relies on error estimation
    # error = d(tf)*V1 + d(tf)*V3 - d(ts)*V2 - d(ts)*V4
    # if error == 0 then...
    # V1 = (-d(tf)*V3 + d(ts)*V2 + d(ts)*V4) / d(tf)
    # V3 = (-d(tf)*V1 + d(ts)*V2 + d(ts)*V4) / d(tf)

    # fZ(V1,V3) = fZ(V2,V4) ...
    # 2*b(tf)*V1 + 2*b(tf)*V3 = 2*b(ts)*V2 + 2*b(ts)*V4
    # 0.7400*V1 + 0.7400*V3 = 0.5518*V2 + 0.5518*V4

    # replaced_by = lambda good_beam : (-d(tf)*good_beam + d(ts)*V2 + d(ts)*V4) / d(tf)
    # replaced_by = lambda good_beam :
    # (-(0.6782/np.sqrt(2))*good_beam + (1.1831/np.sqrt(2))*V2 + (1.1831/np.sqrt(2))*V4) / (0.6782/np.sqrt(2))

    # replaced_by = lambda good_beam : (0.5518*V2 + 0.5518*V4 - 0.7400*good_beam)/0.7400
    def replaced_by(good_beam):
        return (2 * b(ts) * V2 + 2 * b(ts) * V4 - 2 * b(tf) * good_beam) / (2 * b(tf))

    upcasts = ADCP["Pitch"] > 0
    downcasts = ~upcasts

    if options["top_mounted"]:
        _log.info("Assuming ADCP is top mounted")
        V1[downcasts, :] = replaced_by(V3)[
            downcasts, :
        ]  # Use aft beam on downcasts when upward facing
        V3[upcasts, :] = replaced_by(V1)[
            upcasts, :
        ]  # Use fore beam on upcasts when upward facing
    else:
        _log.info("Assuming ADCP is bottom mounted")
        V1[upcasts, :] = replaced_by(V3)[
            upcasts, :
        ]  # Use aft beam on upcasts when downward facing
        V3[downcasts, :] = replaced_by(V1)[
            downcasts, :
        ]  # Use fore beam on downcasts when downward facing

    X = c * a(tf) * V1 - c * a(tf) * V3
    Y = -c * a(ts) * V2 + c * a(ts) * V4  # TODO: sign uncertainty here
    Z = 2 * b(ts) * V2 + 2 * b(ts) * V4

    ADCP["X"] = (["time", "gridded_bin"], X)
    ADCP["Y"] = (["time", "gridded_bin"], Y)
    ADCP["Z"] = (["time", "gridded_bin"], Z)
    plog("Calculating X,Y,Z from isobaric 3-beam measurements.")

    return ADCP


def calcENUfromXYZ(ADCP, options):
    """
    Coordinate transform that converts velocity estimates relative to the glider (X, Y Z) into the earth relative
    reference frame east north up (ENU)

    :param ADCP: xr.DataSet of ADCP data
    :param options: options dictionary
    :return: ADCP with E, N, U velocities
    """

    def M_xyz2enu(heading, pitch, roll):
        hh = np.pi * (heading - 90) / 180
        pp = np.pi * pitch / 180
        rr = np.pi * roll / 180

        GPT = [
            [
                np.cos(hh) * np.cos(pp),
                -np.cos(hh) * np.sin(pp) * np.sin(rr) + np.sin(hh) * np.cos(rr),
                -np.cos(hh) * np.sin(pp) * np.cos(rr) - np.sin(hh) * np.sin(rr),
            ],
            [
                -np.sin(hh) * np.cos(pp),
                np.sin(hh) * np.sin(pp) * np.sin(rr) + np.cos(hh) * np.cos(rr),
                np.sin(hh) * np.sin(pp) * np.cos(rr) - np.cos(hh) * np.sin(rr),
            ],
            [np.sin(pp), np.cos(pp) * np.sin(rr), np.cos(pp) * np.cos(rr)],
        ]
        return GPT

    if options["top_mounted"]:
        direction = 1
    else:
        direction = -1

    H = ADCP["Heading"]
    P = ADCP["Pitch"]
    R = ADCP["Roll"]
    M = M_xyz2enu(H, P, R)

    E = (
        M[0][0] * ADCP["X"]
        + M[0][1] * ADCP["Y"] * direction
        + M[0][2] * ADCP["Z"] * direction
    )
    N = (
        M[1][0] * ADCP["X"]
        + M[1][1] * ADCP["Y"] * direction
        + M[1][2] * ADCP["Z"] * direction
    )
    U = (
        M[2][0] * ADCP["X"]
        + M[2][1] * ADCP["Y"] * direction
        + M[2][2] * ADCP["Z"] * direction
    )

    ADCP["E"] = E
    ADCP["N"] = N
    ADCP["U"] = U

    plog("Converted from XYZ to ENU")

    profiles = ADCP["profile_number"].values
    profile_center = np.nanmean(profiles)
    profile_range = np.nanmax(profiles) - np.nanmin(profiles)
    _gd = (ADCP["Pressure"].values > 5) & (
        np.abs(ADCP["profile_number"] - profile_center) < profile_range * 0.1
    )
    ADCP["Sh_E"] = (
        ["time", "gridded_bin"],
        ADCP["E"].differentiate("gridded_bin").values,
    )
    ADCP["Sh_N"] = (
        ["time", "gridded_bin"],
        ADCP["N"].differentiate("gridded_bin").values,
    )
    ADCP["Sh_U"] = (
        ["time", "gridded_bin"],
        ADCP["U"].differentiate("gridded_bin").values,
    )
    x = np.arange(0, np.shape(ADCP.Sh_N.values)[0], 1)
    SHEm, XI, YI = grid2d(
        np.tile(ADCP.profile_number.values, (len(ADCP.gridded_bin), 1))
        .T[x, :]
        .flatten(),
        ADCP.bin_depth.values[x, :].flatten(),
        ADCP.Sh_N.values[x, :].flatten(),
        xi=1,
        yi=3,
        fn="mean",
    )

    isnan = np.isnan(SHEm)
    SHEm = np.nancumsum(SHEm[::-1, :], 0)[::-1, :]
    SHEm[isnan] = np.nan
    return ADCP


# def get_DAC(ADCP, glider):
#     """
#     Estimate dive averaged horizontal currents for the glider. This function requires estimates of the glider's
#     horizontal and vertical speed through water

#     :param ADCP: xr.DataSet of ADCP data
#     :param glider: pd.DataFrame of glider data
#     :param options: options dictionary
#     :return: glider DataFrame with estimates of dive average current
#     """

#     ## Calculate full x-y dead reckoning during each dive
#     def reset_transport_at_GPS(arr):
#         def ffill(arr):
#             return pd.DataFrame(arr).ffill().values.flatten()

#         ref = np.zeros(np.shape(arr)) * np.nan
#         ref[_gps] = arr[_gps]
#         return arr - ffill(ref)

#     _gps = (glider.dead_reckoning.values < 1) & (glider.nav_resource.values == 116)

#     t = glider.date_float.values * 1e-9
#     heading = interp(
#         ADCP["time"].values.astype("float"),
#         ADCP["Heading"].values,
#         glider.date_float.values,
#     )
#     vg_e = np.nan_to_num(glider["speed_horz"] * np.sin(heading * np.pi / 180))
#     vg_n = np.nan_to_num(glider["speed_horz"] * np.cos(heading * np.pi / 180))

#     glider["speed_e"] = vg_e
#     glider["speed_n"] = vg_n

#     de = np.cumsum(np.append(0, vg_e[1:] * np.diff(t)))
#     dn = np.cumsum(np.append(0, vg_n[1:] * np.diff(t)))

#     de = reset_transport_at_GPS(de)
#     dn = reset_transport_at_GPS(dn)

#     ## Calculate on per dive basis
#     dnum = np.unique(glider.dive_number.values)
#     sidx = np.zeros(np.shape(dnum)) * np.nan
#     didx = np.zeros(np.shape(dnum)) * np.nan

#     for idx, dx in enumerate(dnum):
#         try:
#             sidx[idx] = np.flatnonzero((glider.dive_number.values == dx) & _gps)[0]
#             didx[idx] = np.flatnonzero((glider.dive_number.values == dx) & _gps)[-1]
#         except IndexError:
#             continue

#     _gd = np.isfinite(sidx + didx + dnum)
#     dnum = dnum[_gd]
#     sidx = sidx[_gd]
#     didx = didx[_gd]

#     sidx = sidx.astype(int)
#     didx = didx.astype(int)

#     surf_lat = glider.latitude.values[sidx]
#     surf_lon = glider.longitude.values[sidx]
#     surf_time = t[sidx]

#     dive_lat = glider.latitude.values[didx]
#     dive_lon = glider.longitude.values[didx]
#     dive_time = t[didx]

#     dr_e = np.zeros(np.shape(dnum)) * np.nan
#     dr_n = np.zeros(np.shape(dnum)) * np.nan
#     gps_e = np.zeros(np.shape(dnum)) * np.nan
#     gps_n = np.zeros(np.shape(dnum)) * np.nan
#     dt = np.zeros(np.shape(dnum)) * np.nan
#     meant = np.zeros(np.shape(dnum)) * np.nan

#     def lon2m(x, y):
#         return gsw.distance([x, x + 1], [y, y])

#     def lat2m(x, y):
#         return gsw.distance([x, x], [y, y + 1])

#     for idx, dx in enumerate(dnum):
#         try:
#             dr_e[idx] = de[sidx[idx + 1] - 1]
#             dr_n[idx] = dn[sidx[idx + 1] - 1]
#             gps_e[idx] = (surf_lon[idx + 1] - dive_lon[idx]) * lon2m(
#                 dive_lon[idx], dive_lat[idx]
#             )[0]
#             gps_n[idx] = (surf_lat[idx + 1] - dive_lat[idx]) * lat2m(
#                 dive_lon[idx], dive_lat[idx]
#             )[0]
#             dt[idx] = surf_time[idx + 1] - dive_time[idx]
#             meant[idx] = (surf_time[idx + 1] + dive_time[idx]) / 2
#         except IndexError:
#             _log.info("No final GPS for dive " + str(dx))

#     glider["DAC_E"] = interp(meant, (gps_e - dr_e) / dt, t)
#     glider["DAC_N"] = interp(meant, (gps_n - dr_n) / dt, t)

#     glider["DAC_E"] = glider["DAC_E"].bfill().ffill()
#     glider["DAC_N"] = glider["DAC_N"].bfill().ffill()

#     return glider


# def getSurfaceDrift(glider):
#     _gps = (glider.dead_reckoning.values < 1) & (glider.nav_resource.values == 116)

#     def lon2m(x, y):
#         return gsw.distance([x, x + 1], [y, y])

#     def lat2m(x, y):
#         return gsw.distance([x, x], [y, y + 1])

#     dnum = glider.dive_number.values[_gps]

#     lons = glider.longitude.values[_gps]
#     lats = glider.latitude.values[_gps]

#     dlons = np.gradient(lons)
#     dlats = np.gradient(lats)

#     for idx in range(len(lons)):
#         dlons[idx] = dlons[idx] * lon2m(lons[idx], lats[idx])[0]
#         dlats[idx] = dlats[idx] * lat2m(lons[idx], lats[idx])[0]

#     times = glider.time.values.astype("float")[_gps] / 10**9
#     dtimes = np.gradient(times)

#     dE = np.full(int(np.nanmax(glider.dive_number)), np.nan)
#     dN = np.full(int(np.nanmax(glider.dive_number)), np.nan)
#     dT = np.full(int(np.nanmax(glider.dive_number)), np.nan)

#     for idx in range(len(dE)):
#         _gd = (dtimes < 21) & (dnum == idx + 1)
#         dE[idx] = np.nanmedian(dlons[_gd] / dtimes[_gd])
#         dN[idx] = np.nanmedian(dlats[_gd] / dtimes[_gd])
#         dT[idx] = np.nanmean(times[_gd])

#     dT = dT * 10**9
#     return dE, dN, dT


# def bottom_track(ADCP, adcp_file_path, options):
#     if options["top_mounted"]:
#         plog("ERROR: ADCP is top mounted. Not processing bottom track data")
#         return ADCP

#     def sin(x):
#         return np.sin(np.deg2rad(x))

#     def cos(x):
#         return np.cos(np.deg2rad(x))

#     # BT gives us speed of the glider
#     # If we subtract BT velocity from XYZ
#     # then we get speed of water
#     BT = xr.open_mfdataset(adcp_file_path, group="Data/AverageBT")
#     BT = BT.drop_vars(['MatlabTimeStamp'])
#     BT = BT.where(BT.time < ADCP.time.values[-1]).dropna(dim="time", how="all")

#     thresh = 12

#     ind = (
#         (BT["VelocityBeam1"] > -2)
#         & (BT["VelocityBeam2"] > -2)
#         & (BT["VelocityBeam4"] > -2)
#     )
#     ind2 = (
#         (BT["FOMBeam1"] < thresh)
#         & (BT["FOMBeam2"] < thresh)
#         & (BT["FOMBeam4"] < thresh)
#     )
#     BT = BT.isel(time=ind.values & ind2.values)

#     full_time = ADCP["time"].values.astype("float")
#     BT_time = BT["time"].values.astype("float")
#     matching = []
#     for idx in range(len(BT_time)):
#         matching.append(np.argmin(np.abs(BT_time[idx] - full_time)))

#     ADCP_profile = ADCP["profile_number"].values.copy()
#     ADCP_depth = ADCP["Pressure"].values.copy()
#     for idx in range(np.nanmax(ADCP_profile).astype(int)):
#         _gd = ADCP_profile == idx
#         if np.count_nonzero(_gd) == 0:
#             _log.info("Profile " + str(idx) + " was empty")
#         else:
#             ADCP_depth[_gd] = np.nanmax(ADCP_depth[_gd])
#     ADCP_depth = ADCP_depth[matching]

#     _gd = np.abs(ADCP_depth - BT["Pressure"]).values < 15
#     BT = BT.isel(time=_gd)
#     full_time = ADCP["time"].values.astype("float")
#     BT_time = BT["time"].values.astype("float")
#     matching = []
#     for idx in range(len(BT_time)):
#         matching.append(np.argmin(np.abs(BT_time[idx] - full_time)))

#     C_old = BT["SpeedOfSound"].values
#     C_new = ADCP["SpeedOfSound"].isel(time=matching).values

#     a = 47.5  # Beam 1 and 3 angle from Z
#     b = 25  # Beam 2 and 4 angle from Z
#     xyz2beam_fore = np.array(
#         [[sin(a), 0, cos(a)], [0, -sin(b), cos(b)], [0, sin(b), cos(b)]]
#     )
#     beam2xyz_fore = np.linalg.inv(xyz2beam_fore)

#     BT_X4, BT_Y4, BT_Z4 = beam2xyz_fore @ np.array(
#         [
#             BT["VelocityBeam1"] * (C_new / C_old),
#             BT["VelocityBeam2"] * (C_new / C_old),
#             BT["VelocityBeam4"] * (C_new / C_old),
#         ]
#     )

#     def M_xyz2enu(heading, pitch, roll):
#         hh = np.pi * (heading - 90) / 180
#         pp = np.pi * pitch / 180
#         rr = np.pi * roll / 180

#         _H = np.array(
#             [[np.cos(hh), np.sin(hh), 0], [-np.sin(hh), np.cos(hh), 0], [0, 0, 1]]
#         )
#         _P = np.array(
#             [[np.cos(pp), 0, -np.sin(pp)], [0, 1, 0], [np.sin(pp), 0, np.cos(pp)]]
#         )
#         _R = np.array(
#             [[1, 0, 0], [0, np.cos(rr), -np.sin(rr)], [0, np.sin(rr), np.cos(rr)]]
#         )

#         _M = _H @ _P @ _R
#         return _M

#     H = BT["Heading"].values
#     P = BT["Pitch"].values
#     R = BT["Roll"].values

#     BT_E = np.full_like(H, np.nan)
#     BT_N = np.full_like(H, np.nan)
#     BT_U = np.full_like(H, np.nan)

#     if options["top_mounted"]:
#         direction = 1
#     else:
#         direction = -1

#     n = len(BT_X4)
#     for i in range(n):
#         BT_E[i], BT_N[i], BT_U[i] = M_xyz2enu(H[i], P[i], R[i]) @ [
#             BT_X4[i],
#             BT_Y4[i] * direction,
#             BT_Z4[i] * direction,
#         ]

#     bt_e = np.full_like(full_time, np.nan)
#     bt_e[matching] = BT_E
#     bt_n = np.full_like(full_time, np.nan)
#     bt_n[matching] = BT_N
#     bt_u = np.full_like(full_time, np.nan)
#     bt_u[matching] = BT_U

#     ADCP["BT_E"] = ("time", bt_e)
#     ADCP["BT_N"] = ("time", bt_n)
#     ADCP["BT_U"] = ("time", bt_u)

#     return ADCP


def grid_shear_data(glider): # TODO name and docstring are way off - this is just defining axis
    """
    Grid ENU velocities into standardised vertical bins.

    :param ADCP: xr.DataSet of ADCP data
    :param glider: pd.DataFrame of glider data
    :param options: options dictionary
    :return: ADCP with vertically gridded ENU velocities
    """
    yaxis = np.arange(0, np.nanmax(np.ceil(glider.pressure.values)), y_res)
    xaxis = glider.date_float.groupby(glider.profile_number).agg("mean").index
    taxis = pd.to_datetime(
        glider.date_float.groupby(glider.profile_number).agg("mean").values
    )
    days = np.unique(glider.time.round("D"))
    return xaxis, yaxis, taxis, days


# def reference_shear(ADCP, glider, xaxis, yaxis):
#     """
#     Reference the estimates of vertical shear of horizontal velocity using a per-profile average velocity

#     :param ADCP: xr.DataSet of ADCP data
#     :param options: options dictionary
#     :param dE: eastward surface drift velocity
#     :param dN: northward surface drift velocity
#     :param dT: timestamps surface drift velocity
#     :param xaxis: profile numbers of gridded velocity data
#     :param yaxis: depth bins of gridded velocity data
#     :param taxis: time of gridded velocity data
#     :param options: options dict for processing
#     :return:  xr.DataSet of referenced gridded N and E velocities
#     """
#     out = {}

#     var = ["E", "N"]

#     for pstep in range(len(var)):
#         letter = var[pstep]
#         # Grid shear to average out sensor + zooplankton noise
#         Sh, XI, YI = grid2d(
#             np.tile(ADCP.profile_number.values, (len(ADCP.gridded_bin), 1)).T.flatten(),
#             ADCP.bin_depth.values.flatten(),
#             ADCP["Sh_" + letter].values.flatten(),
#             xi=xaxis,
#             yi=yaxis,
#             fn="mean",
#         )

#         # Integrate shear vertically
#         _bd = ~np.isfinite(
#             Sh
#         )  # Preserve what are originally NaN values to recover later as need conversion to 0 for cumsum-
#         Sh = np.nan_to_num(Sh)  # Replace NaNs with 0 for cumsum
#         V = (
#             np.cumsum(Sh, axis=0) * y_res
#         )  # Cumulative sum of shear to recover velocity profile
#         V[_bd] = np.nan  # Return NaNs to their rightful place.
#         V = V - np.tile(
#             np.nanmean(V, axis=0), (np.shape(V)[0], 1)
#         )  # Make mean of baroclinic profiles equal to 0

#         # Grid DAC
#         DAC, XI, YI = grid2d(
#             glider.profile_number.values,
#             glider.pressure.values,
#             glider["DAC_" + letter].values,
#             xi=xaxis,
#             yi=yaxis,
#             fn="mean",
#         )

#         # Grid vertical speed
#         dPdz, XI, YI = grid2d(
#             glider.profile_number.values,
#             glider.pressure.values,
#             glider["speed_vert"].values,
#             xi=xaxis,
#             yi=yaxis,
#             fn="mean",
#         )

#         # Seconds spent in each depth bin, to weight referencing
#         SpB = y_res / dPdz
#         SpB[np.isinf(SpB)] = 0
#         strictness = 1
#         SpB_std = np.nanstd(SpB.flatten())
#         SpB[np.abs(SpB) > (strictness * SpB_std)] = strictness * SpB_std

#         # Baroclinic velocity, weighted by depth residence time, should be equal to DAC
#         # So the reference to add to a baroclinic profile of mean = 0 is the DAC - the weighted baroclinic velocity.
#         Ref = np.nanmean(DAC, axis=0) - np.nansum(V * SpB, axis=0) / np.nansum(
#             SpB, axis=0
#         )

#         # Now we reference the velocity
#         V = V + np.tile(Ref, (np.shape(V)[0], 1))
#         out["ADCP_" + letter] = V
#     return out


# def _grid_glider_data(glider, out, xaxis, yaxis):
#     exclude_from_grid = [
#         "AD2CP_ALT",
#         "AD2CP_HEADING",
#         "AD2CP_PITCH",
#         "AD2CP_PRESSURE",
#         "AD2CP_ROLL",
#         "AROD_FT_DO",
#         "AROD_FT_TEMP",
#         "Altitude",
#         "AngCmd",
#         "AngPos",
#         "BallastCmd",
#         "BallastPos",
#         "dead_reckoning",
#         "declination",
#         "Depth",
#         "DesiredH",
#         "LinCmd",
#         "LinPos",
#         "NAV_DEPTH",
#         "NAV_LATITUDE",
#         "NAV_LONGITUDE",
#         "nav_resource",
#         "NavState",
#         "Pa",
#         "Pitch",
#         "Roll",
#         "Heading",
#         "SecurityLevel",
#         "Temperature",
#         "time",
#         "Unnamed: 22",
#         "Unnamed: 28",
#         "Voltage",
#         "missionNum",
#         "Lat",
#         "Lon",
#     ]

#     variables = glider.columns
#     variables = [x for x in variables if x not in exclude_from_grid]

#     def grid(name):
#         return grid2d(
#             glider.profile_number.values,
#             glider.pressure.values,
#             glider[name].values,
#             xi=xaxis,
#             yi=yaxis,
#             fn="mean",
#         )[0]

#     for varname in variables:
#         try:
#             out[varname] = grid(varname)
#         except IndexError:
#             _log.info('Variable "' + varname + '" failed to grid.')

#     return out


# def grid_data(ADCP, glider, out, xaxis, yaxis):
#     ADCP_pnum = np.tile(ADCP.profile_number, (len(ADCP.gridded_bin), 1)).T
#     out["Sh_E"] = grid2d(
#         ADCP_pnum.flatten(),
#         ADCP.bin_depth.values.flatten(),
#         ADCP.Sh_E.values.flatten(),
#         xi=xaxis,
#         yi=yaxis,
#         fn="mean",
#     )[0]
#     out["Sh_N"] = grid2d(
#         ADCP_pnum.flatten(),
#         ADCP.bin_depth.values.flatten(),
#         ADCP.Sh_N.values.flatten(),
#         xi=xaxis,
#         yi=yaxis,
#         fn="mean",
#     )[0]
#     out["Sh_U"] = grid2d(
#         ADCP_pnum.flatten(),
#         ADCP.bin_depth.values.flatten(),
#         ADCP.Sh_U.values.flatten(),
#         xi=xaxis,
#         yi=yaxis,
#         fn="mean",
#     )[0]
#     out["Heading"] = grid2d(
#         ADCP.profile_number.values,
#         ADCP.Pressure.values,
#         ADCP["Heading"].values,
#         xi=xaxis,
#         yi=yaxis,
#         fn="mean",
#     )[0]
#     out["Pitch"] = grid2d(
#         ADCP.profile_number.values,
#         ADCP.Pressure.values,
#         ADCP["Pitch"].values,
#         xi=xaxis,
#         yi=yaxis,
#         fn="mean",
#     )[0]
#     out["Roll"] = grid2d(
#         ADCP.profile_number.values,
#         ADCP.Pressure.values,
#         ADCP["Roll"].values,
#         xi=xaxis,
#         yi=yaxis,
#         fn="mean",
#     )[0]
#     out["latitude"] = grid2d(
#         ADCP.profile_number.values,
#         ADCP.Pressure.values,
#         ADCP["Latitude"].values,
#         xi=xaxis,
#         yi=yaxis,
#         fn="mean",
#     )[0]
#     out["longitude"] = grid2d(
#         ADCP.profile_number.values,
#         ADCP.Pressure.values,
#         ADCP["Longitude"].values,
#         xi=xaxis,
#         yi=yaxis,
#         fn="mean",
#     )[0]
#     out["profile_number"] = grid2d(
#         ADCP.profile_number.values,
#         ADCP.Pressure.values,
#         ADCP["profile_number"].values,
#         xi=xaxis,
#         yi=yaxis,
#         fn="mean",
#     )[0]
#     out["Pressure"] = grid2d(
#         ADCP.profile_number.values,
#         ADCP.Pressure.values,
#         ADCP["Pressure"].values,
#         xi=xaxis,
#         yi=yaxis,
#         fn="mean",
#     )[0]

#     out = _grid_glider_data(glider, out, xaxis, yaxis)

#     return out


# def verify_depth_bias(out, yaxis, E="ADCP_E", N="ADCP_N"):
#     north = np.gradient(out["latitude"], axis=1) > 0
#     south = np.gradient(out["latitude"], axis=1) < 0
#     depths = np.linspace(0, np.max(yaxis) - 5, 20)
#     drange = np.mean(np.diff(depths)) / 2
#     bins = np.linspace(-1, 1, 100) * 0.5

#     variables = [E, N]

#     for idx, var in enumerate(variables):

#         for d in depths:
#             depth = np.abs(out["Pressure"] - d) < drange

#             Nvals = out[var][(north & depth)]
#             Svals = out[var][(south & depth)]
#             N, _ = np.histogram(Nvals, bins=bins, density=True)
#             S, _ = np.histogram(Svals, bins=bins, density=True)

#             N[N == 0] = np.nan
#             S[S == 0] = np.nan
#     for d in depths:
#         depth = np.abs(out["Pressure"] - d) < drange

#         Nvals = np.sqrt(out["ADCP_E"] ** 2 + out["ADCP_N"] ** 2)[(north & depth)]
#         Svals = np.sqrt(out["ADCP_E"] ** 2 + out["ADCP_N"] ** 2)[(south & depth)]
#         N, _ = np.histogram(Nvals, bins=bins, density=True)
#         S, _ = np.histogram(Svals, bins=bins, density=True)

#         N[N == 0] = np.nan
#         S[S == 0] = np.nan


# def calc_bias(out, yaxis):
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

#     out["ADCP_E"] = out["ADCP_E"] + get_bias(out["speed_e"], coeff)
#     out["ADCP_N"] = out["ADCP_N"] + get_bias(out["speed_n"], coeff)
#     return out


def make_dataset(out):
    profiles = np.arange(out["Pressure"].shape[1])
    pressure_bins = np.arange(out["Pressure"].shape[0])

    ds_dict = {}
    for key, val in out.items():
        if key in ["pressure_bin", "profile_num"]:
            continue
        ds_dict[key] = (
            (
                "pressure_bin",
                "profile_num",
            ),
            val,
        )
    coords_dict = {
        "profile_num": ("profile_num", profiles),
        "pressure_bin": ("pressure_bin", pressure_bins),
    }
    ds = xr.Dataset(data_vars=ds_dict, coords=coords_dict)
    ds["profile_datetime"] = (
        "profile_num",
        pd.to_datetime(ds.date_float.mean(dim="pressure_bin").values),
    )
    return ds


def shear_from_adcp(adcp_file_path, glider_file_path, options=None):
    if not options:
        options = default_options
    ADCP, data, options = load_adcp_glider_data(
        adcp_file_path, glider_file_path, options
    )
    ADCP = remapADCPdepth(ADCP, options)
    ADCP = correct_heading(ADCP, data, options)
    ADCP = soundspeed_correction(ADCP)
    ADCP = remove_outliers(ADCP, options)
    # ADCP = correct_shear(ADCP, options)
    # ADCP = correct_backscatter(ADCP, data)
    ADCP = regridADCPdata(ADCP, options)
    ADCP = calcXYZfrom3beam(ADCP, options)
    ADCP = calcENUfromXYZ(ADCP, options)
    return ADCP, data


# def grid_shear(ADCP, data):
#     xaxis, yaxis, taxis, days = grid_shear_data(data)
#     out = grid_data(ADCP, data, {}, xaxis, yaxis)
#     ds = make_dataset(out)
#     return ds


# def velocity_from_shear(adcp_file_path, glider_file_path, ADCP, data, options=None):
#     if not options:
#         options = default_options
#     extra_data = pd.read_parquet(glider_file_path)
#     extra_data.index = data.index
#     data["speed_vert"] = extra_data["speed_vert"]
#     data["speed_horz"] = extra_data["speed_horz"]
#     data["dead_reckoning"] = extra_data["dead_reckoning"]
#     data["nav_resource"] = extra_data["nav_resource"]
#     data["dive_number"] = extra_data["dive_number"]
#     xaxis, yaxis, taxis, days = grid_shear_data(data)
#     data = get_DAC(ADCP, data)
#     ADCP = bottom_track(ADCP, adcp_file_path, options)
#     out = reference_shear(ADCP, data, xaxis, yaxis)
#     out = grid_data(ADCP, data, out, xaxis, yaxis)
#     out = calc_bias(out, yaxis)
#     ds = make_dataset(out)
#     return ds
