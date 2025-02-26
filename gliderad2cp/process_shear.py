"""
Calculate shear from a Nortek AD2CP mounted on a glider.
Compatible with any glider as long as the ADCP is returning data for the 4 beams.
TODO: make compatible with the 3-beam configuration.


gliderad2cp.process_shear
--------------------------
process
    The main function which converts beam velocity data shear in the ENU directions.
load_data
    Loads AD2CP netcdf file and creates supplemental variables based on input glider data.
_velocity_soundspeed_correction
    Correct for salinity error in the soundspeed calculation.
_quality_control_velocities
    Removes bad velocity measurements based on velocity, amplitude and correlation criteria defined in options.
_determine_velocity_measurement_depths
    Determines what depth each velocity measurement was taken at for each beam, account for glider attitude.
_regrid_beam_velocities_to_isobars
    Regrids beam velocities onto isobars to avoid shear smearing in the final shear data.
_rotate_BEAMS_to_XYZ
    Transforms velocity measurements from along-beam coordinates to XYZ coordinates.
_rotate_XYZ_to_ENU
    Rotates velocity measurements in XYZ coordinates into ENU coordinates.
    

Notes
-------
.process() runs the following functions in this order

1. load_data
2. (correct_heading - Awaiting publication, contact Bastien Queste for access)
3. _velocity_soundspeed_correction
4. _quality_control_velocities
5. _determine_velocity_measurement_depths 
6. _regrid_beam_velocities_to_isobars
7. _rotate_BEAMS_to_XYZ
8. _rotate_XYZ_to_ENU
    
"""

import warnings
from glob import glob
import pandas as pd
import numpy as np
import datetime
import xarray as xr
import gsw
from scipy.interpolate import interp1d
from .tools import plog, interp, get_options

warnings.filterwarnings(action='ignore', message='Mean of empty slice')
warnings.filterwarnings(action='ignore', message='invalid value encountered in divide')
warnings.filterwarnings(action='ignore', message='invalid value encountered in true_divide')
warnings.filterwarnings(action='ignore', message='Degrees of freedom <= 0 for slice.')


"""
Data loading functions
"""
def load_data(adcp_file_path, glider_file_path, options):
    """
    Loads AD2CP netcdf file and creates supplemental variables based on input glider data.


    Inputs
    ----------
    adcp_file_path : str
        Path to the AD2CP netcdf files created by the Nortek MIDAS software. Can handle wildcards through glob.
    glider_file_path : str
        Path to the glider file containing time, latitude, longitude, soundspeed, and profile_num.
        Can be csv, pandas or xarray. Can also pass pandas or xarray dataframes directly.
    options : dict
        Set of options for gliderAD2CP, created by the gliderad2cp.tools.get_options() function.


    Outputs
    -------
    ADCP : xr.Dataset
        xr.Dataset of the AD2CP netcdf with necessary glider variables interpolated onto the time dimension, with some minor cleaning of coordinates.
    options : dict
        Same as the input options with the very important addition of ADCP mounting direction if the 'auto' option is set for 'ADCP_mounting_direction'.
    glider_data : pd.Dataframe
        pd.Dataframe of the input glider data. Serves no real purpose but may be helpful in debugging. Discarded when called from .process()
    
    """
    def __load_glider_data(glider_file):
        """
        Loads glider data.


        Inputs
        ----------
        glider_file_path : str
            Path to the glider file containing time, latitude, longitude, soundspeed, and profile_num. Can be csv, pandas or xarray.


        Outputs
        -------
        glider_data : pd.Dataframe
            pd.Dataframe of the input glider data. Serves no real purpose but may be helpful in debugging. Discarded when called from .process()

        """
        if type(glider_file) is pd.core.frame.DataFrame:
            data = glider_file
        elif type(glider_file) is xr.core.dataset.Dataset:
            data = glider_file.to_dataframe()
            if 'time' not in data:
                data['time'] = data.index        
            data['date_float'] = data['time'].values.astype('float')
            if 'profile_index' in data:
                data['profile_number'] = data['profile_index']
            if 'profile_num' in data:
                data['profile_number'] = data['profile_num']
        elif str(glider_file)[-4:] == '.csv':
            data = pd.read_table(glider_file, parse_dates=['time'])
        elif str(glider_file)[-4:] == '.pqt':
            data = pd.read_parquet(glider_file)
        elif str(glider_file)[-3:] == '.nc':
            data = xr.open_dataset(glider_file).to_dataframe()
            if 'time' not in data:
                data['time'] = data.index        
            data['date_float'] = data['time'].values.astype('float')
            if 'profile_index' in data:
                data['profile_number'] = data['profile_index']
            if 'profile_num' in data:
                data['profile_number'] = data['profile_num']
        else:
            raise(ValueError(f"could not parse input dataset {glider_file}"))
        sel_cols = [
            'time',
            'temperature',
            'salinity',
            'latitude',
            'longitude',
            'profile_number',
            'pressure',
            'dive_num'
        ]
        valid_cols = list(set(list(data)).intersection(sel_cols))
        data = data[valid_cols]
        time_ms = data.time.values
        if time_ms.dtype != '<M8[ns]':
            divisor = 1e3
            if time_ms.dtype == '<M8[us]':
                divisor = 1e6
            time_float = time_ms.astype('float') / divisor
            base = datetime.datetime(1970, 1, 1)
            time_ms = []
            for seconds in time_float:
                time_ms.append(base + datetime.timedelta(seconds=seconds + 0.0000001))
            time_nanoseconds = pd.to_datetime(time_ms)
            data['time'] = time_nanoseconds
            
            
        data['date_float'] = data['time'].values.astype('float')
        p = data['pressure']
        SA = gsw.conversions.SA_from_SP(data.salinity, p, data.longitude, data.latitude)
        CT = gsw.CT_from_t(SA, data['temperature'], p)
        data['soundspeed'] = gsw.sound_speed(SA, CT, p)
        data.index = data.time
        data.index.name = None
        return data
    
    ## Load glider data
    plog('Loading glider data.')
    glider_data = __load_glider_data(glider_file_path)
    
    ## Load ADCP netcdf file(s)
    plog('Loading ADCP data.')
    ADCP = xr.open_mfdataset(adcp_file_path, group='Data/Average')
    ADCP_settings = xr.open_mfdataset(glob(adcp_file_path)[0], group='Config')
    ADCP.attrs = ADCP_settings.attrs
    
    plog('Merging glider data into ADCP dataset.')
    adcp_time_float = ADCP.time.values.astype('float')
    glider_time_float = glider_data['date_float'].values
    ADCP = ADCP.drop_vars(['MatlabTimeStamp']) # In protest of closed source software.
    
    # Coordinates
    ADCP = ADCP.assign_coords(Latitude  = ('time', interp(glider_time_float, glider_data['latitude'], adcp_time_float)))
    ADCP = ADCP.assign_coords(Longitude = ('time',interp(glider_time_float, glider_data['longitude'], adcp_time_float)))

    # Profile and depth
    ADCP = ADCP.assign_coords(profile_number=('time', np.round(interp(glider_time_float, glider_data['profile_number'], adcp_time_float))))
    ADCP = ADCP.assign_coords(Depth=('time', -gsw.z_from_p(ADCP['Pressure'].values, ADCP['Latitude'].values)) )
    
    ADCP['glider_soundspeed'] = ('time', interp(glider_time_float, glider_data['soundspeed'], adcp_time_float) )

    # Get rid of pointless dimensions and make them coordinates instead
    ADCP = ADCP.assign_coords( bin=('Velocity Range', np.arange(len(ADCP['Velocity Range'].values)) ) )
    ADCP = ADCP.swap_dims({'Velocity Range': 'bin'})

    ADCP = ADCP.assign_coords( bin=('Correlation Range', np.arange(len(ADCP['Correlation Range'].values)) ) )
    ADCP = ADCP.swap_dims({'Correlation Range': 'bin'})

    ADCP = ADCP.assign_coords( bin=('Amplitude Range', np.arange(len(ADCP['Amplitude Range'].values)) ) )
    ADCP = ADCP.swap_dims({'Amplitude Range': 'bin'})

    if options['ADCP_mounting_direction'] == 'auto':
        # Detect if ADCP is top mounted using the accelerometer
        if np.nanmedian(ADCP.AccelerometerZ.values) < 0:
            options['ADCP_mounting_direction'] = 'bottom'
        else:
            options['ADCP_mounting_direction'] = 'top'
        plog(f'ADCP is '+options['ADCP_mounting_direction']+'-mounted according to the accelerometer.')

    return ADCP, options, glider_data


"""
QC functions
"""
def _velocity_soundspeed_correction(ADCP):
    """
    Corrects along beam velocity data for lack of salinity measurements in ADCP default soundspeed settings.


    Inputs
    ----------
    ADCP : xr.Dataset
        GliderAD2CP format xr.Dataset produced by the gliderad2cp.load_data() function.


    Outputs
    -------
    ADCP : xr.Dataset
        Same as input, with correction for sound speed applied to along beam velocities.
    
    """
    if 'NoSal_SpeedOfSound' not in ADCP:
        plog('Performing sound speed correction.')
        ADCP = ADCP.rename({'SpeedOfSound': 'NoSal_SpeedOfSound'})
        ADCP = ADCP.rename({'glider_soundspeed': 'SpeedOfSound'})
        for beam in ['1', '2', '3', '4']:
            # V_new = V_old * (c_new/c_old)
            ADCP['VelocityBeam' + beam] = ADCP['VelocityBeam' + beam] * (
                ADCP['SpeedOfSound'] / ADCP['NoSal_SpeedOfSound']
            )
            plog('    Corrected beam ' + beam + ' velocity for sound speed.')
    else:
        plog('Speed of sound correction already applied')
    return ADCP


def _quality_control_velocities(ADCP, options):
    """
    Removes bad velocity measurements based on velocity, amplitude and correlation criteria defined in options.


    Inputs
    ----------
    ADCP : xr.Dataset
        GliderAD2CP format xr.Dataset produced by the gliderad2cp.load_data() function.
    options : dict
        Set of options for gliderAD2CP, created by the gliderad2cp.tools.get_options() function. Requires 'ADCP_mounting_direction' to be either 'top' or 'bottom'.


    Outputs
    -------
    ADCP : xr.Dataset
        Same as input with bad velocity data replaced with NaNs.
    
    """
    prct = lambda x : '%.2f'%(np.count_nonzero(x) * 100 / n)
    plog('Removing data exceeding correlation, amplitude and velocity thresholds.')
    
    noise_floor = np.nanpercentile(np.concatenate([
        ADCP.AmplitudeBeam1.values,
        ADCP.AmplitudeBeam2.values,
        ADCP.AmplitudeBeam3.values,
        ADCP.AmplitudeBeam4.values
        ]).flatten(),[0.5])
    
    for beam in ['1', '2', '3', '4']:
        
        # Calculate correlation mask
        C = ADCP['CorrelationBeam' + beam].values.copy()
        n = len(C.flatten()) # For percentage calculation above. 
        ind = C < options['QC_correlation_threshold']
        C[ind] = np.nan
        C[np.isfinite(C)] = 1
        plog('    Beam ' + beam + ' correlation: ' + prct(ind) + '% removed')

        # Calculate amplitude mask
        A = ADCP['AmplitudeBeam' + beam].values.copy()
        ind = (A > options['QC_amplitude_threshold']) | (A < noise_floor + options['QC_SNR_threshold'])
        A[ind] = np.nan
        A[np.isfinite(A)] = 1
        plog('    Beam ' + beam + ' amplitude: ' + prct(ind) + '% removed')

        # Calculate velocity mask
        V = ADCP['VelocityBeam' + beam].values.copy()
        ind = np.abs(V) > options['QC_velocity_threshold']
        V[ind] = np.nan
        V[np.isfinite(V)] = 1
        plog('    Beam ' + beam + ' velocity: ' + prct(ind) + '% removed')
        
        # Remove bad data
        ADCP['VelocityBeam' + beam] = ADCP['VelocityBeam' + beam] * C * A * V

    return ADCP


"""
Remapping functions
"""
def _determine_velocity_measurement_depths(ADCP, options):    
    """
    Determines what depth each velocity measurement was taken at for each beam, account for glider attitude.


    Inputs
    ----------
    ADCP : xr.Dataset
        GliderAD2CP format xr.Dataset produced by the gliderad2cp.load_data() function.
    options : dict
        Set of options for gliderAD2CP, created by the gliderad2cp.tools.get_options() function. Requires 'ADCP_mounting_direction' to be either 'top' or 'bottom'.


    Outputs
    -------
    ADCP : xr.Dataset
        Same as input with additional D1, D2, D3 and D4 depth variabiles.
    
    """
    H = ADCP['Heading'] + options['heading_offset']
    P = ADCP['Pitch'] + options['pitch_offset']
    R = ADCP['Roll'] + options['roll_offset']
    
    
    # All *_Range coordinates are distance along beam. Verified with data.
    if options['ADCP_mounting_direction'] == 'top':
        direction = 1
        theta_rad_1 = np.arccos(
            np.cos(np.deg2rad(47.5 - P)) * np.cos(np.deg2rad(R))
        )
        theta_rad_2 = np.arccos(
            np.cos(np.deg2rad(25 - R)) * np.cos(np.deg2rad(P))
        )
        theta_rad_3 = np.arccos(
            np.cos(np.deg2rad(47.5 + P)) * np.cos(np.deg2rad(R))
        )
        theta_rad_4 = np.arccos(
            np.cos(np.deg2rad(25 + R)) * np.cos(np.deg2rad(P))
        )
    else:
        direction = -1
        theta_rad_1 = np.arccos(
            np.cos(np.deg2rad(47.5 + P)) * np.cos(np.deg2rad(R))
        )
        theta_rad_2 = np.arccos(
            np.cos(np.deg2rad(25 + R)) * np.cos(np.deg2rad(P))
        )
        theta_rad_3 = np.arccos(
            np.cos(np.deg2rad(47.5 - P)) * np.cos(np.deg2rad(R))
        )
        theta_rad_4 = np.arccos(
            np.cos(np.deg2rad(25 - R)) * np.cos(np.deg2rad(P))
        )
    
    # The above functions return angles of each beam from the UP direction based on the attitude of the glider.
    # Technically, if the glider is pitched 17.4 degrees, 3 beams are at 30.1 degrees.
    # To remap bins, we need to understand how they are mapped. Is the range equal to distance along beam, or distance along
    # an instrument relative Z-axis.

    # Following pers. comm. with Nortek staff, it has been confirmed that the Velocity Range values are equal to distances
    # from the ADCP along a Z-axis. This means that distance along beam needs to be remapped by the cosine of the beam-to-Z-axis cosine.
    # HOWEVER counterintuitively - they did not use 30.1 degrees. See below correspondence:
    #
    # Sven Nylund<Sven.Nylund@nortekgroup.com> : Sent 27 November 2024 10:28 Swedish time to Bastien Queste.
    # "Sorry for the confusion regarding slant angle yesterday, the correction for slant angle is done at a system level so 1 MHz
    # instruments on the AD2CP platform (like your glider unit) are corrected for a nominal 25-degree slant angle whenever more 
    # than two beams are enabled for measurement. So, in your case, the along beam cell size will be 2.21 m when you set the cell
    # size to 2 m. Like I said yesterday, there is no correction of this based on tilt measurements. There is neither any correction 
    # of the cell size with regards to the sound velocity, we do the time gating based on a nominal sound velocity of 1500 m/s." 
    # 
    # So it is hard-coded for 25 degrees, which doesn't really mesh with the geometry but it's fine as long as we know.
    # It will be important to keep an eye on firmware updates in case they adjust this down the line as the information is not 
    # included in any attributes.
    z_bin_distance = ADCP['Velocity Range'].values/np.cos(np.deg2rad(25))

    ADCP['D1'] = (
        ['time', 'bin'],
        np.tile(ADCP['Depth'], (len(ADCP.bin), 1)).T
        - direction
        * np.tile(z_bin_distance, (len(ADCP.time), 1))
        * np.tile(np.cos(theta_rad_1), (len(ADCP.bin), 1)).T,
    )
    ADCP['D2'] = (
        ['time', 'bin'],
        np.tile(ADCP['Depth'], (len(ADCP.bin), 1)).T
        - direction
        * np.tile(z_bin_distance, (len(ADCP.time), 1))
        * np.tile(np.cos(theta_rad_2), (len(ADCP.bin), 1)).T,
    )
    ADCP['D3'] = (
        ['time', 'bin'],
        np.tile(ADCP['Depth'], (len(ADCP.bin), 1)).T
        - direction
        * np.tile(z_bin_distance, (len(ADCP.time), 1))
        * np.tile(np.cos(theta_rad_3), (len(ADCP.bin), 1)).T,
    )
    ADCP['D4'] = (
        ['time', 'bin'],
        np.tile(ADCP['Depth'], (len(ADCP.bin), 1)).T
        - direction
        * np.tile(z_bin_distance, (len(ADCP.time), 1))
        * np.tile(np.cos(theta_rad_4), (len(ADCP.bin), 1)).T,
    )
    
    return ADCP


def _regrid_beam_velocities_to_isobars(ADCP, options):
    """
    Regrids beam velocities onto isobars to avoid shear smearing in the final shear data.


    Inputs
    ----------
    ADCP : xr.Dataset
        GliderAD2CP format xr.Dataset produced by the gliderad2cp.load_data() function.
    options : dict
        Set of options for gliderAD2CP, created by the gliderad2cp.tools.get_options() function. Requires 'ADCP_mounting_direction' to be either 'top' or 'bottom'.


    Outputs
    -------
    ADCP : xr.Dataset
        Same as input but with the additional gridded_bin dimension (defined as range from the glider), 
        and the V1, V2, V3, V4 velocities which are velocities from beams 1-4 projected onto isobars.
    
    """
    if options['velocity_regridding_distance_from_glider'] == 'auto':
        ## This is to avoid shear smearing because of tilted ADCP
        if options['ADCP_mounting_direction'] == 'top':
            direction = 1
        else:
            direction = -1
        bin_size = ADCP.attrs['avg_cellSize']
        blanking_distance = ADCP.attrs['avg_blankingDistance']
        n_bins = len(ADCP.bin.values)
        max_distance = blanking_distance + n_bins * bin_size + (0.5 * bin_size)
        depth_offsets = np.arange(0, max_distance + bin_size, bin_size / 2) * direction
    else:
        depth_offsets = options['velocity_regridding_distance_from_glider']

    plog('Using the following depth offsets when regridding velocities onto isobars:')
    plog(depth_offsets)
    
    ## Extract to np array for speed
    for beam in ['1', '2', '3', '4']:
        plog(f'    Regridding beam {beam} onto isobars.')

        def __interp1d_np(x, y):
            _gd = np.isfinite(y)
            if np.count_nonzero(_gd) > 1:
                xi = interp1d(x[_gd], y[_gd], bounds_error=False, fill_value=np.nan)(
                    depth_offsets
                )
            else:
                xi = depth_offsets * np.nan
            return xi

        ADCP.load()
        ADCP['V' + beam] = xr.apply_ufunc(
            __interp1d_np,
            ADCP['Depth'] - ADCP['D' + beam],
            ADCP['VelocityBeam' + beam],
            input_core_dims=[['bin'], ['bin']],
            output_core_dims=[['gridded_bin']],
            exclude_dims={'bin'},
            vectorize=True,
            output_sizes={'gridded_bin': len(depth_offsets)},
        )

    ADCP = ADCP.assign_coords({'depth_offset': (['gridded_bin'], depth_offsets)})
    ADCP = ADCP.assign_coords(
        {
            'bin_depth': (
                ['time', 'gridded_bin'],
                np.tile(
                    ADCP['Depth'].values.astype('float'), (len(ADCP.gridded_bin), 1)
                ).T
                - np.tile(depth_offsets, (len(ADCP.time), 1)),
            )
        }
    )
    return ADCP


def _rotate_BEAMS_to_XYZ(ADCP, options):
    """
    Coordinate transform that converts BEAM velocities to glider relative velocities X, Y, Z.
    Those of you with a keen eye will notice we rely on the redundancy of 4 beams to calculate 3 beam which allows you to estimate values of any missing beam from the others.
    We do this for a specific - and I promise you, rather clever - reason.
    You see, if we could calculate from 3 beams, but X would not be aligned to glider pitch.
    So either we rotate the resulting matrix, or we simply calculate from the original 4 beam matrix.
    
    TODO : Bastien Queste : double check a 3-beam algorithm gives the same result.


    Inputs
    ----------
    ADCP : xr.Dataset
        GliderAD2CP format xr.Dataset produced by the gliderad2cp.load_data() function.
    options : dict
        Set of options for gliderAD2CP, created by the gliderad2cp.tools.get_options() function. Requires 'ADCP_mounting_direction' to be either 'top' or 'bottom'.


    Outputs
    -------
    ADCP : xr.Dataset
        Same as input but with the additional X, Y and Z velocities (X aligned with the glider, Z aligned with ADCP direction)
    
    """
    plog('Calculating X,Y,Z from isobaric 3-beam measurements.')
    a = lambda t : 1 / (2 * np.sin(t * np.pi / 180))
    b = lambda t : 1 / (4 * np.cos(t * np.pi / 180))
    c = 1  # for convex transducer head
    tf = 47.5  # theta front - Beam 1 and 3 angle from Z
    ts = 25  # theta side - Beam 2 and 4 angle from Z

    V1 = ADCP['V1'].values
    V2 = ADCP['V2'].values
    V3 = ADCP['V3'].values
    V4 = ADCP['V4'].values
    
    upcasts = (ADCP['Pitch']+options['pitch_offset']) > 0
    downcasts = ~upcasts
    
    replaced_by = lambda good_beam : (2 * b(ts) * V2 + 2 * b(ts) * V4 - 2 * b(tf) * good_beam) / (2 * b(tf))

    if options['ADCP_mounting_direction'] == 'top':
        V1[downcasts, :] = replaced_by(V3)[
            downcasts, :
        ]  # Use aft beam on downcasts when upward facing
        V3[upcasts, :] = replaced_by(V1)[
            upcasts, :
        ]  # Use fore beam on upcasts when upward facing
    else:
        V1[upcasts, :] = replaced_by(V3)[
            upcasts, :
        ]  # Use aft beam on upcasts when downward facing
        V3[downcasts, :] = replaced_by(V1)[
            downcasts, :
        ]  # Use fore beam on downcasts when downward facing

    X = c * a(tf) * V1 - c * a(tf) * V3
    Y = -c * a(ts) * V2 + c * a(ts) * V4  # TODO: sign uncertainty here
    Z = 2 * b(ts) * V2 + 2 * b(ts) * V4

    ADCP['X'] = (['time', 'gridded_bin'], X)
    ADCP['Y'] = (['time', 'gridded_bin'], Y)
    ADCP['Z'] = (['time', 'gridded_bin'], Z)

    return ADCP


def _rotate_XYZ_to_ENU(ADCP, options):
    """
    Coordinate transform that converts velocity estimates relative to the glider (X, Y, Z) into the earth relative reference frame east, north, up (ENU)


    Inputs
    ----------
    ADCP : xr.Dataset
        GliderAD2CP format xr.Dataset produced by the gliderad2cp.load_data() function.
    options : dict
        Set of options for gliderAD2CP, created by the gliderad2cp.tools.get_options() function. Requires 'ADCP_mounting_direction' to be either 'top' or 'bottom'.


    Outputs
    -------
    ADCP : xr.Dataset
        Same as input but with the additional X, Y and Z velocities (X aligned with the glider, Z aligned with ADCP direction)
    
    """
    def __M_xyz2enu(heading, pitch, roll):
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

    if options['ADCP_mounting_direction'] == 'top':
        direction = 1
    else:
        direction = -1

    plog('Calculating E, N, U velocities from X, Y, Z velocities.')
    H = ADCP['Heading'] + options['heading_offset']
    P = ADCP['Pitch'] + options['pitch_offset']
    R = ADCP['Roll'] + options['roll_offset']
    M = __M_xyz2enu(H, P, R)

    E = (
        M[0][0] * ADCP['X']
        + M[0][1] * ADCP['Y'] * direction
        + M[0][2] * ADCP['Z'] * direction
    )
    N = (
        M[1][0] * ADCP['X']
        + M[1][1] * ADCP['Y'] * direction
        + M[1][2] * ADCP['Z'] * direction
    )
    U = (
        M[2][0] * ADCP['X']
        + M[2][1] * ADCP['Y'] * direction
        + M[2][2] * ADCP['Z'] * direction
    )

    ADCP['E'] = E
    ADCP['N'] = N
    ADCP['U'] = U
    
    plog('Calculating shear from velocities.')
    ADCP['Sh_E'] = (
        ['time', 'gridded_bin'],
        ADCP['E'].differentiate('gridded_bin').values,
    )
    ADCP['Sh_N'] = (
        ['time', 'gridded_bin'],
        ADCP['N'].differentiate('gridded_bin').values,
    )
    ADCP['Sh_U'] = (
        ['time', 'gridded_bin'],
        ADCP['U'].differentiate('gridded_bin').values,
    )
    return ADCP


"""
Main
"""
def process(adcp_file_path, glider_file_path, options=None):
    """
    Calculate shear from a Nortek AD2CP mounted on a glider.
    Compatible with any glider as long as the ADCP is returning data for the 4 beams.
    

    Inputs
    ----------
    adcp_file_path : str
        Path to the AD2CP netcdf files created by the Nortek MIDAS software. Can handle wildcards through glob.
    glider_file_path : str
        Path to the glider file containing time, latitude, longitude, soundspeed, and profile_num.
        Can be csv, pandas or xarray. Can also pass pandas or xarray dataframes directly.
    options : dict
        Set of options for gliderAD2CP, created by the gliderad2cp.tools.get_options() function.
    
        
    Outputs
    -------
    ADCP : xr.Dataset
        xr.Dataset of the AD2CP netcdf with many supplemental variables, most importantly X, Y, Z and E, N, U velocities projected onto isobars to avoid shear smearing.


    Notes
    -------
    .process() runs the following functions in this order
    1. load_data - load data from a Nortek AD2CP netCDF file and a glider data file or dataset
    2. (correct_heading - Awaiting publication, contact Bastien Queste for access)
    3. _velocity_soundspeed_correction- Corrects along beam velocity data for lack of salinity measurements in ADCP default soundspeed settings.
    4. _quality_control_velocities -  Removes bad velocity measurements based on velocity, amplitude and correlation criteria defined in options
    5. _determine_velocity_measurement_depths - Determines what depth each velocity measurement was taken at for each beam, account for glider attitude
    6. _regrid_beam_velocities_to_isobars -  Regrids beam velocities onto isobars to avoid shear smearing in the final shear data
    7. _rotate_BEAMS_to_XYZ -  Coordinate transform that converts BEAM velocities to glider relative velocities X, Y, Z
    8. _rotate_XYZ_to_ENU - Coordinate transform that converts velocity estimates relative to the glider (X, Y, Z) into the earth relative reference frame east, north, up (ENU)
    """
    # Load default options if not present.
    if not options:
        options = get_options(verbose=False)
        plog('Using default set of options. See gliderad2cp.tools.get_options() for settings.')
    
    # Load glider and ADCP data.
    # Merge necessary glider variables into the ADCP dataset.
    # Discard glider data as not necessary for other functions
    ADCP, options, _ = load_data(
                    adcp_file_path, 
                    glider_file_path, 
                    options
                    )
    
    # Correct heading based on magnetometer data.
    # Option to be implemented after publication of methods paper. Contact Bastien Queste for details.
    # if options['correct_compass_calibration']:   
    #     ADCP = correct_heading(ADCP, options)
    
    # Correct data for soundspeed and quality control based on correlation, amplitude and velocity.
    ADCP = _velocity_soundspeed_correction(ADCP)
    ADCP = _quality_control_velocities(ADCP, options)
    
    # Determine depths of velocity measurements and regrid onto isobars to prevent shear smearing.
    ADCP = _determine_velocity_measurement_depths(ADCP, options)  
    ADCP = _regrid_beam_velocities_to_isobars(ADCP, options)
    
    # Now that we have clean beam velocities, we can convert the data into XYZ and ENU.
    ADCP = _rotate_BEAMS_to_XYZ(ADCP, options)
    ADCP = _rotate_XYZ_to_ENU(ADCP, options)
    
    return ADCP