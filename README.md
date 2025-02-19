# gliderad2cp

gliderad2cp processes data from the Nortek AD2CP acoustic doppler current profiler (ADCP) mounted on a glider. gliderad2cp takes data from the ADCP unit and the glider and combines them to produce estimates of vertical shear of velocity. The gliderad2cp toolbox greatly simplifies file handling, integration of any glider data to ADCP data, and the complex trigonometry necessary to obtain high quality shear data. In particular, the integration of the Nortek AD2CP varies across glider manufacturers, either using alternating 3-beam configurations between up and down profiles (on the Seaglider or the Spray) or using 4 beams at all times (on the SeaExplorer). This python package allows users to easily load Nortek AD2CP netCDF files and pull the raw data to provide clean shear estimates with consistent processing and quality control independent of which glider they use. Finally, it provides a final referenced velocity profile and corrects for shear bias when the data permits.

Documentation is hosted at [flow-lab.org/gliderad2cp/](https://www.flow-lab.org/gliderad2cp/)

Code is hosted at [https://github.com/bastienqueste/gliderad2cp](https://github.com/bastienqueste/gliderad2cp)

## Installation

gliderad2cp can be installed with pip

`python -m pip install gliderad2cp`

## Usage


gliderad2cp requires 2 inputs:
1. A netcdf (.nc) file of Nortek AD2CP data as output by the Nortek MIDAS post-processing software. This software should be procured from Nortek.
2. Glider data in a timeseries dataset with the following variables:

- time
- pressure
- temperature
- salinity
- latitude
- longitude
- profile_number

### Process shear

The primary functionality of gliderad2cp is to produce a gridded dataset of velocity shear using data from the glider and AD2CP. This is achieved with the function `process_shear.process`:

```python
from gliderad2cp import process_adcp
adcp_file_path = "path/to/adcp/files"
glider_file_path = "path/to/glider/files"
ds_adcp = process_shear.process(adcp_file_path, glider_file_path)
```

### Process currents

Once shear velocities have been calculated using this function, the function `process_currents.process` can be used to estimate absolute velocity profiles:

```python
from gliderad2cp import process_currents
currents, DAC = process_currents.process(ds_adcp, gps_predive, gps_postdive)
```

This requires the output of `process_shear.process` and GPS locations before and after each dive.

### Correct shear bias

Shear bias can be estimated and corrected with the function `process_bais.process`

```python
from gliderad2cp import process_bias
process_bias.process(currents)
```

### Additional options
By default, the following options are used during processing:

```python
    {correct_compass_calibration : [False, 'compass correction algorithm is awaiting publication and will be added upon acceptance. Contact Bastien Queste if you require.'],
    shear_to_velocity_method : ['integrate'],
    ADCP_mounting_direction : ['auto', 'top', 'bottom'],
    QC_correlation_threshold : [80, 'minimum acceptable along-beam correlation value.'],
    QC_amplitude_threshold : [80, 'maximum acceptable along-beam amplitude.'],
    QC_velocity_threshold : [0.8, 'maximum acceptable along-beam velocity in m.s-1.'],
    QC_SNR_threshold : [3, 'minimum acceptable dB above the noise floor.'],
    velocity_regridding_distance_from_glider : ['auto', 'array of depth-offsets from the glider, in m, at which to interpolate beam velocities onto isobars to avoid shear-smearing. Negative for bottom-mounted ADCPs.'],
    xaxis : [1, 'x-axis resolution in number of profiles of the final gridded products.'],
    yaxis : [None, 'If None: ADCP cell size. If int: y-axis resolution in metres of the final gridded products.'],
    weight_shear_bias_regression : [False, True, 'Give greater weight to dives with greater travel distance which can increase signal to noise.'],
    velocity_dependent_shear_bias_correction : [False, True, 'Determine velocity dependent shear-bias correction coefficients rather than constant coefficients.'],
    shear_bias_regression_depth_slice : [(0, 1000), 'A tuple containing the upper and lower depth limits over which to determine shear bias. Helpful to avoid increased noise due to surface variability. For deep diving gliders (500,1000) is good.'],
    pitch_offset : [0, 'value to be added to pitch to correct for transducer-compass misalignment'],
    roll_offset : [0, 'value to be added to roll to correct for transducer-compass misalignment'],}
```

These options can be changed by using the `options` kwarg in the relevant functions.

-------------------------------

# Contributing

gliderad2cp welcomes contributions. Feel free to submit an Issue or Pull Request if you have recommendations or have experienced issues with the package. Please see the community guidelines, CONTRIBUTING.md and make use of the Issue and Pull Request templates.
