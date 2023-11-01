# gliderad2cp

gliderad2cp processes data from the Nortek AD2CP acoustic doppler current profiler (ADCP) mounted on a glider. gliderad2cp takes data from the ADCP unit and the glider and combines them to produce estimates of vertical shear of velocity. It also prodives functionality to integrate these velocity shear profiles into absolute earth relative water vlocities.

## Installation

gliderad2cp can be installed with pip

`python -m pip install gliderad2cp`

## Usage

gliderad2cp requires 2 inputs:
1. A netcdf (.nc) file of Nortek AD2CP data as output by the Nortek MIDAS post-processing software. This software should be procured from Nortek.
2. Glider data in a timeseries csv or parquet file (.pqt) file with the following variables: ["time",
            "temperature",
            "salinity",
            "latitude",
            "pressure",
            "longitude",
            "profile_number",
            "declination",
            ]

**N.B.** if declination data are not present in the datafile, they can be added with the function `process_adcp.get_declination`.

The primary functionality of gliderad2cp is to produce a gridded dataset of velocity shear using data from the glider and AD2CP. This is achieved with two steps:
```python
from gliderad2cp import process_adcp
adcp_file_path = "path/to/adcp/files"
glider_file_path = "path/to/glider/files"
ds_adcp, df_glider = process_adcp.shear_from_adcp(adcp_file_path, glider_file_path)
ds_shear = process_adcp.grid_shear(ds_adcp, df_glider)
```

`shear_from_adcp` Returns an xarray dataset of gridded ADCP data in 3 coordinate systems (beam, XYZ and ENU) and a pandas dataframe of glider data

`grid_shear` Combines the ADCP and glider data into a gridded xarray dataset with velocity shear binned on a per-glider profile basis at 1 m vertical resolution.

Additionally, the function `velocity_from_shear` can be used to estimate absolute velocity profiles:

```python
ds_vel = process_adcp.velocity_from_shear(adcp_file_path, glider_file_path, ds_adcp, df_glider)
```

This requires estimates of speed through water and the following additional variables in the glider file: ["speed_vert",
            "speed_horz",
            "dead_reckoning",
            "nav_resource",
            "dive_number"]

### Additional options
By default, the following options are used during processing:

```python
default_options = {
            'debug_plots': True,
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
            'plots_directory': 'plots',
        }
```

These options can be changed by using the `options` kwarg in the relevant functions.

These utility functions can be controlled with much greater granularity by calling their constittuent functions individually. This process is detailed in the two example notebooks.

### Notebooks

Notebook `01_calculate_velocity_shear.ipynb` demonstrates the core functionality of gliderad2cp. It takes a glider timeseries and a netCDF file created by a Nortek MIDAS from AD2CP data and produces a gridded dataset of velocity shear

Notebook `02_integrate_velocity_shear.ipynb` shows one method of integrating this velocity shear into earth relative absolute velocities.

Both notebooks use data hosted on zenodo, downloaded with [pooch](https://github.com/fatiando/pooch). The datasets can be downloaded here [https://zenodo.org/record/8431329](https://zenodo.org/record/8431329)

-------------------------------

# Contributing

gliderad2cp welcomes contributions. Feel free to submit an Issue or Pull Request if you have recommendations or have experienced issues with the package. Please see the community guidelines, CONTRIBUTING.md and make use of the Issue and Pull Request templates.
