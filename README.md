# gliderad2cp
Processing toolchain for AD2CP data collected from autonomous underwater gliders. Processes from raw data to shear values ready for referencing.

### Usage

3 inputs are required:
1. A netcdf (.nc) file of Nortek AD2CP data as output by the Nortek MIDAS post-processing software
2. glider data in a timeseries parquet file (.pqt) file with the following variables: ["time",
            "temperature",
            "salinity",
            "latitude",
            "pressure",
            "longitude",
            "profile_number",
            "declination",
            ]
3. A dictionary of processing options. For examples, see the notebooks described below

The primary functionality of gliderad2cp is to produce a gridded dataset of velocity shear using data from the glider and AD2CP. This is achieved with two functions:

1. `process_adcp.shear_from_adcp(adcp_path, glider_pqt_path, options)` takes as arguments the three inputs above. It returns processed ADCP and glider data
2. `process_adcp.grid_shear(ADCP, data, options)` takes as arguments the output of `shear_from_adcp` and outputs gridded velocity shear

Additionally, gliderad2cp can integrate these profiles of velocity shear into absolute velocity using `process_adcp.velocity_from_shear(adcp_path, glider_pqt_path, options, data, ADCP)`. This requires estimates of speed through water and the following additional variables in the glider .pqt file: ["speed_vert",
            "speed_horz",
            "dead_reckoning",
            "nav_resource",
            "dive_number"]

These utility functions can be controlled with much greater granularity by calling their constittuent functions individually. This process is detailed in detail in the two example notebooks.

### Notebooks

Notebook `01_calculate_velocity_shear.ipynb` demonstrates the core functionality of gliderad2cp. It takes a glider timeseries and a netCDF file created by a Nortek MIDAS from AD2CP data and produces a gridded dataset of velocity shear

Notebook `02_integrate_velocity_shear.ipynb` shows one method of integrating this velocity shear into earth relative absolute velocities.

Both notebooks use data hosted on zenodo, downloaded with [pooch](https://github.com/fatiando/pooch). The datasets can be downloaded here [https://zenodo.org/record/8431329](https://zenodo.org/record/8431329)

### Contributing

gliderad2cp welcomes contributions. Feel free to submit an Issue or Pull Request if you have recommendations or have experienced issues with the package. Please see the community guidelines, CONTRIBUTING.md and make use of the Issue and Pull Request templates.
