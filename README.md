# gliderad2cp
Processing toolchain for AD2CP data collected from autonomous underwater gliders. Processes from raw data to shear values ready for referencing.

### Notebooks

Notebook `01_calculate_velocity_shear.ipynb` demonstrates the core functionality of gliderad2cp. It takes a glider timeseries and a netCDF file created by a Nortek MIDAS from AD2CP data and produces a gridded dataset of velocity shear

Notebook `02_integrate_velocity_shear.ipynb` shows one method of integrating this velocity shear into earth relative absolute velocities.

Both noteboooks use data hosted on zenodo, downloaded with [pooch](https://github.com/fatiando/pooch). The datasets can be downloaded here [https://zenodo.org/record/8431329](https://zenodo.org/record/8431329)

### Contributing

gliderad2cp welcomes contributions. Feel free to submit an Issue or Pull Request if you have recommendations or have experienced issues with the package. Please see the community guidelines, CONTRIBUTING.md and make use of the Issue and Pull Request templates.
