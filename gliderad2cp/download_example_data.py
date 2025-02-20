"""
gliderad2cp.download_example_data
-----------------------------------
Functions for downloading sample datasets
"""

import pooch
server = "https://erddap.observations.voiceoftheocean.org/examples/gliderad2cp/"
data_source = pooch.create(
    path=pooch.os_cache("gliderad2cp"),
    base_url=server,
    registry={
        "adcp_profiles_160_to_210.nc": "sha256:ee74a61af86547be0ea40d5bc6e09b88fdf89f87c4e48940a8f9fc694c2717b5",
        "glider_profiles_160_to_210.pqt": "sha256:86b63e82a2cbc9128fc015299de09d8090344ff01614d9d9f96a8712a30e63a7",
        "sea055_M82.nc": "sha256:55f033a6c14a7e1882e8ed981e3dc39ec4efade4a269eb16675fbbeeced826fc",
        "sea055_M82.ad2cp.00000.nc": "sha256:8affe4520133aba6336c4db90cbbe5264fa10db49c21887f921f680d0d6cabe4",
    },
)

def load_sample_dataset(dataset_name="SEA055_M82.nc"):
    """Download sample datasets for use with gliderad2cp

    Parameters
    ----------
    dataset_name: str, optional
        Default is "sea045_20230530T0832_delayed.nc".

    Raises
    ------
    ValueError: 
        If the requests dataset is not known, raises a value error

    Returns
    -------
    xarray.Dataset: 
        Requested sample dataset
    """
    if dataset_name in data_source.registry.keys():
        file_path = data_source.fetch(dataset_name)
        return file_path
    else:
        msg = f"Requested sample dataset {dataset_name} not known. Specify one of the following available datasets: {list(data_source.registry.keys())}"
        raise KeyError(msg)
