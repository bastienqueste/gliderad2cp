import pooch

server = "https://zenodo.org/record/10606484/files/"
data_source = pooch.create(
    path=pooch.os_cache("gliderad2cp"),
    base_url=server,
    registry={
        "adcp_profiles_160_to_210.nc": "sha256:ee74a61af86547be0ea40d5bc6e09b88fdf89f87c4e48940a8f9fc694c2717b5",
        "glider_profiles_160_to_210.pqt": "sha256:86b63e82a2cbc9128fc015299de09d8090344ff01614d9d9f96a8712a30e63a7",
        "processed_shear_160_to_210.nc": "sha256:fe4ee1367b4a16916d22f9fc7db7f1349fc66df1e42e68369dc3994b4655eb5f",
        "processed_velocity_160_to_210.nc": "sha256:479f81f5e3c501b9df007ff7ac357653d474f9f847661d05ea1fbf4d6c818528",
    },
)
