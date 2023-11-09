import pooch

server = "https://zenodo.org/record/8431329/files/"
data_source = pooch.create(
    path=pooch.os_cache("gliderad2cp"),
    base_url=server,
    registry={
        "adcp_profiles_160_to_210.nc": "sha256:323ff3cc6402b6c7034a57369ee637c1398af38c2d5f876c0456dbbf9928ab6f",
        "glider_profiles_160_to_210.pqt": "sha256:86b63e82a2cbc9128fc015299de09d8090344ff01614d9d9f96a8712a30e63a7",
        "processed_shear_160_to_210.nc": "sha256:8e629bde110e9163425ba06fce0e173f384b5e3c1ab014807e903f8408f88f78",
        "processed_velocity_160_to_210.nc": "sha256:cb6f0ccd580db111ad6b54da9f8db831632f461740eb78b26141964b6abe97b6",
    },
)
