import pooch
import requests


# If faster server (callumrollo.com) is available, use it, otherwise default to zenodo
if (
    requests.get(
        "https://callumrollo.com/files/processed_velocity_160_to_210.nc"
    ).status_code
    == 200
):
    server = "https://callumrollo.com/files/"
else:
    server = "https://zenodo.org/record/8431329/files/"
data_source = pooch.create(
    path=pooch.os_cache("gliderad2cp"),
    base_url=server,
    registry={
        "adcp_profiles_10_to_1170.nc": "sha256:a759047946176c704b0aa39701ef23c78cbaf41b6eac6253fa280388e3d63c90",
        "adcp_profiles_10_to_210.nc": "sha256:7d82db74164ba38adf97f43d5d447a867d155ccc1bbe1727d8fe0ffcedc4f0db",
        "adcp_profiles_160_to_210.nc": "sha256:323ff3cc6402b6c7034a57369ee637c1398af38c2d5f876c0456dbbf9928ab6f",
        "glider_profiles_10_to_1170.pqt": "sha256:9066ef8e8009a953ca572ae40b2070c6cf2307722c09c30ebb7a3a0ad9a6ae36",
        "glider_profiles_10_to_210.pqt": "sha256:fbea59a95470b69d1c6c2fe9f7db6cf6a4343a1cc7a471892fc654f09d195dd5",
        "glider_profiles_160_to_210.pqt": "sha256:86b63e82a2cbc9128fc015299de09d8090344ff01614d9d9f96a8712a30e63a7",
        "processed_shear_10_to_1170.nc": "sha256:d8da73d95688dd4be7995af5b06ab774782ab314fa268d1f90d3fbc40cc543a9",
        "processed_shear_10_to_210.nc": "sha256:2b92de87276622a4373e122b59171e911063b75f67501ae808e150e96208281b",
        "processed_shear_160_to_210.nc": "sha256:8e629bde110e9163425ba06fce0e173f384b5e3c1ab014807e903f8408f88f78",
        "processed_velocity_10_to_1170.nc": "sha256:76b8f1e14a1896a7c086b769f848d84f7e61fe6e5d455d59b95b9795c042e0ae",
        "processed_velocity_10_to_210.nc": "sha256:e70be7a7d035b86d6bb25959415a751abea4127c9599cb07c555e9e9da77c356",
        "processed_velocity_160_to_210.nc": "sha256:cb6f0ccd580db111ad6b54da9f8db831632f461740eb78b26141964b6abe97b6",
    },
)
