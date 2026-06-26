import datetime
import numpy as np
import ppigrf

def declination_intensity(lon: float, lat: float, height: float, date: datetime.datetime):
    """
    Calculate magnetic field declination (positive eastward from north) and total
    magnetic intensity

    Inputs
    ----------
    lon : float
        longitude
    lat : float
        longitude
    height : float
        height above sea level (km)
    date : datetime.datetime
        date

    Outputs
    -------
    declination : float
        Degrees eastward from north
    total_intensity: float
        Intensity of magnetic field (nT)

    """
    be, bn, bu = ppigrf.igrf(lon, lat, height, date)
    __, declination = ppigrf.get_inclination_declination(be, bn, bu)
    total_intensity = np.sqrt(be ** 2 + bn ** 2 + bu ** 2)
    return declination, total_intensity
