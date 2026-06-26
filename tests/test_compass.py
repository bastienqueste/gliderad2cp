import datetime
import numpy as np
from gliderad2cp import process_compass

def test_declination():
    declination, intensity = process_compass.declination_intensity(10, 60, 0, datetime.datetime(2025, 1, 1))
    assert np.isclose(declination, 4.47926)
    assert np.isclose(intensity, 51511.78821)
