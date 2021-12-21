import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatch
from matplotlib import rc
from astropy import units as u
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord, Angle, get_moon, get_sun

def calc_rahu_ketu_pos(known_eclipse_time,test_date_utc_time,moon_lambda_known_eclipse):
    
    ketu_speed = 360/18.612958 # deg per year

    if (test_date_utc_time.jd > known_eclipse_time.jd):
        time_diff = ((test_date_utc_time - known_eclipse_time).value*u.d).to(u.year)
        ketu_degrees_moved = (ketu_speed*time_diff.value)%360
        print(ketu_degrees_moved)
        ketu_lambda = moon_lambda_known_eclipse - ketu_degrees_moved
        if (ketu_lambda < 0):
            ketu_lambda = 360 - ketu_lambda
    else:
        time_diff = ((known_eclipse_time - test_date_utc_time).value*u.d).to(u.year)
        ketu_degrees_moved = (ketu_speed*time_diff.value)%360
        print(ketu_degrees_moved)
        ketu_lambda = moon_lambda_known_eclipse + ketu_degrees_moved
        if (ketu_degrees_moved > 360):
            ketu_lambda = 360 - ketu_lambda


    if (ketu_lambda < 180):
        rahu_lambda = 180 + ketu_lambda
    else:
        rahu_lambda = (180 + ketu_lambda) - 360

    return rahu_lambda, ketu_lambda




















