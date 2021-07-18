from typing import Callable, List, Tuple, Union
from numbers import Number

import numpy as np
from astropy.time import Time
import batman

import lightkurve as lk


def create_lightcurve_from_batman(
    ma: batman.TransitParams, time: Union[Time, Tuple], noise: Union[float, Callable[[List[Number]], List[Number]]], **kwargs
) -> lk.LightCurve:
    """Create a lightcurve using a Transit Model from `batman <https://lweb.cfa.harvard.edu/~lkreidberg/batman/>`_.

    Examples
    --------

    ma = batman.TransitParams()
    ma.t0 = 3.14  # time of inferior conjunction; first transit is X days after start
    ma.per = 10.123  # orbital period
    ma.rp = 6371 / 696342  # 6371 planet radius (in units of stellar radii)
    ma.a = 19  # semi-major axis (in units of stellar radii)
    ma.inc = 90  # orbital inclination (in degrees)
    ma.ecc = 0 # eccentricity
    ma.w = 90  # longitude of periastron (in degrees)
    ma.u = [0.4, 0.4]  # limb darkening coefficients
    ma.limb_dark = "quadratic"  # limb darkening model
    lc = create_lightcurve_from_batman(ma, (0, 100, "btjd", "tdb", 4800), 50 * 10**-6)
    ax = lc.scatter()

    """
    if not isinstance(time, Time):
        time_raw_start, time_raw_stop_exclusive, format, scale, num_samples = time
        time = Time(
            np.linspace(time_raw_start, time_raw_stop_exclusive, num=num_samples, endpoint=False), format=format, scale=scale
        )

    m = batman.TransitModel(ma, time.value, **kwargs)  # initializes model
    flux_model = m.light_curve(ma)  # calculates light curve

    num_samples = len(time)

    if isinstance(noise, Number):
        # case noise is specified as the standard deviation of Gaussian noises
        flux_noise = np.random.normal(0, noise, num_samples)
        flux_err = np.ones_like(flux_noise) * noise
    else:
        # case noise is created by a user-supplied function
        flux_noise = noise(flux_model)
        flux_err = np.full(len(flux_noise), np.nan)

    flux = flux_model + flux_noise

    lc = lk.LightCurve(time=time, flux=flux, flux_err=flux_err)
    lc["flux_model"] = flux_model
    lc["flux_noise"] = flux_noise

    lc.meta["NORMALIZED"] = True
    lc.meta["LABEL"] = f"Synthetic Batman LC - period={ma.per}, epoch={ma.t0}"
    lc.meta["CREATOR"] = "Batman"
    lc.meta["BATMAN_TRANSIT_PARAMS"] = ma

    return lc
