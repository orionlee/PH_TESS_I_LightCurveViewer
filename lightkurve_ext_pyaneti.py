#
# Lightkurve Extension to interface with Pyaneti Modeling Library
# - https://github.com/oscaribv/pyaneti
#

from collections.abc import Iterable
import os

from astropy.table import Table
from astropy.time import Time
from astropy import units as u
import numpy as np
import lightkurve as lk

#
# Prepare and export `LightCurve` to Pyaneti input data
#


def to_sector_str(sector):
    def format(a_sector):
        return f"{a_sector:02d}"

    if isinstance(sector, Iterable):
        return "_".join([format(i) for i in sector])
    else:
        format(sector)


def _create_dir_if_needed(path):
    basedir = os.path.dirname(path)
    if not os.path.isdir(basedir):
        os.makedirs(basedir)


def truncate_lc_to_around_transits(lc, transit_specs):
    def calc_duration_to_use(spec):
        """Calc a duration for the purpose of masking,
        by adding transit duration with an additional `surround_time`,
        with a default of (`max(transit duration * 2, 1 day)`
        """
        duration = spec["duration_hr"] / 24
        surround_time = spec.get("surround_time", None)
        if surround_time is None:
            surround_time = max(duration * 3, 1)  # duration + 2 * duration
        return surround_time

    def calc_period_to_use(spec):
        period = spec.get("period", None)
        if period is not None and period > 0:
            return period
        # else period is not really specified,
        # for the use case that a single dip is observed with noted
        # use an arbitrary large period as a filler
        return 99999999

    lc = lc.remove_nans().normalize()

    period = [calc_period_to_use(t) for t in transit_specs]
    duration = [calc_duration_to_use(t) for t in transit_specs]
    transit_time = [t["epoch"] for t in transit_specs]
    # a mask to include the transits and their surrounding (out of transit observations)
    mask = lc.create_transit_mask(period=period, duration=duration, transit_time=transit_time)

    return lc[mask]


def to_pyaneti_dat(lc, out_path, transit_specs, return_processed_lc=False):
    "Output lc data to a file readable by Pyaneti, with lc pre-processed to be suitable for Pyaneti modeling"
    if transit_specs is not None:
        lc = truncate_lc_to_around_transits(lc, transit_specs)

    # finally write to output
    #  lc["time", "flux", "flux_err"]  # somehow "time" is causing problem
    lc1 = type(lc)(time=lc.time.copy(), flux=lc.flux, flux_err=lc.flux_err)
    _create_dir_if_needed(out_path)
    lc1.write(out_path, format="ascii.commented_header", overwrite=True)
    if return_processed_lc:
        return lc


def print_input_fit_excerpts(lc, lc_pyaneti_dat_filename, transit_specs, epoch_err, period_err):
    """Output parts of Pyaneti `input_fit.py` based on the specification included"""

    print(
        f"""
#Light curve data file
fname_tr = ['{os.path.basename(lc_pyaneti_dat_filename)}']
"""
    )

    tic = lc.meta.get("TARGETID", None)

    print(
        f"""
#Stellar parameters from:
# https://exofop.ipac.caltech.edu/tess/target.php?id={tic}
# TODO - fill them in:
mstar_mean =  -1
mstar_sigma = -1
rstar_mean = -1
rstar_sigma = -1
tstar_mean = -1
tstar_sigma = -1
"""
    )

    # TODO: handle multiple transits
    # TODO: a better min/max estimate than using 10% of the duration
    if epoch_err is None:
        epoch_err = transit_specs[0]["duration_hr"] * 0.05 / 24  # default to a +- 5% of the duration
    min_t0 = transit_specs[0]["epoch"] - epoch_err
    max_t0 = transit_specs[0]["epoch"] + epoch_err
    time_format = lc.time.format
    print(
        f"""
#Minimum and maximum limits for T0
min_t0  = [{min_t0}]  # {time_format}
max_t0  = [{max_t0}]  # {time_format}
"""
    )

    min_P = transit_specs[0]["period"] - period_err
    max_P = transit_specs[0]["period"] + period_err
    print(
        f"""
#Minimum and maximum limits for P
min_P   = [{min_P}]  # days
max_P   = [{max_P}]  # days
"""
    )


#
# Read Pyaneti lightcurve data as `LightCurve` objects
#


def read_pyaneti_lc_dat(filename, time_format="btjd", time_converter_func=None):
    """Read Pyaneti lightcurve files, e.g.,
    `inpy/.../<starname>.dat`, `outpy/.../<starname>-trdata_lightcurve.txt`."""
    # format="ascii.commented_header" does not work
    tab = Table.read(filename, format="ascii")
    (
        n_time,
        n_flux,
        n_flux_err,
    ) = tab.colnames

    if time_converter_func is not None:
        time = time_converter_func(tab[n_time])
    else:
        time = Time(tab[n_time], format=time_format)

    lc_cls = lk.LightCurve
    if time.format == "btjd":
        lc_cls = lk.TessLightCurve
    return lc_cls(time=time, flux=tab[n_flux], flux_err=tab[n_flux_err])
