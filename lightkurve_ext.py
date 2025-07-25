"""
Convenience helpers for `lightkurve` package.
"""

# so that TransitTimeSpec can be referenced in type annotation in the class itself
# see: https://stackoverflow.com/a/49872353
from __future__ import annotations


import itertools
import os
import logging
import math
import json
import re
import warnings
from collections import OrderedDict
from types import SimpleNamespace

from retry import retry

import astropy
from astropy.io import fits
from astropy import coordinates as coord
from astropy.coordinates import SkyCoord
from astropy.table import QTable, Table
from astropy.time import Time
import astropy.units as u

import numpy as np
from scipy.interpolate import UnivariateSpline

from IPython.display import display, HTML

import lightkurve as lk

import asyncio_compat

log = logging.getLogger(__name__)


def of_lcs(lc_coll, filter_func):
    """Filter a LightCurveCollection using the given filter_func.

    Example
    --------

    Only retain TESS SPOC 2 minute cadence lightcurves
    > of_lcs(lc_coll, lambda lc: lc.author == "SPOC")
    """
    return type(lc_coll)([lc for lc in lc_coll if filter_func(lc)])


def of_sector(lcf_coll, sectorNum):
    res_list = of_sectors(lcf_coll, sectorNum)
    return res_list[0] if len(res_list) > 0 else None


def of_sectors(*args):
    lk_coll_or_sr = args[0]
    if lk_coll_or_sr is None:
        return []
    if len(args) == 1:
        # when no sectors are specified, return entire collection
        # For convenience: when a notebooks is modified such that
        # a user sometimes use a subset of sectors , and sometimes everything
        # the user can can still use of_sectors() wrapper regardless
        return lk_coll_or_sr
    sector_nums = args[1:]

    if hasattr(lk_coll_or_sr, "sector"):
        return lk_coll_or_sr[np.in1d(lk_coll_or_sr.sector, sector_nums)]
    elif hasattr(lk_coll_or_sr, "table") and lk_coll_or_sr.table["sequence_number"] is not None:
        res = lk.SearchResult(lk_coll_or_sr.table[np.in1d(lk_coll_or_sr.table["sequence_number"], sector_nums)])
        res = _sort_chronologically(res)
        return res
    else:
        raise TypeError(f"Unsupported type of collection: {type(lk_coll_or_sr)}")


def of_sector_n_around(lk_coll_or_sr, sector_num, num_additions=8):
    def get_sector_for_lk_coll(lk_coll):
        return lk_coll.sector

    def get_sector_for_sr(sr):
        return sr.table["sequence_number"]

    sector_accessor_func = None
    if hasattr(lk_coll_or_sr, "sector"):
        sector_accessor_func = get_sector_for_lk_coll
    elif hasattr(lk_coll_or_sr, "table") and lk_coll_or_sr.table["sequence_number"] is not None:
        sector_accessor_func = get_sector_for_sr
    else:
        raise TypeError(f"Unsupported type of collection: {type(lk_coll_or_sr)}")

    subset_slice = _get_slice_for_of_sector_n_around(
        lk_coll_or_sr,
        sector_accessor_func,
        sector_num,
        num_additions=num_additions,
    )
    if subset_slice is not None:
        return lk_coll_or_sr[subset_slice]
    else:
        return type(lk_coll_or_sr)([])


def _get_slice_for_of_sector_n_around(coll, sector_accessor_func, sector_num, num_additions):
    if sector_num not in sector_accessor_func(coll):
        return None

    idx = np.where(sector_accessor_func(coll) == sector_num)[0][0]
    # if num_additions is odd number, we add one to older sector
    start = max(idx - math.ceil(num_additions / 2), 0)
    end = min(idx + math.floor(num_additions / 2) + 1, len(coll))

    # case the start:end slice does not fill up the requested num_additions,
    # try to fill it up
    cur_slice_size = end - start - 1
    if cur_slice_size < num_additions:
        num_more_needed = num_additions - cur_slice_size
        if start > 0:
            start = max(start - num_more_needed, 0)
        else:
            end = min(end + num_more_needed, len(coll))

    return slice(start, end)


def of_2min_cadences(lcf_coll):
    """Return LightCurveFiles of short, typically 2-minute cadence, only.
    Primary use case is to filter out 20-second files.
    """
    filtered = [lcf for lcf in lcf_coll if "short" == estimate_cadence_type(lcf)]
    return lk.LightCurveCollection(filtered)


def estimate_cadence(lc, unit=None, round_unit_result=True):
    """Estimate the cadence of a lightcurve by returning the median of a sample"""
    if isinstance(lc, lk.FoldedLightCurve):
        time_vals = lc.time_original.value
        time_vals = np.sort(time_vals)
    else:
        time_vals = lc.time.value
    res = np.nanmedian(np.diff(time_vals[:100]))
    if unit is not None:
        res = (res * u.day).to(unit)  # LATER: handle cases lc.time is not in days
        if round_unit_result:
            res = res.round()
    return res


def map_cadence_type(cadence_in_days):
    long_minimum = 9.9 / 60 / 24  # 10 minutes in days, with some margin of error
    short_minimum = 0.9 / 60 / 24  # 1 minute in days, with some margin of error
    if cadence_in_days is None:
        return None
    if cadence_in_days >= long_minimum:
        return "long"
    if cadence_in_days >= short_minimum:
        return "short"
    return "fast"


def estimate_cadence_type(lc):
    """Estimate the type of cadence to be one of long, short, or fast.
    The definition is the same as ``exptime`` in `lightkurve.search_lightcurve()`.
    """
    return map_cadence_type(estimate_cadence(lc))


def of_tic(lcf_coll, tic):
    """Return LightCurveFiles of the given TIC.

    Useful in case the default MAST result returned nearby targets.
    """
    filtered = [lcf for lcf in lcf_coll if lcf.meta.get("TICID", None) == tic]
    return lk.LightCurveCollection(filtered)


def select(lcf_coll_or_sr, filter_func):
    """Filter the given LightCurveCollection or SearchResult with the filter_func."""
    return type(lcf_coll_or_sr)([obj for obj in lcf_coll_or_sr if filter_func(obj)])


def exclude_range(lc, start, end, column="time"):
    """Exclude the specified range of time from the given lightcurve."""
    tmask = (lc[column].value >= start) & (lc[column].value < end)
    return lc[~tmask]


def get_obs_date_range(lcf_coll):
    """Return the observation date span and the number of days with observation."""
    # the code assumes the time are in all in BTJD, or other consistent format in days
    if isinstance(lcf_coll, lk.LightCurve):
        lcf_coll = lk.LightCurveCollection([lcf_coll])

    # to support folded lightcurve
    time_colname = "time_original" if "time_original" in lcf_coll[0].colnames else "time"

    t_start = lcf_coll[0][time_colname].min().value
    t_end = lcf_coll[-1][time_colname].max().value

    obs_span = t_end - t_start
    obs_actual = len(set(np.concatenate([lc[time_colname].value.astype("int") for lc in lcf_coll])))

    return obs_span, obs_actual


def estimate_object_radius_in_r_jupiter(lc, depth):
    """Return a back of envelope estimate of a companion object's radius."""
    R_JUPITER_IN_R_SUN = 71492 / 695700

    r_star = lc.meta.get("RADIUS")  # assumed to be in R_sun
    if r_star is None or r_star < 0 or depth <= 0:
        return None  # cannot estimate
    r_obj = math.sqrt(r_star * r_star * depth)
    r_obj_in_r_jupiter = r_obj / R_JUPITER_IN_R_SUN
    return r_obj_in_r_jupiter


def _sort_chronologically(sr: lk.SearchResult):
    # Resort the SearchResult rows, because
    # lightkurve v2.4.2 does not honor mission (chronological)
    # the workaround here is to resort with pre-v2.4.2 criteria
    # Note: we must use a copy of the table
    # because if the underlying table is actually a subset of the underlying table,
    # the sorting seems to affect the underlying table in some cases,
    # creating really weird / corrupted results.
    if len(sr) < 1:
        return sr
    res = lk.SearchResult(sr.table.copy())
    res.table.sort(["distance", "year", "mission", "sort_order", "exptime"])
    return res


# BEGIN lightkurve search with retry

LK_SEARCH_NUM_RETRIES = 4
LK_DOWNLOAD_NUM_RETRIES = 2

@retry(IOError, tries=LK_SEARCH_NUM_RETRIES, delay=0.5, backoff=2, jitter=(0, 0.5))
def search_lightcurve(*args, **kwargs):
    return lk.search_lightcurve(*args, **kwargs)


@retry(IOError, tries=LK_SEARCH_NUM_RETRIES, delay=0.5, backoff=2, jitter=(0, 0.5))
def search_targetpixelfile(*args, **kwargs):
    return lk.search_targetpixelfile(*args, **kwargs)


@retry(IOError, tries=LK_SEARCH_NUM_RETRIES, delay=0.5, backoff=2, jitter=(0, 0.5))
def search_tesscut(*args, **kwargs):
    return lk.search_tesscut(*args, **kwargs)


@retry(IOError, tries=LK_DOWNLOAD_NUM_RETRIES, delay=0.5, backoff=2, jitter=(0, 0.5))
def download_all(search_result, *args, **kwargs):
    return search_result.download_all(*args, **kwargs)


# END   lightkurve search with retry


def download_lightcurves_of_tic_with_priority(
    tic, author_priority=["SPOC", "QLP", "TESS-SPOC"], download_filter_func=None, download_dir=None
):
    """For a given TIC, download lightcurves across all sectors.
    For each sector, download one based on pre-set priority.
    """

    sr_unfiltered = search_lightcurve(f"TIC{tic}", mission="TESS")
    if len(sr_unfiltered) < 1:
        print(f"WARNING: no result found for TIC {tic}")
        return None, None, None

    sr_unfiltered = sr_unfiltered[sr_unfiltered.target_name == str(tic)]  # in case we get some other nearby TICs
    sr_unfiltered = _sort_chronologically(sr_unfiltered)

    # filter out HLSPs not supported by lightkurve yet
    sr = sr_unfiltered[sr_unfiltered.author != "DIAMANTE"]
    if len(sr) < len(sr_unfiltered):
        print("Note: there are products not supported by Lightkurve, which are excluded from download.")

    # for each sector, filter based on the given priority.
    # - note: prefer QLP over TESS-SPOC because QLP is detrended, with multiple apertures within 1 file
    sr = filter_by_priority(
        sr,
        author_priority=author_priority,
        exptime_priority=["short", "long", "fast"],
    )
    num_filtered = len(sr_unfiltered) - len(sr)
    num_fast = len(sr_unfiltered[sr_unfiltered.exptime < 60 * u.second])
    if num_filtered > 0:
        msg = f"{num_filtered} rows filtered"
        if num_fast > 0:
            msg = msg + f" ; {num_fast} fast (20secs) products."
        print(msg)

    with astropy.conf.set_temp("max_lines", -1):
        display(sr)

    # let caller to optionally further restrict a subset to be downloaded
    sr_to_download = sr
    if download_filter_func is not None:
        sr_to_download = download_filter_func(sr)
        sr_to_download = _sort_chronologically(sr_to_download)
        if len(sr_to_download) < len(sr):
            display(
                HTML(
                    """<font style="background-color: yellow;">Note</font>:
SearchResult is further filtered - only a subset will be downloaded."""
                )
            )

    lcf_coll = download_all(sr_to_download, download_dir=download_dir)

    print(f"TIC {tic} \t, all available sectors: {abbrev_sector_list(sr)}")
    if lcf_coll is not None and len(lcf_coll) > 0:
        print(f"downloaded #sectors: {len(lcf_coll)} ; {abbrev_sector_list(lcf_coll)}")
        print(
            (
                f"   sector {lcf_coll[-1].meta['SECTOR']}: \t"
                f"camera = {lcf_coll[-1].meta['CAMERA']} ; ccd = {lcf_coll[-1].meta['CCD']}"
            )
        )
    else:
        print(f"TIC {tic}: no data")

    return lcf_coll, sr, sr_unfiltered


def download_lightcurve(
    target,
    mission=("Kepler", "K2", "TESS"),
    exptime="short",
    author="SPOC",
    download_dir=None,
    use_cache="yes",
    display_search_result=True,
):
    """
    Wraps `lightkurve.search_lightcurve()` and the
    subsequent `lightkurve.search.SearchResult.download_all()` calls,
    with the option of caching, so that for a given search,
    if the the result has been downloaded, the cache will be used.

    The parameters all propagate to the underlying `search_lightcurvefile()`
    and `download_all()` calls. The lone exception is `use_cache`.

    Parameters
    ----------
    use_cache : str, must be one of 'yes', or 'no'\n
        OPEN: an option of 'fallback': cache will be used when offline.\n
        OPEN: for now, actual download lightcurve cache will still be used if
        available irrespective of the settings.

    Returns
    -------
    collection : `~lightkurve.collections.Collection` object
        Returns a `~lightkurve.collections.LightCurveCollection`
        containing all lightcurve files that match the criteria
    """

    if use_cache == "no":
        return _search_and_cache(target, mission, exptime, author, download_dir, display_search_result)
    if use_cache == "yes":
        result_file_ids = _load_from_cache_if_any(target, mission, download_dir)
        if result_file_ids is not None:
            result_files = list(map(lambda e: f"{download_dir}/mastDownload/{e}", result_file_ids))
            return lk.collections.LightCurveCollection(list(map(lambda f: lk.read(f), result_files)))
        # else
        return _search_and_cache(target, mission, exptime, author, download_dir, display_search_result)
    # else
    raise ValueError("invalid value for argument use_cache")


# Private helpers for `download_lightcurvefiles`


def _search_and_cache(target, mission, exptime, author, download_dir, display_search_result):
    search_res = lk.search_lightcurve(target=target, mission=mission, exptime=exptime, author=author)
    if len(search_res) < 1:
        return None
    if display_search_result:
        _display_search_result(search_res)
    _cache_search_result_product_identifiers(search_res, download_dir, target, mission)
    return search_res.download_all(quality_bitmask="default", download_dir=download_dir)


def _display_search_result(search_res):
    from IPython.core.display import display

    tab = search_res.table
    # move useful columns to the front
    preferred_cols = ["proposal_id", "target_name", "sequence_number", "t_exptime"]
    colnames_reordered = preferred_cols + [c for c in tab.colnames if c not in preferred_cols]
    display(tab[colnames_reordered])


def _load_from_cache_if_any(target, mission, download_dir):
    key = _get_cache_key(target, mission)
    return _load_search_result_product_identifiers(download_dir, key)


def _cache_search_result_product_identifiers(search_res, download_dir, target, mission):
    key = _get_cache_key(target, mission)
    identifiers = _to_product_identifiers(search_res)
    _save_search_result_product_identifiers(identifiers, download_dir, key)
    return key


def _get_search_result_cache_dir(download_dir):
    # TODO: handle download_dir is None (defaults)
    cache_dir = f"{download_dir}/mastQueries"

    if os.path.isdir(cache_dir):
        return cache_dir

    # else it doesn't exist, make a new cache directory
    try:
        os.mkdir(cache_dir)
    # downloads locally if OS error occurs
    except OSError:
        log.warning(
            "Warning: unable to create {}. "
            "Cache MAST query results to the current "
            "working directory instead.".format(cache_dir)
        )
        cache_dir = "."
    return cache_dir


def _get_cache_key(target, mission):
    # TODO: handle cases the generated key is not a valid filename
    return f"{target}_{mission}_ids"


def _to_product_identifiers(search_res):
    """
    Returns
    -------
    A list of str, constructed from `(obs_collection, obs_id, productFilename)` tuples, that can
    identify cached lightcurve file,s if any.
    """
    return list(
        map(
            lambda e: e["obs_collection"] + "/" + e["obs_id"] + "/" + e["productFilename"],
            search_res.table,
        )
    )


def _save_search_result_product_identifiers(identifiers, download_dir, key):
    resolved_cache_dir = _get_search_result_cache_dir(download_dir)
    filepath = f"{resolved_cache_dir}/{key}.json"
    fp = open(filepath, "w+")
    json.dump(identifiers, fp)
    return filepath


def _load_search_result_product_identifiers(download_dir, key):
    resolved_cache_dir = _get_search_result_cache_dir(download_dir)
    filepath = f"{resolved_cache_dir}/{key}.json"
    try:
        fp = open(filepath, "r")
        return json.load(fp)
    except OSError as err:
        # errno == 2: file not found, typical case of cache miss
        # errno != 2: unexpected error, log a warning
        if err.errno != 2:
            log.warning("Unexpected OSError in retrieving cached search result: {}".format(err))
        return None


def filter_by_priority(
    sr,
    author_priority=["SPOC", "TESS-SPOC", "QLP"],
    exptime_priority=["short", "long", "fast"],
):
    if sr is None or len(sr) < 1:
        return sr

    author_sort_keys = {}
    for idx, author in enumerate(author_priority):
        author_sort_keys[author] = idx + 1

    exptime_sort_keys = {}
    for idx, exptime in enumerate(exptime_priority):
        exptime_sort_keys[exptime] = idx + 1

    def calc_filter_priority(row):
        # Overall priority key is in the form of <author_key><exptime_key>, e.g., 101
        # - "01" is the exptime_key
        # - the leading "1" is the author_key, given it is the primary one
        author_default = max(dict(author_sort_keys).values()) + 1
        author_key = author_sort_keys.get(row["author"], author_default) * 100

        # secondary priority
        exptime_default = max(dict(exptime_sort_keys).values()) + 1
        exptime_key = exptime_sort_keys.get(map_cadence_type(row["exptime"] / 60 / 60 / 24), exptime_default)
        return author_key + exptime_key

    sr.table["_filter_priority"] = [calc_filter_priority(r) for r in sr.table]

    # A temporary table that sorts the table by the priority
    sorted_t = sr.table.copy()
    sorted_t.sort(["mission", "_filter_priority"])

    # create an empty table for results, with the same set of columns
    res_t = sr.table[np.zeros(len(sr), dtype=bool)].copy()

    # for each mission (e.g., TESS Sector 01), select a row based on specified priority
    # - select the first row given the table has been sorted by priority
    uniq_missions = list(OrderedDict.fromkeys(sorted_t["mission"]))
    for m in uniq_missions:
        mission_t = sorted_t[sorted_t["mission"] == m]
        # OPEN: if for a given mission, the only row available is not listed in the priorities,
        # the logic still add a row to the result.
        # We might want it to be an option specified by the user.
        res_t.add_row(mission_t[0])

    sr = lk.SearchResult(table=res_t)
    sr = _sort_chronologically(sr)
    return sr


# Download TPF asynchronously


def search_and_download_tpf(*args, **kwargs):
    """Search and Download a TPFs.

    All parameters are passed on ``search_targetpixelfile()``,
    with the exception of ``download_dir`` and ``quality_bitmask``,
    which are passed to ``download_all()`
    """

    # extract download_all() parameters
    download_dir = kwargs.pop("download_dir", None)
    quality_bitmask = kwargs.pop("quality_bitmask", None)
    sr = search_targetpixelfile(*args, **kwargs)  # pass the rest of the argument to search_targetpixelfile
    sr = _sort_chronologically(sr)
    tpf_coll = download_all(sr, download_dir=download_dir, quality_bitmask=quality_bitmask)
    return tpf_coll, sr


def create_download_tpf_task(*args, **kwargs):
    return asyncio_compat.create_background_task(search_and_download_tpf, *args, **kwargs)


#
# Other misc. extensions
#
def get_bkg_lightcurve(lcf):
    """Returns the background flux, i.e., ``SAP_BKG`` in the file"""
    lc = lcf.copy()
    lc["flux"] = lc["sap_bkg"]
    lc["flux_err"] = lc["sap_bkg_err"]
    lc.label = lc.label + " BKG"
    return lc


def _do_create_quality_issues_mask(quality, flux, flags_included=0b0101001010111111):
    """Returns a boolean array which flags cadences with *issues*.

    The default `flags_included` is a TESS default, based on
    https://outerspace.stsci.edu/display/TESS/2.0+-+Data+Product+Overview#id-2.0DataProductOverview-Table:CadenceQualityFlags
    """
    if np.issubdtype(quality.dtype, np.integer):
        return np.logical_and(quality & flags_included, np.isfinite(flux))
    else:
        # quality column is not an integer, probably a non-standard product
        return np.zeros_like(quality, dtype=bool)


def create_quality_issues_mask(lc, flags_included=0b0101001010111111):
    """Returns a boolean array which flags cadences with *issues*.

    The default `flags_included` is a TESS default, based on
    https://outerspace.stsci.edu/display/TESS/2.0+-+Data+Product+Overview#id-2.0DataProductOverview-Table:CadenceQualityFlags
    """

    # use sap_flux when available (it may not be there in some HLSP)
    # we prefer sap_flux over pdcsap_flux as
    # pdcsap_flux is more likely to be NaN (due to exclusion by quality flags)
    if "sap_flux" in lc.colnames:
        flux = lc["sap_flux"]
    else:
        flux = lc.flux

    return _do_create_quality_issues_mask(lc.quality, flux)


def _get_n_truncate_fits_data(lc, before, after, return_columns, return_mask=False):
    with fits.open(lc.filename) as hdu:
        time = hdu[1].data["TIME"]
        mask = (time >= before) & (time < after)
        res = dict()
        for col in return_columns:
            res[col] = hdu[1].data[col][mask]
        if return_mask:
            return res, mask
        else:
            return res


def list_times_w_quality_issues(lc, include_excluded_cadences=False):
    if not include_excluded_cadences:
        mask = create_quality_issues_mask(lc)
        return lc.time[mask], lc.quality[mask]
    else:
        # case we want cadences that have been excluded in the lc object
        # use the underlying fits file

        flux_colname = lc.meta.get("FLUX_ORIGIN", "sap_flux")
        if flux_colname == "pdcsap_flux":
            flux_colname = "sap_flux"  # pdcsap_flux would not have values in excluded cadences, defeating the purpose

        # data is truncated if the given lc is truncated
        # TODO: cadences excluded before lc.time.min() or after lc.time.max() will still be missing.
        data = _get_n_truncate_fits_data(lc, lc.time.min().value, lc.time.max().value, ["time", "quality", flux_colname])
        mask = _do_create_quality_issues_mask(data["quality"], data[flux_colname])
        return data["time"][mask], data["quality"][mask]


def list_transit_times(t0, period, steps_or_num_transits=range(0, 10), return_string=False):
    """List the transit times based on the supplied transit parameters"""
    if isinstance(steps_or_num_transits, int):
        steps = range(0, steps_or_num_transits)
    else:
        steps = steps_or_num_transits
    times = [t0 + period * i for i in steps]
    if return_string:
        return ",".join(map(str, times))
    else:
        return times


def get_segment_times_idx(times, break_tolerance=5):
    """Segment the input array of times into segments due to data gaps. Return the indices of the segments.

    The minimal gap size is determined by `break_tolerance`.

    The logic is adapted from `LightCurve.flatten`
    """
    if hasattr(times, "value"):  # convert astropy Time to raw values if needed
        times = times.value
    if len(times) < 1:
        return (None, None)
    dt = times[1:] - times[0:-1]
    median_cadence = np.nanmedian(dt)
    # print(f"TRACE median_cadence={median_cadence * 24 * 60} min,  break_tolerance={break_tolerance}")
    with warnings.catch_warnings():  # Ignore warnings due to NaNs
        warnings.simplefilter("ignore", RuntimeWarning)
        cut = np.where(dt > break_tolerance * median_cadence)[0] + 1
    low = np.append([0], cut)
    high = np.append(cut, len(times))
    return (low, high)


def get_segment_times(times, **kwargs):
    if hasattr(times, "value"):  # convert astropy Time to raw values if needed
        times = times.value
    low, high = get_segment_times_idx(times, **kwargs)
    if low is None:  # edge case times is empty
        return []
    else:  # normal case
        # add a small 1e-10 to end so that the end time is exclusive (follow convention in range)
        return [(times[lo], times[hi - 1] + 1e-10) for lo, hi in zip(low, high)]


def get_transit_times_in_range(t0, period, start, end):
    t_start = t0 + math.ceil((start - t0) / period) * period
    num_t = math.ceil((end - t_start) / period)
    return [t_start + period * i for i in range(num_t)]


def get_transit_times_in_lc(lc, t0, period, return_string=False, **kwargs):
    """Get the transit times with observations of the given lightcurve, based on the supplied transit parameters.

    The method will exclude the times where there is no observation due to data gaps.
    """

    lc = lc.remove_nans()  # exclude cadences with no flux.
    # break up the times to exclude times in gap
    times_list = get_segment_times(lc.time, **kwargs)
    # print(f"TRACE  times_list={times_list}")
    transit_times = []
    for start, end in times_list:
        transit_times.extend(get_transit_times_in_range(t0, period, start, end))
    if return_string:
        return ",".join(map(str, transit_times))
    else:
        return transit_times


def _concatenate_list_of_lists(lists):
    # https://stackoverflow.com/a/14807729

    combined = list(itertools.chain.from_iterable(lists))
    # equivalent to using list comprehension, but itertools is supposed to be faster
    # combined = [item for sublist in lists for item in sublist]
    return combined


def get_compatible_periods_of_2dips(
    epoch1, epoch2, lcc: lk.LightCurveCollection, flux_column="flux", min_period=10, verbose=False
):
    """For the case where 2 dips are observed. Return the compatible periods that fit the observations.

    Note: Parameter `lcc` accepts a `LightCurveCollection` instead of a `LightCurve`. This is because
    the logic of segmenting lightcurve used internally is sensitive to the lightcurve cadence,
    and could yield unexpected results if the cadence is mixed.
    E.g., if users have 1 30-minute cadence LC and 2 2-min cadence LCs and stitches them together,
    the median cadence determined would be 2 min.
    The segmenting logic would then (by default) segmenting the lightcurve when there is
    a 10 min gap (2min X 5). As a result, the 30-min cadence portion of the lightcurve will be be broken
    up such that each cadence is its own segment.
    Counting of transit times could miss some values due to rounding issues in those 1-cadence segments.
    """
    num_cycles = 2
    compat_p_list = []
    while True:
        trial_p = abs(epoch2 - epoch1) / num_cycles
        if trial_p <= min_period:
            break
        trial_ttimes = _concatenate_list_of_lists(
            [get_transit_times_in_lc(lc.select_flux(flux_column), epoch1, trial_p) for lc in lcc]
        )
        if len(trial_ttimes) <= 2:  # assuming epoch1, epoch2 is always in the trial_ttimes.
            compat_p_list.append(trial_p)
            if verbose:
                print(f"  Period {trial_p} compatible. {trial_ttimes}")
        else:
            if verbose:
                print(f"  Period {trial_p} is incompatible. Has unexpected transits at times {trial_ttimes}")
        num_cycles += 1
    return compat_p_list


def get_compatible_periods_of_3dips(
    epoch1, epoch2, epoch3, lcc: lk.LightCurveCollection, flux_column="flux", min_period=10, verbose=False
):
    """For the case where 2 dips are observed. Return the compatible periods that fit the observations.

    Note: Parameter `lcc` accepts a `LightCurveCollection` instead of a `LightCurve`. This is because
    the logic of segmenting lightcurve used internally is sensitive to the lightcurve cadence,
    and could yield unexpected results if the cadence is mixed.
    E.g., if users have 1 30-minute cadence LC and 2 2-min cadence LCs and stitches them together,
    the median cadence determined would be 2 min.
    The segmenting logic would then (by default) segmenting the lightcurve when there is
    a 10 min gap (2min X 5). As a result, the 30-min cadence portion of the lightcurve will be be broken
    up such that each cadence is its own segment.
    Counting of transit times could miss some values due to rounding issues in those 1-cadence segments.
    """
    num_cycles = 2
    compat_p_list = []
    while True:
        trial_p = abs(epoch3 - epoch1) / num_cycles
        if trial_p <= min_period:
            break
        # check if trial_p fits epoch2
        num_cycles_for_epoch2 = round((epoch2 - epoch1) / trial_p)
        epoch2_from_trial_p = epoch1 + trial_p * num_cycles_for_epoch2
        epoch2_diff = abs(epoch2 - epoch2_from_trial_p)
        tolerance = 0.05  # observed epoch2 needs to be within tolerance (in days) with predicted epoch 2 time
        if epoch2_diff > tolerance:
            if verbose:
                print(f"  Period {trial_p} is incompatible. Does not fit epoch2 (off by {epoch2_diff:.2f} days).")
            num_cycles += 1
            continue
        trial_ttimes = _concatenate_list_of_lists(
            [get_transit_times_in_lc(lc.select_flux(flux_column), epoch1, trial_p) for lc in lcc]
        )
        if len(trial_ttimes) <= 3:  # assuming epoch1, epoch2, epoch3 is always in the trial_ttimes.
            compat_p_list.append(trial_p)
            if verbose:
                print(f"  Period {trial_p} compatible. {trial_ttimes}")
        else:
            if verbose:
                print(f"  Period {trial_p} is incompatible. Has unexpected transits at times {trial_ttimes}")
        num_cycles += 1
    return compat_p_list


class TransitTimeSpec(dict):

    def __init__(
        self,
        epoch: float = None,
        period: float = None,
        duration_hr: float = None,
        transit_depth_percent: float = None,
        sector: int = None,
        steps_to_show: list = None,
        surround_time: float = None,
        label: str = None,
        defaults: TransitTimeSpec = None,
    ):
        # core parameters
        self["epoch"] = epoch
        self["period"] = period
        self["duration_hr"] = duration_hr
        # depth is used occasionally, for typical plotting, it is not used
        self["transit_depth_percent"] = transit_depth_percent

        # used for plotting
        self["sector"] = sector
        self["steps_to_show"] = steps_to_show
        self["surround_time"] = surround_time
        self["label"] = label

        if defaults is None:
            defaults = {}
        self._defaults = defaults  # put it as a custom attribute

    def __getitem__(self, key):
        res = super().get(key)
        if res is None:
            res = self._defaults.get(key)
        return res

    def get(self, key, default=None):
        res = self.__getitem__(key)
        if res is None:
            res = default
        return res


class TransitTimeSpecList(list):
    def __init__(self, *tt_spec_dict_list, defaults={}):
        self._defaults = TransitTimeSpec(**defaults)
        for tt_spec_dict in tt_spec_dict_list:
            self.append(TransitTimeSpec(**tt_spec_dict, defaults=self._defaults))

    def _spec_property_values(self, property_name):
        return np.array([tt[property_name] for tt in self])

    #
    # The following properties return the specific transit parameters
    # in an array. Together they can be used to create a mask
    # for the transits using ``LightCurve.create_transit_mask()``
    #

    @property
    def epoch(self):
        return self._spec_property_values("epoch")

    @property
    def period(self):
        return self._spec_property_values("period")

    @property
    def duration_hr(self):
        return self._spec_property_values("duration_hr")

    @property
    def duration(self):
        return self.duration_hr / 24

    @property
    def transit_depth_percent(self):
        return self._spec_property_values("transit_depth_percent")

    @property
    def label(self):
        return self._spec_property_values("label")

    def to_table(self, columns=("label", "epoch", "duration_hr", "period")):
        """Convert the specs to an ``astropy.Table``"""
        data = [getattr(self, col) for col in columns]
        return Table(data, names=columns)


def add_sector_like_as_column(lcf_coll, warn_if_failed=False):
    """For each lightcurve, add a column indicating its sequence in a mission,
    e.g., sector for TESS.
    Useful when the lightcurves are stitched together.
    """
    # deduce the header indicating the sequence
    header_name = None
    if isinstance(lcf_coll[0], lk.TessLightCurve):
        header_name = "SECTOR"
    # TODO: handle Kepler / K2
    else:
        if warn_if_failed:
            warnings.warn("Cannot determine the sector-like header.")
        return lcf_coll

    def add_sector_like_to_lc(lc):
        lc = lc.copy()
        lc[header_name.lower()] = lc.meta.get(header_name, -1)
        return lc

    return lk.LightCurveCollection([add_sector_like_to_lc(lc) for lc in lcf_coll])


def stitch(lcf_coll, ignore_incompatible_column_warning=False, to_add_sector_like_as_column=False, **kwargs):
    """Wrapper over native stitch(), and tweak the metadata so that it behaves like a typical single-sector lightcurve."""

    def update_meta_if_exists_in(lc_src, keys):
        for key in keys:
            val = lc_src.meta.get(key, None)
            if val is not None:
                lc_stitched.meta[key] = val

    def safe_del_meta(key):
        if lc_stitched.meta.get(key, None) is not None:
            del lc_stitched.meta[key]

    if to_add_sector_like_as_column:
        lcf_coll = add_sector_like_as_column(lcf_coll, warn_if_failed=True)

    if ignore_incompatible_column_warning:
        with warnings.catch_warnings():
            # suppress useless warning. Use cases: stitching QLP lightcurves with SPOC lightcurves (sap_flux is incompatible)
            warnings.filterwarnings(
                "ignore",
                category=lk.LightkurveWarning,
                message="The following columns will be excluded from stitching because the column types are incompatible:.*",
            )
            lc_stitched = lcf_coll.stitch(**kwargs)
    else:
        lc_stitched = lcf_coll.stitch(**kwargs)

    # now update the metadata

    lc_stitched.meta["STITCHED"] = True

    # update observation start/stop dates
    update_meta_if_exists_in(lcf_coll[0], ("TSTART", "DATE-OBS"))
    update_meta_if_exists_in(lcf_coll[-1], ("TSTOP", "DATE-END"))

    # TODO: recalculate TELAPSE, LIVETIME, DEADC (which is LIVETIME / TELAPSE)

    safe_del_meta("FILENAME")  # don't associate it with a file anymore

    # record the sectors stitched and the associated metadata
    sector_list = [lc.meta.get("SECTOR") for lc in lcf_coll if lc.meta.get("SECTOR") is not None]
    meta_list = [lc.meta for lc in lcf_coll if lc.meta.get("SECTOR") is not None]
    if len(sector_list) > 0:
        lc_stitched.meta["SECTORS"] = sector_list
        meta_dict = dict()
        for sector, meta in zip(sector_list, meta_list):
            meta_dict[sector] = meta.copy()
        lc_stitched.meta["HEADERS_ORIGINAL"] = meta_dict

    return lc_stitched


def stitch_lc_dict(lc_dict: dict, source_colname="source", normalize=True) -> lk.LightCurve:
    lcc = []
    for k in lc_dict:
        lc = lc_dict[k].copy()
        lc[source_colname] = k
        if normalize:
            if lc.flux.unit is u.mag:
                lc = to_normalized_flux_from_mag(lc)
            else:
                lc = lc.normalize()
        lcc.append(lc)
    return stitch(lk.LightCurveCollection(lcc), ignore_incompatible_column_warning=True, corrector_func=lambda lc: lc)


def to_window_length_for_2min_cadence(length_day):
    """Helper for LightCurve.flatten().
    Return a `window_length` for the given number of days, assuming the data has 2-minute cadence."""
    return to_window_length_for_cadence(length_day * u.day, 2 * u.min)


def to_window_length_for_cadence(length, cadence):
    """Helper for LightCurve.flatten().
    Return a `window_length` for the given length and cadence.

    Parameters length and cadence should be ~~astropy.quantity.Quantity~~.
    If they are unitless number, they should be in the same unit.
    """
    res = math.floor(length / cadence)
    if res % 2 == 0:
        res += 1  # savgol_filter window length must be odd number
    return res


# detrend using spline
# Based on:  https://github.com/barentsen/kepler-athenaeum-tutorial/blob/master/how-to-find-a-planet-tutorial.ipynb
def flatten_with_spline_normalized(lc, return_trend=False, **kwargs):
    lc = lc.remove_nans()
    spline = UnivariateSpline(lc.time, lc.flux, **kwargs)
    trend = spline(lc.time)
    # detrended = lc.flux - trend
    detrended_relative = 100 * ((lc.flux / trend) - 1) + 100  # in percentage
    lc_flattened = lc.copy()
    lc_flattened.flux = detrended_relative
    lc_flattened.flux_unit = "percent"
    if not return_trend:
        return lc_flattened
    else:
        lc_trend = lc.copy()
        lc_trend.flux = trend
        return (lc_flattened, lc_trend)


def _lksl_statistics(ts):
    """Compute LKSL Statistics of the given (time-series) values.
    Useful to compare the noises in a folded lightcurve.

    Based on https://arxiv.org/pdf/1901.00009.pdf , equation 4.
    See https://www.aanda.org/articles/aa/pdf/2002/17/aa2208.pdf for more information.
    """
    ts = ts[~np.isnan(ts)]

    vector_length_sq_sum = np.square(np.diff(ts)).sum()
    if len(ts) > 2:
        # to fully utilize the data by including the vector length between the last and first measurement
        # section 2 of Clarke, 2002
        vector_length_sq_sum += np.square(ts[-1] - ts[0])

    diff_from_mean_sq_sum = np.square(ts - np.mean(ts)).sum()
    if diff_from_mean_sq_sum == 0:  # to avoid boundary case that'd cause division by zero
        diff_from_mean_sq_sum = 1e-10

    return (vector_length_sq_sum / diff_from_mean_sq_sum) * (len(ts) - 1) / (2 * len(ts))


def lksl_statistics(lc, column="flux"):
    return _lksl_statistics(lc[column].value)


#
# TODO: util to estimate transit depth
#


def estimate_snr(
    lc,
    signal_depth,
    signal_duration,
    num_signals,
    savgol_to_transit_window_ratio=4,
    cdpp_kwargs=None,
    return_diagnostics=False,
):
    """Estimate Signal-to-Noise Ratio (SNR) of the signals, e.g., transits.
    The estimate assumes:
    - there is no red noises (due to systematics, etc.)
    - noise level does not vary over time.

    References:
    - based on KSCI-19085-001: Planet Detection Metrics
      (section 3:  one-sigma depth function. the basic form quoted is used instead of one-sigma depth function)
      https://exoplanetarchive.ipac.caltech.edu/docs/KSCI-19085-001.pdf
    - Poster with more in-depth treatment on white noises, red noises, etc.
      https://mirasolinstitute.org/kaspar/publications/Window_Functions.pdf

    See also:
    - https://arxiv.org/abs/2506.04136
      on how TESS SPOC estimate SNR using a similar procedure (it uses a more robust CDPP estimate.)
    """

    if not isinstance(signal_duration, u.Quantity):
        signal_duration = signal_duration * u.hour
    if cdpp_kwargs is None:
        cdpp_kwargs = dict()

    cadence = estimate_cadence(lc, unit=u.min)
    transit_window = int(np.ceil((signal_duration / cadence).decompose()))

    savgol_window = transit_window * savgol_to_transit_window_ratio
    if savgol_window % 2 == 0:
        savgol_window += 1

    if cdpp_kwargs.get("transit_duration") is None:
        cdpp_kwargs["transit_duration"] = transit_window
    if cdpp_kwargs.get("savgol_window") is None:
        cdpp_kwargs["savgol_window"] = savgol_window

    cdpp = lc.estimate_cdpp(**cdpp_kwargs)

    snr = np.sqrt(num_signals) * (signal_depth / cdpp).decompose().value

    if not return_diagnostics:
        return snr
    else:
        diagnostics = cdpp_kwargs.copy()
        diagnostics["cdpp"] = cdpp
        diagnostics["cadence"] = cadence
        return snr, diagnostics


def estimate_b(depth, t_full, t_total):
    # equation 1.8 of https://www.astro.ex.ac.uk/people/alapini/Publications/PhD_chap1.pdf,
    # which quotes  Seager & Mallén-Ornelas (2003)
    # limitations include: limb darkening is not taken into account. Circular orbit is assumed
    d, tF, tT = depth, t_full, t_total
    return (((1 - d**0.5) ** 2 - (tF / tT) ** 2 * (1 + d**0.5) ** 2) / (1 - (tF / tT) ** 2)) ** 0.5


def estimate_transit_duration_for_circular_orbit(period, rho, b):
    """Estimate the transit duration for circular orbit, T_circ.
    Usage: if the observed transit duration, T_obs,
    Case 1. T_obs < T_circ significantly: the planet candidate is potentially
    - Case 1a:
      - in a highly eccentric orbit, and
      - the transit occurs near periastron (nearest point to the host star)
        - the planet moves faster, thus the duration is shorter, figure 1 orange line
    - Case 1b:
      - actual impact parameter b is higher than the estimate
        - higher b implies the planet traverses across less of the host star's cross section
    - the uncertainty is manifestation of e-w-b degeneracy.

    Case 2. T_obs > T_circ significantly: the planet candidate is almost definitely
    - in a highly eccentric orbit, and
    - the transit occurs near apastron (farthest point to the host star)
      - the planet moves slower, thus the duration is longer, figure 1 green line

    cf. The TESS–Keck Survey. VI. Two Eccentric Sub-Neptunes Orbiting HIP-97166, MacDougall et. al
    https://ui.adsabs.harvard.edu/abs/2021AJ....162..265M/abstract

    - implementing the equations in section 2.2
    - see figure 1 for the concept, and section 2.2 for interpreting the result
    """
    if isinstance(b, (list, tuple, np.ndarray)):
        if len(b) == 2:  # in (mean, error) form
            b_mean = b[0]
            b_lower = b[0] - b[1]
            b_upper = b[0] + b[1]
        elif len(b) == 3:  # in (mean, lower error, upper error) form
            b_mean = b[0]
            b_lower = b[0] - b[1]
            b_upper = b[0] + b[2]
        else:
            raise ValueError("b must be a scalar, (mean, error), or (mean, lower_error, upper_error)")
        d_mean = estimate_transit_duration_for_circular_orbit(period, rho, b_mean)
        # lower b would lead to higher duration
        d_upper = estimate_transit_duration_for_circular_orbit(period, rho, b_lower)
        d_lower = estimate_transit_duration_for_circular_orbit(period, rho, b_upper)

        return np.array([d_mean.value, d_lower.value, d_upper.value]) * d_mean.unit

    # case b is a single scalar, do the actual calc

    # use default units if the input is not quantity
    if not isinstance(period, u.Quantity):
        period = period * u.day

    if not isinstance(rho, u.Quantity):
        rho = rho * u.gram / u.cm**3

    # implement eq. 1, 2, and 3
    return (
        period ** (1 / 3)
        * rho ** (-1 / 3)
        * (1 - b**2) ** (1 / 2)
        *
        # constants that are not explicitly stated in eq 3, but can be derived from eq 1 and 2
        np.pi ** (-2 / 3)
        * 3 ** (1 / 3)
        * astropy.constants.G ** (-1 / 3)
    ).to(u.hour)


def estimate_period_for_circular_orbit(duration, rho, b):
    """Estimate the period for circular orbit, P
    https://ui.adsabs.harvard.edu/abs/2021AJ....162..265M/abstract

    - implementing the equations in section 2.2
    """
    if isinstance(b, (list, tuple, np.ndarray)):
        if len(b) == 2:  # in (mean, error) form
            b_mean = b[0]
            b_lower = b[0] - b[1]
            b_upper = b[0] + b[1]
        elif len(b) == 3:  # in (mean, lower error, upper error) form
            b_mean = b[0]
            b_lower = b[0] - b[1]
            b_upper = b[0] + b[2]
        else:
            raise ValueError("b must be a scalar, (mean, error), or (mean, lower_error, upper_error)")
        p_mean = estimate_period_for_circular_orbit(duration, rho, b_mean)
        # lower b would lead to shorter duration
        p_lower = estimate_period_for_circular_orbit(duration, rho, b_lower)
        p_upper = estimate_period_for_circular_orbit(duration, rho, b_upper)

        return np.array([p_mean.value, p_lower.value, p_upper.value]) * p_mean.unit

    # case b is a single scalar, do the actual calc

    # use default units if the input is not quantity
    if not isinstance(duration, u.Quantity):
        duration = duration * u.hour

    if not isinstance(rho, u.Quantity):
        rho = rho * u.gram / u.cm**3

    # implement eq. 1, 2, and 3
    #
    # Note: Pyaneti single transit mode uses a slightly different method,
    # based on https://academic.oup.com/mnras/article/457/3/2273/2588921
    # it  takes into the account of Rp/R* (related to transit depth)
    #
    # The method used here does not take Rp/R* into consideration.
    # It should be an approximation that assumes Rp/R* is 0.
    #
    # Pyaneti source:
    # https://github.com/oscaribv/pyaneti/blob/241bb09931737c522dd000e1309bd2c8bfe0e7ba/src/print_values.py#L294-L302
    #
    return (
        (
            duration
            / (
                rho ** (-1 / 3)
                * (1 - b**2) ** (1 / 2)
                *
                # constants that are not explicitly stated in eq 3, but can be derived from eq 1 and 2
                np.pi ** (-2 / 3)
                * 3 ** (1 / 3)
                * astropy.constants.G ** (-1 / 3)
            )
        )
        ** 3
    ).to(u.day)


def _calc_median_flux_around(lc_f: lk.FoldedLightCurve, epoch_phase, flux_window_in_min):
    flux_window = flux_window_in_min / 60 / 24  # in days
    lc_trunc = lc_f.truncate(epoch_phase - flux_window / 2, epoch_phase + flux_window / 2).remove_nans()
    flux_median = np.median(lc_trunc.flux)
    flux_median_sample_size = len(lc_trunc)

    return flux_median, flux_median_sample_size


def calc_flux_at_minimum(lc_f: lk.FoldedLightCurve, flux_window_in_min=10):
    """Return the flux at minimum by calculating the median of the flux at minimum, assumed to be at phase 0"""
    return _calc_median_flux_around(lc_f, 0, flux_window_in_min=flux_window_in_min)


def calc_peak_to_peak(lc_f: lk.FoldedLightCurve, flux_window_in_min=10):
    """Return the flux at minimum by calculating the median of the flux at minimum"""

    argmax_func, argmin_func = "argmax", "argmin"
    if lc_f.flux.unit == u.mag:
        # if the unit is magnitude, reverse max/min
        argmax_func, argmin_func = "argmin", "argmax"

    time_max = lc_f.time[getattr(lc_f.flux, argmax_func)()]
    flux_max, flux_max_sample_size = _calc_median_flux_around(lc_f, time_max.value, flux_window_in_min=flux_window_in_min)

    time_min = lc_f.time[getattr(lc_f.flux, argmin_func)()]
    flux_min, flux_min_sample_size = _calc_median_flux_around(lc_f, time_min.value, flux_window_in_min=flux_window_in_min)

    peak_to_peak = np.abs(flux_max - flux_min)

    return SimpleNamespace(
        peak_to_peak=peak_to_peak,
        time_max=time_max,
        flux_max=flux_max,
        flux_max_sample_size=flux_max_sample_size,
        time_min=time_min,
        flux_min=flux_min,
        flux_min_sample_size=flux_min_sample_size,
    )


def select_flux(lc, flux_cols):
    """Return a Lightcurve object with the named column as the flux column.

    flux_cols: either a column name (string), or a list of prioritized column names
    such that the first one that the lightcurve contains will be used.
    """

    def _to_lc_with_1flux(lc, flux_1col):
        flux_1col = flux_1col.lower()
        if "flux" == flux_1col:
            return lc
        elif flux_1col in lc.colnames:
            return lc.select_flux(flux_1col)
        else:
            return None

    if isinstance(flux_cols, str):
        flux_cols = [flux_cols]

    for flux_1col in flux_cols:
        res = _to_lc_with_1flux(lc, flux_1col)
        if res is not None:
            return res
    raise ValueError(f"'column {flux_cols}' not found")


def correct_crowding(lc_sap: lk.LightCurve, crowdsap=None, flfrcsap=None):
    """Create crowding corrected lightcurve based on SAP_FLUX lightcurve.
    # based on / see: https://heasarc.gsfc.nasa.gov/docs/tess/UnderstandingCrowding.html
    """
    # warn users if lc_sap is likely to be non-SAP
    flux_origin = lc_sap.meta.get("FLUX_ORIGIN")
    if flux_origin is not None and "sap_flux" != flux_origin:
        warnings.warn(
            f"correct_crowding(): supplied lightcurve is likely not SAP based, with FLUX_ORIGIN: {flux_origin}. "
            "Corrected lightcurve returned might not be valid."
        )

    if crowdsap is None:
        crowdsap = lc_sap.meta["CROWDSAP"]
    if flfrcsap is None:
        flfrcsap = lc_sap.meta["FLFRCSAP"]

    median_flux = np.nanmedian(lc_sap.flux)
    excess_flux = (1 - crowdsap) * median_flux
    # note: a **constant** amount of flux is removed, so a dip would be proportionally deeper
    flux_excess_removed = lc_sap.flux - excess_flux
    flux_corr = flux_excess_removed / flfrcsap

    # Calculate the new uncertainties
    flux_err_corr = lc_sap.flux_err / flfrcsap

    # OPEN: should we copy the entire lc_sap instead?
    lc_corr = type(lc_sap)(time=lc_sap.time, flux=flux_corr, flux_err=flux_err_corr)
    lc_corr.meta.update(lc_sap.meta)
    lc_corr.meta["CROWDSAP"] = crowdsap
    lc_corr.meta["FLFRCSAP"] = flfrcsap
    return lc_corr


def deblend_mag(lc, contaminant_mag):
    """Deblend a blended lightcurve by subtracting the (constant)
    contaminant magnitude.

    The purpose is similar to ``correct_crowding``, but the approach
    here is used by VSX, in the form of a spreadsheet:
    https://www.aavso.org/vsx/_images/COMBLEND.xlsx
    listed in the FAQ:
    https://www.aavso.org/vsx/index.php?view=about.faq
    """

    if lc.flux.unit != u.mag:
        raise ValueError(f"Given lightcurve's flux must be in mag. Actual {lc.flux.unit}")

    blended_mag = lc.flux.value
    target_mag = contaminant_mag - 2.5 * np.log10(10 ** ((contaminant_mag - blended_mag) * 0.4) - 1)

    target_lc = lk.LightCurve(time=lc.time, flux=target_mag * u.mag)
    target_lc.meta.update(lc.meta)
    return target_lc


def combine_magnitudes(mag1, mag2):
    """Return the combined magnitude of 2 stars."""
    # https://www.astro.keele.ac.uk/jkt/pubs/JKTeq-fluxsum.pdf
    return -2.5 * np.log10(10 ** (-0.4 * mag1) + 10 ** (-0.4 * mag2))


def normalized_flux_val_to_mag(flux_val, base_mag):
    return base_mag + 2.5 * np.log10(1 / flux_val)


def to_flux_in_mag_by_normalization(lc, base_mag_header_name=None, base_mag=None):
    """Convert the a lightcurve's flux to magnitude via a normalized lightcurve with a known average / base magnitude."""
    if lc.flux.unit is u.mag:
        return lc

    lc = lc.copy()

    if base_mag is not None:
        pass
    elif base_mag_header_name is not None:
        base_mag = lc.meta.get(base_mag_header_name)
    else:
        base_mag = lc.meta.get("TESSMAG", lc.meta.get("KEPMAG"))  # Try TESS or Kepler as base
    if base_mag is None:
        raise ValueError(f"The given lightcurve does not have base magnitude in {base_mag_header_name} header ")

    lc_norm = lc.normalize()
    flux_mag = (base_mag + 2.5 * np.log10(1 / lc_norm.flux)) * u.mag
    flux_err_mag = (1.086 * lc_norm.flux_err / lc_norm.flux) * u.mag
    lc.flux = flux_mag
    lc.flux_err = flux_err_mag
    lc.meta["NORMALIZED"] = False
    return lc


def to_normalized_flux_from_mag(lc):
    """Convert the a lightcurve's flux from magnitude to normalized flux."""
    if lc.flux.unit is not u.mag:
        raise ValueError("The flux must be in magnitude")

    lc = lc.copy()

    median_flux_mag = np.nanmedian(lc.flux.value)

    flux_delta_mag = lc.flux.value - median_flux_mag

    flux_norm = 1 / np.power(10, flux_delta_mag / 2.5)
    flux_err_norm = np.abs(1 - 1 / np.power(10, lc.flux_err.value / 2.5))

    lc.flux = flux_norm * u.dimensionless_unscaled
    lc.flux_err = flux_err_norm * u.dimensionless_unscaled
    lc.meta["NORMALIZED"] = True
    return lc


def to_normalized_flux_from_mag_vals(vals, errs):
    """Convert values from magnitude to normalized values."""
    if vals.unit is not u.mag:
        raise ValueError("The values must be in magnitude")

    median_flux_mag = np.nanmedian(vals.value)

    flux_delta_mag = vals.value - median_flux_mag
    vals_norm = 1 / np.power(10, flux_delta_mag / 2.5)

    if errs is not None:
        errs_norm = np.abs(1 - 1 / np.power(10, errs.value / 2.5))
    else:
        errs_norm = None

    return vals_norm, errs_norm


def ratio_to_mag(val_in_ratio):
    """Convert normalized transit depth to magnitude."""
    return 2.5 * np.log10(1 / (1 - val_in_ratio))


def to_hjd_utc(t_obj: Time, sky_coord: SkyCoord) -> Time:
    # Based on astropy documentation
    # https://docs.astropy.org/en/stable/time/#barycentric-and-heliocentric-light-travel-time-corrections

    t_jd_tdb = t_obj.copy("jd").tdb

    # 1. convert the given time in JD TDB to local time UTC
    greenwich = coord.EarthLocation.of_site("greenwich")
    ltt_bary = t_jd_tdb.light_travel_time(sky_coord, location=greenwich, kind="barycentric")
    t_local_jd_utc = (t_jd_tdb - ltt_bary).utc

    # 2. convert local time UTC to HJD UTC
    ltt_helio = t_local_jd_utc.light_travel_time(sky_coord, location=greenwich, kind="heliocentric")
    t_hjd_utc = t_local_jd_utc + ltt_helio

    return t_hjd_utc


def convert_lc_time_to_hjd_utc(lc, target_coord, cache_dir=".", cache_key=None, cache_key_suffix=None):
    """Convert the lightcurve's time to HJD UTC, with the HJD time cached in disk.

    `cache_key_suffix`: optional, applicable only when `cache_key` is `None`. It lets user to
    identify different variants of the same LC (which may have slightly different sets of cadence),
     e.g., one with outliers removed, detrended, etc.
    """

    # a reasonable default for TESS / Kepler / K2 based lc
    def get_cache_key(lc):
        target_str = lc.meta.get("LABEL", "").replace(" ", "_")
        author = lc.meta.get("AUTHOR", None)
        if author is not None:
            target_str = f"{target_str}_{author}"
        sectors_str = "_".join([str(s) for s in lc.meta.get("SECTORS", [])])
        if sectors_str == "":
            sectors_str = lc.meta.get("SECTOR", "na")
        if cache_key_suffix is None:
            return f"hjd_{target_str}_{sectors_str}.txt"
        else:
            return f"hjd_{target_str}_{sectors_str}_{cache_key_suffix}.txt"

    def time_in_hjd_utc(lc, target_coord, cache_dir):
        # actual time conversion to HJD UTC
        # it can be take a while (close to 1 minute for 160K)
        # so we support caching the result persistently
        hjd_time_val = None
        nonlocal cache_key
        if cache_key is None:
            cache_key = get_cache_key(lc)
        cache_file = f"{cache_dir}/hjd/{cache_key}"
        if (os.path.exists(cache_file)) and (os.path.getsize(cache_file) > 0):
            hjd_time_val = np.genfromtxt(cache_file)
            #  ensure the cached time is compatible with the input LC
            if len(lc) == len(hjd_time_val):
                return hjd_time_val
            else:
                print(f"The cached HJD time in {cache_file} has different length. discard it.")
                hjd_time_val = None

        # else cache miss
        if not isinstance(target_coord, SkyCoord):
            target_coord = SkyCoord(target_coord["ra"], target_coord["dec"], unit=(u.deg, u.deg), frame="icrs")
        hjd_time_val = to_hjd_utc(lc.time, target_coord).value

        # write to cache
        os.makedirs(os.path.dirname(cache_file), exist_ok=True)
        np.savetxt(cache_file, hjd_time_val, fmt="%f", header="hjd")
        return hjd_time_val

    # the main logic
    hjd_time_val = time_in_hjd_utc(lc, target_coord, cache_dir)

    # create a new LC object with the HJD time.
    # It does not modify the time of an existing LC object, because it is extremely slow.
    data = QTable(data=lc)
    data.remove_column("time")
    lc_hjd = lc.__class__(time=Time(hjd_time_val, format="jd", scale="utc"), data=data)
    lc_hjd.meta.update(lc.meta)

    return lc_hjd


# fast_nanmean(): an optimized nanmean for binning
# https://github.com/astropy/astropy/pull/17574
# (planned for astropy 7.1.0)


def nanmean_reduceat(data, indices):
    mask = np.isnan(data)

    if mask.any():  # If there are NaNs
        # Create a writeable copy and mask NaNs
        data = data.copy()
        data[mask] = 0
        count_data = np.add.reduceat(~mask, indices)
        # Avoid division by zero warnings
        count_data = count_data.astype(data.dtype)
        count_data[count_data == 0] = np.nan
    else:
        # Derive counts from indices
        count_data = np.diff(indices, append=len(data))
        count_data[count_data <= 0] = 1

    sum_data = np.add.reduceat(data, indices)
    return sum_data / count_data


def fast_nanmean(a, **kwargs):
    return np.nanmean(a, **kwargs)


fast_nanmean.reduceat = nanmean_reduceat


HAS_BOTTLENECK = False
try:
    import bottleneck

    HAS_BOTTLENECK = True
except:
    HAS_BOTTLENECK = False


def parse_aggregate_func(cenfunc):
    # based on Astropy SigmaClip
    # https://github.com/astropy/astropy/blob/326435449ad8d859f1abf36800c3fb88d49c27ea/astropy/stats/sigma_clipping.py#L263

    # cenfunc should really be aggfunc, but I keep the name from SigmaClip
    if cenfunc is None:
        cenfunc = fast_nanmean

    if isinstance(cenfunc, str):
        if cenfunc == "median":
            if HAS_BOTTLENECK:
                cenfunc = bottleneck.nanmedian  # SigmaClip has more robust version
            else:
                cenfunc = np.nanmedian  # pragma: no cover

        elif cenfunc == "mean":
            if HAS_BOTTLENECK:
                cenfunc = bottleneck.nanmean  # SigmaClip has more robust version
            else:
                cenfunc = np.nanmean  # pragma: no cover

        else:
            raise ValueError(f"{cenfunc} is an invalid cenfunc.")

    return cenfunc


def bin_flux(lc, columns=["flux", "flux_err"], **kwargs):
    """Helper to bin() more efficiently."""
    # Note: the biggest slowdown comes from astropy regression
    # that this impl cannot address:
    # https://github.com/astropy/astropy/issues/13058

    aggregate_func = parse_aggregate_func(kwargs.get("aggregate_func"))
    kwargs["aggregate_func"] = aggregate_func

    # construct a lc_subset that only has a subset of columns,
    # to minimize the number of columns that need to be binned
    # see: https://github.com/lightkurve/lightkurve/issues/1191

    # lc_subset = lc['time', 'flux', 'flux_err'] does not work
    # due to https://github.com/lightkurve/lightkurve/issues/1194
    lc_subset = type(lc)(time=lc.time.copy())
    lc_subset.meta.update(lc.meta)
    for c in columns:
        if c in lc.colnames:
            col = lc[c]
            # convert astropy Masked to regular Column / Quantity
            # needed for the default custom fast_nanmean above
            # see: https://github.com/astropy/astropy/pull/17875
            # OPEN: consider to always convert to Column / Quantity
            # even without using the fast_nanmean, binning with regular Column / Quantity
            # is several times faster than with Masked (in astropy 6.0.1)
            if aggregate_func is fast_nanmean and isinstance(col, astropy.utils.masked.Masked):
                col = col.filled(np.nan)
            lc_subset[c] = col
        else:
            warnings.warn(f"bin_flux(): column {c} cannot be found in lightcurve. It is ignored.")

    return lc_subset.bin(**kwargs)


def normalize(lc, **kwargs):
    if "TASOC" == lc.meta.get("AUTHOR") and str(lc.flux.unit) == "ppm":
        # convert TASOC's flux (e.g. corr_flux), in ppm and median at 0,
        # to the more typical unscaled one with median at 1
        # so that it works nice with lightkurve's normalize()
        lc = lc.copy()
        lc.flux = lc.flux.to(u.dimensionless_unscaled) + 1
        lc.flux_err = lc.flux_err.to(u.dimensionless_unscaled)

    return lc.normalize(**kwargs)


def abbrev_sector_list(lcc_or_sectors_or_lc_or_sr):
    """Abbreviate a list of sectors, e.g., `1,2,3, 9` becomes `1-3, 9`."""

    # OPEN: consider to handle Kepler quarter / K2 campaign.

    sectors = lcc_or_sectors_or_lc_or_sr
    if sectors is None:
        sectors = []
    elif isinstance(sectors, lk.collections.Collection):  # LC / TPF collection
        sectors = [lc.meta.get("SECTOR") for lc in lcc_or_sectors_or_lc_or_sr]
    elif isinstance(sectors, lk.LightCurve):
        sectors = sectors.meta.get("SECTORS", [sectors.meta.get("SECTOR")])  # case my custom stitched lc that has SECTORS meta
    elif isinstance(sectors, lk.SearchResult):
        sectors = [s for s in sectors.table["sequence_number"]]
    # else it's assumed to be a list of sector number

    sectors = sectors.copy()
    sectors.sort()

    if len(sectors) < 1:
        return ""

    def a_range_to_str(range):
        if range[0] == range[1]:
            return str(range[0])
        else:
            return f"{range[0]}-{range[1]}"

    ranges = []
    cur_range = [sectors[0], sectors[0]]
    for s in sectors[1:]:
        if s == cur_range[1] + 1:
            cur_range = [cur_range[0], s]
        else:
            ranges.append(cur_range)
            cur_range = [s, s]
    ranges.append(cur_range)

    return ", ".join([a_range_to_str(r) for r in ranges])


#
# Misc. generic astronomy helpers
#

RHO_sun = astropy.constants.M_sun / (4 / 3 * np.pi * astropy.constants.R_sun**3)


def estimate_rho(mass, radius, return_unit=None):
    """Estimate density based on mass and radius.
    mass / radius is assumed to be in the units of M_sun, R_sun respectively,
    if they are not of `Quantity` type.
    """
    if not isinstance(mass, u.Quantity):
        mass = mass * astropy.constants.M_sun
    if not isinstance(radius, u.Quantity):
        radius = radius * astropy.constants.R_sun

    rho = mass / (4 / 3 * np.pi * radius**3)
    if return_unit is None:  # return in the unit of solar density
        rho = rho / RHO_sun
    else:
        rho = rho.to(return_unit)
    return rho


#
# TargetPixelFile helpers
#


def truncate(tpf: lk.targetpixelfile.TargetPixelFile, before: float, after: float) -> lk.targetpixelfile.TargetPixelFile:
    return tpf[(tpf.time.value >= before) & (tpf.time.value <= after)]


def to_lightcurve_with_custom_aperture(tpf, aperture_mask, background_mask):
    """Create a lightcurve from the given TargetPixelFile, with suppplied aperture and background"""
    n_aperture_pixels = aperture_mask.sum()
    aperture_lc = tpf.to_lightcurve(aperture_mask=aperture_mask)  # aperture + background

    n_background_pixels = background_mask.sum()
    background_lc_per_pixel = tpf.to_lightcurve(aperture_mask=background_mask) / n_background_pixels
    background_lc = background_lc_per_pixel * n_aperture_pixels

    corrected_lc = aperture_lc - background_lc
    # OPEN: consider filling in lc metadata, label, etc.

    return corrected_lc, aperture_lc, background_lc


def to_mask_in_pixel_coordinate(lc_or_tpf, mask=None):
    """Convert the given aperture mask to the form of CCD pixel coordinate array.
    The CDD pixel coordinate array form is used by some programs, such as triceratops.
    """

    def get_mask_coord_ref_tpf(tpf, mask):
        if mask is None:  # default to pipeline mask
            mask = tpf.pipeline_mask

        row_base, col_base = tpf.row, tpf.column
        return mask, row_base, col_base

    def get_mask_coord_ref_lc(lc, mask):
        # TODO: the logic probably only works for TESS SPOC LC fits, maybe Kepler or QLP too.
        with fits.open(lc.filename) as hdul:
            row_base, col_base = hdul[2].header.get("CRVAL2P "), hdul[2].header.get("CRVAL1P ")
            if mask is None:
                # the data encodes the aperture pixel, background pixels, etc.
                # convert it into aperture mask
                pixels = hdul[2].data
                mask = pixels == np.max(pixels)

        return mask, row_base, col_base

    if isinstance(lc_or_tpf, lk.targetpixelfile.TargetPixelFile):
        mask, row_base, col_base = get_mask_coord_ref_tpf(lc_or_tpf, mask)
    elif isinstance(lc_or_tpf, lk.LightCurve):
        mask, row_base, col_base = get_mask_coord_ref_lc(lc_or_tpf, mask)
    else:
        raise ValueError(f"Only LightCurve or TPF is supported. {type(lc_or_tpf)}")

    num_rows, num_cols = mask.shape
    res = []
    for y in range(0, num_rows):
        for x in range(0, num_cols):
            if mask[y][x]:
                res.append([col_base + x, row_base + y])
    return np.array(res)


#
# Astropy extension
#


# Specify / display time delta in hours
# useful for specifying / displaying planet transits
class TimeDeltaHour(astropy.time.TimeDeltaNumeric):
    """Time delta in hours (3600 SI seconds)"""

    import erfa

    name = "hour"
    unit = 3600.0 / erfa.DAYSEC  # for quantity input


#
# Others
#


def coordinate_like_id_to_coordinate(id, style="decimal"):
    """Given an identifier that describes the target's coordinate, return the coordInate in string or `SkyCoord` object.
    The type of IDs supported are in the form of '<prefix-of-catalog> J<ra><dec>'.
    Catalogs that adopts such form of ids include PTF1, ASASSN-V, WISE, WISEA, 1SWASP, etc.

    """
    # e.g.,
    # PTF1 J2219+3135 (the shortened form often found in papers)
    # PTF1 J221910.09+313523.1  (the actual id, with higher precision)
    result = re.findall(r"^.{2,}\s*J(\d{2})(\d{2})([0-9.]*)([+-])(\d{2})(\d{2})([0-9.]*)\s*$", id)
    if len(result) < 1:
        return None
    [ra_hh, ra_mm, ra_ss, dec_sign, dec_deg, dec_min, dec_ss] = result[0]
    coord = SkyCoord(f"{ra_hh} {ra_mm} {ra_ss} {dec_sign} {dec_deg} {dec_min} {dec_ss}", unit=(u.hourangle, u.deg))
    if style is None:
        return coord
    else:
        return coord.to_string(style=style)


from astropy.coordinates import SkyCoord, Angle
from astroquery.vizier import Vizier


def calc_separation(
    target_coord, result, ra_colname="RAJ2000", dec_colname="DEJ2000", sep_colname="_r", pos_angle_colname="_p"
):
    """Calculate angular separation from the target coordinate.
    If `target_coord` == `"first_row"`, use the first row as the target coordinate.
    """
    if len(result) < 1:
        return result

    coord_unit = result[ra_colname].unit or u.deg
    if target_coord == "first_row":
        target_coord = SkyCoord(result[0][ra_colname], result[0][dec_colname], unit=coord_unit)
    result_coord = SkyCoord(result[ra_colname], result[dec_colname], unit=coord_unit)
    separation = target_coord.separation(result_coord)
    result[sep_colname] = separation.arcsec
    result[sep_colname].unit = u.arcsec
    result[sep_colname].info.format = "{:6.3f}"
    pos_angle = target_coord.position_angle(result_coord)
    result[pos_angle_colname] = pos_angle.deg
    result[pos_angle_colname].unit = u.deg
    result[pos_angle_colname].info.format = "{:5.1f}"

    return result


def search_nearby(
    ra,
    dec,
    equinox="J2000.0",
    catalog_name="I/355/gaiadr3",
    radius_arcsec=60,
    magnitude_limit_column=None,
    magnitude_lower_limit=None,
    magnitude_upper_limit=None,
    pmra=None,
    pmdec=None,
    pmra_lower=None,
    pmra_upper=None,
    pmdec_lower=None,
    pmdec_upper=None,
    pmra_limit_column="pmRA",
    pmdec_limit_column="pmDE",
    calc_separation_from_first_row=False,
    include_gaiadr3_astrophysical=False,
    warn_if_all_filtered=True,  # warn if PM/Mag filter filters outs all query result
):
    """Stars around the given coordinate from Gaia DR2/EDR3, etc."""

    c1 = SkyCoord(ra, dec, equinox=equinox, frame="icrs", unit="deg")
    Vizier.ROW_LIMIT = -1
    # include and sorted by angular separation. "_r" is the generic Vizier-calculated column for it.
    # "_p": Position angle (E of N), also a generic Vizier-calculated column.
    columns = ["*", "+_r", "_p"]
    if catalog_name == "I/350/gaiaedr3":
        columns += ["epsi", "sepsi"]  # add astrometric excess noise to the output (see if a star wobbles)
    if catalog_name == "I/355/gaiadr3":
        columns += ["epsi", "sepsi", "VarFlag", "IPDfmp", "Dup", "GRVSmag"]  # also add variability and Duplicate flag

    if catalog_name == "I/355/gaiadr3" and include_gaiadr3_astrophysical:
        catalog_names_in_query = ["I/355/gaiadr3", "I/355/paramp"]
        columns += [
            "SpType-ELS",  # spectral type
            # emission line star class (beStar, etc). See CLASSLABEL_ESPELS, 20.2.1 astrophysical_parameters, Gaia DR3 doc
            "ClassELS",
            "f_ClassELS",
            "Evol",  # stellar evolution stage
            "Flags-Flame",
            # include errors for mass / radius
            "b_Rad",
            "B_Rad",
            "b_Rad-Flame",
            "B_Rad-Flame",
            "b_Mass-Flame",
            "B_Mass-Flame",
        ]
    else:
        catalog_names_in_query = catalog_name

    with warnings.catch_warnings():
        # suppress useless warning.  https://github.com/astropy/astroquery/issues/2352
        warnings.filterwarnings(
            "ignore", category=astropy.units.UnitsWarning, message="Unit 'e' not supported by the VOUnit standard"
        )
        result_all = Vizier(columns=columns).query_region(
            c1,
            catalog=catalog_names_in_query,
            radius=Angle(radius_arcsec, "arcsec"),
        )
    if len(result_all) < 1:  # handle no search result case
        return None
    result = result_all[catalog_name]

    # Convert Gaia DR3 mag to Vmag
    if catalog_name in ["I/355/gaiadr3", "I/350/gaiaedr3"] and all([c in result.columns for c in ["Gmag", "BP-RP"]]):
        result["Vmag"] = gaia_dr3_mag_to_vmag(result["Gmag"], result["BP-RP"])
        result["Vmag"].info.unit = result["Gmag"].info.unit
        result["Vmag"].info.description = "V-band mean magnitude, derived from Gmag and BP-RP."

    result_pre_filter = result

    if magnitude_lower_limit is not None:
        result = result[(magnitude_lower_limit <= result[magnitude_limit_column]) | result[magnitude_limit_column].mask]

    if magnitude_upper_limit is not None:
        result = result[(result[magnitude_limit_column] <= magnitude_upper_limit) | result[magnitude_limit_column].mask]

    if pmra is not None and pmra_upper is not None and pmra_lower is not None:
        result = result[(pmra_lower <= result[pmra_limit_column]) | result[pmra_limit_column].mask]
        result = result[(result[pmra_limit_column] <= pmra_upper) | result[pmra_limit_column].mask]
    if pmdec is not None and pmdec_upper is not None and pmdec_lower is not None:
        result = result[(pmdec_lower <= result[pmdec_limit_column]) | result[pmdec_limit_column].mask]
        result = result[(result[pmdec_limit_column] <= pmdec_upper) | result[pmdec_limit_column].mask]

    # tweak default format to make magnitudes and separation more succinct
    for col in ["separation", "RPmag", "Gmag", "BPmag", "BP-RP", "GRVSmag", "Vmag"]:
        if col in result.colnames:
            result[col].info.format = ".3f"
    # fill in the missing unit (Vizier / astroquery does not provide)
    result["_r"].unit = u.arcsec
    result["_p"].unit = u.deg
    if calc_separation_from_first_row:
        result = calc_separation("first_row", result, ra_colname="RA_ICRS", dec_colname="DE_ICRS")

    if warn_if_all_filtered and len(result) == 0 and len(result_pre_filter) > 0:
        warnings.warn(f"All query results filtered due to mag/PM filter. Num. of entries pre-filter: {len(result_pre_filter)}")

    if catalog_name == "I/355/gaiadr3" and include_gaiadr3_astrophysical:
        result_paramp = result_all["I/355/paramp"]

        # only Sources in the filtered main result
        result_paramp = result_paramp[np.isin(result_paramp["Source"], result["Source"])]

        for col in ["Pstar", "Pbin", "PWD"]:  # tweak default format
            if col in result_paramp.colnames:
                result_paramp[col].info.format = ".2f"
        # fill in the missing unit (Vizier / astroquery does not provide)
        result_paramp["_r"].unit = u.arcsec
        result_paramp["_p"].unit = u.deg
        # Note for case calc_separation_from_first_row is True
        # for now we don't re-calculate it, as getting the result consistent is a bit tricky
        # Reason:
        # - result_paramp do not have coordiante columns of RAJ2000, DEJ2000
        #   one could use functionally identical _RAJ2000, _DEJ2000 (calculated by Vizier),
        #   but the values are not identical.
        #   So users might sitll still small discrepancy at times.
        #   To keep things simple, we don't do the recalcuation for now

        return result, result_paramp
    else:
        return result


def gaia_dr3_mag_to_vmag(gmag, b_minus_r):
    """
    Convert Gaia DR3 magnitude (`Gmag` and `BP-RP`) to Johnson V.

    Based on Table 5.9 in Gaia DR3 Documentation, Section 5.5.1 Photometric relationships with other photometric systems:
    https://gea.esac.esa.int/archive/documentation/GDR3/Data_processing/chap_cu5pho/cu5pho_sec_photSystem/cu5pho_ssec_photRelations.html

    The applicable range of `BP-RP` is from Table 5.10 in the documentation.

    Applicable to Gaia EDR3 data as well: The formula is identical.
    """
    if isinstance(gmag, list):
        gmag = np.asarray(gmag)
    g_minus_v = -0.02704 + 0.01424 * b_minus_r - 0.2156 * b_minus_r**2 + 0.01426 * b_minus_r**3

    # check BP-RP is in appiicale range,
    # we check them one by one so that teh warning could be specific.
    if isinstance(b_minus_r, (float, int)):
        b_minus_r_ary = [b_minus_r]
    else:
        b_minus_r_ary = b_minus_r
    for a_b_r in b_minus_r_ary:
        if not (-0.5 < a_b_r < 5.0):
            warnings.warn(
                f"gaia_dr3_mag_to_vmag(): b_minus_r value ({a_b_r}) is outside the applicable range for the transformation. "
                "The result is probably not reliable."
            )

    return gmag - g_minus_v


def are_stars_bound_by_plx_pm(
    rs,
    plx_diff_sig_threshold=3,
    cpm_index_threshold=10,
    parllax_col_name="Plx",
    parllax_err_col_name="e_Plx",
    pmra_col_name="pmRA",
    pmde_col_name="pmDE",
    verbose=False,
):
    """Deterine if two stars are bound by parllax and proper motion.
    The two stars are the first two rows of the table,
    with the default table format / expected columns from Gaia DR3 Main in Vizier (I/355).

    The methdology and the default cutoff thresholds are based on 2025arXiv250501470M (Section 3)
    https://ui.adsabs.harvard.edu/abs/2025arXiv250501470M/abstract

    """

    # Note: in the paper search radius for potential companion is dependent on target parllax,
    # so that separation is <= 10000 AU
    #   - eq 1: r [arcsec] = 10 * pi * plx [mas]

    def calc_plx_diff(rs, parllax_col_name="Plx", parllax_err_col_name="e_Plx"):
        # OPEN: the paper suggests epsi is also used in considering signifiance but I don't know how.
        # e.g. for TOI-4661, the paper reports plx_diff_sig=2.0, while the procedure here reports 2.25
        #      it is as if epsi is somehow added to Gaia DR3's plx_diff_err (as the denominator)
        # implication: plx_diff_sig calculated here is possibly slightly larger than that of the paper (table 4)
        plx_diff = abs(rs[parllax_col_name][0] - rs[parllax_col_name][1])
        plx_diff_err = max(rs[parllax_err_col_name][0], rs[parllax_err_col_name][1])
        plx_diff_sig = plx_diff / plx_diff_err
        return plx_diff, plx_diff_err, plx_diff_sig

    def calc_pm_diff(rs, pmra_col_name="pmRA", pmde_col_name="pmDE"):
        # i.e., differential proper motion 𝜇~rel~ relative to the target (row 0). See paper Section 3
        pmra_d = rs[pmra_col_name][1] - rs[pmra_col_name][0]
        pmde_d = rs[pmde_col_name][1] - rs[pmde_col_name][0]
        pm_diff = np.sqrt(pmra_d**2 + pmde_d**2)
        pmra_err_col_name = "e_pmRA"
        pmde_err_col_name = "e_pmDE"
        pm0_err = np.sqrt(rs[pmra_err_col_name][0] ** 2 + rs[pmde_err_col_name][0] ** 2)
        pm1_err = np.sqrt(rs[pmra_err_col_name][1] ** 2 + rs[pmde_err_col_name][1] ** 2)
        pm_diff_err = (pm0_err + pm1_err) / 2
        return pm_diff, pm_diff_err

    def calc_pm_diff_n_cpm_index(rs, pmra_col_name="pmRA", pmde_col_name="pmDE"):
        # an index characterizing the degree of common proper motion. eq (2) in the paper
        # combined PM
        pmra_c = rs[pmra_col_name][1] + rs[pmra_col_name][0]
        pmde_c = rs[pmde_col_name][1] + rs[pmde_col_name][0]
        pm_c = np.sqrt(pmra_c**2 + pmde_c**2)

        pm_diff, pm_diff_err = calc_pm_diff(rs, pmra_col_name=pmra_col_name, pmde_col_name=pmde_col_name)
        pm_diff_sig = pm_diff / pm_diff_err

        cpm_index = pm_c / pm_diff

        return pm_diff, pm_diff_err, pm_diff_sig, cpm_index

    plx_diff, plx_diff_err, plx_diff_sig = calc_plx_diff(
        rs, parllax_col_name=parllax_col_name, parllax_err_col_name=parllax_err_col_name
    )
    pm_diff, pm_diff_err, pm_diff_sig, cpm_index = calc_pm_diff_n_cpm_index(
        rs, pmra_col_name=pmra_col_name, pmde_col_name=pmde_col_name
    )

    bound = plx_diff_sig < plx_diff_sig_threshold and cpm_index >= cpm_index_threshold

    if verbose:
        print(
            f"are_stars_bound_by_plx_pm(): bound={bound}, "
            f"plx_diff={plx_diff:.3f}, plx_diff_err={plx_diff_err:.3f}, plx_diff_sig={plx_diff_sig:.2f}, "
            f"pm_diff={pm_diff:.3f}, pm_diff_err={pm_diff_err:.3f}, pm_diff_sig={pm_diff_sig:.2f}, cpm_index={cpm_index:.1f}"
        )

    return bound


def to_separation_in_au(
    rs,  # result set table
):
    plx = rs["Plx"][0]  # assumed to be in mas
    sep_angular = rs["_r"][1]  # assumed to be in arcsec
    return 1000 / plx * sep_angular
