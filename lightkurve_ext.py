"""
Convenience helpers for `lightkurve` package.
"""

import os
import logging
import math
import json
import warnings
from collections import OrderedDict

import astropy.units as u
import numpy as np
from scipy.interpolate import UnivariateSpline

from IPython.display import display, HTML

import lightkurve as lk
from lightkurve.search import SearchResult

import asyncio_compat

log = logging.getLogger(__name__)


def of_sector(lcf_coll, sectorNum):
    for lcf in lcf_coll:
        if lcf.meta["SECTOR"] == sectorNum:
            return lcf
    return None


def of_sectors(*args):
    lcf_coll = args[0]
    if len(args) == 1:
        # when no sectors are specified, return entire collection
        # For convenience: when a notebooks is modified such that
        # a user sometimes use a subset of sectors , and sometimes everything
        # the user can can still use of_sectors() wrapper regardless
        return lcf_coll
    sector_nums = args[1:]
    return lcf_coll[np.in1d(lcf_coll.sector, sector_nums)]


def of_sector_n_around(lk_coll_or_sr, sector_num, num_additions=8):
    def do_for_lk_coll():
        subset_slice = _get_slice_for_of_sector_n_around(
            lk_coll_or_sr,
            lambda coll: coll.sector,
            sector_num,
            num_additions=num_additions,
        )
        if subset_slice is not None:
            # workaround bug that lcf_coll[start:end] returns a list only
            return lk.LightCurveCollection(lk_coll_or_sr[subset_slice])
        else:
            return lk.LightCurveCollection([])

    def do_for_sr():
        subset_slice = _get_slice_for_of_sector_n_around(
            lk_coll_or_sr,
            lambda sr: sr.table["sequence_number"],
            sector_num,
            num_additions=num_additions,
        )
        if subset_slice is not None:
            return lk_coll_or_sr[subset_slice]
        else:
            return SearchResult()

    if hasattr(lk_coll_or_sr, "sector"):
        return do_for_lk_coll()
    elif hasattr(lk_coll_or_sr, "table") and lk_coll_or_sr.table["sequence_number"] is not None:
        return do_for_sr()
    else:
        raise TypeError(f"Unsupported type of collection: {type(lk_coll_or_sr)}")


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


def estimate_cadence(lc):
    """Estimate the cadence of a lightcurve by returning the median of a sample"""
    return np.nanmedian(np.diff(lc.time[:100].value))


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


def estimate_object_radius_in_r_jupiter(lc, depth):
    """Return a back of envelope estimate of a companion object's radius."""
    R_JUPITER_IN_R_SUN = 71492 / 695700

    r_star = lc.meta.get("RADIUS")  # assumed to be in R_sun
    if r_star is None or depth <= 0:
        return None  # cannot estimate
    r_obj = math.sqrt(r_star * r_star * depth)
    r_obj_in_r_jupiter = r_obj / R_JUPITER_IN_R_SUN
    return r_obj_in_r_jupiter


def download_lightcurves_of_tic_with_priority(tic, download_filter_func=None, download_dir=None):
    """For a given TIC, download lightcurves across all sectors.
    For each sector, download one based on pre-set priority.
    """

    sr_unfiltered = lk.search_lightcurve(f"TIC{tic}", mission="TESS")
    if len(sr_unfiltered) < 1:
        print(f"WARNING: no result found for TIC {tic}")
        return None, None, None

    sr_unfiltered = sr_unfiltered[sr_unfiltered.target_name == str(tic)]  # in case we get some other nearby TICs

    # filter out HLSPs not supported by lightkurve yet
    sr = sr_unfiltered[sr_unfiltered.author != "DIAMANTE"]
    if len(sr) < len(sr_unfiltered):
        print("Note: there are products not supported by Lightkurve, which are excluded from download.")

    # for each sector, filter based on the given priority.
    # - note: prefer QLP over TESS-SPOC because QLP is detrended, with multiple apertures within 1 file
    sr = filter_by_priority(
        sr,
        author_priority=["SPOC", "QLP", "TESS-SPOC"],
        exptime_priority=["short", "long", "fast"],
    )
    num_filtered = len(sr_unfiltered) - len(sr)
    num_fast = len(sr_unfiltered[sr_unfiltered.exptime < 60 * u.second])
    if num_filtered > 0:
        msg = f"{num_filtered} rows filtered"
        if num_fast > 0:
            msg = msg + f" ; {num_fast} fast (20secs) products."
        print(msg)

    display(sr)

    # let caller to optionally further restrict a subset to be downloaded
    sr_to_download = sr
    if download_filter_func is not None:
        sr_to_download = download_filter_func(sr)
        if len(sr_to_download) < len(sr):
            display(
                HTML(
                    f"""<font style="background-color: yellow;">Note</font>:
SearchResult is further filtered - only a subset will be downloaded."""
                )
            )

    lcf_coll = sr_to_download.download_all(download_dir=download_dir)

    if lcf_coll is not None and len(lcf_coll) > 0:
        print(f"TIC {tic} \t#sectors: {len(lcf_coll)} ; {lcf_coll[0].meta['SECTOR']} - {lcf_coll[-1].meta['SECTOR']}")
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

    return lk.SearchResult(table=res_t)


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
    sr = lk.search_targetpixelfile(*args, **kwargs)  # pass the rest of the argument to search_targetpixelfile
    tpf_coll = sr.download_all(download_dir=download_dir, quality_bitmask=quality_bitmask)
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


def create_quality_issues_mask(lc, flags_included=0b0101001010111111):
    """Returns a boolean array which flags cadences with *issues*.

    The default `flags_included` is a TESS default, based on
    https://outerspace.stsci.edu/display/TESS/2.0+-+Data+Product+Overview#id-2.0DataProductOverview-Table:CadenceQualityFlags
    """
    if np.issubdtype(lc["quality"].dtype, np.integer):
        return np.logical_and(lc.quality & flags_included, np.isfinite(lc.flux))
    else:
        # quality column is not an integer, probably a non-standard product
        return np.zeros_like(lc.flux, dtype=bool)


def list_times_w_quality_issues(lc):
    mask = create_quality_issues_mask(lc)
    return lc.time[mask], lc.quality[mask]


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
    dt = times[1:] - times[0:-1]
    with warnings.catch_warnings():  # Ignore warnings due to NaNs
        warnings.simplefilter("ignore", RuntimeWarning)
        cut = np.where(dt > break_tolerance * np.nanmedian(dt))[0] + 1
    low = np.append([0], cut)
    high = np.append(cut, len(times))
    return (low, high)


def get_segment_times(times, **kwargs):
    if hasattr(times, "value"):  # convert astropy Time to raw values if needed
        times = times.value
    low, high = get_segment_times_idx(times, **kwargs)
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

    # break up the times to exclude times in gap
    times_list = get_segment_times(lc.time)
    transit_times = []
    for start, end in times_list:
        transit_times.extend(get_transit_times_in_range(t0, period, start, end))
    if return_string:
        return ",".join(map(str, transit_times))
    else:
        return transit_times


def to_window_length_for_2min_cadence(length_day):
    """Helper for LightCurve.flatten().
    Return a `window_length` for the given number of days, assuming the data has 2-minute cadence."""
    res = math.floor(720 * length_day)
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
