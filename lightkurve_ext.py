"""
Convenience helpers for `lightkurve` package.
"""

import os
import logging
import math
import json
import warnings

import numpy as np
from scipy.interpolate import UnivariateSpline

import lightkurve as lk

log = logging.getLogger(__name__)

def of_sector(lcf_coll, sectorNum):
    for lcf in lcf_coll:
        if lcf.meta['SECTOR'.lower()] == sectorNum:
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
    sectorNums = args[1:]
    res = []
    for lcf in lcf_coll:
        if lcf.meta['SECTOR'.lower()] in sectorNums:
            res.append(lcf)
    return lk.LightCurveCollection(res)


def download_lightcurvefiles(target, mission=('Kepler', 'K2', 'TESS'),
                             download_dir=None, use_cache='yes'):
    '''
    Wraps `lightkurve.search.search_lightcurvefile()` and the
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
    '''

    if use_cache == 'no':
        return _search_and_cache(target, mission, download_dir)
    if use_cache == 'yes':
        result_file_ids = _load_from_cache_if_any(target, mission, download_dir)
        if result_file_ids is not None:
            result_files = list(map(lambda e: f"{download_dir}/mastDownload/{e}",
                                    result_file_ids))
            return lk.collections.LightCurveCollection(
                list(map(lambda f: lk.open(f), result_files)))
        # else
        return _search_and_cache(target, mission, download_dir)
    # else
    raise ValueError('invalid value for argument use_cache')

# Private helpers for `download_lightcurvefiles`

def _search_and_cache(target, mission, download_dir):
    search_res = lk.search.search_lightcurvefile(target=target, mission=mission)
    _cache_search_result_product_identifiers(search_res, download_dir, target, mission)
    return search_res.download_all(quality_bitmask='default', download_dir=download_dir)

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
        log.warning('Warning: unable to create {}. '
                    'Cache MAST query results to the current '
                    'working directory instead.'.format(cache_dir))
        cache_dir = '.'
    return cache_dir

def _get_cache_key(target, mission):
    #TODO: handle cases the generated key is not a valid filename
    return f"{target}_{mission}_ids"


def _to_product_identifiers(search_res):
    '''
    Returns
    -------
    A list of str, constructed from `(obs_collection, obs_id, productFilename)` tuples, that can
    identify cached lightcurve file,s if any.
    '''
    return list(map(
        lambda e: e['obs_collection'] +  '/' + e['obs_id'] + '/' + e['productFilename'],
        search_res.table))

def _save_search_result_product_identifiers(identifiers, download_dir, key):
    resolved_cache_dir = _get_search_result_cache_dir(download_dir)
    filepath = f"{resolved_cache_dir}/{key}.json"
    fp = open(filepath, 'w+')
    json.dump(identifiers, fp)
    return filepath

def _load_search_result_product_identifiers(download_dir, key):
    resolved_cache_dir = _get_search_result_cache_dir(download_dir)
    filepath = f"{resolved_cache_dir}/{key}.json"
    try:
        fp = open(filepath, 'r')
        return json.load(fp)
    except OSError as err:
        # errno == 2: file not found, typical case of cache miss
        # errno != 2: unexpected error, log a warning
        if err.errno != 2:
            log.warning('Unexpected OSError in retrieving cached search result: {}'.format(err))
        return None


# Other misc. extensions

def create_quality_issues_mask(lc, flags_included=0b0101001010111111):
    """Returns a boolean array which flags cadences with *issues*.

    The default `flags_included` is a TESS default, based on https://outerspace.stsci.edu/display/TESS/2.0+-+Data+Product+Overview#id-2.0DataProductOverview-Table:CadenceQualityFlags
    """

    return np.nonzero(np.logical_and(lc.quality.value & flags_included, np.isfinite(lc.flux)))

def list_times_w_quality_issues(lc):
    mask = create_quality_issues_mask(lc)
    return lc.time.value[mask], lc.quality[mask]


def list_transit_times(t0, period, steps_or_num_transits=range(0, 10), return_string=False):
    """List the transit times based on the supplied transit parameters"""
    if isinstance(steps_or_num_transits, int):
        steps = range(0, steps_or_num_transits)
    else:
        steps = steps_or_num_transits
    times = [t0 + period * i for i in steps]
    if return_string:
        return ','.join(map(str, times))
    else:
        return times

def get_segment_times_idx(times, break_tolerance=5):
    """Segment the input array of times into segments due to data gaps. Return the indices of the segments.

    The minimal gap size is determined by `break_tolerance`.

    The logic is adapted from `LightCurve.flatten`
    """
    dt = times[1:] - times[0:-1]
    with warnings.catch_warnings():  # Ignore warnings due to NaNs
        warnings.simplefilter("ignore", RuntimeWarning)
        cut = np.where(dt > break_tolerance * np.nanmedian(dt))[0] + 1
    low = np.append([0], cut)
    high = np.append(cut, len(times))
    return (low, high)

def get_segment_times(times, **kwargs):
    low, high = get_segment_times_idx(times, **kwargs)
    # add a small 1e-10 to end so that the end time is exclusive (follow convention in range)
    return [(times[l], times[h-1] + 1e-10) for l, h in zip(low, high)]


def get_transit_times_in_range(t0, period, start, end):
    t_start = t0 + math.ceil((start - t0) / period) * period
    num_t = math.ceil((end - t_start) / period)
    return [t_start + period * i for i in range(num_t)]


def get_transit_times_in_lc(lc, t0, period, return_string=False, **kwargs):
    """Get the transit times with observations of the given lightcurve, based on the supplied transit parameters.

       The method will exclude the times where there is no observation due to data gaps.
    """

    # break up the times to exclude times in gap
    times_list = get_segment_times(lc.time.value)
    transit_times = []
    for start, end in times_list:
        transit_times.extend(get_transit_times_in_range(t0, period, start, end))
    if return_string:
        return ','.join(map(str, transit_times))
    else:
        return transit_times



def to_window_length_for_2min_cadence(length_day):
    """Helper for LightCurve.flatten(). Return a `window_length` for the given number of days, assuming the data has 2-minute cadence."""
    res = math.floor(720 * length_day)
    if res % 2 == 0:
        res += 1 # savgol_filter window length must be odd number
    return res


# detrend using spline
# Based on:  https://github.com/barentsen/kepler-athenaeum-tutorial/blob/master/how-to-find-a-planet-tutorial.ipynb
def flatten_with_spline_normalized(lc, return_trend=False, **kwargs):
    lc = lc.remove_nans()
    spline = UnivariateSpline(lc.time.value, lc.flux.value, **kwargs)
    trend = spline(lc.time.value)
    detrended = lc.flux.value - trend
    detrended_relative = 100 * ((lc.flux.value / trend) - 1) + 100 # in percentage
    lc_flattened = lc.copy()
    lc_flattened.flux = detrended_relative
    lc_flattened.flux_unit = 'percent'
    if not return_trend:
        return lc_flattened
    else:
        lc_trend = lc.copy()
        lc_trend.flux = trend
        return (lc_flattened, lc_trend)

