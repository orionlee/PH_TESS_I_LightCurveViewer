"""
Convenience helpers for `lightkurve` package.
"""

import os
import logging
import json

import lightkurve as lk

log = logging.getLogger(__name__)

def of_sector(lcf_coll, sectorNum):
    for lcf in lcf_coll:
        if lcf.get_header()['SECTOR'] == sectorNum:
            return lcf
    return None

def of_sectors(*args):
    lcf_coll = args[0]
    sectorNums = args[1:]
    res = []
    for lcf in lcf_coll:
        if lcf.get_header()['SECTOR'] in sectorNums:
            res.append(lcf)
    return lk.LightCurveFileCollection(res)


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
        Returns a `~lightkurve.collections.LightCurveFileCollection`
        containing all lightcurve files that match the criteria
    '''

    if use_cache == 'no':
        return _search_and_cache(target, mission, download_dir)
    if use_cache == 'yes':
        result_file_ids = _load_from_cache_if_any(target, mission, download_dir)
        if result_file_ids is not None:
            result_files = list(map(lambda e: f"{download_dir}/mastDownload/{e}",
                                    result_file_ids))
            return lk.collections.LightCurveFileCollection(
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
