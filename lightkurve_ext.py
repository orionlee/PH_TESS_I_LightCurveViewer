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

import asyncio_compat

log = logging.getLogger(__name__)

def of_sector(lcf_coll, sectorNum):
    for lcf in lcf_coll:
        if lcf.meta['SECTOR'] == sectorNum:
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


def of_2min_cadences(lcf_coll):
    """Return LightCurveFiles of 2-minute cadence only.

       Primary use case is to filter out 20-second files.
    """
    # accept those in the range of 1.5 - 2.5 minutes
    filtered = [lcf for lcf in lcf_coll if 150 / 86400 >= np.median(np.diff(lcf.time[:100])) >= 90 / 86400]
    return lk.LightCurveCollection(filtered)


def of_tic(lcf_coll, tic):
    """Return LightCurveFiles of the given TIC.

    Useful in case the default MAST result returned nearby targets.
    """
    filtered = [lcf for lcf in lcf_coll if lcf.meta.get('TICID', None) == tic]
    return lk.LightCurveCollection(filtered)


def download_lightcurve(target, mission=('Kepler', 'K2', 'TESS'),
                             exptime='short',
                             author='SPOC',
                             download_dir=None, use_cache='yes',
                             display_search_result=True):
    '''
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
    '''

    if use_cache == 'no':
        return _search_and_cache(target, mission, exptime, author, download_dir, display_search_result)
    if use_cache == 'yes':
        result_file_ids = _load_from_cache_if_any(target, mission, download_dir)
        if result_file_ids is not None:
            result_files = list(map(lambda e: f"{download_dir}/mastDownload/{e}",
                                    result_file_ids))
            return lk.collections.LightCurveCollection(
                list(map(lambda f: lk.read(f), result_files)))
        # else
        return _search_and_cache(target, mission, exptime, author, download_dir, display_search_result)
    # else
    raise ValueError('invalid value for argument use_cache')

# Private helpers for `download_lightcurvefiles`


def _search_and_cache(target, mission, exptime, author, download_dir, display_search_result):
    search_res = lk.search_lightcurve(target=target, mission=mission, exptime=exptime, author=author)
    if len(search_res) < 1:
        return None
    if display_search_result:
        _display_search_result(search_res)
    _cache_search_result_product_identifiers(search_res, download_dir, target, mission)
    return search_res.download_all(quality_bitmask='default', download_dir=download_dir)


def _display_search_result(search_res):
    from IPython.core.display import display
    tab = search_res.table
    # move useful columns to the front
    preferred_cols = ['proposal_id', 'target_name', 'sequence_number', 't_exptime']
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
    # - still somewhat blocking in jupyter: the cell completes, but running in subsequent cells are stalled)
    # - try to use specific event loop to see if it helps, so far not promising.
    #    https://stackoverflow.com/questions/47518874/how-do-i-run-python-asyncio-code-in-a-jupyter-notebook
    return asyncio_compat.create_task(asyncio_compat.to_thread(search_and_download_tpf, *args, **kwargs))


#
# Other misc. extensions
#
def get_bkg_lightcurve(lcf):
    """Returns the background flux, i.e., ``SAP_BKG`` in the file"""
    lc = lcf.copy()
    lc['flux'] = lc['sap_bkg']
    lc['flux_err'] = lc['sap_bkg_err']
    lc.label = lc.label + ' BKG'
    return lc


def create_quality_issues_mask(lc, flags_included=0b0101001010111111):
    """Returns a boolean array which flags cadences with *issues*.

    The default `flags_included` is a TESS default, based on https://outerspace.stsci.edu/display/TESS/2.0+-+Data+Product+Overview#id-2.0DataProductOverview-Table:CadenceQualityFlags
    """

    return np.logical_and(lc.quality & flags_included, np.isfinite(lc.flux))

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
        return ','.join(map(str, times))
    else:
        return times

def get_segment_times_idx(times, break_tolerance=5):
    """Segment the input array of times into segments due to data gaps. Return the indices of the segments.

    The minimal gap size is determined by `break_tolerance`.

    The logic is adapted from `LightCurve.flatten`
    """
    if hasattr(times, 'value'):  # convert astropy Time to raw values if needed
        times = times.value
    dt = times[1:] - times[0:-1]
    with warnings.catch_warnings():  # Ignore warnings due to NaNs
        warnings.simplefilter("ignore", RuntimeWarning)
        cut = np.where(dt > break_tolerance * np.nanmedian(dt))[0] + 1
    low = np.append([0], cut)
    high = np.append(cut, len(times))
    return (low, high)

def get_segment_times(times, **kwargs):
    if hasattr(times, 'value'):  # convert astropy Time to raw values if needed
        times = times.value
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
    times_list = get_segment_times(lc.time)
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
    spline = UnivariateSpline(lc.time, lc.flux, **kwargs)
    trend = spline(lc.time)
    detrended = lc.flux - trend
    detrended_relative = 100 * ((lc.flux / trend) - 1) + 100 # in percentage
    lc_flattened = lc.copy()
    lc_flattened.flux = detrended_relative
    lc_flattened.flux_unit = 'percent'
    if not return_trend:
        return lc_flattened
    else:
        lc_trend = lc.copy()
        lc_trend.flux = trend
        return (lc_flattened, lc_trend)


# Backport tpf.plot_pixels() from lightkurve2 for use in lightkurve1
PATCH_LK = False
if PATCH_LK:
    # adapted from https://github.com/KeplerGO/lightkurve/blob/v2.0b3/lightkurve/targetpixelfile.py
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    from lightkurve import MPLSTYLE
    from lightkurve import LightkurveWarning
    def _plot_pixels(self, ax=None, periodogram=False, aperture_mask=None,
                    show_flux=False, corrector_func=None, style='lightkurve',
                    title=None, **kwargs):
        """Show the light curve of each pixel in a single plot.

        Note that all values are autoscaled and axis labels are not provided.
        This utility is designed for by-eye inspection of signal morphology.

        Parameters
        ----------
        ax : `~matplotlib.axes.Axes`
            A matplotlib axes object to plot into. If no axes is provided,
            a new one will be generated.
        periodogram : bool
            Default: False; if True, periodograms will be plotted, using normalized light curves.
            Note that this keyword overrides normalized.
        aperture_mask : ndarray or str
            Highlight pixels selected by aperture_mask.
            Only `pipeline`, `threshold`, or custom masks will be plotted.
            `all` and None masks will be ignored.
        show_flux : bool
            Default: False; if True, shade pixels with frame 0 flux colour
            Inspired by https://github.com/noraeisner/LATTE
        corrector_func : function
            Function that accepts and returns a `~lightkurve.lightcurve.LightCurve`.
            This function is applied to each light curve in the collection
            prior to stitching. The default is to normalize each light curve.
        style : str
            Path or URL to a matplotlib style file, or name of one of
            matplotlib's built-in stylesheets (e.g. 'ggplot').
            Lightkurve's custom stylesheet is used by default.
        kwargs : dict
            e.g. extra parameters to be passed to `lc.to_periodogram`.
        """
        if style == 'lightkurve' or style is None:
            style = MPLSTYLE
        if title is None:
            title = f'Target ID: {self.targetid}'
        if corrector_func is None:
            corrector_func = lambda x: x.remove_outliers()
        if show_flux:
            cmap = plt.get_cmap()
            norm = plt.Normalize(vmin=np.nanmin(self.flux[0]),
                                 vmax=np.nanmax(self.flux[0]))
        mask = self._parse_aperture_mask(aperture_mask)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=(RuntimeWarning, LightkurveWarning))

            # get an aperture mask for each pixel
            masks = np.zeros((self.shape[1]*self.shape[2], self.shape[1], self.shape[2]),
                             dtype='bool')
            for i in range(self.shape[1]*self.shape[2]):
                masks[i][np.unravel_index(i, (self.shape[1], self.shape[2]))] = True

            pixel_list = []
            for j in range(self.shape[1]*self.shape[2]):
                lc = self.to_lightcurve(aperture_mask=masks[j])
                lc = corrector_func(lc)

                if periodogram:
                    try:
                        pixel_list.append(lc.to_periodogram(**kwargs))
                    except IndexError:
                        pixel_list.append(None)
                else:
                    if len(lc.remove_nans().flux) == 0:
                        pixel_list.append(None)
                    else:
                        pixel_list.append(lc)

        with plt.style.context(style):
            fig = plt.figure()
            if ax is None:  # Configure axes if none is given
                ax = plt.gca()
                ax.get_xaxis().set_ticks([])
                ax.get_yaxis().set_ticks([])
                if periodogram:
                    ax.set(title=title, xlabel='Frequency', ylabel='Power')
                else:
                    ax.set(title=title, xlabel='Time', ylabel='Flux')

            gs = gridspec.GridSpec(self.shape[1], self.shape[2], wspace=0.01, hspace=0.01)

            for k in range(self.shape[1]*self.shape[2]):
                if pixel_list[k]:
                    x, y = np.unravel_index(k, (self.shape[1], self.shape[2]))

                    # Highlight aperture mask in red
                    if aperture_mask is not None and mask[x,y]:
                        rc = {"axes.linewidth": 2, "axes.edgecolor": 'red'}
                    else:
                        rc = {"axes.linewidth": 1}
                    with plt.rc_context(rc=rc):
                        gax = fig.add_subplot(gs[self.shape[1] - x - 1, y])

                    # Determine background and foreground color
                    if show_flux:
                        gax.set_facecolor(cmap(norm(self.flux[0,x,y])))
                        markercolor = "white"
                    else:
                        markercolor = "black"

                    # Plot flux or periodogram
                    if periodogram:
                        gax.plot(pixel_list[k].frequency.value,
                                 pixel_list[k].power.value,
                                 marker='None', color=markercolor, lw=0.5)
                    else:
                        gax.plot(pixel_list[k].time,
                                 pixel_list[k].flux,
                                 marker='.', color=markercolor, ms=0.5, lw=0)

                    gax.margins(y=.1, x=0)
                    gax.set_xticklabels('')
                    gax.set_yticklabels('')
                    gax.set_xticks([])
                    gax.set_yticks([])

            fig.set_size_inches((y*1.5, x*1.5))

        return ax


    lk.KeplerTargetPixelFile.plot_pixels = _plot_pixels
    lk.TessTargetPixelFile.plot_pixels = _plot_pixels


def plot_pixels(tpf, pixel_size_inch=1, **kwargs):
    """A thin wrapper over ``plot_pixels()`` that provides preferred UI tweak: figure size, etc."""
    ax = tpf.plot_pixels(**kwargs)

    y, x = tpf.shape[1], tpf.shape[2]
    ax.get_figure().set_size_inches((x * pixel_size_inch, y * pixel_size_inch))
    return ax
