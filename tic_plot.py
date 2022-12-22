# -*- coding: utf-8 -*-
"""
Helpers to plot the lightcurve of a TESS subject, given a
LightCurveCollection
"""
# so that type hint Optional(LC_Ylim_Func_Type) can be used
# otherwise Python complains TypeError: Cannot instantiate typing.Optional
from __future__ import annotations

import inspect
import warnings
import re
from types import SimpleNamespace


import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, AutoMinorLocator
import matplotlib.animation as animation
import numpy as np
import pandas as pd

from astropy.io import fits
from astropy import units as u
from astropy.table import Table

import IPython
from IPython.display import display, HTML, Audio
from ipywidgets import interactive, interactive_output, fixed
import ipywidgets as widgets

import lightkurve as lk
from lightkurve import LightCurve, LightCurveCollection, LightkurveWarning, FoldedLightCurve
from lightkurve.utils import TessQualityFlags
from lightkurve_ext import of_sectors
import lightkurve_ext as lke
import lightkurve_ext_tess as lket

# typing
from typing import Callable, Optional, Tuple

LC_Ylim_Func_Type = Callable[[LightCurve], Tuple[float, float]]


def get_tic_meta_in_html(lc, a_subject_id=None, download_dir=None):
    # import locally so that if it fails (due to missing dependency)
    # it will not break the rest of the codes
    import lightkurve_ext_tess as lke_tess

    return lke_tess.get_tic_meta_in_html(lc, a_subject_id=a_subject_id, download_dir=download_dir)


def beep():
    """Emits a beep sound. It works only in IPython / Jupyter environment only"""
    # a beep to remind the users that the data has been downloaded
    # css tweak to hide beep
    display(
        HTML(
            """<script>
function tweakCSS() {
  if (document.getElementById("hide-beep-css")) {
      return;
  }
  document.head.insertAdjacentHTML('beforeend', `<style id="hide-beep-css" type="text/css">
  #beep { /* hide the audio control for the beep, generated from tplt.beep() */
    width: 1px;
    height: 1px;
  }
</style>`);
}
tweakCSS();
</script>
"""
        )
    )
    # the actual beep
    ## somehow ssl error
    ## beep_url = "https://upload.wikimedia.org/wikipedia/commons/f/fb/NEC_PC-9801VX_ITF_beep_sound.ogg"
    beep_url = "beep_sound.ogg"
    if int(re.sub(r"[.].+", "", IPython.__version__)) < 7:
        # compatibility with older older IPython (e.g., google colab)
        ## audio = Audio(url=beep_url, autoplay=True, embed=True)
        audio = Audio(filename=beep_url, autoplay=True, embed=True)
    else:
        ## audio = Audio(url=beep_url, autoplay=True, embed=True, element_id="beep")
        audio = Audio(filename=beep_url, autoplay=True, embed=True, element_id="beep")
    display(audio)


def _normalize_to_percent_quiet(lc):
    # Some product are in normalized flux, e.g., as 1, we still want to normalize them to percentage
    # for consistency
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=LightkurveWarning, message=".*in relative units.*")
        return lc.normalize(unit="percent")


# Plot the flux changes (not flux themselves) to get a sense of the rate of changes, not too helpful yet.
def plot_lcf_flux_delta(lcf, ax, xmin=None, xmax=None, moving_avg_window="30min"):

    # possible input arguments
    lc = _normalize_to_percent_quiet(lcf)

    # Basic scatter of the observation
    #    ax = lc.scatter(ax=ax)

    # convert to dataframe to add moving average
    df = lc.to_pandas()
    df["time_ts"] = [pd.Timestamp(x, unit="D") for x in df.index]
    # the timestamp above is good for relative time.
    # if we want the timestamp to reflect the actual time, we need to convert the BTJD in time to timetamp, e.g.
    # pd.Timestamp(astropy.time.Time(x + 2457000, format='jd', scale='tdb').datetime.timestamp(), unit='s')
    df["flux_mavg"] = df.rolling(moving_avg_window, on="time_ts")["flux"].mean()
    #    ax.plot(lc.time.value, df['flux_mavg'], c='black', label=f"Moving average ({moving_avg_window})")

    df["flux_delta"] = df.rolling(moving_avg_window, on="time_ts")["flux_mavg"].apply(
        lambda vary: vary[-1] - vary[0], raw=True
    )
    ax.plot(
        lc.time.value,
        df["flux_delta"],
        c="blue",
        label=f"Flux delta ({moving_avg_window})",
    )

    ax.set_xlim(xmin, xmax)

    return ax


def lk_ax(*args, **kwargs):
    """Create a matplotlib figure, and return its Axes object (`gca()`) with Lightkurve style."""
    with plt.style.context(lk.MPLSTYLE):
        # MUST return the Axes object, rather than just the Figure object
        # if only the Figure object is returned, it will have Lightkurve's default figsize
        # but its style won't be used in actual plot
        return plt.figure(*args, **kwargs).gca()


def flux_near(lc, time):
    if time is None or lc is None:
        return None
    else:
        idx = (np.abs(lc.time - time)).argmin()
        return lc.flux[idx]


def flux_mavg_near(df, time):
    if time is None or df is None:
        return None
    else:
        idx = (np.abs(df.index.values - time)).argmin()
        # must use df.iloc[idx]['flux_mavg'], rather than df['flux_mavg'][idx]
        # because dataframe from lightkurve is indexed by time (rather than regular 0-based index)
        # df.iloc[] ensures we can still access the value by 0-based index
        return df.iloc[idx]["flux_mavg"]


def _to_unitless(n):
    if hasattr(n, "value"):
        return n.value
    else:
        return n


def as_4decimal(float_num):
    if float_num is None:
        return None
    elif isinstance(float_num, tuple) or isinstance(float_num, list):
        return [float("{0:.4f}".format(_to_unitless(n))) for n in float_num]
    else:
        return float("{0:.4f}".format(_to_unitless(float_num)))


def scatter(lc, **kwargs):
    """lc.scatter() with the proper support of plotting flux in magnitudes"""
    ax = lc.scatter(**kwargs)
    y_column = "flux"
    if lc[y_column].unit is u.mag:
        ax.invert_yaxis()
    return ax


def add_flux_moving_average(lc, moving_avg_window):
    df = lc.to_pandas()
    begin_t = df.index[0]
    df["time_ts"] = [pd.Timestamp(t - begin_t, unit="D") for t in df.index]
    # the timestamp above is good for relative time.
    # 1. we subtract the time with the timestamp because for some products, e.g., CDIPS, the time value itself
    #    is so large that creating pd.Timestamp with it causes Overflow error
    # 2. if we want the timestamp to reflect the actual time, we need to convert the BTJD in time to timetamp, e.g.
    #      pd.Timestamp(astropy.time.Time(x + 2457000, format='jd', scale='tdb').datetime.timestamp(), unit='s')
    df["flux_mavg"] = df.rolling(moving_avg_window, on="time_ts")["flux"].mean()
    return df


def add_relative_time(lc, lcf):
    t_start = lcf.meta.get("TSTART")
    if t_start is None:
        return False
    lc["time_rel"] = lc.time - t_start
    return True


def mask_gap(x, y, min_x_diff):
    """
    Help to plot graphs with gaps in the data, so that straight line won't be draw to fill the gap.
    Return a masked y that can be passed to pyplot.plot() that can show the gap.
    """
    # handle case that x is a astropy Time object, rather than simple float array
    x = _to_unitless(x)

    x_diff = np.diff(x, prepend=-min_x_diff)
    return np.ma.masked_where(x_diff > min_x_diff, y)


def normalize_percent(lc):
    """
    A syntactic surgar for lambda for normalize as percentage.
    Useful when calling ``lc.fold()``, ``tpf.interact()``, etc.
    """
    return lc.normalize(unit="percent")


def _add_flux_origin_to_ylabel(ax, lc):
    # e.g., QLP uses sap_flux as the standard.
    # it needs to support other products too
    standard_flux_col_map = {
        "SPOC": "pdcsap_flux",
        "TESS-SPOC": "pdcsap_flux",
        "QLP": "sap_flux",
        "TASOC": "flux_raw",
        "CDIPS": "irm1",
        "PATHOS": "psf_flux_cor",
    }

    def make_italic(text):
        # convert the text to latex italic expression
        return r"$\it{" + text.replace("_", r"\_") + "}$"

    flux_origin = lc.meta.get("FLUX_ORIGIN", None)
    if flux_origin is not None and flux_origin != standard_flux_col_map.get(lc.meta.get("AUTHOR", None), None):
        ax.yaxis.set_label_text(ax.yaxis.get_label_text().replace("Flux", make_italic(lc.flux_origin)))


_cache_plot_n_annotate_lcf = dict(lcf=None, flux_col=None, normalize=None, lc=None)


def plot_n_annotate_lcf(
    lcf,
    ax,
    flux_col="flux",
    xmin=None,
    xmax=None,
    t0=None,
    t_start=None,
    t_end=None,
    moving_avg_window="30min",
    t0mark_ymax=0.3,
    mark_momentum_dumps=True,
    set_title=True,
    show_r_obj_estimate=True,
    title_fontsize=18,
    t0_label_suffix=None,
    normalize=True,
    lc_tweak_fn=None,
    ax_tweak_fn=None,
    legend_kwargs=dict(),
):
    if lcf is None:
        print("Warning: lcf is None. Plot skipped")
        return

    # cache lc to speed up plots repeatedly over the same lcf
    global _cache_plot_n_annotate_lcf
    if (
        lcf is _cache_plot_n_annotate_lcf["lcf"]
        and flux_col == _cache_plot_n_annotate_lcf["flux_col"]
        and normalize == _cache_plot_n_annotate_lcf["normalize"]
    ):
        lc = _cache_plot_n_annotate_lcf["lc"]
    else:
        lc = lke.select_flux(lcf, flux_col)
        if normalize:
            lc = _normalize_to_percent_quiet(lc)
        _cache_plot_n_annotate_lcf["lcf"] = lcf
        _cache_plot_n_annotate_lcf["flux_col"] = flux_col
        _cache_plot_n_annotate_lcf["normalize"] = normalize
        _cache_plot_n_annotate_lcf["lc"] = lc

    if xmin is None and t_start is not None:
        xmin = t_start - 0.5
    if xmax is None and t_end is not None:
        xmax = t_end + 0.5

    # implement xmin / xmax by limiting the LC itself, rather than using ax.set_xlim after the plot
    # - the Y-scale will then automatically scaled to the specified time range, rather than over entire lightcurve
    # - make plotting faster (fewer data points)
    if xmin is not None:
        lc = lc[lc.time.value >= xmin]
    if xmax is not None:
        lc = lc[lc.time.value <= xmax]

    if lc_tweak_fn is not None:
        lc = lc_tweak_fn(lc)

    lcfh = lcf.meta

    # Basic scatter of the observation
    if "long" == lke.estimate_cadence_type(lc):
        # long cadence has more spare data, use a larger "x" to represent them
        # "x" is also useful to distinguish it from moving average,
        # which will likely overlap with the points given the sparse data
        ax = lc.scatter(ax=ax, s=36, marker="x")
    else:
        ax = lc.scatter(ax=ax)

    if len(lc) < 1:
        print(
            (
                "Warning: specified (xmin, xmax) is out of the range of the lightcurve "
                f"{lc.label} sector {lcfh['SECTOR']}. Nothing to plot"
            )
        )
        return ax

    # convert to dataframe to add moving average
    if moving_avg_window is not None:
        df = add_flux_moving_average(lc, moving_avg_window)
        # mask_gap: if there is a gap larger than 2 hours,
        # show the gap rather than trying to fill the gap with a straight line.
        ax.plot(
            lc.time.value,
            mask_gap(lc.time, df["flux_mavg"], 2 / 24),
            c="#3AF",
            label=f"Moving average ({moving_avg_window})",
        )
    else:
        df = add_flux_moving_average(lc, "10min")  # still needed for some subsequent calc, but don't plot it

    # annotate the graph
    if t_start is not None:
        ax.axvline(t_start)
    if t_end is not None:
        ax.axvline(t_end)
    if t0 is not None:
        t_lc_start = lcf.meta.get("TSTART", None)
        t0_rel_text = ""
        if t_lc_start is not None:
            t0_rel = t0 - t_lc_start
            t0_rel_text = f" ({t0_rel:.3f})"
        label_vline = f"t0 ~= {t0:.3f}{t0_rel_text}"
        if t0_label_suffix is not None:
            label_vline = f"{label_vline}\n{t0_label_suffix}"
        ax.axvline(
            t0,
            ymin=0,
            ymax=t0mark_ymax,
            color="black",
            linewidth=3,
            linestyle="--",
            label=label_vline,
        )

    if mark_momentum_dumps:
        plot_momentum_dumps(lc, ax)

    if set_title:
        title_text = lc.label
        if len(lcfh.get("SECTORS", [])) > 1:
            sector_text = lke.abbrev_sector_list(lcfh.get("SECTORS", []))
        else:
            sector_text = lcfh.get("SECTOR", None)
        if sector_text is not None:
            title_text += f", sector {sector_text}"

        author = lc.meta.get("AUTHOR", None)
        if author is not None and author != "SPOC":
            title_text += f", by {author}"
        if t0 is not None:
            transit_duration_msg = ""
            if t_start is not None and t_end is not None:
                transit_duration_msg = f"\ntransit duration ~= {as_4decimal(24 * (t_end - t_start))}h"
            flux_t0 = flux_mavg_near(df, t0)
            if flux_t0 is not None:
                flux_begin = max(flux_mavg_near(df, t_start), flux_mavg_near(df, t_end))
                flux_dip = flux_begin - flux_t0
                r_obj_msg = ""
                r_obj = lke.estimate_object_radius_in_r_jupiter(lc, flux_dip / 100)  # convert flux_dip in percent to fractions
                if show_r_obj_estimate and r_obj is not None:
                    r_obj_msg = f", R_p ~= {r_obj:0.2f} R_j"
                title_text += (
                    f" \nflux@$t_0$ ~= {as_4decimal(flux_t0)}%, "
                    f"dip ~= {as_4decimal(flux_dip)}%{r_obj_msg}{transit_duration_msg}"
                )
        ax.set_title(title_text, {"fontsize": title_fontsize})
    ax.legend(**legend_kwargs)

    _add_flux_origin_to_ylabel(ax, lc)

    ax.xaxis.label.set_size(18)
    ax.yaxis.label.set_size(18)

    # to avoid occasional formating in scientific notations
    ax.xaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))

    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(axis="x", which="minor", length=4)
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(axis="y", which="minor", length=4)

    if ax_tweak_fn is not None:
        ax_tweak_fn(ax)

    return ax


def plot_transit(lcf, ax, t0, duration, surround_time, **kwargs):
    return plot_n_annotate_lcf(
        lcf,
        ax=ax,
        t0=t0 if duration > 0 else None,
        t_start=t0 - duration / 2 if duration > 0 else None,
        t_end=t0 + duration / 2 if duration > 0 else None,
        xmin=t0 - (duration + surround_time) / 2,
        xmax=t0 + (duration + surround_time) / 2,
        **kwargs,
    )


def plot_transits(lcf_coll, transit_specs, ax_fn=lambda: lk_ax(), **kwargs):
    """Helper to plot transits zoomed-in."""
    flux_col = kwargs.get("flux_col", "flux")
    if not isinstance(flux_col, str) or flux_col.lower() not in ["flux", "pdcsap_flux"]:
        display(HTML(f"""<span style="background-color: yellow"> Note: </span> Not standard flux is plotted: {flux_col}"""))

    axs = []
    for spec in transit_specs:
        for lcf in of_sectors(lcf_coll, spec["sector"]):  # in case we have multiple lcf per sector
            #  process the supplied spec and apply defaults
            t0 = spec["epoch"]
            duration = spec["duration_hr"] / 24
            period = spec["period"]
            steps_to_show = spec["steps_to_show"]
            surround_time = spec.get("surround_time", 1.5)  # a hardcoded last resort default

            # TODO: warn if period is 0, but steps to show is not [0]

            for i in steps_to_show:
                cur_t0 = t0 + period * i
                t0_label_suffix = None
                if spec.get("label", "") not in ["", "dip", "dips"]:
                    t0_label_suffix = spec.get("label")
                ax = plot_transit(lcf, ax_fn(), cur_t0, duration, surround_time, t0_label_suffix=t0_label_suffix, **kwargs)
                axs.append(ax)
    return axs


def print_data_range(lcf_coll):
    """Print the data range for the given LightCurveCollection

    For each LightCurveFile:
    * sector start/stop time
    * first / last observation time
    * camera used
    """
    html = '<pre style="line-height: 1.1;">\n'
    html += "<summary>Sectors: " + lke.abbrev_sector_list(lcf_coll) + f" ({len(lcf_coll)})" + "\n"
    html += "Observation period range / data range:" + "\n"
    html += "<details>"
    for lc in lcf_coll:
        html += f"  Sector {lc.meta.get('SECTOR')}: {lc.meta.get('TSTART')} - {lc.meta.get('TSTOP')}" + "\n"
        html += f"   (cam.ccd {lc.meta.get('CAMERA')}.{lc.meta.get('CCD')})   {lc.time.min()} - {lc.time.max()}" + "\n"
    html += "</details></summary></pre>"
    display(HTML(html))


def get_momentum_dump_times(lcf):

    # momentum_dump_times_from_file() is no longer useful
    # `lket.MomentumDumpsAccessor.get_in_range()` is more flexible (e.g., better support for stitched lc better), and faster
    # but I keep it here as reference
    def momentum_dump_times_from_file():
        # Note: momentum_dump signals are by default masked out in LightCurve objects.
        # To access times marked as such, I need to access the raw LightCurveFile directly.
        filename = lcf.meta.get("FILENAME", None)
        if filename is None:
            warnings.warn("get_momentum_dump_times(): No-Op, because there is the LightCurve object has no backing FITS file.")
            return np.array([])
        with fits.open(filename) as hdu:
            if "TIME" not in hdu[1].columns.names:
                # case the file has no TIME column, typically non SPOC-produced ones, e.g., CDIPS,
                # the logic of finding momentum dump would not apply to such files anyway.
                return np.array([])

            # normal flow
            time = hdu[1].data["TIME"]
            mom_dumps_mask = np.bitwise_and(hdu[1].data["QUALITY"], TessQualityFlags.Desat) >= 1
            time_mom_dumps = time[mom_dumps_mask]
            return time_mom_dumps

    # main logic
    time_mom_dumps = lcf.meta.get("momentum_dumps", None)
    if time_mom_dumps is None:
        time_mom_dumps = lket.MomentumDumpsAccessor.get_in_range(lcf)
        lcf.meta["momentum_dumps"] = time_mom_dumps

    # in case the lcf has been truncated, we preserve the truncation
    return time_mom_dumps[(lcf.time.min().value <= time_mom_dumps) & (time_mom_dumps <= lcf.time.max().value)]


def vlines_y_in_axes_coord(ax, x, ymin, ymax, **kwargs):
    """Wrapper over `Axes.vlines()` for cases where ymin/ymax are in axes coordinates.
    It is to workaround bug: https://github.com/matplotlib/matplotlib/issues/23171
    """
    ybottom, ytop = ax.get_ylim()  # saved for workaround

    trans = ax.get_xaxis_transform()  # for ymin/ymax in axes coordinates
    if kwargs.get("transform", None) is not None and kwargs["transform"] is not trans:
        raise ValueError("_vlines_y_in_axes_coord() does not accept transform() parameter. It uses its own")
    kwargs["transform"] = trans
    res = ax.vlines(x, ymin, ymax, **kwargs)

    # Applying workaround: the bug scaled the ax's y-axis incorrectly, we compensate it by
    # rescaling it using the saved ylim
    #
    # edge case: if ymax > 1, we need to rescale the ytop
    if ymax > 1:
        ytop = ytop * ymax
    ax.set_ylim(ybottom, ytop)

    return res


def plot_momentum_dumps(lcf, ax, use_relative_time=False, mark_height_scale=0.15, color="red"):
    """Mark  momentum dumps on the given plot."""
    time_mom_dumps = get_momentum_dump_times(lcf)
    if len(time_mom_dumps) < 1:
        return ax
    # case have data to plot
    if use_relative_time:
        t_start = lcf.meta.get("TSTART")
        time_mom_dumps = time_mom_dumps - t_start

    vlines_y_in_axes_coord(
        ax,
        time_mom_dumps,
        ymin=0,
        ymax=mark_height_scale,
        color=color,
        linewidth=1,
        linestyle="-.",
        label="Momentum dumps",
    )

    return ax


# Do the actual plots
def plot_all(
    lcf_coll,
    flux_col="flux",
    moving_avg_window=None,
    lc_tweak_fn=None,
    ax_fn=None,
    use_relative_time=False,
    mark_quality_issues=True,
    mark_momentum_dumps=True,
    set_title=True,
    ax_tweak_fn=None,
):
    """Plot the given LightCurveFile collection, one graph for each LightCurve

    Returns
    -------
    axs : the list of plots in `matplotlib.Axes`
    """
    # choice 1: use the built-in plot method
    #    ax_all = plt.figure(figsize=(30, 15)).gca()
    #     lcf_coll.PDCSAP_FLUX.plot(ax=ax_all) # Or lcf_coll.SAP_FLUX.plot()

    # choice 2: stitch lightcurves of the collection together, and then use more flexible methods, e.g., scatter
    #   Note: pass lambda x: x to stitch() so that the code won't normalize the flux value sector by sector
    #     lc_all = lcf_coll.PDCSAP_FLUX.stitch(lambda x: x)
    #     lc_all.scatter(ax=ax_all, normalize=True)

    # choice 3: plot the lightcurve sector by sector: each sector has its own color
    #     for i in range(0, len(lcf_coll)):
    #         lcf_coll[i].PDCSAP_FLUX.scatter(ax=ax_all)

    #     ax_all.set_title((f"TIC {lcf_coll[0].PDCSAP_FLUX.label}, "
    #                       f"sectors {list(map(lambda lcf: lcf.meta.get('SECTOR'), lcf_coll))}"))
    #     return ax_all

    # choice 4: plot the lightcurve sector by sector: each sector in its own graph
    axs = []
    for i in range(0, len(lcf_coll)):
        if ax_fn is None:
            ax = lk_ax()
        else:
            ax = ax_fn()
        lcf = lcf_coll[i]
        lc = lke.select_flux(lcf, flux_col)

        lc = _normalize_to_percent_quiet(lc)
        if lc_tweak_fn is not None:
            lc = lc_tweak_fn(lc)

        # temporarily change time to a relative one if specified
        if use_relative_time:
            rel_time_added = add_relative_time(lc, lcf)
            if rel_time_added:
                lc["time_orig"] = lc.time
                lc.time = lc.time_rel
            else:
                # the file has no observation start time, so we cannot add it
                use_relative_time = False

        # tweak label to include sector if any
        sector = lcf_coll[i].meta.get("SECTOR", None)
        label_long = lc.label
        if sector is not None:
            lc.label += f", s.{sector}"
            label_long += f", sector {sector}"
        if lc.author is not None and lc.author != "SPOC":
            label_long += f", by {lc.author}"

        if "long" == lke.estimate_cadence_type(lc):
            # long cadence has more spare data, use a larger "x" to represent them
            # "x" is also useful to distinguish it from moving average,
            # which will likely overlap with the points given the sparse data
            ax = lc.scatter(ax=ax, s=16, marker="x")
        else:
            # for typical short cadence data, make dots smaller than the default (s=4)
            # so that the output doesn't look overly dense
            ax = lc.scatter(ax=ax, s=0.5)

        # convert to dataframe to add moving average
        if moving_avg_window is not None:
            df = add_flux_moving_average(lc, moving_avg_window)
            # mask_gap: if there is a gap larger than 2 hours,
            # show the gap rather than trying to fill the gap with a straight line.
            ax.plot(
                lc.time.value,
                mask_gap(lc.time, df["flux_mavg"], 2 / 24),
                c="#3AF",
                lw=0.4,
                label=f"Moving average ({moving_avg_window})",
            )

        title_extras = ""
        if lc_tweak_fn is not None:
            title_extras = "\nLC tweaked, e.g., outliers removed"

        if set_title:
            ax.set_title(f"{label_long} {title_extras}")  # {"fontsize": 18}
        if use_relative_time:
            ax.xaxis.set_label_text("Time - relative")
            # restore original time after plot is done
            lc.time = lc.time_orig
        else:
            t_start = lc.meta.get("TSTART")
            if t_start is not None:
                ax.xaxis.set_label_text(ax.xaxis.label.get_text() + f", TSTART={t_start:0.2f}")

        _add_flux_origin_to_ylabel(ax, lc)

        # to avoid occasional formating in scientific notations
        ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))

        # minor tick, 1 day interval in practice
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.tick_params(axis="x", which="minor", length=4)
        # ax.xaxis.grid(True, which='minor') # too noisy to be there by default

        ax.xaxis.label.set_size(fontsize=18)
        ax.yaxis.label.set_size(fontsize=18)
        if ax_tweak_fn is not None:
            ax_tweak_fn(ax)

        # mark quality issue is applied after ax_tweak_fn, in case users use ax_tweak_fn and change the graph's ylim
        if mark_quality_issues:
            # the time where flux might have potential issues, using the suggested starting quality flag mask
            time = lc.time if not use_relative_time else lc.time_rel
            time_w_quality_issues = time[lke.create_quality_issues_mask(lc)]
            if len(time_w_quality_issues) > 0:
                # add marks as vertical lines at bottom 10% of the plot
                vlines_y_in_axes_coord(
                    ax,
                    time_w_quality_issues.value,
                    ymin=0,
                    ymax=0.1,
                    color="red",
                    linewidth=1,
                    linestyle="--",
                    label="potential quality issue",
                )

        if mark_momentum_dumps:
            plot_momentum_dumps(lcf, ax, use_relative_time=use_relative_time)

        ax.legend()
        axs.append(ax)
    return axs


_lcf_4_plot_interactive = None


def _update_plot_lcf_interactive(figsize, flux_col, xrange, moving_avg_window, ymin, ymax, widget_out2):
    # use global to accept lct
    global _lcf_4_plot_interactive
    lcf = _lcf_4_plot_interactive

    ax = lk_ax(figsize=figsize)
    plot_n_annotate_lcf(
        lcf,
        ax,
        flux_col=flux_col,
        xmin=xrange[0],
        xmax=xrange[1],
        moving_avg_window=moving_avg_window,
    )
    codes_text = f"ax.set_xlim({xrange[0]}, {xrange[1]})"
    ymin_to_use = ymin if ymin >= 0 else None
    ymax_to_use = ymax if ymax >= 0 else None
    if (ymin_to_use is not None) or (ymax_to_use is not None):
        ax.set_ylim(ymin_to_use, ymax_to_use)
        codes_text += f"\n\nax.set_ylim({ymin_to_use}, {ymax_to_use})"

    widget_out2.clear_output()
    with widget_out2:
        print(codes_text)

    return None


def plot_lcf_interactive(lcf, figsize=(15, 8), flux_col="flux"):
    desc_style = {"description_width": "25ch"}
    slider_style = {"description_width": "25ch"}
    slider_layout = {"width": "100ch"}
    t_start = lcf.meta.get("TSTART")
    t_stop = lcf.meta.get("TSTOP")
    # Add a second output for textual
    widget_out2 = widgets.Output()

    # pass lcf with a global rather than the slow fixed(lcf) with lkv2
    #
    # import warnings
    # with warnings.catch_warnings():
    #     # lkv2 workaround: to suppress astropy table warning, stating that the semantics of == will be changed in the future.
    #     warnings.filterwarnings("ignore", category=FutureWarning)
    #     fixed_lcf = fixed(lcf)
    global _lcf_4_plot_interactive
    _lcf_4_plot_interactive = lcf

    w = interactive(
        _update_plot_lcf_interactive,
        figsize=fixed(figsize),
        # lcf = fixed_lcf,
        flux_col=fixed(flux_col),
        xrange=widgets.FloatRangeSlider(
            min=t_start,
            max=t_stop,
            step=0.1,
            value=(t_start, t_stop),
            description="Time",
            continuous_update=False,
            readout_format=".1f",
            layout=slider_layout,
            style=slider_style,
        ),
        moving_avg_window=widgets.Dropdown(
            options=[
                ("None", None),
                ("10 min", "20min"),
                ("20 min", "20min"),
                ("30 min", "30min"),
                ("1 hour", "1h"),
                ("2 hours", "2h"),
                ("4 hours", "4h"),
            ],
            value="30min",
            description="Moving average window",
            style=desc_style,
        ),
        ymin=widgets.FloatText(value=-1, description="Flux min, -1 for default", style=desc_style),
        ymax=widgets.FloatText(value=-1, description="Flux max, -1 for default", style=desc_style),
        widget_out2=fixed(widget_out2),
    )
    w.layout.border = "1px solid lightgray"
    w.layout.padding = "1em 0px"

    widget_out2.layout.padding = "1em"
    w.children = w.children + (widget_out2,)

    display(w)
    return w


def plot_transit_interactive(lcf, figsize=(15, 8), flux_col="flux", defaults=None):
    def _update_plot_transit_interactive(
        flux_col,
        t0,
        duration_hr,
        period,
        step,
        surround_time,
        moving_avg_window,
        t0mark_ymax,
        ymin,
        ymax,
        widget_out2,
    ):
        # for typical inline matplotlib backend, the figure needs to be recreated every time.
        ax = lk_ax(figsize=figsize)
        codes_text = "# Snippets to generate the plot"
        moving_avg_window_for_codes = "None" if moving_avg_window is None else f"'{moving_avg_window}'"
        if t0 < 0:
            plot_n_annotate_lcf(lcf, ax, flux_col=flux_col, moving_avg_window=moving_avg_window)
            codes_text += f"\nplot_n_annotate_lcf(lcf, ax, moving_avg_window={moving_avg_window_for_codes})"
        else:
            t0_to_use = t0 + step * period
            plot_transit(
                lcf,
                ax,
                t0_to_use,
                duration_hr / 24,
                surround_time,
                flux_col=flux_col,
                moving_avg_window=moving_avg_window,
                t0mark_ymax=t0mark_ymax,
                # fix legend to upper left to avoid clashing with the notebook nav at the upper right
                legend_kwargs=dict(loc="upper left"),
            )
            codes_text += f"""
#   transit parameters - t0: BTJD {t0}, duration: {duration_hr} hours, period: {period} days
plot_transit(lcf, ax, {t0_to_use}, {duration_hr} / 24, {surround_time}, \
moving_avg_window={moving_avg_window_for_codes}, t0mark_ymax={t0mark_ymax})

# transit_specs for calling plot_transits()
transit_specs = TransitTimeSpecList(
    dict(epoch={t0}, duration_hr={duration_hr}, period={period}, label="dip",
         sector={lcf.meta.get('SECTOR')}, steps_to_show=[{step}],
        ),
    defaults=dict(surround_time={surround_time})
)
"""

        ymin_to_use = ymin if ymin >= 0 else None
        ymax_to_use = ymax if ymax >= 0 else None
        if (ymin_to_use is not None) or (ymax_to_use is not None):
            ax.set_ylim(ymin_to_use, ymax_to_use)
            codes_text += f"""
# Zoom in on flux
ax.set_ylim({ymin_to_use}, {ymax_to_use})
"""
        widget_out2.clear_output()
        with widget_out2:
            print(codes_text)
        return None

    # Construct the UI

    desc_style = {"description_width": "25ch"}

    # Add a second output for textual
    widget_out2 = widgets.Output()

    if defaults is None:
        defaults = {}

    t0 = widgets.FloatText(
        value=defaults.get("epoch", -1),
        step=0.01,
        description=r"$t_{epoch}$, -1 for unspecified",
        style=desc_style,
    )
    duration_hr = widgets.FloatText(
        value=defaults.get("duration_hr", 1), step=0.1, description="duration (hours)", style=desc_style
    )
    period = widgets.FloatText(value=defaults.get("period", 999), step=0.01, description="period (days)", style=desc_style)
    step = widgets.IntText(
        value=defaults.get("step", 0), description=r"cycle (0 for transit at $t_{epoch}$)", style=desc_style
    )
    surround_time = widgets.FloatText(
        value=defaults.get("surround_time", 7), step=0.5, description="padding (days)", style=desc_style
    )
    moving_avg_window = widgets.Dropdown(
        options=[
            ("None", None),
            ("10 min", "10min"),
            ("20 min", "20min"),
            ("30 min", "30min"),
            ("1 hour", "1h"),
            ("2 hours", "2h"),
            ("4 hours", "4h"),
        ],
        value=defaults.get("moving_avg_window", "20min"),
        description="moving average window",
        style=desc_style,
    )
    ymin = widgets.FloatText(value=-1, step=0.1, description="flux min, -1 for default", style=desc_style)
    ymax = widgets.FloatText(value=-1, step=0.1, description="flux max, -1 for default", style=desc_style)
    t0mark_ymax = widgets.BoundedFloatText(
        value=0.05,
        step=0.05,
        min=0.0,
        max=1.0,
        description=r"$t_{epoch}$ mark height",
        style=desc_style,
    )
    VB = widgets.VBox
    HB = widgets.HBox
    ui = VB(
        [
            HB([t0, duration_hr, period]),
            HB([step, surround_time, moving_avg_window]),
            HB([ymin, ymax, t0mark_ymax]),
        ]
    )

    w = interactive_output(
        _update_plot_transit_interactive,
        dict(
            flux_col=fixed(flux_col),
            t0=t0,
            duration_hr=duration_hr,
            period=period,
            step=step,
            surround_time=surround_time,
            moving_avg_window=moving_avg_window,
            t0mark_ymax=t0mark_ymax,
            ymin=ymin,
            ymax=ymax,
            widget_out2=fixed(widget_out2),
        ),
    )

    w.layout.border = "1px solid lightgray"
    w.layout.padding = "1em 0px"

    widget_out2.layout.padding = "1em"

    display(ui, w, widget_out2)
    return w


def plot_flux_sap_flux_comparison(lc, sap_col="sap_flux", ax=None, offset=None, **kwargs):
    """Plot flux (typically PDCSAP_FLUX) and sap_flux together,
    to spot any anomaly in processed lightcurve."""
    lc_sap = lc.copy()
    lc_sap["flux"] = lc[sap_col]
    if sap_col + "_err" in lc.colnames:
        lc_sap["flux_err"] = lc[sap_col + "_err"]
    else:  # some products, e.g., QLP, does not give err
        # Hit a bug - ValueError: TessLightCurve object is invalid - expected 'time' as the first columns but found 'time'
        # lc_sap.remove_column('flux_err')
        # zero out the column as a workaround
        lc_sap["flux_err"] = np.zeros_like(lc_sap["flux_err"])

    if offset is None:
        # auto offset: move lc_sap curve so that
        # - its median is at about 10 percentile of main flux
        # - move farther down by a factor of the 40% amplitude of the main flux, so that it is most below the main flux
        #    without much gap, but not overlapping too much either.
        offset = (
            np.nanpercentile(lc.flux.value, 10)
            - np.nanmedian(lc_sap.flux.value)
            - (np.nanmedian(lc.flux.value) - np.nanpercentile(lc.flux.value, 10)) * 3
        )

    lc_sap.label += f" {sap_col} + {offset:.0f}"

    ax = LightCurveCollection([lc, lc_sap]).plot(ax=ax, offset=offset, **kwargs)
    ax.set_title(f"{lc.label}, sector {lc.sector} - flux vs {sap_col}")
    return ax


def mark_transit_times(
    lc,
    tt_specs,
    axvline_kwargs_specs=None,
    skip_no_transit_plot=False,
    mark_data_gap=False,
    legend_loc=None,
    lc_plot_func_name="scatter",
    ax=None,
):
    """Plot the given LC, and mark the transit times based on `tt_specs`."""
    break_tolerance = 5 if not mark_data_gap else 1e7
    tt_list = [
        lke.get_transit_times_in_lc(lc, a_spec["epoch"], a_spec["period"], break_tolerance=break_tolerance)
        for a_spec in tt_specs
    ]

    # skip if no transit found
    # (tt_list is a list of list, so it needs to be flattend for counting)
    if skip_no_transit_plot and len(np.array(tt_list, dtype=object).flatten()) < 1:
        print(f"{lc._repr_simple_()} is skipped - no matching transits.")
        return None, None

    # base plot
    #
    if ax is None:
        ax = lk_ax(figsize=(12, 4))
    ax = getattr(lc, lc_plot_func_name)(ax=ax, color="black", label=f"{lc.label} s.{getattr(lc, 'sector', 'N/A')}")
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(axis="x", which="minor", length=4)
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(axis="y", which="minor", length=4)

    # pre-process axvline_kwargs
    #
    if axvline_kwargs_specs is None:
        axvline_kwargs_specs = [dict(label="dip", linestyle="--", color="red")]

    # use the label in tt_specs if not specified in axvline_kwargs
    for (a_spec, an_axvline_kwargs, idx_0_based) in zip(tt_specs, axvline_kwargs_specs, range(len(axvline_kwargs_specs))):
        if an_axvline_kwargs.get("label") is None:
            an_axvline_kwargs["label"] = a_spec.get("label", f"dip {idx_0_based + 1}")

    # Mark transit times on the base plot
    #
    # a hack: mark the first line for each tt set, then set legend
    # so that each tt set will have 1 legend
    # if we simply set legend at the end, each dip will have its own legend!
    # TODO: use `vlines_y_in_axes_coord()` helper instead.
    for (transit_times, axvline_kwargs) in zip(tt_list, axvline_kwargs_specs):
        if len(transit_times) > 0 and axvline_kwargs is not None:
            axvline_kwargs = axvline_kwargs.copy()  # as we might need to modify them locally
            ymin, ymax = axvline_kwargs.pop("ymin", 0), axvline_kwargs.pop("ymax", 0.1)
            ax.axvline(transit_times[0], ymin, ymax, **axvline_kwargs)
    ax.legend(loc=legend_loc)

    for (transit_times, axvline_kwargs) in zip(tt_list, axvline_kwargs_specs):
        if axvline_kwargs is not None:
            axvline_kwargs = axvline_kwargs.copy()  # as we might need to modify them locally
            ymin, ymax = axvline_kwargs.pop("ymin", 0), axvline_kwargs.pop("ymax", 0.1)
            for tt in transit_times:
                ax.axvline(tt, ymin, ymax, **axvline_kwargs)

    return ax, tt_list


def normalize_lc_across_collection(lc, lcf_coll, **kwargs):
    """Normalize the given lightcurve across the collection of lightcurves."""

    # TODO: check the flux have the same unit, reject or convert if they don't have the same unit

    # calling np.nanmedain() directly over an array of fluxes does not work
    # np.nanmedain requires all flux object has the same length (thus forming a normal 2D array)
    median_all = np.nanmedian(np.concatenate([a_lc.flux.value for a_lc in lcf_coll]))
    median_lc = np.nanmedian(lc.flux.value)

    res = lc.normalize(**kwargs)
    res.flux = res.flux * (median_lc / median_all)
    res.flux_err = res.flux_err * (median_lc / median_all)

    return res


def break_vals_to_intervals(vals, min_v_diff):
    v_diff = np.diff(vals, prepend=vals[0])
    break_idxs = np.where(v_diff > min_v_diff)[0]
    # add begin to the indices
    break_idxs = np.concatenate(([0], break_idxs))

    def to_interval(i):
        i_start = break_idxs[i]
        if i < len(break_idxs) - 1:
            i_stop = break_idxs[i + 1]
        else:
            i_stop = len(vals)
        return (vals[i_start], vals[i_stop - 1])

    return np.array([to_interval(i) for i in np.arange(0, len(break_idxs))])


def intervals_to_ratio(intervals):
    def width(a_pair):
        res = a_pair[1] - a_pair[0]
        return res if res > 0 else 0.1  # ensure width is positive

    return np.array([width(a_pair) for a_pair in intervals])


def plot_skip_data_gap(lc, wspace=0.2, figsize=(16, 4), data_gap_min_days=10, **kwargs):
    """Plot a lightcurve, but compress large time gap.
    Use cases: for a lc stitched from multiple sectors, there are often huge data gap.
    This plot compresses the data gap.
    """

    # derived from broken axis example
    # https://matplotlib.org/stable/gallery/subplots_axes_and_figures/broken_axis.html

    # break the lc into a a few where there is a data gap between them
    # each is represented by an interval [time_start, time_stop_inclusive]
    intervals = break_vals_to_intervals(lc.time.value, min_v_diff=data_gap_min_days)
    width_ratio = intervals_to_ratio(intervals)

    num_plots = len(intervals)
    with plt.style.context(lk.MPLSTYLE):
        figsize = kwargs.pop("figsize", figsize)
        fig, axs = plt.subplots(1, num_plots, sharey=True, figsize=figsize, gridspec_kw={"width_ratios": width_ratio})
        if isinstance(axs, matplotlib.axes.Axes):
            # case it returns a single Axes object
            axs = [axs]

    # compress vertical spaces between plots
    fig.subplots_adjust(wspace=wspace)

    def create_tstart_sector_list(lc):
        res = []
        header_dict = lc.meta.get("HEADERS_ORIGINAL", {})  # a sector: meta dict
        for sector in header_dict:
            meta = header_dict[sector]
            tstart = meta.get("TSTART", None)
            if tstart is not None:
                res.append((tstart, sector))
        return res

    # to be consumed (and eventually emptied) in the for loop
    tstart_sector_list = create_tstart_sector_list(lc)

    xlabel_text = ""
    for i, (ax, interval) in enumerate(zip(axs, intervals)):
        lc.truncate(interval[0], interval[1] + 0.0001).scatter(ax=ax, **kwargs)

        # add sector start marker(s) to the current ax when in range
        while len(tstart_sector_list) > 0:
            tstart = tstart_sector_list[0][0]
            if interval[0] - 1 <= tstart <= interval[1]:
                # the current chunk is in range, mark tstart on graph
                ax.axvline(tstart, c="gray", alpha=0.4, linestyle="--")
                ax.text(  # annotate the line
                    tstart,
                    0.99,
                    f" S{tstart_sector_list[0][1]}",
                    transform=ax.get_xaxis_transform(),
                    va="top",
                    c="gray",
                )
                del tstart_sector_list[0]
            else:
                break

        # hide the spines between axs
        if 0 < i < num_plots - 1:
            ax.spines["left"].set_visible(False)
            ax.spines["right"].set_visible(False)
        elif i == 0:
            ax.spines["right"].set_visible(False)
        else:
            ax.spines["left"].set_visible(False)

        # hide the axis/and ticks
        if i > 0:
            ax.yaxis.set_visible(False)

        # Now, let's turn towards the cut-out slanted lines.
        # We create line objects in axes coordinates, in which (0,0), (0,1),
        # (1,0), and (1,1) are the four corners of the axes.
        # The slanted lines themselves are markers at those locations, such that the
        # lines keep their angle and position, independent of the axes size or scale
        # Finally, we need to disable clipping.

        d = 0.5  # proportion of vertical to horizontal extent of the slanted line
        marker_kwargs = dict(
            marker=[(-1, -d), (1, d)],
            markersize=12,
            linestyle="none",
            color="darkgray",
            alpha=0.4,
            mec="k",
            mew=1,
            clip_on=False,
        )
        if i < num_plots - 1:
            ax.plot(1, 0, transform=ax.transAxes, **marker_kwargs)  # the right bottom mark
            ax.plot(1, 1, transform=ax.transAxes, **marker_kwargs)  # the right top mark
        if i > 0:
            ax.plot(0, 0, transform=ax.transAxes, **marker_kwargs)  # the left bottom mark
            ax.plot(0, 1, transform=ax.transAxes, **marker_kwargs)  # the left top mark

        # 1 xlabel for entire plot that will be constructed later
        xlabel_text = ax.get_xlabel()
        ax.set_xlabel("")

        # 1 legend for entire plot
        if i < num_plots - 1:
            ax.legend().remove()
        else:
            ax.legend(loc="upper right")

    # use a special text object for the single xlabel
    # so that it can be positioned at the center
    custom_xlabel_text_obj = fig.text(0.5, 0, xlabel_text, ha="center")
    # propagagate the text properties from the native xlabel to our custom oe
    for prop in ["color", "fontsize", "fontproperties", "fontfamily", "fontname", "fontstyle", "fontweight", "fontvariant"]:
        # generalized version of: custom_xlabel_text_obj.set_fontsize(axs[0].xaxis.get_label().get_fontsize())
        getattr(custom_xlabel_text_obj, f"set_{prop}")(getattr(axs[0].xaxis.get_label(), f"get_{prop}")())

    # OPEN: create a wrapper over the axs, say, `HorizontalAxesCollection`, such that users have convenience plot methods
    # as if it is a single plot object, e.g.,
    # - axvline(): find the right ax to forward the call
    # - get/set_xlabel(): forward to the custom xlabel text object `custom_xlabel_text_obj`
    # - get/set_ylabel(): forward to axs[0]
    # - set_title(): use figure

    return axs


def plot_centroids(lc, transit_time=None, window=None):

    # based on:
    # https://github.com/noraeisner/PH_Coffee_Chat/blob/4a723030ed80eeabdfb0eca49d948b97d61e35f6/False%20Positive/False%20positives%20-%20(4)%20centroid-position.ipynb

    # bin the data, skip it to speed up the plot due to astropy performance regression
    # lc_bin = lc.bin(7/60/24)
    lc_bin = lc

    # generate a mask so that we only see the times around the transit event
    if (window is None) or (transit_time is None):
        transit_mask = lc_bin.time.value > 0
    else:
        transit_mask = (lc_bin.time.value > transit_time - window) & (lc_bin.time.value < transit_time + window)

    # make a plot with three panels so that we can see the lightcurve and the centroid positions
    with plt.style.context(lk.MPLSTYLE):
        fig, ax = plt.subplots(3, 1, figsize=(8, 5), sharex=True)

    # plot the lightcurve in the top panel (in orange)
    ax[0].plot(
        lc_bin.time.value[transit_mask], lc_bin.sap_flux.value[transit_mask], color="darkorange", lw=0, marker=".", ms=3
    )

    # plot the centroid motions in the column direction in the middle panel
    ax[1].plot(
        lc_bin.time.value[transit_mask],
        lc_bin.mom_centr1.value[transit_mask] - np.nanmean(lc_bin.mom_centr1.value[transit_mask]),
        color="black",
        lw=0,
        marker=".",
        ms=2,
        alpha=0.5,
    )
    ax[1].plot(
        lc_bin.time.value[transit_mask],
        lc_bin.pos_corr1.value[transit_mask] - np.nanmean(lc_bin.pos_corr1.value[transit_mask]),
        color="red",
        lw=0,
        marker=".",
        ms=2,
        alpha=0.5,
    )

    # plot the centroid motions in the row direction in the middle panel
    ax[2].plot(
        lc_bin.time.value[transit_mask],
        lc_bin.mom_centr2.value[transit_mask] - np.nanmedian(lc_bin.mom_centr2.value[transit_mask]),
        color="black",
        lw=0,
        marker=".",
        ms=2,
        alpha=0.5,
        label="Brightness motion",
    )
    ax[2].plot(
        lc_bin.time.value[transit_mask],
        lc_bin.pos_corr2.value[transit_mask] - np.nanmedian(lc_bin.pos_corr2.value[transit_mask]),
        color="red",
        lw=0,
        marker=".",
        ms=2,
        alpha=0.5,
        label="Satellite motion",
    )

    if transit_time != None:
        # draw a vertical line at the time of the transit event
        ax[0].axvline(transit_time, color="grey", zorder=-1)
        ax[1].axvline(transit_time, color="grey", zorder=-1)
        ax[2].axvline(transit_time, color="grey", zorder=-1)

    # label the axes
    ax[0].set_ylabel("Flux")
    ax[1].set_ylabel("Column")
    ax[2].set_ylabel("Row")
    plt.xlabel(f"Time ({lc.time.format.upper()})")

    ax[2].xaxis.set_major_formatter(FormatStrFormatter("%.2f"))

    if (window is not None) and (transit_time is not None):
        plt.xlim(transit_time - window, transit_time + window)

    plt.tight_layout()
    plt.legend()
    plt.show()

    return ax


def scatter_centroids(
    lcf,
    fig=None,
    highlight_time_range=None,
    time_range=None,
    c="blue",
    c_highlight="red",
):
    """
    Scatter centroids, and highlight the specific time range
    """

    if fig is None:
        fig = plt.figure(figsize=(12, 12))

    lc = _normalize_to_percent_quiet(lcf)
    sector = lcf.meta.get("SECTOR")

    if time_range is not None:
        lc = lc.truncate(time_range[0], time_range[1])

    fig.gca().yaxis.set_major_formatter(FormatStrFormatter("%.3f"))  # avoid scientific notations
    fig.gca().scatter(lc.centroid_col.value, lc.centroid_row.value, c=c, label=f"TIC {lc.targetid}")

    if highlight_time_range is not None:
        lc_highlight = lc.truncate(highlight_time_range[0], highlight_time_range[1])
        if len(lc_highlight) < 1:
            print("WARNING: scatter_centroids() no observations in highlight_time_range")
        fig.gca().scatter(
            lc_highlight.centroid_col.value,
            lc_highlight.centroid_row.value,
            c=c_highlight,
            label="highlights",
        )

    title = f"TIC {lc.targetid} Centroids, sector {sector}"
    if time_range is not None:
        title += f"\n{as_4decimal(time_range)}"
    if highlight_time_range is not None:
        title += f"\nHighlights:{as_4decimal(highlight_time_range)}"
    fig.gca().set_title(title)
    fig.legend()
    return fig


def _update_anim(n, ax, lc, label, num_centroids_to_show, use_relative_time, c):
    ax.cla()
    # fix the x/y scale to ensure it doesn't change over the animation
    c_col, c_row = _to_unitless(lc.centroid_col), _to_unitless(lc.centroid_row)
    ax.set_xlim(np.nanmin(c_col), np.nanmax(c_col))
    ax.set_ylim(np.nanmin(c_row), np.nanmax(c_row))

    # avoid scientific notation for y-axis
    # x-axis might need scientific notation so that the labels won't get too cramped with long decimals
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.3f"))

    if num_centroids_to_show is None:
        col = lc.centroid_col[:n]
        row = lc.centroid_row[:n]
        time_label = f"{as_4decimal(lc.time[n])}"
        if use_relative_time:
            time_label = time_label + f" ({as_4decimal(lc.time_rel[n])})"
    else:
        n_start = max(0, n - num_centroids_to_show)
        col = lc.centroid_col[n_start:n]
        row = lc.centroid_row[n_start:n]
        time_label = f"{as_4decimal(lc.time[n_start])} - {as_4decimal(lc.time[n])}"
        if use_relative_time:
            time_label = time_label + f" ({as_4decimal(lc.time_rel[n_start])} - {as_4decimal(lc.time_rel[n])})"

    ax.set_title(f"TIC {lc.targetid} Centroids, {label}\nday: {time_label}")
    ax.scatter(col, row, c=c)


def animate_centroids(
    lcf,
    fig=None,
    frames=None,
    num_obs_per_frame=240,
    interval=250,
    use_relative_time=False,
    time_range=None,
    accumulative=True,
    c=None,
    display=True,
):
    """
    Animate centroids to visualize changes over time.

    """
    lc = lcf
    label = f"sector {lcf.meta.get('SECTOR')}"

    # Zoom to a particular time range if specified
    if time_range is not None:
        # use pandas to zoom to a particular time_range
        df = _normalize_to_percent_quiet(lc).to_pandas(columns=["time", "flux", "centroid_row", "centroid_col"])
        df = df[(df.time >= time_range[0]) & (df.time <= time_range[1])]
        if len(df) < 1:
            raise Exception(f"Zoomed lightcurve has no observation. time_range={time_range}")

        lc_z = SimpleNamespace()  # zoomed-in lightcurve-like object for the purpose of animation
        setattr(lc_z, "time", df.time.values)
        setattr(lc_z, "flux", df.flux.values)
        setattr(lc_z, "centroid_row", df.centroid_row.values)
        setattr(lc_z, "centroid_col", df.centroid_col.values)
        setattr(lc_z, "targetid", lc.targetid)
        lc = lc_z

    if fig is None:
        fig = plt.figure(figsize=(12, 12))
    if frames is None:
        num_obs = len(lc.centroid_row)
        num_frames = int(num_obs / num_obs_per_frame)  # default 240 is about every 8 hours, given 2-minute intervals
        ary_n = np.linspace(1, num_obs, num=num_frames, endpoint=False)
        ary_n[0] = np.ceil(ary_n[1] / 2)
        ary_n = list(map(lambda n: int(n), ary_n))
    else:
        ary_n = frames
        num_obs_per_frame = frames[-1] - frames[-2]  # assume the spacing of input is linear

    num_centroids_to_show = num_obs_per_frame
    if accumulative:
        num_centroids_to_show = None
    #     print(f'Steps: {ary_n}')
    if use_relative_time:
        rel_time_added = add_relative_time(lc, lcf)
        if not rel_time_added:
            use_relative_time = False
    anim = animation.FuncAnimation(
        fig,
        _update_anim,
        frames=ary_n,
        fargs=(fig.gca(), lc, label, num_centroids_to_show, use_relative_time, c),
        interval=interval,
        blit=False,
    )
    if display:
        # for inline display in jupyter
        try:
            from IPython.display import HTML
            from IPython.display import display as iDisplay

            return iDisplay(HTML(anim.to_jshtml(default_mode="once")))
        except ImportError:
            print("WARNING: animate_centroids() - inline display not possible Not in IPython environment.")
            return anim


def markTimes(ax, times, **kwargs):
    """Helper to mark specifics time as vertical lines on a plot"""
    axvline_kwargs = kwargs.copy()
    # apply defaults
    axvline_kwargs.setdefault("c", "gray")
    axvline_kwargs.setdefault("linewidth", 1)
    axvline_kwargs.setdefault("linestyle", "--")
    for t in times:
        ax.axvline(t, **axvline_kwargs)


def plot_n_annotate_folded(lc_f, figsize=(10, 5), annotate=True, also_plot_zoom_transit=False, duration=None):
    ax = lk_ax(figsize=figsize)
    scatter(lc_f, ax=ax, s=1)

    if annotate:
        epoch = lc_f.meta.get("EPOCH_TIME")  # a time object
        period = lc_f.meta.get("PERIOD").to(u.day).value
        obs_span_days, obs_actual_days = lke.get_obs_date_range(lc_f)
        obs_span_cycles = obs_span_days / period
        obs_actual_cycles = len(set(lc_f.cycle))

        ax.set_title(
            f"""{lc_f.label} folded, period={period:.4f}d, epoch={epoch.format.upper()} {epoch.value:.3f}
    {lc_f.time_original.format.upper()} {lc_f.time_original.min().value:.2f} - {lc_f.time_original.max().value:.2f}
    time span: {obs_span_days:.2f}d / {obs_span_cycles:.0f} cycles
    observation: {obs_actual_days}d / {obs_actual_cycles} cycles
    """
        )

    if not also_plot_zoom_transit:
        return ax

    ax_z = lk_ax(figsize=figsize)
    scatter(lc_f.truncate(-duration * 1.5, duration * 1.5), ax=ax_z)
    ax_z.axvline(x=0, c="red")
    ax_z.axvspan(0 - duration / 2, 0 + duration / 2, color="red", alpha=0.2)

    return ax, ax_z


def fold_and_plot(lc, period, epoch_time, flux_in_mag=False, **kwargs):
    lc_f = lc.fold(period=period, epoch_time=epoch_time, epoch_phase=0)
    if flux_in_mag:
        lc_f = lke.lc_to_flux_in_mag_by_normalization(lc_f)
    return plot_n_annotate_folded(lc_f, **kwargs), lc_f


def fold_and_plot_odd_even(lc, period, epoch_time, figsize=(10, 5), title_extra=""):
    lc_folded = lc.fold(period=period, epoch_time=epoch_time, epoch_phase=0)

    ax = lk_ax(figsize=figsize)
    lc_f_odd = lc_folded[lc_folded.odd_mask]
    lc_f_odd.scatter(ax=ax, c="r", label="odd", marker=".", s=4)
    lc_f_even = lc_folded[lc_folded.even_mask]
    lc_f_even.scatter(ax=ax, c="b", label="even", marker="x", s=4)

    pct01_odd = np.nanpercentile(lc_f_odd.flux, 0.1)
    pct01_even = np.nanpercentile(lc_f_even.flux, 0.1)

    ax.axhline(
        pct01_odd * 100,
        c="r",
        linestyle="--",
        label=f"odd 0.1 pctile {pct01_odd:0.4f}",
    )
    ax.axhline(
        pct01_even * 100,
        c="b",
        linestyle="dotted",
        label=f"even 0.1 pctile {pct01_even:0.4f}",
    )

    ax.legend()
    plt.title(f"{lc.label} folded {title_extra}\nperiod={period:.4f} d")

    print("odd  0.1 percentile: ", pct01_odd)
    print("even 0.1 percentile: ", pct01_even)
    return ax, lc_folded


def fold_2x_periods_and_plot(lc, period, epoch_time, figsize=(10, 5), title_extra=""):
    lc_folded = lc.fold(period=period * 2, epoch_time=epoch_time, epoch_phase=period / 2)

    ax = lk_ax(figsize=figsize)
    lc_folded.scatter(ax=ax)

    ax.legend()
    ax.xaxis.set_label_text(
        ax.xaxis.get_label_text() + f" , {lc.time.format.upper()} {lc.time.min().value:.2f} - {lc.time.max().value:.2f}"
    )
    plt.title(f"{lc.label} folded at 2X periods {title_extra}\nperiod={period:.4f} d")

    return ax, lc_folded


def calc_cycles(lc: FoldedLightCurve):
    cycle_epoch_start = lc.epoch_time - lc.period / 2
    cycles = np.asarray(np.floor(((lc.time_original - cycle_epoch_start) / lc.period).value), dtype=int)
    # the cycle where epoch is set as 0, adjust it so that the first cycle is 0
    cycles = cycles - cycles.min()

    return cycles


def animate_folded_lightcurve(lc: FoldedLightCurve, ax=None, num_frames=10, interval=1000, plot_kwargs={}):
    def _update_anim(frame, ax, lc, cycle_column, cycle_list, num_cycles_per_plot, plot_kwargs):
        cycle_list_subset = cycle_list[num_cycles_per_plot * frame : num_cycles_per_plot * (frame + 1)]

        ax.cla()
        lc_subset = lc[np.in1d(cycle_column, cycle_list_subset)]
        lc_subset.scatter(ax=ax, **plot_kwargs)

        # ensure all plots have the same scale
        ax.set_xlim(lc.time.min().value, lc.time.max().value)
        ax.set_ylim(lc.flux.min().value, lc.flux.max().value)

        # set title
        if len(cycle_list_subset) < 2:
            cycle_list_subset_label = cycle_list_subset[0]
        else:
            cycle_list_subset_label = f"{np.min(cycle_list_subset)} - {np.max(cycle_list_subset)}"
        ax.set_title(
            f"""{lc.label} , cycle {cycle_list_subset_label}, {lc_subset.time_original.format.upper()} \
{lc_subset.time_original.min().value:.4f} - {lc_subset.time_original.max().value:.4f}"""
        )

    lc = lc.remove_nans()  # remove any potential empty frames with no flux
    if ax is None:
        ax = lk_ax()

    cycle_column = calc_cycles(lc)
    cycle_list = np.unique(cycle_column)
    cycle_list.sort()

    num_cycles_per_plot = int(np.ceil(len(cycle_list) / num_frames))
    num_frames = int(np.ceil(len(cycle_list) / num_cycles_per_plot))

    cycle_idxs = list(range(0, num_frames))
    anim = animation.FuncAnimation(
        ax.get_figure(),
        _update_anim,
        frames=cycle_idxs,
        fargs=(ax, lc, cycle_column, cycle_list, num_cycles_per_plot, plot_kwargs),
        interval=interval,
        blit=False,
    )

    # for inline display in jupyter
    try:
        from IPython.display import HTML
        from IPython.display import display as iDisplay

        return iDisplay(HTML(anim.to_jshtml(default_mode="once")))
    except ImportError:
        print("WARNING: animate_folded_lightcurve() - display not possible in non-IPython environment.")
        return anim


from bokeh.plotting import ColumnDataSource


class MarkListUiModel(object):
    def __init__(self, mark_list, plot_height, upper_fraction=0.2, lower_fraction=0, data_source_attr_list=["time"]):
        self._mark_list = mark_list
        self._data_source_attr_list = data_source_attr_list
        self.plot_height = plot_height
        self.upper_fraction = upper_fraction
        self.lower_fraction = lower_fraction
        # it needs to be initialized last
        self._source = ColumnDataSource(data=self._to_dict())

    def _to_dict(self, attr_list=None):
        if attr_list is None:
            attr_list = self._data_source_attr_list
        res = dict()
        for attr in attr_list:
            res[attr] = [m.get(attr, None) for m in self._mark_list]

        # upper/lower is needed for UI, in screen unit
        upper = self.plot_height * self.upper_fraction
        lower = self.plot_height * self.lower_fraction
        res["upper"] = [upper for i in range(0, len(self._mark_list))]
        res["lower"] = [lower for i in range(0, len(self._mark_list))]

        return res

    def append(self, mark):
        self._mark_list.append(mark)
        # could possibly be optimized by using stream()
        self._source.data = self._to_dict()

    def del_by_idx(self, idx):
        del self._mark_list[idx]
        self._source.data = self._to_dict()

    @property
    def source(self):
        return self._source


def interact(
    lc: LightCurve,
    ylim_func: Optional(LC_Ylim_Func_Type) = None,
    plot_height: float = 490,
    plot_width: float = 900,
    show_line: bool = False,
    notebook_url: str = "localhost:8888",
) -> SimpleNamespace:
    from bokeh.plotting import ColumnDataSource, output_notebook, show
    from bokeh.layouts import layout, row, column, Spacer
    from bokeh.models import Button, Div, Whisker, TextInput, Dropdown
    from bokeh.models.tools import BoxZoomTool, WheelZoomTool, UndoTool, RedoTool
    from astropy.table import Table

    mark_list = []  # to be returned, so it needs to be deined at the top

    def get_tool_of_class(toolbar_or_fig, cls):
        if hasattr(toolbar_or_fig, "toolbar"):
            tools = toolbar_or_fig.toolbar.tools
        else:
            tools = toolbar_or_fig.tools
        for t in tools:
            if isinstance(t, cls):
                return t
        return None

    def create_interact_ui(doc):
        lc_source = lk.interact.prepare_lightcurve_datasource(lc)
        fig_lc, vertical_line = lk.interact.make_lightcurve_figure_elements(lc, lc_source, ylim_func=ylim_func)
        fig_lc.output_backend = "webgl"  # use GPU accelerated graphics when possible

        # customize the plot to make it more suitable for our purpose
        # hack: assume the renderers are in specific order
        #       can be avoided if the renderers have name when they are created.
        # r_lc_step = [r for r in fig_lc.renderers if r.name == "lc_step"][0]
        r_lc_step = fig_lc.renderers[0]
        r_lc_step.visible = show_line

        # r_lc_circle = [r for r in fig_lc.renderers if r.name == "lc_circle"][0]
        r_lc_circle = fig_lc.renderers[1]
        r_lc_circle.glyph.fill_color = "gray"
        r_lc_circle.glyph.fill_alpha = 1.0
        r_lc_circle.nonselection_glyph.fill_color = "gray"
        r_lc_circle.nonselection_glyph.fill_alpha = 1.0

        fig_lc.plot_height = plot_height
        fig_lc.plot_width = plot_width
        fig_lc.toolbar.active_drag = get_tool_of_class(fig_lc, BoxZoomTool)
        fig_lc.toolbar.active_scroll = get_tool_of_class(fig_lc, WheelZoomTool)
        fig_lc.toolbar.active_inspect = None
        fig_lc.add_tools(UndoTool(), RedoTool())

        #
        # UI to select a specific point (with vertical_line as the marker)
        #
        vertical_line.visible = False

        def jump_to_lightcurve_position(attr, old, new):
            if new == []:
                return

            time = lc_source.data["time"][new[0]]
            vertical_line.update(location=time)
            vertical_line.visible = True
            lc_source.selected.indices = [new[0]]

            # pan/zoom to the selected mark if it's out of range
            if time < fig_lc.x_range.start or time > fig_lc.x_range.end:
                fig_lc.x_range.start = time - 1.5
                fig_lc.x_range.end = time + 1.5

            flux = lc_source.data["flux"][new[0]]
            current_point_info_div.text = f"({time:.3f}, {flux:.3f})"

            if len(mark_list) > 0:
                duration = abs(time - mark_list[-1]["time"])
                flux_delta = flux - mark_list[-1]["flux"]
                delta_info_text = f"({duration:.3f}d, {flux_delta:.3f})"
            else:
                delta_info_text = ""

            delta_info_div.text = delta_info_text
            delta_label_div.visible = delta_info_text != ""

        lc_source.selected.on_change("indices", jump_to_lightcurve_position)

        #
        # Mark-related Glyphs, widgets and callbacks
        #
        select_mark_dropdown = Dropdown(label="Select mark", disabled=True, width=200)
        out_div_default_text = "Marked times to be shown here"
        out_div = Div(text=out_div_default_text)

        current_point_info_div = Div(style={"font-family": "monospace"}, width=150)
        delta_label_div = Div(text="Delta from last mark:", visible=False)
        delta_info_div = Div(style={"font-family": "monospace"})

        mark_btn = Button(label="Mark", button_type="default", width=100)
        mark_label_input = TextInput(width=200, placeholder="Optional label for the mark")

        # UI(Glyph) and UI Model for marks
        marks_ui_model = MarkListUiModel(mark_list, fig_lc.plot_height, lower_fraction=0, upper_fraction=0.2)
        marks_whisker = Whisker(
            base="time",
            upper="upper",
            lower="lower",
            upper_units="screen",
            lower_units="screen",
            dimension="height",
            source=marks_ui_model.source,
            level="annotation",
            upper_head=None,
            lower_head=None,
            line_color="red",
            line_width=4,
            line_dash="dashed",
        )
        fig_lc.add_layout(marks_whisker)

        # define callbacks (where a bokeh server is needed, for Python callback)

        def update_mark_btn_ui(attr, old, new):
            # don't care about attr, old, new, but needed to work as bokeh callback
            idx_in_mark_list = idx_of_selected_in_mark_list()
            if idx_in_mark_list < 0:
                mark_btn.label = "Mark"
                mark_label_input.disabled = False
            else:
                mark_btn.label = "Un-mark"
                mark_label_input.disabled = True

        lc_source.selected.on_change("indices", update_mark_btn_ui)
        marks_ui_model.source.on_change("data", update_mark_btn_ui)

        def update_select_mark_dropdown_ui(attr, old, new):
            # don't care about attr, old, new, but needed to work as bokeh callback
            time = [mark["time"] for mark in mark_list]
            select_mark_dropdown.menu = [(f"{t:.3f}", str(i)) for i, t in enumerate(time)]
            select_mark_dropdown.disabled = True if len(time) < 1 else False

        marks_ui_model.source.on_change("data", update_select_mark_dropdown_ui)

        def on_select_mark_from_dropdown(event):
            idx_of_in_mark_list = int(event.item)
            mark = mark_list[idx_of_in_mark_list]
            idx_of_mark_in_lc_source = mark["idx"]
            lc_source.selected.indices = [idx_of_mark_in_lc_source]

        select_mark_dropdown.on_click(on_select_mark_from_dropdown)

        def idx_of_selected_in_mark_list():
            """Return the index of the selected point in `mark_list`, -1 if it is not in the list.
            Used by `toggle_selected_mark()`
            """
            if len(lc_source.selected.indices) < 1:
                return -1
            idx = lc_source.selected.indices[0]
            selected_time = lc_source.data["time"][idx]
            match_indices = [i for i, mark in enumerate(mark_list) if mark["time"] == selected_time]
            if len(match_indices) < 1:
                return -1
            return match_indices[0]

        def toggle_selected_mark():
            if len(lc_source.selected.indices) < 1:
                return
            idx_in_mark_list = idx_of_selected_in_mark_list()
            if idx_in_mark_list < 0:
                idx = lc_source.selected.indices[0]
                mark = dict(time=lc_source.data["time"][idx], flux=lc_source.data["flux"][idx], idx=idx)
                label = mark_label_input.value
                if label != "":
                    mark["label"] = label

                marks_ui_model.append(mark)
            else:
                marks_ui_model.del_by_idx(idx_in_mark_list)

            out_div.text = Table(mark_list)._repr_html_() if len(mark_list) > 0 else out_div_default_text

        mark_btn.on_click(toggle_selected_mark)

        #
        # Pan-related widgets / callbacks
        #
        rr_button = Button(label=">>", button_type="default", width=30)
        ll_button = Button(label="<<", button_type="default", width=30)
        pan_amount_input = TextInput(width=100, placeholder="default: plot width")

        def get_pan_width():
            try:
                return float(pan_amount_input.value)
            except ValueError:
                # case empty string (or any invalid value for now)
                return fig_lc.x_range.end - fig_lc.x_range.start

        def pan_left():
            width = get_pan_width()
            fig_lc.x_range.start = fig_lc.x_range.start - width
            fig_lc.x_range.end = fig_lc.x_range.end - width

        ll_button.on_click(pan_left)

        def pan_right():
            width = get_pan_width()
            fig_lc.x_range.start = fig_lc.x_range.start + width
            fig_lc.x_range.end = fig_lc.x_range.end + width

        rr_button.on_click(pan_right)

        #
        # Overall layout
        #
        doc_layout = row(
            column(
                fig_lc,
                row(
                    current_point_info_div,
                    delta_label_div,
                    delta_info_div,
                    Spacer(width=30),
                    ll_button,
                    rr_button,
                    pan_amount_input,
                    Spacer(width=30),
                    mark_btn,
                    mark_label_input,
                    align="end",
                ),
            ),
            Spacer(width=30),
            column(select_mark_dropdown, out_div),
        )

        doc.add_root(doc_layout)

    # main logic
    output_notebook(verbose=False, hide_banner=True)
    show(create_interact_ui, notebook_url=notebook_url)
    # Use a namespace rather than just returning a mark_list, to give flexibility in case more information is to be returned
    return SimpleNamespace(mark_list=mark_list)


def plot_spoc_aperture(lc, ax=None):
    """Plot the aperture associated with the SPOC lightcurve"""
    return _plot_lc_aperture(lc, "SPOC", 2, ax)


def plot_tasoc_aperture(lc, ax=None):
    """Plot the aperture associated with the TASOC lightcurve"""
    return _plot_lc_aperture(lc, "TASOC", 3, ax)


def _plot_lc_aperture(lc, supported_author, aperture_hdu_idx, ax=None):
    """Plot the aperture associated with the lightcurve.
    The specifics depends on the pipeline produced the lightcurve."""
    if lc.meta.get("AUTHOR") != supported_author:
        warnings.warn(f"The given lightcurve is not a {supported_author} lightcurve. No-op")
        return ax

    with fits.open(lc.filename) as hdu:
        # it'd be great is I could show the actual CCD row/column, but the data is not there.
        ax = lk.utils.plot_image(hdu[aperture_hdu_idx].data, title=lc.meta.get("LABEL"), ax=ax)
        return ax


def plot_tasoc_pixels_n_aperture(lc, ax=None):
    """Plot the time-average pixels, and the aperture associated with the TASOC lightcurve"""

    def draw_aperture_mask(aperture_mask, ax, mask_color):
        # adapted from TargetPixelFile.plot()
        from matplotlib import patches

        # Overlay the aperture mask if given
        if aperture_mask is None:
            return ax
        for i in range(aperture_mask.shape[0]):
            for j in range(aperture_mask.shape[1]):
                if aperture_mask[i, j]:
                    rect = patches.Rectangle(
                        xy=(j - 0.5, i - 0.5),
                        width=1,
                        height=1,
                        color=mask_color,
                        fill=False,
                        hatch="//",
                    )
                    ax.add_patch(rect)
        return ax

    if lc.meta.get("AUTHOR") != "TASOC":
        warnings.warn("The given lightcurve is not a TASOC lightcurve. No-op")
        return ax

    # for TASOC lightcurve
    # - image of time-average pixels is in hdu[2]
    # - aperture is in hdu[3]
    # https://archive.stsci.edu/hlsp/tasoc
    with fits.open(lc.filename) as hdu:
        # it'd be great is I could show the actual CCD row/column, but the data is not there.
        ax = lk.utils.plot_image(hdu[2].data, title=lc.meta.get("LABEL"), ax=ax)
        aperture = hdu[3].data

        # TASOC aperture has multiple values, the exact meaning is unclear. It seems that:
        # - the one with the highest value is the aperture
        # - the one with 2nd highest is the background
        unique_vals = np.unique(aperture)
        unique_vals[::-1].sort()

        val = unique_vals[0]
        draw_aperture_mask((aperture == val), ax, "red")
        #         val = unique_vals[1]  # background?
        #         draw_aperture_mask((aperture == val), ax, 'white')

        return ax


def scatter_partition_by(lc, partition_by_column, ax=None, **kwargs):
    """Generate a scatter plot of the given lightcurve, with flux
    partitioned by the given column.

    Use cases: provide a plot where flux comes from multiple bands / cameras.
    """
    if ax is None:
        ax = lk_ax()

    for val in np.unique(lc[partition_by_column]):
        lc[lc[partition_by_column] == val].scatter(ax=ax, label=f"{partition_by_column} {val}", **kwargs)
    ax.set_title(lc.label)
    return ax


#
# TargetPixelFile helpers
#


def show_tpf_orientation(tpf):
    """ "Helper to visualize the TPF's orientation in the sky. Requires IPython.
    Long arm is north, short arm with arrow is east.
    """
    coord_bottom_left = tpf.wcs.pixel_to_world(0, 0)
    coord_upper_right = tpf.wcs.pixel_to_world(tpf.shape[2] - 1, tpf.shape[1] - 1)
    coord_upper_left = tpf.wcs.pixel_to_world(0, tpf.shape[2] - 1)
    deg_from_north = coord_bottom_left.position_angle(coord_upper_left).to(u.deg).value

    display(
        HTML(
            f"""<div style="position: relative; margin-left: 16px;height: 64px;">
    <div title="Long arm: North; Short arm with arrow: East"
         style="float: left; max-width: 64px;font-size: 32px;margin: 16px;\
transform: rotate({-deg_from_north}deg);transform-origin: left; cursor:pointer;">↳</div>
        <div style="font-family: monospace;">Upper right offset from bottom left - <br>
        RA: {(coord_upper_right.ra - coord_bottom_left.ra).to(u.arcmin):0.6},
        Dec: {(coord_upper_right.dec - coord_bottom_left.dec).to(u.arcmin):0.6}
        </div>
    </div>"""
        )
    )


def interact_sky(tpf, notebook_url="localhost:8888", aperture_mask="empty", magnitude_limit=18):
    """tpf.interact_sky wrapper to handle different lightkurve versions."""
    if "aperture_mask" in inspect.getfullargspec(tpf.interact_sky).args:
        # case using a pre-release lightkurve that supports aperture_mask
        return tpf.interact_sky(notebook_url=notebook_url, aperture_mask=aperture_mask, magnitude_limit=magnitude_limit)
    else:
        # using release lightkurve that not yet supports aperture_mask
        return tpf.interact_sky(notebook_url=notebook_url, magnitude_limit=magnitude_limit)


def show_nearby_tic_summary_form():
    """Display a form that create a 1-line summary of a nearby TIC from the first three rows of the selected TIC table"""
    display(
        HTML(
            r"""
First 3 rows of the TIC info table:<br>
<textarea id="inStarInfo" style="width: 40ch; height: 6em;" placeholder="TIC \t12345678\nTESS Mag \t10.123\nSeparation ...">
</textarea><br>
<button id="ctlStarInfo">Create Nearby TIC summary</button>
<input id="outStarInfo" style="width: 40ch;" value="" readonly>

<script>

function convertToMultiLinePlaceholder(elem) { // create multiline placeholder
  elem.placeholder = elem.placeholder.replace(/\\n/g, '\n');
  elem.placeholder = elem.placeholder.replace(/\\t/g, '\t');
}
convertToMultiLinePlaceholder(document.querySelector('#inStarInfo'));


function createNearbyTicSummary(text) {
const toCells = (line) => {
    return line.split("\t");
}

const lines = text.split(/[\r\n]+/);


const ticId = toCells(lines[0])[1].replace(/^\s*(\d+).*$/, '$1');
const tessMag = toCells(lines[1])[1];
const separation = toCells(lines[2])[1];

return `TIC ${ticId} (TESS magnitude ${tessMag}, ${separation} arcsec away)`;
}

document.querySelector('#ctlStarInfo').onclick = (evt) => {
const summary = createNearbyTicSummary(document.querySelector('#inStarInfo').value);
document.querySelector('#outStarInfo').value = summary;
};
</script>
"""
        )
    )


def plot_with_aperture_n_background(
    tpf: lk.targetpixelfile.TargetPixelFile, aperture_mask=None, background_mask=None, ax=None, title_extra=""
) -> matplotlib.axes.Axes:
    """Plot the given targetpixel file, with  aperture mask and background_mask."""
    if aperture_mask is None:
        aperture_mask = tpf.pipeline_mask

    if background_mask is None:
        background_mask = tpf.background_mask

    ax = tpf.plot(ax=ax, aperture_mask=aperture_mask)
    ax = tpf.plot(ax=ax, aperture_mask=background_mask, mask_color="white", show_colorbar=False)
    ax.set_title(f"{getattr(tpf, 'targetid', '')}{title_extra}")
    return ax


def plot_in_out_diff(tpf, epoch, transit_half_duration=0.25, oot_outer_relative=0.5, oot_inner_relative=0.3, plot_lc=True):
    """
    Plot the in transit average flux and the out of transit average flux and compare the two (difference image).
    """

    # based on plot_in_out_TPF() in
    # https://github.com/noraeisner/PH_Coffee_Chat/blob/4a723030ed80eeabdfb0eca49d948b97d61e35f6/False%20Positive/False%20positives%20-%20(2)%20in%20out%20transit%20flux.ipynb

    tpf_list = [tpf.flux.value]
    t_list = [tpf.time.value]
    T0_list = [epoch]

    plt.figure(figsize=(9, 2.5 * len(T0_list)))

    plt.tight_layout()

    # keep track of how many images have been plotted to that they appear on a subgrid of plots which has three columns
    count = 0

    # loop through all of the list of PCA corrected flux vs time arrays for each marked transit-event
    for idx, tpf_filt in enumerate(tpf_list):  # idx is for each maked transit-event

        T0 = T0_list[idx]  # the time of the transit-like event
        t = t_list[idx]  # the time array

        intr = abs(T0 - t) < transit_half_duration  # create a mask of the in transit times
        oot = (abs(T0 - t) < oot_outer_relative) * (
            abs(T0 - t) < oot_inner_relative
        )  # create a mask of the out of transit times
        img_intr = tpf_filt[intr, :, :].sum(axis=0) / float(intr.sum())  # apply the masks and normalize the flux
        img_oot = tpf_filt[oot, :, :].sum(axis=0) / float(oot.sum())
        img_diff = img_oot - img_intr  # calculate the difference image (out of transit minus in-transit)

        # ---- PLOT -------

        # in transit
        count += 1  # add to the count before each plot
        plt.subplot(len(T0_list), 3, count)
        plt.axis("off")
        plt.imshow(img_intr, cmap=plt.cm.viridis, origin="lower")
        plt.colorbar()
        plt.title("t = {} days \n In Transit Flux (e-/candence)".format(T0), fontsize=9)

        # out of transit
        count += 1
        plt.subplot(len(T0_list), 3, count)
        plt.axis("off")
        plt.imshow(img_oot, cmap=plt.cm.viridis, origin="lower")
        plt.colorbar()
        plt.title("Out of Transit Flux (e-/candence)", fontsize=9)

        # out of transit minus in-transit
        count += 1
        plt.subplot(len(T0_list), 3, count)
        plt.axis("off")
        plt.imshow(img_diff, cmap=plt.cm.viridis, origin="lower")
        plt.colorbar()
        plt.title("Difference Flux (e-/candence)", fontsize=9)

    plt.subplots_adjust(wspace=0)
    plt.tight_layout()

    if not plot_lc:
        return None

    # ---- additional lightcurve plot to help visualization of the time span measured -------
    lc = tpf.to_lightcurve().remove_nans()
    lc = lc.truncate(T0 - oot_outer_relative * 1.25, T0 + oot_outer_relative * 1.25)
    with plt.style.context(lk.MPLSTYLE):
        ax = plt.figure(figsize=(6, 3)).gca()
        ax = lc.scatter(ax=ax)
        ax.axvline(T0, color="red", ymax=0.15, linewidth=1, linestyle="--", label="epoch")
        ax.axvspan(T0 - transit_half_duration, T0 + transit_half_duration, facecolor="red", alpha=0.3, label="In Transit")
        ax.axvspan(T0 - oot_outer_relative, T0 - oot_inner_relative, facecolor="green", alpha=0.3, label="Out of Transit")
        # no label to avoid double legend
        ax.axvspan(T0 + oot_inner_relative, T0 + oot_outer_relative, facecolor="green", alpha=0.3)
        ax.legend(loc="upper right", fontsize="small")


def plot_pixel_level_LC(tpf, epoch, transit_half_duration=0.25, oot_outer_relative=0.5, oot_inner_relative=0.3):
    """
    Plot the LC for each pixel around the time of the transit like event.
    Each LC is fitted with a 3 order polynomial in order to flatten.

    Returns
    -------
        Plot of the normalised LC for each pixel around the time of the transit like event.
        The pixel background colour represents the average flux.
        The time of the transit is highlighted in red/gold for each pixel LC.
    """

    # based on plot_pixel_level_LC() in
    # https://github.com/noraeisner/PH_Coffee_Chat/blob/4a723030ed80eeabdfb0eca49d948b97d61e35f6/False%20Positive/False%20positives%20-%20(3)%20pixel-level-lightcurve-plot.ipynb

    transit_list = [epoch]
    t_list = [tpf.time.value]
    tpf_list = [tpf.flux.value]
    bkg_list = [np.nanmean(tpf.flux.value, axis=0)]
    arrshape_list = [tpf.flux.shape]

    # loop through the transits and make plot for each ( only the first is currently displayed in the pdf report)
    plot_half_width = oot_outer_relative * 1.25
    for idx, X1_original in enumerate(tpf_list):

        bkg = np.flip(bkg_list[idx], axis=0)
        arrshape = arrshape_list[idx]
        peak = transit_list[idx]
        tpf = tpf_list[idx]

        s = X1_original.shape
        X1 = X1_original.reshape(s[0], s[1] * s[2])

        T0 = transit_list[idx]  # the time of the transit-like event
        t = t_list[idx]  # the time array

        intr = abs(T0 - t) < transit_half_duration  # create a mask of the in transit times
        oot = (abs(T0 - t) < oot_outer_relative) * (
            abs(T0 - t) < oot_inner_relative
        )  # create a mask of the out of transit times

        fig, ax = plt.subplots(
            arrshape[1], arrshape[2], sharex=True, sharey=False, gridspec_kw={"hspace": 0, "wspace": 0}, figsize=(5.5, 5.5)
        )

        plt.tight_layout()

        # see if the background of this plot can be the average pixel flux (if there are too many nans this will fail and the background will just be black which is also okay)
        try:
            color = plt.cm.viridis(np.linspace(0, 1, int(np.nanmax(bkg)) - int(np.nanmin(bkg)) + 1))
            simplebkg = False
        except:
            simplebkg = True

        for i in range(0, arrshape[1]):
            ii = arrshape[1] - 1 - i  # we want to plot this such that the pixels increase from left to right and bottom to top

            for j in range(0, arrshape[2]):

                apmask = np.zeros(arrshape[1:], dtype=np.int)
                apmask[i, j] = 1
                apmask = apmask.astype(bool)

                flux = X1[:, apmask.flatten()].sum(axis=1)

                m = np.nanmedian(flux[oot])

                normalizedflux = flux / m

                # bin the data
                f1 = normalizedflux
                time = t

                binfac = 7

                N = len(time)
                n = int(np.floor(N / binfac) * binfac)
                X = np.zeros((2, n))
                X[0, :] = time[:n]
                X[1, :] = f1[:n]
                Xb = rebin(X, (2, int(n / binfac)))

                # binned data
                time_binned = np.array(Xb[0])
                flux_binned = np.array(Xb[1])

                # create a mask that only looks at the times cut around the transit-event
                timemask = (time_binned < peak + plot_half_width) & (time_binned > peak - plot_half_width)

                time_binned = time_binned[timemask]
                flux_binned = flux_binned[timemask]

                # ----------
                # fit a spline to the cut-out of each pixel LC in order to flatten it
                p = np.poly1d(np.polyfit(time_binned, flux_binned, 3))
                flux_binned = flux_binned / p(time_binned)
                # ----------

                intr = abs(peak - time_binned) < 0.1

                if simplebkg == True:
                    ax[ii, j].set_facecolor(color="k")
                    linecolor = "w"
                    transitcolor = "gold"
                else:
                    ax[ii, j].set_facecolor(color=color[int(bkg[ii, j]) - int(np.nanmin(bkg))])

                    if int(bkg[ii, j]) - abs(int(np.nanmin(bkg))) > ((np.nanmax(bkg)) - abs(int(np.nanmin(bkg)))) / 2:
                        linecolor = "k"
                        transitcolor = "orangered"
                    else:
                        linecolor = "w"
                        transitcolor = "gold"

                ax[ii, j].plot(time_binned, flux_binned, color=linecolor, marker=".", markersize=1, lw=0)
                ax[ii, j].plot(time_binned[intr], flux_binned[intr], color=transitcolor, marker=".", markersize=1, lw=0)

                # get rid of ticks and ticklabels
                ax[ii, j].set_yticklabels([])
                ax[ii, j].set_xticklabels([])
                ax[ii, j].set_xticks([])
                ax[ii, j].set_yticks([])

        # ------------------

        # print("done.\n")
        # ------------------

        # label the pixels

        fig.text(0.5, 0.01, "column (pixel)", ha="center", fontsize=13)
        fig.text(0.01, 0.5, "row (pixel)", va="center", rotation="vertical", fontsize=13)

        # - - - - - - - - - -

        plt.subplots_adjust(top=0.95, right=0.99, bottom=0.04, left=0.04)

        plt.suptitle(r"T0 = {} $\pm$ {:.2f} d".format(peak, plot_half_width), y=0.98, fontsize=12)
        plt.xlim(peak - plot_half_width, peak + plot_half_width)
        plt.show()


def rebin(arr, new_shape):
    """
    function used to rebin the data
    """
    shape = (new_shape[0], arr.shape[0] // new_shape[0], new_shape[1], arr.shape[1] // new_shape[1])
    return arr.reshape(shape).mean(-1).mean(1)
