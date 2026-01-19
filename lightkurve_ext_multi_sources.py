#
# Lightkurve extension for use cases for combining data from multiple sources, e.g., for VSX submission.
#

from copy import deepcopy

from astropy.time import Time, TimeDelta
from astropy import units as u
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

import tic_plot as tplt
import lightkurve_ext as lke


def shift_flux(lc, lc_ref, inplace=False):
    shift = np.nanmedian(lc_ref.flux) - np.nanmedian(lc.flux)
    if not inplace:
        lc = lc.copy()
    lc.flux += shift
    lc.meta["FLUX_SHIFT"] = shift
    if inplace:
        return
    else:
        return lc


def combine_multi_bands_and_shift(lc_dict, split_lc_to_multi_bands=True, shift_to=None):
    def split_some_lcs_to_multi_bands(lc_dict):
        """Helper to split single ASAS-SN (and possibly other) lc to multiple ones by bands"""
        new_dict = {}
        for band, lc in lc_dict.items():
            if band == "ASAS-SN":
                for filter in np.unique(lc.filter):
                    lc_of_filter = lc[lc.filter == filter]
                    new_dict[f"ASAS-SN {filter}"] = lc_of_filter
            else:
                new_dict[band] = lc
        return new_dict

    if split_lc_to_multi_bands:
        lc_dict = split_some_lcs_to_multi_bands(lc_dict)

    res = {}
    lc_ref = lc_dict[shift_to] if shift_to is not None else None
    for band, lc in lc_dict.items():
        lc = (
            lc.copy()
        )  # ensure users who further modify the result won't affect the source
        if lc_ref is not None and band != shift_to:
            shift_flux(lc, lc_ref=lc_ref, inplace=True)
        lc.meta["BAND"] = band
        res[band] = lc
    return res


def combine_tess_n_k2(lc_tess, lc_k2, shift_to_tess=True):
    lc_tess = (
        lc_tess.copy()
    )  # ensure users who further modify the result won't affect the source
    lc_k2 = lc_k2.copy()

    # OPEN: shift the data to a common band?!
    if shift_to_tess:
        shift_flux(lc_k2, lc_ref=lc_tess, inplace=True)

    return {"TESS": lc_tess, "K2": lc_k2}


def get_label_of_source(lc_dict, source, mag_shift_precision=3):
    def get_rounded_str(val, precision):
        # Can't use the naive version: `f"{np.round(val, precision)}"`
        #
        # for some reason, in a specific case
        # (in my notebook TIC_75722876_EA.ipynb),
        # in rounding the value 1.1312150955200195 mag (TESS mag shift, in MaskedQuantity),
        # the naive version returns 1.1299999952316284 mag
        # I cannot reproduce the problem in isolated codes.
        #
        # This function is a workaround, which requires:
        # 1. strip the unit (handle it separately), and
        # 2. use str() on the rounded value (f-string would produce the incorrect value)

        if isinstance(val, u.Quantity):
            val_unit = f" {val.unit}"  # space prepended
            val = val.value
        else:
            val_unit = ""
        val_rounded = np.round(val, precision)
        return str(val_rounded), val_unit  # must use str(), cannot use f-string

    lc = lc_dict[source]
    mag_shift = lc.meta.get("FLUX_SHIFT", None)
    if mag_shift is not None and mag_shift != 0:
        sign_str = "+" if mag_shift > 0 else ""
        mag_shift_rounded, mag_unit = get_rounded_str(mag_shift, mag_shift_precision)
        if mag_shift_rounded != "0":
            return f"{source} {sign_str}{mag_shift_rounded}{mag_unit}"
        else:
            return source
    else:
        return source


# the default plot options support the pattern that
# the first lc is from high cadences plot such as TESS / Kepler,
# and the remaining LCs are from ground-based observations,
# with relatively sparse data and larger errors.
DEFAULT_MULTI_BANDS_PLOT_OPTIONS = [
    (
        "scatter",
        dict(
            c="#3AF",
            s=0.1,
            alpha=1.0,
        ),
    ),
    (
        "errorbar",
        dict(
            marker=".",
            c="green",
            linewidth=0.5,
            ls="none",
        ),
    ),
    (
        "errorbar",
        dict(
            marker=".",
            c="pink",
            linewidth=0.5,
            ls="none",
        ),
    ),
    (
        "errorbar",
        dict(
            marker=".",
            c="violet",
            linewidth=0.5,
            ls="none",
        ),
    ),
    (
        "errorbar",
        dict(
            marker=".",
            c="orange",
            linewidth=0.5,
            ls="none",
        ),
    ),
]


def get_default_plot_multi_bands_options_copy():
    """Return a copy of the default plot_multi_bands(), so that
    callers could customize it for specific usage.
    """
    return deepcopy(DEFAULT_MULTI_BANDS_PLOT_OPTIONS)


def _flip_yaxis_for_mag(ax, lc, plot_kwargs):
    y_column = plot_kwargs.get("column", "flux")
    # invert y-axis only when it hasn't been inverted
    # to support multiple scatter/plot/errorbar calls on the same ax object
    if lc[y_column].unit == u.mag and ax.get_ylim()[1] > ax.get_ylim()[0]:
        ax.invert_yaxis()
    return ax


def plot_multi_bands(
    lc_combined_dict,
    figsize,
    target_name,
    ax=None,
    plot_options=None,
    mag_shift_precision=2,
    # parameters used by fold_n_plot_multi_bands() use cases
    include_labels=True,
    time_shift_func=None,
    show_colorbar_if_applicable=False,  # <-- needed to avoid showing colorbar twice in doing a 2X folded plot
):
    if ax is None:
        ax = tplt.lk_ax(figsize=figsize)

    if plot_options is None:
        plot_options = DEFAULT_MULTI_BANDS_PLOT_OPTIONS
    plot_options = deepcopy(plot_options)  # avoid modifying the original copy

    for i, band in enumerate(lc_combined_dict):
        lc = lc_combined_dict[band]
        plot_funcname, plot_kwargs = plot_options[i]
        show_colorbar = False
        if plot_funcname == "errorbar":
            plot_kwargs["yerr"] = lc.flux_err.value
        plot_func = getattr(ax, plot_funcname)
        if plot_kwargs.get("c") == "_time_":
            plot_kwargs["c"] = lc.time_original.value
            if show_colorbar_if_applicable:
                show_colorbar = True
        if include_labels:
            plot_label = get_label_of_source(
                lc_combined_dict, band, mag_shift_precision
            )
        else:
            plot_label = None
        x_vals = lc.time.value
        if time_shift_func is not None:
            x_vals = time_shift_func(x_vals)
        path_collection = plot_func(
            x_vals, lc.flux.value, label=plot_label, **plot_kwargs
        )
        if show_colorbar:
            cbar = plt.colorbar(path_collection, ax=ax)
            cbar.set_label("Time")
        _flip_yaxis_for_mag(ax, lc, plot_kwargs)
    ax.legend()

    if isinstance(lc.time, Time):
        xlabel = f"Time [{lc.time.format.upper()}]"
    elif isinstance(lc.time, TimeDelta):
        xlabel = f"Phase [{lc.time.format.upper()}]"
    else:  # normalized phase
        xlabel = "Phase"
    ax.set_xlabel(xlabel)

    ax.set_ylabel("Magnitude")
    ax.set_title(f"""{target_name}""")

    if lc.time.min().value > 10000:  # probably JD / MJD
        ax.xaxis.set_major_formatter(
            FuncFormatter(lambda x, p: format(int(x), ","))
        )  # thousands separator

    return ax


def plot_tess_n_ztf(lc_combined_dict, figsize, target_name, mag_shift_precision=3):
    ax = tplt.lk_ax(figsize=figsize)
    ax.invert_yaxis()

    # scatter plot for the dense TESS data, error is relatively small
    lc = lc_combined_dict["TESS"]
    ax.scatter(
        lc.time.value,
        lc.flux.value,
        c="#3AF",
        s=0.1,
        alpha=1.0,
        label=get_label_of_source(lc_combined_dict, "TESS", mag_shift_precision),
    )

    # scatter plot for ZTF data
    # TODO: do not assume ZTF data must be in ZTF g
    lc = lc_combined_dict["ZTF g"]
    # ax.scatter(lc.time.value, lc.flux.value, c="green", s=1.0, alpha=1.0, label=get_label_of_source(lc_combined_dict, "ZTF g", mag_shift_precision))
    ax.errorbar(
        x=lc.time.value,
        y=lc.flux.value,
        yerr=lc.flux_err.value,
        marker=".",
        c="green",
        linewidth=0.5,
        ls="none",
        label=get_label_of_source(lc_combined_dict, "ZTF g", mag_shift_precision),
    )

    ax.legend()
    ax.set_xlabel("Time [HJD]")
    ax.set_ylabel("Magnitude")
    ax.set_title(f"""{target_name}""")
    ax.xaxis.set_major_formatter(
        FuncFormatter(lambda x, p: format(int(x), ","))
    )  # thousands separator
    return ax


def fold_n_plot_multi_bands(
    lc_combined_dict,
    period,
    epoch: Time,
    phase_scale=2,
    target_name=None,
    duration_hr=None,  # used for plotting purpose only
    # for plotting only, the midpoint of which duration_hr is based on
    # typically to draw lines for secondary eclipses / transits
    duration_midpoint_phase=0,
    mag_shift_precision=2,
    figsize=(8, 4),
    ax=None,
    plot_options=None,
):
    if phase_scale not in [1, 2]:
        raise ValueError("phase_scale must be 1 (plotted once) or 2 (plotted twice).")

    lc_f_combined_dict = {}
    for band, lc in lc_combined_dict.items():
        # Note: the folded lc_f is the same regardless of phase_scale, which would
        # only affect the plot result below.
        lc_f = lc.fold(epoch_time=epoch, period=period, normalize_phase=True)
        lc_f_combined_dict[band] = lc_f

    ax = plot_multi_bands(
        lc_f_combined_dict,
        figsize=figsize,
        target_name=target_name,
        ax=ax,
        plot_options=plot_options,
        mag_shift_precision=mag_shift_precision,
        show_colorbar_if_applicable=True,  # show colorbar only in 1st plot for case phase_scale == 2
    )
    if phase_scale == 2:
        # phase [-0.5, +0.5] has been plotted above
        # now plot phases [0.5, 1.0] and [-1.0, -0.5]
        def time_shit_func(x):
            # for x in [-0.5, 0.0], shift to  phase [0.5, 1.0]
            res1 = x + 1
            res1[x > 0] = 0

            # for x in [0.0, +0.5], shift to phase [-1.0, -0.5]
            res2 = x - 1
            res2[x <= 0] = 0

            return res1 + res2

        ax = plot_multi_bands(
            lc_f_combined_dict,
            figsize=figsize,
            target_name=target_name,
            ax=ax,
            plot_options=plot_options,
            mag_shift_precision=mag_shift_precision,
            include_labels=False,
            time_shift_func=time_shit_func,
        )

    if duration_hr is not None:
        duration_phase = duration_hr / 24 / period
        ax.axvline(
            duration_midpoint_phase - duration_phase / 2, linestyle="--", c="blue"
        )
        ax.axvline(
            duration_midpoint_phase + duration_phase / 2, linestyle="--", c="blue"
        )

    # Set phase plot specific title
    time_all = np.array([])
    for lc in lc_combined_dict.values():
        time_all = np.concatenate([time_all, lc.time.to_value("jd")])
    plot_time_span = time_all.max() - time_all.min()

    title = f"""{target_name}
    period: {period} d"""
    if period < 1 / 24:
        period_min = round(period * 24 * 60, 3)
        title += f" ({period_min} m)"
    title += f", epoch={epoch.format.upper()} {epoch.value}, time span: {plot_time_span:.0f}d"
    ax.set_title(title)

    return ax, lc_f_combined_dict


def fold_n_plot_tess_n_ztf(
    lc_combined_dict,
    period,
    epoch: Time,
    phase_scale,
    target_coord=None,
    figsize=(8, 4),
    target_name=None,
    mag_shift_precision=3,
    ax=None,
):
    def fold_at_scale(lc, **kwargs):
        if lc is None:
            return None
        kwargs = kwargs.copy()
        kwargs["period"] = kwargs["period"] * phase_scale
        lc_f = lc.fold(**kwargs)
        return lc_f

    from astropy.coordinates import SkyCoord

    if epoch.format == "jd" and epoch.scale == "utc":
        # if already in HJD UTC, avoid unnecessary conversion which would also have undesirable effect on precision
        epoch_hjd = epoch
    else:
        epoch_hjd = lke.to_hjd_utc(
            epoch,
            SkyCoord(
                target_coord["ra"],
                target_coord["dec"],
                unit=(u.deg, u.deg),
                frame="icrs",
            ),
        )

    lc_tess_f = fold_at_scale(
        lc_combined_dict.get("TESS"),
        epoch_time=epoch_hjd,
        period=period,
        normalize_phase=True,
    )
    # TODO: do not assume ZTF data must be in ZTF g
    lc_ztf_f = fold_at_scale(
        lc_combined_dict.get("ZTF g"),
        epoch_time=epoch_hjd,
        period=period,
        normalize_phase=True,
    )

    if ax is None:
        ax = tplt.lk_ax(figsize=figsize)

    ax.invert_yaxis()

    if lc_tess_f is not None and len(lc_tess_f) > 0:
        ax.scatter(
            lc_tess_f.time * phase_scale,
            lc_tess_f.flux.value,
            c="#3AF",
            s=0.1,
            alpha=1.0,
            label=get_label_of_source(lc_combined_dict, "TESS", mag_shift_precision),
        )

    if lc_ztf_f is not None and len(lc_ztf_f) > 0:
        ax.errorbar(
            x=lc_ztf_f.time * phase_scale,
            y=lc_ztf_f.flux.value,
            yerr=lc_ztf_f.flux_err.value,
            marker=".",
            c="green",
            linewidth=0.5,
            ls="none",
            label=get_label_of_source(lc_combined_dict, "ZTF g", mag_shift_precision),
        )

    ax.legend()
    ax.set_xlabel("Phase")
    ax.set_ylabel("Magnitude")
    time_all = np.array([])
    if lc_tess_f is not None:
        time_all = np.concatenate([time_all, lc_tess_f.time_original.to_value("jd")])
    if lc_ztf_f is not None:
        time_all = np.concatenate([time_all, lc_ztf_f.time_original.to_value("jd")])
    plot_time_span = time_all.max() - time_all.min()
    ax.set_title(
        f"""{target_name}
    period: {period}d, epoch={epoch_hjd.value}, time span: {plot_time_span:.0f}d
    """
    )

    return ax, {"TESS": lc_tess_f, "ZTF g": lc_ztf_f}


def plot_tess_n_k2(lc_combined_dict, figsize, target_name):
    ax = tplt.lk_ax(figsize=figsize)
    ax.invert_yaxis()

    # scatter plot for the dense TESS data, error is relatively small
    lc = lc_combined_dict["TESS"]
    ax.scatter(
        lc.time.value,
        lc.flux.value,
        c="#3AF",
        s=0.1,
        alpha=1.0,
        label=get_label_of_source(lc_combined_dict, "TESS"),
    )

    # scatter plot for K2 data
    lc = lc_combined_dict["K2"]
    ax.scatter(
        lc.time.value,
        lc.flux.value,
        c="green",
        s=0.5,
        alpha=1.0,
        label=get_label_of_source(lc_combined_dict, "K2"),
    )

    ax.legend()
    ax.set_xlabel("Time [HJD]")
    ax.set_ylabel("Magnitude")
    ax.set_title(f"""{target_name}""")
    ax.xaxis.set_major_formatter(
        FuncFormatter(lambda x, p: format(int(x), ","))
    )  # thousands separator
    return ax


def fold_n_plot_tess_n_k2(
    lc_combined_dict,
    period,
    epoch: Time,
    phase_scale,
    target_coord=None,
    figsize=(8, 4),
    target_name=None,
    ax=None,
):
    def fold_at_scale(lc, **kwargs):
        if lc is None:
            return None
        kwargs = kwargs.copy()
        kwargs["period"] = kwargs["period"] * phase_scale
        lc_f = lc.fold(**kwargs)
        return lc_f

    from astropy.coordinates import SkyCoord

    if epoch.format == "jd" and epoch.scale == "utc":
        # if already in HJD UTC, avoid unnecessary conversion which would also have undesirable effect on precision
        epoch_hjd = epoch
    else:
        epoch_hjd = lke.to_hjd_utc(
            epoch,
            SkyCoord(
                target_coord["ra"],
                target_coord["dec"],
                unit=(u.deg, u.deg),
                frame="icrs",
            ),
        )

    lc_tess_f = fold_at_scale(
        lc_combined_dict.get("TESS"),
        epoch_time=epoch_hjd,
        period=period,
        normalize_phase=True,
    )
    lc_k2_f = fold_at_scale(
        lc_combined_dict.get("K2"),
        epoch_time=epoch_hjd,
        period=period,
        normalize_phase=True,
    )

    if ax is None:
        ax = tplt.lk_ax(figsize=figsize)

    ax.invert_yaxis()

    if lc_tess_f is not None and len(lc_tess_f) > 0:
        ax.scatter(
            lc_tess_f.time * phase_scale,
            lc_tess_f.flux.value,
            c="#3AF",
            s=0.1,
            alpha=1.0,
            label=get_label_of_source(lc_combined_dict, "TESS"),
        )

    if lc_k2_f is not None and len(lc_k2_f) > 0:
        ax.scatter(
            lc_k2_f.time * phase_scale,
            lc_k2_f.flux.value,
            c="green",
            s=0.5,
            alpha=1.0,
            label=get_label_of_source(lc_combined_dict, "K2"),
        )

    ax.legend()
    ax.set_xlabel("Phase")
    ax.set_ylabel("Magnitude")
    time_all = np.array([])
    if lc_tess_f is not None:
        time_all = np.concatenate([time_all, lc_tess_f.time_original.to_value("jd")])
    if lc_k2_f is not None:
        time_all = np.concatenate([time_all, lc_k2_f.time_original.to_value("jd")])
    plot_time_span = time_all.max() - time_all.min()
    ax.set_title(
        f"""{target_name}
    period: {period}d, epoch={epoch_hjd.value}, time span: {plot_time_span:.0f}d
    """
    )

    return ax, {"TESS": lc_tess_f, "K2": lc_k2_f}
