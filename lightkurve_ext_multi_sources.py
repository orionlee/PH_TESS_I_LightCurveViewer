#
# Lightkurve extension for use cases for combining data from multiple sources, e.g., for VSX submission.
#

from astropy.time import Time
from astropy import units as u
import numpy as np

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


def combine_tess_n_k2(lc_tess, lc_k2, shift_to_tess=True):
    lc_tess = lc_tess.copy()  # ensure users who further modify the result won't affect the source
    lc_k2 = lc_k2.copy()

    # OPEN: shift the data to a common band?!
    if shift_to_tess:
        shift_flux(lc_k2, lc_ref=lc_tess, inplace=True)

    return {"TESS": lc_tess, "K2": lc_k2}


def get_label_of_source(lc_dict, source):
    lc = lc_dict[source]
    mag_shift = lc.meta.get("FLUX_SHIFT", None)
    if mag_shift is not None and mag_shift != 0:
        sign_str = "+" if mag_shift > 0 else ""
        return f"{source} {sign_str}{mag_shift:.4f}"
    else:
        return source


def plot_tess_n_k2(lc_combined_dict, figsize, target_name):
    ax = tplt.lk_ax(figsize=figsize)
    ax.invert_yaxis()

    # scatter plot for the dense TESS data, error is relatively small
    lc = lc_combined_dict["TESS"]
    ax.scatter(lc.time.value, lc.flux.value, c="#3AF", s=0.1, alpha=1.0, label=get_label_of_source(lc_combined_dict, "TESS"))

    # scatter plot for K2 data
    lc = lc_combined_dict["K2"]
    ax.scatter(lc.time.value, lc.flux.value, c="green", s=0.5, alpha=1.0, label=get_label_of_source(lc_combined_dict, "K2"))

    ax.legend()
    ax.set_xlabel("Time [HJD]")
    ax.set_ylabel("Magnitude")
    ax.set_title(f"""{target_name}""")
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, p: format(int(x), ",")))  # thousands separator
    return ax


def fold_n_plot_tess_n_k2(
    lc_combined_dict, period, epoch: Time, phase_scale, target_coord=None, figsize=(8, 4), target_name=None, ax=None
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
        epoch_hjd = lke.to_hjd_utc(epoch, SkyCoord(target_coord["ra"], target_coord["dec"], unit=(u.deg, u.deg), frame="icrs"))

    lc_tess_f = fold_at_scale(lc_combined_dict.get("TESS"), epoch_time=epoch_hjd, period=period, normalize_phase=True)
    lc_k2_f = fold_at_scale(lc_combined_dict.get("K2"), epoch_time=epoch_hjd, period=period, normalize_phase=True)

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
