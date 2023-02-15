# Various helpers to work with Periodogram
#
import logging
import warnings

from types import SimpleNamespace
from memoization import cached

from astropy.time import Time
from astropy import units as u
import numpy as np

import lightkurve as lk
from lightkurve.periodogram import BoxLeastSquaresPeriodogram

from IPython.display import display, HTML

# for type annotation
from numbers import Number

log = logging.getLogger(__name__)


def plot_pg_n_mark_max(pg, ax=None, max_period_factor=None):
    ax = pg.plot(ax=ax, view="period")
    ax.axvline(pg.period_at_max_power.value, c="blue", alpha=0.4)
    period_val_text = f"{pg.period_at_max_power:.4f}"
    if pg.period_at_max_power < 0.2 * u.day and pg.period_at_max_power.unit == u.day:
        # for short period in days, also show the period in hours
        period_val_text += f" ({pg.period_at_max_power.to(u.hour):.4f})"
    max_text = f"Power: {pg.max_power:.2f}, Period: {period_val_text}"
    if hasattr(pg, "depth_at_max_power"):  # to accommodate LombScarglePeriodogram
        max_text += f", Depth: {pg.depth_at_max_power:.6f}"

    if max_period_factor is not None and ax.get_xlim()[1] > pg.period_at_max_power.value * max_period_factor:
        ax.set_xlim(0 if ax.get_xlim()[0] < 0 else None, pg.period_at_max_power.value * max_period_factor)
        display(HTML("""<span style="background-color: yellow;">Note:</span>Long potential periods truncated from the plot"""))

    x, y = pg.period_at_max_power.value, pg.max_power.value * 0.9
    x_min, x_max = ax.get_xlim()
    x_mid = (x_max - x_min) / 2 + x_min
    horizontalalignment = "left" if x < x_mid else "right"
    ax.text(x, y, " " + max_text + " ", c="blue", horizontalalignment=horizontalalignment)
    return ax


def model(pg, lc, **kwargs):
    if hasattr(pg, "get_transit_model"):  # case BLS pg
        return _bls_model(pg, lc, **kwargs)
    else:  # case LS pg
        return pg.model(lc.time, pg.frequency_at_max_power)


def _bls_model(pg, lc, time=None, period=None, duration=None, transit_time=None):
    # BLS pg.get_transit_model() does not support user-supplied time
    if time is None:
        time = lc.time  # defaulted to lc.time rather than pg.time
    if period is None:
        period = pg.period_at_max_power
        # log.warning("No period specified. Using period at max power")
    if duration is None:
        duration = pg.duration_at_max_power
        # log.warning("No duration specified. Using duration at max power")
    if transit_time is None:
        transit_time = pg.transit_time_at_max_power
        # log.warning("No transit time specified. Using transit time at max power")
    if not isinstance(transit_time, Time):
        transit_time = Time(transit_time, format=pg.time.format, scale=pg.time.scale)

    model_flux = pg._BLS_object.model(
        time,
        u.Quantity(period, "d").value,
        u.Quantity(duration, "d").value,
        transit_time,
    )
    model = lk.LightCurve(time=time, flux=model_flux, label="Transit Model Flux")
    return model


def plot_lc_with_model(lc, pg, plot_lc=True, plot_model=True, plot_folded_model=True, also_return_lcs=False):
    lc_model = model(pg, lc)
    ax1 = None
    if plot_lc:
        ax1 = lc.scatter()
        if hasattr(pg, "get_transit_mask"):
            lc[pg.get_transit_mask()].scatter(ax=ax1, c="orange", marker="x", s=9, label="in transits")

    ax2 = None
    if plot_model:
        ax2 = lc.scatter(alpha=0.5)
        lc_model.plot(ax=ax2, c="red", linewidth=2)

    # folded, zoom -in
    ax_f = None
    if plot_folded_model:
        period = pg.period_at_max_power
        if hasattr(pg, "transit_time_at_max_power"):
            epoch_time = pg.transit_time_at_max_power
        else:
            epoch_time = None
        lc_f = lc.fold(epoch_time=epoch_time, period=period)
        lc_model_f = lc_model.fold(epoch_time=epoch_time, period=period)

        ax_f = lc_f.scatter(alpha=0.3)
        lc_model_f.plot(ax=ax_f, c="red", linewidth=4)
        if hasattr(pg, "duration_at_max_power"):
            # zoom in for BLS model:
            ax_f.set_xlim(-pg.duration_at_max_power.value, pg.duration_at_max_power.value)

    if not also_return_lcs:
        return ax1, ax2, ax_f
    else:
        lcs = SimpleNamespace(lc=lc, lc_f=lc_f, lc_model=lc_model, lc_model_f=lc_model_f)
        return [ax1, ax2, ax_f], lcs


def errorbar_transit_depth(pg):
    lc_model = pg.get_transit_model()
    ax = lc_model.plot(c="red", alpha=0.2)

    time = [lc_model.time.min().value, lc_model.time.max().value]
    ax.plot(
        time, [1 - pg.depth_at_max_power, 1 - pg.depth_at_max_power], linestyle="-.", c="red", alpha=0.5, label="Model depth"
    )

    transit_depth_mean = np.nanmean(
        pg.transit_depth
    )  # note pg._TLS_result.depth_mean[0] would be more accurate, it also has the error
    ax.plot(
        time,
        [1 - transit_depth_mean, 1 - transit_depth_mean],
        linestyle="--",
        c="black",
        alpha=0.8,
        label="Actual transit depth mean",
    )

    ax.errorbar(
        pg.transit_time.value,
        1 - pg.transit_depth,
        pg.transit_depth_err,
        fmt="o",
        color="black",
        label="Actual transit mid points",
    )

    ax.legend()
    return ax


def find_peaks(pg, powerlimit=None):
    """Lists in descending order the peaks found in the periodogram with
    the find_peaks function, which calculates the Prominence of the peak,
    and Lower & Upper Half Width at Half Maximum (HWHM) of this Prominence.
    For more information on this SciPy signal function,
    please follow the link to their documentation:
    https://docs.scipy.org/doc//scipy/reference/generated/scipy.signal.find_peaks.html
    Additionally final column of the list shows the ratio of the peak with
    the top peak, where deviation from neat fractions (like 2/1, 1/2, 3/2,
    2/3, 4/3, 3/4, etc.) is an indication of another candidate.
    Parameters
    ----------
    powerlimit : number or ndarray or sequence, optional
        Required power of peaks. Either a number, None, an array matching x
        or a 2-element sequence of the former. The first element is always
        interpreted as the minimal power and the second, if supplied, as the
        maximal required power. By default or None, 5% of the periodogram's
        max_power value will be used as the number for the powerlimit.
    Returns
    -------
    Table : `astropy table` object
        Returns a Table object extracted from the periodogram.
    """
    # based on https://github.com/lightkurve/lightkurve/pull/1255/
    # - added FWHM
    # - rename the label x from Periodicity to Period
    from scipy.signal import find_peaks as scipy_find_peaks
    from astropy.table import Table
    from astropy import units as u

    if pg.default_view == "period":
        view = pg.period
        x = "period"
        y = ".1f"
    elif pg.default_view == "frequency":
        view = pg.frequency
        x = "frequency"
        y = None
    if powerlimit is None:
        powerlimit = pg.max_power / 20
        if hasattr(powerlimit, "value"):  # powerlimit must be unit-less for scipy.find_peaks()
            powerlimit = powerlimit.value
    peaks, stats = scipy_find_peaks(pg.power, height=powerlimit, width=1)
    lhwhm_int_down = view[np.floor(stats["left_ips"]).astype(int)]
    lhwhm_int_up = view[np.ceil(stats["left_ips"]).astype(int)]
    lhwhm_int_remainder = stats["left_ips"] - np.floor(stats["left_ips"])
    lhwhm_period = lhwhm_int_down + lhwhm_int_remainder * (lhwhm_int_up - lhwhm_int_down)
    lhwhm = lhwhm_period - view[peaks]
    uhwhm_int_down = view[np.floor(stats["right_ips"]).astype(int)]
    uhwhm_int_up = view[np.ceil(stats["right_ips"]).astype(int)]
    uhwhm_int_remainder = stats["right_ips"] - np.floor(stats["right_ips"])
    uhwhm_period = uhwhm_int_down + uhwhm_int_remainder * (uhwhm_int_up - uhwhm_int_down)
    uhwhm = uhwhm_period - view[peaks]
    fwhm = np.abs(uhwhm) + np.abs(lhwhm)  # make it immune to the sign of lhwhm, uhwhm
    result = Table(
        data=[stats["peak_heights"], view[peaks], stats["prominences"], lhwhm, uhwhm, fwhm],
        names=("power", x, "prominence", "lower_hwhm", "upper_hwhm", "fwhm"),
    )
    result.sort("prominence", reverse=True)
    result[x + "_ratio"] = result[x][0] / result[x]
    result["power"].format = y
    result["power"].unit = pg.power.unit
    result["prominence"].format = y
    result["lower_hwhm"].format = ".5g"
    result["upper_hwhm"].format = "+.5g"
    result["fwhm"].format = ".5g"
    result[x + "_ratio"].format = ".3f"
    result[x + "_ratio"].unit = u.dimensionless_unscaled
    return result


@cached
def idx_bin_peaks_by_half_power(x, y):
    """Bin a power spectrum `(x, y)`, such that for each peak,
    the surrounding values which  are more than half the peak are removed.
    """
    idx_kept = np.ones_like(y, dtype=int) * -1
    num_kept = 0

    idx_sorted = np.flip(np.argsort(y))

    y_working_copy = np.array(y, dtype=float)
    for i in idx_sorted:
        cur_peak = y_working_copy[i]
        #         print('X1', i, cur_peak, ' -- ', y_working_copy)
        if np.isnan(cur_peak):
            continue

        #         print('X1a keeping ', i, cur_peak)

        # case keep the peak
        idx_kept[num_kept] = i
        num_kept += 1
        y_working_copy[i] = np.nan  # mark is as counted

        # now remove adjacent signals, if the value (y) is more than half of current peak
        threshold_for_removal = cur_peak / 2
        j = i + 1
        while True:
            if j >= len(y_working_copy) or j < 0:
                break
            if np.isnan(y_working_copy[j]) or y_working_copy[j] < threshold_for_removal:
                #                 print('X2', cur_peak, y_working_copy[j],  threshold_for_removal)
                break
            y_working_copy[j] = np.nan
            j = j + 1
        j = i - 1
        while True:
            if j >= len(y_working_copy) or j < 0:
                break
            if np.isnan(y_working_copy[j]) or y_working_copy[j] < threshold_for_removal:
                break
            y_working_copy[j] = np.nan
            j = j - 1

    #     print('X3', num_kept, idx_kept)
    idx_kept = idx_kept[:num_kept]

    # idx_kept is used to select x / y, without changing the order, so I need to sort it again
    return np.sort(idx_kept)


def plot_pg_n_top_frequencies(pg, num_top, freq_filter_func=lambda f: True):
    """Plot a periodogram and its top frequencies"""
    # condense values surrounding peaks to peaks themselves,
    # so that the top frequencies returned would be more meaningful
    # (they won't be just those surrounding the global peak)
    idx_peaks = idx_bin_peaks_by_half_power(pg.frequency.value, pg.power.value)
    # OPEN: the peaks identified are still too frequent in many cases
    # consider to further reduce them, e.g., consolidating a peak with those within 3 sigma
    # with half width to half power treated as 1 sigma

    frequency_b, power_b = pg.frequency[idx_peaks], pg.power[idx_peaks]

    ax = pg.plot(view="frequency")
    #     ax = lk_ax()
    #     ax.scatter(frequency_b.value, power_b.value, marker=".", s=4)

    # plot top pulsating frequency (set a cutoff point to be 450 upon the periodogram)
    idx_sorted = np.flip(np.argsort(power_b))
    top_frequencies = [f.value for f in frequency_b[idx_sorted] if freq_filter_func(f)][:num_top]
    ymin, ymax = ax.get_ylim()
    ax.vlines(top_frequencies, ymin=ymin, ymax=ymax, linestyle="dashed", alpha=0.5, linewidth=0.5)

    return ax, top_frequencies


def sde(pg: BoxLeastSquaresPeriodogram) -> Number:
    return (pg.max_power - np.mean(pg.power)) / np.std(pg.power)


def validate_bls_n_report(pg, to_display=True):
    # See https://docs.astropy.org/en/stable/timeseries/bls.html
    # and
    # https://docs.astropy.org/en/stable/api/astropy.timeseries.BoxLeastSquares.html#astropy.timeseries.BoxLeastSquares.compute_stats

    def metrics(name, value, flag_func=None):
        to_flag = flag_func(value) if flag_func is not None else False
        return SimpleNamespace(name=name, value=value, flag=to_flag)

    def output_metrics(metrics_list):
        def html(m):
            flag_style = "color: red; font-weight: bold" if m.flag else ""
            name = m.name.replace("\n", "<br>")
            return f"""<tr><td>{name}</td><td style="{flag_style}">{m.value}</td></tr>"""

        html_body = "\n".join([html(m) for m in metrics_list])
        return f"""
<table>
<tr><th>Metrics</th><th>Value</th</tr>
{html_body}
</table>
        """

    # SDE, Signal Detection Efficiency, eq 6 of the paper
    # assuming the power of the BLS is proportional to SR (Signal Residue, defined in eq 5, which is a form of chi-square statistics)
    # chi-square stats: stats based on sum of squares of random variable
    # in this context, the random variable is the difference between the model flux and the observed flux
    # it seems to be the case, based on https://docs.astropy.org/en/stable/timeseries/bls.html
    # - model flux are the y_in and y_out in the equation (dervied from observation, for a given period + epoch + duration)
    # - observed flux is the y_n in the equation
    #  (the power and SR might differ by some constant, but for SDE calculation constant does not matter)
    sde = (pg.max_power - np.mean(pg.power)) / np.std(pg.power)

    # From TLS paper, it seems that the SDE to FAP mapping can be used for BLS results too
    # (it's based on empirical fits on synthetic data with white noises)
    fap = None
    try:
        from transitleastsquares.stats import FAP as fap_from_tls

        fap = fap_from_tls(sde)
    except ImportError:
        warnings.warn(
            "validate_bls_n_report() cannot calculate FAP, because the depdent package transitleastsquares is not installed"
        )

    stats = pg.compute_stats()

    # adapted from
    # https://github.com/hippke/tls/blob/1cb4e599812a181ea5b95ee5386e8966fd47577b/build/lib/transitleastsquares/stats.py#L80
    def period_uncertainty(periods, power):
        # Determine estimate for uncertainty in period
        # Method: Full width at half maximum
        try:
            # Upper limit
            index_highest_power = np.argmax(power)
            idx = index_highest_power
            while True:
                idx += 1
                if power[idx] <= 0.5 * power[index_highest_power]:
                    idx_upper = idx
                    break
            # Lower limit
            idx = index_highest_power
            while True:
                idx -= 1
                if power[idx] <= 0.5 * power[index_highest_power]:
                    idx_lower = idx
                    break
            period_uncertainty = 0.5 * (periods[idx_upper] - periods[idx_lower])
        except:
            period_uncertainty = float("inf")
        return period_uncertainty

    period_at_max_power_err = period_uncertainty(pg.period.value, pg.power.value)

    def calc_empty_transit_count(stats):
        count = 0
        for ptc in stats["per_transit_count"]:
            if ptc < 1:
                count += 1
        return count

    empty_transit_count = calc_empty_transit_count(stats)

    def calc_odd_even_mismatch(stats):
        depth_mean_odd, depth_mean_odd_std = stats["depth_odd"]
        depth_mean_even, depth_mean_even_std = stats["depth_even"]

        # Odd even mismatch in standard deviations
        odd_even_difference = abs(depth_mean_odd - depth_mean_even)
        odd_even_std_sum = depth_mean_odd_std + depth_mean_even_std
        odd_even_mismatch = odd_even_difference / odd_even_std_sum
        return odd_even_mismatch

    odd_even_mismatch = calc_odd_even_mismatch(stats)

    # Use False Alarm Probability (on white noise only run)
    # threshold at 0.1%, or SDE 8.3
    # OPEN: consider to use the more stringent threshold 0.01% / SDE 9.1
    # (the value computed is optimistic as it's based on white noise data, actual data has some correlation varation (poink noise))
    # See:
    # https://transitleastsquares.readthedocs.io/en/latest/FAQ.html#false-alarm-probability

    def flag_fap(false_alarm_probability):
        # None means it is not available because the dependent package is not in the environment
        return false_alarm_probability is not None and (np.isnan(false_alarm_probability) or false_alarm_probability > 0.001)

    def flag_sde(sde):
        return sde < 8.3

    def flag_snr(snr):
        return snr < 10  # arbitrary for now

    def flag_odd_even(odd_even_mismatch):
        return odd_even_mismatch > 3  # i.e, > 3 sigma

    def flag_empty_transit_count(empty_transit_count):
        return empty_transit_count > 0

    def flag_harmonic_delta_log_likelihood(harmonic_delta_log_likelihood):
        return harmonic_delta_log_likelihood > 0

    metrics_list = [
        metrics("Period", pg.period_at_max_power),
        metrics("- Period error\n(half width at half max)", period_at_max_power_err),
        metrics(f"Epoch ({pg.transit_time_at_max_power.format.upper()})", pg.transit_time_at_max_power),
        metrics("Duration", pg.duration_at_max_power),
        metrics("Depth (Model)", pg.depth_at_max_power),
        metrics("Depth (Mean)", stats["depth"][0]),  # mean of all in-transit flux
        metrics("- Depth (Mean) error", stats["depth"][1]),
        metrics(
            "Period grid size", len(pg.period)
        ),  # larger period grid enerally results in more accurate period and better SDE
        metrics("FAP (white noises)", fap, flag_fap),
        metrics("SDE", sde, flag_sde),
        metrics("Power (Log likelihood)", pg.max_power),
        metrics("SNR", np.max(pg.snr), flag_snr),
        metrics("SNR per transit", "N/A"),
        metrics("SNR (pink)\nconerns if greatly different from SNR", "N/A"),
        metrics("Empty Transits\nactual period might be double", empty_transit_count, flag_empty_transit_count),  # todo
        metrics("Odd-even mismatch (sigma)", odd_even_mismatch, flag_odd_even),
        metrics(
            "Sine-like\n+ve => sine model better fit",
            stats["harmonic_delta_log_likelihood"],
            flag_harmonic_delta_log_likelihood,
        ),
        metrics("Elapsed time (ms)", getattr(pg, "elapsed_time", "N/A")),
    ]

    html = output_metrics(metrics_list)
    if to_display:
        return display(HTML(html))
    else:
        return html


def validate_tls_n_report(pg, to_display=True):
    # See https://transitleastsquares.readthedocs.io/en/latest/Python%20interface.html#return-values
    # for TLS return values

    def metrics(name, value, flag_func=None):
        to_flag = flag_func(value) if flag_func is not None else False
        return SimpleNamespace(name=name, value=value, flag=to_flag)

    def output_metrics(metrics_list):
        def html(m):
            flag_style = "color: red; font-weight: bold" if m.flag else ""
            name = m.name.replace("\n", "<br>")
            return f"""<tr><td>{name}</td><td style="{flag_style}">{m.value}</td></tr>"""

        html_body = "\n".join([html(m) for m in metrics_list])
        return f"""
<table>
<tr><th>Metrics</th><th>Value</th</tr>
{html_body}
</table>
        """

    # Use False Alarm Probability (on white noise only run)
    # threshold at 0.1%, or SDE 8.3
    # OPEN: consider to use the more stringent threshold 0.01% / SDE 9.1
    # (the value computed is optimistic as it's based on white noise data, actual data has some correlation varation (poink noise))
    # See:
    # https://transitleastsquares.readthedocs.io/en/latest/FAQ.html#false-alarm-probability

    def flag_fap(false_alarm_probability):
        return np.isnan(false_alarm_probability) or false_alarm_probability > 0.001

    def flag_sde(sde):
        return sde < 8.3

    def flag_snr(snr):
        return snr < 10  # arbitrary for now

    def flag_odd_even(odd_even_mismatch):
        return odd_even_mismatch > 3  # i.e, > 3 sigma

    def flag_empty_transit_count(empty_transit_count):
        return empty_transit_count > 0

    metrics_list = [
        metrics("Period", pg.period_at_max_power),
        metrics("- Period error\n(half width at half max)", pg.period_at_max_power_err),
        metrics(f"Epoch ({pg.transit_time_at_max_power.format.upper()})", pg.transit_time_at_max_power),
        metrics("Duration", pg.duration_at_max_power),
        metrics("Depth (Model)", pg.depth_at_max_power),
        metrics("Depth (Mean)", 1 - pg._TLS_result.depth_mean[0]),  # mean of all in-transit flux, need 1 - TLS_result
        metrics("- Depth (Mean) error", pg._TLS_result.depth_mean[1]),
        metrics(
            "Period grid size", len(pg.period)
        ),  # larger period grid enerally results in more accurate period and better SDE
        metrics("FAP (white noises)", pg.false_alarm_probability, flag_fap),
        metrics("SDE (Power)", pg.max_power, flag_sde),
        metrics("SR (Log likelihood)", np.max(pg.sr)),  # TODO not too useful in normaized form, it's practically always 1
        metrics("SNR", pg.snr_at_max_power, flag_snr),
        metrics("SNR per transit", pg._TLS_result.snr_per_transit),
        metrics("SNR (pink)\nconerns if greatly different from SNR", pg._TLS_result.snr_pink_per_transit),
        metrics("Empty Transits\nactual period might be double", pg._TLS_result.empty_transit_count, flag_empty_transit_count),
        metrics("Odd-even mismatch (sigma)", pg._TLS_result.odd_even_mismatch, flag_odd_even),
        metrics("Elapsed time (ms)", getattr(pg, "elapsed_time", "N/A")),
    ]

    html = output_metrics(metrics_list)
    if to_display:
        return display(HTML(html))
    else:
        return html


def iterative_bls(
    lc,
    num_iterations,
    pg_kwargs=dict(),
    duration_factor_for_mask=2,
    plot_pg=True,
    plot_lc_model=dict(plot_lc=False, plot_model=True, plot_folded_model=False),
):
    """Iteratively run BLS, to find multiple sets of transit/eclipse like signals."""

    import lightkurve_ext_pg_runner as lke_pg_runner

    result_list = []

    lc_in = lc.copy()
    for i in range(1, num_iterations + 1):
        lc_in.meta["LABEL"] = f"{lc.meta.get('LABEL')}, #{i}"

        display(HTML(f"<h3>Iteration {i}<h3>"))
        result = lke_pg_runner.run_bls(lc_in, pg_kwargs, plot_pg=plot_pg, plot_lc_model=plot_lc_model)

        result_list.append(result)

        # remove identified dips from the LC, then fit it to the next iteration
        t0 = result.pg.transit_time_at_max_power
        period = result.pg.period_at_max_power
        # duration_factor_for_mask: a factor to mask out extra time surrounding dips identified by the model
        # useful when the model's duration it too short, which would leave residual dips that could
        # confuse subsequent BLS runs.
        duration = result.pg.duration_at_max_power * duration_factor_for_mask
        tmask = lc_in.create_transit_mask(period=period, transit_time=t0, duration=duration)
        lc_in = lc_in[~tmask]

    return result_list


def iterative_sine_fit(
    lc,
    num_iterations,
    mask_for_model=None,
    pg_kwargs=dict(),
    plot_kwargs=dict(figsize=(30, 5), s=4, alpha=0.5),
    plot_diagnostics=False,
):
    """Remove sine-wave like periodic signals using iterative sine fitting

    Based on:
    https://docs.lightkurve.org/tutorials/3-science-examples/periodograms-measuring-a-rotation-period.html#5.-Removing-Periodic-Signals-Using-Iterative-Sine-Fitting
    """
    import tic_plot as tplt  # for plot

    lc = lc.normalize()  # use normalized as the base so that it can compute with the model lcs easily later on
    lc = lc["time", "flux", "flux_err"]  # reduce the input size to iterative_sine_fit

    if plot_diagnostics:
        axs = tplt.plot_skip_data_gap(lc, **plot_kwargs)
        axs[0].set_title("Input LC")

    pgs, lc_models, lc_residuals = [], [], []
    lc_in = lc
    for i in range(1, num_iterations + 1):
        lc_4_pg = lc_in
        if mask_for_model is not None:  # the optional mask is to exclude cadence that could skew the periodogram calculation
            lc_4_pg = lc_in[~mask_for_model]
        pg = lc_4_pg.to_periodogram(method="lombscargle", **pg_kwargs)

        lc_model = pg.model(lc_in.time, pg.frequency_at_max_power)
        lc_model.meta["LS_MODEL_ITERATION"] = i
        lc_residual = lc_in.copy()
        lc_residual.flux = lc_in.flux / lc_model.flux
        lc_residual.meta["LS_RESIDUAL_ITERATION"] = i

        if plot_diagnostics:
            axs = tplt.plot_skip_data_gap(lc_residual, label=f"lc_residual{i}", **plot_kwargs)
            axs[0].set_title(f"Iteration {i}; signals removed: period={pg.period_at_max_power}, power={pg.max_power}")
        #             axs = tplt.plot_skip_data_gap(lc_model, label=f"lc_model{i}", **plot_kwargs);

        # accumulate output
        pgs.append(pg)
        lc_models.append(lc_model)
        lc_residuals.append(lc_residual)

        # Send the residual to the next iteration
        lc_in = lc_residual

    return dict(pgs=pg, lc_models=lc_models, lc_residuals=lc_residuals, lc_input=lc)
