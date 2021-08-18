# Various helpers to work with Periodogram
#
import warnings

from types import SimpleNamespace
from memoization import cached

import numpy as np

from lightkurve.periodogram import BoxLeastSquaresPeriodogram

from IPython.display import display, HTML

# for type annotation
from numbers import Number


def plot_pg_n_mark_max(pg, ax=None, max_period_factor=None):
    ax = pg.plot(ax=ax, view="period")
    ax.axvline(pg.period_at_max_power.value, c="blue", alpha=0.4)
    max_text = f"Power: {pg.max_power:.2f}, Period: {pg.period_at_max_power:.4f}"
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


def _model(pg, lc):
    if hasattr(pg, "get_transit_model"):  # case BLS pg
        return pg.get_transit_model()
    else:  # case LS pg
        return pg.model(lc.time, pg.frequency_at_max_power)


def plot_lc_with_model(lc, pg):
    # with plt.style.context(lk.MPLSTYLE):
    #     ax = plt.figure(figsize=(15, 6)).gca()
    ax1 = lc.scatter()
    if hasattr(pg, "get_transit_mask"):
        lc[pg.get_transit_mask()].scatter(ax=ax1, c="orange", marker="x", s=9, label="in transits")

    lc_model = _model(pg, lc)

    ax2 = lc.scatter()
    lc_model.plot(ax=ax2, c="red", alpha=0.9, linewidth=2)

    # folded, zoom -in
    period = pg.period_at_max_power
    if hasattr(pg, "transit_time_at_max_power"):
        epoch_time = pg.transit_time_at_max_power
    else:
        epoch_time = None
    lc_f = lc.fold(epoch_time=epoch_time, period=period)
    lc_model_f = lc_model.fold(epoch_time=epoch_time, period=period)

    ax_f = lc_f.scatter()
    lc_model_f.scatter(ax=ax_f, c="red")
    if hasattr(pg, "duration_at_max_power"):
        # zoom in for BLS model:
        ax_f.set_xlim(-pg.duration_at_max_power.value, pg.duration_at_max_power.value)

    return ax1, ax2, ax_f


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
