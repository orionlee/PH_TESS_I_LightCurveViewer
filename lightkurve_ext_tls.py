from __future__ import annotations
import logging
from time import time_ns
import warnings

from astropy.time import Time
import astropy.units as u
import numpy as np

from lightkurve import LightCurve, FoldedLightCurve, LightkurveWarning
from lightkurve.periodogram import Periodogram, BoxLeastSquaresPeriodogram

# for type annotation
from numbers import Number
from typing import Tuple, Union

# for the pattern float or Quantity
QuantityLike = Union[u.Quantity, Number]

log = logging.getLogger(__name__)

# Derived from https://github.com/lightkurve/lightkurve/pull/476


def _set_if_not_exists(dict_obj, key_value_pairs, value_converter_func=None):
    def _isnan(value):
        res = np.isnan(value)
        return res if np.isscalar(res) else np.all(res)

    for key, value in key_value_pairs:
        if value is not None and (not _isnan(value)):
            if value_converter_func is not None:
                value = value_converter_func(value)
            if dict_obj.get(key, None) is None:
                dict_obj[key] = value
            else:
                warnings.warn(
                    f"The argument {key} is specified twice. Use value {dict_obj.get(key)}", LightkurveWarning, stacklevel=2
                )


def _set_min_max_if_needed(dict_obj, key, min_default, max_default):
    # provide some min/max if not available from catalog
    value = dict_obj.get(key)
    if value is not None:
        if dict_obj.get(key + "_min") is None:
            min = value * 0.5
            if min > min_default:
                min = min_default
            dict_obj[key + "_min"] = min
        if dict_obj.get(key + "_max") is None:
            max = value * 2
            if max < max_default:
                max = max_default
            dict_obj[key + "_max"] = max


def _time_like(time_template, val, val2=None):
    return Time(
        val,
        val2=val2,
        format=time_template.format,
        scale=time_template.scale,
        precision=time_template.precision,
        in_subfmt=time_template.in_subfmt,
        out_subfmt=time_template.out_subfmt,
        location=time_template.location,
    )


def _current_time_millis():
    return time_ns() / 1000000


def _catalog_info(lc: LightCurve):
    try:
        from transitleastsquares import catalog_info
    except ImportError:
        raise Exception(
            "This feature requires the `transitleastsquares` package. "
            "You can install it using `pip install transitleastsquares`."
        )

    ab, mass, mass_min, mass_max, radius, radius_min, radius_max = None, None, None, None, None, None, None

    mission_to_ci_arg_name = {"TESS": "TIC_ID", "Kepler": "KIC_ID", "K2": "EPIC_ID"}
    target_id = lc.meta.get("TARGETID")
    id_key_name = mission_to_ci_arg_name.get(lc.meta.get("MISSION"))
    if id_key_name is not None and target_id is not None:
        ci_kwargs = {}
        ci_kwargs[id_key_name] = target_id
        time_b = _current_time_millis()
        ab, mass, mass_min, mass_max, radius, radius_min, radius_max = catalog_info(**ci_kwargs)
        time_e = _current_time_millis()
        log.debug(f"catalog_info() elapsed time: {time_e - time_b}ms")
        # return value is weird ndarray with no dimension (shape ())
        mass, mass_min, mass_max, radius, radius_min, radius_max = [
            mass.flatten()[0],
            mass_min.flatten()[0],
            mass_max.flatten()[0],
            radius.flatten()[0],
            radius_min.flatten()[0],
            radius_max.flatten()[0],
        ]
        # the returned min/max is actually error, at least for TESS
        if mass_min is not None:
            mass_min = mass - mass_min
        if mass_max is not None:
            mass_max = mass + mass_max
        if radius_min is not None:
            radius_min = radius - radius_min
        if radius_max is not None:
            radius_max = radius + radius_max
    else:
        raise ValueError(f"No supported Target ID: {id_key_name} {target_id}")

    return ab, mass, mass_min, mass_max, radius, radius_min, radius_max


class TransitLeastSquaresPeriodogram(Periodogram):
    """Subclass of :class:`Periodogram <lightkurve.periodogram.Periodogram>`
    representing a power spectrum generated using the Transit Least Squares (TLS) method.
    """

    def _init_from_kwargs(self, kwargs, keys):
        for key in keys:
            value = kwargs.pop(key, None)
            setattr(self, key, value)

    def __init__(self, *args, **kwargs):
        self._TLS_result = kwargs.pop("tls_result", None)
        self._TLS_object = kwargs.pop("tls_obj", None)
        self._init_from_kwargs(
            kwargs,
            [
                # attributes that are the same as BLS
                "transit_time",
                "transit_time_at_max_power",
                "time",
                "time_unit",
                "flux",
                # additional per-period/frequency info
                "sr",
                "chi2",
                "chi2red",
                # additional summary info
                # note: TLS impl, unlike BLS impl, does not provide per-period depth, duration, and snr
                "depth_at_max_power",
                "duration_at_max_power",
                "snr_at_max_power",
                "period_at_max_power_err",
                "false_alarm_probability",  # use the term false_alarm_probability, based on astropy's LombScargle
                # TODO: additional info, e.g., odd_even_mismatch
                # TODO: consider to move some of the statistics to a separate compute_stats() function
                "transit_depth",
                "transit_depth_err",
                # TODO: additional per-transit info
            ],
        )
        super(TransitLeastSquaresPeriodogram, self).__init__(*args, **kwargs)

    def __repr__(self):
        return "TransitLeastSquaresPeriodogram(ID: {})".format(self.targetid)

    @staticmethod
    def from_lightcurve(lc: LightCurve, **kwargs) -> TransitLeastSquaresPeriodogram:
        """Creates a Periodogram from a LightCurve using the TLS method."""
        try:
            from transitleastsquares import transitleastsquares, catalog_info, tls_constants
        except ImportError:
            raise Exception(
                "This feature requires the `transitleastsquares` package. "
                "You can install it using `pip install transitleastsquares`."
            )

        lc = lc.remove_nans()

        log.debug(f"TLS.from_lightcurve() args: {kwargs}")
        if np.isfinite(lc.flux_err).all():
            flux_err = lc.flux_err
        else:
            flux_err = None
        tls = transitleastsquares(lc.time.value, lc.flux.value, flux_err)

        # TODO: convert lk Periodgram standard arguments to tls.power() arguments,
        def _to_unitless_day(val):
            if hasattr(val, "value"):
                return val.to(u.day).value
            else:
                return val

        _set_if_not_exists(
            kwargs,
            [
                ("period_min", kwargs.pop("minimum_period", None)),
                ("period_max", kwargs.pop("maximum_period", None)),
            ],
            value_converter_func=_to_unitless_day,
        )

        try:
            derive_stellar_priors = kwargs.pop("derive_stellar_priors", True)
            if derive_stellar_priors:
                ab, mass, mass_min, mass_max, radius, radius_min, radius_max = _catalog_info(lc)
                _set_if_not_exists(
                    kwargs,
                    [
                        ("u", ab),
                        ("R_star", radius),
                        ("R_star_min", radius_min),
                        ("R_star_max", radius_max),
                        ("M_star", mass),
                        ("M_star_min", mass_min),
                        ("M_star_max", mass_max),
                    ],
                )
                _set_min_max_if_needed(kwargs, "R_star", tls_constants.R_STAR_MIN, tls_constants.R_STAR_MAX)
                _set_min_max_if_needed(kwargs, "M_star", tls_constants.M_STAR_MIN, tls_constants.M_STAR_MAX)

            log.debug(f"TLS.from_lightcurve() - tls.power() args: {kwargs}")
        except Exception as e:
            warnings.warn(f"TLS.from_lightcurve(): cannot derive stellar priors. Use defaults. Reason: {e}", LightkurveWarning)

        result = tls.power(**kwargs)
        if not isinstance(result.period, u.quantity.Quantity):
            result.periods = u.Quantity(result.periods, u.day)
        if not isinstance(result.power, u.quantity.Quantity):
            result.power = result.power * u.dimensionless_unscaled

        return TransitLeastSquaresPeriodogram(
            # attributes that are the same as BLS
            default_view="period",
            label=lc.meta.get("LABEL"),
            targetid=lc.meta.get("TARGETID"),
            frequency=1.0 / result.periods,
            power=result.power,
            transit_time=_time_like(lc.time, result.transit_times),
            time=lc.time,
            time_unit="day",
            flux=lc.flux,
            # TLS-specific per-period / frequency info
            sr=result.SR,
            chi2=result.chi2,
            chi2red=result.chi2red,
            # TLS-specific summary info
            duration_at_max_power=result.duration * u.day,
            depth_at_max_power=1 - result.depth,  # lk convention: depth is the dip's depth.
            transit_time_at_max_power=_time_like(lc.time, result.T0),
            snr_at_max_power=result.snr,
            period_at_max_power_err=result.period_uncertainty * u.day,
            false_alarm_probability=result.FAP,
            # TLS-specific per-transit info
            transit_depth=1 - result.transit_depths,
            transit_depth_err=result.transit_depths_uncertainties,
            # TLS impl
            tls_result=result,
            tls_obj=tls,
        )

    def get_transit_model(self):
        return LightCurve(
            time=_time_like(self.time, self._TLS_result.model_lightcurve_time),
            flux=self._TLS_result.model_lightcurve_model,
            meta=dict(LABEL=f"{self.label} Transit Model Flux"),
        )

    def get_transit_mask(self, period: QuantityLike = None, duration: QuantityLike = None, transit_time: QuantityLike = None):
        from transitleastsquares import transit_mask

        if period is None:
            period = self.period_at_max_power.to(u.day).value
        if duration is None:
            duration = self.duration_at_max_power.to(u.day).value
        if transit_time is None:
            transit_time = self.transit_time_at_max_power.value

        return transit_mask(self.time.value, period, duration, transit_time)

    def plot(self, **kwargs):
        # TODO: support a new parameter, spectrum="power", (to plot sr, chi2, etc.)
        ax = super(TransitLeastSquaresPeriodogram, self).plot(**kwargs)
        if "ylabel" not in kwargs:
            ax.set_ylabel("TLS SDE")
        return ax

    def fold(
        self, lc: LightCurve, period: QuantityLike = None, transit_time: QuantityLike = None, **kwargs
    ) -> Tuple[FoldedLightCurve, FoldedLightCurve]:
        if period is None:
            period = self.period_at_max_power
        if transit_time is None:
            transit_time = self.transit_time_at_max_power

        return (
            lc.fold(period=period, epoch_time=transit_time, **kwargs),
            self.get_transit_model().fold(period=period, epoch_time=transit_time, **kwargs),
        )


def create_bls_pg_with_stellar_specific_search_grid(lc: LightCurve, **kwargs) -> BoxLeastSquaresPeriodogram:
    """Run BLS using the stellar-specific grid from TLS implementation."""
    # See: https://github.com/hippke/tls/blob/master/tutorials/09%20Optimal%20period%20grid%20and%20optimal%20duration%20grid.ipynb
    # for background
    try:
        from transitleastsquares import period_grid, duration_grid, tls_constants
    except ImportError:
        raise Exception(
            "This feature requires the `transitleastsquares` package. "
            "You can install it using `pip install transitleastsquares`."
        )

    def _to_absolute(duration_in_fraction, period, log_step=tls_constants.DURATION_GRID_STEP):
        duration_min, duration_max = duration_in_fraction[0] * period[0], duration_in_fraction[-1] * period[-1]
        # redoing what TLS duration_grid() does to create the grid given min, max
        # - essentially creating a geometric space, with the log_step as the multiple
        durations = [duration_min]
        current_depth = duration_min
        while current_depth * log_step < duration_max:
            current_depth = current_depth * log_step
            durations.append(current_depth)
        durations.append(duration_max)  # Append endpoint. Not perfectly spaced.
        return durations

    ab, mass, mass_min, mass_max, radius, radius_min, radius_max = _catalog_info(lc)
    if mass is not None:
        # TODO: handle optional parameters, in particular, minimum_period and maximum_period
        period = period_grid(radius, mass, (lc.time.max() - lc.time.min()).value)
        duration_in_fraction = duration_grid(period, shortest=None)  # shortest not used by implementation
        # convert the duration grid, in fraction of period, to one with absolute value.
        # It is less accurate, because for a given triad period, the duration grid is the same
        # (rather than specific to the period).
        # But astropy BLS does not support expressing duration in fractions of periods
        duration = _to_absolute(duration_in_fraction, period)
        log.debug(
            f"""\
To Run BLS with grid: period: ({len(period)}){period[0]} - {period[-1]}  ; \
duration: ({len(duration)}){duration[0]} - {duration[-1]}"""
        )
    else:
        period, duration = None, None

    return lc.to_periodogram(method="bls", period=period, duration=duration, **kwargs)
