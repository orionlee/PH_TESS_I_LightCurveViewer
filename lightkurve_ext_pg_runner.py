import logging
from time import time_ns
from types import SimpleNamespace
import warnings

import ipywidgets as widgets

from IPython.display import display

import lightkurve_ext as lke
import lightkurve_ext_tls as lke_tls
import lightkurve_ext_pg as lke_pg


def _current_time_millis():
    return time_ns() / 1000000


def _flatten(lc, flatten_kwargs):
    if flatten_kwargs is None:
        return lc
    flatten_kwargs = flatten_kwargs.copy()
    window_length_in_days = flatten_kwargs.pop("window_length_in_days", None)
    if window_length_in_days is not None:
        window_length = lke.to_window_length_for_2min_cadence(window_length_in_days)
        flatten_kwargs["window_length"] = window_length
    return lc.flatten(**flatten_kwargs)


def run_tls(lc, pg_kwargs={}, flatten_kwargs=None, plot_pg=True, plot_lc_model=True, plot_transit_depth=True):
    lc = lc.remove_nans().normalize()
    lc = _flatten(lc, flatten_kwargs)
    time_b = _current_time_millis()
    pg = lke_tls.TransitLeastSquaresPeriodogram.from_lightcurve(lc, **pg_kwargs)
    time_e = _current_time_millis()
    pg.elapsed_time = time_e - time_b

    lke_pg.validate_tls_n_report(pg)

    ax_pg = None
    if plot_pg:
        ax_pg = lke_pg.plot_pg_n_mark_max(pg)

    ax_lc_model_1, ax_lc_model_2, ax_lc_model_f = None, None, None
    if plot_lc_model:
        ax_lc_model_1, ax_lc_model_2, ax_lc_model_f = lke_pg.plot_lc_with_model(lc, pg)

    ax_tt_depth = None
    if plot_transit_depth:
        ax_tt_depth = lke_pg.errorbar_transit_depth(pg)

    return SimpleNamespace(
        pg=pg,
        lc=lc,
        ax_pg=ax_pg,
        ax_lc_model_1=ax_lc_model_1,
        ax_lc_model_2=ax_lc_model_2,
        ax_lc_model_f=ax_lc_model_f,
        ax_tt_depth=ax_tt_depth,
    )


def run_bls(lc, pg_kwargs={}, flatten_kwargs=None, plot_pg=True, plot_lc_model=True):
    lc = lc.remove_nans().normalize()
    lc = _flatten(lc, flatten_kwargs)
    time_b = _current_time_millis()
    pg = lc.to_periodogram(method="bls", **pg_kwargs)
    time_e = _current_time_millis()
    pg.elapsed_time = time_e - time_b

    lke_pg.validate_bls_n_report(pg)

    ax_pg = None
    if plot_pg:
        ax_pg = lke_pg.plot_pg_n_mark_max(pg)
    ax_lc_model_1, ax_lc_model_2, ax_lc_model_f = None, None, None

    if plot_lc_model:
        with warnings.catch_warnings():
            # avoid warnings about using max power values
            warnings.filterwarnings("ignore", message=".*Using.*")
            logger = logging.getLogger("lightkurve.periodogram")
            logger.setLevel(logging.ERROR)
            ax_lc_model_1, ax_lc_model_2, ax_lc_model_f = lke_pg.plot_lc_with_model(lc, pg)

    ax_tt_depth = None
    # ax_tt_depth = lke_pg.errorbar_transit_depth(pg) # bls has no info directly

    return SimpleNamespace(
        pg=pg,
        lc=lc,
        ax_pg=ax_pg,
        ax_lc_model_1=ax_lc_model_1,
        ax_lc_model_2=ax_lc_model_2,
        ax_lc_model_f=ax_lc_model_f,
        ax_tt_depth=ax_tt_depth,
    )


def run_bls_n_tls(lc, plot_pg=True, plot_lc_model=True, plot_transit_depth=True, bls_pg_kwargs={}, tls_pg_kwargs={}):
    # Run TLS and BLS and have their results displayed side-by-side.
    #
    # For the matplotlib figures to be displayed inside the respective boxes in Jupyter, magic
    # %matplotlib widget
    # is needed (requiring ipympl package)
    #
    # sometimes it crashes the browsers (possibly too many interactive figures?!)

    out_bls = widgets.Output(layout={"border": "1px solid lightgray"})
    out_tls = widgets.Output(layout={"border": "1px solid lightgray"})
    ctr = widgets.GridBox(
        children=[out_bls, out_tls],
        layout=widgets.Layout(width="auto", grid_template_rows="auto", grid_template_columns="50% 50%", grid_gap="5px 10px"),
    )

    with out_bls:
        run_bls(lc, bls_pg_kwargs, plot_pg=plot_pg, plot_lc_model=plot_lc_model)
    with out_tls:
        run_tls(lc, tls_pg_kwargs, plot_pg=plot_pg, plot_lc_model=plot_lc_model, plot_transit_depth=plot_transit_depth)
    return display(ctr)
