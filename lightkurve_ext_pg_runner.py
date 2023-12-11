import contextlib
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


def _remove_fig_title(*ax_args):
    # Used to remove the extra title in %matplotlib widget mode
    # alternative would be disable them globally, see
    # https://github.com/matplotlib/ipympl/issues/229#issuecomment-633430427
    for ax in ax_args:
        if ax is not None:
            ax.get_figure().canvas.header_visible = False
            ax.get_figure().canvas.footer_visible = False
            # ax.get_figure().canvas.toolbar_visible = False
            # ax.get_figure().canvas.resizable = False


def run_tls(
    lc, pg_kwargs={}, flatten_kwargs=None, plot_pg=True, plot_lc_model=True, plot_transit_depth=True, display_context=None
):
    if display_context is None:
        # note : nullcontext() requires Python 3.7
        ctx_validate, ctx_plot = contextlib.nullcontext(), contextlib.nullcontext()
    else:
        ctx_validate, ctx_plot = display_context["validate"], display_context["plot"]

    with ctx_validate:
        lc = lc.remove_nans().normalize()
        lc = _flatten(lc, flatten_kwargs)
        time_b = _current_time_millis()
        pg = lke_tls.TransitLeastSquaresPeriodogram.from_lightcurve(lc, **pg_kwargs)
        time_e = _current_time_millis()
        pg.elapsed_time = time_e - time_b

        lke_pg.validate_tls_n_report(pg)

    with ctx_plot:
        ax_pg = None
        if plot_pg:
            ax_pg = lke_pg.plot_pg_n_mark_max(pg)
            _remove_fig_title(ax_pg)

        ax_lc_model_1, ax_lc_model_2, ax_lc_model_f = None, None, None
        if plot_lc_model:
            ax_lc_model_1, ax_lc_model_2, ax_lc_model_f = lke_pg.plot_lc_with_model(lc, pg)
            _remove_fig_title(ax_lc_model_1, ax_lc_model_2, ax_lc_model_f)

        ax_tt_depth = None
        if plot_transit_depth:
            ax_tt_depth = lke_pg.errorbar_transit_depth(pg)
            _remove_fig_title(ax_tt_depth)

    return SimpleNamespace(
        pg=pg,
        lc=lc,
        ax_pg=ax_pg,
        ax_lc_model_1=ax_lc_model_1,
        ax_lc_model_2=ax_lc_model_2,
        ax_lc_model_f=ax_lc_model_f,
        ax_tt_depth=ax_tt_depth,
    )


def run_bls(
    lc,
    use_stellar_specific_search_grid=False,
    pg_kwargs={},
    flatten_kwargs=None,
    plot_pg=True,
    plot_lc_model=True,
    display_context=None,
):
    if display_context is None:
        ctx_validate, ctx_plot = contextlib.nullcontext(), contextlib.nullcontext()
    else:
        ctx_validate, ctx_plot = display_context["validate"], display_context["plot"]

    with ctx_validate:
        lc = lc.remove_nans().normalize()
        lc = _flatten(lc, flatten_kwargs)
        time_b = _current_time_millis()
        if use_stellar_specific_search_grid:
            pg = lke_tls.create_bls_pg_with_stellar_specific_search_grid(lc, **pg_kwargs)
        else:
            pg = lc.to_periodogram(method="bls", **pg_kwargs)
        time_e = _current_time_millis()
        pg.elapsed_time = time_e - time_b

        lke_pg.validate_bls_n_report(pg)

    with ctx_plot:
        ax_pg = None
        if plot_pg:
            ax_pg = lke_pg.plot_pg_n_mark_max(pg)
            _remove_fig_title(ax_pg)

        ax_lc_model_1, ax_lc_model_2, ax_lc_model_f, lcs = None, None, None, None
        if plot_lc_model is not False:
            # convert plot_lc_model to kwargs
            if plot_lc_model is True:
                plot_lc_model_kwargs = dict(plot_lc=True, plot_model=True, plot_folded_model=True)
            elif isinstance(plot_lc_model, dict):
                plot_lc_model_kwargs = plot_lc_model
            else:
                raise TypeError("Argument plot_lc_model is not of supported Type (boolean or dict)")

            with warnings.catch_warnings():
                # avoid warnings about using max power values
                warnings.filterwarnings("ignore", message=".*Using.*")
                logger = logging.getLogger("lightkurve.periodogram")
                logger.setLevel(logging.ERROR)
                [ax_lc_model_1, ax_lc_model_2, ax_lc_model_f], lcs = lke_pg.plot_lc_with_model(
                    lc, pg, also_return_lcs=True, **plot_lc_model_kwargs
                )
                _remove_fig_title(ax_lc_model_1, ax_lc_model_2, ax_lc_model_f)
        else:
            # case not plotting LC model, create a dummy to satisfy the return expression below
            lcs = SimpleNamespace(lc=lc, lc_f=None, lc_model=None, lc_model_f=None)

        ax_tt_depth = None
        # ax_tt_depth = lke_pg.errorbar_transit_depth(pg) # bls has no info directly

    return SimpleNamespace(
        pg=pg,
        lc=lc,
        lc_f=lcs.lc_f,
        lc_model=lcs.lc_model,
        lc_model_f=lcs.lc_model_f,
        ax_pg=ax_pg,
        ax_lc_model_1=ax_lc_model_1,
        ax_lc_model_2=ax_lc_model_2,
        ax_lc_model_f=ax_lc_model_f,
        ax_tt_depth=ax_tt_depth,
    )


def run_bls_n_tls(
    lc,
    use_stellar_specific_search_grid_for_bls=False,
    plot_pg=True,
    plot_lc_model=True,
    plot_transit_depth=True,
    bls_pg_kwargs={},
    tls_pg_kwargs={},
):
    # Run TLS and BLS and have their results displayed side-by-side.
    #
    # For the matplotlib figures to be displayed inside the respective boxes in Jupyter, magic
    # %matplotlib widget
    # is needed (requiring ipympl package)
    #
    # sometimes it crashes the browsers (possibly too many interactive figures?!)

    out_bls_validate = widgets.Output(layout={"border": "0px solid lightgray"})
    out_bls_plot = widgets.Output(layout={"border": "0px solid lightgray"})
    out_tls_validate = widgets.Output(layout={"border": "0px solid lightgray"})
    out_tls_plot = widgets.Output(layout={"border": "0px solid lightgray"})
    ctr = widgets.GridBox(
        children=[out_bls_validate, out_tls_validate, out_bls_plot, out_tls_plot],
        layout=widgets.Layout(width="auto", grid_template_rows="auto", grid_template_columns="50% 50%", grid_gap="5px 10px"),
    )

    run_bls(
        lc,
        use_stellar_specific_search_grid=use_stellar_specific_search_grid_for_bls,
        pg_kwargs=bls_pg_kwargs,
        plot_pg=plot_pg,
        plot_lc_model=plot_lc_model,
        display_context=dict(validate=out_bls_validate, plot=out_bls_plot),
    )

    run_tls(
        lc,
        tls_pg_kwargs,
        plot_pg=plot_pg,
        plot_lc_model=plot_lc_model,
        plot_transit_depth=plot_transit_depth,
        display_context=dict(validate=out_tls_validate, plot=out_tls_plot),
    )

    # with out_bls:
    #     run_bls(lc, bls_pg_kwargs, plot_pg=plot_pg, plot_lc_model=plot_lc_model)
    # with out_tls:
    #     run_tls(lc, tls_pg_kwargs, plot_pg=plot_pg, plot_lc_model=plot_lc_model, plot_transit_depth=plot_transit_depth)
    return display(ctr)
