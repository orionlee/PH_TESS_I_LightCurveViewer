# ETV based on trapezoid model of eclipses

import numpy as np

from astropy.modeling.functional_models import Trapezoid1D

from etv_functions import plot_initial_guess_of_model, fit_each_eclipse_of_model, run_mcmc_initial_fit_of_model


def trapezoid_model_fit(x, val_constant, amplitude, t0, bottom_width, slope):
    # Note: unlike cosh model, amplitude here must be positive
    return val_constant - Trapezoid1D.evaluate(x, amplitude, t0, bottom_width, slope)


def plot_initial_guess(data, ph_binned, flux_binned, err_binned, *start_vals, ax=None, **kwargs):
    return plot_initial_guess_of_model(
        trapezoid_model_fit, data, ph_binned, flux_binned, err_binned, *start_vals, ax=ax, **kwargs
    )


def log_prior(theta):
    # val_constant, amplitude, t0, bottom_width, slope = theta
    return 0.0
    # return -np.inf  # case reject the params


def log_likelihood(theta, x, y, yerr):
    val_constant, amplitude, t0, bottom_width, slope = theta
    # essentially inlining trapezoid_model_fit()
    model = val_constant - Trapezoid1D.evaluate(x, amplitude, t0, bottom_width, slope)

    return -0.5 * np.sum((y - model) ** 2 / (yerr**2))


def log_probability(theta, x, y, yerr):
    # check that the priors are satisfied
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, x, y, yerr)


def run_mcmc_initial_fit(
    data,
    start_vals,
    nruns=1000,
    discard=600,
    thin=15,
    pool=None,
    plot_chains=False,
    plot=True,
    also_return_stats=False,
    **kwargs,
):
    val_constant, amplitude, t0, bottom_width, slope = start_vals
    start_vals_dict = dict(
        val_constant=val_constant,
        amplitude=amplitude,
        t0=t0,
        bottom_width=bottom_width,
        slope=slope,
    )
    return run_mcmc_initial_fit_of_model(
        log_probability,
        trapezoid_model_fit,
        data,
        start_vals_dict,
        nruns=nruns,
        discard=discard,
        thin=thin,
        pool=pool,
        plot_chains=plot_chains,
        plot=plot,
        also_return_stats=also_return_stats,
        **kwargs,
    )


#
# for fitting individual eclipses
#


def log_prior_fitting(theta):
    return 0.0
    # return -np.inf  # case reject the params


def log_likelihood_fitting(theta, x, y, yerr, mean_bottom_width, mean_slope):
    val_constant, amplitude, t0 = theta
    # essentially inlining trapezoid_model_fit()
    # TODO: consider to replace Trapezoid1D.evaluate() with some more optimized codes,
    # as it appears to be slower than the cosh gaussian model
    # (especially for fitting individual eclipses here)
    model = val_constant - Trapezoid1D.evaluate(x, amplitude, t0, mean_bottom_width, mean_slope)

    return -0.5 * np.sum((y - model) ** 2 / (yerr**2))


def log_probability_fitting(theta, x, y, yerr, mean_bottom_width, mean_slope):
    # check that the priors are satisfied
    lp = log_prior_fitting(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood_fitting(theta, x, y, yerr, mean_bottom_width, mean_slope)


def fit_each_eclipse(
    data,
    n_transits,
    t0,
    period,
    mean_val_constant,
    mean_amplitude,
    mean_t0,
    mean_bottom_width,
    mean_slope,
    outfile_path,
    pool=None,
    min_number_data=20,
):
    start_vals_dict = dict(
        val_constant=mean_val_constant,
        amplitude=mean_amplitude,
        # note: mean_t0 is not used
    )
    fixed_vals_dict = dict(bottom_width=mean_bottom_width, slope=mean_slope)
    return fit_each_eclipse_of_model(
        log_probability_fitting,  # the version for individual eclipse
        log_prior_fitting,  # the version for individual eclipse
        trapezoid_model_fit,
        data,
        n_transits,
        t0,
        period,
        start_vals_dict,
        fixed_vals_dict,
        outfile_path,
        pool=pool,
        min_number_data=min_number_data,
    )
