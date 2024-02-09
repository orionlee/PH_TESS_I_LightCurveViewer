#
# A Cosh Gaussian model including period `p` as the parameters
#
# The function names followed the convention of having a `_p` suffix,
# comparing to their counterpart in `etv_functions.py`
#

import logging
import os
import multiprocessing
from multiprocessing import Pool

import numpy as np
import matplotlib.pyplot as plt
import emcee


from etv_functions import phase_data, get_starting_positions

log = logging.getLogger(__name__)


def _parse_pool_param(pool):
    """
    Parse the `pool` parameter shared by various mcmc functions.
    It allows caller to pass a `Pool` instance or various shorthands
    """
    if pool is None:
        return None, False  # <-- is pool_from_caller

    if isinstance(pool, multiprocessing.pool.Pool):
        # multiprocessing.pool.Pool is the class, while
        # multiprocessing.Pool is a factory function to create the instances
        is_pool_from_caller = True
    else:
        is_pool_from_caller = False
        if pool == "default":
            # defaulted to use about 80% of the CPUs
            num_processes_to_use = max(1, int(os.cpu_count() * 0.8))
            log.info(f"emcee parallel enabled, defaulted to use {num_processes_to_use} CPUs.")
            pool = Pool(num_processes_to_use)
        elif pool == "all":
            num_processes_to_use = os.cpu_count()
            log.info(f"emcee parallel enabled, use all {num_processes_to_use} CPUs.")
            pool = Pool(num_processes_to_use)
        elif isinstance(pool, int):
            num_processes_to_use = pool
            log.info(f"emcee parallel enabled, use {num_processes_to_use} CPUs.")
            pool = Pool(num_processes_to_use)
        else:
            raise TypeError('`pool` must be None, "default", int (for num processes), or a `Pool` instance.')
        return pool, is_pool_from_caller


def run_mcmc_initial_fit_p(
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
    figsize = kwargs.get("figsize", (8, 4))
    figsize_chains = kwargs.get("figsize_chains", (8, 12))

    pool, is_pool_from_caller = _parse_pool_param(pool)

    try:
        if pool is not None:
            # case parallel emcee is enabled,
            # turn off numpy parallel operations to avoid possible conflicts:
            # see: https://emcee.readthedocs.io/en/stable/tutorials/parallel/#parallel
            os.environ["OMP_NUM_THREADS"] = "1"

        pos = list(get_starting_positions(start_vals, nwalkers=128))[0]

        nwalkers = 128
        ndim = len(start_vals)

        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability_p, args=(data.time, data.flux, data.err), pool=pool)

        sampler.run_mcmc(pos, nruns, progress=True, store=True)

        tau = sampler.get_autocorr_time(tol=0)

        samples = sampler.get_chain()
        labels = ["alpha0", "alpha1", "t0", "d", "Tau", "p"]

        if plot_chains == True:

            fig, axes = plt.subplots(ndim, figsize=figsize_chains, sharex=True)

            for i in range(ndim):
                ax = axes[i]
                ax.plot(samples[:, :, i], "k", alpha=0.3)
                ax.set_xlim(0, len(samples))
                ax.set_ylabel(labels[i])
                ax.yaxis.set_label_coords(-0.1, 0.5)

            axes[-1].set_xlabel("step number")
            plt.show()

        flat_samples = sampler.get_chain(discard=discard, thin=thin, flat=True)

        mean_alpha0 = np.median(flat_samples[:, 0])
        mean_alpha1 = np.median(flat_samples[:, 1])
        mean_t0 = np.median(flat_samples[:, 2])
        mean_d = np.median(flat_samples[:, 3])
        mean_Tau = np.median(flat_samples[:, 4])
        mean_p = np.median(flat_samples[:, 5])

        inds = np.random.randint(len(flat_samples), size=nruns)

        if plot == True:
            fig, axes = plt.subplots(figsize=figsize, sharex=True)

            for ind in inds:
                sample = flat_samples[ind]
                plt.plot(
                    data.phase,
                    coshgauss_model_fit_p(data.time, sample[0], sample[1], sample[2], sample[3], sample[4], sample[5]),
                    "C1",
                    lw=0,
                    marker=".",
                    markersize=0.1,
                    alpha=0.1,
                    zorder=1,
                )

            plt.errorbar(data.phase, data.flux, yerr=data.err, fmt=".k", capsize=0, zorder=-2)

            plt.plot(
                data.phase,
                coshgauss_model_fit_p(data.time, mean_alpha0, mean_alpha1, mean_t0, mean_d, mean_Tau, mean_p),
                lw=0,
                marker=".",
                markersize=0.5,
                alpha=1,
                zorder=2,
                color="red",
            )

            plt.xlabel("x")
            plt.ylabel("y")
            plt.show()

        if not also_return_stats:
            return mean_alpha0, mean_alpha1, mean_t0, mean_d, mean_Tau, mean_p
        else:
            stats = dict(
                std_alpha0=np.std(flat_samples[:, 0]),
                std_alpha1=np.std(flat_samples[:, 1]),
                std_t0=np.std(flat_samples[:, 2]),
                std_d=np.std(flat_samples[:, 3]),
                std_Tau=np.std(flat_samples[:, 4]),
                std_p=np.std(flat_samples[:, 5]),
            )
            return mean_alpha0, mean_alpha1, mean_t0, mean_d, mean_Tau, mean_p, stats
    finally:
        if pool is not None:
            del os.environ["OMP_NUM_THREADS"]  # TODO: should restore the original value
        if pool is not None and not is_pool_from_caller:
            pool.close()


def log_prior_p(theta):
    alpha0, alpha1, t0, d, Tau, p = theta

    # no t0 check, as t0 is time, not phase
    if (0 < alpha0 < 10) and (-10 < alpha1 < 0) and (0 < d < 10) and (0.5 < Tau < 50):
        return 0.0
    return -np.inf


def log_likelihood_p(theta, x, y, yerr):

    alpha0, alpha1, t0, d, Tau, p = theta

    x_phase = phase_data(x, t0, p)
    t0_phase = phase_data([t0], t0, p)[0]
    cosh_term = np.cosh((x_phase - t0_phase) / d)
    exp_term = np.exp(1 - cosh_term)
    pow_term = pow((1 - exp_term), Tau)

    psi = 1 - pow_term

    model = alpha0 + (alpha1 * psi)

    return -0.5 * np.sum((y - model) ** 2 / (yerr**2))


def log_probability_p(theta, x, y, yerr):

    # check that the priors are satisfied
    lp = log_prior_p(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood_p(theta, x, y, yerr)


def coshgauss_model_fit_p(x, alpha0, alpha1, t0, d, Tau, p):
    """
    Calculates the coshgauss model fit for the given parameters.

    Parameters:
    x (float): The input value (time, NOT phase)
    alpha0 (float): The coefficient for the constant term.
    alpha1 (float): The coefficient for the psi term.
    t0 (float): The center of the cosh term. (time, NOT phase)
    d (float): The width of the cosh term.
    Tau (float): The exponent for the pow term.
    p (float): the period

    Returns:
    float: The model fit value.
    """
    x_phase = phase_data(x, t0, p)
    t0_phase = phase_data([t0], t0, p)[0]
    cosh_term = np.cosh((x_phase - t0_phase) / d)
    exp_term = np.exp(1 - cosh_term)
    pow_term = pow((1 - exp_term), Tau)

    psi = 1 - pow_term

    model = alpha0 + (alpha1 * psi)

    return model
