#
# A Cosh Gaussian model including period `p` as the parameters
#
# The function names followed the convention of having a `_p` suffix,
# comparing to their counterpart in `etv_functions.py`
#

import warnings

import numpy as np

from etv_functions import run_mcmc_initial_fit_of_model, phase_data, do_plot_autocorrelation

# the period-aware coshgauss model and log likelihood functions


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


def run_mcmc_initial_fit_p(
    data,
    start_vals,
    nruns=1000,
    discard=600,
    thin=15,
    log_probability_func=log_probability_p,
    autocorr_time_kwargs=None,
    pool=None,
    plot_chains=False,
    plot_autocorrelation=False,
    plot=True,
    also_return_stats=False,
    **kwargs,
):
    import matplotlib.pyplot as plt
    import emcee
    from etv_functions import phase_data, get_starting_positions, EmceePoolContext, _parse_pool_param

    figsize = kwargs.get("figsize", (8, 4))
    figsize_chains = kwargs.get("figsize_chains", (8, 12))

    pool, is_pool_from_caller = _parse_pool_param(pool)

    with EmceePoolContext(pool, auto_close=not is_pool_from_caller):
        pos = list(get_starting_positions(start_vals, nwalkers=128))[0]

        nwalkers = 128
        ndim = len(start_vals)

        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability_func, args=(data.time, data.flux, data.err), pool=pool)

        sampler.run_mcmc(pos, nruns, progress=True, store=True)

        # integrated autocorrelation time estimate
        # (used for convergence check, and number of independent samples estimate)
        autocorr_time = sampler.get_autocorr_time(tol=0)  # tol=0 to provide estimate without rasing errors

        # issue a warning if the chain is possibly too short for the above estimate,
        # the threshold is primarily tuned with `tol` parameter
        if autocorr_time_kwargs is None:
            autocorr_time_kwargs = dict(tol=50)
        autocorr_time_kwargs["quiet"] = True  # to issue a warning instead of raising an error
        sampler.get_autocorr_time(**autocorr_time_kwargs)

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

        if plot_autocorrelation:
            ax = do_plot_autocorrelation(samples, labels)

            def to_int(t):
                if np.isfinite(t):
                    return int(round(t, 0))
                else:  # handle nan, etc.
                    return t

            ax.set_title(f"Integrated autocorrelation time estimate:\n{[to_int(t) for t in autocorr_time]}")

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
            stats["autocorr_time"] = autocorr_time
            stats["sampler"] = sampler  # the underlying sampler
            return mean_alpha0, mean_alpha1, mean_t0, mean_d, mean_Tau, mean_p, stats


# run_mcmc_initial_fit_p() implementation using generic function
# it is not working yet, as the new implementation behaves differently and converges poortly
def x_new_broken_run_mcmc_initial_fit_p(
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
    alpha0, alpha1, t0, d, Tau, p = start_vals
    # Note: t0 here is the time in Lightcurve (e.g., BTJD), not normalized phase.
    start_vals_dict = dict(alpha0=alpha0, alpha1=alpha1, t0=t0, d=d, Tau=Tau, p=p)
    return run_mcmc_initial_fit_of_model(
        log_probability_p,
        coshgauss_model_fit_p,
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
