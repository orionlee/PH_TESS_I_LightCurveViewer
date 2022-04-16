# Input file for {alias}
# template: {template_type}

# Light curve data file
fname_tr = ["{fname_tr}"]

# MCMC controls
thin_factor = {mcmc_thin_factor}
niter = {mcmc_niter}
nchains = {mcmc_nchains}

# method to use 'mcmc' sample the posteriors, 'plot' creates the plots using posteriors from a previous run
method = "mcmc"

# We want the planet values in Earth units
unit_mass = "earth"

# Stellar parameters, typically from MAST / Gaia. See ExoFOP parameters:
# https://exofop.ipac.caltech.edu/tess/target.php?id={tic}
mstar_mean = {mass}
mstar_sigma = {e_mass}
rstar_mean = {rad}
rstar_sigma = {e_rad}
tstar_mean = {Teff}
tstar_sigma = {e_Teff}

# Since we want to fit transit data, we need to include the next line
fit_tr = [True]

# THIS IS THE MOST IMPORTANT LINE TO ADD IN A SINGLE TRANSIT MODEL
# If set to `True`, we tell pyaneti that we will be fitting only one transit
# This will tell to the code how to deal with the period and the dummy semi-amplitude
# see: https://github.com/oscaribv/pyaneti/blob/master/inpy/example_single/input_fit.py
# - note: single transit mode does not appear to be compatible with
#   rho fitting mode (`sample_stellar_density = True`). Don't use them together.
is_single_transit = {is_single_transit}


# Specify which kind of priors we want
# For a transit fit of a single planet with eccentric orbit
# we fit only for time of minimum conjunction T0, period P, impact parameter b, scaled semi-major axis (a)
# scaled planet radius rp, and limb darkening coefficients q1 and q2
# We let MCMC pick eccentricity e, and angle of periastron w (with uniform priors on full range)
fit_t0 = ["{type_epoch}"]
fit_P = ["{type_period}"]

# use ew1 and ew2, polar form of e and w, to parametrize them.
# they'd be better sampled for small eccentricities
# see: https://github.com/oscaribv/pyaneti/wiki/Parametrizations#limb-darkening-coefficients
# for case both ew1 and ew2 are both uniform priors, they are equivalent to uniform priors for e
is_ew = True
fit_ew1 = ["{type_ew}"]
fit_ew2 = ["{type_ew}"]

fit_b = ["u"]  # TODO:
# sampling stellar density: if True
# - fit_a would fit rho* instead of a/R* (aka matched_rho),
# - an useful constraint if there are reliable, independent rho* available
# - see:  https://github.com/oscaribv/pyaneti/wiki/Parametrizations#stellar-density
sample_stellar_density = {sample_stellar_density}
fit_a = ["{type_a}"]  # {comment_a}
fit_rp = ["{type_rp}"]  # Rp/R*
fit_q1 = "{type_q1}"
fit_q2 = "{type_q2}"

# Set the prior ranges
# here we use a rough range based on the estimated T0 , P

# T0 prior type: {type_epoch}
min_t0 = [{val1_epoch}]  # {time_format}
max_t0 = [{val2_epoch}]  # {time_format}

# P prior type: {type_period}
min_P = [{val1_period}]  # days
max_P = [{val2_period}]  # days


# ew prior type: {type_ew}
min_ew1 = [{val1_ew1}]
max_ew1 = [{val2_ew1}]
min_ew2 = [{val1_ew2}]
max_ew2 = [{val2_ew2}]

# Minimum and maximum limits for b
min_b = [0.0]
max_b = [1.0]
# {comment_a} prior type: {type_a}
min_a = [{val1_a}]
max_a = [{val2_a}]
# Rp/R* prior type: {type_rp}
min_rp = [{val1_rp}]
max_rp = [{val2_rp}]
# q1 prior type: {type_q1}
min_q1 = {val1_q1}
max_q1 = {val2_q1}
# q2 prior type: {type_q2}
min_q2 = {val1_q2}
max_q2 = {val2_q2}

#
# Plot Controls
#

# If True it creates a correlation plot
is_plot_correlations = True

plot_binned_data = True

tr_xlabel = "{lc_time_label}"


# End of the input file for {alias}
