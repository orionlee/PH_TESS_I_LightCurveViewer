# Input file for {alias}
# template: {template_filename}

# Light curve data file
fname_tr = ["{fname_tr}"]

# MCMC controls
thin_factor = 1
niter = 500
nchains = 100

# method to use 'mcmc' sample the posteriors, 'plot' creates the plots using posteriors from a previous run
method = "mcmc"

# We want the planet values in Earth units
unit_mass = "earth"

# Stellar parameters from MAST, which is the basis of ExoFOP data at:
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
# We need to tell pyaneti that we will be fitting only one transit
# This will tell to the code how to deal with the period and the dummy semi-amplitude
# see: https://github.com/oscaribv/pyaneti/blob/master/inpy/example_single/input_fit.py
# - note: single transit mode does not appear to be compatible with
#   matching rho mode (`_sample_stellar_density = True`). Don't use them together.
is_single_transit = True

# Specify which kind of priors we want
# For a transit fit of a single planet with eccentric orbit
# we fit only for time of minimum conjunction T0, period P, impact parameter b, scaled semi-major axis (a)
# scaled planet radius rp, and limb darkening coefficients q1 and q2
# We let MCMC pick eccentricity e, and angle of periastron w (with uniform priors on full range)
fit_t0 = ["u"]  # uniform prior for T0
fit_P = ["u"]  # uniform prior for P

# use ew1 and ew2, polar form of e and w, to parametrize them.
# they'd be better sampled for small eccentricities
# see: https://github.com/oscaribv/pyaneti/wiki/Parametrizations#limb-darkening-coefficients
is_ew = True
fit_ew1 = ["u"]  # uniform priors for ew1
fit_ew2 = ["u"]  # uniform priors for ew2; together implying uniform priors for e

fit_b = ["u"]  # uniform prior for b
fit_a = ["u"]  # uniform prior for a/R*
fit_rp = ["u"]  # uniform prior for Rp/R*
fit_q1 = "g"  # gaussian prior for q1
fit_q2 = "g"  # gaussian prior for q2

# Set the prior ranges
# here we use a rough range based on the estimated T0 , P

# Minimum and maximum limits for T0")
min_t0 = [{epoch_min}]  # btjd
max_t0 = [{epoch_max}]  # btjd

# Minimum and maximum limits for P
min_P = [{period_min}]  # days
max_P = [{period_max}]  # days


# sample full range of e/w
min_ew1 = [-1.0]
max_ew1 = [1.0]
min_ew2 = [-1.0]
max_ew2 = [1.0]
# Minimum and maximum limits for b
min_b = [0.0]
max_b = [1.0]
# Minimum and maximum limits for a/R*
min_a = [{a_min}]
max_a = [{a_max}]
# Minimum and maximum limits for Rp/R*
min_rp = [{r_planet_in_r_star_min}]
max_rp = [{r_planet_in_r_star_max}]
# Gaussian priors on q1
min_q1 = {q1}  # Mean
max_q1 = {e_q1}  # Standard deviation
# Gaussian priors on q2
min_q2 = {q2}  # Mean
max_q2 = {e_q2}  # Standard deviation

#
# Plot Controls
#

# If True it creates a correlation plot
is_plot_correlations = True

plot_binned_data = True

tr_xlabel = "{lc_time_label}"


# End of the input file for {alias}
