# ETV functions

import csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator
from scipy.interpolate import interp1d
from os.path import basename, exists
import sys  
# sys.path.insert(1, '/Users/neisner/Documents/code/utils/')
import filters
import norm
import emcee



def phase_data(time, t0, period):

	'''
	Calualate the phase given a t0 and period.
	'''
	return np.array([-0.5+( ( t - t0-0.5*period) % period) / period for t in time])


def sort_on_x(x,y):
    zipped = list(zip(x,y))
    zipped.sort(key=lambda x:x[0])
    x,y = list(zip(*zipped))
    return x,y


def run_binning(x,y,yerr=None,nbins=100,alias=False):
    # Binning function -- takes into account aliasing and error
    # propogation on errorbins
    bwidth = 1./nbins
    if alias==True:
        phStart,phStop = -0.6,0.6
    else:
        phStart,phStop = -0.5,0.5

    bins      = np.arange(phStart,phStop+bwidth,bwidth)
    bin_means = ( np.histogram(x,bins,weights=y)[0] / np.histogram(x,bins)[0] )
    if yerr is not None:
        bin_errs = ( np.histogram(x,bins,weights=yerr)[0] / np.histogram(x,bins)[0] )
    else:
        bin_errs = None

    return bwidth,bins,bin_means,bin_errs



def run_binning_phased(x,y,yerr=None,nbins=100):
    # Binning function -- takes into account aliasing and error
    # propogation on errorbins
    
    bwidth = 1./nbins
    phStart,phStop = np.nanmin(x), np.nanmax(x)

    bins      = np.arange(phStart,phStop+bwidth,bwidth)
    bin_means = ( np.histogram(x,bins,weights=y)[0] / np.histogram(x,bins)[0] )
    if yerr is not None:
        bin_errs = ( np.histogram(x,bins,weights=yerr)[0] / np.histogram(x,bins)[0] )
    else:
        bin_errs = None
    
    bins =  bins[:-1] + 0.5*bwidth
    
    return bins,bin_means,bin_errs



def time_to_ph(time, period=1., t0=0., pshift=0.):
    '''
    converts time to phase from input ephemeris
    DOES NOT ACCOUNT FOR BARYCENTRIC OR HELIOCENTRIC CORRECTION

    input: time (float or array)  --> time point or array
    input: period (float)
    input: t0 (float)
    input: pshift (float) --> phase shift
    ------
    output: phase (float or array)
    '''
    time = np.array(time)
    ph = np.mod((time-t0)/period, 1.0) + pshift
    ph[ph < 0.0] += 1.
    ph[ph > 0.5] -= 1.

    return ph


def interpolate_signal(x_indep,x_model,y_model,nbins):

    
    bin_width,bin_edges,bin_means,_ = run_binning(x_model,y_model,nbins=nbins)
    
    mask = np.isfinite(bin_means) 

    bin_means = bin_means[mask]
    
    bins = bin_edges[:-1]+0.5*bin_width
    
    bins = bins[mask]
    interpolated_model = PchipInterpolator(bins,bin_means,extrapolate=False)
    interpolated_y = interpolated_model(x_indep)

    return interpolated_y,interpolated_model


def model_signal(period,t0,x,y,x1=None,y1=None,nbins=100):

    ph_x = time_to_ph(x,period,t0)
    ph_x_sorted,y_sorted = sort_on_x(ph_x,y)
    y_interpolated, model = interpolate_signal(ph_x,ph_x_sorted,y_sorted,nbins)

    y_model_subtracted = y-y_interpolated

    return y_model_subtracted, y_interpolated, model


def trend_removal_interact(period_a, t0_a, data, **kwargs):

    figsize = kwargs.get('figsize', (7,8))

    times_original = data['time']
    flux_original = data['flux']
    time_cut = data['time']
    flux_cut = data['flux']
    
    flux_cut_sub_period_a, flux_cut_period_a_model,\
    flux_cut_period_a_function = model_signal(period_a, t0_a,time_cut,flux_cut)
    
    ph_original = time_to_ph(times_original,period_a,t0_a)
    ph_cut = time_to_ph(time_cut,period_a,t0_a)
    
    flux_original_period_a_model = flux_cut_period_a_function(ph_original)
    flux_original_sub_period_a = flux_original - np.array(flux_original_period_a_model)
    
    fig,axes = plt.subplots(2,1, figsize = figsize, sharex = True)
    fig.subplots_adjust(bottom=0.2)

    axes[0].plot(ph_original,flux_original,marker = '.', lw = 0, color = 'grey')
    axes[0].plot(ph_cut,flux_cut,marker = '.', lw = 0, color = 'k')
    
    axes[0].plot(ph_original,flux_original_period_a_model,marker = '.', lw = 0, color = 'darkorange', ms = 2)
    
    
    axes[1].plot(ph_original,flux_original_sub_period_a,'k.', alpha = 0.5)
    
    axes[0].set_ylabel(r'$Flux$',fontsize=18)
    axes[1].set_ylabel(r'$Flux$',fontsize=18)
    axes[1].set_xlabel(r'$\Phi$',fontsize=18)
    
    
    axes[0].tick_params(direction='in', length = 3, which ='minor', colors='grey', labelsize=13)
    axes[0].tick_params(axis="y",direction="inout") #, pad= -20)
    axes[0].tick_params(axis="x",direction="inout") #, pad= -17)   
    axes[0].tick_params(axis='both', length = 5, left='on', top='on', right='on', bottom='on')

    axes[1].tick_params(direction='in', length = 3, which ='minor', colors='grey', labelsize=13)
    axes[1].tick_params(axis="y",direction="inout") #, pad= -20)
    axes[1].tick_params(axis="x",direction="inout") #, pad= -17)   
    axes[1].tick_params(axis='both', length = 5, left='on', top='on', right='on', bottom='on')


    #plt.savefig('/Users/Nora/Documents/research/projects/fluffy/figs/sig_removal.png', dpi = 300)
    plt.show()  
    
    data['flux_sub_binary'] = flux_original_sub_period_a
    data['binary_model'] = flux_original_period_a_model

    return data


def log_prior(theta):
    alpha0, alpha1, t0, d, Tau = theta
    
    if  (0 < alpha0 < 10) and ( -10 < alpha1 < 0) and (-0.5 < t0 < 0.5) and (0 < d < 10) and (0.5 < Tau < 50):
        return 0.0
    return -np.inf


def log_likelihood(theta, x, y, yerr):
    
    alpha0, alpha1, t0, d, Tau = theta
    
    cosh_term = np.cosh( (x - t0) /d )
    exp_term = np.exp(1 - cosh_term)
    pow_term = pow((1 - exp_term), Tau)
    
    psi = 1 - pow_term
    
    model = alpha0 + (alpha1 * psi)
    
    return -0.5 * np.sum((y - model) ** 2 / (yerr ** 2) )


def log_probability(theta, x, y, yerr):
    
    # check that the priors are satisfied
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, x, y, yerr)


def get_starting_positions(start_vals, nwalkers=128):

    p0 = np.array( [ [ val+1e-5*np.random.randn()
                       for jj,val in enumerate(start_vals) ]
                      for i in range(nwalkers) ] )
    
    yield p0



def detrend_suz(data, T_dur = 10, plot = True, **kwargs):

    # remove the overall trend and then re do the fitting to remoe the signal to see whether we can improve on the trend removal

     # define the eclipse duration 

    figsize = kwargs.get('figsize', (8,4))
    
    nmed = int(30*3*T_dur)
    nmed = 2*int(nmed/2)+1

    ff = filters.NIF(np.array(data.flux_sub_binary),nmed,10,fill=True,verbose=False)
    # first number is three times transit durations, the second quite small (10,20 )
    l = np.isfinite(ff)
    g = interp1d(data.time[l],ff[l],bounds_error=False,fill_value=np.nan)
    ff = g(data.time)

    data['flux_detrended'] = data.flux - (ff)
    
    if plot == True:
        fig, ax = plt.subplots(figsize=figsize)

        # plot the original data
        plt.plot(data.time,data.flux,'.', color = 'k', ms = 1) 

        # plot the data minus the binary signal
        plt.plot(data.time,data.flux_sub_binary+1,'.')

        # plot the model that was fit to the binary model
        plt.plot(data.time,((ff*1.2)+1),'.', color = 'orange')

        # make a second figure with the detrended data
        fig, ax = plt.subplots(figsize=figsize)

        #plt.plot(alltime, allflux,'.k', ms = 1)
        plt.plot(data.time, data.flux_detrended,'.r', alpha = 1, ms = 1)

    return data


def coshgauss_model_fit(x, alpha0, alpha1, t0, d, Tau):
    """
    Calculates the coshgauss model fit for the given parameters.

    Parameters:
    x (float): The input value.
    alpha0 (float): The coefficient for the constant term.
    alpha1 (float): The coefficient for the psi term.
    t0 (float): The center of the cosh term.
    d (float): The width of the cosh term.
    Tau (float): The exponent for the pow term.

    Returns:
    float: The model fit value.
    """
    cosh_term = np.cosh((x - t0) / d)
    exp_term = np.exp(1 - cosh_term)
    pow_term = pow((1 - exp_term), Tau)
    
    psi = 1 - pow_term
    
    model = alpha0 + (alpha1 * psi)
    
    return model


def plot_initial_guess(data, ph_binned, flux_binned, err_binned, *start_vals, **kwargs):
    """
    Plots the initial guess for the data using the given parameters.

    Parameters:
    data (DataFrame): The data to be plotted.
    ph_binned (array-like): The binned phase values.
    flux_binned (array-like): The binned flux values.
    err_binned (array-like): The binned error values.
    *start_vals: Variable length argument list of start values.
    **kwargs: Additional keyword arguments.

    Returns:
    None
    """
    figsize = kwargs.get('figsize', (8,4))
    
    fig = plt.subplots(figsize=figsize, sharex=True)
    plt.errorbar(ph_binned, flux_binned, yerr=err_binned, fmt=".k", capsize=0, zorder = 2)
    plt.scatter(data.phase, data.flux, zorder = -2)
    plt.plot(data.phase, coshgauss_model_fit(data.phase, *start_vals), lw = 0, marker = '.', markersize = 0.5, alpha=1, zorder = 2, color= 'red')


# for the mcmc 

def run_mcmc_initial_fit(data, start_vals, nruns = 1000, plot_chains = False, plot = True, **kwargs):

    figsize = kwargs.get('figsize', (8,4))
    
    pos = list(get_starting_positions(start_vals,nwalkers=128))[0]

    nwalkers = 128
    ndim = len(start_vals)

    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=(data.phase, data.flux, data.err))

    sampler.run_mcmc(pos, nruns, progress=True, store = True)

    tau = sampler.get_autocorr_time(tol=0)

    samples = sampler.get_chain()
    labels = ["alpha0", "alpha1", "t0", "d", "Tau"]
    
    if plot_chains == True:
        
        fig, axes = plt.subplots(5, figsize=figsize, sharex=True)
        
        for i in range(ndim):
            ax = axes[i]
            ax.plot(samples[:, :, i], "k", alpha=0.3)
            ax.set_xlim(0, len(samples))
            ax.set_ylabel(labels[i])
            ax.yaxis.set_label_coords(-0.1, 0.5)

        axes[-1].set_xlabel("step number");
        plt.show()
        
    flat_samples = sampler.get_chain(discard=600, thin=15, flat=True)
    
    mean_alpha0 = np.median(flat_samples[:,0])
    mean_alpha1 = np.median(flat_samples[:,1])
    mean_t0 = np.median(flat_samples[:,2])
    mean_d = np.median(flat_samples[:,3])
    mean_Tau = np.median(flat_samples[:,4])

    inds = np.random.randint(len(flat_samples), size=1000)

    if plot == True:
        fig, axes = plt.subplots(figsize=figsize, sharex=True)
        
        for ind in inds:
            sample = flat_samples[ind]
            plt.plot(data.phase, coshgauss_model_fit(data.phase, sample[0], sample[1], sample[2], sample[3], sample[4]), "C1", lw = 0, marker = '.',alpha=0.1, zorder = 2)

        plt.errorbar(data.phase, data.flux, yerr=data.err, fmt=".k", capsize=0, zorder = -2)

        plt.plot(data.phase, coshgauss_model_fit(data.phase, mean_alpha0, mean_alpha1,mean_t0, mean_d, mean_Tau), lw = 0, marker = '.', markersize = 0.5, alpha=1, zorder = 2, color= 'red')

        plt.xlabel("x")
        plt.ylabel("y");
        plt.show()
        
    return mean_alpha0, mean_alpha1, mean_t0, mean_d, mean_Tau



# now that we have the best fit model, we fit this model to each individual eclipse using mcmc
# the free parameters are t0, alpha0 and alpha1, the rest are fixed by phase folded model

# May have to change the priors!
def log_prior_fitting(theta):
    alpha0, alpha1, t0 = theta
    
    if  (0.95 < alpha0 < 1.1) and ( -0.16 < alpha1 < 0) and (-0.5 < t0 < 0.5):
        return 0.0
    
    return -np.inf


def log_likelihood_fitting(theta, x, y, yerr, mean_d, mean_Tau):
    
    alpha0, alpha1, t0 = theta
        
    cosh_term = np.cosh( (x - t0) /mean_d )
    exp_term = np.exp(1 - cosh_term)
    pow_term = pow((1 - exp_term), mean_Tau)
    
    psi = 1 - pow_term
    
    model = alpha0 + (alpha1 * psi)
    
    return -0.5 * np.sum((y - model) ** 2 / (yerr ** 2) )

def log_probability_fitting(theta, x, y, yerr, mean_d, mean_Tau):
    
    # check that the priors are satisfied
    lp = log_prior_fitting(theta)
    if not np.isfinite(lp):
        return -np.inf
    lnp = lp + log_likelihood_fitting(theta, x, y, yerr, mean_d, mean_Tau)

    return lnp


def get_starting_positions_fitting(start_vals, nwalkers=128):

    p0 = np.array( [ [ val+1e-5*np.random.randn()
                       for jj,val in enumerate(start_vals) ]
                      for i in range(nwalkers) ] )
    
    yield p0

    
    
def fit_each_eclipse(data, n_transits, t0, period, mean_alpha0, mean_alpha1, mean_t0, mean_d, mean_Tau, outfile_path, min_number_data = 20):
    
    if exists('{}'.format(outfile_path)):
        print("Existing manifest file found, will skip previously processed LCs and append to end of manifest file")
        sys.stdout.flush()
    else:
        print("Creating new manifest file")
        sys.stdout.flush()
        metadata_header = ["number", 'epoch', 't0', 'stdv_t0', 'alpha1', 'alpha2', 'd', 'Tau']
        with open('{}'.format(outfile_path), 'w') as f: # save in the photometry folder
            writer = csv.writer(f, delimiter=',')
            writer.writerow(metadata_header)
    
    manifest_table = pd.read_csv('{}'.format(outfile_path))
    number_done = manifest_table['number']

    tr_index = range(0,n_transits)
    
    for i in tr_index:
        if not np.isin(i,number_done):
            
            transit_time = t0+(period*i)

            x = np.array(data.time)
            y = np.array(data.flux) 
            yerr = np.array(data.err)
            
            mask = (x > (transit_time - (0.2*period))) & (x < (transit_time + (0.2*period))) 
            
            x = np.array(x[mask])
            y = np.array(y[mask]) 
            yerr = np.array(yerr[mask])

            # convert x from time (JD-like) time to normalized phase, with t0 as midpoint (0)
            x = np.array([-0.5+( ( t - t0-0.5*period) % period) / period for t in x])
            
            if len(x) > min_number_data:

                print (transit_time, mean_alpha0, mean_alpha1, mean_t0)
                
                start_vals = [mean_alpha0, mean_alpha1, 0]  # use 0 instead of mean_t0, as x is shifted to be 0 for t0
                
                pos = list(get_starting_positions(start_vals,nwalkers=64))[0]
                
                nwalkers = 64
                ndim = len(start_vals)
                
                # start the mcmc fitting
                sampler2 = emcee.EnsembleSampler(nwalkers, ndim, log_probability_fitting, args=(x, y, yerr, mean_d, mean_Tau))
                
                sampler2.run_mcmc(pos,10000, progress=True)
        
                flat_samples2 = sampler2.get_chain(discard=400, thin=15, flat=True)
                
                mean_alpha0_fit = np.nanmedian(flat_samples2[:,0])
                mean_alpha1_fit = np.nanmedian(flat_samples2[:,1])
                mean_t0_fit = np.nanmedian(flat_samples2[:,2])
                stdv_t0_fit = np.nanstd(flat_samples2[:,2])
                
                fig = plt.subplots(figsize=(10, 3), sharex=True)
                
                plt.errorbar(x, y, yerr=yerr, fmt=".k", capsize=0, zorder = -2)
                plt.plot(x, coshgauss_model_fit(x, mean_alpha0_fit, mean_alpha1_fit,mean_t0_fit, mean_d, mean_Tau), lw = 1, marker = '.', markersize = 0.5, alpha=1, zorder = 2, color= 'red')
        
                plt.show()
                
                with open('{}'.format(outfile_path), 'a') as f: # save in the photometry folder
                    writer = csv.writer(f, delimiter=',')
                    writer.writerow([i, transit_time, mean_t0_fit, stdv_t0_fit, mean_alpha0_fit, mean_alpha1_fit, mean_d, mean_Tau])
            else:  
                continue
        else:
            print ("Number {} has already been completed -- skip".format(i))
        