import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import scipy.stats
from ultranest import ReactiveNestedSampler
from scipy.special import lambertw
import os
from scipy.optimize import curve_fit
from ultranest import ReactiveNestedSampler
import scipy.stats
import numpy
import warnings
warnings.filterwarnings("ignore")


# Load in the data, get rid of NaNs
lc_data = pd.read_csv('/Users/lalvopis/SNe_ICE/Photometry/geza_typeII/2020nny.csv', sep = ',')
cnd = np.isnan(lc_data['f'])

lc_data = lc_data[~cnd]

# Separate the data into a list according to separate bands
bands = ['c', 'o', 'ztfg', 'ztfr']

lcs = []
for i in range (len(bands)):
    lcs.append(lc_data[lc_data['band'] == bands[i]])

# Cut the relevant parts. Since the fitting function is simple, we can only use the initial rise, the first peak/plateu and the non-detections
cut_lcs = []

c_lim = [58950, 59035]
o_lim = [58950, 59040]
g_lim = [58980, 59040]
r_lim = [58980, 59038]

for i in range (len(lcs)):
    if lcs[i]['band'].values[0] == 'c':
        lim_cond = (lcs[i]['mjd'] > c_lim[0]) & (lcs[i]['mjd'] < c_lim[1])
    elif lcs[i]['band'].values[0] == 'o':
        lim_cond = (lcs[i]['mjd'] > o_lim[0]) & (lcs[i]['mjd'] < o_lim[1])
    elif lcs[i]['band'].values[0] == 'ztfg':
        lim_cond = (lcs[i]['mjd'] > g_lim[0]) & (lcs[i]['mjd'] < g_lim[1])
    else:
        lim_cond = (lcs[i]['mjd'] > r_lim[0]) & (lcs[i]['mjd'] < r_lim[1])
        
    cut_lcs.append(lcs[i][lim_cond])
    
# Separate actual detections and non-detections by putting in a limit. This is a very by-hand step; I just made this to follow the prescription of the past analysis
lcs_nd = []
photo_vls = []
for i in range (len(cut_lcs)):
    if len(cut_lcs[i]) > 0:
        cond = cut_lcs[i]['mjd'] < 59010
        lcs_nd.append(cut_lcs[i][cond].values[:,:3].astype(float))
        photo_vls.append(cut_lcs[i][~cond].values[:,:3].astype(float))
    else:
        lcs_nd.append([])
        photo_vls.append([])

# The initial t0 guess for 2020nny in MJD
t0_guess = 59000 

# Normalize the flux values, subtract initial t0 guess from times.
for i in range (len(photo_vls)):
    if len(photo_vls[i]) > 0:
        norm_fact = max(photo_vls[i][:,1])
        lcs_nd[i][:,2] /= norm_fact
        lcs_nd[i][:,1] /= norm_fact
        lcs_nd[i][:,0] -= t0_guess
        
        photo_vls[i][:,2] /= norm_fact
        photo_vls[i][:,1] /= norm_fact
        photo_vls[i][:,0] -= t0_guess
        

# Fitting procedure using ultranest

# For ultranest one has to set up the parameters and the priors in a reasonable range
parameters = ['t0', 'fmc', 'fmo', 'fmg', 'fmr', 'tec', 'teo', 'teg', 'ter' ,'f']

def prior_transform(cube):
    # the argument, cube, consists of values from 0 to 1
    # we have to convert them to physical scales

    params = cube.copy()
    params[0] = cube[0] * 10 + 5      # parameter for t0; t0 guess has already been subtracted from the data
    params[1] = cube[1] * 0.6 + 0.4   # the plateou flux values; since the fluxes are normalized, this should be close to 1
    params[2] = cube[2] * 0.6 + 0.4
    params[3] = cube[3] * 0.6 + 0.4
    params[4] = cube[4] * 0.6 + 0.4
    params[5] = cube[5] * 5           # the charachteristic rise time; this should be on the order of days
    params[6] = cube[6] * 5
    params[7] = cube[7] * 5
    params[8] = cube[8] * 5
    params[9] = cube[9] * 3           # error inflation parameter; in case uncertainties are underestimated
    return params
    

# The fitting function; inverse exponential
def flux_evol(x, fm, t0, te):    
    return fm * (1 - np.exp(-(x - t0)/te))

def flux_evol_mod(x, fm, t0, te, offs):
    cond = x < t0
    result = np.zeros(x.shape)
    
    result[cond] = offs
    result[~cond] = fm * (1 - np.exp(-(x[~cond] - t0)/te)) + offs
    
    return result
    
    
# Likelihood for the Bayesian fitting: this is where the different band light curves get connected through t0 and the rise time condition
def log_likelihood(params):
    
    t0, fmc, fmo, fmg, fmr, tec, teo, teg, ter, f = params
    fms = [fmc, fmo, fmg, fmr]
    tes = [tec, teo, teg, ter]

    # Check for the models in each band if the model light curve goes ABOVE the non-detections; if yes, that is a bad model
    for i in range (len(photo_vls)):   # I have put in here the last g_nondet for simplicity
        try:
            non_dec_cond = lcs_nd[i][-1,1] < flux_evol(lcs_nd[i][-1,0], fm=fms[i], t0=t0, te=tes[i])  # It is enough to check this for the last non-det, if it is deep enough

            if non_dec_cond:
                return -9999
        except (IndexError,AttributeError):  # It can happen that we have no non-detections, then we cannot check the condition, this is here to catch that
            pass
    
    # check if te.c > te.o; the ATLAS C band light curve should rise faster than the O band (redder bands rise slower)
    if tes[0] > tes[1]:
        return -9999
    # same for ztf g and r
    if tes[2] > tes[3]:
        return -9999

    # if any of the above, return -9999, force the sampler away
    
    # if we are still in here, do the actual loglike calculation
    loglike = 0
    for i in range (len(photo_vls)):
        y_model = flux_evol(photo_vls[i][:,0], fm=fms[i], t0=t0, te=tes[i])
    
        sigma2 = photo_vls[i][:,2]**2 + y_model **2 * f**2    # the error inflation comes in here
        loglike += -0.5 * np.sum((y_model - photo_vls[i][:,1])**2 / sigma2 + np.log(2*np.pi*sigma2))

    return loglike

    
# Do the fitting; this code can throw an AssortionError, if so, check if there are NaNs in the data. If no, just run the above def scripts again, and try once more; sometimes it starts from the wrong intial parameters
sampler = ReactiveNestedSampler(parameters, log_likelihood, prior_transform,
    wrapped_params=np.full(len(parameters), False)
)

result = sampler.run(min_num_live_points=400, dKL=np.inf, min_ess=100)  # 400 does the trick, but if you need sharp posteriors, might want to increase it; computationally expensive!



# Plot the results
f = plt.figure(figsize = (12,8))
pdf=PdfPages('/Users/lalvopis/SNe_ICE/Photometry/geza_typeII/2020nny_test.pdf')
post_mean = sampler.results['posterior']['mean']

cls_s = ['blue', 'orange', 'green','red']

cls_line = ['royalblue', 'orange', 'limegreen','coral']


for i in range (len(bands)):
    ax = f.add_subplot(2,2,i+1)
    
    if len(photo_vls[i]) > 0:
        plt.errorbar(photo_vls[i][:,0], photo_vls[i][:,1], np.sqrt(photo_vls[i][:,2]**2 + photo_vls[i][:,1]**2 * post_mean[-1]**2), fmt = '.', color = 'royalblue')
        plt.scatter(photo_vls[i][:,0], photo_vls[i][:,1], color = 'blue', edgecolor = 'black', zorder = 10)
        plt.scatter(lcs_nd[i][:,0], lcs_nd[i][:,1], color = 'red', marker = 'v')

    
        t_grid = np.linspace(-10, 45, 700)

        from ultranest.plot import PredictionBand
        band = PredictionBand(t_grid)

        # go through the solutions
        for t0, fmc, fmo, fmg, fmr, tec, teo, teg, ter, fadd in sampler.results['samples']:
            # compute for each time the y value
            fms = [fmc, fmo, fmg, fmr]
            tes = [tec, teo, teg, ter]
            band.add(flux_evol(t_grid, fm=fms[i], t0=t0, te=tes[i]))

        y_mod = flux_evol(t_grid, fm=post_mean[1+i], t0=post_mean[0], te=post_mean[i+5])
        cond = y_mod > 0
        band.shade(color=cls_line[i], alpha=0.5, label = '1$\sigma$')
        band.shade(q=0.48, color=cls_line[i], alpha=0.3, label = '2$\sigma$')
        plt.plot(t_grid[cond], y_mod[cond], color = 'grey', lw = 1.5, label = 'Median fit')
        
        plt.ylim(-0.2, 1.2)
        
        plt.xlim(-35, 55)
        
        if i > 1:
            plt.xlabel('MJD - ' + str(np.round(t0_guess, 2)), fontsize = 14)

    pdf.savefig(f)
    plt.close(f)
    pdf.close()