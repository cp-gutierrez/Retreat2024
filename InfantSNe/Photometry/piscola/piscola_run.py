import piscola
import pdb
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table


version = piscola.__version__
print(f'PISCOLA version: v{version}')

#pdb.set_trace()
sn_name='2020jgl'

#lc = Table.read('/home/lara/ICE/Photometry/piscola/new/'+sn_name+'.dat', format='csv',comment='#',delimiter=' ')

sn=piscola.call_sn('/home/lara/ICE/Photometry/piscola/new/'+sn_name+'_atlas.dat')

#Filters, light curves and SED found in:
#filters=sn.filters
#lcs=sn.lcs
#sed=sn.sed

sn.set_sed_template('conley09f')

print(sn)
print(sn.sed)
print(f'Observed bands: { sn.bands}')

#Create .txt file to save parameters of the lc of each band
results=open('/home/lara/ICE/Photometry/gaussian_processes/'+sn_name+'_params_ATLAS.txt','w+')
results.write('band mmax emmax tmax etmax dm15 edm15 \n')


#lcs=sn.lcs
#lcs.get_lc_params()

#sn.rest_lcs_fits.ZTF_g.get_max()

#To plot light curves
#sn.plot_lcs()

#To perform GP and plot the light curves with the respective fit
sn.fit()

sn.plot_fits()

#Fitted light curves --> to estimate tmax
lcs_fit=sn.lc_fits

#Obtain estimated parameters from the fit and not from the observed light curve
for band in lcs_fit.bands:
    lc=lcs_fit[band]
    #pdb.set_trace()
    lc.get_max()
    lc.get_dm15()
    
    print('####'+band+'####')
    print(lc.tmax, lc.tmax_err)
    
    
    #pdb.set_trace()
    if (np.isnan(lc.mmax)==False)&(np.isnan(lc.tmax)==False):
        results.write(band+' '+str(lc.mmax)+' '+str(lc.mmax_err)+' '+str(lc.tmax)+' '+str(lc.tmax_err)+' '+str(lc.dm15)+' '+str(lc.dm15_err)+'\n')
    else:
        print(band)
        results.write(band+' nan nan nan nan nan nan\n')
    
    #fit=piscola.gaussian_process.gp_lc_fit(lc.time,lc.flux)
    #fig,ax=plt.subplots(1,1,figsize=(10,8))
    #mask=(fit[0]>=59070)&(fit[0]<=59170)
    #plt.plot(fit[0][mask],fit[1][mask])
    #plt.show()
    
results.write('\n'+str(sn.lc_parameters))
results.close()

#print(sn.lcs)
#print(sn.lcs.Megacam_g)
#print(sn.lcs.Megacam_g.__dict__)


#sn.rest_lcs_fit()


#sn.export_fits()


#x_pred,y_interp,sigma_interp=piscola.gaussian_process.gp_lc_fit(lc.time,lc.flux,yerr=lc.flux_err) #to obtain fit in order to estimate t0

#pdb.set_trace()

print(sn.lc_parameters)


