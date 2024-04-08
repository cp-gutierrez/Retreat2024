import os, sys
import numpy as np
import operator
from astropy.io import fits
from astropy.io import ascii
import matplotlib.pyplot as plt
import pdb
import numpy.ma
import astropy
from astropy.table import Table
import pandas as pd
import linecache
import math

def mag_to_flux(mag,emag):
    #AB system
    flux_ab=10**(0.4*(8.90-mag))
    eflux_ab=np.abs(-0.4*np.log(10)*flux_ab*emag)
    
    #in microJansky
    flux=flux_ab*10**6
    eflux=eflux_ab*10**6
    
    return flux,eflux

sn_name='2020jgl'
    
#Open .dat file where we will write the format read by Piscola
table=open('/home/lara/ICE/Photometry/piscola/new/'+sn_name+'_delete.dat','w+')
table.write('name z ra dec \n')

dataset=Table.read('/home/lara/ICE/Photometry/piscola/dataset.csv',format='csv')

name_all,ra_all,dec_all,z_all=list(dataset['name']),list(dataset['ra']),list(dataset['dec']),list(dataset['redshift'])
ind=name_all.index(sn_name)
ra,dec,z=ra_all[ind],dec_all[ind],z_all[ind]
table.write(sn_name+' '+str(z)+' '+str(ra)+' '+str(dec)+'\n')

table.write('time flux flux_err zp band mag_sys \n')

lc = pd.read_csv('/home/lara/ICE/Photometry/'+sn_name+'.csv', comment='#')

#For JGL
flux,eflux=mag_to_flux(lc.m,lc.em)

for i in range(len(lc.f)):
    if (lc.det[i]==True)&(lc.band[i]!='ztfi'):#&(lc.mjd[i]<=59100):
        if lc.band[i]=='o':
            band='ATLAS_o'
        elif lc.band[i]=='c': 
            band='ATLAS_c'
        elif lc.band[i]=='ztfr': 
            band='ZTF_r'
        else: 
            band='ZTF_g'
 
        #table.write(str(lc.mjd[i])+' '+str(lc.f_counts[i])+' '+str(lc.ef_counts[i])+' '+str(lc.zp[i])+' '+band+' '+'AB '+'\n')
        table.write(str(lc.mjd[i])+' '+str(flux[i])+' '+str(eflux[i])+' '+str(lc.zp[i])+' '+band+' '+'AB '+'\n')
'''        
####### for SN 2020jgl

lc2 = Table.read('/home/lara/ICE/Photometry/OTHER/sn'+sn_name+'.txt', format='ascii')

mag2,emag2=np.array(lc2['mag']),np.array(lc2['emag'])
flux2,eflux2=mag_to_flux(mag2,emag2)

#pdb.set_trace()

for i in range(len(lc2)):
    if lc2['band'][i]=='U':
        band2='Standard_U'
        mag_sys='BD17'
    elif lc2['band'][i]=='B':
        band2='Standard_B'
        mag_sys='BD17'
    elif lc2['band'][i]=='V':
        band2='Standard_V'
        mag_sys='BD17'
    elif lc2['band'][i]=='g':
        band2='sdss_g'
        mag_sys='AB'
    elif lc2['band'][i]=='r':
        band2='sdss_r'
        mag_sys='AB'
    else:
        band2='sdss_i'
        mag_sys='AB'
        
    table.write(str(lc2['mjd'][i])+' '+str(flux2[i])+' '+str(eflux2[i])+' '+str(8.90)+' '+band2+' '+mag_sys+' \n')
'''
  
table.close()

