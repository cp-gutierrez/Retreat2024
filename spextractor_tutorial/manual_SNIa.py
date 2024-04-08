#!/usr/bin/env python
# coding: utf-8

import matplotlib.pyplot as plt
from spextractor import Spextractor
import pandas as pd
import numpy as np
from astropy.table import Table
import pandas as pd
from astropy.io import ascii
import glob


####### read spectrum 
 
spectra=glob.glob('spectra/SN2023vjh/*')
spectra.sort()

i=0

print(spectra[i])
f = open(spectra[i], 'r')



# Name of the supernova
string = spectra[i]
SName=string.split('/')[-1].split('.')[0]
print(SName)


# spetra to a numpy array

tab = pd.read_csv(spectra[i], comment='#',header=13,delim_whitespace=True,names=('x','f'))
#,dtype={'x': np.float64, 'f': np.float64, 'e': np.float64})
#if type(tab['e']) != float: tab['e']=tab['f']*0. + 1.25*np.std(tab['f'])
tab2=tab.to_numpy()
tab2

# Name of the supernova

string = spectra[i]
SName=string.split('/')[-1].split('.')[0]
SName


#spextractor
sn_type='Ia'
z=None
spex = Spextractor(tab2, sn_type,z=z, manual_range=True, plot=True)

spex.create_model(downsampling=2.0,sigma_outliers=8.5,model_uncertainty=True,optimize_noise=False)
spex.process()



# CA II H&K

ca = 'Ca II H&K'
vca = spex.vel[ca]
vca_err = spex.vel_err[ca]
print(f'vca = {vca:.3f} +- {vca_err:.3f}')

pca = spex.pew[ca]
pca_err = spex.pew_err[ca]
print(f'pca = {pca:.3f} +- {pca_err:.3f}')

try: 
    dca = spex.depth[ca][0]
except:
    dca = np.nan
    
try: 
    dca_err = spex.depth_err[ca][0]
except:
    dca_err = np.nan
    
print(f'dca = {dca:.3f} +- {dca_err:.3f}')

# SI II 4000A

si4 = 'Si II 4000A'
vsi4 = spex.vel[si4]
vsi4_err = spex.vel_err[si4]
print(f'vsi4 = {vsi4:.3f} +- {vsi4_err:.3f}')

psi4 = spex.pew[si4]
psi4_err = spex.pew_err[si4]
print(f'psi4 = {psi4:.3f} +- {psi4_err:.3f}')

try:
    dsi4 = spex.depth[si4][0]
except: 
    dsi4 = np.nan

try: 
    dsi4_err = spex.depth_err[si4][0]
except:
    dsi4_err = np.nan
    
print(f'dsi4 = {dsi4:.3f} +- {dsi4_err:.3f}')

# MG II 4300A

mg = 'Mg II 4300A'
vmg = spex.vel[mg]
vmg_err = spex.vel_err[mg]
print(f'vmg = {vmg:.3f} +- {vmg_err:.3f}')

pmg = spex.pew[mg]
pmg_err = spex.pew_err[mg]
print(f'pmg = {pmg:.3f} +- {pmg_err:.3f}')

try: 
    dmg = spex.depth[mg][0]
except:
    dmg = np.nan
    
try: 
    dmg_err = spex.depth_err[mg][0]
except:
    dmg_err = np.nan
    
print(f'dmg = {dmg:.3f} +- {dmg_err:.3f}')

# FE II 4800A

fe = 'Fe II 4800A'
vfe = spex.vel[fe]
vfe_err = spex.vel_err[fe]
print(f'vfe = {vfe:.3f} +- {vfe_err:.3f}')

pfe = spex.pew[fe]
pfe_err = spex.pew_err[fe]
print(f'pfe = {pfe:.3f} +- {pfe_err:.3f}')

try: 
    dfe = spex.depth[fe][0]
except:
    dfe = np.nan
    
try: 
    dfe_err = spex.depth_err[fe][0]
except:
    dfe_err = np.nan
    
print(f'dfe = {dfe:.3f} +- {dfe_err:.3f}')

# S II 5500A

s = 'S II 5500A'
vs = spex.vel[s]
vs_err = spex.vel_err[s]
print(f'vs = {vs:.3f} +- {vs_err:.3f}')

ps = spex.pew[s]
ps_err = spex.pew_err[s]
print(f'ps = {ps:.3f} +- {ps_err:.3f}')

try: 
    ds = spex.depth[s][0]
except:
    ds = np.nan
    
try: 
    ds_err = spex.depth_err[s][0]
except:
    ds_err = np.nan

print(f'ds = {ds:.3f} +- {ds_err:.3f}')

# SI II 5800A

si5 = 'Si II 5800A'
vsi5 = spex.vel[si5]
vsi5_err = spex.vel_err[si5]
print(f'vsi = {vsi5:.3f} +- {vsi5_err:.3f}')

psi5 = spex.pew[si5]
psi5_err = spex.pew_err[si5]
print(f'psi5 = {psi5:.3f} +- {psi5_err:.3f}')

try: 
    dsi5 = spex.depth[si5][0]
except:
    dsi5 = np.nan
    
try: 
    dsi5_err = spex.depth_err[si5][0]
except:
    dsi5_err = np.nan

print(f'dsi5 = {dsi5:.3f} +- {dsi5_err:.3f}')

# SI II 6150A

si = 'Si II 6150A'
vsi = spex.vel[si]
vsi_err = spex.vel_err[si]
print(f'vsi = {vsi:.3f} +- {vsi_err:.3f}')

psi = spex.pew[si]
psi_err = spex.pew_err[si]
print(f'psi = {psi:.3f} +- {psi_err:.3f}')

try: 
    dsi = spex.depth[si][0]
except:
    dsi = np.nan
    
try: 
    dsi_err = spex.depth_err[si][0]
except:
    dsi_err = np.nan

print(f'dsi = {dsi:.3f} +- {dsi_err:.3f}')


# O I 7500A

o = 'O I 7500A'
vo = spex.vel[o]
vo_err = spex.vel_err[o]
print(f'vo = {vo:.3f} +- {vo_err:.3f}')

po = spex.pew[o]
po_err = spex.pew_err[o]
print(f'po = {po:.3f} +- {po_err:.3f}')

try: 
    do = spex.depth[o][0]
except:
    do = np.nan
    
try: 
    do_err = spex.depth_err[o][0]
except:
    do_err = np.nan
    
print(f'do = {do:.3f} +- {do_err:.3f}')


# Ca II IR

caIR = 'Ca II 8542A'
vcaIR = spex.vel[caIR]
vcaIR_err = spex.vel_err[caIR]
print(f'vcaIR = {vcaIR:.3f} +- {vcaIR_err:.3f}')

pcaIR = spex.pew[caIR]
pcaIR_err = spex.pew_err[caIR]
print(f'pcaIR = {pcaIR:.3f} +- {pcaIR_err:.3f}')

try: 
    dcaIR = spex.depth[caIR][0]
except:
    dcaIR = np.nan
    
try: 
    dcaIR_err = spex.depth_err[caIR][0]
except:
    dcaIR_err = np.nan
    
print(f'dcaIR = {dcaIR:.3f} +- {dcaIR_err:.3f}')



#Create plot 
fig, ax = spex.plot

ax.set_title(SName)

plt.tight_layout()
fig.savefig('Images/SN2023vjh/'+SName, dpi=300)

plt.close('all')


from tabulate import tabulate

datav = [['%.3f +- %.3f'%(vca,vca_err), '%.3f +- %.3f'%(vsi4, vsi4_err), '%.3f +- %.3f'%(vmg, vmg_err), '%.3f +- %.3f'%(vfe, vfe_err), '%.3f +- %.3f'%(vs, vs_err), '%.3f +- %.3f'%(vsi5, vsi5_err), '%.3f +- %.3f'%(vsi, vsi_err), '%.3f +- %.3f'%(vo, vo_err), '%.3f +- %.3f'%(vcaIR, vcaIR_err)]]
print (tabulate(datav, headers=["Ca II H&K", "Si II 4000A", "Mg II 4300A", "Fe II 4800A", "S II 5500A", "Si II 5800A", "Si II 6150A","O I 7500A","Ca II 8542A"]))

datap = [['%.3f +- %.3f'%(pca, pca_err), '%.3f +- %.3f'%(psi4,psi4_err), '%.3f +- %.3f'%(pmg, pmg_err), '%.3f +- %.3f'%(pfe, pfe_err), '%.3f +- %.3f'%(ps, ps_err), '%.3f +- %.3f'%(psi5,psi5_err), '%.3f +- %.3f'%(psi,psi_err),'%.3f +- %.3f'%(po,po_err),'%.3f +- %.3f'%(pcaIR,pcaIR_err)]]
print (tabulate(datap, headers=["Ca II H&K", "Si II 4000A", "Mg II 4300A", "Fe II 4800A", "S II 5500A", "Si II 5800A", "Si II 6150A","O I 7500A","Ca II 8542A"]))

datad = [['%.3f +- %.3f'%(dca, dca_err), '%.3f +- %.3f'%(dsi4, dsi4_err), '%.3f +- %.3f'%(dmg, dmg_err), '%.3f +- %.3f'%(dfe, dfe_err), '%.3f +- %.3f'%(ds, ds_err), '%.3f +- %.3f'%(dsi5,dsi5_err), '%.3f +- %.3f'%(dsi,dsi_err),'%.3f +- %.3f'%(do,do_err),'%.3f +- %.3f'%(dcaIR,dcaIR_err)]]
print (tabulate(datad, headers=["Ca II H&K", "Si II 4000A", "Mg II 4300A", "Fe II 4800A", "S II 5500A", "Si II 5800A", "Si II 6150A","O I 7500A","Ca II 8542A"]))



f = open('results/Results_SN2023vjh.txt', 'a')
f.write('%10s %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f \n'%(SName, vca, vca_err, pca, pca_err, dca, dca_err, vsi4, vsi4_err, psi4, psi4_err, dsi4, dsi4_err, vmg, vmg_err, pmg, pmg_err, dmg, dmg_err, vfe, vfe_err, pfe, pfe_err, dfe, dfe_err, vs, vs_err, ps, ps_err, ds, ds_err, vsi5, vsi5_err, psi5, psi5_err, dsi5, dsi5_err, vsi, vsi_err, psi, psi_err, dsi, dsi_err))
f.close()

