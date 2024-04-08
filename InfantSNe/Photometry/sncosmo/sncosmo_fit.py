import sncosmo
import pdb
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import ascii

file='2020jgl_standard_sdss.dat'
data=ascii.read(file)

redshift=[0.006758,0.01543852,0.083562,0.0319161,0.039761]

if file[4:7]=='jgl': k=0
elif file[4:7]=='jhf': k=1
elif file[4:7]=='kku': k=2
elif file[4:7]=='kyx': k=3
else: k=4

# create a model
model = sncosmo.Model(source='salt2')

#mask=(data['time']<=59060)

model.set(z=redshift[k])

# run the fit
result, fitted_model = sncosmo.fit_lc(
    data, model,
    ['t0', 'x0', 'x1', 'c'])#,  # parameters of model to vary
    #bounds={'z':(0, 0.05)})  # bounds on parameters (if any)

fig=sncosmo.plot_lc(data, model=fitted_model, errors=result.errors)
plt.show()

#pdb.set_trace()
