import matplotlib.pyplot as plt
import numpy as np
from colossus.cosmology import cosmology
from colossus.lss import mass_function
from scipy import integrate
from astropy import constants as const
from astropy.cosmology import Planck15

def log_int(function, x0, xf, nstep = 1000):
    integral = 0
    domain = np.logspace(x0, xf, nstep)
    for i in range(len(domain)-1): 
        integral += function(domain[i])*(np.log(domain[i+1]-domain[i]))
    return integral

def zint(mass, zmin, zmax, nstep = 1000):
    z = np.linspace(zmin, zmax, nstep)
    
    integral = 0
    
    if mass > 100:
        #Shitty check if mass was given as just the exponant
        mass = np.log(mass)

    for i in range(len(z)-1):
        integrand = lambda m: mass_function.massFunction(m, z[i], mdef = '200m', model = 'tinker08', q_out = 'dndlnM')
        integrated_density = log_int(integrand, mass, 20.5)
        #Comoving volume computed over full sky i.e. 4pi steradians. Scale accordingly
        comov_vol = Planck15.comoving_volume(z[i]).value
        integral += integrated_density*comov_vol*(z[i+1]-z[i])
    return integral

cosmology.setCosmology('planck15')

print(zint(15, 0.1, 0.2, nstep = 10000))



