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
    domain = np.log10(domain)
    for i in range(len(domain)-1): 
        #print(domain[i])
        #mtemp = (domain[i+1]-domain[i])/2
        #print(function(10**domain[i]))
        #print(mass_function.massFunction(10**domain[i], 0.15, mdef = '200m', model = 'tinker08', q_out = 'dndlnM'))
        integral += function(10**domain[i])*(domain[i+1] - domain[i])
    return integral

def zint(mass, zmin, zmax, nstep = 1000):
    z = np.linspace(zmin, zmax, nstep)
    
    integral = 0
    
    if mass > 100:
        #Shitty check if mass was given as just the exponant
        mass = np.log(mass)

    for i in range(len(z)-1):
        ztemp = (z[i+1]+z[i])/2
        print(i)
        integrand = lambda m: mass_function.massFunction(m, ztemp, mdef = '200m', model = 'tinker08', q_out = 'dndlnM')
        integrated_density = log_int(integrand, mass, 20.5)
        #Comoving volume computed over full sky i.e. 4pi steradians. Scale accordingly
        comov_vol = Planck15.comoving_volume(ztemp).value
        integral += integrated_density*comov_vol*(z[i+1]-z[i])
    return integral

cosmology.setCosmology('planck15')

print(zint(14, 0.1, 0.2, nstep = 1000))



