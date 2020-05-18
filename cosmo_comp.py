import matplotlib.pyplot as plt
import numpy as np
from colossus.cosmology import cosmology
from colossus.lss import mass_function
from scipy import integrate
from astropy import constants as const

cosmology.setCosmology('planck18')

nstep = 5

z = np.linspace(.1,.2, nstep)

for i in range(len(z)):

    integrand = lambda m: mass_function.massFunction(m, z[i], mdef = '200m', model = 'tinker08', q_out = 'dndlnM')
    integrated_density = integrate.quad(integrand, 1e15, np.inf)[0]
    comov_vol = 4*np.pi*const.c.value
    print(comov_vol)
    print("Integrated density at z = {}: {}".format(z[i], integrated_density))

