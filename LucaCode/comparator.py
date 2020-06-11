import numpy  as np
import matplotlib.pyplot as plt 
from colossus.cosmology import cosmology 
from colossus.lss import mass_function 
from scipy import integrate 
from astropy import constants as const 
from astropy.cosmology import Planck15 
import pickle as pk 

# Functions for performing the theory integrals
######################################################################################################################
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
        mass = np.log10(mass)
        #print(mass)
    for i in range(len(z)-1):
        ztemp = (z[i+1]+z[i])/2
        #print(i)
        integrand = lambda m: mass_function.massFunction(m, ztemp, mdef = '500m', model = 'tinker08', q_out = 'dndlnM')
        integrated_density = log_int(integrand, mass, 20.5)
        #Comoving volume computed over full sky i.e. 4pi steradians. Scale accordingly
        comov_vol = Planck15.comoving_volume(ztemp).value
        integral += integrated_density*comov_vol*(z[i+1]-z[i])
    return integral

def arr_count(array, low, up):
    n = 0
    for z in array:
        if low<= z < up: n += 1
    return n

# Set cosmology here
cosmology.setCosmology('planck15')



mad_dict = pk.load(open("mad_dict.p", 'rb'))

z_step = 0.25

mass_cut = 4

carma_z = mad_dict['carma'][0][np.where(mad_dict['carma'][1]>mass_cut)]
print(mad_dict['carma'][1])
print(carma_z)

datum = []
theorys = []

for z in np.arange(0.75, 1.75, z_step):
    data = arr_count(carma_z, z, z + z_step)
    datum.append(data)
    print("Data: {}".format(data))
    theory = 0.6*zint(mass_cut*10**14, z, z + z_step)
    theorys.append(theory)
    print("Theory: {}".format(theory))

zs = [0.875, 1.125, 1.375, 1.625]

plt.scatter(zs, datum, label = 'Data')
plt.scatter(zs, theorys, label = 'Theory')
plt.title('Carma scaling relation for M_500 > 4*10^14 vs theory assuming 60% sky coverage')
plt.legend()
plt.savefig('madcows_v_cosmo.png')
plt.show()
