import matplotlib
from colossus.cosmology import cosmology
from colossus.lss import mass_function
from scipy import integrate
from astropy import constants as const
from astropy.cosmology import Planck15
import pickle as pk

import platform; plat = platform.system()
if plat=='Darwin': matplotlib.use('MacOSX')
else: matplotlib.use('TkAgg')

matplotlib.rc('font',**{'family':'serif','sans-serif':['Times'],'size':11})
matplotlib.rc('text',usetex=True)

from matplotlib import gridspec
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.pyplot as plt
plt.rc('text.latex', preamble=r'\usepackage{newtxtext}')
plt.rc('text.latex', preamble=r'\usepackage{newtxmath}')
plt.rc('font',family='serif',size=11)
plt.rc('text',usetex=True)

import scipy.special
import scipy.stats

import numpy  as np
import pandas as pd

import dill
import os

from masszeta_tool import *

cmap = {'blue'   : ['#08519C','#D2DCEA'],
        'green'  : ['#448C8A','#DDE8E7'],
        'red'    : ['#D53E4F','#F3DBDD'],
        'yellow' : ['#9C0000','#9C0000']}

perc = {'1sigma': [0.50-0.50*scipy.special.erf(1.00/np.sqrt(2.00)),0.50+0.50*scipy.special.erf(1.00/np.sqrt(2.00))],
        '2sigma': [0.50-0.50*scipy.special.erf(2.00/np.sqrt(2.00)),0.50+0.50*scipy.special.erf(2.00/np.sqrt(2.00))],
        '3sigma': [0.50-0.50*scipy.special.erf(3.00/np.sqrt(2.00)),0.50+0.50*scipy.special.erf(3.00/np.sqrt(2.00))]}

dpi = 72.27*390.00/504.00
factorx = 1.50 # 0.75
factory = 0.75 # 0.50
figsize = (factorx*504.00/dpi,factory*504.00/dpi)

# Functions for performing the theory integrals
######################################################################################################################
def log_int(function, x0, xf, nstep = 1000):
    integral = 0
    domain = np.logspace(x0, xf, nstep)
    for i in range(len(domain)-1):
        mtemp = (domain[i+1]-domain[i])/2
        integral += function(domain[i])*(np.log(domain[i+1])-np.log(domain[i]))
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
        integrand = lambda m: mass_function.massFunction(m, ztemp, mdef = '500m', model = 'tinker08', q_out = 'dndlnM')
        integrated_density = log_int(integrand, mass, 20.5)
        #Comoving volume computed over full sky i.e. 4pi steradians. Scale accordingly
        comov_vol = Planck15.comoving_volume(ztemp).value
        integral += integrated_density*comov_vol*(z[i+1]-z[i])
    return integral

# Set cosmology here
cosmology.setCosmology('planck15')


# Open spreadsheets with cluster properties
#######################################################################################################################

actlis = pd.read_excel (r'massrich_act.xlsx')
moolis = pd.read_excel (r'massrich_dat.xlsx')

# Select the mass-richness relations to use
#######################################################################################################################
# 'do'     : plot the data points if True
# 'name'   : reference name of the data used to derive the mass-richness relation
# 'merger' : True/False for using the mass-richness scaling copmuted including/excluding marger clusters
# 'nosign' : same as 'merger', but for the non detection
# 'color'  : color of the data points in the final plot

fitlis = [{'do':  True, 'name': 'carma', 'merger': False, 'nosign': False, 'color':   cmap['blue']},
          {'do': False, 'name': 'must2', 'merger': False, 'nosign':  True, 'color':    cmap['red']},
          {'do':  True, 'name': 'total', 'merger': False, 'nosign':  True, 'color': cmap['yellow']}]

# Generate mass-redshift distributions
#######################################################################################################################

fig = plt.figure(figsize=figsize)
gs = fig.add_gridspec(1,1,hspace=0,height_ratios=[1])
ax = fig.add_subplot(gs[0])

totrich = np.array([totrich1[0][totrich1[0]!=0.00],totrich1[1][totrich1[0]!=0.00],totrich1[2][totrich1[0]!=0.00]])
totzeta = np.array([totphot1[0][totrich1[0]!=0.00],totphot1[1][totrich1[0]!=0.00],totphot1[2][totrich1[0]!=0.00]])

mad_dict = {}

for f, fit in enumerate(fitlis):
  totm500 = np.zeros((3,totrich.shape[1]))

  if fit['do']:
    strm = 'wt-merger' if fit['merger'] else 'no-merger'
    strs = 'wt-nosign' if fit['nosign'] else 'no-nosign'

    if fit['name']=='advact':
      chain = dill.load(open('support/lm_chain_1_total_no-merger_no-nosign_wt-advact.pickle','rb'))
    else:
      chain = dill.load(open('support/lm.chain.1.{0}.{1}.{2}.pickle'.format(fit['name'],strm,strs),'rb'))

    alpha = np.copy(chain['alpha'][1000:-1:10])
    beta  = np.copy(chain['beta'][1000:-1:10])

    alpha = np.broadcast_to(alpha[...,np.newaxis],(alpha.shape[0],totrich.shape[-1]))
    beta  = np.broadcast_to(beta[...,np.newaxis],alpha.shape)

    linrich = scipy.stats.truncnorm.rvs(loc=totrich[0],scale=totrich[1],a=-totrich[0]/totrich[1],b=np.inf,size=alpha.shape)

    linm500 = alpha+np.log10(linrich)*beta
    linm500 = np.percentile(linm500,[100*perc['1sigma'][0],50,100*perc['1sigma'][1]],axis=0)

    totm500[0] = 10**(linm500[1])
    totm500[1] = 10**(linm500[1])*np.log(10.00)*np.abs(linm500[0]-linm500[1])
    totm500[2] = 10**(linm500[1])*np.log(10.00)*np.abs(linm500[2]-linm500[1])

    ax.errorbar(totzeta[0],totm500[0],yerr=totm500[1:],fmt='.',color=fit['color'][0])
    
    mad_dict[fit['name']] = [totzeta[0],totm500[0]]

    print(fit['name'])
    print(totm500[0])
# Plot other data points
#######################################################################################################################

pk.dump(mad_dict, open("mad_dict.p", "wb"))

ax.plot(plczeta[0],plcm500[0],'x',color='#63A6A0',markersize=2,label=r'\textit{Planck}')
ax.plot(sptzeta[0],sptm500[0],'D',color='#B4D9CC',markersize=2,label=r'SPT')
ax.plot(actzeta[0],actm500[0],'v',color='#89C0B6',markersize=2,label=r'ACT')

ax.set_xlim(0.00E+00,2.00E+00)
ax.set_ylim(0.00E+00,1.50E+01)

ax.set_xlabel(r'$z$')
ax.set_ylabel(r'$M_{500}~\mathrm{[10^{14}M_{\odot}]}$')
ax.tick_params(axis='both',which='both',direction='in',top=True,left=True,right=True)
ax.legend()

plt.savefig('halo_mass.png')
plt.show()
