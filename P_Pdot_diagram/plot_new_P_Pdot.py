#!/usr/bin/env python

import matplotlib 
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import rc, rcParams
#import adjustText
#from adjustText import adjust_text
import os
import math
import numpy as np
from scipy import constants
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap
from matplotlib import colors
import scipy
import scipy.constants as sc
import astropy
import astropy.constants as astro_cons

cm = 1/2.54  # centimeters in inches

# Nature requires sans-serif fonts
plt.rcParams.update({
    "font.size": 7,
    "font.sans-serif": ["Helvetica"],
    "figure.figsize" : (17.8*cm,11.9*cm)})

inside_text_size = 5
my_cmap = "viridis"
# Move text slightly to the right of scatter points
offpad = 1.3

electron_mass_cgs = astro_cons.m_e.cgs.value
electron_charge_cgs = astro_cons.e.gauss.value
solar_mass_cgs = astro_cons.M_sun.cgs.value
planck_bar_cgs = astro_cons.hbar.cgs.value
light_c_cgs = astro_cons.c.cgs.value

df_sources_data = pd.read_csv('PPdotB_LumkTage_allclasses_04042023.csv')
#atnf_data = pd.read_csv('atnf-catalog-10-02-2023.csv', delimiter=';')
#atnf_data_clean = atnf_data[atnf_data['BSURF']!='*']
#df_models_data = pd.read_csv('cooling_curves_1.4Msun_ppdot_all_1e6yr.csv')
df_atnf_data = pd.read_csv('ATNF_catalog_Pgr0.01_25062021.csv', delimiter=';')

# Colorp MAP source data with P Pdot and B
texts = []
my_norm = colors.LogNorm(vmin=0.0001, vmax=50)

# Dictionaries to reuse
magnetargs = {'norm' : my_norm, 'marker' : '*', 'alpha' : 0.9, 'edgecolors' : 'black', 'linewidths' : 0.5, 's' : 60}
ulpmargs = { 'norm' : my_norm, 'alpha' : 1, 'marker' : 'o', 'edgecolors' : 'black', 's' : 40, 'linewidths' : 0.5 , 'color' : 'grey'}
rcwargs = { 'norm' : my_norm, 'alpha' : 0.9, 'marker' : '*', 'edgecolors' : 'black', 's' : 60, 'linewidths' : 0.5 , 'color' : 'grey'}
xdinargs = { 'norm' : my_norm, 'alpha' : 0.8, 'marker' : 's', 'edgecolors' : 'black', 's' : 30, 'linewidths' : 0.5 , 'zorder' : 3}
ccoargs = { 'norm' : my_norm, 'alpha' : 0.9, 'marker' : '^', 'edgecolors' : 'black', 's' : 40, 'linewidths' : 0.5 }
pulsargs = { 'norm' : my_norm, 'alpha' : 0.9, 'marker' : '.', 'edgecolors' : 'black', 's' : 30, 'linewidth' : 0, 'zorder' : -3}
radiomagnetargs = {'facecolors' : 'none', 'alpha' : 0.7, 'edgecolors' : 'black', 's' : 60, 'linewidths' : 0.5 }
arrowargs = { 'zorder' : -1, 'alpha' : 0.7, 'arrowprops' : dict(arrowstyle="->", lw=0.5, color='grey') }

# Magnetars 
selection = (df_sources_data['class']=='Magnetar') 
plt.scatter(x=df_sources_data[selection]['period(s)'], 
            y=df_sources_data[selection]['pdot(1e-11s/s)'], 
            c=df_sources_data[selection]['Bp(e14G)'], cmap=my_cmap,
            label='Magnetar-like emission',
            **magnetargs)

for index, row in df_sources_data[selection].iterrows():
  # Plot radio magnetars
    if row['radio']==True:
        plt.scatter(row['period(s)'], row['pdot(1e-11s/s)'], **radiomagnetargs)
              
    # Plot Pdot upperlimits
    if row['pdot_ul']==True:
        plt.annotate("", xy=(row['period(s)'], row['pdot(1e-11s/s)']/10), 
                    xytext=(row['period(s)'], row['pdot(1e-11s/s)']), **arrowargs)
    

sm = plt.cm.ScalarMappable(cmap=my_cmap, norm=my_norm)
cbar = plt.colorbar(sm, pad=0.022, aspect=30)
cbar.ax.set_ylabel(r'Surface dipolar B-field at pole ($10^{14}$ Gauss)')

#XDINS 
selection = (df_sources_data['class']=='XDINS')
plt.scatter(x=df_sources_data[selection]['period(s)'], 
            y=df_sources_data[selection]['pdot(1e-11s/s)'], 
            c=df_sources_data[selection]['Bp(e14G)'], cmap=my_cmap, 
             label='Thermally-emitting (XDINSs)',
             **xdinargs)

#CCOs
selection = ((df_sources_data['class']=='CCO') & (~np.isnan(df_sources_data['Edot(e33erg/s)'])))
plt.scatter(x=df_sources_data[selection]['period(s)'], 
            y=df_sources_data[selection]['pdot(1e-11s/s)'], 
            c=df_sources_data[selection]['Bp(e14G)'], cmap=my_cmap, 
            label='Central Compact Objects (CCOs)',
            **ccoargs)

#RCW103
selection = ((df_sources_data['source']=='1E161348-5055'))
plt.scatter(x=df_sources_data[selection]['period(s)'], 
            y=df_sources_data[selection]['pdot(1e-11s/s)'], 
            **rcwargs)

for index, row in df_sources_data[selection].iterrows():
    texts.append(plt.text(offpad*row['period(s)'], row['pdot(1e-11s/s)'], 
                         row['source'], ha='left', va='center', color='black', size=inside_text_size))
    # Plot Pdot upperlimits
    if row['pdot_ul']==True:
        plt.annotate("", xy=(row['period(s)'], row['pdot(1e-11s/s)']/100), 
                    xytext=(row['period(s)'], row['pdot(1e-11s/s)']), **arrowargs)
                    
        
#Long-period pulsars
selection = ((df_sources_data['source']=='PSRJ0250+5854') | (df_sources_data['source']=='PSRJ0901-4046'))
plt.scatter(x=df_sources_data[selection]['period(s)'], 
            y=df_sources_data[selection]['pdot(1e-11s/s)'], 
            c=df_sources_data[selection]['Bp(e14G)'], cmap=my_cmap, 
            **pulsargs)

#GLEAM-XJ162759-523504
selection = ((df_sources_data['source']=='GLEAM-XJ162759-523504'))
plt.scatter(x=df_sources_data[selection]['period(s)'], 
            y=df_sources_data[selection]['pdot(1e-11s/s)'], 
            #c=df_sources_data[selection]['Bp(e14G)'], cmap=my_cmap, 
            **ulpmargs)
            
for index, row in df_sources_data[selection].iterrows():
    texts.append(plt.text(offpad*row['period(s)'], row['pdot(1e-11s/s)'], 
                         row['source'].replace("J", " J"), ha='left', va='center', color='black', size=inside_text_size))
 
    # Plot Pdot upperlimits
    if row['pdot_ul']==True:
        plt.annotate("", xy=(row['period(s)'], row['pdot(1e-11s/s)']/100), 
                    xytext=(row['period(s)'], row['pdot(1e-11s/s)']), **arrowargs)
        
#GPM J1839-10
selection = ((df_sources_data['source']=='GPMJ1839-10'))
plt.scatter(x=df_sources_data[selection]['period(s)'], 
            y=df_sources_data[selection]['pdot(1e-11s/s)'], 
            **ulpmargs)

for index, row in df_sources_data[selection].iterrows():
    texts.append(plt.text(offpad*row['period(s)'], row['pdot(1e-11s/s)'], 
                         row['source'].replace("J", " J"), ha='left', va='center', color='black', size=inside_text_size))

    if row['pdot_ul']==True:
        plt.annotate("", xy=(row['period(s)'], row['pdot(1e-11s/s)']/100), 
                    xytext=(row['period(s)'], row['pdot(1e-11s/s)']), **arrowargs)


# ATNF pulsars
plt.scatter(x=df_atnf_data['P0'], 
            y=df_atnf_data['P1']/1.e-11, 
            c=df_atnf_data['BSURF']/0.5e14, cmap=my_cmap, 
             label='Radio Pulsars',
             **pulsargs)

plt.xlabel(r'Spin Period (s)')#, size=26)
plt.ylabel(r'Period derivative ($10^{-11}$ s s$^{-1}$)')#, size=26)

plt.yscale('log')
plt.xscale('log')

plt.legend(loc='lower right')

# Magnetic field lines
x = np.linspace(1e-4,1e6,3000)
def yplace(Gpow, x):
    return ((10**Gpow)**2)/(x * 1.e-11*(6.4*10**19)**2) 

for G in range (11, 17):
    plt.plot(x, yplace(G, x), 'grey', alpha=0.3, lw=1)

# Label magnetic field
xl=20
xm = 2.e3
xr=1.e5

Gtargs = {'color' : 'grey', 'alpha' : 0.9, 'size' : inside_text_size}

# Plot these in the whitespace around the legend
for G in range(11, 13):
    plt.text(xl, yplace(G, xl), r'$10^{{{G}}}$G'.format(G=G), **Gtargs)
G = 13
plt.text(xm, yplace(G, xm), r'$10^{{{G}}}$G'.format(G=G), **Gtargs)
for G in range(14, 17):
    plt.text(xr, yplace(G, xr), r'$10^{{{G}}}$G'.format(G=G), **Gtargs)

# DEATH LINES

R_NS6 = 1.2
R_NS = 1.2*1.e6
Msolar = 2*1.e33
M_NS = 1.4*Msolar
K = 4 * (3 * light_c_cgs**3) * 2 / (5 * 8 * math.pi**2)

x = np.linspace(1e-3,1e7,6000)

# Chen 1 Pure Dipole
chen1_ppdot_pd = 1/K * ((2.2*1.e12 * (R_NS6)**(-19/8) * x**(15/8))**2) * (R_NS)**4 / M_NS  / x
plt.plot(x,chen1_ppdot_pd/1.e-11, linestyle='dashed', color='black', markersize=12, alpha=1, lw=0.5)
#x11=2.2e3
#y11 = 1/K * ((2.2*1.e12 * (R_NS6)**(-19/8) * x11**(15/8))**2) * (R_NS)**4 / M_NS  / x11 / 1.e-11
#plt.text(x11, y11*1.5, 'Pure dipole', color='black', alpha=1, rotation=39, size=34)


# Chen 2 Twisted dipole
chen2_ppdot_twd = 1/K * ((2.7*1.e11 * (R_NS6)**(-17/8) * x**(13/8))**2) * (R_NS)**4 / M_NS  / x
plt.plot(x,chen2_ppdot_twd/1.e-11, linestyle='dotted', color='black', markersize=12, alpha=1, lw=0.5)


# Chen 4 Twisted Multipolar Spot beta=10
beta = 10
chen4_ppdot_tw_multi = 1/K * ((9.2*1.e10 * beta**(-1/4) * (R_NS6)**(-2) * x**(3/2))**2) * (R_NS)**4 / M_NS  / x
plt.plot(x,chen4_ppdot_tw_multi/1.e-11, linestyle='solid', color='black', markersize=12, alpha=1, lw=0.5)



# Zhang III 
zhangIII_ppdot = 1/K * ((9.2*1.e25 * (R_NS)**(-9/4) * x**(7/4))**2) * (R_NS)**4 / M_NS  / x
plt.plot(x,zhangIII_ppdot/1.e-11, linestyle='dashed', color='black', markersize=12, alpha=1, lw=0.5)


# Zhang III prime
zhangIIIprime_ppdot = 1/K * ((3.5*1.e23 * (R_NS)**(-2) * x**(3/2))**2) * (R_NS)**4 / M_NS  / x
plt.plot(x,zhangIIIprime_ppdot/1.e-11, linestyle='dotted', color='black', markersize=12, alpha=1, lw=0.5)
plt.fill_between(x, chen1_ppdot_pd/1.e-11, chen4_ppdot_tw_multi/1.e-11, color='grey', alpha=0.2)

#
plt.xlim(8e-3, 1.e6)
plt.ylim(1.001e-9, 2000)
plt.clim(0.0001, 50)
plt.savefig('P_Pdot_deathlines.pdf', bbox_inches="tight")
