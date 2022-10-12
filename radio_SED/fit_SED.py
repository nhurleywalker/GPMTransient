import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from matplotlib.ticker import FormatStrFormatter

import matplotlib.font_manager
from matplotlib import rc
# Nature requires sans-serif fonts
plt.rcParams.update({
    "text.usetex": True,
    "font.size": 7,
    "font.sans-serif": ["Helvetica"]})

#plt.rcParams.update({
#    "text.usetex": True,
#    "font.family": "serif",
#    "font.size": 14})

cm = 1/2.54

# TOPCSV='GLEAM-X J065228.65-255000.34.csv'
#TOPCSV = 'three_point_spectrum.csv'

markers = {
    'MWA':{'marker':'o', 'color':'black', 'label':'GLEAM-X', 'markeredgewidth':0.1, 'elinewidth':0.2, 'markersize':1},
    'MeerKAT':{'marker':'s', 'color':'green', 'label':'MeerKAT', 'markeredgewidth':0.1, 'elinewidth':0.2, 'markersize':1},
    'ASKAP':{'marker':'x', 'color':'blue', 'label':'ASKAP', 'markersize':0.5},
    'Parkes':{'marker':'+', 'color':'red', 'label':'Parkes'}
}

ref_nu = 1000

def pl(nu, norm, alpha):
    return norm * nu **alpha 

def curved_law(nu, s_nu, alpha, q):
    spec_nu = nu / ref_nu
        
    return s_nu * spec_nu ** alpha * \
            np.exp(q * np.log(spec_nu)**2)

def curved_law_integral(nu, s_nu, alpha, q, nu_min, nu_max):
    pass

def curved_law_luminosity(nu, s_nu, alpha, q, nu_min, nu_max, P, Pdot, d):
    pass

def make_ax1(ax1, nu, df):
    
    mwa_mask = df.freq < 4000

    ax1.errorbar(
        df.freq[mwa_mask],
        df.flux[mwa_mask],
        yerr=df.fluxerr[mwa_mask],
        linestyle='None',
        **markers['MWA']
    )
# MeerKAT
    corr = 0.085
    df3 = pd.read_csv("radio_SED3.csv")
    df3['flux'] /= corr
    df3['fluxerr'] /= corr
    ax1.errorbar(
        df3.freq,
        df3.flux,
        yerr=df3.fluxerr,
        linestyle='None',
        **markers['MeerKAT']
    )
# MeerKAT


# Good initials for two_comp_model
#    fit_p0 = [np.median(df.flux), -0.7, 1.5, 2]
# Initialisation for cureved_law
    fit_p0 = [np.median(df.flux), -0.7, 0.0]
    fit_res = curve_fit(
        curved_law,
        df.freq,
        df.flux,
        fit_p0,
        sigma=df.fluxerr,
        absolute_sigma=True
    )

    best_p = fit_res[0]

    no_samps = 1000
    samps = np.random.multivariate_normal(
        fit_res[0], fit_res[1], size=no_samps
    ).swapaxes(0,1)
    
    #freq = nu
    models = curved_law(
        nu[:, None],
        *samps
    )

    q16, q50, q84 = np.percentile(models, [16, 50, 84], axis=1)

    alpha = best_p[1]
    q = best_p[2]
    S1GHz = curved_law(1000, *best_p)
    
    ax1.plot(
        nu,
        q50,
        lw=0.5,
    )
    ax1.fill_between(
        nu,
        q16, q84, alpha=0.3
    )
    #print("fit res=", fit_res)

    covar = fit_res[1]
    err_p = np.sqrt(np.diag(covar))
    print(err_p)
    dof = len(df.freq) - 2
    chi2 = np.sum(
        ((df.flux - curved_law(df.freq, *best_p)) / df.fluxerr)**2
        )
    rchi2 = chi2 / dof

# Power-law (just for comparison, not used)

    model = (pl, (np.median(df.flux), -0.7), 'Power Law')

    nu = np.geomspace(70, 3000, 100)
    fit_func = model[0]
    fit_p0 = model[1]
    fit_res = curve_fit(
        fit_func,
        df.freq,
        df.flux,
        fit_p0,
        sigma=df.fluxerr,
        absolute_sigma=True
    )

    no_samps = 1000
    samps = np.random.multivariate_normal(
        fit_res[0], fit_res[1], size=no_samps
    ).swapaxes(0,1)
    
    freq = nu
    models = fit_func(
        nu[:, None],
        *samps
    )

    q16, q50, q84 = np.percentile(models, [16, 50, 84], axis=1)
    
    ax1.plot(
        nu,
            q50,
            lw=0.5,
        )
    ax1.fill_between(
        nu,
        q16, q84, alpha=0.3
        )

    ax1.set(
        xscale='log',
        yscale='log',
        xlabel='Frequency (MHz)',
        ylabel='Flux density (Jy)',
    )
    ax1.xaxis.set_major_formatter(FormatStrFormatter('%3.0f'))
    ax1.yaxis.set_major_formatter(FormatStrFormatter('%3.2f'))
# For ATCA points
#    ax1.yaxis.set_major_formatter(FormatStrFormatter('%3.5f'))
    return S1GHz, alpha, q, rchi2

def make_sed_figure(df, output='radio_SED.pdf'):

# For ATCA points, eventually?
#    example_nu_large = np.geomspace(70, 9000, 100)

    example_nu_large = np.geomspace(70,2000, 100)

    ax1_loc = (0.1, 0.1, 0.8, 0.8)

    fig = plt.figure(figsize=(7*cm, 7*cm))

    ax1 = fig.add_axes(ax1_loc)
    S1GHz, alpha, q, rchi2 = make_ax1(ax1, example_nu_large, df)

    fig.savefig(output, bbox_inches="tight")

    return S1GHz, alpha, q, rchi2

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Make some nice SED figure')

    parser.add_argument('-o','--output', type=str, default='radio_SED.pdf', help='Name of the output file')

    args = parser.parse_args()

# MWA - MWA
    df1 = pd.read_csv("radio_SED1.csv")
# MWA - Parkes
    df2 = pd.read_csv("radio_SED2.csv")
# Scale df2 to match df1
    mg = pd.merge(df1, df2, how='inner', on=['freq'])
    corr = np.average(mg['flux_x']/mg['flux_y'])
    df1['flux'] /= corr
    df1['fluxerr'] /= corr

# ASKAP
    corr = 0.95
    df4 = pd.read_csv("radio_SED4.csv")
    df4['flux'] /= corr
    df4['fluxerr'] /= corr

    df = pd.concat([df1, df2, df4])

    #print(df)
    
    S1GHz, alpha, q, rchi2 = make_sed_figure(df,
        output=args.output
    )


    print(S1GHz, alpha, q, rchi2)
    # From Wolfram Alpha
    # integrate (Power[\(40)Divide[x,Power[10,9]]\(41),-1.08])*exp(-0.60*ln(Power[\(40)Divide[x,Power[10,9]]\(41),2])) over x = Power[10,7] to Power[10,14]

    intr = 3.06e11

    # From YMW 2017 calculator
    # http://119.78.162.254/dmodel/index.php?mode=Gal&gl=22.15264754&gb=-2.06298163&dm=275&DM_Host=&ndir=1
    d = 5.8
    kpc = 3.086e+19
    Jy2Wm = 1.e-26
    Wm2ergs = 1.e7

    # Very small opening angle
    rho = 0.2
    print("Radio luminosity {0:2.2e} erg/s for rho={1}deg".format(S1GHz * intr * Jy2Wm * Wm2ergs *(2*np.pi*((d * kpc)**2)/(0.06))*(1-np.cos(np.radians(rho))),rho))

    # More reasonable, maybe?
    rho = 2
    print("Radio luminosity {0:2.2e} erg/s for rho={1}deg".format(S1GHz * intr * Jy2Wm * Wm2ergs *(2*np.pi*((d * kpc)**2)/(0.06))*(1-np.cos(np.radians(rho))),rho))

    # Spin-down luminosity
    P = 1318.19578
    Pdot = 1.e-13

    print("Spin-down luminosity < {0:2.2e} erg/s".format((4*np.pi*np.pi*1.e45*Pdot)/(P**3)))


