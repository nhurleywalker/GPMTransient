import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.special import erfi, erf

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

def curved_law_integral(nu_min, nu_max, s_nu, alpha, q, print_terms=False):
    '''
    Definite integral for the function defined in curved_law()
    '''
    Q = 2*np.sqrt(np.abs(q))
    R = (alpha + 1)/Q
    if q > 0:
        exp_term = np.exp(-R**2)
        erf_min_term = erfi(R + 2*q*np.log(nu_min/ref_nu)/Q)
        erf_max_term = erfi(R + 2*q*np.log(nu_max/ref_nu)/Q)
    else:
        exp_term = np.exp(R**2)
        erf_min_term = erf(-(R + 2*q*np.log(nu_min/ref_nu)/Q))
        erf_max_term = erf(-(R + 2*q*np.log(nu_max/ref_nu)/Q))

    frontmatter = s_nu * (ref_nu*1e6) * np.sqrt(np.pi) / Q

    if print_terms:
        print(f"s_nu ν0 (1/2) √(π/q)    = {frontmatter}")
        print(f"exp(|R|^2)              = {exp_term}")
        print(f"-R                      = {-R}")
        print(f"-2q/Q                   = {-2*q/Q}")
    return frontmatter * exp_term * (erf_max_term - erf_min_term)

def curved_law_numeric_integral(nu_min, nu_max, s_nu, alpha, q, nsteps):
    '''
    Numerical integration of curved law, used to check the analytic expression for errors
    '''

    # Area divided up into NSTEPS slivers in log-space, with the top of each sliver
    # being approximated with a power law: S = Kν^α
    nus = np.logspace(np.log10(nu_min), np.log10(nu_max), nsteps+1)
    Ss  = curved_law(nus, s_nu, alpha, q) # curved_law does the division by reference frequency, so don't do it here.

    nus_lo = nus[:-1]
    nus_hi = nus[1:]
    Ss_lo  = Ss[:-1]
    Ss_hi  = Ss[1:]

    alphas = np.log(Ss_hi/Ss_lo) / np.log(nus_hi/nus_lo)
    Ks = Ss_hi / (nus_hi/ref_nu)**alphas

    # The area of each sliver is the integral of a power law:
    # ∫ Kν^α dν = K ν^(α+1) / (α+1)
    areas = (Ks / (alphas + 1)) * ((nus_hi/ref_nu)**(alphas + 1) - (nus_lo/ref_nu)**(alphas + 1))
    # but dν needs to be converted from GHz to Hz as well, so the area is stretched by
    areas *= 1e9

    # Return the total area
    return np.sum(areas)

def f(P, Pdot):
    return 4.68e-3 * (Pdot/1e-15)**0.07 * P**(-0.7)

def curved_law_luminosity_Speak(nu_min, nu_max, s_nu, alpha, q, P, Pdot, d, print_terms=False):
    '''
    See Derivation #3 in README.md
    '''

    beta = -0.26
    integral_value = curved_law_integral(nu_min, nu_max, s_nu, alpha + beta, q, print_terms=print_terms)

    if print_terms:
        print("integral = ", integral_value)

    nsteps = 10000
    numeric_integral_value = curved_law_numeric_integral(nu_min, nu_max, s_nu, alpha + beta, q, nsteps)

    if print_terms:
        print("numeric integral for ", nsteps, " steps = ", numeric_integral_value)

    if print_terms:
        print(f"4πd²        = {4*np.pi*d**2}")
        print(f"f(P, Pdot)  = {f(P, Pdot)}")

    return 4 * np.pi * d**2 * f(P, Pdot) * integral_value

def curved_law_luminosity_Smean(nu_min, nu_max, s_nu, alpha, q, P, Pdot, d, sin_mag_incl):
    '''
    See Derivation #4 in README.md
    '''

    beta = -0.26
    return 2 * np.pi**2 * d**2 * sin_mag_incl * np.sqrt(f(P, Pdot)) * curved_law_integral(nu_min, nu_max, s_nu, alpha + beta/2, q)

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
    corr = 0.072
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


# Broken power-law -- two different power laws

#    model = (pl, (np.median(df.flux[df.freq < 500]), -0.7), 'Power Law')
#
#    nu = np.geomspace(70, 500, 100)
#    fit_func = model[0]
#    fit_p0 = model[1]
#    fit_res = curve_fit(
#        fit_func,
#        df.freq[df.freq < 500],
#        df.flux[df.freq < 500],
#        fit_p0,
#        sigma=df.fluxerr[df.freq < 500],
#        absolute_sigma=True
#    )
#
#    no_samps = 1000
#    samps = np.random.multivariate_normal(
#        fit_res[0], fit_res[1], size=no_samps
#    ).swapaxes(0,1)
#    
#    freq = nu
#    models = fit_func(
#        nu[:, None],
#        *samps
#    )
#
#    q16, q50, q84 = np.percentile(models, [16, 50, 84], axis=1)
#    
#    ax1.plot(
#        nu,
#            q50,
#            lw=0.5,
#            color='green'
#        )
#    ax1.fill_between(
#        nu,
#        q16, q84, alpha=0.3,
#        color='green'
#        )
#
#
#    model = (pl, (np.median(df.flux[df.freq > 500]), -3), 'Power Law')
#
#    nu = np.geomspace(500, 2000, 100)
#    fit_func = model[0]
#    fit_p0 = model[1]
#    fit_res = curve_fit(
#        fit_func,
#        df.freq[df.freq > 500],
#        df.flux[df.freq > 500],
#        fit_p0,
#        sigma=df.fluxerr[df.freq > 500],
#        absolute_sigma=True
#    )
#
#    no_samps = 1000
#    samps = np.random.multivariate_normal(
#        fit_res[0], fit_res[1], size=no_samps
#    ).swapaxes(0,1)
#    
#    freq = nu
#    models = fit_func(
#        nu[:, None],
#        *samps
#    )
#
#    q16, q50, q84 = np.percentile(models, [16, 50, 84], axis=1)
#    
#    ax1.plot(
#        nu,
#            q50,
#            lw=0.5,
#            color='green'
#        )
#    ax1.fill_between(
#        nu,
#        q16, q84, alpha=0.3,
#        color='green'
#        )

# Power-law (just for comparison)

    model = (pl, (np.median(df.flux), -0.7), 'Power Law')

    nu = np.geomspace(70, 2000, 100)
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
    print(fit_res)

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

    best_p = fit_res[0]
    pla = best_p[1]
    plS = pl(1000, *best_p)
    return S1GHz, alpha, q, rchi2, plS, pla

def make_sed_figure(df, output='radio_SED.pdf'):

# For ATCA points, eventually?
#    example_nu_large = np.geomspace(70, 9000, 100)

    example_nu_large = np.geomspace(70,2000, 100)

    ax1_loc = (0.1, 0.1, 0.8, 0.8)

    fig = plt.figure(figsize=(8.9*cm, 8*cm))

    ax1 = fig.add_axes(ax1_loc)
    S1GHz, alpha, q, rchi2, plS, pla = make_ax1(ax1, example_nu_large, df)

    fig.savefig(output, bbox_inches="tight", dpi=300)

    return S1GHz, alpha, q, rchi2, plS, pla

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
    corr = 0.7
    df4 = pd.read_csv("radio_SED4.csv")
    df4['flux'] /= corr
    df4['fluxerr'] /= corr

    df = pd.concat([df1, df2, df4])

    #print(df)
    
    S1GHz, alpha, q, rchi2, plS, pla = make_sed_figure(df,
        output=args.output
    )


    print("Curved spectrum fit parameters: S at 1 GHz = {0:3.2f}mJy, alpha = {1:3.2f}, q = {2:3.2f}, reduced chi^2 = {3:3.2f}".format(1.e3*S1GHz, alpha, q, rchi2))
    print("Power-law fit parameters: S at 1 GHz = {0:3.2f}mJy, alpha = {1:3.2f}".format(1.e3*plS, pla))


# From Wolfram Alpha
    # integrate (Power[\(40)Divide[x,Power[10,9]]\(41),-1.08])*exp(-0.60*ln(Power[\(40)Divide[x,Power[10,9]]\(41),2])) over x = Power[10,7] to Power[10,14]

    intr = 3.06e11

    # From YMW 2017 calculator
    # http://119.78.162.254/dmodel/index.php?mode=Gal&gl=22.15264754&gb=-2.06298163&dm=275&DM_Host=&ndir=1
    d = 5.8
    kpc = 3.086e+19
    Jy2Wm = 1.e-26
    Wm2ergs = 1.e7
    P = 1318.19578
    Pdot = 1.e-13

    # Formula for prefactor of rho (in degrees)
    rhop = 1.24 * ( 40 * (Pdot/1.e-15)**0.07 * P**0.3 )**(0.5) * (P**-0.5)

    print(rhop * 0.1**(-0.26/2))
    
    # From Wolfram alpha, integral of rho^2 S as a function of frequency:
    # integrate ((Power[v,b])Power[\(40)v\(41),a])*exp(q*ln(Power[\(40)v\(41),2]))
    def integral(numin, numax, q, alpha, beta=-0.26):
        return (numax ** (2*q + alpha + beta + 1)) / (2*q + alpha + beta + 1) - \
               (numin ** (2*q + alpha + beta + 1)) / (2*q + alpha + beta + 1)

    # Doing the integral properly
    print("Radio luminosity {0:2.2e} erg/s for frequency-dependent rho".format(Jy2Wm * Wm2ergs * curved_law_luminosity_Speak(1.e7/1.e6, 1.e15/1.e6, S1GHz, alpha, q, P, Pdot, d * kpc, print_terms=True)))
    print(f"Unit factors: Jy2Wm * Wm2ergs    = {Jy2Wm * Wm2ergs}")

    # Doing the integral properly
    print("Radio luminosity {0:2.2e} erg/s, using Smean, for frequency-dependent rho and duty-cycle".format(Jy2Wm * Wm2ergs * curved_law_luminosity_Smean(1.e7/1.e6, 1.e15/1.e6, S1GHz, alpha, q, P, Pdot, d * kpc, 1)))

    # Treating it as a power law
    #print("Radio luminosity {0:2.2e} erg/s for frequency-dependent rho".format(S1GHz * Jy2Wm * Wm2ergs *(np.pi*((d * kpc)**2)) * rhop**2 * integral(1.e7/1.e9, 1.e15/1.e9, q, alpha)))

    rho = 0.2 # degrees
    # Vs. Lorimer & Kramer 2012 single power-law
    print("Radio luminosity {0:2.2e} erg/s for simple power-law fit".format(Jy2Wm * Wm2ergs * 2 * np.pi * (d * kpc)**2 * (1- np.cos(np.radians(rho))) * plS * ( (1.e9)** (-pla) / (pla + 1) ) * (2.e9 ** (pla + 1) - 72.e6 ** (pla + 1) ) ) )
    


    # Spin-down luminosity

    print("Spin-down luminosity < {0:2.2e} erg/s".format((4*np.pi*np.pi*1.e45*Pdot)/(P**3)))


