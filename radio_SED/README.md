# Radio SEDs and Luminosity

## Luminosity calculations along the lines of [Erkut (2022)](https://ui.adsabs.harvard.edu/abs/2022MNRAS.514L..41E%2F/abstract)

In their Equation (15), the general expression for the luminosity over the frequency range $\nu_\text{min} \le \nu \le \nu_\text{max}$ for a beamed pulsar-like object is

$$ L = 4\pi d^2 \int_{\nu_\text{min}}^{\nu_\text{max}} \sin^2 \left[ \frac{\rho(\nu)}{2} \right] S_p(\nu) \text{d}\nu,$$

where $d$ is the distance, $\rho$ is the opening angle of the beam, and $S_p(\nu)$ is the peak flux density.
This assumes that the beam is uniformly illuminated with flux density $S_p(\nu)$.

In the case of GLEAM-X J162759.5–523504.3, using the peak flux density as a proxy for the emission beam is problematic, owing to the paucity of pulses and the large uncertainty of the average pulse shape.
Hence, it is arguably more reliable (and practical, for continuum imaging!) to use the mean flux density, $S_m(\nu)$.
In this case, it is assumed that the emission beam flux density is just the mean flux density averaged over the whole rotation divided by the duty cycle, $S_p(\nu) = S_m(\nu)/\delta$.
The duty cycle is dependent on the opening angle and the viewing geometry (see section below), so the above integral in general becomes

$$ L = 4\pi d^2 \int_{\nu_\text{min}}^{\nu_\text{max}} \sin^2 \left[ \frac{\rho(\nu)}{2} \right] \frac{S_m(\nu)}{\delta(\nu)} \text{d}\nu,$$

## Relation between duty cycle and opening angle

The general relationship between the opening angle and the duty cycle, from spherical geometry, is

$$ \sin^2\left(\frac{\pi}{2}\delta\right) = \frac{\sin^2(\rho/2) - \sin^2(\beta/2)}{\sin\alpha \sin(\alpha + \beta)},$$

where $\alpha$ is the star's magnetic inclination angle (relative to the spin axis), and $\beta$ is the impact angle (the closest approach of the magnetic axis to the line of sight).
Since Erkut's argument is that the luminosity has been overestimated, we can assume a "worst case scenario" that the impact angle is very small, since as the impact angle increases, so too must the opening angle in order to match the observed duty cycle, which would give rise to even larger total luminosities.
Thus, setting $\beta = 0$,

$$ \sin\left(\frac{\pi}{2}\delta\right) = \frac{\sin(\rho/2)}{\sin\alpha},$$

which implies a beaming fraction

$$ \sin^2\left(\frac{\rho}{2}\right) = \sin^2\left(\frac{\pi}{2}\delta\right) \sin^2\alpha. $$

### The consequences of $\rho$'s dependence on $\nu$

The most straightforward justification for the assumption of small $\rho$ is that the extremely large period implies a large light cylinder radius, which in turn implies a very small polar cap region.
In this limit, the beam opening angle is

$$ \sin^2\left(\frac{\rho}{2}\right) \approx \frac{\rho^2}{4}, $$

giving the luminosity

$$ L = \pi d^2 \int_{\nu_\text{min}}^{\nu_\text{max}} \rho(\nu)^2 \frac{S_m(\nu)}{\delta(\nu)} \text{d}\nu. $$

Moreover, the duty cycle can be estimated as

$$ \delta \approx \frac{\rho}{\pi\sin\alpha}. $$

with the key relation being $\delta \propto \rho$.

However, it appears that Erkut may have taken $\delta$ to be constant (which would be ironic, because the whole point of the paper is about including the beam opening angle in the luminosity estimation), so I will derive the luminosity formula below both ways, to see what the difference is.
In both cases, the opening angle is assumed to have the form

$$ \rho = 1.^\circ24 \left( \frac{r}{10\text{km}} \right)^{1/2} \left( \frac{P}{\text{s}} \right)^{-1/2}, $$

(see Erkut's equation (13)), where

$$ \frac{r}{10\text{km}} \approx 40 \left( \frac{\nu}{\nu_s} \right)^\beta \left( \frac{\dot{P}}{10^{-15}} \right)^{0.07} \left( \frac{P}{\text{s}} \right)^{0.30}, $$

(see Erkut's equation (14)), where $\nu_s = 1$ GHz, and $\beta \approx -0.26$ is the radius-to-frequency index, not to be confused with the impact angle above.
Hence,

$$ \rho(\nu) \approx 1.^\circ24 \sqrt{40} \left( \frac{\nu}{\nu_s} \right)^{\beta/2} \left( \frac{\dot{P}}{10^{-15}} \right)^{0.035} \left( \frac{P}{\text{s}} \right)^{-0.35}, $$

$$ \rho^2(\nu) \approx 0.0187 \left( \frac{\nu}{\nu_s} \right)^\beta \left( \frac{\dot{P}}{10^{-15}} \right)^{0.07} \left( \frac{P}{\text{s}} \right)^{-0.7}, $$

or, separating out all the factors that do not depend on frequency, and keeping the factor of 4 with the $\rho^2$,

$$ \frac{\rho(\nu)^2}{4} \approx f(P,\dot{P}) \left( \frac{\nu}{\nu_s} \right)^\beta, $$

$$ f(P,\dot{P}) = 4.68\times10^{-3} \left( \frac{\dot{P}}{10^{-15}} \right)^{0.07} \left( \frac{P}{\text{s}} \right)^{-0.7}. $$

#### Derivation #1: Constant duty cycle, single power law

Unfortunately, $\alpha$ is overloaded, being both the magnetic inclination (as used above), and the usual symbol for the spectral index.
Therefore, in this derivation, I'll use $\gamma$ for the spectral index, making the observed mean flux

$$ S_m(\nu) = S_m(\nu_0) \left( \frac{\nu}{\nu_0} \right)^\gamma. $$

The luminosity is then

$$ L = 4 \pi d^2 \int_{\nu_\text{min}}^{\nu_\text{max}} \frac{\rho(\nu)^2}{4} \frac{S_m(\nu)}{\delta} \text{d}\nu $$

$$ = 4 \pi d^2 f(P,\dot{P}) \frac{S_m(\nu_0)}{\delta} \int_{\nu_\text{min}}^{\nu_\text{max}} \left( \frac{\nu^{\gamma + \beta}}{\nu_s^\beta \nu_0^\gamma} \right) \text{d}\nu $$

$$ = 4 \pi d^2 f(P,\dot{P}) \frac{\nu_\text{max}^{\gamma + \beta + 1} - \nu_\text{min}^{\gamma + \beta + 1}}{\gamma + \beta + 1} \frac{S_m(\nu_0)}{\delta \nu_s^\beta \nu_0^\gamma}, $$

which is Erkut's Equation (17).

#### Derivation #2: Assumed profile evolution

With the estimated dependence of $\delta$ on $\rho$ (see above), the luminosity integral becomes

$$ L = 2 \pi^2 d^2 \sin\alpha \int_{\nu_\text{min}}^{\nu_\text{max}} \frac{\rho(\nu)}{2} S_m(\nu) \text{d}\nu $$

$$ = 2 \pi^2 d^2 \sin\alpha \sqrt{f(P,\dot{P})} \int_{\nu_\text{min}}^{\nu_\text{max}} \left( \frac{\nu^{\gamma + \beta/2}}{\nu_s^{\beta/2} \nu_0^\gamma} \right) \text{d}\nu $$

which evaluates to

$$ = 2 \pi^2 d^2 \sqrt{f(P,\dot{P})} \frac{\nu_\text{max}^{\gamma + \beta/2 + 1} - \nu_\text{min}^{\gamma + \beta/2 + 1}}{\gamma + \beta/2 + 1} \frac{S_m(\nu_0)}{\delta \nu_s^{\beta/2} \nu_0^\gamma}. $$

#### Derivation #3: Using a log-parabolic spectrum

The proposed log-parabolic spectrum is

$$ S_p(\nu) = S_p(\nu_0) \left( \frac{\nu}{\nu_0} \right)^\alpha e^{q \left[ \ln \left( \frac{\nu}{\nu_0} \right) \right]^2}, $$

where **$\alpha$ is now a spectral index** (not the magnetic inclination angle, as before), and the peak fluxes are used here (the case for $S_m(\nu)$ is treated later).
Since we're currently using $S_p(\nu)$, we don't need to factor the duty cycle into the luminosity calculation.
However, it should be understood that this will *overestimate* the (time-averaged) luminosity, since clearly the emission does not maintain the peak luminosity over the whole pulse window.

Still assuming a small $\rho$ (because of the large light cylinder and necessarily small polar cap), and also assuming the other empirical relations between $\rho$ and $\nu$ determined for normal pulsars (stated above), the luminosity is

$$ L = 4 \pi d^2 \int_{\nu_\text{min}}^{\nu_\text{max}} \frac{\rho(\nu)^2}{4} S_p(\nu) \text{d}\nu $$

$$ = 4 \pi d^2 f(P,\dot{P}) S_p(\nu_0) \left( \frac{\nu_0}{\nu_s} \right)^\beta \int_{\nu_\text{min}}^{\nu_\text{max}} \left( \frac{\nu}{\nu_0} \right)^{\alpha+\beta} e^{q \left[ \ln \left( \frac{\nu}{\nu_0} \right) \right]^2} \text{d}\nu $$

$$ = 4 \pi d^2 f(P,\dot{P}) S_p(\nu_0) \left( \frac{\nu_0}{\nu_s} \right)^\beta \left[ \frac12 \sqrt{ \frac{\pi}{q} } e^{-\frac{(\alpha+\beta+1)^2}{4q}} \text{erfi} \left( \frac{a + 2q \ln \left( \frac{\nu}{\nu_0} \right) + 1}{2\sqrt{q}} \right) \right]_{\nu_\text{min}}^{\nu_\text{max}} $$

$$ = 4 \pi d^2 f(P,\dot{P}) S_p(\nu_0) \left( \frac{\nu_0}{\nu_s} \right)^\beta \left[ \frac12 \sqrt{ \frac{\pi}{q} } e^{-\frac{(\alpha+\beta+1)^2}{4q}} \text{erfi} \left( \frac{a + 2q \ln \left( \frac{\nu}{\nu_0} \right) + 1}{2\sqrt{q}} \right) \right]^{\nu_\text{max}}\_{\nu_\text{min}} $$

