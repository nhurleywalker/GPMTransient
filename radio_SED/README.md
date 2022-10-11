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

### Small opening angle approximation

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

$$ \frac{r}{10\text{km}} \approx 40 \left( \frac{\nu}{\text{GHz}} \right)^\beta \left( \frac{\dot{P}}{10^{-15}} \right)^{0.07} \left( \frac{P}{\text{s}} \right)^{0.30}, $$

(see Erkut's equation (14)), where $\beta \approx -0.26$ is the radius-to-frequency index, not to be confused with the impact angle above.
Hence,

$$ \rho \approx 1.^\circ24 \sqrt{40} 

#### Derivation #1: Constant duty cycle, single power law

Let

$$ f(P,\dot{P}) \approx $$
