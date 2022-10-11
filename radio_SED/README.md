# Radio SEDs and Luminosity

## Luminosity calculations following [Erkut (2022)](https://ui.adsabs.harvard.edu/abs/2022MNRAS.514L..41E%2F/abstract)

In their Equation (15), the general expression for the luminosity over the frequency range $\nu_\text{min} \le \nu \le \nu_\text{max}$ for a beamed pulsar-like object is

$$ L = 4\pi d^2 \int_{\nu_\text{min}}^{\nu_\text{max}} \sin^2 \left[ \frac{\rho(\nu)}{2} \right] S_p(\nu) \text{d}\nu,$$

where $d$ is the distance, $\rho$ is the opening angle of the beam, and $S_p(\nu)$ is the peak flux density.
This assumes that the beam is uniformly illuminated with flux density $S_p(\nu)$.

In the case of GLEAM-X J162759.5–523504.3, using the peak flux density as a proxy for the emission beam is problematic, owing to the paucity of pulses and the large uncertainty of the average pulse shape.
Hence, it is arguably more reliable (and practical, for continuum imaging!) to use the mean flux density, $S_m(\nu)$, but doing so requires knowing (or assuming) the viewing geometry.

The general relationship between the opening angle and the duty cycle, from spherical geometry, is

$$ \sin^2\left(\frac{\pi}{2}\delta\right) = \frac{\sin^2(\rho/2) - \sin^2(\beta/2)}{\sin\alpha \sin(\alpha + \beta)},$$

where $\alpha$ is the star's magnetic inclination angle (relative to the spin axis), and $\beta$ is the impact angle (the closest approach of the magnetic axis to the line of sight).
Since Erkut's argument is that the luminosity has been overestimated, we can assume a "worst case scenario" that the impact angle is very small, since as the impact angle increases, so too must the opening angle in order to match the observed duty cycle.
Thus, setting $\beta = 0$,

$$ \sin\left(\frac{\pi}{2}\delta\right) = \frac{\sin(\rho/2)}{\sin\alpha},$$

Erkut's stated relationship between the mean flux density and the peak, $S_p = S_m/\delta$, where $\delta$ is the duty cycle, implicitly assumes a very special viewing geometry.
