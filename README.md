# HurwitzFanoFit
Python program for fitting Hurwitz-Fano lineshapes to experimental data

## Theoretical background
The python script can be used to fit different lineshapes that occur in the theory of STM spectroscopy (STS) of Kondo systems to STS data.

### Frota-Fano lineshape
First, the Frota lineshape yields an accurate description of the Kondo resonance in the spectral function of a Kondo system. 
Taking into account quantum interference by a Fano phase factor $\phi$, the Kondo resonance can be well described by the
following Frota-Fano lineshape:
```math
A(\omega) = A_0 \, {\rm Re} \left[ e^{i\phi} \sqrt{\frac{i\Delta_{\rm K}}{\omega+i\Delta_{\rm K}}} \right] 
```
where $A_0$ is the amplitude of the Kondo resonance, $\phi\in[0,2\pi]$ is the Fano phase and $\Delta_{\rm K}$ is the Frota width parameter
which yields the halfwidth of the Kondo resonance as $\Gamma_{\rm K}=2.542\Delta_{\rm K}$.
When Fermi-Dirac smearing at the STM tip can be neglected ($T\rightarrow0$), the measured differential conductance ($dI/dV$) is directly proportional to $A(\omega$):
$\mathcal{G}(V) = dI/dV \propto A(eV)$.

### Hurwitz-Fano lineshape
Usually, smearing of the Fermi surface due to finite temperature $T$ at the STM tip cannot be neglected. 
In this case the lineshape of the Kondo resoance in the $dI/dV$ is given by a convolution with the derivative of the Fermi-Dirac distribution.
As was shown recently by us [1] the resulting lineshape can be described in terms of a Hurwitz $\zeta$-function:
```math
\mathcal{G}(V) = \frac{dI}{dV} \propto \int d\omega \, \left(-f^\prime(\omega)\right) \, A(\omega+eV) =
A_0 \, \sqrt{\frac{\beta \Delta_{\rm K}}{8\pi} } \,{\rm Re} \left[
    e^{i\phi} \zeta\left( \frac{3}{2}, \frac{\beta\Delta_{\rm K}}{2\pi} + \frac{1}{2} +i\frac{\beta eV}{2\pi} \right) 
\right]
```
where $\beta=1/kT$ the inverse temperature, $\Delta_{\rm K}$ is the Frota width parameter of the *underlying* Frota lineshape in the spectral function $A(\omega)$,
and $\phi$ the Fano phase factor describing quantum interference.
$\zeta(s,a)$ is the Hurwitz $\zeta$-function, a generalization of the Riemann $\zeta$-function,
defined by the infinite series:
```math
\zeta(s,a) = \sum_{n=0}^{\infty} \frac{1}{(n+a)^s}
```
The Riemann $\zeta$-function is recovered for $a=1$.

### Lock-in modulation
The effect of lock-modulation on the spectra in the STS experiment can be taken into account by an additional numerical convolution:
```math
\tilde{\mathcal{G}}(V) = \int dV^\prime \, \chi_{\rm LI}(V^\prime) \, \mathcal{G}(V+V^\prime)
```
where the lock-in modulation $\chi_{\rm LI}$ is given by
```math
\chi_{\rm LI}(V) = \left\{ \begin{array}{cc} 2\sqrt{2V_{\rm rms}^2 -V^2}/\pi V_{\rm rms}^2 & \mbox{ for } |V|\le \sqrt2 V_{\rm rms} \\
 0 & \mbox{ for } |V| > \sqrt2 V_{\rm rms} \end{array} \right.
```

## Installation

## Program usage

### Linux

### Windows

### MacOS
