# HurwitzFanoFit
Python program for fitting Hurwitz-Fano lineshapes to experimental data

## Theoretical background
The python script can be used to fit different lineshapes that occur in the theory of STM spectroscopy (STS) of Kondo systems to STS data.

### Frota-Fano lineshape
First, the Frota lineshape yields an accurate description of the Kondo resonance in the spectral function of a Kondo system.
```math
A(\omega) = A_0 \, {\rm Re} \left[ e^{i\phi} \sqrt{\frac{i\Delta_{\rm K}}{\omega+i\Delta_{\rm K}}} \right] 
```
where $A_0$ is the amplitude of the Kondo resonance and $\Delta_{\rm K}$ is the Frota width parameter
which yields the halfwidth of the Kondo resonance as $\Gamma_{\rm K}=2.542\Delta_{\rm K}$.
When Fermi-Dirac smearing at the STM tip can be neglected ($T\rightarrow0$), the measured differential conductance ($dI/dV$) is directly proportional to $A(\omega$):
$\mathcal{G}(V) = dI/dV \propto A(eV)$.

### Hurwitz-Fano lineshape
Usually, smearing of the Fermi surface due to finite temperature at the STM tip cannot be neglected. 
In this case the lineshape of the Kondo resoance in the $dI/dV$ is given by a convolution with the derivative of the Fermi-Dirac distribution.
As was shown recently by us [1] the resulting lineshape can be described in terms of a Hurwitz $\zeta$-function:
```math
\mathcal{G}(V) = \frac{dI}{dV} \propto \int d\omega \, \left(-f^\prime(\omega)\right) \, A(\omega+eV) =
A_0 \, \sqrt{\frac{\beta \Delta_{\rm K}}{8\pi} } \,{\rm Re} \left[
    e^{i\phi} \zeta\left( \frac{3}{2}, \frac{\beta\Delta_{\rm K}}{2\pi} + \frac{1}{2} +i\frac{\beta eV}{2\pi} \right) 
\right]
```

### Lock-in modulation

## Installation

## Program usage

### Linux

### Windows

### MacOS
