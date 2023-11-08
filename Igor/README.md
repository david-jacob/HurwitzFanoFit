# Igor program - KondoFunctions.ipf
Igor routines for fitting Hurwitz-Fano lineshapes to experimental data. 

## Installation
The Igor program can be installed by simply downloading the procedure file `KondoFunctions.ipf` from GitHub and copying it into the `Igor Procedures` or `User Procedures` folder.

For a single experiment the procedure file could also be loaded via drag and drop into the experiment or via the Menu `File` $\rightarrow$ `Open File` $\rightarrow$ `Procedure...`.

## System requirements
WaveMetrics Igor, tested with version 7, 8 and 9.

## Program usage
The following functions are included in the `KondoFunctions.ipf`

### Frota-Fano Lineshape

* `FrotaFanoLineshape(x,HWHM,Phi)`
Returns for a given `x` the value of a Frota-Fano lineshape with $\Gamma=2.542\Delta$ (`HWHM`) and Fano phase factor `Phi` 

* `FrotaFit(w,x) : FitFunc`
  Fit function for Frota-Fano lineshape with linear background. It can be used via the GUI in `Analysis` $\rightarrow$ `Curve Fitting...` or via the `FuncFit` command.

### Hurwitz-Fano Lineshape

* `HurwitzFanoLineshape(x,HWHM,Phi,T)`
  Returns for a given `x` the value of a Hurwitz-Fano lineshape with $\Gamma=2.542\Delta$ (`HWHM`), Fano phase factor `Phi` and temperature `T`.

* `HurwitzFit(w,x) : FitFunc`
  Fit function for Hurwitz-Fano lineshape with linear background. It can be used via the GUI in `Analysis` $\rightarrow$ `Curve Fitting...` or via the `FuncFit` command.
  The temperature `T` must be set to `Hold` by checking the corresponding checkbox in the GUI or using `FuncFit` with the `/H="0000001"` flag.

* `HurwitzLockinFit(pw, yw, xw) : FitFunc`
  All-in-one fit function for a lock-in-broadened Hurwitz-Fano lineshape with linear background. It can be used via the GUI in `Analysis` $\rightarrow$ `Curve Fitting...` or via the `FuncFit` command.
  The input parameters for temperature `T`, Lock-In modulation `Vmod` and convolution density `dX` must be set to `Hold` by checking the corresponding checkbox in the GUI or using `FuncFit` with the `/H="000000111"` flag.
  The convolution density `dX` should be chosen small enough for good results, e.g. $dX \leq V_\rm{mod}/10$.

  For Igor Version 7 or earlier: In order to use an All-in-one fit function in the `Curve Fitting` dialog, a coefficient wave with correct length has to be created beforehand, *e.g.*
  ```code
  Make/O/N=(9) W_coef
  ```
  for the `HurwitzLockinFit` fit function. Later versions of Igor will create this coefficient wave automatically.

### Kondo Temperature

* `Gamma2TK(HWHM,T)`
  Returns the Kondo temperature according to Wilson's definition for a halfwidth `HWHM` of a Kondo resonance at temperature `T` (see equation (5) in [1]).

*  `TK2Gamma(HWHM,T)`
  Returns the halfwidth of a Kondo resonance with Kondo temperature `TK` (Wilson's definition) at a given temperature `T` (see equation (1) in [1]).

* `GammaT_Jacob(w,T) : FitFunc`
  Fit function wrapper for the `TK2Gamma` function, to be used for fitting halfwidth versus temperature data sets (*e.g.* Fig3b in [1]).

## Example

The `example.ibw` file contains a Kondo resonance measured at T=4.45K and Vmod=0.4mV.
To fit the data one can either use the `Curve fitting` GUI or the commandline, *e.g.* with the Lockin-broadened Hurwitz-Fano function
```code
• Make/O W_coef = {0,0.0035,0,40,10,0,4.45,0.0004,10E-6}
• FuncFit/H="000000111"/TBOX=768 HurwitzLockinFit W_coef example /D
  Fit converged properly
  Duplicate/O fit_example,WMCF_TempAutoXWave
  HurwitzLockinFit(W_coef,fit_example,WMCF_TempAutoXWave)
  KillWaves/Z WMCF_TempAutoXWave
  W_coef={0.00011517,0.0037781,0.15616,41.37,9.0139,-15.521,4.45,0.0004,1e-05}
  V_chisq= 95.1896;V_npnts= 321;V_numNaNs= 0;V_numINFs= 0;
  V_startRow= 0;V_endRow= 320;
  W_sigma={2.6e-05,5.71e-05,0.00891,0.185,0.0966,5.86,0,0,0}
  Coefficient values ± one standard deviation
  	pw_0	= 0.00011517 ± 2.6e-05
  	pw_1	= 0.0037781 ± 5.71e-05
  	pw_2	= 0.15616 ± 0.00891
  	pw_3	= 41.37 ± 0.185
  	pw_4	= 9.0139 ± 0.0966
  	pw_5	= -15.521 ± 5.86
  	pw_6	= 4.45 ± 0
  	pw_7	= 0.0004 ± 0
  	pw_8	= 1e-05 ± 0
```
<p align="center"><img class="marginauto" src="Figure_example_HurwitzLockinFit.png" width="400"></p>

The obtained halfwidth is given by the parameter `pw_1` with 3.7781mV. Using the `Gamma2TK` function yields directly the Kondo temperature of the system:
```code
• print Gamma2TK(3.7781E-3,4.45)
  9.69676
```


## References
[1] E. Turco, M. Aapro, S. C. Ganguli, N. Krane, R. Drost, N. Sobrino, A. Bernhardt, M. Juríček, R. Fasel, P. Ruffieux, P. Liljeroth, D. Jacob, 
"Accurate Kondo temperature determination of spin-1/2 magnetic impurities", arXiv:2310.09326 (2023). https://doi.org/10.48550/arXiv.2310.09326

