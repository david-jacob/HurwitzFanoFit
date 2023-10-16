#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3				// Use modern global access method and strict wave access
#pragma DefaultTab={3,20,4}		// Set default tab width in Igor Pro 9 and later

// ######################################################################
// #																	#
// #		Fit Functions for Kondo Resonances in STS					#
// #					(c) by Nils Krane								#
// #				Empa Dübendorf, Switzerland							#
// #																	#
// ######################################################################


// Version 1.0 (2023-10-16)


// The Following Functions are included in the KondoFunctions.ipf

//		FrotaFanoLineshape: Returns Frota-Fano lineshape
//		FrotaFit: Fit Function for Frota-Fano lineshape with linear background

//		HurwitzFanoLineshape: Returns Hurwitz-Fano lineshape (aka Sobrino lineshape)
//		HurwitzFit: Fit Function for Hurwitz-Fano lineshape with linear background
//		HurwitzLockinFit: Fit Function for lock-in-broadened Hurwitz-Fano lineshape

//		GammaT_Jacob: Fit Function for the temperature dependent broadening of intrinsic Kondo halfwidht "Γ" (aka "HWHM") 
//		TK2Gamma: Returns Kondo halfwidth "Γ" (aka "HWHM") for given experimental temperature T and Kondo temperature TK
//		Gamma2TK: Returns Kondo temperature TK for given Kondo halfwidht "Γ" (aka "HWHM") and experimental temperature T


// ################################################# KONDO LINESHAPES #################################################


// Return y-Value of Frota-Fano lineshape at given x-value (xi)
Function FrotaFanoLineshape(xi,HWHM,phi)
	Variable xi		// given x-value
	Variable HWHM	// linewidth (aka Γ = 2.542Δ)
	Variable phi	// Phano Phase
	
	Variable Delta=HWHM/2.542 // Kondo width parameter 
	
	return real(exp(cmplx(0,phi))*sqrt((cmplx(0,Delta))/(xi+cmplx(0,Delta))))
End


// -------------------------------------------------------------------------------------------------------


// Fit Function for Frota-Fano lineshape
// including linear background
Function FrotaFit(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = Amplitude*FrotaFanoLineshape(x-x_offset],HWHM,Phi)+y_offset+slope*x
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 6
	//CurveFitDialog/ w[0] = x_offset
	//CurveFitDialog/ w[1] = HWHM
	//CurveFitDialog/ w[2] = Phi
	//CurveFitDialog/ w[3] = Amplitude
	//CurveFitDialog/ w[4] = y_offset
	//CurveFitDialog/ w[5] = slope

	return w[3]*FrotaFanoLineshape(x-w[0],w[1],w[2])+w[4]+w[5]*x
End


// ==================================================================================================================


// Return y-value of Hurwitz-Fano lineshape at given x-value and Temperature T
// as described by Eq. (3) in E. Turco et al. (arXiv:TBA)
// The Hurwitz-Fano (or Sobrino) lineshape corresponds to a Frota-Fano lineshape
// broadened by Fermi-Dirac smearing 
Function HurwitzFanoLineshape(xi,HWHM,Phi,T)
	Variable xi		// given x-value
	Variable HWHM	// linewidth (aka Γ = 2.542Δ)
	Variable Phi	// Phano Phase
	Variable T		// Temperature
	
	Variable tau=pi*T*86.173E-6	// temperature parameter
	Variable delta=HWHM/2.542	// Kondo width parameter
	
	return sqrt(delta/tau/8) * real(exp(cmplx(0,phi))*zeta(3/2,1/2+delta/tau/2+cmplx(0,xi/tau/2)))

End


// -------------------------------------------------------------------------------------------------------


// Fit Function for Hurwitz-Fano lineshape
// including linear background
// Temperature T is the known experimental Temperature and a fixed parameter ("Hold": yes)
Function HurwitzFit(w,x) : FitFunc
	Wave w
	Variable x

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(x) = Amplitude*HurwitzFanoLineshape(x-x_offset,Phi,HWHM,T)+y_offset+slope*x
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ x
	//CurveFitDialog/ Coefficients 7
	//CurveFitDialog/ w[0] = x_offset
	//CurveFitDialog/ w[1] = HWHM
	//CurveFitDialog/ w[2] = Phi
	//CurveFitDialog/ w[3] = Amplitude
	//CurveFitDialog/ w[4] = y_offset
	//CurveFitDialog/ w[5] = slope
	//CurveFitDialog/ w[6] = T

	return w[3]*HurwitzFanoLineshape(x-w[0],w[1],w[2],w[6])+w[4]+w[5]*x
End


// -------------------------------------------------------------------------------------------------------


// Fit Function for lock-in-broadened Hurwitz-Fano lineshape
// including linear background
// Temperature T, Lock-In Modulation Vmod are fixed parameter ("Hold": yes)
// dX is the point density for convolution and...
//		... should be smaller than point density of spectrum (e.g. x1/10)
//		... is fixed parameter ("Hold": yes)
Function HurwitzLockinFit(pw, yw, xw) : FitFunc
	Wave pw, yw, xw
	
	// pw[0] = x0		x-offset
	// pw[1] = gam		HWHM of FrotaFunction	
	// pw[2] = phi		Fano-Phase in units of pi
	// pw[3] = A		Amplitude of Peak
	// pw[4] = yoff		y-offset
	// pw[5] = slope	background slope
	// pw[6] = T		Temperature
	// pw[7] = Vmod		modulation amplitude of LockIn
	// pw[8] = dX		point density for convolution (e.g. 10E-6 V)
	

	Variable Vmod=abs(pw[7])
	Variable dX=abs(pw[8])

	// xw seems to be of increasing order
	Variable xmin=xw[0]
	Variable xrange=xw[numpnts(xw)-1]-xmin

	// Determine how many points the lockin wave needs
	Variable VmodPnts= Vmod ? ceil(Vmod/dX) : 0	// how many points per amplitude?
	Variable nYPnts=round(xrange/dX)+2*VmodPnts	// #points for spectrum + 2*amplitude
	Make/FREE/D/N=(nYPnts) yWave
	SetScale/P x, xmin-VmodPnts*dX, dX, yWave
	
	yWave=HurwitzFanoLineshape(x-pw[0],pw[1],pw[2],pw[6])	// HF lineshape without lockin

	// if Vmod>0: convolute HF lineshape with lockIn broadening
	if(VmodPnts)
		// create wave for lockin broadening function chi	
		Make/FREE/D/N=(VmodPnts*2+1) HFLockIn
		SetScale/P x, -VmodPnts*dX, dX, HFLockIn
		HFLockIn= abs(x)<Vmod ? dX*2*sqrt(Vmod^2-x^2)/pi/Vmod^2 : 0
		
		Convolve/A HFLockIn, yWave
	endif 
	
	// add linear background and return
	yw=ywave(xw[p])*pw[3]+pw[4]+pw[5]*xw[p]

End


// ################################################# KONDO TEMPERATURE #################################################


// Fit Function for temperature-dependence of intrinsic Kondo halfwidth Γ(T)
// according to Eq. (10) from David Jacob (arXiv:2306.13136)
// or Eq. (1) form E. Turco et al. (arXiv: TBA) 
// using Wilson's Definition of TK: Δ_0 = 1.542*TK*kB and Γ_0 = 3.920*TK*kB
//
// TK: Wilson’s thermodynamic definition of the Kondo temperature [1], corrected by Wiegman and Tsvelick [2]
// [1]: K. G. Wilson, Rev. Mod. Phys. 47, 773 (1975)
// [2]: P. B. Wiegmann and A. M. Tsvelick, Journal of Physics C: Solid State Physics 16, 2281 (1983)
Function GammaT_Jacob(w,T) : FitFunc
	Wave w
	Variable T

	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(T)=TK2Gamma(TK,T)
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ T
	//CurveFitDialog/ Coefficients 1
	//CurveFitDialog/ w[0] = TK

	return TK2Gamma(w[0],T)
End


// -------------------------------------------------------------------------------------------------------


// Return intrinsic Kondo halfwidth "Γ" (aka "HWHM")
// for a Kondo system with Kondo temperature "TK" at finite temperature "T"
// according to Eq. (10) from David Jacob (arXiv:2306.13136)
// or Eq. (1) from E. Turco (arXiv: TBA) 
// using Wilson's definition of TK: Δ_0 = 1.542*TK*kB and Γ_0 = 3.920*TK*kB
//
// TK: Wilson’s thermodynamic definition of the Kondo temperature [1], corrected by Wiegman and Tsvelick [2]
// [1]: K. G. Wilson, Rev. Mod. Phys. 47, 773 (1975)
// [2]: P. B. Wiegmann and A. M. Tsvelick, Journal of Physics C: Solid State Physics 16, 2281 (1983)
Function TK2Gamma(TK,T)
	Variable TK	// Kondo Temperature (TK = 0.255*HWHM/kB)
	Variable T	// Experimental Temperature
	
	Variable kB=86.17E-6 // Boltzmann constant in eV/K
	Variable a=1+sqrt(3)
	Variable b=2+sqrt(3)
	Variable c=sqrt(3)/2
	
	return TK*1.542*kB * sqrt(a+b*sqrt(1+(pi*T/(TK*1.542))^2)+c*(pi*T/(TK*1.542))^2)
End


// -------------------------------------------------------------------------------------------------------


// Return Kondo temperature "TK"
// for a Kondo system with intrinsic Kondo halfwidth "Γ" (aka "HWHM")
// and experimental temperature "T"
// according to Eq. (5) from E. Turco et al. (arXiv:TBA) 
// using Wilson's definition of TK: Δ_0 = 1.542*TK*kB and Γ_0 = 3.920*TK*kB
//
// TK: Wilson’s thermodynamic definition of the Kondo temperature [1], corrected by Wiegman and Tsvelick [2]
// [1]: K. G. Wilson, Rev. Mod. Phys. 47, 773 (1975)
// [2]: P. B. Wiegmann and A. M. Tsvelick, Journal of Physics C: Solid State Physics 16, 2281 (1983)
Function Gamma2TK(HWHM,T)
	Variable HWHM	// intrinsic Kondo halfwidth "Γ"
	Variable T		// Experimental Temperature
	
	Variable a=(2+sqrt(3))/6
	Variable b=1-1/sqrt(12)
	Variable c=1-1/sqrt(3)
	Variable d=pi*T*86.17E-6

	return sqrt(sqrt(HWHM^4/3+d^2*HWHM^2/3+a*d^4)-b*d^2-c*HWHM^2)/1.542/86.17E-6
End


// #####################################################################################################################
