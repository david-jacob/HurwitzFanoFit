# *************************
# ***                   ***
# ***   hurwitzfit.py   ***
# ***                   ***
# *************************
#
# Python program for fitting Hurwitz-Fano / Frota-Fano lineshapes to experimental data 
# 
# (c) 2023 by David Jacob, Universidad del Pais Vasco
#

import sys
import os
import time

from math import sqrt, pi
from cmath import exp
from cmath import sqrt as csqrt

import numpy as np
from scipy.optimize import curve_fit
import scipy.integrate as integrate
from scipy.special import bernoulli, poch
from math import factorial
import matplotlib.pyplot as plt

# ----------------------------------------
# Class implementing Hurwitz zeta-function
# ----------------------------------------
# - makes use of Euler-Maclaurin formula
class hurwitz_zeta:
	def __init__(self,n=20,r=20):
		self.n = n
		self.r = r
		self.B = bernoulli(r)
	# Hurwitz zeta function takes two arguments: s,a
	def __call__(self,s,a):
		n=self.n
		r=self.r
		sum = ((n+a)**(1-s))/(s-1)
		for k in range(n+1):
			sum += 1.0/(k+a)**s
		sum -= 0.5/(a+n)**s
		for kk in range(1,r//2):
			add = poch(2*kk-2+s,2*kk-1)*self.B[2*kk]/(factorial(2*kk)*(a+n)**(2*kk-1+s))
			if( abs(add/sum) < 1e-6 ): break
			sum += add
		return sum
	
# --------------------------------------------
# Function implementing Hurwitz-Fano lineshape
# --------------------------------------------
# Arguments:
# - V (real):      bias voltage
# - A0 (real):     amplitude
# - DeltaK (real): Frota width parameter
# - phi (real):    Fano phase
#
def hurwitz_fano( V, A0, DeltaK, phi ):
	return A0 * 0.5*sqrt(0.5*beta*DeltaK/pi) * ( exp( 1.j*phi ) * zeta(1.5, 0.5*(beta*DeltaK/pi+1.0+1.j*beta*V/pi))).real

# ------------------------------------------
# Function implementing Frota-Fano lineshape
# ------------------------------------------
# Arguments:
# - V (real):      bias voltage
# - A0 (real):     amplitude
# - DeltaK (real): Frota width parameter
# - phi (real):    Fano phase
#
def frota_fano( V, A0, DeltaK, phi ):
	return A0 * ( exp(1.j*phi)/np.sqrt( 1.0 + 1.j*V/DeltaK ) ).real

# lock-in modulation function
def Xi_mod( V ):
	if abs(V) < np.sqrt(2)*Vrms:
		return np.sqrt(2*Vrms**2-V**2)/pi/Vrms**2
	return 0.0

# Convolution of Frota-Fano with lockin modulation
def ff_w_lockin( V, A0, DeltaK, phi ):
	integr = lambda Vp, A0, DeltaK, phi, V : Xi_mod(Vp) * frota_fano( V+Vp, A0, DeltaK, phi )
	conv, _ = integrate.quad(integr, -np.inf, np.inf, args=(A0,DeltaK,phi,V))
	return conv

# Convolution of Hurwitz-Fano with lockin modulation
def hf_w_lockin( V, A0, DeltaK, phi ):
	integr = lambda Vp, A0, DeltaK, phi, V : Xi_mod(Vp) * hurwitz_fano( V+Vp, A0, DeltaK, phi )
	conv, _ = integrate.quad(integr, -np.inf, np.inf, args=(A0,DeltaK,phi,V))
	return conv

def fit_func( V, V0, A0, DeltaK, phi, a, b ):
	return vec_func( V-V0, A0, DeltaK, phi )+a*V+b


### MAIN PROGRAM ###

print()
print( "**************************************************************" )
print( "***                                                        ***" )
print( "***    hurwitzfit.py - fit to Hurwitz / Frota lineshape    ***" ) 
print( "***                                                        ***" )
print( "***  (c) 2023 by David Jacob, Universidad del Pais Vasco   ***" )
print( "***                                                        ***" )
print( "**************************************************************" )
print()

if len(sys.argv[1:]) < 4 :
	print( "Missing arguments.", file=sys.stderr )
	print( "Usage: hurwitzfit.py <fname> <col> <temp> <func> [Vrms] [show]", file=sys.stderr )
	print( "  <fname> : name of data file (string)" )
	print( "  <col>   : column for dI/dV in data file (0, 1, 2, ...)" )
	print( "  <temp>  : tip temperature in K (real number)" )
	print( "  <func>  : lineshape to fit (string) - 'hurwitz' or 'frota'" )
	print( "  [Vrms]  : root-mean-square bias of lock-in modulation in mV (float) [optional]" )
	print( "  [show]  : whether to show plot after fit [optional]." )
	print()
	exit(-1)

fname_full = sys.argv[1]
ncol = int(sys.argv[2])
T_str = sys.argv[3]
func_name = sys.argv[4]

w_lockin = False
if len(sys.argv[1:]) >= 5 :
	Vrms = float(sys.argv[5])
	if abs(Vrms)>1e-6:
		w_lockin = True

show_plot = False
if len(sys.argv[1:]) >= 6 :
	if sys.argv[6] == 'show': show_plot = True
		   
if func_name == "hurwitz":
	zeta = hurwitz_zeta()
	# vectorize function for use in _curve_fit
	if w_lockin:
		vec_func = np.vectorize( hf_w_lockin, excluded=['A0','DeltaK','phi'] )
	else:
		vec_func = np.vectorize( hurwitz_fano, excluded=['A0','DeltaK','phi'] )
elif func_name == "frota":
	if w_lockin:
		vec_func = np.vectorize( ff_w_lockin, excluded=['A0','DeltaK','phi'] ) 
	else:
		vec_func = np.vectorize( frota_fano, excluded=['A0','DeltaK','phi'] ) 
else:
	print( "Error: unknown function name: ", func_name, file=sys.stderr )
	print( " Possible choices: 'frota', 'hurwitz`", file=sys.stderr )
	exit(-1)

fname = os.path.basename(fname_full)

print()
print( "*** Fitting dI/dV data to Hurwitz zeta-function ***" )
print()
print( " Data file: ", fname_full )
print( " dI/dV data in col #", ncol )
print( " T = ", T_str, "K" )
print()

kT = float(T_str) * 8.61732814974493E-02 # Temperature in meV
beta = 1.0/kT

data = np.loadtxt(fname_full)
ndata=np.size(data,0)
print( " Number of data points: ", ndata )
print()

bias_data = np.zeros( ndata )
dIdV_data = np.zeros( ndata )

for i in range(ndata):
	bias_data[i] = data[i][0]
	dIdV_data[i] = data[i][ncol]

#
# Curve fitting with initial guess and lower and upper boundaries for fit parameters
#
# V0     : peak position (0.01mV) [-20mV:+20mV]
# A0     : Peak amplitude (1.0) [0:100]
# DeltaK : Kondo peak width parameter (0.5mV) [0:200] - HWHM ~ 2.542 * DeltaK
# phi    : Fano phase (0.1) [-pi:+pi] 
# a,b    : linear background ~ a*V + b (0.01,0.1) [-infty:+infty]
#
start = time.time()
popt, pcov = curve_fit(fit_func, bias_data, dIdV_data, p0=[0.01,1.0,0.5,0.1,0.01,0.1],
					   bounds=([-20.,0.0,0.0,-pi,-np.inf,-np.inf],[20.,100.,200.,pi,np.inf,np.inf]))
end = time.time()
print( "*** Fit converged ***" )
print()

print( " Elapsed time: ", end-start, "s." ) 

print()
print( " Fit parameters:" )
print( " V0 = ", popt[0] )
print( " A0 = ", popt[1] )
print( " DeltaK = ", popt[2] )
print( " phi = ", popt[3] )
print( " a =", popt[4] )
print( " b =", popt[5] )
print()
GammaK = 2.542 * popt[2]
print( " GammaK = ", GammaK )
print()


if show_plot:
	plt.plot( bias_data, dIdV_data, 'bo' )
	plt.plot( bias_data, fit_func( bias_data, *popt ), 'r-' )
	plt.show()

print( " Writing temperature and all fit parameters to ",  "optparams.dat" ) 
f = open( "optparams.dat", "a" )
print( "%8s" % T_str,file=f, end="  ")
print( 6*"%15.6e" % tuple(popt), file=f, end="  #" )
print( fname, file=f )
f.close()

print( " Writing temperature and GammaK to ", "GammaK_vs_temp.dat" ) 
f = open( "GammaK_vs_temp.dat", "a" )
print( "%8s" % T_str, file=f, end="  ")
print( "%15.6e" % GammaK, file=f, end="  #" )
print( fname, file=f )
f.close()

# Output file name for fitted lineshape:

if func_name == 'hurwitz':
	ofname = "hf" # Hurwitz-Fano
else:
	ofname = "ff" # Frota-Fano
if w_lockin:
	ofname = ofname + "wl" # with-lockin
	
ofname = ofname + "_" + fname
	
print( " Write fitted Hurwitz zeta-function to ",  ofname )
f = open( ofname, "w" )
for i in range(ndata):
	V = bias_data[i]
	print( 3*"%15.6e" % (V,fit_func(V,*popt),dIdV_data[i]), file=f )
f.close()
print()
print( "Done." )


