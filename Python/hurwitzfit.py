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

from math import sqrt, pi,nan,isnan
from cmath import exp
from cmath import sqrt as csqrt

import numpy as np
from scipy.optimize import curve_fit
import scipy.integrate as integrate
from scipy.special import bernoulli, poch
from math import factorial
import matplotlib.pyplot as plt

INFTY = 1e+100
sum_acc = 1e-6

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
			if( abs(add/sum) < sum_acc ): break
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
	if abs(V) < Vm:
		return np.sqrt(2*Vrms**2-V**2)/pi/Vrms**2
	return 0.0

# Convolution of Frota-Fano with lockin modulation
def ff_w_lockin( V, A0, DeltaK, phi ):
	integr = lambda Vp, A0, DeltaK, phi, V : Xi_mod(Vp) * frota_fano( V+Vp, A0, DeltaK, phi )
	# Changed Nov-30, 2023:
	# - Use Gaussian quadrature to speed-up convolution;
	# - also larger default tolerance (eps=1e-4) speeds-up integral considerably while still accurate enough
	#conv, _ = integrate.quad(integr, -Vm, Vm, args=(A0,DeltaK,phi,V))
	conv, _ = integrate.quadrature(integr, -Vm, Vm, args=(A0,DeltaK,phi,V), vec_func=False, tol=eps, rtol=eps)
	return conv

# Convolution of Hurwitz-Fano with lockin modulation
def hf_w_lockin( V, A0, DeltaK, phi ):
	integr = lambda Vp, A0, DeltaK, phi, V : Xi_mod(Vp) * hurwitz_fano( V+Vp, A0, DeltaK, phi )
	# Changed Nov-30, 2023:
	# - Use Gaussian quadrature to speed-up convolution;
	# - also larger default tolerance (eps=1e-4) speeds-up integral considerably while still accurate enough
	#conv, _ = integrate.quad(integr, Vm, Vm, args=(A0,DeltaK,phi,V))
	conv, _ = integrate.quadrature(integr, -Vm, Vm, args=(A0,DeltaK,phi,V), vec_func=False, tol=eps, rtol=eps)
	return conv

def fit_func( V, V0, A0, DeltaK, phi, a, b ):
	return vec_func( V-V0, A0, DeltaK, phi )+a*V+b


######################################
###                                ###
###          MAIN PROGRAM          ###
###                                ###
######################################

print()
print( "**************************************************************" )
print( "***                                                        ***" )
print( "***      hurwitzfit.py - fit to Hurwitz-Fano lineshape     ***" ) 
print( "***                                                        ***" )
print( "***  (c) 2023 by David Jacob, Universidad del Pais Vasco   ***" )
print( "***                                                        ***" )
print( "**************************************************************" )
print()

if len(sys.argv[1:]) < 3 :
	print( "Missing arguments.", file=sys.stderr )
	print( "Usage: hurwitzfit.py <fname> <col> <temp> [OPTIONS]", file=sys.stderr )
	print( "Mandatory arguments:", file=sys.stderr )
	print( "  <fname> : name of data file (string)", file=sys.stderr )
	print( "  <col>   : column number (starting at 1) for dI/dV in data file.", file=sys.stderr )
	print( "  <temp>  : tip temperature in K (real number)", file=sys.stderr )
	print( "Optional arguments:", file=sys.stderr )
	print( "  --frota          : Use Frota-Fano lineshape instead of Hurwitz-Fano.", file=sys.stderr )
	print( "  --lock-in=<Vrms> : Take into account lock-in modulation by numerical convolution.", file=sys.stderr )
	print( "                     <Vrms> is the root-mean-square bias of LI modulation in mV (float).", file=sys.stderr )
	print( "  --range=V1,V2    : Bias range for fitting with data: V1<=V<=V2.", file=sys.stderr )
	print( "                     If not specfied the entire data range will be used.", file=sys.stderr )
	print( "  --show           : whether to show plot after fit (uses matplotlib).", file=sys.stderr )
	print( "  --acc=<eps>      : Set accuracy for evalutaion of convolution integral to <eps>.", file=sys.stderr )
	print( "                     Default value (eps=1e-4) is usually accurate enough.", file=sys.stderr )
	print()
	exit(-1)

fname_full = sys.argv[1]
ncol = int(sys.argv[2])
T_str = sys.argv[3]

w_lockin = False
show_plot = False
use_frota = False
V1 = None
V2 = None

# NEW Nov-30:
# - Default tolerance for convolution integral:
#   can be much larger than the default values of quadrature routines (1.49e-8)
# - With this value convilution is much faster (by a factor of 25) but still accurate enough 
eps = 1e-4

# Iterate through remaining list of arguments
# to find optional arguments
for arg in sys.argv[4:]:
	#print( arg )
	if arg == '--frota':
		use_frota = True
	elif arg[0:10] == '--lock-in=':
		w_lockin = True
		Vrms = float(arg[10:])
	elif arg[0:8] == '--range=':
		V1,V2 = tuple(float(x) for x in arg[8:].split(","))
	elif arg == '--show':
		show_plot = True
	elif arg[0:6] == '--acc=':
		eps = float(arg[6:])
	else:
		print( "Unknown option: ", arg, " Abort.", file=sys.stderr )
		exit(-1)
		   
if use_frota:
	if w_lockin:
		vec_func = np.vectorize( ff_w_lockin, excluded=['A0','DeltaK','phi'] ) 
	else:
		vec_func = np.vectorize( frota_fano, excluded=['A0','DeltaK','phi'] )
else:
	zeta = hurwitz_zeta()
	if w_lockin:
		vec_func = np.vectorize( hf_w_lockin, excluded=['A0','DeltaK','phi'] )
	else:
		vec_func = np.vectorize( hurwitz_fano, excluded=['A0','DeltaK','phi'] )

fname = os.path.basename(fname_full)

print()

if not use_frota:
	print( "*** Fitting Hurwitz-Fano lineshape to dI/dV data ***" )
else:
	print( "*** Fitting Frota-Fano lineshape to dI/dV data ***" )
print()
print( " Data file: ", fname_full )
print( " dI/dV data in col #", ncol )
print( " T = ", T_str, "K" )
if V1 is not None and V2 is not None:
	print( " Fit range: = [", V1, ",", V2, "]" )
if w_lockin:
	print( " With lock-in modulation: Vrms = ", Vrms )
	Vm = np.sqrt(2)*Vrms
print()

kT = float(T_str) * 8.61732814974493E-02 # Temperature in meV
beta = 1.0/kT

data = np.loadtxt(fname_full)
ndata=np.size(data,0)
print( " Number of data points: ", ndata )

bias_data = np.zeros( ndata )
dIdV_data = np.zeros( ndata )

## if fitting range is specified, only take data within that range ###

if V1 is not None and V2 is not None:
	for i in range(ndata):	
		if data[i][0] >= V1 and data[i][0] <= V2:
			bias_data[i] = data[i][0]
			dIdV_data[i] = data[i][ncol-1]
		else:
			bias_data[i] = INFTY
			dIdV_data[i] = INFTY
	bias_data = np.extract(bias_data<INFTY,bias_data)
	dIdV_data = np.extract(dIdV_data<INFTY,dIdV_data)
	ndata = len(bias_data)
	print( " Number of data points in fit range: ", ndata )
else:
	for i in range(ndata):	
		bias_data[i] = data[i][0]
		dIdV_data[i] = data[i][ncol-1]

print()

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
	plt.plot( bias_data, dIdV_data, 'o', markersize=4 )
	plt.plot( bias_data, fit_func( bias_data, *popt ), '-', linewidth=2 )
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

if use_frota:
	ofname = "ff" # Frota-Fano
else:
	ofname = "hf" # Hurwitz-Fano

if w_lockin:
	ofname = ofname + "li" # with-lockin
	
ofname = ofname + "_" + fname
	
print( " Write fitted lineshape and dI/dV data to ",  ofname )
f = open( ofname, "w" )
for i in range(ndata):
	V = bias_data[i]
	print( 3*"%15.6e" % (V,fit_func(V,*popt),dIdV_data[i]), file=f )
f.close()
print()
print( "Done." )


