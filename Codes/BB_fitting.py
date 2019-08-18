import numpy as np # numerical tools
from scipy import *
import matplotlib.pyplot as plt 
from matplotlib.pyplot import figure, show, rcParams
from scipy.optimize import curve_fit

################################################################################################
################################################################################################
################################################################################################
#
#          BB temperature from spectrum
#          Thomas de Jaeger
#          2019
#
################################################################################################
################################################################################################
################################################################################################

########################################## Constantes #######################
class C:
	'''
	A container class that holds a bunch of constants.
	'''
	##############################
	# constants
	##############################
	c = 2.998E10 #cm s^-1
	sig_B = 5.67E-5 #erg cm^-2 s^-1 K^-4
	a_B = 7.56E-15 #erg cm^-3 K^-4
	k = 1.38E-16 #erg K^-1
	wein_lam = .29 #cm K^-1
	wein_nu = 5.88E10 #Hz K^-1
	h = 6.6260755E-27 #erg s
	h_bar = 1.05457266E-27 #erg s
	G = 6.674E-8 #cm^3 g^-1 s^-2
	sig_T = 6.65E-25 #cm^2
	pi = 3.141592653589793 #none
	H0_h = 3.2408E-18 #h * s^-1
	H0 = 2E-18 #s^-1
	T_cmb = 2.725 #K

	##############################
	# properties
	##############################
	m_p = 1.67E-24 #g
	m_e = 9.11E-28 #g
	M_sun = 1.988E33 #g
	M_earth = 5.972E27 #g
	R_sun = 6.955E10 #cm
	R_earth = 6.3675E8 #cm
	L_sun = 3.846E33 #erg s^-1
	e = 4.803E-10 #statC
	T_sun = 5778. # K, surface

	##############################
	# conversions
	##############################
	eV2erg = 1.602E-12 #ergs eV^-1
	year2sec = 3.154E7 #seconds yr^-1
	pc2cm = 3.086E18 #cm parsec^-1

rcParams['legend.numpoints']=1

##############################
# functions BB
##############################

def black_body_nu(nu, T):
	''' blackbody curve as a function of frequency (array) and T (float) '''
	B = (2.*C.h*nu**3)/(C.c**2) * ( np.exp((C.h*nu)/(C.k*T)) - 1.)**-1
	return B

def black_body_lam(lam, T):
	''' blackbody curve as function of wavelength (angstroms, array) and T (float) '''
	# convert lambda (in angstroms) to cm
	lam = 1E-8*lam
	B = (2.*C.h*C.c)/(lam**5) * ( np.exp((C.h*C.c)/(lam*C.k*T)) -1.)**-1
	return B

def bb_fit( x, T, Factor ):
	'''used in fit_blackbody, this function returns a blackbody of
	temperature T, evaluated at wavelengths x (Angstroms), and
	rescaled by a factor Factor
	'''
	return Factor * black_body_lam( x, T )

######################################################################3

##############################
# functions fit BB
##############################

def fit_blackbody(file, interactive=True, guessT=10000, plot=True):
	'''
	Fit a blackbody curve to an input spectrum.
	file= text file with 2 columns
	lam: wavelength  in Angs
	flux: flux in erg/s/cm2/Angs
	'''
	spectra_SN=np.loadtxt(file).transpose()
	lam=spectra_SN[0] 
	flux=spectra_SN[1] 
	x = np.array(lam)
	y = np.array(flux)

	if interactive:
		plt.ion()
		plt.figure()
		plt.show()
		plt.clf()
		plt.plot( x, y, alpha=.75 )
		plt.title('Click to define points')
		plt.draw()
		print ("Click on points to which you'd like to fit the BB curve.")
		print (" Left click: choose a point")
		print (" Right click: remove most recent point")
		print (" Middle click: done")
		points = plt.ginput(n=0, timeout=0)
		xPoints, yPoints = map( np.array, zip(*points) )

		# pick a first-guess scale factor
		guessA = yPoints[0]/black_body_lam( xPoints[0], guessT )

		[T,A], pcov = curve_fit(bb_fit, xPoints, yPoints, p0=[guessT, guessA], maxfev=100000)

		plt.ioff()
		plt.close()

	else:
		# use an average value to get a best-guess scale factor
		l = len(x)
		xavg = np.mean(x[int(.25*l):int(.25*l)+10])
		yavg = np.mean(y[int(.25*l):int(.25*l)+10])
		guessA = yavg/black_body_lam( xavg, guessT )

		[T,A], pcov = curve_fit(bb_fit, x, y, p0=[guessT, guessA], maxfev=100000)


	if plot:
		plt.plot( x, y,'k',label='Spectrum' )
		plt.plot( x, bb_fit(x, T, A), 'r',label='BB fit %s K'%int(T))
		plt.title( 'Blackbody fitting' )
		plt.ylabel('$F_\lambda$ [erg s$^{-1}$ cm$^{-2}$ $\AA$$^{-1}$]')
		plt.xlabel('Wavelength in [$\AA$]')
		plt.legend(loc=0,markerscale=0.5,prop={'size':8},ncol=1)
		plt.show()
	print ('Best-fit T:', int(T))

	return T, A


#Example
file='Data/SN2016esw_spec.dat'#input('give the file name: ')
fit_blackbody(file, interactive=True, guessT=10000, plot=True)


