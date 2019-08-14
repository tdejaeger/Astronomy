import re # use regular patters
import sys # system commands
import string as string # string functions4
import math
import numpy as np # numerical tools
from scipy import *
from pylab import *
import pyfits
import os
from scipy import integrate
from scipy import interpolate
import itertools
from scipy.optimize import leastsq
import math as maths
from matplotlib.pyplot import figure, show, rc
from scipy.optimize import curve_fit


################################################################################################
################################################################################################
################################################################################################
#
#          Permet de calculer la temperature BB
#          Thomas de Jaeger
#          2013
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
c=299792.458#in km/s
c_m=299792458.0#in km/s
c_mario=2.99792458*1.0e10#in cm/s
c_AA=299792458*1.0e10#in AA/s
kB=1.38066e-23; # Boltzmanns constant in J/K
h_planck = 6.63e-34 #Plancks constant (J.s)
h_erg = 6.63e-27 #Plancks constant (erg.s)
h=0.71
H0 = 100.0 * h  #Hubble constant in Km/s/Mpc
omega_m=0.27   # Dark matter density for standard model
omega_k=0.     #Universe curvature for standard model
omega_lambda=0.73 #Dark energy density for standard model
d_h=c*1.0/H0 #Hubble distance in Mpc

### Definition fonctions utiles

def find_closest(A, target):
    #A must be sorted
    idx = A.searchsorted(target)
    idx = np.clip(idx, 1, len(A)-1)
    left = A[idx-1]
    right = A[idx]
    idx -= target - left < right - target
    return idx

def date_to_jd(year,month,day):

    if month == 1 or month == 2:
        yearp = year - 1
        monthp = month + 12
    else:
        yearp = year
        monthp = month
    
    # this checks where we are in relation to October 15, 1582, the beginning
    # of the Gregorian calendar.
    if ((year < 1582) or
        (year == 1582 and month < 10) or
        (year == 1582 and month == 10 and day < 15)):
        # before start of Gregorian calendar
        B = 0
    else:
        # after start of Gregorian calendar
        A = math.trunc(yearp / 100.)
        B = 2 - A + math.trunc(A / 4.)
        
    if yearp < 0:
        C = math.trunc((365.25 * yearp) - 0.75)
    else:
        C = math.trunc(365.25 * yearp)
        
    D = math.trunc(30.6001 * (monthp + 1))
    
    jd = B + C + D + day + 1720994.5
    
    return jd

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

def fit_blackbody( lam, flux, interactive=True, guessT=10000, plot=True):
    '''Fit a blackbody curve to an input spectrum.
       lam: wavelength (A)
       flux: flux (flam)
    '''
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
        print "Click on points to which you'd like to fit the BB curve."
        print " Left click: choose a point"
        print " Right click: remove most recent point"
        print " Middle click: done"
        points = plt.ginput(n=0, timeout=0)
        xPoints, yPoints = map( np.array, zip(*points) )

        # pick a first-guess scale factor
        guessA = yPoints[0]/black_body_lam( xPoints[0], guessT )

        [T,A], pcov = curve_fit(bb_fit, xPoints, yPoints, p0=[guessT, guessA], maxfev=100000)

        print 'Best-fit T:', round(T, 2)

        plt.plot( x, y,'k',label='Spectrum' )
        plt.plot( x, bb_fit(x, T, A), 'r',label='BB fit %s'%round(T, 4))
        plt.title( 'Blackbody fitting %s'%spectra[i] )
	ylabel('$F_\lambda$ [erg s$^{-1}$ cm$^{-2}$ $\AA$$^{-1}$]')
        xlabel('Wavelength in [$\AA$]')
	legend(loc=0,markerscale=0.5,prop={'size':8},ncol=1)
	savefig('Figures/%s_BB_fitting_%s.png'%(spectra[i],round(T, 2)))
	savefig('Figures/%s_BB_fitting_%s.eps'%(spectra[i],round(T, 2)), format='eps', dpi=1000)

        plt.draw()
        print ' (click plot to dismiss)'
        tmp = plt.ginput(n=1, timeout=120)
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
            plt.plot( x, bb_fit(x, T, A), 'r',label='BB fit %s'%T)
            plt.title( 'Blackbody fitting' )
	    ylabel('$F_\lambda$ [erg s$^{-1}$ cm$^{-2}$ $\AA$$^{-1}$]')
            xlabel('Wavelength in [$\AA$]')
	    legend(loc=0,markerscale=0.5,prop={'size':8},ncol=1)
            plt.show()

    return T, A


######## SN Data ##############################

R_V= 3.1
z_SN=0.017105
A_SN=0.077   # Milky Way extinction #

###### Inventory of all the files #######
dir='/home/thomas/Bureau/BlackBody_fitting/spectres/' # Folder with files
files=os.listdir(dir)  #liste les fichiers
files.sort()  #tri alphabetique
nfiles=size(files)

spectra=[]
for i in range(nfiles):
	if files[i].endswith('.dat'):  #only text files in the folder
		spectra.append(files[i])
nspectra=size(spectra)

###### Inventory of all the files #######
dir='/home/thomas/Bureau/BlackBody_fitting/spectres/' # Folder with files
files=os.listdir(dir)  #liste les fichiers
files.sort()  #tri alphabetique
nfiles=size(files)

spectra=[]
for i in range(nfiles):
	if files[i].endswith('.dat'):  #only text files in the folder
		spectra.append(files[i])
nspectra=size(spectra)


## Here we do a loop for all the spectra


interactive=input("If you want to choose points on the spectrum for BB fitting, 1 for yes ")

for i in range(nspectra):
	f,(ax1)=subplots(1,sharex=True,sharey=False)
	f.subplots_adjust(hspace=0.0)
	os.chdir("/home/thomas/Bureau/BlackBody_fitting/spectres/") #where spectra are

	spectra_SN=np.loadtxt(spectra[i]).transpose()
	lam=spectra_SN[0] #observed wavelength in Angs
	flux=spectra_SN[1] #observed flux in erg/s/cm2/ang
	name_SN=string.split(spectra[i],'_')[0]

	os.chdir("/home/thomas/Bureau/BlackBody_fitting/") #where spectra are

	if interactive==1:

		fit_blackbody( lam, flux, interactive=True, guessT=10000, plot=True)

	else:
	
		guessT=10000
      		guessA = flux[0]/black_body_lam(lam[0], guessT )

        	[T,A], pcov = curve_fit(bb_fit, lam, flux, p0=[guessT, guessA], maxfev=100000)

        	print 'Best-fit T:', round(T, 2)
        	ax1.plot( lam, flux,'k',label='Spectrum' )
        	ax1.plot( lam, bb_fit(lam, T, A), 'r',label='BB fit %s'%round(T,2))
        	title( 'Blackbody fitting %s'%spectra[i] )
		ax1.set_ylabel('$F_\lambda$ [erg s$^{-1}$ cm$^{-2}$ $\AA$$^{-1}$]')
        	ax1.set_xlabel('Wavelength in [$\AA$]')
		legend(loc=0,markerscale=0.5,prop={'size':8},ncol=1)
		savefig('Figures/%s_BB_fitting_%s.png'%(spectra[i],round(T, 2)))
		savefig('Figures/%s_BB_fitting_%s.eps'%(spectra[i],round(T, 2)), format='eps', dpi=1000)
        	show()
		close()
