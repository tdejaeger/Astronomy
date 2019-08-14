##############################################################
#
#	ccm_unred: Deredden a flux vector using the CCM 1989 parameterization
#
#	Cardelli_coeff: Calculate a,b and a+b/Rv for the Cardelli dust 
#	law given a wavelength lam in angstroms
#			
#	calc_Av_from_Balmer_decrement: derive extinction using Balmer decrement (Ha/Hb)
#
#
##############################################################

def ccm_unred(wave, flux, av, **kwargs):

    """
    NAME:
     CCM_UNRED
    PURPOSE:
     Deredden a flux vector using the CCM 1989 parameterization
    EXPLANATION:
     The reddening curve is that of Cardelli, Clayton, and Mathis (1989 ApJ.
     345, 245), including the update for the near-UV given by O'Donnell 
     (1994, ApJ, 422, 158).   Parameterization is valid from the IR to the 
     far-UV (3.5 microns to 0.1 microns).    

     Users might wish to consider using the alternate procedure FM_UNRED
     which uses the extinction curve of Fitzpatrick (1999).
    
    CALLING SEQUENCE:
     ccm_unred(wave, flux, ebv [, R_V = ])      
     
    INPUT:
     WAVE - wavelength vector (Angstroms)
     FLUX - calibrated flux vector, same number of elements as WAVE
             If only 3 parameters are supplied, then this vector will
             updated on output to contain the dereddened flux.
     EBV  - color excess E(B-V), scalar.  If a negative EBV is supplied,
             then fluxes will be reddened rather than deredenned.

    OUTPUT:
     FUNRED - unreddened flux vector, same units and number of elements
             as FLUX

    OPTIONAL INPUT KEYWORD
     R_V - scalar specifying the ratio of total selective extinction
             R(V) = A(V) / E(B - V). If not specified, then R_V = 3.1
             Extreme values of R(V) range from 2.75 to 5.3

    EXAMPLE:
     Determine how a flat spectrum (in wavelength) between 1200 A and 3200 A
     is altered by a reddening of E(B-V) = 0.1.  Assume an "average"
     reddening for the diffuse interstellar medium (R(V) = 3.1)

       >>> w = 1200 + arange(40)*50       #Create a wavelength vector
       >>> f = w*0 + 1                    #Create a "flat" flux vector
       >>> fnew = ccm_unred(w, f, -0.1)   #Redden (negative E(B-V)) flux vector
       >>> plot(w,fnew)                   

    NOTES:
     (1) The CCM curve shows good agreement with the Savage & Mathis (1979)
             ultraviolet curve shortward of 1400 A, but is probably
             preferable between 1200 and 1400 A.
     (2)  Many sightlines with peculiar ultraviolet interstellar extinction 
             can be represented with a CCM curve, if the proper value of 
             R(V) is supplied.
     (3)  Curve is extrapolated between 912 and 1000 A as suggested by
             Longo et al. (1989, ApJ, 339,474)
     (4)  Use the 4 parameter calling sequence if you wish to save the 
             original flux vector.
     (5)  Valencic et al. (2004, ApJ, 616, 912) revise the ultraviolet CCM
             curve (3.3 -- 8.0 um-1).    But since their revised curve does
             not connect smoothly with longer and shorter wavelengths, it is
             not included here.

    REQUIRED MODULES:
     scipy, numpy
    REVISION HISTORY:
       Written   W. Landsman        Hughes/STX   January, 1992
       Extrapolate curve for wavelengths between 900 and 1000 A   Dec. 1993
       Use updated coefficients for near-UV from O'Donnell   Feb 1994
       Allow 3 parameter calling sequence      April 1998
       Converted to IDLV5.0                    April 1998
       Ported to Python        C. Theissen    August 2012
    """

    # Import modules   
    import numpy as n

    # Set defaults   
    R_V = 3.1

    for key in kwargs:
        if key.lower() == 'r_v':
            R_V = kwargs[key]
            
    if isinstance(wave, int) or isinstance(wave, float):
        x = 10000. / n.array([wave])              # Convert to inverse microns
    else:
        x = 10000. / n.array(wave)                # Convert to inverse microns 

    npts = len( x )

    a = n.zeros((npts))  
    b = n.zeros((npts))
    
    ###############################

    good = n.where( (x > 0.3) & (x < 1.1) )     # Infrared
    Ngood = len(x[good])


    if Ngood > 0:
        a[good] =  0.574 * x[good]**(1.61)
        b[good] = -0.527 * x[good]**(1.61)

    ###############################

    good = n.where( (x >= 1.1) & (x < 3.3) )     # Optical/NIR
    Ngood = len(good[0])

    if Ngood > 0:                               # Use new constants from O'Donnell (1994)
        y = x[good] - 1.82
        #c1 = n.array([ 0.32999, -0.77530, 0.01979, 0.72085,        # Original
        #               -0.02427,  -0.50447, 0.17699, 1. ])         # coefficients              
        #c2 = n.array([ -2.09002, 5.30260, -0.62251, -5.38434,       # from CCM89
        #               1.07233, 2.28305, 1.41338, 0. ]) 
        c1 = n.array([ -0.505 , 1.647, -0.827, -1.718,              # New coefficients
                       1.137, 0.701, -0.609, 0.104, 1. ])           # from O'Donnell
        c2 = n.array([ 3.347,  -10.805, 5.491, 11.102,              # (1994)
                       -7.985, -3.989, 2.908, 1.952, 0. ])

        a[good] = n.polyval(c1, y)
        b[good] = n.polyval(c2, y)

    ###############################

    good = n.where( (x >= 3.3) & (x < 8) )                # Mid-UV
    Ngood = len(x[good])

    if Ngood > 0:
        y = x[good]
        F_a = n.zeros((Ngood))
        F_b = n.zeros((Ngood))
        good1 = n.where( (y > 5.9) )
        Ngood1 = len(y[good1])
        
        if Ngood1 > 0:
            y1 = y[good1] - 5.9
            F_a[good1] = -0.04473 * y1**2 - 0.009779 * y1**3
            F_b[good1] = 0.2130 * y1**2  +  0.1207 * y1**3
    
        a[good] = 1.752 - 0.316*y - (0.104 / ( (y-4.67)**2 + 0.341 )) + F_a
        b[good] = -3.090 + 1.825*y + (1.206 / ( (y-4.62)**2 + 0.263 )) + F_b

    ###############################

    good = n.where( (x >= 8) & (x < 11) )              #Far-UV
    Ngood = len(x[good])

    if Ngood > 0:
        y = x[good] - 8.
        c1 = [ -0.070, 0.137, -0.628, -1.073 ]
        c2 = [ 0.374, -0.420, 4.257, 13.670 ]
        
        a[good] = n.polyval(c1, y)
        b[good] = n.polyval(c2, y)

    ###############################

    # Now apply extinction correction to input flux vector

    A_V = av
    A_lambda = A_V * (a + b / R_V)
    return flux * 10.**(0.4 * A_lambda)       

"""
Calculate a,b and a+b/Rv for the Cardelli dust law given a wavelength lam in angstroms
"""
def Cardelli_coeff(lamb,Rv):
	import numpy as np
	scalar=np.isscalar(lamb)
	x=1e4/np.array(lamb,ndmin=1) #CCM x is 1/microns
	a,b=np.ndarray(x.shape,x.dtype),np.ndarray(x.shape,x.dtype)

	if any((x<0.3)|(10<x)):
	    raise ValueError('some wavelengths outside CCM 89 extinction curve range')

	irs=(0.3 <= x) & (x <= 1.1)
	opts = (1.1 <= x) & (x <= 3.3)
	nuv1s = (3.3 <= x) & (x <= 5.9)
	nuv2s = (5.9 <= x) & (x <= 8)
	fuvs = (8 <= x) & (x <= 10)

	#TODO:pre-compute polys

	#CCM Infrared
	a[irs]=.574*x[irs]**1.61
	b[irs]=-0.527*x[irs]**1.61

	#CCM NIR/optical
	a[opts]=np.polyval((.32999,-.7753,.01979,.72085,-.02427,-.50447,.17699,1),x[opts]-1.82)
	b[opts]=np.polyval((-2.09002,5.3026,-.62251,-5.38434,1.07233,2.28305,1.41338,0),x[opts]-1.82)

	#CCM NUV
	a[nuv1s]=1.752-.316*x[nuv1s]-0.104/((x[nuv1s]-4.67)**2+.341)
	b[nuv1s]=-3.09+1.825*x[nuv1s]+1.206/((x[nuv1s]-4.62)**2+.263)

	y=x[nuv2s]-5.9
	Fa=-.04473*y**2-.009779*y**3
	Fb=-.2130*y**2-.1207*y**3
	a[nuv2s]=1.752-.316*x[nuv2s]-0.104/((x[nuv2s]-4.67)**2+.341)+Fa
	b[nuv2s]=-3.09+1.825*x[nuv2s]+1.206/((x[nuv2s]-4.62)**2+.263)+Fb

	#CCM FUV
	a[fuvs]=np.polyval((-.070,.137,-.628,-1.073),x[fuvs]-8)
	b[fuvs]=np.polyval((.374,-.42,4.257,13.67),x[fuvs]-8)

	AloAv = a+b/Rv

	if scalar:
	    return a[0],b[0],AloAv[0]
	else:
	    return a,b,AloAv


def calc_Av_from_Balmer_decrement(halpha,hbeta,halpha_err=None,hbeta_err=None,recom_ratio=2.86,R_v=3.1):
	import numpy as np

	halpha_ext_ratio = cardelli(6562.81, R_v)
	num = recom_ratio * hbeta / halpha
	try:
		A_v = 2.5 * math.log10(num) / (halpha_ext_ratio - hbeta_ext_ratio)    
	except:
		A_v = 999.9
	if (halpha_err != None) and (hbeta_err != None):
		dA_v = abs(2.5 / np.log(10) / (halpha_ext_ratio - hbeta_ext_ratio) * ( (halpha_err/halpha)**2. + (hbeta_err/hbeta)**2.) )**0.5
	else:
		dA_v = None
	return A_v, dA_v
