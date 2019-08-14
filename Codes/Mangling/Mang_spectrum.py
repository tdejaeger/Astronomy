import re # use regular patters
import sys # system commands
import os 
import string as string # string functions4
import math
import numpy as num # numerical tools
from scipy import *
from pylab import *
import glob 
import scipy.interpolate
from scipy.optimize import leastsq
from scipy.optimize import minimize
import numdifftools
from scipy import integrate
c_AA=299792458*1.0e10#in AA/s
h_erg = 6.63e-27 #Plancks constant (erg.s)
c_light=299792.458

sn_name=str('sn2016hnk')

## Import filter transmissions ###
exec(open("Programmes/filter_trans.py").read())
'''->filt_B_func,lambda_B,s_B

dict_ZP['filter']=ZP, lamda_eff,f_transmi,lambda_filter
'''

exec(open("Programmes/get_sn.py").read())
'''-> MJD[band],mags[band],emags[band],bands'''

## Plot SN ###
exec(open("Programmes/plot_sn.py").read())

## Grouped data by observation day ###
exec(open("Programmes/group_mag.py").read())
'''->ret_data['MJD'],ret_data['B_swope'],ret_data['e_B_swope']'''

## load spectra ###
exec(open("Programmes/spectra.py").read())
'''->epoch_spectra,wl_spectra,flux_spectra'''

#Mangling
exec(open("Programmes/mangling.py").read())


