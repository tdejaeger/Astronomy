################################################################################################
################################################################################################
#
#		Fit Supernova light curve using Gaussian process 
#		Thomas de Jaeger
#		2017
#
#	    
################################################################################################

import re # use regular patters
import sys # system commands
import string as string # string functions4
import numpy as np # numerical tools
from scipy import *
import os
import math as maths
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, show, rcParams
import shutil

########################################## Constantes #######################
rcParams['legend.numpoints']=1
c_light=299792.458#in km/s

### Import Data ###
'''
Import SN data from a datafile in the following format:
line 1:     name z ra decl
line 2:     filter {filter name}
line 3:     Date   magnitude  error
...
line N:     filter {filter name}
line N+1:   Date   magnitue   error
....
'''
file='Data/SN2005J.txt'

data_sn=open(file)
lines = data_sn.readlines()
fields = lines[0].split()
if len(fields) != 5: 
	print('first line of %s must have 5 fields:  name, redshift, AvG,Texp (MJD), Survey'%file)

else:

	name_sn = fields[0]
	z_sn=float(fields[1])
	AvG=float(fields[2])
	Texp=float(fields[3])
	survey=fields[4]

lines = lines[1:]
this_filter = None
bands=[]
MJD = {}
mags = {}
emags = {}

for line in lines:
	if line[0] == "#":  continue
	if line.find('filter') >= 0:
		this_filter = line.split()[1]
		bands.append(this_filter)
		MJD[this_filter] = []
		mags[this_filter] = []
		emags[this_filter] = []
	elif this_filter is not None:
		try:
			t,m,em = map(float, str.split(str.strip(line)))
		except:
			print('Bad format in line:\n %s' % (line))
		MJD[this_filter].append(t)
		mags[this_filter].append(m)
		emags[this_filter].append(em)

for f in MJD:
	MJD[f] = np.array(MJD[f])
	mags[f] = np.array(mags[f])
	emags[f] = np.array(emags[f])

#Remove Texp

MJD.update((x, y-Texp) for x, y in MJD.items())

#Interactive mode
interactive = 1

#import Gaussian process library http://dfm.io/george/current/
import george

## Predicted epochs
nepochs =  200+10
epochs = np.linspace(-9,200,nepochs)

#Kernel george: the main ingredient of a gaussian process is the kernel.
#It describes how correlated each point is with every other.
kernel = george.kernels.ExpSquaredKernel(0.5)
model = george.GP(kernel)

# Initialisation new mag and error of the GP
mu, std = np.zeros((len(bands),nepochs),dtype='f'),np.zeros((len(bands),nepochs),dtype='f')
sfig, ax = plt.subplots(int(round(size(bands)/2,0)+size(bands)%2),2,sharex=True,sharey=True)
fig, ax1 = plt.subplots(facecolor='w', edgecolor='k')  # 
#ref epoch
refepoch, errfac = np.zeros(size(bands)), np.zeros(size(bands))
for i,j in enumerate(bands):
	#transform magnitude to flux
	ep,fl,efl=MJD[j],10**(-0.4*(mags[j])),10**(-0.4*(mags[j]))*(emags[j]*log(10)*1.0/(2.5)) 

	fmax=max(fl)
	ep,fl,efl = ep,fl/fmax,efl/fmax
	## --interactive GP and plot
	repeat = 1
	trefepoch, terrfac = -5, 1.5 #refepoch taken before the explosion. terrfac, if noisy increase value   
	plmaxep = 140
	while repeat == 1:
		## --george
		kernelGP = np.var(fl) * kernel
		modelGP = george.GP(kernelGP)

		lgep = np.arcsinh(ep-trefepoch)
		lgepochs = np.arcsinh(epochs-trefepoch)#arcsinh

		## -- Pre-compute the factorization of the matrix.
		modelGP.compute(lgep, efl*terrfac)

		## -- Compute the log likelihood.
		modelGP.lnlikelihood(fl)

		## -- optimize
		modelGP.optimize(lgep,fl,efl*terrfac)

		## -- Prediction
		tmu, tvar = modelGP.predict(fl, lgepochs)
		tstd = np.sqrt(np.diag(tvar))
		mu[i,:],std[i,:] = np.reshape(tmu,nepochs),np.reshape(tstd,nepochs)
		## --interactive plot
		repeat = 0
		if interactive:
			print("  Current refepoch: %5.2f and errfac: %5.2f " %(trefepoch,terrfac))
			ifig,iax = plt.subplots(1,figsize=(8,5))
			iax.errorbar(ep,fl,efl,label='%s'%bands[i],marker="o",linestyle='None')
			iax.plot(epochs,mu[i,:])
			iax.fill_between(epochs,mu[i,:]-std[i,:],mu[i,:]+std[i,:],alpha=0.5)
			iax.set_xlim(-9,plmaxep)
			iax.set_ylim(-0.25,1.5)
			iax.legend(fontsize=8)
			iax.set_title('%s'%str.split(name_sn,'.')[0])
			ifig.show()
			readvar = input("Enter new set of refepoch,errfac separated by comma "+\
				"(or 'n' for none): ")
			if readvar != 'n':
				trefepoch, terrfac = np.array(readvar.split(','),dtype=float)
				repeat = 1
			plt.close(ifig)
		del modelGP
		
	## -- final values
	refepoch[i], errfac[i] = trefepoch, terrfac

	## -- plot all bands together in flux
	tdataplot = plt.errorbar(ep,fl,efl,label=j,marker="o",linestyle='None')
	ax1.plot(epochs,mu[i,:],color=tdataplot[0].get_color())
	ax1.fill_between(epochs,mu[i,:]-std[i,:],mu[i,:]+std[i,:],alpha=0.5,color=tdataplot[0].get_color())

	## -- plot filter one by one
	p1,p2 = int(i/2),int(i%2)
	ax[p1,p2].errorbar(ep,fl,efl,label=bands[i],marker="o",linestyle='None')
	ax[p1,p2].plot(epochs,mu[i,:])
	ax[p1,p2].fill_between(epochs,mu[i,:]-std[i,:],mu[i,:]+std[i,:],alpha=0.5)
	ax[p1,p2].legend(fontsize=8)
	ax[p1,p2].set_xlim(-20,plmaxep)
	ax[p1,p2].set_ylim(-0.25,1.5)

	mu[i,:], std[i,:] = mu[i,:]*fmax, std[i,:]*fmax

plt.close('all')
## -- plot all bands together but magnitude
fig, ax1 = plt.subplots(facecolor='w', edgecolor='k')  # ax1 is apparent magnitude

for i,j in enumerate(bands):
	mag_GP=np.array(-2.5*log10(mu[i]))
	std_GP=std[i]*2.5/(mu[i]*log(10))

	tdataplot = ax1.errorbar(MJD[j], mags[j],yerr=emags[j],label=j,marker="o",ms=5,linestyle='None')
	ax1.plot(epochs,mag_GP,color=tdataplot[0].get_color())
	ax1.fill_between(epochs,mag_GP-std_GP,mag_GP+std_GP,alpha=0.5,color=tdataplot[0].get_color())
ax1.set_ylim([min(mags['i'])-1,max(mags['u'])+1])
ax1.set_xlim([min(np.array(ep))-7,120])
ax1.tick_params(axis='both', which='major', labelsize=20)
ax1.set_xlabel('Epoch since explosion [days]',fontsize=15,fontweight='bold')
ax1.set_ylabel('Magnitude [mag]',fontsize=15,fontweight='bold')
ax1.legend(loc=0,title='',ncol=2,prop={'size':12})
plt.minorticks_on()
plt.gca().invert_yaxis()
plt.show()
#fig.savefig('%s_mag.png'%str.split(name_sn,'.')[0])
plt.close("all")

save_file=input('Do you want to same the GP in a file? y or n: ')
if save_file=='y':
	## --print information into one file per SN
	for i,j in enumerate(bands):
		matrice=np.zeros((size(epochs),3))
		matrice[:,0]=np.array(epochs) #epo
		matrice[:,1]=np.array(-2.5*log10(mu[i,:])) #mag
		matrice[:,2]=np.array(std[i,:]*2.5/(mu[i,:]*log(10))) #err mag
		np.savetxt('%s_%s.txt'%(name_sn,j),matrice)

