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
from kapteyn import kmpfit
import shutil
import george

################################################################################################
################################################################################################
################################################################################################
#
#          Gaussian process fit 
#          Thomas de Jaeger
#          2017
#
################################################################################################
################################################################################################
################################################################################################
def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

########################################## Constantes #######################
rcParams['legend.numpoints']=1
c_light=299792.458#in km/s
##
data_sn=np.loadtxt('Data_all_SNe_CSP.dat',usecols=[1,2,3,4,5,6,7,8,9,10,11]).transpose()
A_V=data_sn[0]
redshift=data_sn[1]*1.0/c_light
err_z_hel=data_sn[2]*1.0/c_light
JD_explosion=data_sn[3]
err_JD_explosion=data_sn[4]
optd=data_sn[5]
s2=data_sn[6]
err_s2=data_sn[7]
l_plateau=data_sn[8]
z_CMB=data_sn[9]*1.0/c_light
err_z_CMB=data_sn[10]*1.0/c_light
JD_ref=2453000

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

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

####################################################################################
#
#               On va charger les magnitudes
###################################################################################
###### On fait l'inventaire de tous les fichiers #######
dir='/home/thomas/Recherche/BERKELEY/CSP_photo/Raw_opt'
files=os.listdir(dir)  #liste les fichiers
files.sort()  #tri alphabetique
nfiles=size(files)
table=[]
name_SN=[]
for i in range(nfiles):
	if files[i].endswith('.dat'):  #on selectionne
		table.append(files[i])
		name_SN.append(str.split(files[i],'.')[0])
ntable=size(table)

SN_epoch=[]
u=[]
err_u=[]
g=[]
err_g=[]
r=[]
err_r=[]
iprime=[]
err_i=[]
B=[]
err_B=[]
V=[]
err_V=[]
for i in range(ntable):

	SN_photo=np.loadtxt('/home/thomas/Recherche/BERKELEY/CSP_photo/Raw_opt/%s'%table[i]).transpose()
	SN_epoch.append(SN_photo[0]-JD_explosion[i])
	u.append(SN_photo[1])
	err_u.append(SN_photo[2])
	g.append(SN_photo[3])
	err_g.append(SN_photo[4])
	r.append(SN_photo[5])
	err_r.append(SN_photo[6])
	iprime.append(SN_photo[7])
	err_i.append(SN_photo[8])
	B.append(SN_photo[9])
	err_B.append(SN_photo[10])
	V.append(SN_photo[11])
	err_V.append(SN_photo[12])

####################### Magnitude sans les 999.00 #################################

phase_u_final=[]
phase_g_final=[]
phase_r_final=[]
phase_i_final=[]
phase_B_final=[]
phase_V_final=[]

u_final=[]
g_final=[]
r_final=[]
i_final=[]
B_final=[]
V_final=[]

err_u_final=[]
err_g_final=[]
err_r_final=[]
err_i_final=[]
err_B_final=[]
err_V_final=[]

for i in range(ntable):


	phase_u_final.append(np.array(SN_epoch[i][np.array(u[i])!=0.0]))
	phase_g_final.append(np.array(SN_epoch[i][np.array(g[i])!=0.0]))
	phase_r_final.append(np.array(SN_epoch[i][np.array(r[i])!=0.0]))
	phase_i_final.append(np.array(SN_epoch[i][np.array(iprime[i])!=0.0]))
	phase_B_final.append(np.array(SN_epoch[i][np.array(B[i])!=0.0]))
	phase_V_final.append(np.array(SN_epoch[i][np.array(V[i])!=0.0]))

	u_final.append(np.array(u[i][np.array(u[i])!=0.0]))
	g_final.append(np.array(g[i][np.array(g[i])!=0.0]))
	r_final.append(np.array(r[i][np.array(r[i])!=0.0]))
	i_final.append(np.array(iprime[i][np.array(iprime[i])!=0.0]))
	B_final.append(np.array(B[i][np.array(B[i])!=0.0]))
	V_final.append(np.array(V[i][np.array(V[i])!=0.0]))

	err_u_final.append(np.array(err_u[i][np.array(u[i])!=0.0]))
	err_g_final.append(np.array(err_g[i][np.array(g[i])!=0.0]))
	err_r_final.append(np.array(err_r[i][np.array(r[i])!=0.0]))
	err_i_final.append(np.array(err_i[i][np.array(iprime[i])!=0.0]))
	err_V_final.append(np.array(err_V[i][np.array(V[i])!=0.0]))
	err_B_final.append(np.array(err_B[i][np.array(B[i])!=0.0]))

'''
On a ttes les magnitudes observees, phase_X_final,X_final,err_X_final 
'''


'''
All the raw photometry 
'''

# name filters
filters=['B','g','V','r','i']
#Interactive mode
interactive = 1
## Predicted epochs
nepochs =  200+10
epochs = np.linspace(-9,200,nepochs)

#Kernel george
kernel = george.kernels.ExpSquaredKernel(0.5)
model = george.GP(kernel)

for i in range(ntable):
	print(i,table[i])


#SN=int(input('Select SN number: '))
ZP_tot=[14.328,15.111,14.439,14.902,14.545]# ZP http://csp.obs.carnegiescience.edu/data/filters
for SN in range(ntable):
	SN=SN+41
	mag_TOT=[np.array(B_final[SN]),np.array(g_final[SN]),np.array(V_final[SN]),np.array(r_final[SN]),np.array(i_final[SN])]

	err_mag_TOT=[np.array(err_B_final[SN]),np.array(err_g_final[SN]),np.array(err_V_final[SN]),np.array(err_r_final[SN]),np.array(err_i_final[SN])]

	epoch_TOT=[np.array(phase_B_final[SN]),np.array(phase_g_final[SN]),np.array(phase_V_final[SN]),np.array(phase_r_final[SN]),np.array(phase_i_final[SN])]

	# Initialisation new mag and error
	mu, std = np.zeros((len(mag_TOT),nepochs),dtype='f'),np.zeros((len(mag_TOT),nepochs),dtype='f')
	sfig, ax = plt.subplots(int(round(size(filters)/2,0)+size(filters)%2),2,sharex=True,sharey=True)
	fig, ax1 = plt.subplots(facecolor='w', edgecolor='k')  # ax1 is apparent magnitude
	#ref epoch
	refepoch, errfac = np.zeros(size(filters)), np.zeros(size(filters))

	for i in range(size(filters)):
		ep,fl,efl=epoch_TOT[i],10**(-0.4*(np.array(mag_TOT[i])-ZP_tot[i])),10**(-0.4*(np.array(mag_TOT[i])-ZP_tot[i]))*(np.array(err_mag_TOT[i])*log(10)*1.0/(2.5)) 
	     
		fmax = np.max(10**(-0.4*(np.array(mag_TOT[i])-ZP_tot[i])))
		ep,fl,efl = np.array(ep),np.array(fl)/fmax,np.array(efl)/fmax
		## --interactive GP and plot
		repeat = 1
		trefepoch, terrfac = -5, 1.5 #trefepoch epoque ou il y aura le plus de variations. terrfac, si grand grande erreur   
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

			## -- inference

			## -- Prediction
			tmu, tvar = modelGP.predict(fl, lgepochs)
			tstd = np.sqrt(np.diag(tvar))
			mu[i,:],std[i,:] = np.reshape(tmu,nepochs),np.reshape(tstd,nepochs)
			## --interactive plot
			repeat = 0
			if interactive:
				print("  Current refepoch: %5.2f and errfac: %5.2f " %(trefepoch,terrfac))
				ifig,iax = plt.subplots(1,figsize=(8,5))
				iax.errorbar(ep,fl,efl,label='%s'%filters[i],marker="o",linestyle='None')
				iax.plot(epochs,mu[i,:])
				iax.fill_between(epochs,mu[i,:]-std[i,:],mu[i,:]+std[i,:],alpha=0.5)
				iax.set_xlim(-9,plmaxep)
				iax.set_ylim(-0.25,1.5)
				iax.legend(fontsize=8)
				iax.set_title('%s'%str.split(table[SN],'.')[0])
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

		## -- plot all filters together
		tdataplot = plt.errorbar(ep*1.0/(1+redshift[SN]),fl,efl,label=filters[i],marker="o",linestyle='None')
		ax1.plot(epochs*1.0/(1+redshift[SN]),mu[i,:],color=tdataplot[0].get_color())
		ax1.fill_between(epochs*1.0/(1+redshift[SN]),mu[i,:]-std[i,:],mu[i,:]+std[i,:],alpha=0.5,color=tdataplot[0].get_color())

		## -- plot filter one by one
		p1,p2 = int(i/2),int(i%2)
		ax[p1,p2].errorbar(ep,fl,efl,label=filters[i],marker="o",linestyle='None')
		ax[p1,p2].plot(epochs,mu[i,:])
		ax[p1,p2].fill_between(epochs,mu[i,:]-std[i,:],mu[i,:]+std[i,:],alpha=0.5)
		ax[p1,p2].legend(fontsize=8)
		ax[p1,p2].set_xlim(-20,plmaxep)
		ax[p1,p2].set_ylim(-0.25,1.5)

		mu[i,:], std[i,:] = mu[i,:]*fmax, std[i,:]*fmax

	## --plot filters one by one
	sfig.subplots_adjust(hspace=0,wspace=0)
	sfig.text(0.5, 0.04, 'Epoch(d)', ha='center',fontsize=12)
	sfig.text(0.04, 0.5, 'Flux (arbitrary units)', va='center', rotation='vertical',fontsize=12)
	sfig.text(0.5,0.90,'%s'%str.split(table[SN],'.')[0],fontsize=12)
	sfig.savefig('Figures/%s-sep.png'%str.split(table[SN],'.')[0])
	sfig.gca().cla()

	## --plot all filters together
	ax1.set_xlabel('Epoch since explosion rest-frame [days]',fontsize=15,fontweight='bold')
	ax1.set_ylabel('Flux (arbitrary units)',fontsize=15,fontweight='bold')
	plt.xlim(-20,plmaxep)
	plt.ylim(-0.25,1.5)
	plt.legend(fontsize=6)
	plt.title('%s'%str.split(table[SN],'.')[0])
	fig.savefig('Figures/%s_flux.png'%str.split(table[SN],'.')[0])
	fig.gca().cla()
	plt.close("all")
	## -- plot all filters together but magnitude
	fig, ax1 = plt.subplots(facecolor='w', edgecolor='k')  # ax1 is apparent magnitude

	for i in range(size(filters)):
		mag_GP=np.array(-2.5*log10(mu[i]))+median(ZP_tot[i])
		std_GP=std[i]*2.5/(mu[i]*log(10))

		tdataplot = ax1.errorbar(np.array(epoch_TOT[i])*1.0/(1+redshift[SN]), np.array(mag_TOT[i]),yerr=err_mag_TOT[i],label=filters[i],marker="o",ms=5,linestyle='None')
		ax1.plot(np.array(epochs)*1.0/(1+redshift[SN]),mag_GP,color=tdataplot[0].get_color())
		ax1.fill_between(epochs*1.0/(1+redshift[SN]),mag_GP-std_GP,mag_GP+std_GP,alpha=0.5,color=tdataplot[0].get_color())
	ax1.set_ylim([min(mag_TOT[4])-1,max(mag_TOT[0])+1])
	ax1.set_xlim([min(np.array(ep)*1.0/(1+redshift[SN]))-5,120])
	ax1.tick_params(axis='both', which='major', labelsize=20)
	ax1.set_xlabel('Epoch since explosion rest-frame [days]',fontsize=15,fontweight='bold')
	ax1.set_ylabel('Magnitude [mag]',fontsize=15,fontweight='bold')
	ax1.legend(loc=0,title='',ncol=1,prop={'size':16})
	plt.minorticks_on()
	plt.gca().invert_yaxis()
	fig.savefig('Figures/%s_mag.png'%str.split(table[SN],'.')[0])
	close("all")


	## --print information into one file per SN
	for i in range(0,size(filters)):
		matrice=np.zeros((size(epochs),5))
		matrice[:,0]=np.array(epochs)
		matrice[:,1]=np.array(mu[i,:]) #flux
		matrice[:,2]=np.array(std[i,:]) #err flux
		matrice[:,3]=np.array(-2.5*log10(mu[i,:]))+median(ZP_tot[i]) #mag
		matrice[:,4]=np.array(std[i,:]*2.5/(mu[i,:]*log(10))) #err mag
		np.savetxt('GP/%s_%s.txt'%(str.split(table[SN],'.')[0],filters[i]),matrice)

