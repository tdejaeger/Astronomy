#!/usr/bin/env python
from scipy.optimize import minimize
import sys,string,os
import warnings
import numdifftools
import scipy as sc
from scipy import integrate
import numpy as np
import matplotlib.pyplot as plt
import math

c_AA=299792458*1.0e10#in AA/s
h_erg = 6.63e-27 #Plancks constant (erg.s)
c_light=299792.458
dir=os.getcwd()
def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

class SN(object):   
	def import_lc(file):   
		'''
		Import SN data from a datafile in the following format:
		line 1:     name
		line 2:     filter {filter name}
		line 3:     Date   magnitude  error
		...
		line N:     filter {filter name}
		line N+1:   Date   magnitue   error
		....
		'''
		data_sn=open(dir+'/Photometry/%s.txt'%file)
		lines = data_sn.readlines()
		fields = lines[0].split()
		if len(fields) != 1: 
			raise RuntimeError('first line of %s must have 1 field:  name'%file)
		sn_name = fields[0]
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

		return MJD,mags,emags,bands,sn_name


	def get_sn(object_name, **kw):
		'''Attempt to get a sn object from  an existing file name using import_lc()'''

		try:
			s = SN.import_lc(object_name)
		except RuntimeError:
			raise RuntimeError("Not file with this name %s") % object_name

		if np.size(list(set(s[3]) & set( SN.get_filter()[1])))==0:
			raise RuntimeError("Photometric system unknown. Add filters and ZP") 
		return s

	def plot_sn(object_name,show_plot='yes',save_plot='yes'):
		'''Plot and save the light curve in a nice panel figure'''
		MJD,mags,emags,bands,sn_name=SN.get_sn(object_name)
		n_subplot=np.size(bands)
		cols_subplot=2
		# Compute Rows required
		rows_subplot = n_subplot // cols_subplot 
		rows_subplot += n_subplot % cols_subplot
		# Create a Position index
		pos_subplot = range(1,n_subplot + 1)
		# Create main figure
		f= plt.figure(1,figsize=(8, 8), facecolor='w', edgecolor='k')
		f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.0, hspace=0.0)
		for i in range(n_subplot):
			# add every single subplot to the figure with a for loop
			ax = f.add_subplot(rows_subplot,cols_subplot,pos_subplot[i])
			ax.tick_params(axis='both', which='major', labelsize=15,direction='out')
			ax.minorticks_on()
			ax.tick_params(axis='both', which='minor', labelsize=15,direction='out')
			ax.tick_params(which='minor', length=4,width=2, color='k')
			ax.tick_params(which='major',width=2, color='k')
			if pos_subplot[i] % 2 == 0:
				ax.yaxis.tick_right()
			ax.text(0.8, 0.80,'%s'%bands[i],horizontalalignment='center',verticalalignment='center',fontsize=15,transform = ax.transAxes)
			ax.errorbar(MJD['%s'%bands[i]]-MJD['%s'%bands[i]][0],mags['%s'%bands[i]],yerr=emags['%s'%bands[i]],marker='o',color='b',linestyle='None')
			ax.invert_yaxis()
			f.suptitle('%s'%sn_name, fontsize=16)
			f.text(0.5, 0.01, 'Date since first epoch [days]', ha='center',fontsize=15)
			f.text(0.01, 0.5, 'Magnitude [mags]', va='center', rotation='vertical',fontsize=15)
		if show_plot=='yes':
			plt.show()
		else:
			plt.close('all')
		if save_plot=='yes':
			f.savefig(dir+'/Figures/%s.png'%sn_name)

	def date_to_mjd(year,month,day):
		'''Calculate MJD date from date'''
		if month == 1 or month == 2:
			yearp = year - 1
			monthp = month + 12
		else:
			yearp = year
			monthp = month
	    
		# this checks where we are in relation to October 15, 1582, the beginning
		# of the Gregorian calendar.
		if ((year < 1582) or (year == 1582 and month < 10) or (year == 1582 and month == 10 and day < 15)):
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

		return jd- 2400000.5



	def get_spec(file):   
		'''
		Import SN spectrum from a datafile in the following format:
		snXXXXX_yyyymmdd.txt

		line N:    wavelength (AA) flux(erg/s/cm2/A)
		....
		'''
		data_spec=np.loadtxt(dir+'/Spectra/%s.txt'%file).transpose()
		date_spectra=file.split('_')[1]	
		sn_name=file.split('_')[0]
		if len(date_spectra)==8:
			date_spectra=int(date_spectra)
		else:
			raise RuntimeError("File should be snXXXXX_yyyymmdd.txt")
		wl_spec=data_spec[0]
		fl_spec=data_spec[1]

		mjd_spec=SN.date_to_mjd(int(str(date_spectra)[0:4]),int(str(date_spectra)[4:6]),int(str(date_spectra)[6:8]))

		return wl_spec,fl_spec,mjd_spec,sn_name

	def plot_spec(spec_date,show_plot='yes',save_plot='yes'):
		'''Plot and save the light curve in a nice panel figure'''
		wl_spec,fl_spec,ep_spec,sn_name=SN.get_spec(spec_date)
		# Create main figure
		f= plt.figure(1,figsize=(8, 8), facecolor='w', edgecolor='k')
		plt.plot(wl_spec,fl_spec,'k-')

		f.suptitle('%s'%sn_name, fontsize=16)
		f.text(0.5, 0.01, 'Observed wavelength [AA]', ha='center',fontsize=15)
		f.text(0.01, 0.5, 'Flux [erg/s/cm2/A]', va='center', rotation='vertical',fontsize=15)
		if save_plot=='yes':
			f.savefig(dir+'/Figures/%s_spec.png'%sn_name)

		if show_plot=='yes':
			plt.show()
		else:
			plt.close('all')
	def filter_function(filt):
		data_filter=np.loadtxt(dir+'/Filters/%s.txt'%filt).transpose()
		lambda_filt,s_filt=data_filter[0],data_filter[1]
		filt_func=sc.interpolate.interp1d(lambda_filt,s_filt)
		return lambda_filt,s_filt,filt_func


	def get_filter():
		'''
		Download all the transmission function for each available filters
		return a dictionnary with ZP,l_eff,S_lambda(interpolation),lambda_filter		
		
		'''
		ZP_csp=[12.986,14.328,15.111,14.439,14.902,14.545]# ZP http://csp.obs.carnegiescience.edu/data/filters
		l_csp=[3639.0,4350.6,4765.1,5375.2,6223.3,7609.2]
		bands_csp=['u_swope','B_swope','g_swope','V_swope','r_swope','i_swope']
		ZP_kait2=[15.379,14.928,14.992,14.467]
		l_kait2=[4364,5389,6297,8077]
		bands_kait2=['B_kait2','V_kait2','R_kait2','I_kait2']
		ZP_kait3=[15.352,14.935,15.025,14.470]
		l_kait3=[4398,5397,6323,8076]
		bands_kait3=['B_kait3','V_kait3','R_kait3','I_kait3']
		ZP_kait4=[15.266,14.937,14.989,14.446]
		l_kait4=[4445,5389,6273,8061]
		bands_kait4=['B_kait4','V_kait4','R_kait4','I_kait4']
		ZP_nickel1=[15.241,14.843,14.945,14.736]
		l_nickel1=[4369,5329,6259,8125]
		bands_nickel1=['B_nickel1','V_nickel1','R_nickel1','I_nickel1']
		ZP_nickel2=[15.241,14.843,14.945,14.736]
		l_nickel2=[4369,5329,6259,8175]
		bands_nickel2=['B_nickel2','V_nickel2','R_nickel2','I_nickel2']
		ZP_sdss=[13.177,14.663,14.475,14.157,13.363] #ugriz cross checked using 5 standard stars
		l_sdss=[3594.9,4640.4,6122.3,7439.5,9987.1]#ugriz
		bands_sdss=['u_sdss','g_sdss','r_sdss','i_sdss','z_sdss']
		ZP_snls=[15.485,14.817,14.569,14.070] #griz cross checked using BD17 and Vega from Regnault et al. 09
		l_snls=[4870,6300,7760,9250]#griz
		bands_snls=['g_snls','r_snls','i_snls','z_snls']
		ZP_des=[13.5969,-13.7707,-14.0583,-14.501,-15.6223]
		ZP_des=[5.658,5.483,5.196,4.753,3.632]
		l_des=[4842,6439,7760,9172]#griz
		bands_des=['g_des','r_des','i_des','z_des']

		bands_tot=np.concatenate((bands_csp,bands_kait2,bands_kait3,bands_kait4,bands_nickel1,bands_nickel2,bands_sdss,bands_snls,bands_des))
		ZP_tot=np.concatenate((ZP_csp,ZP_kait2,ZP_kait3,ZP_kait4,ZP_nickel1,ZP_nickel2,ZP_sdss,ZP_snls,ZP_des))
		l_eff=np.concatenate((l_csp,l_kait2,l_kait3,l_kait4,l_nickel1,l_nickel2,l_sdss,l_snls,l_des))
		dict_filter=dict()
		for i in range(np.size(bands_tot)):
			dict_filter["%s"%bands_tot[i]] = [ZP_tot[i],l_eff[i],SN.filter_function(bands_tot[i])[2],SN.filter_function(bands_tot[i])[0]]
		return dict_filter,bands_tot


	def mangling(spec_date,object_name):
		wl_spec,fl_spec,mjd_spec=SN.get_spec(spec_date)[0:3]
		MJD,mags,emags,bands,sn_name=SN.get_sn(object_name)
		dict_filter=SN.get_filter()[0]
		#Filter available, which filters are totally inside the spectrum range and for which we have photometry
		filter_available=[]
		for filt in bands:
			#Inside the spectrum range
			if (min(wl_spec)<min(dict_filter[filt][3])) and (max(wl_spec)>max(dict_filter[filt][3])): 
				# with photometry points before and after spectrum epoch
				if (mjd_spec<max(MJD[filt])) and (mjd_spec>min(MJD[filt])):
					filter_available.append(filt)

		if np.size(filter_available)==0:
			raise RuntimeError("Need at least 1 spectrum for the mangling")

		else:
			#1D interpolation of the spectrum
			F_spec_func=sc.interpolate.interp1d(wl_spec,fl_spec)
			f_filter=np.zeros(np.size(filter_available))
			l_eff_X=np.zeros(np.size(filter_available))

			mag_tot_obs=[]
			for i,filt in enumerate(filter_available):
				#observed magnitudes
				X_func=sc.interpolate.interp1d(MJD[filt],mags[filt]) 
				mag_X_obs=X_func(mjd_spec)*1
				mag_tot_obs.append(mag_X_obs)

				#Synthetic magnitude from the spectrum m=-2.5log10[int (F*S*lam)/ch)]+ZP

				m_X_syn=-2.5*np.log10(1.0/(h_erg*c_AA)*integrate.simps(F_spec_func(dict_filter[filt][3])*dict_filter[filt][2](dict_filter[filt][3])*dict_filter[filt][3],dict_filter[filt][3]))+dict_filter[filt][0]
				#Diffrence between the two mags
				delta_X=mag_X_obs-m_X_syn
				f_X=10**(-0.4*delta_X)
				#points of the warping function: coefficient vs lam eff
				f_filter[i]=f_X
				l_eff_X[i]=dict_filter[filt][1]

			#sort by wavelength
			x_mang,y_mang,filter_available=zip(*sorted(zip(l_eff_X,f_filter,filter_available)))
			#If you use different photometric system with similar filter, you need to remove some points
			# 500 Ang min
			ind_good=[ i for i in np.arange(0,np.size(x_mang)-1,1) if x_mang[i+1]-x_mang[i]>500]
			ind_good.append(np.size(x_mang)-1)
			mag_tot_obs=np.array(mag_tot_obs)[ind_good] 
			x_mang,y_mang=np.array(x_mang)[ind_good],np.array(y_mang)[ind_good]
			x_mang,y_mang=list(x_mang),list(y_mang)

			#First and last SED wl will have the same values of the first and the last filter respectively
			# the SED wavelength range is much larger than the filter range, that is why we do that
			y_mang.insert(0,y_mang[0]) 
			x_mang.insert(0,min(wl_spec))
			y_mang.insert(np.size(y_mang),y_mang[np.size(y_mang)-1])
			x_mang.insert(np.size(y_mang),max(wl_spec))

			#Spline the warping function
			if np.size(x_mang)==3:
				tck = sc.interpolate.splrep(x_mang, y_mang, k=2,s=0)
			elif np.size(x_mang)>3:
				tck = sc.interpolate.splrep(x_mang, y_mang, k=3,s=0)
			elif np.size(x_mang)==2:
				tck = sc.interpolate.splrep(x_mang, y_mang, k=1,s=0)
			coeff_flux_mang= sc.interpolate.splev(wl_spec, tck, der=0) # all the coeff of each wavelength
			#Mangling the template SED
			fl_spec_mang=(coeff_flux_mang)*(fl_spec)
			fl_spec_mang_func=sc.interpolate.interp1d(wl_spec,fl_spec_mang)

			# Mangling check. We use the warping function and recalculate the synthetic magnitudes
			mag_tot_mang=[]
			for i,filt in enumerate(filter_available):
				#observed magnitudes
				X_func=sc.interpolate.interp1d(MJD[filt],mags[filt]) 
				mag_X_obs=X_func(mjd_spec)
				#Synthetic magnitude from the spectrum m=-2.5log10[int (F*S*lam)/ch)]+ZP

				m_X_syn=-2.5*np.log10(1.0/(h_erg*c_AA)*integrate.simps(fl_spec_mang_func(dict_filter[filt][3])*dict_filter[filt][2](dict_filter[filt][3])*dict_filter[filt][3],dict_filter[filt][3]))+dict_filter[filt][0]

				mag_tot_mang.append(m_X_syn)
			mag_tot_mang=np.array(mag_tot_mang)[ind_good] 
			maxiter=11
			while False in (abs(np.array(mag_tot_obs)-np.array(mag_tot_mang))<0.07):
				maxiter = maxiter - 1
				if maxiter <= 0:
					print ("Did not converge.")
					break
				mag_tot_mang=[]
				for i,filt in enumerate(filter_available):
					#observed magnitudes
					X_func=sc.interpolate.interp1d(MJD[filt],mags[filt]) 
					mag_X_obs=X_func(mjd_spec)
					#Synthetic magnitude from the spectrum m=-2.5log10[int (F*S*lam)/ch)]+ZP
					m_X_syn=-2.5*np.log10(1.0/(h_erg*c_AA)*integrate.simps(F_spec_func(dict_filter[filt][3])*dict_filter[filt][2](dict_filter[filt][3])*dict_filter[filt][3],dict_filter[filt][3]))+dict_filter[filt][0]
					mag_tot_mang.append(m_X_syn)
					#Diffrence between the two mags
					delta_X=mag_X_obs-m_X_syn
					f_X=10**(-0.4*delta_X)
					#points of the warping function: coefficient vs lam eff
					f_filter[i]=f_X
					l_eff_X[i]=dict_filter[filt][1]
				mag_tot_mang=np.array(mag_tot_mang)[ind_good]
				#sort by wavelength
				x_mang,y_mang,filter_available=zip(*sorted(zip(l_eff_X,f_filter,filter_available)))
				#If you use different photometric system with similar filter, you need to remove some points
				# 500 Ang min
				ind_good=[ i for i in np.arange(0,np.size(x_mang)-1,1) if x_mang[i+1]-x_mang[i]>500]
				ind_good.append(np.size(x_mang)-1)
				x_mang,y_mang=np.array(x_mang)[ind_good],np.array(y_mang)[ind_good]
				x_mang,y_mang=list(x_mang),list(y_mang)

				#Spline the warping function
				if np.size(x_mang)==3:
					tck = sc.interpolate.splrep(x_mang, y_mang, k=2,s=0)
				elif np.size(x_mang)>3:
					tck = sc.interpolate.splrep(x_mang, y_mang, k=3,s=0)
				elif np.size(x_mang)==2:
					tck = sc.interpolate.splrep(x_mang, y_mang, k=1,s=0)
				coeff_flux_mang= sc.interpolate.splev(wl_spec, tck, der=0) # all the coeff of each wavelength
				#Mangling the template SED
				fl_spec_mang=(coeff_flux_mang)*(fl_spec)
				fl_spec_mang_func=sc.interpolate.interp1d(wl_spec,fl_spec_mang)	
		print('')	
		print('### Filters used ###')	
		print(np.array(filter_available)[ind_good])	
		print('### Observed magnitude ###')	
		print(np.around(mag_tot_obs,2))
		print('### Synthetic magnitude after mangling ###')
		print(np.around(mag_tot_mang,2))			
		return fl_spec_mang

	def plot_spec_warp(object_name,spec_date,show_plot='yes',save_plot='yes'):
		spec=object_name+'_'+spec_date
		fl_spec_mang=SN.mangling(spec,object_name)
		wl_spec,fl_spec,ep_spec,sn_name=SN.get_spec(spec)
		# Create main figure
		f= plt.figure(1,figsize=(8, 8), facecolor='w', edgecolor='k')
		plt.plot(wl_spec,fl_spec,'k-',alpha=0.8,label='Original')
		plt.plot(wl_spec,fl_spec_mang,'r-',label='Warped')

		f.suptitle('%s'%sn_name, fontsize=16)
		f.text(0.5, 0.01, 'Observed wavelength [AA]', ha='center',fontsize=15)
		f.text(0.01, 0.5, 'Flux [erg/s/cm2/A]', va='center', rotation='vertical',fontsize=15)
		plt.legend(loc=0,title='',markerscale=0.8,ncol=6,prop={'size':10})

		if save_plot=='yes':
			f.savefig(dir+'/Figures/%s_spec.png'%sn_name)
		if show_plot=='yes':
			plt.show()
		else:
			plt.close('all')
