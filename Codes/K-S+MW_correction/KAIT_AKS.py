import math
import numpy as np # numerical tools
import os
from scipy import integrate
from scipy import interpolate
import scipy as sc
import math as maths
from matplotlib.pyplot import figure, show, rcParams
import matplotlib.pyplot as plt
from function import *
import pandas as pd
################################################################################################
################################################################################################
################################################################################################
#
# 
#			K-S correction
#
#
################################################################################################
################################################################################################
################################################################################################

########################################## Constantes #######################
rcParams['legend.numpoints']=1
c_light=299792.458#in km/s
c_AA=299792458*1.0e10#in AA/s
kB=1.38066e-23; # Boltzmanns constant in J/K
h_erg = 6.63e-27 #Plancks constant (erg.s)
################		SN 		#############################		
SN_name=input('Give the SN name: sn2006ee, sn2007od, sn2008aw, sn2009N: ')
#get SN info:
A_V,z_hel,err_z_hel,JD_explo,err_JD_explo,z_cmb,err_z_cmb=get_SN_info(SN_name)
do_AKS_corr=input('Do you want to do AKS correction for %s, y or n: '%(SN_name))
if do_AKS_corr==str('y'):

	#get SN photometry
	MJD,mags,emags,tel=get_sn(SN_name)

	#Transform MJD to epoch since the explosion
	MJD.update((x, y-JD_explo) for x, y in MJD.items())

	# On va interpoler les magnitudes!!
	B_band,B_band_plus,V_band,V_band_plus,R_band,R_band_plus,I_band,I_band_plus=inter_mag(MJD,mags,emags)

	#Filter CSP
	exec(open("filters_CSP.py").read())
	#Filter KAIT
	exec(open("filters_KAIT.py").read())

	#Define dict for Zero Point
	ZP={}
	ZP['csp']=[]
	ZP['csp'].append([14.328,14.439,14.902,14.545])#  #BVri
	ZP['kait2']=[]
	ZP['kait2'].append([15.379,14.928,14.992,14.467])	#BVRI
	ZP['kait3']=[]
	ZP['kait3'].append([15.352,14.935,15.025,14.470])
	ZP['kait4']=[]
	ZP['kait4'].append([15.266,14.937,14.989,14.446])
	ZP['nickel1']=[]
	ZP['nickel1'].append([15.241,14.843,14.945,14.736])
	ZP['nickel2']=[]
	ZP['nickel2'].append([14.873,14.686,14.806,15.120])

	#Get theoritcal SED from Dessart et al. 2013
	model=pd.read_table('Models/m15_du_sch_mlt3_FeC_mix0p4_list',delim_whitespace=True,header=None)
	Model=model[0].values
	epoch_mod=model[1].values

	#Let's sampling the filter transmission
	sampling_model=10000
	lambda_B_kait2=np.linspace(min(lambda_B_kait2),max(lambda_B_kait2),sampling_model)
	lambda_V_kait2=np.linspace(min(lambda_V_kait2),max(lambda_V_kait2),sampling_model) 
	lambda_R_kait2=np.linspace(min(lambda_R_kait2),max(lambda_R_kait2),sampling_model) 
	lambda_I_kait2=np.linspace(min(lambda_I_kait2),max(lambda_I_kait2),sampling_model) 
	dlambda_B_kait2=lambda_B_kait2[1]-lambda_B_kait2[0]
	dlambda_V_kait2=lambda_V_kait2[1]-lambda_V_kait2[0]
	dlambda_R_kait2=lambda_R_kait2[1]-lambda_R_kait2[0] 
	dlambda_I_kait2=lambda_I_kait2[1]-lambda_I_kait2[0]
	lambda_B_kait3=np.linspace(min(lambda_B_kait3),max(lambda_B_kait3),sampling_model)
	lambda_V_kait3=np.linspace(min(lambda_V_kait3),max(lambda_V_kait3),sampling_model) 
	lambda_R_kait3=np.linspace(min(lambda_R_kait3),max(lambda_R_kait3),sampling_model) 
	lambda_I_kait3=np.linspace(min(lambda_I_kait3),max(lambda_I_kait3),sampling_model) 
	dlambda_B_kait3=lambda_B_kait3[1]-lambda_B_kait3[0]
	dlambda_V_kait3=lambda_V_kait3[1]-lambda_V_kait3[0]
	dlambda_R_kait3=lambda_R_kait3[1]-lambda_R_kait3[0] 
	dlambda_I_kait3=lambda_I_kait3[1]-lambda_I_kait3[0]
	lambda_B_kait4=np.linspace(min(lambda_B_kait4),max(lambda_B_kait4),sampling_model)
	lambda_V_kait4=np.linspace(min(lambda_V_kait4),max(lambda_V_kait4),sampling_model) 
	lambda_R_kait4=np.linspace(min(lambda_R_kait4),max(lambda_R_kait4),sampling_model) 
	lambda_I_kait4=np.linspace(min(lambda_I_kait4),max(lambda_I_kait4),sampling_model) 
	dlambda_B_kait4=lambda_B_kait4[1]-lambda_B_kait4[0]
	dlambda_V_kait4=lambda_V_kait4[1]-lambda_V_kait4[0]
	dlambda_R_kait4=lambda_R_kait4[1]-lambda_R_kait4[0] 
	dlambda_I_kait4=lambda_I_kait4[1]-lambda_I_kait4[0]
	lambda_B_nickel1=np.linspace(min(lambda_B_nickel1),max(lambda_B_nickel1),sampling_model)
	lambda_V_nickel1=np.linspace(min(lambda_V_nickel1),max(lambda_V_nickel1),sampling_model) 
	lambda_R_nickel1=np.linspace(min(lambda_R_nickel1),max(lambda_R_nickel1),sampling_model) 
	lambda_I_nickel1=np.linspace(min(lambda_I_nickel1),max(lambda_I_nickel1),sampling_model) 
	dlambda_B_nickel1=lambda_B_nickel1[1]-lambda_B_nickel1[0]
	dlambda_V_nickel1=lambda_V_nickel1[1]-lambda_V_nickel1[0]
	dlambda_R_nickel1=lambda_R_nickel1[1]-lambda_R_nickel1[0] 
	dlambda_I_nickel1=lambda_I_nickel1[1]-lambda_I_nickel1[0]
	lambda_B_nickel2=np.linspace(min(lambda_B_nickel2),max(lambda_B_nickel2),sampling_model)
	lambda_V_nickel2=np.linspace(min(lambda_V_nickel2),max(lambda_V_nickel2),sampling_model) 
	lambda_R_nickel2=np.linspace(min(lambda_R_nickel2),max(lambda_R_nickel2),sampling_model) 
	lambda_I_nickel2=np.linspace(min(lambda_I_nickel2),max(lambda_I_nickel2),sampling_model) 
	dlambda_B_nickel2=lambda_B_nickel2[1]-lambda_B_nickel2[0]
	dlambda_V_nickel2=lambda_V_nickel2[1]-lambda_V_nickel2[0]
	dlambda_R_nickel2=lambda_R_nickel2[1]-lambda_R_nickel2[0] 
	dlambda_I_nickel2=lambda_I_nickel2[1]-lambda_I_nickel2[0]

	lambda_B_CSP=np.linspace(min(lambda_B_CSP),max(lambda_B_CSP),sampling_model)
	lambda_V_CSP=np.linspace(min(lambda_V_CSP),max(lambda_V_CSP),sampling_model) 
	lambda_r_CSP=np.linspace(min(lambda_r_CSP),max(lambda_r_CSP),sampling_model) 
	lambda_i_CSP=np.linspace(min(lambda_i_CSP),max(lambda_i_CSP),sampling_model)

	#Regroup all the transmission function
	F_X_func_kait2_tot=[F_B_kait2_func,F_V_kait2_func,F_R_kait2_func,F_I_kait2_func]
	lambda_X_kait2_tot=[lambda_B_kait2,lambda_V_kait2,lambda_R_kait2,lambda_I_kait2]
	F_X_func_kait3_tot=[F_B_kait3_func,F_V_kait3_func,F_R_kait3_func,F_I_kait3_func]
	lambda_X_kait3_tot=[lambda_B_kait3,lambda_V_kait3,lambda_R_kait3,lambda_I_kait3]
	F_X_func_kait4_tot=[F_B_kait4_func,F_V_kait4_func,F_R_kait4_func,F_I_kait4_func]
	lambda_X_kait4_tot=[lambda_B_kait4,lambda_V_kait4,lambda_R_kait4,lambda_I_kait4]

	F_X_func_nickel1_tot=[F_B_nickel1_func,F_V_nickel1_func,F_R_nickel1_func,F_I_nickel1_func]
	lambda_X_nickel1_tot=[lambda_B_nickel1,lambda_V_nickel1,lambda_R_nickel1,lambda_I_nickel1]
	F_X_func_nickel2_tot=[F_B_nickel2_func,F_V_nickel2_func,F_R_nickel2_func,F_I_nickel2_func]
	lambda_X_nickel2_tot=[lambda_B_nickel2,lambda_V_nickel2,lambda_R_nickel2,lambda_I_nickel2]

	F_X_func_tot=[F_X_func_kait2_tot,F_X_func_kait3_tot,F_X_func_kait4_tot,F_X_func_nickel1_tot,F_X_func_nickel2_tot]
	lambda_X_tot=[lambda_X_kait2_tot,lambda_X_kait3_tot,lambda_X_kait4_tot,lambda_X_nickel1_tot,lambda_X_nickel2_tot]

	#Mangling the models at all the epochs for which we have photometry points
	exec(open("mangling.py").read())
	'''
	For each epoch we have a warp model in the obs frame
	wl_model_mang
	fl_model_mang
	'''

	F_X_func_CSP_tot=[F_B_func_CSP,F_V_func_CSP,F_r_func_CSP,F_i_func_CSP]
	lambda_X_CSP_tot=[lambda_B_CSP,lambda_V_CSP,lambda_r_CSP,lambda_i_CSP]

	filtre_obs=['B_kait','V_kait','R_kait','I_kait']
	filtre_rest=['B_swope','V_swope','r_swope','i_swope']
	#Now for each filter we will calculate the AKS: is the difference of the synthetic magnitudes obtained using the 
	#Kait/Nickel photometric system and using the CSP photometric system
	#For Kait we integrate the warped spectrum in the observed frame, for CSP for integrate the warped spectrum in the rest
	#frame and corrected for MW extinction.
	Filter_rest_used=[]
	for i,filter_name in enumerate(MJD):
		
		matrice=np.zeros((np.size(MJD[filter_name]),5))
		x_final=mags[filter_name]	#photo
		err_x_final=emags[filter_name]	#err photo
		phase_x_final=MJD[filter_name]	#epoque
	
		for epo in range(np.size(MJD[filter_name])):
			#Take the good warped spectrum
			ind_warp=np.where(np.array(epoch)==MJD[filter_name][epo])[0][0]
			wavelength_model_obs=wl_model_mang[ind_warp]
			flux_model_mang=fl_model_mang[ind_warp]

	
			#Now we calculate F_spec_model_z(KAIT) 
			#First which system was used (kait2,kait3,kait4,nickel1,nickel2..)
			ind=np.where(MJD[filter_name]==MJD[filter_name][epo])[0][0]
			photo_sys=tel[filter_name][ind]
			ZP_X_kait=ZP[photo_sys][0]
			ind_sys=np.where(np.array(KAIT_sys)==photo_sys)[0][0]
	
			#Synthetic magnitude of the warped spetrum in the observed frame using KAIT system
			F_spec_model_z=interpolate.interp1d(wavelength_model_obs,flux_model_mang)
			
			m_kait=-2.5*np.log10(1.0/(h_erg*c_AA)*integrate.simps(F_spec_model_z(lambda_X_tot[ind_sys][i])*F_X_func_tot[ind_sys][i](lambda_X_tot[ind_sys][i])*lambda_X_tot[ind_sys][i],lambda_X_tot[ind_sys][i]))+ZP_X_kait[i]


			# We put the warped spectrum in the restframe
			wavelength_model_rest=wavelength_model_obs*1.0/(1+z_hel)
			flux_model_rest=flux_model_mang*(1+z_hel)

			# We remove the Milky Way extinction
			flux_model_obs_AvG=ccm_unred(wavelength_model_obs,flux_model_mang,A_V)

			#wavelength_model_obs,flux_model_mang
			wavelength_model_rest=wavelength_model_obs*1.0/(1+z_hel)
			flux_model_rest_AvG=flux_model_obs_AvG*(1+z_hel)

			# Interpolate the restframe spectrum after mang and AvG correction
			F_spec_model_0=interpolate.interp1d(wavelength_model_rest,flux_model_rest_AvG)

			'At low redshift is ok, a photon emitted in the B band will be received in the same band. For example, here, a photon received in the V-band (KAIT) have been emitted in the V rest (CSP). In any case, here we write the general case, I try to find the good band.'

			#First we need to derive the effective wl of the restframe filters, i.e., CSP.

			lambda_eff_B_CSP=effective_wavelength_csp(wavelength_model_rest,flux_model_rest_AvG,'B')
			lambda_eff_V_CSP=effective_wavelength_csp(wavelength_model_rest,flux_model_rest_AvG,'V')
			lambda_eff_r_CSP=effective_wavelength_csp(wavelength_model_rest,flux_model_rest_AvG,'r')
			lambda_eff_i_CSP=effective_wavelength_csp(wavelength_model_rest,flux_model_rest_AvG,'i')
			lambda_eff_csp=[lambda_eff_B_CSP,lambda_eff_V_CSP,lambda_eff_r_CSP,lambda_eff_i_CSP]

			#Look at the observed filt, divide by (1+z) to put it in the rest frame
			new_lambda=lambda_eff_tot_kait[i]*1.0/(1+z_hel) 
			#Look at which CSP filter this corresponds
			ind_lambda_eff_csp=np.where(lambda_eff_csp==find_nearest(lambda_eff_csp,new_lambda))[0][0]
			filtre_rest_choisi=filtre_rest[ind_lambda_eff_csp]

			print ('Obs filter: '+str(filtre_obs[i])+' corresponds to  '+'Restframe filter: '+str(filtre_rest_choisi))
			#if the filter is already taken, we choose the next one
			if filtre_rest_choisi in Filter_rest_used:
				ind_lambda_eff_csp=ind_lambda_eff_csp+1
				filtre_rest_choisi=filtre_rest[ind_lambda_eff_csp]
			else:
				ind_lambda_eff_csp=ind_lambda_eff_csp
				filtre_rest_choisi=filtre_rest[ind_lambda_eff_csp]
			print ('Obs filter: '+str(filtre_obs[i])+' corresponds to  '+'Restframe filter: '+str(filtre_rest_choisi))	

			F_X_func_CSP=F_X_func_CSP_tot[ind_lambda_eff_csp]
			lambda_X_CSP=lambda_X_CSP_tot[ind_lambda_eff_csp]
			ZP_X=ZP['csp'][0][ind_lambda_eff_csp]

			## Synthetic magnitude of the warped spectrum in the rest frame and corrected for AvG using CSP
			m_csp=-2.5*np.log10(1.0/(h_erg*c_AA)*integrate.simps(F_spec_model_0(lambda_X_CSP)*F_X_func_CSP(lambda_X_CSP)*lambda_X_CSP,lambda_X_CSP))+ZP_X

			#AKS correction:
			AKS=m_kait-m_csp
	
			### We save the new photometry in a txt file: SN_filtre.txt
			matrice[epo,0]=phase_x_final[epo]+JD_explo # MJD
			matrice[epo,1]=x_final[epo]       #Obs mag
			matrice[epo,2]=err_x_final[epo]   #err obs mag
			matrice[epo,3]=round(x_final[epo]-AKS,2)	#Mag corrected for AKS
			matrice[epo,4]=AKS	#value AKS
		Filter_rest_used.append(filtre_rest_choisi)
		np.savetxt('AKS_mag/%s_%s_to_%s.dat' %(SN_name,filtre_obs[i],filtre_rest_choisi),matrice ,fmt='%s')

'''

Here, for 4 SNe, we have data obtained using KAIT and CSP. We can easily check if our technique works by comparing
the magnitudes obtained after corrected the KAIT photometry and those obtained with CSP.
The initial KAIT magnitudes should be different than CSP (and mainly I where the filters are very different) and the final AKS
 KAIT photometry should be more similar to CSP
'''



#Photometry from CSP for the 4 SNe in common with KAIT
exec(open("photo_CSP_AKS.py").read())

#Photometry from KAIT after AKS
exec(open("photo_KAIT_AKS.py").read())

# Figure to see the AKS correction
exec(open("figure.py").read())


