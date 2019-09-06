import numpy as np # numerical tools
from scipy import integrate
from scipy import interpolate
c_light=299792.458#in km/s

#Find nearest value
def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]
#### DATA SN
def get_SN_info(targetname):
	data_sn=np.loadtxt('Info_SNe_KAIT.txt',usecols=[1,2,3,4,5,6,7]).transpose()
	name_SN_kait=np.array(np.genfromtxt('Info_SNe_KAIT.txt',usecols=[0],dtype='str'))

	ind_SN=np.where(np.array(name_SN_kait)==targetname)[0][0]

	A_V=data_sn[0][ind_SN]
	z_hel=data_sn[1][ind_SN]*1.0/c_light
	err_z_hel=data_sn[2][ind_SN]*1.0/c_light
	JD_explo=data_sn[3][ind_SN]
	err_JD_explo=data_sn[4][ind_SN]
	z_cmb=data_sn[5][ind_SN]*1.0/c_light
	err_z_cmb=data_sn[6][ind_SN]*1.0/c_light

	return A_V,z_hel,err_z_hel,JD_explo,err_JD_explo,z_cmb,err_z_cmb

#Get SN photometry
def get_sn(targetname):
	data_sn=open('Nat_KAIT/%s.txt'%targetname,'r')
	lines = data_sn.readlines()
	fields = lines[0].split()
	ind_B=np.where(np.array(fields)=='B')[0][0]
	ind_V=np.where(np.array(fields)=='V')[0][0]
	ind_R=np.where(np.array(fields)=='R')[0][0]
	ind_I=np.where(np.array(fields)=='I')[0][0]
	MJD = {}
	mags = {}
	emags = {}
	tel = {}
	for i in range(4):
		this_filter = ['B','V','R','I']
		MJD[this_filter[i]] = []
		mags[this_filter[i]] = []
		emags[this_filter[i]] = []
		tel[this_filter[i]] = []
	for j in range(np.size(lines)):
		if (j!=0):

			if ((lines[j].split()[ind_B+1])<'0.8') and ((lines[j].split()[ind_B+1])!='NaN') and ((lines[j].split()[0][0])!='#'): 
				mags['B'].append(float(lines[j].split()[ind_B]))
				emags['B'].append(float(lines[j].split()[ind_B+1]))
				MJD['B'].append(float(lines[j].split()[1]))
				tel['B'].append(lines[j].split()[3])
			if ((lines[j].split()[ind_V+1])<'0.8') and ((lines[j].split()[ind_V+1])!='NaN')and ((lines[j].split()[0][0])!='#'): 
				mags['V'].append(float(lines[j].split()[ind_V]))
				emags['V'].append(float(lines[j].split()[ind_V+1]))
				MJD['V'].append(float(lines[j].split()[1]))
				tel['V'].append(lines[j].split()[3])
			if ((lines[j].split()[ind_R+1])<'0.8') and ((lines[j].split()[ind_R+1])!='NaN') and ((lines[j].split()[0][0])!='#'): 
				mags['R'].append(float(lines[j].split()[ind_R]))
				emags['R'].append(float(lines[j].split()[ind_R+1]))
				MJD['R'].append(float(lines[j].split()[1]))
				tel['R'].append(lines[j].split()[3])

			if ((lines[j].split()[ind_I+1])<'0.8') and ((lines[j].split()[ind_I+1])!='NaN') and ((lines[j].split()[0][0])!='#'): 
				mags['I'].append(float(lines[j].split()[ind_I]))
				emags['I'].append(float(lines[j].split()[ind_I+1]))
				MJD['I'].append(float(lines[j].split()[1]))
				tel['I'].append(lines[j].split()[3])
	for f in MJD:
		MJD[f],mags[f],emags[f],tel[f]=zip(*sorted(zip(MJD[f],mags[f],emags[f],tel[f])))	
		MJD[f] = np.array(MJD[f])
		mags[f] = np.array(mags[f])
		emags[f] = np.array(emags[f])
		tel[f] = np.array(tel[f])

	return MJD,mags,emags,tel

#Linear interpolation of the magnitude
def inter_mag(MJD,mags,emags):

	B_band=interpolate.interp1d(MJD['B'],mags['B'])
	B_band_plus=interpolate.interp1d(MJD['B'],mags['B']+emags['B'])
	V_band=interpolate.interp1d(MJD['V'],mags['V'])
	V_band_plus=interpolate.interp1d(MJD['V'],mags['V']+emags['V'])
	if np.size(MJD['R'])>0:
		R_band=interpolate.interp1d(MJD['R'],mags['R'])
		R_band_plus=interpolate.interp1d(MJD['R'],mags['R']+emags['R'])
	else:
		R_band=[]
		R_band_plus=[]
	I_band=interpolate.interp1d(MJD['I'],mags['I'])
	I_band_plus=interpolate.interp1d(MJD['I'],mags['I']+emags['I'])

	return B_band,B_band_plus,V_band,V_band_plus,R_band,R_band_plus,I_band,I_band_plus


#Derive for each CSP filter the effective wavelength
def effective_wavelength_csp(lam_spec,flux_spec,filter_name):
	### Each transmission function ###########
	trans_u=np.loadtxt('Filters/CSP/u_swope.txt')
	lambda_u=trans_u[:,0]
	s_u=trans_u[:,1]
	trans_g=np.loadtxt('Filters/CSP/g_swope.txt')
	lambda_g=trans_g[:,0]
	s_g=trans_g[:,1]
	trans_r=np.loadtxt('Filters/CSP/r_swope.txt')
	lambda_r=trans_r[:,0]
	s_r=trans_r[:,1]
	trans_i=np.loadtxt('Filters/CSP/i_swope.txt')
	lambda_i=trans_i[:,0]
	s_i=trans_i[:,1]
	trans_V=np.loadtxt('Filters/CSP/V_swope.txt')
	lambda_V=trans_V[:,0]
	s_V=trans_V[:,1]
	trans_B=np.loadtxt('Filters/CSP/B_swope.txt')
	lambda_B=trans_B[:,0]
	s_B=trans_B[:,1]

	F_u_func=interpolate.interp1d(lambda_u,s_u) #interpolation Filtre u
	F_B_func=interpolate.interp1d(lambda_B,s_B) #interpolation Filtre B
	F_V_func=interpolate.interp1d(lambda_V,s_V) #interpolation Filtre V
	F_g_func=interpolate.interp1d(lambda_g,s_g) 
	F_r_func=interpolate.interp1d(lambda_r,s_r) #interpolation Filtre t
	F_i_func=interpolate.interp1d(lambda_i,s_i) #interpolation Filtre i
	N_pt=3000
	lambda_u=np.linspace(min(lambda_u),max(lambda_u),N_pt) 
	lambda_B=np.linspace(min(lambda_B),max(lambda_B),N_pt) 
	lambda_V=np.linspace(min(lambda_V),max(lambda_V),N_pt) 
	lambda_g=np.linspace(min(lambda_g),max(lambda_g),N_pt) 
	lambda_r=np.linspace(min(lambda_r),max(lambda_r),N_pt) 
	lambda_i=np.linspace(min(lambda_i),max(lambda_i),N_pt)
	if filter_name==str('u'):
		F_filter_func=interpolate.interp1d(lambda_u,F_u_func(lambda_u)) #interpolation Filtre B
		lam_filter=lambda_u
	if filter_name==str('B'):
		F_filter_func=interpolate.interp1d(lambda_B,F_B_func(lambda_B)) #interpolation Filtre B
		lam_filter=lambda_B
	if filter_name==str('g'):	
		F_filter_func=interpolate.interp1d(lambda_g,F_g_func(lambda_g)) 
		lam_filter=lambda_g
	if filter_name==str('V'):	
		F_filter_func=interpolate.interp1d(lambda_V,F_V_func(lambda_V)) #interpolation Filtre V
		lam_filter=lambda_V
	if filter_name==str('r'):	
		F_filter_func=interpolate.interp1d(lambda_r,F_r_func(lambda_r)) #interpolation Filtre r
		lam_filter=lambda_r
	if filter_name==str('i'):	
		F_filter_func=interpolate.interp1d(lambda_i,F_i_func(lambda_i)) #interpolation Filtre i
		lam_filter=lambda_i
	# interpolation spectre
	F_spec=interpolate.interp1d(lam_spec,flux_spec)
	# New wavelength vector with wavelength of filter + spectrum
	wavelength_to_interpolate=np.concatenate([lam_spec,lam_filter])
	# Sort the wavelength
	wavelength_to_interpolate.sort()
	# We select only the wavelenght in the filter
	wavelength_to_interpolate_2=wavelength_to_interpolate[(wavelength_to_interpolate>min(lam_filter)) & (wavelength_to_interpolate<max(lam_filter))]
	# We calculate the filter response
	interpolate_filter_response=F_filter_func(wavelength_to_interpolate_2)
	# We calculate SEDter
	SED_inside_filter=F_spec(wavelength_to_interpolate_2)
	# num=f*s*lambda
	num=SED_inside_filter*interpolate_filter_response*wavelength_to_interpolate_2*wavelength_to_interpolate_2
	# num=f*s
	dem=SED_inside_filter*interpolate_filter_response*wavelength_to_interpolate_2
	# integral de num / integral de dem
	lambda_eff_filter=np.trapz(num)*1.0/np.trapz(dem)	
	return lambda_eff_filter

def effective_wavelength_KAIT(lam_spec,flux_spec,filter_name):

	### KAIT 2 ###########
	trans_B_kait2=np.loadtxt('Filters/KAIT_NICKEL/B_kait2.txt')
	lambda_B_kait2=trans_B_kait2[:,0]
	s_B_kait2=trans_B_kait2[:,1]
	trans_V_kait2=np.loadtxt('Filters/KAIT_NICKEL/V_kait2.txt')
	lambda_V_kait2=trans_V_kait2[:,0]
	s_V_kait2=trans_V_kait2[:,1]
	trans_R_kait2=np.loadtxt('Filters/KAIT_NICKEL/R_kait2.txt')
	lambda_R_kait2=trans_R_kait2[:,0]
	s_R_kait2=trans_R_kait2[:,1]
	trans_I_kait2=np.loadtxt('Filters/KAIT_NICKEL/I_kait2.txt')
	lambda_I_kait2=trans_I_kait2[:,0]
	s_I_kait2=trans_I_kait2[:,1]
	dlambda_B_kait2=lambda_B_kait2[1]-lambda_B_kait2[0]
	dlambda_V_kait2=lambda_V_kait2[1]-lambda_V_kait2[0]
	dlambda_R_kait2=lambda_R_kait2[1]-lambda_R_kait2[0]
	dlambda_I_kait2=lambda_I_kait2[1]-lambda_I_kait2[0]

	### KAIT 3 ###########
	trans_B_kait3=np.loadtxt('Filters/KAIT_NICKEL/B_kait3.txt')
	lambda_B_kait3=trans_B_kait3[:,0]
	s_B_kait3=trans_B_kait3[:,1]
	trans_V_kait3=np.loadtxt('Filters/KAIT_NICKEL/V_kait3.txt')
	lambda_V_kait3=trans_V_kait3[:,0]
	s_V_kait3=trans_V_kait3[:,1]
	trans_R_kait3=np.loadtxt('Filters/KAIT_NICKEL/R_kait3.txt')
	lambda_R_kait3=trans_R_kait3[:,0]
	s_R_kait3=trans_R_kait3[:,1]
	trans_I_kait3=np.loadtxt('Filters/KAIT_NICKEL/I_kait3.txt')
	lambda_I_kait3=trans_I_kait3[:,0]
	s_I_kait3=trans_I_kait3[:,1]

	dlambda_B_kait3=lambda_B_kait3[1]-lambda_B_kait3[0]
	dlambda_V_kait3=lambda_V_kait3[1]-lambda_V_kait3[0]
	dlambda_R_kait3=lambda_R_kait3[1]-lambda_R_kait3[0]
	dlambda_I_kait3=lambda_I_kait3[1]-lambda_I_kait3[0]

	### KAIT 4 ###########
	trans_B_kait4=np.loadtxt('Filters/KAIT_NICKEL/B_kait4.txt')
	lambda_B_kait4=trans_B_kait4[:,0]
	s_B_kait4=trans_B_kait4[:,1]
	trans_V_kait4=np.loadtxt('Filters/KAIT_NICKEL/V_kait4.txt')
	lambda_V_kait4=trans_V_kait4[:,0]
	s_V_kait4=trans_V_kait4[:,1]
	trans_R_kait4=np.loadtxt('Filters/KAIT_NICKEL/R_kait4.txt')
	lambda_R_kait4=trans_R_kait4[:,0]
	s_R_kait4=trans_R_kait4[:,1]
	trans_I_kait4=np.loadtxt('Filters/KAIT_NICKEL/I_kait4.txt')
	lambda_I_kait4=trans_I_kait4[:,0]
	s_I_kait4=trans_I_kait4[:,1]

	dlambda_B_kait4=lambda_B_kait4[1]-lambda_B_kait4[0]
	dlambda_V_kait4=lambda_V_kait4[1]-lambda_V_kait4[0]
	dlambda_R_kait4=lambda_R_kait4[1]-lambda_R_kait4[0]
	dlambda_I_kait4=lambda_I_kait4[1]-lambda_I_kait4[0]


	### Nickel 1 ###########
	trans_B_nickel1=np.loadtxt('Filters/KAIT_NICKEL/B_nickel1.txt')
	lambda_B_nickel1=trans_B_nickel1[:,0]
	s_B_nickel1=trans_B_nickel1[:,1]
	trans_V_nickel1=np.loadtxt('Filters/KAIT_NICKEL/V_nickel1.txt')
	lambda_V_nickel1=trans_V_nickel1[:,0]
	s_V_nickel1=trans_V_nickel1[:,1]
	trans_R_nickel1=np.loadtxt('Filters/KAIT_NICKEL/R_nickel1.txt')
	lambda_R_nickel1=trans_R_nickel1[:,0]
	s_R_nickel1=trans_R_nickel1[:,1]
	trans_I_nickel1=np.loadtxt('Filters/KAIT_NICKEL/I_nickel1.txt')
	lambda_I_nickel1=trans_I_nickel1[:,0]
	s_I_nickel1=trans_I_nickel1[:,1]

	dlambda_B_nicke1l=lambda_B_nickel1[1]-lambda_B_nickel1[0]
	dlambda_V_nickel1=lambda_V_nickel1[1]-lambda_V_nickel1[0]
	dlambda_R_nickel1=lambda_R_nickel1[1]-lambda_R_nickel1[0]
	dlambda_I_nickel1=lambda_I_nickel1[1]-lambda_I_nickel1[0]

	### Nickel 2 ###########
	trans_B_nickel2=np.loadtxt('Filters/KAIT_NICKEL/B_nickel2.txt')
	lambda_B_nickel2=trans_B_nickel2[:,0]
	s_B_nickel2=trans_B_nickel2[:,1]
	trans_V_nickel2=np.loadtxt('Filters/KAIT_NICKEL/V_nickel2.txt')
	lambda_V_nickel2=trans_V_nickel2[:,0]
	s_V_nickel2=trans_V_nickel2[:,1]
	trans_R_nickel2=np.loadtxt('Filters/KAIT_NICKEL/R_nickel2.txt')
	lambda_R_nickel2=trans_R_nickel2[:,0]
	s_R_nickel2=trans_R_nickel2[:,1]
	trans_I_nickel2=np.loadtxt('Filters/KAIT_NICKEL/I_nickel2.txt')
	lambda_I_nickel2=trans_I_nickel2[:,0]
	s_I_nickel2=trans_I_nickel2[:,1]

	dlambda_B_nickel2=lambda_B_nickel2[1]-lambda_B_nickel2[0]
	dlambda_V_nickel2=lambda_V_nickel2[1]-lambda_V_nickel2[0]
	dlambda_R_nickel2=lambda_R_nickel2[1]-lambda_R_nickel2[0]
	dlambda_I_nickel2=lambda_I_nickel2[1]-lambda_I_nickel2[0]

	F_B_kait2_func=interpolate.interp1d(lambda_B_kait2,s_B_kait2) 
	F_V_kait2_func=interpolate.interp1d(lambda_V_kait2,s_V_kait2) 
	F_R_kait2_func=interpolate.interp1d(lambda_R_kait2,s_R_kait2) 
	F_I_kait2_func=interpolate.interp1d(lambda_I_kait2,s_I_kait2) 
	F_B_kait3_func=interpolate.interp1d(lambda_B_kait3,s_B_kait3) 
	F_V_kait3_func=interpolate.interp1d(lambda_V_kait3,s_V_kait3) 
	F_R_kait3_func=interpolate.interp1d(lambda_R_kait3,s_R_kait3) 
	F_I_kait3_func=interpolate.interp1d(lambda_I_kait3,s_I_kait3) 
	F_B_kait4_func=interpolate.interp1d(lambda_B_kait4,s_B_kait4) 
	F_V_kait4_func=interpolate.interp1d(lambda_V_kait4,s_V_kait4) 
	F_R_kait4_func=interpolate.interp1d(lambda_R_kait4,s_R_kait4) 
	F_I_kait4_func=interpolate.interp1d(lambda_I_kait4,s_I_kait4) 
	F_B_nickel1_func=interpolate.interp1d(lambda_B_nickel1,s_B_nickel1) 
	F_V_nickel1_func=interpolate.interp1d(lambda_V_nickel1,s_V_nickel1) 
	F_R_nickel1_func=interpolate.interp1d(lambda_R_nickel1,s_R_nickel1) 
	F_I_nickel1_func=interpolate.interp1d(lambda_I_nickel1,s_I_nickel1) 
	F_B_nickel2_func=interpolate.interp1d(lambda_B_nickel2,s_B_nickel2) 
	F_V_nickel2_func=interpolate.interp1d(lambda_V_nickel2,s_V_nickel2) 
	F_R_nickel2_func=interpolate.interp1d(lambda_R_nickel2,s_R_nickel2) 
	F_I_nickel2_func=interpolate.interp1d(lambda_I_nickel2,s_I_nickel2) 	
	N_pt=5000 
	lambda_B_kait2=np.linspace(min(lambda_B_kait2),max(lambda_B_kait2),N_pt) 
	lambda_V_kait2=np.linspace(min(lambda_V_kait2),max(lambda_V_kait2),N_pt) 
	lambda_R_kait2=np.linspace(min(lambda_R_kait2),max(lambda_R_kait2),N_pt) 
	lambda_I_kait2=np.linspace(min(lambda_I_kait2),max(lambda_I_kait2),N_pt) 
	lambda_B_kait3=np.linspace(min(lambda_B_kait3),max(lambda_B_kait3),N_pt) 
	lambda_V_kait3=np.linspace(min(lambda_V_kait3),max(lambda_V_kait3),N_pt) 
	lambda_R_kait3=np.linspace(min(lambda_R_kait3),max(lambda_R_kait3),N_pt) 
	lambda_I_kait3=np.linspace(min(lambda_I_kait3),max(lambda_I_kait3),N_pt) 
	lambda_B_kait4=np.linspace(min(lambda_B_kait4),max(lambda_B_kait4),N_pt) 
	lambda_V_kait4=np.linspace(min(lambda_V_kait4),max(lambda_V_kait4),N_pt) 
	lambda_R_kait4=np.linspace(min(lambda_R_kait4),max(lambda_R_kait4),N_pt) 
	lambda_I_kait4=np.linspace(min(lambda_I_kait4),max(lambda_I_kait4),N_pt) 
	lambda_B_nickel1=np.linspace(min(lambda_B_nickel1),max(lambda_B_nickel1),N_pt) 
	lambda_V_nickel1=np.linspace(min(lambda_V_nickel1),max(lambda_V_nickel1),N_pt) 
	lambda_R_nickel1=np.linspace(min(lambda_R_nickel1),max(lambda_R_nickel1),N_pt) 
	lambda_I_nickel1=np.linspace(min(lambda_I_nickel1),max(lambda_I_nickel1),N_pt) 
	lambda_B_nickel2=np.linspace(min(lambda_B_nickel2),max(lambda_B_nickel2),N_pt) 
	lambda_V_nickel2=np.linspace(min(lambda_V_nickel2),max(lambda_V_nickel2),N_pt) 
	lambda_R_nickel2=np.linspace(min(lambda_R_nickel2),max(lambda_R_nickel2),N_pt) 
	lambda_I_nickel2=np.linspace(min(lambda_I_nickel2),max(lambda_I_nickel2),N_pt) 
	
	if filter_name==str('Bkait2'):	
		F_filter_func=interpolate.interp1d(lambda_B_kait2,F_B_kait2_func(lambda_B_kait2)) 
		lam_filter=lambda_B_kait2
	if filter_name==str('Vkait2'):	
		F_filter_func=interpolate.interp1d(lambda_V_kait2,F_V_kait2_func(lambda_V_kait2)) 
		lam_filter=lambda_V_kait2
	if filter_name==str('Rkait2'):	
		F_filter_func=interpolate.interp1d(lambda_R_kait2,F_R_kait2_func(lambda_R_kait2)) 
		lam_filter=lambda_R_kait2
	if filter_name==str('Ikait2'):	
		F_filter_func=interpolate.interp1d(lambda_I_kait2,F_I_kait2_func(lambda_I_kait2)) 
		lam_filter=lambda_I_kait2

	if filter_name==str('Bkait3'):	
		F_filter_func=interpolate.interp1d(lambda_B_kait3,F_B_kait3_func(lambda_B_kait3)) 
		lam_filter=lambda_B_kait3
	if filter_name==str('Vkait3'):	
		F_filter_func=interpolate.interp1d(lambda_V_kait3,F_V_kait3_func(lambda_V_kait3)) 
		lam_filter=lambda_V_kait3
	if filter_name==str('Rkait3'):	
		F_filter_func=interpolate.interp1d(lambda_R_kait3,F_R_kait3_func(lambda_R_kait3)) 
		lam_filter=lambda_R_kait3
	if filter_name==str('Ikait3'):	
		F_filter_func=interpolate.interp1d(lambda_I_kait3,F_I_kait3_func(lambda_I_kait3)) 
		lam_filter=lambda_I_kait3

	if filter_name==str('Bkait4'):	
		F_filter_func=interpolate.interp1d(lambda_B_kait4,F_B_kait4_func(lambda_B_kait4)) 
		lam_filter=lambda_B_kait4
	if filter_name==str('Vkait4'):	
		F_filter_func=interpolate.interp1d(lambda_V_kait4,F_V_kait4_func(lambda_V_kait4)) 
		lam_filter=lambda_V_kait4
	if filter_name==str('Rkait4'):	
		F_filter_func=interpolate.interp1d(lambda_R_kait4,F_R_kait4_func(lambda_R_kait4)) 
		lam_filter=lambda_R_kait4
	if filter_name==str('Ikait4'):	
		F_filter_func=interpolate.interp1d(lambda_I_kait4,F_I_kait4_func(lambda_I_kait4)) 
		lam_filter=lambda_I_kait4
	if filter_name==str('Bnickel1'):	
		F_filter_func=interpolate.interp1d(lambda_B_nickel1,F_B_nickel1_func(lambda_B_nickel1)) 
		lam_filter=lambda_B_nickel1
	if filter_name==str('Vnickel1'):	
		F_filter_func=interpolate.interp1d(lambda_V_nickel1,F_V_nickel1_func(lambda_V_nickel1)) 
		lam_filter=lambda_V_nickel1
	if filter_name==str('Rnickel1'):	
		F_filter_func=interpolate.interp1d(lambda_R_nickel1,F_R_nickel1_func(lambda_R_nickel1)) 
		lam_filter=lambda_R_nickel1
	if filter_name==str('Inickel1'):	
		F_filter_func=interpolate.interp1d(lambda_I_nickel1,F_I_nickel1_func(lambda_I_nickel1)) 
		lam_filter=lambda_I_nickel1
	if filter_name==str('Bnickel2'):	
		F_filter_func=interpolate.interp1d(lambda_B_nickel2,F_B_nickel2_func(lambda_B_nickel2)) 
		lam_filter=lambda_B_nickel2
	if filter_name==str('Vnickel2'):	
		F_filter_func=interpolate.interp1d(lambda_V_nickel2,F_V_nickel2_func(lambda_V_nickel2)) 
		lam_filter=lambda_V_nickel2
	if filter_name==str('Rnickel2'):	
		F_filter_func=interpolate.interp1d(lambda_R_nickel2,F_R_nickel2_func(lambda_R_nickel2)) 
		lam_filter=lambda_R_nickel2
	if filter_name==str('Inickel2'):	
		F_filter_func=interpolate.interp1d(lambda_I_nickel2,F_I_nickel2_func(lambda_I_nickel2)) 
		lam_filter=lambda_I_nickel2


	# interpolation spectre
	F_spec=interpolate.interp1d(lam_spec,flux_spec)
	# New wavelength vector with wavelength of filter + spectrum
	wavelength_to_interpolate=np.concatenate([lam_spec,lam_filter])
	# Sort the wavelength
	wavelength_to_interpolate.sort()
	# We select only the wavelenght in the filter
	wavelength_to_interpolate_2=wavelength_to_interpolate[(wavelength_to_interpolate>min(lam_filter)) & (wavelength_to_interpolate<max(lam_filter))]
	# We calculate the filter response
	interpolate_filter_response=F_filter_func(wavelength_to_interpolate_2)
	# We calculate SEDter
	SED_inside_filter=F_spec(wavelength_to_interpolate_2)
	# num=f*s*lambda
	num=SED_inside_filter*interpolate_filter_response*wavelength_to_interpolate_2*wavelength_to_interpolate_2
	# num=f*s
	dem=SED_inside_filter*interpolate_filter_response*wavelength_to_interpolate_2
	# integral de num / integral de dem
	lambda_eff_filter=np.trapz(num)*1.0/np.trapz(dem)	

	return lambda_eff_filter


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
