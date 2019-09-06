# on charge la table avec toutes les donnees de SNe
data_sn=np.loadtxt('CSP_photo/Data_all_SNe_CSP.dat',usecols=[1,2,3,4,5,6,7,8,9,10,11]).transpose()

name_SN_CSP=np.array(np.genfromtxt('CSP_photo/Data_all_SNe_CSP.dat',usecols=[0],dtype='str'))

ind_SN=np.where(np.array(name_SN_CSP)==SN_name)[0][0]

c_light=299792.458#in km/s
z_hel_CSP=data_sn[1][ind_SN]*1.0/c_light
err_z_hel_CSP=data_sn[2][ind_SN]*1.0/c_light
JD_explo_CSP=data_sn[3][ind_SN]
err_JD_explo_CSP=data_sn[4][ind_SN]

SN_name=SN_name
###### For a previous work, we have already done the A and K correction for these SNe
###### Here we download the data
filtre=['B','V','r','i']
phase_X_final_csp=[] #final phase corrected for z
X_final_csp=[]	     #mag raw
err_X_final_csp=[]   #err mag raw
X_AKS_final_csp=[]  #mag AKS
AKS_X_csp=[]		#AKS correction
for filtre_name in range(np.size(filtre)):
	mag_txt_csp=np.loadtxt('CSP_photo/%s_%s.dat' %(SN_name,filtre[filtre_name])).transpose()

	phase_X_final_csp.append(np.around(mag_txt_csp[0][mag_txt_csp[1]!=999.9]*1.0/(1+z_hel_CSP), decimals=2))  
	X_final_csp.append(mag_txt_csp[4][mag_txt_csp[1]!=999.9])	       
	err_X_final_csp.append(mag_txt_csp[2][mag_txt_csp[1]!=999.9])    
	X_AKS_final_csp.append(mag_txt_csp[1][mag_txt_csp[1]!=999.9])   
	AKS_X_csp.append(mag_txt_csp[3][mag_txt_csp[1]!=999.9])         
