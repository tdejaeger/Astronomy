###### On charge la photometrie CSP apres correction AvG,KC, KS #######
filtre=['B_swope','V_swope','r_swope','i_swope']
filtre_kait=['B_kait','V_kait','R_kait','I_kait']
phase_X_final_kait=[] #phase finale apres explosion et corrigee par z
X_final_kait=[]  #mag 
X_AKS_final_kait=[]  #mag AKS
err_X_AKS_final_kait=[] #err_mag AKS
AKS_X_kait=[]		#AKS correction
for filtre_name in range(np.size(filtre)):
	# Probleme mang KC pour 0 40, 46 bk sub
	mag_txt_kait=np.loadtxt('AKS_mag/%s_%s_to_%s.dat' %(SN_name,filtre_kait[filtre_name],filtre[filtre_name])).transpose()
	
	phase_X_final_kait.append(np.around((mag_txt_kait[0][mag_txt_kait[3]!=0.0]-JD_explo)*1.0/(1+z_hel), decimals=2))  
	X_final_kait.append(mag_txt_kait[1][mag_txt_kait[3]!=0.0])
	X_AKS_final_kait.append(mag_txt_kait[3][mag_txt_kait[3]!=0.0])   
	err_X_AKS_final_kait.append(mag_txt_kait[2][mag_txt_kait[3]!=0.0]) 
	AKS_X_kait.append(mag_txt_kait[3][mag_txt_kait[3]!=0.0])        



