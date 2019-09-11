# Define unique epoch
epoch = np.array(list(set().union(MJD['B'],MJD['V'],MJD['R'],MJD['I'])))
epoch.sort()

wl_model_mang=[]
fl_model_mang=[]
ep_model_mang=[]
KAIT_sys=['kait2','kait3','kait4','nickel1','nickel2']

for ep in range(np.size(epoch)):
	#Which filters are available at this epoch?
	epoque_photo=epoch[ep]
	print('%s/%s'%(ep,np.size(epoch)))
	if (epoque_photo>=min(MJD['B'])-1) and (epoque_photo<=max(MJD['B'])+1):
		filter_B=str('B')
	else:
		filter_B=str('')
	if (epoque_photo>=min(MJD['V'])-1) and (epoque_photo<=max(MJD['V'])+1):
		filter_V=str('V')
	else:
		filter_V=str('')
	if (np.size(MJD['R']) >1) and (epoque_photo>=min(MJD['R'])-1) and (epoque_photo<=max(MJD['R'])+1):
		filter_R=str('R')
	else:
		filter_R=str('')
	if (np.size(MJD['I']) >1) and (epoque_photo>=min(MJD['I'])-1) and (epoque_photo<=max(MJD['I'])+1):
		filter_I=str('I')
	else:
		filter_I=str('')
	filter_available=filter_B+filter_V+filter_R+filter_I

	#Which photometric system was used: Kait2, Kait3, Kait4, Nickel1 Nickel2
	if epoque_photo in MJD['V']:
		ind=np.where(MJD['V']==epoque_photo)[0][0]
		photo_sys=tel['V'][ind]
	elif epoque_photo in MJD['R']:
		ind=np.where(MJD['R']==epoque_photo)[0][0]
		photo_sys=tel['R'][ind]
	elif epoque_photo in MJD['I']:
		ind=np.where(MJD['I']==epoque_photo)[0][0]
		photo_sys=tel['I'][ind]

	ind_sys=np.where(np.array(KAIT_sys)==photo_sys)[0][0]
	ZP_X_kait=ZP[photo_sys][0] #Zero point of the photo sys used

	#We select the closest model in epoch (divided by (1+z))
	ind_ep_mod=np.where(epoch_mod==find_nearest(np.array(epoch_mod),epoque_photo*1.0/(1+z_hel)))[0][0]
	txt=np.loadtxt('Models/%s.fl'%Model[ind_ep_mod]).transpose()
	wl_mod=(txt[0][(txt[0]<16000) & (txt[0]>1500)]) #no need of the IR part of the spectrum
	fl_mod=(txt[1][(txt[0]<16000) & (txt[0]>1500)])
	#We put the rest frame model into the observed frame 
	wl_mod_obs=(wl_mod)*(1+z_hel)
	fl_mod_obs=(fl_mod)*1.0/(1+z_hel)
	F_mod_func_obs=interpolate.interp1d(wl_mod_obs,fl_mod_obs)

	#We derive the effective wavelength for KAIT
	lambda_eff_B_kait=effective_wavelength_KAIT(wl_mod,fl_mod,'B%s'%photo_sys)
	lambda_eff_V_kait=effective_wavelength_KAIT(wl_mod,fl_mod,'V%s'%photo_sys)
	lambda_eff_R_kait=effective_wavelength_KAIT(wl_mod,fl_mod,'R%s'%photo_sys)
	lambda_eff_I_kait=effective_wavelength_KAIT(wl_mod,fl_mod,'I%s'%photo_sys)
	lambda_eff_tot_kait=[lambda_eff_B_kait,lambda_eff_V_kait,lambda_eff_R_kait,lambda_eff_I_kait]

	#We normalise the model using the V band. It is always better to first normalise using one band, therefore, the other correction will be smallest.
	m_V_kait=-2.5*np.log10(1.0/(h_erg*c_AA)*integrate.simps(F_mod_func_obs(lambda_X_tot[ind_sys][1])*F_X_func_tot[ind_sys][1](lambda_X_tot[ind_sys][1])*lambda_X_tot[ind_sys][1],lambda_X_tot[ind_sys][1]))+ZP_X_kait[1]

	if str('V') in filter_available :
		if (epoque_photo<max(MJD['V'])) and (epoque_photo>min(MJD['V'])):
			mag_V_photo=V_band(epoque_photo)
		elif (epoque_photo>=max(MJD['V'])): 
			mag_V_photo=V_band(max(MJD['V'])-0.0001)
		elif  (epoque_photo<=min(MJD['V'])): 
			mag_V_photo=V_band(min(MJD['V'])+0.0001)

	coeff_norm=(10**(-0.4*mag_V_photo))*1.0/(10**(-0.4*m_V_kait)) #difference in mag converted in diff in flux
	fl_mod_obs=fl_mod_obs*coeff_norm #Apply the normalisation to all the spectrum
	F_mod_func_obs=interpolate.interp1d(wl_mod_obs,fl_mod_obs)

	# We can now derive the synthetic magnitude after normalisation
	# B band KAIT
	m_B_kait=-2.5*np.log10(1.0/(h_erg*c_AA)*integrate.simps(F_mod_func_obs(lambda_X_tot[ind_sys][0])*F_X_func_tot[ind_sys][0](lambda_X_tot[ind_sys][0])*lambda_X_tot[ind_sys][0],lambda_X_tot[ind_sys][0]))+ZP_X_kait[0]
	# V band KAIT
	m_V_kait=-2.5*np.log10(1.0/(h_erg*c_AA)*integrate.simps(F_mod_func_obs(lambda_X_tot[ind_sys][1])*F_X_func_tot[ind_sys][1](lambda_X_tot[ind_sys][1])*lambda_X_tot[ind_sys][1],lambda_X_tot[ind_sys][1]))+ZP_X_kait[1]
	# R band KAIT
	m_R_kait=-2.5*np.log10(1.0/(h_erg*c_AA)*integrate.simps(F_mod_func_obs(lambda_X_tot[ind_sys][2])*F_X_func_tot[ind_sys][2](lambda_X_tot[ind_sys][2])*lambda_X_tot[ind_sys][2],lambda_X_tot[ind_sys][2]))+ZP_X_kait[2]
	# I band KAIT
	m_I_kait=-2.5*np.log10(1.0/(h_erg*c_AA)*integrate.simps(F_mod_func_obs(lambda_X_tot[ind_sys][3])*F_X_func_tot[ind_sys][3](lambda_X_tot[ind_sys][3])*lambda_X_tot[ind_sys][3],lambda_X_tot[ind_sys][3]))+ZP_X_kait[3]

	#For each filter we look at the observed photometry and compare with the synthetic to obtain the offset
	f_filter=np.zeros(4)
	mag_available=np.zeros(4)
	if str('B') in filter_available:
		if (epoque_photo<max(MJD['B'])) and (epoque_photo>min(MJD['B'])):
			mag_B_photo=B_band(epoque_photo)
		elif (epoque_photo>=max(MJD['B'])): 
			mag_B_photo=B_band(max(MJD['B'])-0.0001)
		elif  (epoque_photo<=min(MJD['B'])): 
			mag_B_photo=B_band(min(MJD['B'])+0.0001)
		delta_M_B=mag_B_photo-m_B_kait
		f_B=10**(-0.4*delta_M_B)
		f_filter[0]=f_B
		mag_available[0]=mag_B_photo
	else:
		f_filter[0]=0
		mag_available[0]=0

	if str('V') in filter_available :
		if (epoque_photo<=max(MJD['V'])) and (epoque_photo>=min(MJD['V'])):
			mag_V_photo=V_band(epoque_photo)
		elif (epoque_photo>max(MJD['V'])): 
			mag_V_photo=V_band(max(MJD['V'])-0.0001)
		elif  (epoque_photo<min(MJD['V'])): 
			mag_V_photo=V_band(min(MJD['V'])+0.0001)
		delta_M_V=mag_V_photo-m_V_kait
		f_V=10**(-0.4*delta_M_V)
		f_filter[1]=f_V
		mag_available[1]=mag_V_photo
	else:
		f_filter[1]=0
		mag_available[1]=0
	if str('R') in filter_available :
		if (epoque_photo<=max(MJD['R'])) and (epoque_photo>=min(MJD['R'])):
			mag_R_photo=R_band(epoque_photo)
		elif (epoque_photo>max(MJD['R'])): 
			mag_R_photo=R_band(max(MJD['R'])-0.0001)
		elif  (epoque_photo<min(MJD['R'])): 
			mag_R_photo=R_band(min(MJD['R'])+0.0001)
		delta_M_R=mag_R_photo-m_R_kait
		f_R=10**(-0.4*delta_M_R)
		f_filter[2]=f_R
		mag_available[2]=mag_R_photo
	else:
		f_filter[2]=0
		mag_available[2]=0
	if str('I') in filter_available :
		#si epoque la plus proche est entre min et max phase i
		if (epoque_photo<=max(MJD['I'])) and (epoque_photo>=min(MJD['I'])):
			mag_I_photo=I_band(epoque_photo)
		elif (epoque_photo>max(MJD['I'])): 
			mag_I_photo=I_band(max(MJD['I'])-0.0001)
		elif  (epoque_photo<min(MJD['I'])): 
			mag_I_photo=I_band(min(MJD['I'])+0.0001)
		delta_M_I=mag_I_photo-m_I_kait
		f_I=10**(-0.4*delta_M_I)
		f_filter[3]=f_I
		mag_available[3]=mag_I_photo
	else:
		f_filter[3]=0
		mag_available[3]=0

	#We create a warp function for all the spectrum. We need to interpolate between the filter and extrapolate in the blue
	#and red part of the spectrum.
	x_mang=np.array(lambda_eff_tot_kait)[f_filter!=0]
	y_mang=f_filter[f_filter!=0]
	if np.size(y_mang)>=4:
		cubic_interp=interpolate.interp1d(x_mang, y_mang,kind='cubic')
		fine_x=np.linspace(min(x_mang), max(x_mang), 150)
		y_interp = cubic_interp(fine_x)
	if (np.size(y_mang)<4) & (np.size(y_mang)>2):
		quadratic_interp=interpolate.interp1d(x_mang, y_mang,kind='quadratic')
		fine_x=np.linspace(min(x_mang), max(x_mang), 150)
		y_interp = quadratic_interp(fine_x)
	if np.size(y_mang)==2:
		linear_interp=interpolate.interp1d(x_mang, y_mang,kind='linear')
		fine_x=np.linspace(min(x_mang), max(x_mang), 100)
		y_interp = linear_interp(fine_x)

	if np.size(y_mang)==1:
		fine_x=x_mang
		y_interp = y_mang

	################### We do a flat extrapolation ###############
	x_list=list(fine_x)
	y_list=list(y_interp)
	#We create a fake filter with the same value as the first filter to cover the blue part of the spectrum
	y_list.insert(0,y_list[0]) 
	x_list.insert(0,min(wl_mod_obs))
	#We create a fake filter with the same value as the last filter to cover the red part of the spectrum
	y_list.insert( np.size(y_list),y_list[np.size(y_list)-1]) 
	x_list.insert( np.size(y_list),max(wl_mod_obs))
	#We create the mangling function
	F_mang_func=interpolate.interp1d(np.array(x_list),np.array(y_list))
	coeff_flux=F_mang_func(wl_mod_obs) #wavelength spectre
	fl_mod_mang=abs((coeff_flux)*fl_mod_obs)
	#We interpolate the warp model
	F_mod_mang=interpolate.interp1d(wl_mod_obs,fl_mod_mang) 
	
	'''
	We have the spectrum after mangling, let's check is now the colours match the observed colours.
	'''
	
	############ Mangling verification ################ 
	if str('B') in filter_available :	
		# B band KAIT
		m_B_kait_mang=-2.5*np.log10(1.0/(h_erg*c_AA)*integrate.simps(F_mod_mang(lambda_X_tot[ind_sys][0])*F_X_func_tot[ind_sys][0](lambda_X_tot[ind_sys][0])*lambda_X_tot[ind_sys][0],lambda_X_tot[ind_sys][0]))+ZP_X_kait[0]
	else:
		m_B_kait_mang=0
	# V band KAIT
	if str('V') in filter_available :	
		m_V_kait_mang=-2.5*np.log10(1.0/(h_erg*c_AA)*integrate.simps(F_mod_mang(lambda_X_tot[ind_sys][1])*F_X_func_tot[ind_sys][1](lambda_X_tot[ind_sys][1])*lambda_X_tot[ind_sys][1],lambda_X_tot[ind_sys][1]))+ZP_X_kait[1]
	else:
		m_V_kait_mang=0
	# R band KAIT
	if str('R') in filter_available :	
		m_R_kait_mang=-2.5*np.log10(1.0/(h_erg*c_AA)*integrate.simps(F_mod_mang(lambda_X_tot[ind_sys][2])*F_X_func_tot[ind_sys][2](lambda_X_tot[ind_sys][2])*lambda_X_tot[ind_sys][2],lambda_X_tot[ind_sys][2]))+ZP_X_kait[2]

	else:
		m_R_kait_mang=0
	# I band KAIT
	if str('I') in filter_available :	
		m_I_kait_mang=-2.5*np.log10(1.0/(h_erg*c_AA)*integrate.simps(F_mod_mang(lambda_X_tot[ind_sys][3])*F_X_func_tot[ind_sys][3](lambda_X_tot[ind_sys][3])*lambda_X_tot[ind_sys][3],lambda_X_tot[ind_sys][3]))+ZP_X_kait[3]

	else:
		m_I_kait_mang=0

	m_tot_mang=[m_B_kait_mang,m_V_kait_mang,m_R_kait_mang,m_I_kait_mang]

	#We compare the new magnitude with the true photometry, we need less than 0.05mag
	maxiter=10
	while False in (abs(mag_available-m_tot_mang)<0.05):
		maxiter = maxiter - 1
		if maxiter <= 0:
			print ("Did not converge.")
			break

		# on calcule les magnitudes synthetiques
		if str('B') in filter_available:
			delta_M_B=mag_B_photo-m_B_kait_mang
			f_B=10**(-0.4*delta_M_B)
			f_filter[0]=f_B
		else:
			f_filter[0]=0

		if str('V') in filter_available :
			delta_M_V=mag_V_photo-m_V_kait_mang
			f_V=10**(-0.4*delta_M_V)
			f_filter[1]=f_V
		else:
			f_filter[1]=0
		if str('R') in filter_available :
			delta_M_R=mag_R_photo-m_R_kait_mang
			f_R=10**(-0.4*delta_M_R)
			f_filter[2]=f_R
		else:
			f_filter[2]=0

		if str('I') in filter_available :
			delta_M_I=mag_I_photo-m_I_kait_mang
			f_I=10**(-0.4*delta_M_I)
			f_filter[3]=f_I
		else:
			f_filter[3]=0

		x_mang=np.array(lambda_eff_tot_kait)[f_filter!=0]
		y_mang=f_filter[f_filter!=0]
		if np.size(y_mang)>=4:
			cubic_interp=interpolate.interp1d(x_mang, y_mang,kind='cubic')
			fine_x=np.linspace(min(x_mang), max(x_mang), 150)
			y_interp = cubic_interp(fine_x)
		if ( np.size(y_mang)<4) & ( np.size(y_mang)>2):
			quadratic_interp=interpolate.interp1d(x_mang, y_mang,kind='quadratic')
			fine_x=np.linspace(min(x_mang), max(x_mang), 150)
			y_interp = quadratic_interp(fine_x)
		if np.size(y_mang)==2:
			linear_interp=interpolate.interp1d(x_mang, y_mang,kind='linear')
			fine_x=np.linspace(min(x_mang), max(x_mang), 100)
			y_interp = linear_interp(fine_x)

		if np.size(y_mang)==1:
			fine_x=x_mang
			y_interp = y_mang


		################### On fait un extrapolation plane ###############
		x_list=list(fine_x)
		y_list=list(y_interp)
		#on rajoute la meme valeur que le premier filtre
		y_list.insert(0,y_list[0]) 
		x_list.insert(0,min(wl_mod_obs))
		y_list.insert( np.size(y_list),y_list[np.size(y_list)-1]) # meme valeur que le dernier filtre utilise
		x_list.insert( np.size(y_list),max(wl_mod_obs))
		F_mang_func=interpolate.interp1d(np.array(x_list),np.array(y_list)) # on cree la fonction mangling
		coeff_flux=F_mang_func(wl_mod_obs) #wavelength spectre
		fl_mod_mang=abs((coeff_flux)*(fl_mod_mang))
		F_model_func_mang=interpolate.interp1d(wl_mod_obs,fl_mod_mang) # on cree la fonction spectre

		wl_mod_rest=wl_mod_obs*1.0/(1+z_hel)
		fl_mod_rest=fl_mod_mang*(1+z_hel)

		if str('B') in filter_available :	
			# B band KAIT
			m_B_kait_mang=-2.5*np.log10(1.0/(h_erg*c_AA)*integrate.simps(F_mod_mang(lambda_X_tot[ind_sys][0])*F_X_func_tot[ind_sys][0](lambda_X_tot[ind_sys][0])*lambda_X_tot[ind_sys][0],lambda_X_tot[ind_sys][0]))+ZP_X_kait[0]
		else:
			m_B_kait_mang=0
		# V band KAIT
		if str('V') in filter_available :	
			m_V_kait_mang=-2.5*np.log10(1.0/(h_erg*c_AA)*integrate.simps(F_mod_mang(lambda_X_tot[ind_sys][1])*F_X_func_tot[ind_sys][1](lambda_X_tot[ind_sys][1])*lambda_X_tot[ind_sys][1],lambda_X_tot[ind_sys][1]))+ZP_X_kait[1]
		else:
			m_V_kait_mang=0
		# R band KAIT
		if str('R') in filter_available :	
			m_R_kait_mang=-2.5*np.log10(1.0/(h_erg*c_AA)*integrate.simps(F_mod_mang(lambda_X_tot[ind_sys][2])*F_X_func_tot[ind_sys][2](lambda_X_tot[ind_sys][2])*lambda_X_tot[ind_sys][2],lambda_X_tot[ind_sys][2]))+ZP_X_kait[2]

		else:
			m_R_kait_mang=0
		# I band KAIT
		if str('I') in filter_available :	
			m_I_kait_mang=-2.5*np.log10(1.0/(h_erg*c_AA)*integrate.simps(F_mod_mang(lambda_X_tot[ind_sys][3])*F_X_func_tot[ind_sys][3](lambda_X_tot[ind_sys][3])*lambda_X_tot[ind_sys][3],lambda_X_tot[ind_sys][3]))+ZP_X_kait[3]

		else:
			m_I_kait_mang=0
	
	#we store for each epoch the model after mangling.
	wl_model_mang.append(wl_mod_obs)
	fl_model_mang.append(fl_mod_mang)
	ep_model_mang.append(epoch[ep])

	print('Photo observed')
	print(mag_available)
	print ('Photo after mangling:')
	print(sc.around(m_tot_mang,2))

