def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]
# Loops for all the spectra

# MANGLING
for i in range(n_spectra):

	#Spectra available:
	filter_available_photo=[]
	#filter available at this epoch but we remove orange
	for filt in range (size(bands)):
		if (min(MJD[bands[filt]])<spectra_MJD[i]) and (max(MJD[bands[filt]])>spectra_MJD[i]) and (bands[filt]!=str('orange')): 
			filter_available_photo.append(bands[filt])
	filter_available=[]
	for filt in range (size(filter_available_photo)):
		#From the filter available we look at those which are totally inside the spectrum range
		if (min(wl_spectra[i])<min(dict_ZP[filter_available_photo[filt]][2])) and (max(wl_spectra[i])>max(dict_ZP[filter_available_photo[filt]][2])): 
			filter_available.append(filter_available_photo[filt])

	#If at least there is 1 filter, we do a mangling
	if size(filter_available)!=0:

		#Mag syn for available filter
		F_spec_func=interpolate.interp1d(wl_spectra[i],flux_spectra[i])

		f_filter=np.zeros(size(bands))
		lambda_eff_filter=np.zeros(size(bands))
	
		for filt in range (size(filter_available)):

			#observed magnitudes
			X_band_func=interpolate.interp1d(MJD[filter_available[filt]],mags[filter_available[filt]]) 

			mag_X_photo=X_band_func(spectra_MJD[i])

			#Synthetic magnitude from the spectrum
			dem_X=integrate.simps(dict_ZP[filter_available[filt]][1](dict_ZP[filter_available[filt]][2])*1.0/dict_ZP[filter_available[filt]][2],dict_ZP[filter_available[filt]][2])
			num_X=integrate.simps(F_spec_func(dict_ZP[filter_available[filt]][2])*dict_ZP[filter_available[filt]][1](dict_ZP[filter_available[filt]][2])*dict_ZP[filter_available[filt]][2],dict_ZP[filter_available[filt]][2])   
			m_X_syn=-2.5*(np.log10(num_X/dem_X))-48.60+2.5*(np.log10(c_AA))

			#Diffrence between the two mags
			delta_M_X=mag_X_photo-m_X_syn
			f_X=10**(-0.4*delta_M_X)

			pos_filter=bands.index(filter_available[filt])
			#points of the warping function: coefficient vs lam eff
			f_filter[pos_filter]=f_X
			lambda_eff_filter[pos_filter]=dict_ZP[filter_available[filt]][0]


		#Remove the 0, i.e, no filters	
		x_mang=list(lambda_eff_filter[f_filter!=0])
		y_mang=list(f_filter[f_filter!=0])

		#sort by wavelength
		x_mang,y_mang=zip(*sorted(zip(x_mang,y_mang)))
		x_mang,y_mang=list(x_mang),list(y_mang)

		#First and last SED wavelength will have the same values of the first and the last filter respectively
		# the SED wavelength range is much larger than the filter range, that is why we do that
		y_mang.insert(0,y_mang[0]) 
		x_mang.insert(0,min(wl_spectra[i]))
		y_mang.insert(size(y_mang),y_mang[size(y_mang)-1])
		x_mang.insert(size(y_mang),max(wl_spectra[i]))

		#Spline the warping function
		if size(x_mang)==3:
			tck = interpolate.splrep(x_mang, y_mang, k=2,s=0)
		elif size(x_mang)>3:
			tck = interpolate.splrep(x_mang, y_mang, k=3,s=0)
		elif size(x_mang)==2:
			tck = interpolate.splrep(x_mang, y_mang, k=1,s=0)
		coeff_flux_mang= interpolate.splev(wl_spectra[i], tck, der=0) # all the coeff of each wavelength
		#Mangling the template SED
		flux_spec_mang=(coeff_flux_mang)*(flux_spectra[i])
		F_spec_func_mang=interpolate.interp1d(wl_spectra[i],flux_spec_mang) # on cree la fonction spectre


		# Mangling check. We use the warping function and recalculate the magnitudes
		f_filter=np.zeros(size(bands))
		lambda_eff_filter=np.zeros(size(bands))
		mag_obs=np.zeros(size(bands))
		mag_syn_tot=np.zeros(size(bands))
		for filt in range (size(filter_available)):
			#obs
			X_band_func=interpolate.interp1d(MJD[filter_available[filt]],mags[filter_available[filt]]) #lineal interpo magnitude
			mag_X_photo=X_band_func(spectra_MJD[i])
			#syn
			dem_X=integrate.simps(dict_ZP[filter_available[filt]][1](dict_ZP[filter_available[filt]][2])*1.0/dict_ZP[filter_available[filt]][2],dict_ZP[filter_available[filt]][2])
			num_X=integrate.simps(F_spec_func_mang(dict_ZP[filter_available[filt]][2])*dict_ZP[filter_available[filt]][1](dict_ZP[filter_available[filt]][2])*dict_ZP[filter_available[filt]][2],dict_ZP[filter_available[filt]][2])   
			m_X_syn_mang=-2.5*(np.log10(num_X/dem_X))-48.60+2.5*(np.log10(c_AA))

			#diff
			delta_M_X=mag_X_photo-m_X_syn_mang
			f_X=10**(-0.4*delta_M_X)

			pos_filter=bands.index(filter_available[filt])
			f_filter[pos_filter]=f_X
			lambda_eff_filter[pos_filter]=dict_ZP[filter_available[filt]][0]

			mag_obs[pos_filter]=mag_X_photo
			mag_syn_tot[pos_filter]=m_X_syn_mang

		print(mag_syn_tot-mag_obs)
		print(spectra_name[i])

		# Spectrum after mangling
		fig, ax1 = plt.subplots()
		wavelength_spec_mang=wl_spectra[i]
		ax1.plot(wl_spectra[i],flux_spec_mang,'r',alpha=0.7,label='mang')
		ax1.plot(wl_spectra[i],flux_spectra[i],'k',alpha=0.7,label='old')

		for filt in range(size(filter_available)):
			norm_filt=max(flux_spectra[i])/max((dict_ZP[filter_available[filt]][1](dict_ZP[filter_available[filt]][2])))
			tdataplot=plot(dict_ZP[filter_available[filt]][2],(dict_ZP[filter_available[filt]][1](dict_ZP[filter_available[filt]][2]))*norm_filt,linestyle='-')
			ax1.plot(dict_ZP[filter_available[filt]][2],(dict_ZP[filter_available[filt]][1](dict_ZP[filter_available[filt]][2]))*norm_filt,color=tdataplot[0].get_color(),linestyle='-',alpha=0.5,label=filter_available[filt])


		ax1.set_xlabel('Flux')
		ax1.legend()

		ax2 = ax1.twinx()
		ax2.plot(wl_spectra[i], coeff_flux_mang,'b',alpha=0.5)
		ax2.plot(x_mang,y_mang, 'sb')
		ax2.set_ylabel('Warping function', color='b')
		ax2.tick_params('y', colors='b')

		fig.tight_layout()

		fig.savefig('Figures/%s_mang.png'%spectra_name[i])
		close()
		matrice=np.zeros((size(flux_spec_mang),3))
		matrice[:,0]=wavelength_spec_mang 
		matrice[:,1]=flux_spec_mang
		matrice[:,2]=err_flux_spectra[i]
		savetxt('Spectra_mang/%s_mang.dat' %(spectra_name[i]),matrice ,fmt='%g')
	else:
		# Spectrum after mangling
		plot(wl_spectra[i],flux_spectra[i],'r',label='mang')
		plot(wl_spectra[i],flux_spectra[i],'k',label='old')
		legend()
		savefig('Figures/%s_mang.png'%spectra_name[i])
		close()
		matrice=np.zeros((size(flux_spectra[i]),3))
		matrice[:,0]=wl_spectra[i]
		matrice[:,1]=flux_spectra[i]
		matrice[:,2]=err_flux_spectra[i]
		savetxt('Spectra_mang/%s_mang.dat' %(spectra_name[i]),matrice ,fmt='%g')
