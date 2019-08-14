##List of spectra ###
import glob
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


spectra_files=glob.glob("Spectra/*.txt")
n_spectra=size(spectra_files)
spectra_files.sort()
date_spectra=[]
spectra_MJD=[]
wl_spectra=[]
flux_spectra=[]
err_flux_spectra=[]
spectra_name=[]
for i in range(n_spectra):
	data_spec=np.loadtxt('%s'%(spectra_files[i])).transpose()
	date_spectra.append(int(spectra_files[i].split('/')[1].split('_')[1].split('.')[0]))
	wl_spectra.append(data_spec[0])
	flux_spectra.append(data_spec[1])
	err_flux_spectra.append(data_spec[2])
	spectra_MJD.append(date_to_jd(int(str(date_spectra[i])[0:4]),int(str(date_spectra[i])[4:6]),int(str(date_spectra[i])[6:8]))- 2400000.5)
	spectra_name.append(spectra_files[i].split('/')[1].split('.')[0])

