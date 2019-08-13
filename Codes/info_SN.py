##############################################################
#
#	Get SN information from TNS and NED
#
#		Input: SN name, e.g SN1999em
#			
#		Output: Galaxy, ra,dec,AvG,z hel, zcmb
#
#
##############################################################

#convert hour min sec to degree
def HMS2deg(ra='', dec=''):
	RA, DEC, rs, ds = '', '', 1, 1
	if dec:
		D, M, S = [float(i) for i in dec.split()]
		if str(D)[0] == '-':
			ds, D = -1, abs(D)
	deg = D + (M/60) + (S/3600)
	DEC = '{0}'.format(deg*ds)
  
	if ra:
		H, M, S = [float(i) for i in ra.split()]
		if str(H)[0] == '-':
			rs, H = -1, abs(H)
		deg = (H*15) + (M/4) + (S/240)
		RA = '{0}'.format(deg*rs)

	if ra and dec:
		return (RA, DEC)
	else:
		return RA or DEC


import numpy as np
import re # use regular patters
from urllib.request import urlopen
from xml.dom.minidom import parse

c_light=c_light=299792.458#in km/s
def SN_infos(name):
	url_TNS="https://wis-tns.weizmann.ac.il/object/%s"

	try:
		page=urlopen(url_TNS %(name[2:10]))
	except:
		print('Unknown SN/format, should be SNXXXX')

	url_TNS="https://wis-tns.weizmann.ac.il/object/%s"
	page=urlopen(url_TNS %(name[2:10]))
	lines=page.readlines()
	ind_lines=np.zeros(1)
	ind_lines_z=np.zeros(1)
	ind_lines_ra_dec=np.zeros(1)

	for line in range(np.size(lines)):
		if str(lines[line]).find('Host Name')!=-1:
			ind_lines[0]=line
		if str(lines[line]).find('="cell-redshift">0')!=-1:
			ind_lines_z[0]=line
		if str(lines[line]).find('RA/DEC')!=-1:
			ind_lines_ra_dec[0]=line

	#RA & DEC in degrees
	
	ra=str(lines[int(ind_lines_ra_dec[0])]).split('ra=')[1].split()[0].split('&')[0]
	dec=str(lines[int(ind_lines_ra_dec[0])]).split('dec=')[1].split()[0].split('&')[0]
	ra_SN=float(ra)
	dec_SN=float(dec)

	
	if int(ind_lines[0])!=0:	

		lines_cut=lines[int(ind_lines[0])].split(b'<b>')
		name_host=lines_cut[6].split(b'Host Name</span><div class="value"><b>')[0].split()[0]+lines_cut[6].split(b'Host Name</span><div class="value"><b>')[0].split()[1]
		name_host=name_host.decode("utf-8")

		url_zhel='https://ned.ipac.caltech.edu/cgi-bin/objsearch?objname=%s&extend=no&hconst=73&omegam=0.27&omegav=0.73&corr_z=1&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=RA+or+Longitude&of=pre_text&zv_breaker=30000.0&list_limit=5&img_stamp=YES'%name_host

		u_zhel = urlopen(url_zhel)
		lines = u_zhel.readlines()

		ind_lines=np.zeros(1)
		for line in range(np.size(lines)):
			if str(lines[line]).find('Velocity/Redshift')!=-1:
				ind_lines[0]=line
				
		z_hel=float(lines[int(ind_lines[0])+2].split()[9])			

	elif (int(ind_lines[0])==0) and (int(ind_lines_z[0])!=0): 
		lines_cut=lines[int(ind_lines_z[0])].split(b'<b>')
		z_hel=float(lines[int(ind_lines_z[0])].split(b'cell-redshift">')[1].split(b'</td><td')[0])
		name_host='None'
	
	elif (int(ind_lines[0])==0) and (int(ind_lines_z[0])==0):

		print(' No host galaxy name')
		print(' No redshift information')
		z_hel=999.9
		name_host='None'


	#AvG
	if name_host!='none':

		url_SN="https://ned.ipac.caltech.edu/cgi-bin/objsearch?objname=%s&extend=no&hconst=70&omegam=0.30&omegav=0.70&corr_z=1&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=RA+or+Longitude&of=pre_text&zv_breaker=30000.0&list_limit=5&img_stamp=YES"
		ind_AvG=np.zeros(1)
		u=urlopen(url_SN %(name_host))
		lines=u.readlines()
		for line in range(np.size(lines)):
			if str(lines[line]).find('>V</td><td>(0.54)')!=-1:
				ind_AvG[0]=line

	else:
		url_SN="https://ned.ipac.caltech.edu/cgi-bin/objsearch?objname=%s&extend=no&hconst=70&omegam=0.30&omegav=0.70&corr_z=1&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=RA+or+Longitude&of=pre_text&zv_breaker=30000.0&list_limit=5&img_stamp=YES"
		ind_AvG=np.zeros(1)
		u=urlopen(url_SN %(name))
		lines=u.readlines()
		for line in range(np.size(lines)):
			if str(lines[line]).find('>V</td><td>(0.54)')!=-1:
				ind_AvG[0]=line
	if int(ind_AvG[0])!=0:
		AvG_SN=float(str(lines[int(ind_AvG[0])].split()[1])[2:7])

	if (z_hel!=999.9):

		# Convert z hel to z cmb
		url_zcmb = "http://ned.ipac.caltech.edu/cgi-bin/velc?lon=%.6fd&in_csys=Equatorial&lat=%.6fd&in_equinox=J2000.0&vel=%.3f&vfrom=Heliocentric&vto=3K&alon=&a_csys=Equatorial&alat=&a_equinox=J2000.0&avel=0.0"
		pat = re.compile(r'([0-9]+(\.[0-9]*)?)')
	
		v_hel = z_hel*c_light
		u = urlopen(url_zcmb % (ra_SN, dec_SN, v_hel))
		lines = u.readlines()
		vcmb = None
		for line in lines:
			if str(line).find('Output:') >= 0:
				vcmb = float(line.split()[1][9:30])
		z_cmb=round(vcmb/c_light,5)

	
		if (int(ind_AvG[0])!=0):
			print('Galaxy ra dec AvG zhel zcmb')
			return name_host,float("{:.4f}".format(ra_SN)),float("{:.4f}".format(dec_SN)),float("{:.3f}".format(AvG_SN)),float("{:.6f}".format(z_hel)),float("{:.6f}".format(z_cmb))


		elif (int(ind_AvG[0])==0):

			print('Galaxy ra dec AvG zhel zcmb')
			return name_host,float("{:.4f}".format(ra_SN)),float("{:.4f}".format(dec_SN)),'None',float("{:.6f}".format(z_hel)),float("{:.6f}".format(z_cmb))


	elif (z_hel==999.9):

		if (int(ind_AvG[0])!=0):
			print('Galaxy ra dec AvG zhel zcmb')
			return name_host,float("{:.4f}".format(ra_SN)),float("{:.4f}".format(dec_SN)),float("{:.3f}".format(AvG_SN))

		elif (int(ind_AvG[0])==0):
			print('Galaxy ra dec')
			return name_host,float("{:.3f}".format(ra_SN)),float("{:.3f}".format(dec_SN))
