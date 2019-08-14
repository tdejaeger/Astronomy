import re # use regular patters
import sys # system commands
import string as string # string functions4
import math
import numpy as np # numerical tools
from scipy import *
from pylab import *
import glob 

c_light=299792.458
#SNname	zhel	ra(deg)	dec(deg)
name_SN=str('sn2016hnk')
zhel_SN=0.016
ra_SN=33.31929
dec_SN=-7.66133

'''
Module for computing the CMB redshift from the heliocentric redshift and
coordinates. For now, uses NED calculatro, but maybe later will do all the 
math, so we don't need the interweb'''

from urllib.request import urlopen
from xml.dom.minidom import parse
url_temp = "http://ned.ipac.caltech.edu/cgi-bin/velc?lon=%.6fd&in_csys=Equatorial&lat=%.6fd&in_equinox=J2000.0&vel=%.3f&vfrom=Heliocentric&vto=3K&alon=&a_csys=Equatorial&alat=&a_equinox=J2000.0&avel=0.0"
pat = re.compile(r'([0-9]+(\.[0-9]*)?)')


url_temp="https://ned.ipac.caltech.edu/cgi-bin/velc?lon=%0.6fd&in_csys=Equatorial&lat=%.6f&in_equinox=J2000.0&vel=%.3f&vfrom=Heliocentric&vto=3K&alon=&a_csys=Equatorial&alat=&a_equinox=J2000.0&avel=0.0"

def z_cmb(z_hel, ra, dec):
	v_hel = z_hel*c_light
	u = urlopen(url_temp % (ra, dec, v_hel))
	lines = u.readlines()
	vcmb = None
	for line in lines:
		if str(line).find('Output:') >= 0:
			vcmb = float(line.split()[1][9:30])
	return round(vcmb/c_light,7)

z_cmb_SN=z_cmb(zhel_SN,ra_SN,dec_SN)


#### Load the photometry ###
#JD orange cyan u g r i z B V
MJD_SN_B=[]
B_SN=[]
err_B_SN=[]
MJD_SN_V=[]
V_SN=[]
err_V_SN=[]
MJD_SN_u=[]
u_SN=[]
err_u_SN=[]
MJD_SN_g=[]
g_SN=[]
err_g_SN=[]
MJD_SN_r=[]
r_SN=[]
err_r_SN=[]
MJD_SN_i=[]
i_SN=[]
err_i_SN=[]
MJD_SN_z=[]
z_SN=[]
err_z_SN=[]

MJD_SN_cyan=[]
cyan_SN=[]
err_cyan_SN=[]
MJD_SN_orange=[]
orange_SN=[]
err_orange_SN=[]

data_B=open('Photometry/%s.out_B'%(name_SN),'r')
lines_B=data_B.readlines()
data_V=open('Photometry/%s.out_V'%(name_SN),'r')
lines_V=data_V.readlines()
data_u=open('Photometry/%s.out_u_2.5m'%(name_SN),'r')
lines_u=data_u.readlines()
data_g=open('Photometry/%s.out_g_2.5m'%(name_SN),'r')
lines_g=data_g.readlines()
data_r=open('Photometry/%s.out_r_2.5m'%(name_SN),'r')
lines_r=data_r.readlines()
data_i=open('Photometry/%s.out_i_2.5m'%(name_SN),'r')
lines_i=data_i.readlines()
data_z=open('Photometry/%s.out_z_2.5m'%(name_SN),'r')
lines_z=data_z.readlines()
data_c=open('Photometry/%s.out_c'%(name_SN),'r')
lines_c=data_c.readlines()
data_o=open('Photometry/%s.out_o'%(name_SN),'r')
lines_o=data_o.readlines()

for j in range(size(lines_u)):
	if (lines_u[j].split()[0][0])==str('5'):
		epoch=float(lines_u[j].split()[0])	
		mag=float(lines_u[j].split()[1])
		err_mag=float(lines_u[j].split()[2])

		MJD_SN_u.append(epoch)
		u_SN.append(mag)
		err_u_SN.append(err_mag)

for j in range(size(lines_g)):
	if (lines_g[j].split()[0][0])==str('5'):
		epoch=float(lines_g[j].split()[0])	
		mag=float(lines_g[j].split()[1])
		err_mag=float(lines_g[j].split()[2])

		MJD_SN_g.append(epoch)
		g_SN.append(mag)
		err_g_SN.append(err_mag)
for j in range(size(lines_r)):
	if (lines_r[j].split()[0][0])==str('5'):
		epoch=float(lines_r[j].split()[0])	
		mag=float(lines_r[j].split()[1])
		err_mag=float(lines_r[j].split()[2])

		MJD_SN_r.append(epoch)
		r_SN.append(mag)
		err_r_SN.append(err_mag)
for j in range(size(lines_i)):
	if (lines_i[j].split()[0][0])==str('5'):
		epoch=float(lines_i[j].split()[0])	
		mag=float(lines_i[j].split()[1])
		err_mag=float(lines_i[j].split()[2])

		MJD_SN_i.append(epoch)
		i_SN.append(mag)
		err_i_SN.append(err_mag)
for j in range(size(lines_z)):
	if (lines_z[j].split()[0][0])==str('5'):
		epoch=float(lines_z[j].split()[0])	
		mag=float(lines_z[j].split()[1])
		err_mag=float(lines_z[j].split()[2])

		MJD_SN_z.append(epoch)
		z_SN.append(mag)
		err_z_SN.append(err_mag)
for j in range(size(lines_B)):
	if (lines_B[j].split()[0][0])==str('5'):
		epoch=float(lines_B[j].split()[0])	
		mag=float(lines_B[j].split()[1])
		err_mag=float(lines_B[j].split()[2])

		MJD_SN_B.append(epoch)
		B_SN.append(mag)
		err_B_SN.append(err_mag)
for j in range(size(lines_V)):
	if (lines_V[j].split()[0][0])==str('5'):
		epoch=float(lines_V[j].split()[0])	
		mag=float(lines_V[j].split()[1])
		err_mag=float(lines_V[j].split()[2])

		MJD_SN_V.append(epoch)
		V_SN.append(mag)
		err_V_SN.append(err_mag)
for j in range(size(lines_c)):
	if (lines_c[j].split()[0][0])==str('5'):
		epoch=float(lines_c[j].split()[0])	
		mag=float(lines_c[j].split()[1])
		err_mag=float(lines_c[j].split()[2])

		MJD_SN_cyan.append(round(epoch,1))
		cyan_SN.append(mag)
		err_cyan_SN.append(err_mag)
for j in range(size(lines_o)):
	if (lines_o[j].split()[0][0])==str('5'):
		epoch=float(lines_o[j].split()[0])	
		mag=float(lines_o[j].split()[1])
		err_mag=float(lines_o[j].split()[2])

		MJD_SN_orange.append(round(epoch,1))
		orange_SN.append(mag)
		err_orange_SN.append(err_mag)

#Average orange/cyan filters

epoch_com_orange=list(set([x for x in MJD_SN_orange if MJD_SN_orange.count(x) > 1]))
epoch_com_orange.sort()
MJD_final_orange=[]
mag_final_orange=[]
err_mag_final_orange=[]
for i in range(size(epoch_com_orange)):

	MJD_final_orange.append(epoch_com_orange[i])
	mag_final_orange.append(mean(np.array(orange_SN)[where(np.array(MJD_SN_orange)==epoch_com_orange[i])[0]]))
	err_mag_final_orange.append(mean(np.array(err_orange_SN)[where(np.array(MJD_SN_orange)==epoch_com_orange[i])[0]]))


epoch_com_cyan=list(set([x for x in MJD_SN_cyan if MJD_SN_cyan.count(x) > 1]))
epoch_com_cyan.sort()
MJD_final_cyan=[]
mag_final_cyan=[]
err_mag_final_cyan=[]
for i in range(size(epoch_com_cyan)):

	MJD_final_cyan.append(epoch_com_cyan[i])
	mag_final_cyan.append(mean(np.array(cyan_SN)[where(np.array(MJD_SN_cyan)==epoch_com_cyan[i])[0]]))
	err_mag_final_cyan.append(mean(np.array(err_cyan_SN)[where(np.array(MJD_SN_cyan)==epoch_com_cyan[i])[0]]))


# Same in Snoopy format
#SN1981D 0.005871 50.65992 -37.23272
#filter Bs
#674.8593   12.94   0.11

filtre_sn=['u_2.5m','g_2.5m','r_2.5m','i_2.5m','z_2.5m','B','V','cyan','orange']
photo_tot=[u_SN,g_SN,r_SN,i_SN,z_SN,B_SN,V_SN,mag_final_cyan,mag_final_orange]
err_photo_tot=[err_u_SN,err_g_SN,err_r_SN,err_i_SN,err_z_SN,err_B_SN,err_V_SN,err_mag_final_cyan,err_mag_final_orange]
MJD_tot=[MJD_SN_u,MJD_SN_g,MJD_SN_r,MJD_SN_i,MJD_SN_z,MJD_SN_B,MJD_SN_V,MJD_final_cyan,MJD_final_orange]

sndata=open('%s.txt'%name_SN,"w")
sndata.write('%s %s %s %s %s\n'%(name_SN,zhel_SN,z_cmb_SN,ra_SN,dec_SN))
	
for fil in range(size(filtre_sn)):		
	if size(np.array(photo_tot[fil])[np.array(photo_tot[fil])!=0])>1:
		sndata.write('filter %s\n'%filtre_sn[fil])
		for j in range(size(MJD_tot[fil])):
			if photo_tot[fil][j]!=0.0:
				sndata.write('%s\t%s\t%s\n'%(round(MJD_tot[fil][j],1),photo_tot[fil][j],err_photo_tot[fil][j]))
sndata.close()


