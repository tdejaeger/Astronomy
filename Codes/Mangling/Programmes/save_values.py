#Search for MW extinction using coordinates
from urllib.request import urlopen
from xml.dom.minidom import parse
BASE_URL = "http://irsa.ipac.caltech.edu/cgi-bin/DUST/nph-dust?locstr=%.5f+%.5f"
u =urlopen(BASE_URL % (ra,decl))
dom = parse(u)
u.close()
tag = {'SF11':'meanValueSandF','SFD98':'meanValueSFD'}
EBVstr = dom.getElementsByTagName(tag['SF11'])[0].childNodes[0].data
EBV_sn = float(EBVstr.strip().split()[0])
AvG_sn=EBV_sn*3.1

#Save Tmax, a_tot, pca1 pca2
sndata=open('Results/Values/%s.txt'%sn_name,"w")
sndata.write('#SN zhel z_cmb EBV_MW ra dec Survey\n')
sndata.write('%s %s %s %s %s %s %s\n'%(sn_name,z_sn,z_sn_cmb,EBV_sn,ra,decl,name_survey))

if emcee_fit==str('True'):
	sndata.write('######### EMCEE values ############\n')
	sndata.write('#Tmax B MJD:\n%s +%s -%s\n%s\n'%(Tmax,err_Tmax_plus,err_Tmax_moins,result_mcmc[size(bands)]))

	for i in range(size(a_tot)):
		sndata.write('#factor %s:\n%s +%s -%s\n'%(bands[i],a_tot[i],err_a_tot_plus[i],err_a_tot_moins[i]))

	if alpha1:
		sndata.write('#factor pca1:\n%s +%s -%s\n'%(alpha1,err_alpha1_plus,err_alpha1_moins))
	if beta1:
		sndata.write('#factor pca2:\n%s +%s -%s\n'%(beta1,err_beta1_plus,err_beta1_moins))

else:
	sndata.write('######### Likelihood ############ \n')
	sndata.write('#Tmax B MJD:\n %s\n %s\n'%Tmax,Tmax-min(epoch_tot))

	for i in range(size(a_tot)):
		sndata.write('#factor %s:\n%s\n'%(bands[i],a_tot[i]))

	if alpha1:
		sndata.write('#factor pca1:\n%s\n'%(alpha1))
	if beta1:
		sndata.write('#factor pca2:\n%s\n'%(beta1))

sndata.close()
