import numpy as np
import matplotlib.pyplot as plt


#PArams figures ###
plt.rcParams['legend.numpoints']=1
plt.rcParams['xtick.major.size'] = 11
plt.rcParams['xtick.minor.size'] = 5
plt.rcParams['ytick.major.size'] = 11
plt.rcParams['ytick.minor.size'] = 5
plt.rcParams['xtick.minor.visible']=True #See minor tick
plt.rcParams['text.usetex']=True #use Latex
plt.rcParams['axes.linewidth']=2 #width axes
plt.rcParams['axes.labelsize']=21 #
plt.rcParams['axes.labelpad']=10 #
plt.rcParams['ytick.labelsize']=15 #fontsize of tick labels
plt.rcParams['xtick.labelsize']=15 #fontsize of tick labels
plt.rcParams['ytick.direction']='inout' ## direction: in, out, or inout
plt.rcParams['xtick.direction']='inout' ## direction: in, out, or inout
plt.rcParams['xtick.major.top']=True #draw x axis top major ticks
plt.rcParams['xtick.major.bottom']=True #draw x axis bottom major ticks
plt.rcParams['xtick.minor.top']=True ## draw x axis top minor ticks
plt.rcParams['xtick.minor.bottom']=True #draw x axis bottom minor ticks
plt.rcParams['ytick.major.left']=True #draw y axis left major ticks
plt.rcParams['ytick.major.right']=True #draw y axis right major ticks
plt.rcParams['ytick.minor.left']=True ## draw y axis left minor ticks
plt.rcParams['ytick.minor.right']=True #draw y axis right minor ticks
plt.rcParams['font.weight']='bold'
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titlesize']=22
plt.rcParams['figure.titlesize']=22

plt.rcParams['text.latex.preamble']=[r'\boldmath']


#From Riess et al. 16
Ab_R16 = 0.71273
err_Ab_R16 = 0.00176

sn_host_R16 = ["744", "1241", "1371", "1794", "2102", "2635","2916", "2992", "5549", "6057", "6558", "7147","8921", "12779", "12781", "12898", "12950","13005", "13038", "13044", "13736", "13796","14024", "14108", "14871", "15234", "16021","16619", "16641", "17240", "17366", "17497","18241", "18298", "18602", "18697", "18809","18835", "18855", "19775", "19899", "19953","20084", "20470", "21502", "722", "774", "1032","2308", "2561", "3592", "5395", "5751", "6406","8719", "10028", "10434", "10805", "11067","11300", "12860", "12928", "12930", "13135","13894", "14437", "15171", "15508", "16069","16185", "16259", "17208", "17258", "17280","17605", "17629", "17745", "18415", "18612","18650", "19543", "19940", "19968", "20064","20625", "21034", "21510", "22075", "06D2fb","2004ef", "2005eq", "2005hc", "2005hf", "2005hj","2005iq", "2005lz", "2005mc", "2005ms", "2005na","2006ac", "2006ak", "2006al", "2006bb", "2006bu","2006cj", "2006cq", "2006cz", "2006en", "2006ev","2006gr", "2006je", "2006mo", "2006oa", "2006on","2006qo", "2006sr", "2006te", "2007ae", "2007ai","2007ba", "2007bd", "2007co", "2007cq", "2007F","2007O", "2008af", "2008bf", "2010dt", "1997dg","1998dx", "1998eg", "1999cc", "1999ef", "1999X","2000cf", "2000cn", "2001ah", "2001az", "2002bf","2002bz", "2002ck", "2002de", "2002G", "2002hd","2003cq", "2003it", "2004as", "2004L", "2002he","2003fa", "2003ic", "2003iv", "2003U", "2007aj","2007cb", "2007cc", "2007is", "2007kh", "2007kk","2007nq", "2007ob", "2007su", "2007ux", "2008Y","2008ac", "2008ar", "2008at", "2008bw", "2008by","2008bz", "2008cf", "2008fr", "2008gb", "2008gl","2009D", "2009ad", "2008050", "2008051", "2004gu","2005ag", "2005be", "2005ir", "2006eq", "2006lu","2006ob", "2006py", "2007hx", "2008bq", "1993ac","1994M", "1994Q", "1996ab", "1996bl", "1996C","1990af", "1990O", "1990T", "1990Y", "1991S","1991U", "1992ae", "1992aq", "1992bg", "1992bh","1992bk", "1992bl", "1992bp", "1992bs", "1992J","1992P", "1993ag", "1993H", "1993O", "2000bh","2000bk", "2000ca", "2001ba"]

n_SN_R16=np.size(sn_host_R16)

#first let's download the data from Scolnic et al. 2018
#Remove all the duplicate. The majority have the name however, a handful of SDSS that have duplicates (sdss=['16314','16392','16333','14318','17186','17784','7876’],truname=['2006oa','2006ob','2006on','2006py','2007hx','2007jg','2005ir’]).
sdss_dup=['16314','16392','16333','14318','17186','17784','7876']

name_SN=[]
zHD=[]
zHD_err=[]
z =[]
z_err =[]
Host_Logmass=[]
Host_Logmass_err=[]
Peak_MJD=[]
Peak_MJD_err=[]
x1=[]
x1_err=[] 
c =[]
c_err=[] 
mB=[]
mB_err=[]
x0 =[]
x0_err=[]
n_skip=15 # skip the first 15 lines
f=open('supercal_vH0.fitres.txt')
for i, l in enumerate(f):
	vals=l.split()
	if (i > n_skip - 1):
		if (vals[1] not in sdss_dup) and (vals[1] not in name_SN) and (float(vals[6]) <0.4) and (float(vals[6]) >0.01):
			name_SN.append(vals[1])
			zHD.append(float(vals[6]))
			zHD_err.append(float(vals[7]))
			z.append(float(vals[8])) 
			z_err.append(float(vals[9]))
			Host_Logmass.append(float(vals[10]))
			Host_Logmass_err.append(float(vals[11]))
			Peak_MJD.append(float(vals[15]))
			Peak_MJD_err.append(float(vals[16]))
			x1.append(float(vals[17]))
			x1_err.append(float(vals[18]))
			c.append(float(vals[19])) 
			c_err.append(float(vals[20]))
			
			if float(vals[10]) > 10:
				mB.append(float(vals[21])+0.03)
			else:
				mB.append(float(vals[21])-0.03)
			mB_err.append(float(vals[22])) 
			x0.append(float(vals[23]))
			x0_err.append(float(vals[24]))

# Now we make some cuts following R16.
z_min=0.0233
z_max=0.15

name_sne_cut=[]
zHD_cut=[]
zHD_err_cut=[]
mB_cut=[]
mB_err_cut=[]
c_cut=[]
c_err_cut=[]
x1_cut=[]
x1_err_cut=[]
#Riess et al.
for i in range(np.size(zHD)):

	if (zHD[i]>z_min and zHD[i]<z_max and name_SN[i] in sn_host_R16):

		name_sne_cut.append(name_SN[i])
		zHD_cut.append(zHD[i])
		zHD_err_cut.append(zHD_err[i])
		mB_cut.append(mB[i])
		mB_err_cut.append(mB_err[i])
		c_cut.append(c[i])
		c_err_cut.append(c_err[i])
		x1_cut.append(x1[i])
		x1_err_cut.append(x1_err[i])

zHD_cut=np.array(zHD_cut)
zHD_err_cut=np.array(zHD_err_cut)
mB_cut=np.array(mB_cut)
mB_err_cut=np.array(mB_err_cut)
c_cut=np.array(c_cut)
c_err_cut=np.array(c_err_cut)
x1_cut=np.array(x1_cut)
x1_err_cut=np.array(x1_err_cut)


c_light=299792.458
q0 = -0.5575 # Betoule et al. 2014
q0_err =  0.0510 # Betoule et al. 2014
j0=1

#
def SN_mag_z(z,q0,j0):
	return c_light*z*(1.0+0.5*(1-q0)*z-(1/6)*(1-q0-3*q0**2+j0)*z**2)

from scipy.odr import *
# Orthogonal distance regression
# Define a function to fit the data with.

#Formula (5) in Riess et al. 16
#ax=log10(c_light*z*(1.0+0.5*(1-q0)*z-(1/6)*(1-q0-3*q0**2+j0)*z**2)
def linear_func(p, x):
    b= p
    return np.array(x) - b

#All the SNe Ia
x_plot=np.log10(SN_mag_z(np.array(zHD),q0,j0))
y_plot=0.2*np.array(mB)
err_y_plot=0.2*np.array(mB_err)

#Only those used in Riess et al. 16 (~217 SNe Ia)
x_plot_cut=np.log10(SN_mag_z(zHD_cut,q0,j0))
y_plot_cut=0.2*mB_cut
err_y_plot_cut=0.2*mB_err_cut

#Intrinsic dispersion for Ia
sig_sys_Ia=0.1

linear_model = Model(linear_func) # Create a model for fitting
data_fit = RealData(x_plot_cut,y_plot_cut,sy=np.sqrt((err_y_plot_cut)**2+(0.2*sig_sys_Ia)**2))
odr_data = ODR(data_fit, linear_model, beta0=[0.5])
odr_data.set_job(fit_type=2) #fit_type=0 full ODR, 2 least squares optimisation

res_odr = odr_data.run() # Run the regression.
print('Using ODR' )
print(f'aB = {res_odr.beta[0]:.6f} +/- {res_odr.sd_beta[0]:.6f}')

plot_fit=input('plot figure yes or no: ')

if plot_fit=='yes':

	fig = plt.figure(figsize=(8, 6))
	ax0 = plt.subplot2grid((3, 3), (0, 0),colspan=3,rowspan=2)	

	ax0.errorbar(x_plot,np.array(y_plot),yerr=np.array(err_y_plot),marker='o',color='k',markersize=3,linestyle='None')
	ax0.errorbar(x_plot_cut,y_plot_cut,yerr=err_y_plot_cut,marker='s',color='r',markersize=3,linestyle='None')
	ax0.plot(x_plot,np.array(x_plot)-res_odr.beta[0],'r--')
	ax0.set_yticks([2.5,3,3.5,4.0,4.5])

	ax1 = plt.subplot2grid((3, 1), (2, 0),colspan=3,sharex=ax0)
	ax1.plot(x_plot,np.array(y_plot)-(np.array(x_plot)-res_odr.beta[0]),'ko')
	ax1.plot(x_plot,np.array(y_plot)-np.array(y_plot),'r--')
	ax1.set_ylim([-0.2,0.2])
	ax1.set_xticks([3.5,4.0,4.5,5.0])
	fig.subplots_adjust(hspace=0)   
	for ax in [ax0]:
	    plt.setp(ax.get_xticklabels(), visible=False)
	    # The y-ticks will overlap with "hspace=0", so we'll hide the bottom tick
	    ax.set_yticks(ax.get_yticks()[1:])

	ax1.set_xlabel(r'\textbf{log(cz[1+0.5(1-q$_{0}$)z-(1/6)(1-q$_{0}$-3q$_{0}^{2}$)z$^{2}$])}')
	ax1.set_ylabel(r'\textbf{$\Delta$ 0.2m$_{B}$]')
	ax0.set_ylabel(r'\textbf{0.2m$_{B}$ [mag]}')
	plt.show()


#Method using matrix 
n_s = len(zHD_cut)
n_par = 1
y_vec = np.zeros(n_s)
l_mat = np.zeros((n_s, n_par))
c_mat_inv = np.zeros((n_s, n_s))
# loop through SNe
for i in range(0, n_s):
	y_vec[i] = np.log10(SN_mag_z(zHD_cut[i],q0,j0)) - 0.2 * mB_cut[i]
	l_mat[i, 0] = 1.0
	c_mat_inv[i, i] = 1.0 / (0.2 ** 2 *(mB_err_cut[i] ** 2 + sig_sys_Ia ** 2))
	#k += 1

# fit, calculate residuals in useable form and return
ltci = np.dot(l_mat.transpose(), c_mat_inv)
coVar_aB_mat = np.linalg.inv(np.dot(ltci, l_mat))
aB_mat = np.dot(np.dot(coVar_aB_mat, ltci), y_vec)
err_aB_mat=np.sqrt(np.diag(coVar_aB_mat))
residual = y_vec - np.dot(l_mat, aB_mat)

print('Using Matrix')
print(f'aB = {aB_mat[0]:.6f} +/- {err_aB_mat[0]:.6f}')
print('Value found by Riess et al. 2016')
print(f'aB = {Ab_R16:.6f} +/- {err_Ab_R16:.6f}')

