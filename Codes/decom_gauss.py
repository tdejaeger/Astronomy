#Module: Decomposition of H_alpha line with 3 Gaussians
#              A non linear least-squares script solves
#               the best fit multi Gaussian decomposition
#                
#Author: de JAEGER Thomas
###########################################################################


import numpy as np
from scipy.optimize import leastsq
import matplotlib.pyplot as plt

plt.rcParams['legend.numpoints']=1
plt.rcParams['xtick.major.size'] = 11
plt.rcParams['xtick.minor.size'] = 5
plt.rcParams['ytick.major.size'] = 11
plt.rcParams['ytick.minor.size'] = 5
plt.rcParams['xtick.minor.visible']=True #See minor tick
plt.rcParams['text.usetex']=True #use Latex
plt.rcParams['axes.linewidth']=2 #width axes
plt.rcParams['axes.labelsize']=25 #
plt.rcParams['ytick.labelsize']=25 #fontsize of tick labels
plt.rcParams['xtick.labelsize']=25 #fontsize of tick labels
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

plt.rcParams['text.latex.preamble']=[r'\boldmath']

text_font = {'size':'27', 'color':'black', 'weight':'heavy',
              'verticalalignment':'center','horizontalalignment':'center'}
axis_font = {'size':'30','weight':'heavy'}

######################################
# Setting up test data
def gauss(x,x0, FWHM, A):
	norm = []
	sd=FWHM*1.0/(2*np.sqrt(2*np.log(2)))
	for i in range(x.size):
		norm += [A*np.exp(-(x[i] - x0)**2/(2*sd**2))]
	return np.array(norm)

c=299792.458 #light speed in m/s
spect=np.loadtxt('Data/16jan2011_NOT.norm.dat').transpose()

x=spect[0]
y_real=spect[1]-1

#Selection only Halpha
x_min=6550
x_max=6700

ind_cut=np.where((x>x_min) & (x<x_max))[0]
x=x[ind_cut]
y_real=y_real[ind_cut]

######################################
# Solving
m1, fwhm1, A1, m2, fwhm2, A2,m3, fwhm3,A3 = [6620, 5, 10, 6625, 3, 3, 6625,5,2]#m=x0, sd1=FWHM
p = [m1, fwhm1, A1, m2, fwhm2, A2, m3, fwhm3, A3] # Initial guesses for leastsq
# Number of Gaussian
y_init = gauss(x,m1, fwhm1, A1) + gauss(x,m2, fwhm2, A2)+gauss(x,m3, fwhm3, A3) # For final comparison plot

#### fit with solution
def res(p, y, x):
	m1, fwhm1, A1, m2, fwhm2, A2,m3, fwhm3, A3 = p
	y_fit = gauss(x,m1, fwhm1, A1) + gauss(x,m2, fwhm2, A2)+gauss(x,m3, fwhm3, A3)
	err = y - y_fit
	return err

residual=res(p,y_real,x)
plsq = leastsq(res, p, args = (y_real, x))

#Velocity of each component
velocity_narrow=(plsq[0][1]*1.0/plsq[0][0])*c
velocity_inter=(plsq[0][4]*1.0/plsq[0][3])*c
velocity_broad=(plsq[0][7]*1.0/plsq[0][6])*c

xx=np.arange(min(x),max(x),1)
y_est = gauss(xx, plsq[0][0], plsq[0][1],plsq[0][2]) + gauss(xx, plsq[0][3], plsq[0][4],plsq[0][5])+gauss(xx, plsq[0][6], plsq[0][7],plsq[0][8])

fig, ax1 = plt.subplots(figsize=(8,8), facecolor='w', edgecolor='k')
ax1.plot(x, y_real,'k',label=r'\textbf{Data}')
ax1.plot(xx, gauss(xx, plsq[0][0], plsq[0][1],plsq[0][2]),'b:',alpha=0.8,label=r'\textbf{Narrow %s km/s}'%int(velocity_narrow))
ax1.plot(xx,gauss(xx, plsq[0][3], plsq[0][4],plsq[0][5]), 'r.',alpha=0.8,label=r'\textbf{Inter %s km/s}'%int(velocity_inter))
ax1.plot(xx, gauss(xx, plsq[0][6], plsq[0][7],plsq[0][8]), 'g-.',alpha=0.8,label=r'\textbf{Broad %s km/s}'%int(velocity_broad))

ax1.plot(xx, y_est, 'm--', label=r'\textbf{Sum of components}')
ax1.set_xlabel(r'\textbf{Observed wavelength [\AA]}')
ax1.set_ylabel(r'\textbf{Normalized flux]}')
ax1.legend(loc=0,title='',markerscale=.5,ncol=1,prop={'size':15})
plt.show()
print ('FWHM narrow, intermediate, broad in Angs')
print (plsq[0][1],plsq[0][4],plsq[0][7])

print ('Velocity narrow component %s km/s'%int(velocity_narrow))
print ('Velocity intermediate component %s km/s'%int(velocity_inter))
print ('Velocity broad component %s km/s'%int(velocity_broad))


