import os
import matplotlib.pyplot as p
from numpy import *
import numpy as np
from operator import itemgetter
import mpfit3
import math
from scipy.stats import f as scipyf
import subprocess
import sys
from scipy.interpolate import interp1d
from scipy.stats import chisquare as scipy_chi2


########## FUNCTIONS TO FIT ##########
def oneline(x,p2):
	i = 0
	v1 = np.zeros(len(x))
	for i in range(len(x)):
		v1[i] = p2[0]*x[i]+p2[1]
	return v1

def onelinee(x,p1,p2):
	return p1*x+p2

def twoline(x,p1):
	v = np.zeros(len(x))
	i = 0
	for i in range(len(x)):
		if x[i]<p1[3]:
			v[i] = p1[0]*x[i] + p1[2]
		elif x[i]>=p1[3]:
			v[i] = p1[1]*x[i] + p1[3]*(p1[0]-p1[1]) + p1[2]
	return v

def twolinee(x,p1,p2,p3,p4):
	if x<p4:
		return p1*x+p3
	elif x>=p4:
		return p2*x+p4*(p1-p2)+p3

def olivares1(x,pt):
	i = 0
	f = np.zeros(len(x))
	for i in range(len(x)):
		fd = -pt[0]/(1.+np.exp((x[i]-pt[1])/pt[2]))
		l = pt[3]*(x[i]-pt[1])+pt[4]
		f[i] = fd+l
	return f

def olivares2(x,pt):
	i = 0
	f = np.zeros(len(x))
	for i in range(len(x)):
		g = -pt[5]*math.exp(-(((x[i]-pt[6])/pt[7])**2))
		fd = -pt[0]/(1.+math.exp((x[i]-pt[1])/pt[2]))
		l = pt[3]*(x[i]-pt[1])+pt[4]
		f[i] = g+fd+l
	return f

def olivares1es(x,pt):
	fd = -pt[0]/(1.+math.exp((x-pt[1])/pt[2]))
	l = pt[3]*(x-pt[1])+pt[4]
	f = fd+l
	#f = 10.**(-0.4*y)
	return f

def olivares2es(x,pt):
	g = -pt[5]*math.exp(-(((x-pt[6])/pt[7])**2))
	fd = -pt[0]/(1.+np.exp((x-pt[1])/pt[2]))
	l = pt[3]*(x-pt[1])+pt[4]
	f = fd+l+g
	#f = 10.**(-0.4*y)
	return f


######### FUNCTIONS TO FIT AND USE MPFIT ########
def myfunct(p1, fjac=None, x=None, y=None, err=None):
	model = twoline(x,p1)
	status=0
	return ([status,(y-model)/err])

def myfuncte0(p1, fjac=None, x=None, y=None):
	model = twoline(x,p1)
	status=0
	return ([status,y-model])

def myfunct1(p2, fjac=None, x=None, y=None, err=None):
	model = oneline(x,p2)
	status=0
	return ([status,(y-model)/err])

def myfunct1e0(p2, fjac=None, x=None, y=None):
	model = oneline(x,p2)
	status=0
	return ([status,y-model])

def myfunctoliv1(pt, fjac=None, x=None, y=None, err=None):
	model = olivares1(x,pt)
	status=0
	return([status,(y-model)/err])

def myfunctoliv2(pt, fjac=None, x=None, y=None, err=None):
	model = olivares2(x,pt)
	status=0
	return([status,(y-model)/err])

def chi2(Data,p1,p2):
	s = 0.
	for x,y,yerr in Data:
		teoria = onelinee(x,p1,p2)
		if yerr==0.:
			s += (teoria-y)**2
		else:
			s += (teoria-y)**2/yerr**2
	return s

def chi22(Data,p1,p2,p3,p4):
	s = 0.
	for x,y,yerr in Data:
		teoria = twolinee(x,p1,p2,p3,p4)
		if yerr==0.:
			s += (teoria-y)**2
		else:
			s += (teoria-y)**2/yerr**2
	return s

####### data!! ########

data=
name=
filter_name=
output
'''

data= epo, mag and err mag
name= supernova name
filter_name= name filter B,V,g,r,...
output= output file name

'''

table = open('%s.txt'%output,'a')
table.write('#SN     s1      err_s1 s2 err_s2 s err_s s3 err_s3   OPTd  err_OPTd \n')


days,magnitude,error=data
name_SN=name

days,magnitude,error=zip(*sorted(zip(days, magnitude, error)))


## Figure :
print('Left click if you want to do the fit')
f,(ax1)=p.subplots(1,sharex=True,sharey=False)
f.subplots_adjust(hspace=0.2)

ea = max(error)
ymax = max(magnitude) + ea + 0.05
ymin = min(magnitude) - ea - 0.05

ax1.errorbar(days,magnitude,yerr = error,fmt='.', label = '%s'%(filter_name))
ax1.set_ylabel('magnitude')
ax1.set_xlabel('days')
ax1.set_ylim((ymax,ymin))
ax1.set_xlim((np.min(days)-5,np.max(days)+5))
ax1.set_title(name_SN)
ax1.legend()
def onclick(event):
	global do_fit
	if event.button==1:
		do_fit=str('y')
		p.close()
	if event.button==3:
		do_fit=str('n')
		p.close()
cid = f.canvas.mpl_connect('button_press_event', onclick)
p.show()

if do_fit=='y':
	
	f,(ax1)=p.subplots(1,sharex=True,sharey=False)
	f.subplots_adjust(hspace=0.2)

	ea = max(error)
	ymax = max(magnitude) + ea + 0.05
	ymin = min(magnitude) - ea - 0.05

	ax1.errorbar(days,magnitude,yerr = error,fmt='.', label = '%s'%(filter_name))
	ax1.set_ylabel('magnitude')
	ax1.set_xlabel('days')
	ax1.set_ylim((ymax,ymin))
	ax1.set_xlim((np.min(days)-5,np.max(days)+5))
	ax1.set_title(name_SN)
	ax1.legend()
	print("Click Tmax, Toptd, Tstart s3, Tfinal s3")
	print("Left click put point, right click 0 for each time, click right and left remove last value")
	epoch=[]
	def onclick(event):
		global epoch
		if event.button==1:
			epoch.append(event.xdata)
			print('clicked!'+str(event.xdata))
		if event.button==2:
			del epoch[-1]
		if event.button==3:
			epoch.append(0)
			epoch.append(0)
			epoch.append(0)
			epoch.append(0)	

	cid = f.canvas.mpl_connect('button_press_event', onclick)
	p.show()
	p.close()
	ts1 = float(epoch[0])
	ts2 = float(epoch[1])
	ts31 = float(epoch[2])
	ts32 = float(epoch[3])
	print(round(ts1,3),round(ts2,3),round(ts31,3),round(ts32,3))


	print('Doing s1 and s2')
	#if input('s1,s2? (|yes(y)|,no(n)) ') !='n':

	p.show(block=True)

	Data = []
	xx = []
	nn = 0

	for o,x in enumerate(days):
		x1 = x
		xx.append(x1)
		if x>=ts1 and x<=ts2:
			nn += 1

	X = np.zeros(nn)
	Y = np.zeros(nn)
	Ye = np.zeros(nn)
	N = 0
	cero = 0

	for o,x in enumerate(days):
		if x>=ts1 and x<=ts2:
			X[N] = x
			Y[N] = magnitude[o]
			if error[o]==0:
				cero = 1
			else:
				Ye[N] = error[o]
			N+=1
			Data.append((x,magnitude[o],error[o]))

	if cero==0:

		p00 = [0.02,16.]
		fa = {'x':X,'y':Y,'err':Ye}
		fab = {'x':X,'y':Y,'err':Ye}

		mm = mpfit3.mpfit(myfunct1, p00, functkw=fa)
		print('parameters= ',mm.params)

		st = mm.params[0]
		b = mm.params[1]

		p0 = [st, st, b, (ts2-ts1)*0.2]

		m = mpfit3.mpfit(myfunct, p0, functkw=fab)
		print('parameters= ',m.params)
		ttrans = m.params[3]
		print('ttrans= ',ttrans)

		tiem = linspace(ts1,ts2,1000)

		p.errorbar(xx,magnitude,yerr = error,fmt='k.', label = filter_name)
		p.plot(tiem,twoline(tiem,m.params),'-r',label='Fit dos')
		p.plot(tiem,oneline(tiem,mm.params),'-b',label='Fit una')
		p.ylim((ymax,ymin))
		p.legend()

	else:

		p00 = [0.02,16.]
		fa = {'x':X,'y':Y}
		fab = {'x':X,'y':Y}

		mm = mpfit3.mpfit(myfunct1e0, p00, functkw=fa)
		print('parameters= ',mm.params)

		st = mm.params[0]
		b = mm.params[1]

		p0 = [st, st, b, (ts2-ts1)*0.2]

		m = mpfit3.mpfit(myfuncte0, p0, functkw=fab)
		print('parameters= ',m.params)

		tiem = linspace(ts1,ts2,1000)

		p.errorbar(xx,magnitude,yerr = error,fmt='k.', label = filter_name)
		p.plot(tiem,twoline(tiem,m.params),'-r',label='Fit dos')
		p.plot(tiem,oneline(tiem,mm.params),'-b',label='Fit una')
		p.ylim((ymax,ymin))
		p.legend()
		    
	p.show()

	N = nn # number of data points
	SSR1 = chi2(Data,mm.params[0],mm.params[1]) # chi-square of oneline
	SSR2 = chi22(Data,m.params[0],m.params[1],m.params[2],m.params[3]) # chi-square of twoline

	ftest = ((SSR1-SSR2)/2.)/(SSR2/(N-4.))
	print('oneline s2 = ',mm.params[0]*100.,'err s2= ',mm.perror[0]*100., 'b= ',mm.params[1])
	if nn>4:
		print('twoline s1= ', m.params[0]*100.,'err s1= ',m.perror[0]*100.,'twoline s2= ',m.params[1]*100.,'err s2= ',m.perror[1]*100., 'b= ', m.params[2])
		print('ttrans= ',m.params[3],'err ttrans= ', m.perror[3])
	else:
		print('twoline s1= ', m.params[0]*100.,'twoline s2= ',m.params[1]*100.,'err s2= ', 'b= ', m.params[2])
		print('ttrans= ',m.params[3])
	print('F-test= ',ftest, 'ddf= ',nn-4., 'p-value= ',scipyf.cdf(ftest, 2, nn-4))
	if scipyf.cdf(ftest, 2, nn-4) > 0.95:
	#if ftest>scipyf.ppf(0.95,2, nn-4):
		print('2-slopes better')
		good_slope=2
	else:
		print('1-slope better')
		good_slope=1

	if good_slope == 1:
		s1 = 999.9
		errs1 = 999.9
		s2 = mm.params[0]*100.
		errs2 = mm.perror[0]*100.
		b1 = mm.params[1]
		s = mm.params[0]*100.
		errs = mm.perror[0]*100.
		b = mm.params[1]
		ttrans = 999.9
		errttrans = 999.9
	elif good_slope == 2:
		s1 = m.params[0]*100.
		errs1 = m.perror[0]*100.
		b1 = m.params[2]
		s2 = m.params[1]*100.
		errs2 = m.perror[1]*100.
		s = mm.params[0]*100.
		errs = mm.perror[0]*100.
		b = mm.params[1]
		ttrans = m.params[3]
		errttrans = m.perror[3]
	else:
		print('Not an option')

	
	#OPTD
	do_optd=input('Doing optd ? yes or no (y or n) ')
	if do_optd=='y':
		#Extrapoloation s2 to last point
		if nn>4:
			epoch_s2=np.arange(epoch[1]-10,max(xx),0.1)
			if good_slope==2:	
				s2_value=twoline(epoch_s2,m.params)
				s2_value_err=twoline(epoch_s2,m.params+m.perror)
			else:
				s2_value=oneline(epoch_s2,mm.params)
				s2_value_err=oneline(epoch_s2,mm.params+mm.perror)

			f_mag = interp1d(xx, magnitude, kind='linear')

			f_s2 = interp1d(epoch_s2, s2_value, kind='linear')

			diff_mag=abs(s2_value-f_mag(epoch_s2))
			diff_mag_err=abs(s2_value_err-f_mag(epoch_s2))

			if size(where(diff_mag>0.1)[0])>0:

				OPTd=epoch_s2[where(diff_mag>0.1)[0][0]]
				err_OPTd=sqrt((OPTd-epoch_s2[where(diff_mag_err>0.1)[0][0]])**2+(err_JD_explosion_kait[SN])**2)

				f,(ax1)=p.subplots(1,sharex=True,sharey=False)
				f.subplots_adjust(hspace=0.2)

				ea = max(error)
				ymax = max(magnitude) + ea + 0.05
				ymin = min(magnitude) - ea - 0.05

				ax1.errorbar(np.array(xx),magnitude,yerr = error,fmt='.', label = '%s'%(filter_name))
				ax1.plot(epoch_s2, s2_value,'r--')
				ax1.axvline(x=OPTd, color='g', linestyle='--')
				ax1.set_ylabel('magnitude')
				ax1.set_xlabel('days')
				ax1.set_ylim((ymax,ymin))
				ax1.set_xlim((np.min(days)-5,np.max(days)+5))
				ax1.set_title(name_SN)
				ax1.legend()
				p.show()
			else:
				OPTd=999.9
				err_OPTd=999.9
		else:
			OPTd=999.9
			err_OPTd=999.9
	else:
		OPTd=999.9
		err_OPTd=999.9
	do_s3=input('Doing s3 ? yes or no (y or n) ')
	if do_s3 ==str('y'):
		nnn = 0
		cero=0
		xx = []
	    
		for o,x in enumerate(days):
			x1 = x-ts31
			xx.append(x1)
			if x>=ts31 and x<=ts32:
				nnn += 1

		X_s3 = np.zeros(nnn)
		Y_s3 = np.zeros(nnn)
		Ye_s3 = np.zeros(nnn)
		N = 0

		for o,x in enumerate(days):
			if x>=ts31 and x<=ts32:
				X_s3[N] = x-ts31
				Y_s3[N] = magnitude[o]
				if error[o]==0:
					cero=1
				else:
					Ye_s3[N] = error[o]
				N+=1

		if cero==0:
		    
			fb = {'x':X_s3,'y':Y_s3,'err':Ye_s3}
			p01 = [0.01, 13.5]
			m3 = mpfit3.mpfit(myfunct1, p01, functkw=fb)

			tiem = linspace(0,ts32-ts31,1000)

			p.errorbar(xx,magnitude,yerr = error,fmt='k.', label = filter_name)
			p.plot(tiem,oneline(tiem,m3.params),'-b',label='Fit s3')
			p.ylim((ymax,ymin))
			p.legend()
			p.show()


		else:
		    
			fb = {'x':X_s3,'y':Y_s3}
			p01 = [0.01, 13.5]
			m3 = mpfit3.mpfit(myfunct1e0, p01, functkw=fb)

			tiem = linspace(0,ts32-ts31,1000)

			p.errorbar(xx,magnitude,yerr = error,fmt='k.', label = filter_name)
			p.plot(tiem,oneline(tiem,m3.params),'-b',label='Fit s3')
			p.ylim((ymax,ymin))
			p.legend()
			p.show()

		print('s3= ',m3.params[0]*100.,'err s3= ',m3.perror[0]*100., 'b3= ', m3.params[1])
		s3 = m3.params[0]*100.
		errs3 = m3.perror[0]*100.
		b3 = m3.params[1]
	else:
		s3=999.9
		errs3=999.9

	### CALCULAR TPT, MEND, MTAIL
	do_tpt=input('Doing Tpt ? yes or no (y or n) ')
	if do_tpt==str('y'):
		n = len(days)
		days1 = np.zeros(n)
		magnitude1 = np.zeros(n)
		error1 = np.zeros(n)
		for k in range(len(days)):
			days1[k] = days[k]
			magnitude1[k] = magnitude[k]
			error1[k] = error[k]
			k+=1

		#print(days, magnitude, error)
		#print(days1, magnitude1, error1)
		fbb = {'x':days1,'y':magnitude1,'err':error1}
		p0 = {'value':3., 'fixed':0,'limited':[1,0],'limits':[0.001,10.]} #step FD (a_0)
		p1 = {'value':90., 'fixed':0,'limited':[0,0],'limits':[40.,150.]} #middle of transition phase FD (t_PT)
		p2 = {'value':6., 'fixed':0,'limited':[1,0],'limits':[1.,12.]} #width of the transition phase FD (w_0)
		p3 = {'value':1.000, 'fixed':0,'limited':[1,1],'limits':[0.,3.]} #slope of radioactive decay  
		p4 = {'value':4., 'fixed':0,'limited':[1,0],'limits':[0.,0.]} #zero point at t=t_pt
		p5 = {'value':1.3, 'fixed':0,'limited':[1,1],'limits':[0.001,2.]} #gaussian height
		p6 = {'value':0., 'fixed':0,'limited':[1,1],'limits':[-15,10.]}#center of gaussian
		p7 = {'value':20., 'fixed':0,'limited':[1,0],'limits':[0.00001,3.]} #gaussian width
		pp = [p0,p1,p2,p3,p4,p5,p6,p7]
		vector = linspace(0,days1[len(days1)-1],1000)
		vector2 = linspace(min(magnitude1),max(magnitude1),100)
		f = mpfit3.mpfit(myfunctoliv2, functkw=fbb, parinfo=pp)
		##
		siga = f.perror[0]/(1+math.exp(-30./f.params[2]))
		sigtpt = 0.#f.params[0]*math.exp(-30./f.params[2])*f.perror[1]/(f.params[2]*(math.exp((f.params[1]-30)/f.params[2])+math.exp(-30./f.params[2])))
		sigw =f.params[0]*-30.*math.exp(-30./f.params[2])*f.perror[2]/((f.params[2]*(math.exp(-30./f.params[2])+1))**2)
		##
		tpt = f.params[1]
		errtpt = f.perror[1]
		w0 = f.params[2]
		ew0 = f.perror[2]
		a0 = f.params[0]
		ea0 = f.perror[0]
		mtail = olivares2es(f.params[1]+3.0*w0,f.params)
		mtailerr =math.sqrt(((f.params[1]+3.0*w0)*f.perror[3])**2 + f.perror[4]**2)
		mend = olivares2es(f.params[1]-3.0*w0,f.params)
		menderr= math.sqrt(siga**2 + sigtpt**2 + sigw**2)
		vectortpt = np.zeros(len(vector2))
		for m in range(len(vector2)):vectortpt[m]=f.params[1]
		##
		print('tpt= ', tpt, 'errtpt= ', errtpt)
		print('mend= ', mend, 'errmend= ', menderr, 'tend= ', tpt-3.0*w0)
		print('mtail= ', mtail, 'errmatil= ', mtailerr, 'ttail= ', tpt+3.0*w0)
		print('w0= ', f.params[2], 'w0err= ', f.perror[2])
		print('a0= ', f.params[0], 'a0err= ', f.perror[0])
		print('s3= ', f.params[3], 's3err= ', f.perror[3], '\n')

		p.errorbar(days1,magnitude1,yerr = error1,fmt='.k',label =filter_name)
		p.plot(vector,olivares2(vector,f.params),'-b',label= 'Fit Tpt')
		p.plot(vectortpt,vector2,'-r',label='Tpt')
		p.plot(f.params[1]-3.0*w0,mend,'og',label='Mend')
		p.plot(f.params[1]+3.0*w0,mtail,'og',label='Mtail')
		p.ylim((max(magnitude1)+max(error1)+0.5,min(magnitude1)-max(error1)-0.5))
		p.legend()
		p.title(name_SN)
		p.show()
		optd=tpt-3.0*w0
		err_optd=sqrt(errtpt**2+err_JD_explosion_kait[SN]**2)

		s3_fit=f.params[3]*100
		errs3_fit=f.perror[3]*100

		if input('is it good?( |yes(y)|,no(n)) ')!='n':
			pass
		else:
			p.errorbar(days1,magnitude1,yerr = error1,fmt='.k',label =filter_name)
			p.ylim((max(magnitude1)+max(error1)+0.5,min(magnitude1)-max(error1)-0.5))
			p.legend()
			p.title(name_SN)
			p.show(block=False)
			tinicio = float(input('initial time? '))
			tfinal = float(input('final time? '))
			p.show(block=True)
			uop = 0
			for k in range(len(days1)):
				if days1[k] >= tinicio and days1[k]<=tfinal:
					uop+=1

			days2 = np.zeros(uop)
			magnitude2 = np.zeros(uop)
			error2 = np.zeros(uop)

			u=0
			for k in range(len(days1)):
				if days1[k] >= tinicio and days1[k]<=tfinal:
					days2[u] = days1[k]
					magnitude2[u] = magnitude1[k]
					error2[u] = error1[k]
					u+=1

			fbb = {'x':days2,'y':magnitude2,'err':error2}
			p0 = {'value':3., 'fixed':0,'limited':[1,0],'limits':[0.001,10.]} #peldagno fd
			p1 = {'value':90., 'fixed':0,'limited':[1,0],'limits':[40.,150.]} #tiempo de transicion fd (tpt)
			p2 = {'value':6., 'fixed':0,'limited':[1,0],'limits':[1.,12.]} #ancho de la transicion FD
			p3 = {'value':1.0, 'fixed':0,'limited':[1,1],'limits':[0.,5.]} #decaimiento tau
			p4 = {'value':4., 'fixed':0,'limited':[1,0],'limits':[0.,0.]} #pendiente decaimiento
			pp = [p0,p1,p2,p3,p4]
			vector = linspace(0,days1[len(days1)-1],1000)
			vector2 = linspace(min(magnitude1),max(magnitude1),100)
			f = mpfit3.mpfit(myfunctoliv1, functkw=fbb, parinfo=pp)
			siga = f.perror[0]/(1+math.exp(-30./f.params[2]))
			sigtpt = f.params[0]*math.exp(-30./f.params[2])*f.perror[1]/(f.params[2]*(math.exp((f.params[1]-30)/f.params[2])+math.exp(-30./f.params[2])))
			sigw =f.params[0]*-30.*math.exp(-30./f.params[2])*f.perror[2]/((f.params[2]*(math.exp(-30./f.params[2])+1))**2)

			tpt = f.params[1]
			w0 = f.params[2]
			ew0 = f.perror[2]
			a0 = f.params[0]
			ea0 = f.perror[0]
			errtpt = f.perror[1]
			mtail = olivares1es(f.params[1]+3.0*w0,f.params)
			mtailerr =math.sqrt(((f.params[1]+3.0*w0)*f.perror[3])**2 + f.perror[4]**2)
			mend = olivares1es(f.params[1]-3.0*w0,f.params)
			menderr= math.sqrt(siga**2 + sigtpt**2 + sigw**2)
			vectortpt = np.zeros(len(vector2))
			for m in range(len(vector2)):vectortpt[m]=f.params[1]

			print('tpt= ', tpt, 'errtpt= ', errtpt)
			print('mend= ', mend, 'errmend= ', menderr, 'tend= ', tpt-3.0*w0)
			print('mtail= ', mtail, 'errmatil= ', mtailerr, 'ttail= ', tpt+3.*w0)
			print('w0= ', f.params[2], 'w0err= ', f.perror[2])
			print('a0= ', f.params[0], 'a0err= ', f.perror[0])
			print('s3= ', f.params[3]*100., 's3err= ', f.perror[3]*100.,'\n')
			#print('parameters and errors = ', '%e' %f.params[0], '%.2f' %f.perror[0], '%e' %f.params[1], '%.2f' %f.perror[1], '%e' %f.params[2], '%.2f' %f.perror[2], '%e' %f.params[3], '%.2f' %f.perror[3], '%e' %f.params[4],  '%.2f' %f.perror[4])

			p.errorbar(days1,magnitude1,yerr = error1,fmt='.k',label =filter_name)
			p.plot(vector,olivares1(vector,f.params),'-b',label= 'Fit Tpt')
			p.plot(vectortpt,vector2,'-r',label='Tpt')
			p.plot(f.params[1]-3.0*w0,mend,'og',label='Mend')
			p.plot(f.params[1]+3.0*w0,mtail,'og',label='Mtail')
			p.ylim((max(magnitude1)+max(error1)+0.5,min(magnitude1)-max(error1)-0.5))
			p.legend()
			p.title(name_SN)
			p.show()
			s3_fit=f.params[3]*100
			errs3_fit=f.perror[3]*100
			if input('is it good?( |yes(y)|,no(n)) ')!='n':
				optd=tpt-3.0*w0
				err_optd=sqrt(errtpt**2+err_JD_explosion_kait[SN]**2)
				s3_fit=f.params[3]*100
				errs3_fit=f.perror[3]*100
			else:
				s3_fit=999.9
				errs3_fit=999.9
				optd=999.9
				err_optd=999.9

	else:
		optd=999.9
		err_optd=999.9
		s3_fit=999.9
		errs3_fit=999.9
	f,(ax1)=p.subplots(1,sharex=True,sharey=False)
	f.subplots_adjust(hspace=0.2)

	ea = max(error)
	ymax = max(magnitude) + ea + 0.05
	ymin = min(magnitude) - ea - 0.05

	ax1.errorbar(days,magnitude,yerr = error,fmt='.', label = '%s'%(filter_name))
	ax1.set_ylabel('magnitude')
	ax1.set_xlabel('days')
	ax1.set_ylim((ymax,ymin))
	ax1.set_xlim((np.min(days)-5,np.max(days)+5))
	ax1.set_title(name_SN)

	if good_slope==2:
		ax1.plot(X,twoline(X,[s1/100,s2/100,b1,ttrans]),'-r',label='Fit dos')
	else:
		ax1.plot(X,oneline(X,[s2/100,b1]),'-b',label='Fit una')

	if do_s3=='y':
		ax1.plot(X_s3+ts31,oneline(X_s3,m3.params),'-m',label='Fit s3')
	if optd!=999.9:
		print('fit O10 OPTd=%s +-/ %s '%(round(optd,2),round(err_optd,2)))
		ax1.axvline(x=optd, color='brown', linestyle='--',label='OPTd')
	print('fit old OPTd=%s +-/ %s '%(round(OPTd,2),round(err_OPTd,2)))

	ax1.axvline(x=OPTd, color='g', linestyle='--',label='OPTd old')
	ax1.legend()
	p.show()
	optd_val=input('Valeur OPTD old or fit: (o) or (f)')
	if optd_val==str('f'):
		optd=optd
		err_optd=err_optd
	else:
		optd=OPTd
		err_optd=err_OPTd

	if do_tpt==str('y'):
		print('s3 fit= ', s3_fit, 's3err= ',errs3_fit)
		print('s3= ', s3, 's3err= ', errs3)
	else:
		print('s3= ', s3, 's3err= ', errs3)
	s3_val=input('Valeur s3 old or fit: (o) or (f)')
	if s3_val==str('f'):
		s3=s3_fit
		errs3=errs3_fit
	else:
		s3=s3
		errs3=errs3



else:
	s1=999.9
	errs1=999.9
	s2=999.9
	errs2=999.9
	s=999.9
	errs=999.9
	s3=999.9
	errs3=999.9
	optd=999.9
	err_optd=999.9


########## SAVE #############
ps1 = '%.2f' %s1
perrs1 = '%.2f' %errs1
ps2 = '%.2f' %s2
perrs2 = '%.2f' %errs2
ps = '%.2f' %s
perrs = '%.2f' %errs
ps3 = '%.2f' %s3
perrs3 = '%.2f' %errs3
poptd = '%.2f' %optd
perroptd = '%.2f' %err_optd
print(name_SN.ljust(10,' ')+ps1.ljust(7,' ')+perrs1.ljust(7,' ')+ps2.ljust(7,' ')+perrs2.ljust(7,' ')+ps.ljust(7,' ')+perrs.ljust(7,' ')+ps3.ljust(7,' ')+perrs3.ljust(7,' ')+poptd.ljust(10,' ')+perroptd.ljust(7,' '))
table.write(name_SN.ljust(10,' ')+ps1.ljust(7,' ')+perrs1.ljust(7,' ')+ps2.ljust(7,' ')+perrs2.ljust(7,' ')+ps.ljust(7,' ')+perrs.ljust(7,' ')+ps3.ljust(7,' ')+perrs3.ljust(7,' ')+poptd.ljust(10,' ')+perroptd.ljust(7,' ')+'\n')
table.close()

