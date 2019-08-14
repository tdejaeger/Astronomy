#!/usr/bin/python
import os
import matplotlib.pyplot as p
from numpy import *
import numpy as np
from operator import itemgetter
import math
########## FUNCTIONS TO FIT ##########

def fit(tiempos):
    
	t1 =tiempos[0]#float(input('lim time1: '))
	t2 =tiempos[1]#float(input('lim time2: '))
    
	print('Close the plot')
	p.show(block=True)
	ptos=0

	for i,a in enumerate(days):
		if a>=t1 and a<=t2:
			ptos+=1
                
	at = np.zeros(ptos)
	att = []
	am = np.zeros(ptos)
	er = np.zeros(ptos)
	j = 0
    
	for i,a in enumerate(days):
		if a>=t1 and a<=t2:
			att.append((days[i],magnitude[i],error[i]))

	att.sort(key=itemgetter(0))
	for i in range(len(att)):
		tupla = att[i]
		at[i] = tupla[0]
		am[j] = tupla[1]
		er[j] = tupla[2]
		j+=1
    
	k=0
	color_fit=['r','b','g','brown']
	for j in range(2,6):
		if ptos>=j+2:
			try:
				aj = polyfit(at,am,j)
				po = poly1d(aj)
				eam = max(error)
				ymaxh = max(am) + eam + 0.05
				yminh = min(am) - eam - 0.05
				x = linspace(at[0]-5,at[len(at)-1]+5,100)
				p.plot(at,am,'.',x,po(x),'-',color=color_fit[j-2],label='p'+str(j))
				#p.gca().invert_yaxis()
				#p.set_ylim((ymaxh,yminh))
				ax = p.gca()
				ax.set_ylim((ymaxh,yminh))
			except ValueError:
				print('It is not possible to make a polynomial grade'+str(j))
				k = k+1
	p.plot(at,am,'ko')

	if k==3:
		print('Problem with fit')
	else:
		p.legend()
		p.show(block=False)
		maxi = []
		good=0
		des = input('polynomial you choose? p2, p3,etc ')
		if des == 'p2':
			aj = polyfit(at,am,2) 
			der = polyder(aj)
			r = roots(der)
			for raiz in r:
				val = polyval(aj,raiz)
				maxi.append([raiz,val])
		elif des == 'p3':
			aj = polyfit(at,am,3)
			der = polyder(aj)
			r = roots(der)
			for raiz in r:
				val = polyval(aj,raiz)
				maxi.append([raiz,val])
			good=1
		elif des == 'p4':
			aj = polyfit(at,am,4)
			der = polyder(aj)
			r = roots(der)
			for raiz in r:
				val = polyval(aj,raiz)
				maxi.append([raiz,val])
			good=1
		elif des == 'p5':
			aj = polyfit(at,am,5)
			der = polyder(aj)
			r = roots(der)
			for raiz in r:
				val = polyval(aj,raiz)
				maxi.append([raiz,val])
			good=1
		else:
			print('Invalid')

	print(maxi)
	s=0.
	for err in er:
		s+=err
	erro = s/(len(er))
	ind_max=where(np.real(maxi)[:,1]==min(np.real(maxi)[:,1]))[0][0]
	print('Time, Maximum magnitudee, Error:')
	print(np.real(maxi)[:,0][ind_max],np.real(maxi)[:,1][ind_max],erro)
	p.show(block=True)

	return [np.real(maxi)[:,1][ind_max],np.real(maxi)[:,0][ind_max],erro]
'''

data=[[epo],[mag],[err]]

'''

def fit_mmax(data,name,filt):

	tabla = open('%s_%s_mmax.txt'%(name,filt),'a')
	days,magnitude,error=data
	name_SN=name

	days,magnitude,error=zip(*sorted(zip(days, magnitude, error)))

	ea = max(error)
	ymax = max(magnitude) + ea + 0.05
	ymin = min(magnitude) - ea - 0.05


	magmax = magnitude[0]
	tmagmax = days[0]
	magmaxe = error[0]
	banda=filt
	print(tmagmax, magmax, magmaxe)

	if tmagmax<20: #Maximum should be inferior a 20 days after the explosion

		f,(ax1)=p.subplots(1,sharex=True,sharey=False)
		f.subplots_adjust(hspace=0.2)
		ax1.errorbar(np.array(days),magnitude,yerr = error,fmt='.', label = '%s'%(banda))
		ax1.set_ylabel('magnitude')
		ax1.set_xlabel('days')
		ax1.set_ylim((ymax,ymin))
		ax1.set_xlim((np.min(np.array(days))-5,np.max(np.array(days))+5))
		ax1.set_title(name_SN)
		ax1.legend()
		ax1.scatter(tmagmax,magmax,color='red',s=60)#,fmt='.')
		p.show()

		if input('is it good?( |yes(y)|,no(n)) ') =='n':

			if input('polynomial fit? ( |yes(y)|,no(n)) ') =='y':

				f,(ax1)=p.subplots(1,sharex=True,sharey=False)
				f.subplots_adjust(hspace=0.2)
				ax1.errorbar(days,magnitude,yerr = error,fmt='.', label = '%s'%(banda))
				ax1.set_ylabel('magnitude')
				ax1.set_xlabel('days')
				ax1.set_ylim((ymax,ymin))
				ax1.set_xlim((np.min(days)-5,np.min(days)+60))
				ax1.set_title(name_SN)
				ax1.legend()
				ax1.scatter(tmagmax,magmax,color='red',s=60)
				tiempos=[]
				print("define limit fit maximum: Left click put point, right click 0 for each time, click right and left remove last value")
				def onclick(event):
					global tiempos
	
					if event.button==1:
						tiempos.append(event.xdata)
						print(str(event.xdata))
					if event.button==2:
						del tiempos[-1]
					if event.button==3:
						tiempos.append(0)
						tiempos.append(0)
				cid = f.canvas.mpl_connect('button_press_event', onclick)
				p.show()
				res = fit(tiempos)
				if res == 0:
					print('Cannot be fit with polyfit')
					magmax = nan
					tmagmax = nan
					magmaxe = nan
				else:
					magmax = res[0]
					tmagmax = res[1]
					magmaxe = res[2]

				if input('is it good?( |yes(y)|,no(n)) ') =='n':

					print("Click on the maximum")

					f,(ax1)=p.subplots(1,sharex=True,sharey=False)
					f.subplots_adjust(hspace=0.2)
					ax1.errorbar(days,magnitude,yerr = error,fmt='.', label = '%s'%(banda))
					ax1.set_ylabel('magnitude')
					ax1.set_xlabel('days')
					ax1.set_ylim((ymax,ymin))
					ax1.set_xlim((np.min(days)-5,np.min(days)+160))
					ax1.set_title(name_SN)
					ax1.legend()
					ax1.scatter(tmagmax,magmax,color='green',s=60)#,fmt='.')  
					tiempos=[]
					def onclick(event):
						global tiempos
						if event.button==1:
							tiempos=(event.xdata,event.ydata)
							print(tiempos)
					cid = f.canvas.mpl_connect('button_press_event', onclick)
					p.show()
					tmagmax = tiempos[0]
					magmax = tiempos[1]
					magmaxe = 0.1

	
			else:
				print("Click on the maximum")

				f,(ax1)=p.subplots(1,sharex=True,sharey=False)
				f.subplots_adjust(hspace=0.2)
				ax1.errorbar(days,magnitude,yerr = error,fmt='.', label = '%s'%(banda))
				ax1.set_ylabel('magnitude')
				ax1.set_xlabel('days')
				ax1.set_ylim((ymax,ymin))
				ax1.set_xlim((np.min(days)-5,np.min(days)+160))
				ax1.set_title(name_SN)
				ax1.legend()
				ax1.scatter(tmagmax,magmax,color='green',s=60)#,fmt='.')  
				tiempos=[]
				def onclick(event):
					global tiempos
					if event.button==1:
						tiempos=(event.xdata,event.ydata)
						print(tiempos)
				cid = f.canvas.mpl_connect('button_press_event', onclick)
				p.show()
				tmagmax = tiempos[0]
				magmax = tiempos[1]
				magmaxe = 0.1

		print('Time, Maximum magnitudee, Error, Mabs:')
		print(round(float(tmagmax),2),round(float(magmax),3),round(float(magmaxe),3))

	else:
		print ('First photometric point to far from explosion date')
		tmagmax=999.9
		magmax=999.9
		magmaxe=999.9
	########## SAVE #############
	Tmax = '%.2f' %tmagmax #T max
	mmax = '%.2f' %magmax #mmax
	err_mmax = '%.2f' %magmaxe #err_max

	#if input('Save?(y,|n|) ')=='y':
	tabla.write(name_SN.ljust(10,' ')+Tmax.ljust(7,' ')+mmax.ljust(7,' ')+err_mmax.ljust(7,' '))
	tabla.close()
	p.close()
