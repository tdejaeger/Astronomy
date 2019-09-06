import re # use regular patters
import sys # system commands
import string as string # string functions4
import math
import numpy as np # numerical tools
from scipy import *
from pylab import *
import os
from scipy import integrate
from scipy import interpolate
import itertools
import math as maths
from matplotlib.pyplot import figure, show, rcParams
import shutil

plt.rcParams['legend.numpoints']=1
plt.rcParams['xtick.major.size'] = 11
plt.rcParams['xtick.minor.size'] = 5
plt.rcParams['ytick.major.size'] = 11
plt.rcParams['ytick.minor.size'] = 5
plt.rcParams['xtick.minor.visible']=True #See minor tick
plt.rcParams['text.usetex']=True #use Latex
plt.rcParams['axes.linewidth']=2 #width axes
plt.rcParams['axes.labelsize']=25 #
plt.rcParams['ytick.labelsize']=22 #fontsize of tick labels
plt.rcParams['xtick.labelsize']=22 #fontsize of tick labels
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


exec(open("filters_CSP.py").read())
#Filter KAIT
exec(open("filters_KAIT.py").read())

fig,(ax1,ax2)=subplots(2,figsize=(6, 6),sharex=True,sharey=False ,facecolor='w', edgecolor='k')
fig.subplots_adjust(hspace=0)

ax1.plot(lambda_B_CSP,s_B_CSP/max(s_B_CSP),color='brown',ls=':',lw=2,label=r'\textbf{CSP}')
ax1.plot(lambda_r_CSP,s_r_CSP/max(s_r_CSP),color='brown',ls=':',lw=2)
ax2.plot(lambda_V_CSP,s_V_CSP/max(s_V_CSP),color='brown',ls=':',lw=2)
ax2.plot(lambda_i_CSP,s_i_CSP/max(s_i_CSP),color='brown',ls=':',lw=2)

ax1.plot(lambda_B_kait2,s_B_kait2/max(s_B_kait2),color='m',ls=(0, (3, 5, 1, 5, 1, 5)),lw=2,label=r'\textbf{KAIT 2}')
ax1.plot(lambda_R_kait2,s_R_kait2/max(s_V_kait2),color='m',ls=(0, (3, 5, 1, 5, 1, 5)),lw=2)
ax2.plot(lambda_V_kait2,s_V_kait2/max(s_R_kait2),color='m',ls=(0, (3, 5, 1, 5, 1, 5)),lw=2)
ax2.plot(lambda_I_kait2,s_I_kait2/max(s_I_kait2),color='m',ls=(0, (3, 5, 1, 5, 1, 5)),lw=2)

ax1.plot(lambda_B_kait3,s_B_kait3/max(s_B_kait3),color='darkgreen',ls='dashdot',lw=2,label=r'\textbf{KAIT 3}')
ax1.plot(lambda_R_kait3,s_R_kait3/max(s_V_kait3),color='darkgreen',ls='dashdot',lw=2)
ax2.plot(lambda_V_kait3,s_V_kait3/max(s_R_kait3),color='darkgreen',ls='dashdot',lw=2)
ax2.plot(lambda_I_kait3,s_I_kait3/max(s_I_kait3),color='darkgreen',ls='dashdot',lw=2)

ax1.plot(lambda_B_kait4,s_B_kait4/max(s_B_kait4),'k',alpha=0.7,ls='solid',lw=2,label=r'\textbf{KAIT 4}')
ax1.plot(lambda_R_kait4,s_R_kait4/max(s_V_kait4),'k',alpha=0.7,ls='solid',lw=2)
ax2.plot(lambda_V_kait4,s_V_kait4/max(s_R_kait4),'k',alpha=0.7,ls='solid',lw=2)
ax2.plot(lambda_I_kait4,s_I_kait4/max(s_I_kait4),'k',alpha=0.7,ls='solid',lw=2)

ax1.plot(lambda_B_nickel1,s_B_nickel1/max(s_B_nickel1),'r',ls='dashed',lw=2,label=r'\textbf{Nickel 1}')
ax1.plot(lambda_R_nickel1,s_R_nickel1/max(s_V_nickel1),'r',ls='dashed',lw=2)
ax2.plot(lambda_V_nickel1,s_V_nickel1/max(s_R_nickel1),'r',ls='dashed',lw=2)
ax2.plot(lambda_I_nickel1,s_I_nickel1/max(s_I_nickel1),'r',ls='dashed',lw=2)

ax1.plot(lambda_B_nickel2,s_B_nickel2/max(s_B_nickel2),'b',ls='dotted',lw=2,label=r'\textbf{Nickel 2}')
ax1.plot(lambda_R_nickel2,s_R_nickel2/max(s_V_nickel2),'b',ls='dotted',lw=2)
ax2.plot(lambda_V_nickel2,s_V_nickel2/max(s_R_nickel2),'b',ls='dotted',lw=2)
ax2.plot(lambda_I_nickel2,s_I_nickel2/max(s_I_nickel2),'b',ls='dotted',lw=2)

ax1.set_xlim([3000,12500])
ax1.set_ylim([0,1.2])
ax2.set_ylim([0,1.2])
ax2.set_xlabel(r'\textbf{Wavelength [\AA]}')
ax1.legend(loc=0,markerscale=0.5,prop={'size':13},ncol=1)
ax1.text(0.15, 0.15,r'\textbf{B}',horizontalalignment='center',verticalalignment='center',color='k',fontsize=18,transform = ax1.transAxes)
ax1.text(0.37, 0.15,r'\textbf{R/r}',horizontalalignment='center',verticalalignment='center',color='k',fontsize=18,transform = ax1.transAxes)
ax2.text(0.25, 0.15,r'\textbf{V}',horizontalalignment='center',verticalalignment='center',color='k',fontsize=18,transform = ax2.transAxes)
ax2.text(0.52, 0.15,r'\textbf{I/i}',horizontalalignment='center',verticalalignment='center',color='k',fontsize=18,transform = ax2.transAxes)
ax1.text(-0.15, 0.00,r'\textbf{Normalised Transmission}',horizontalalignment='center',rotation=90,verticalalignment='center',color='k',fontsize=22,transform = ax1.transAxes)

plt.savefig('Figures/transmission_curve.png',bbox_inches='tight')
close('all')
