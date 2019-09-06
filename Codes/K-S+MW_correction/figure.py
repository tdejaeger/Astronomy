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



fig, ax1 = plt.subplots(figsize=(8,8), facecolor='w', edgecolor='k')  # ax1 is apparent magnitude
fig.subplots_adjust(hspace=0.0)

#CSP photometry
ax1.errorbar(phase_X_final_csp[0],X_AKS_final_csp[0],yerr=err_X_final_csp[0],marker='o',color='b',markersize=10,linestyle='None',label=r'\textbf{B CSP}')
ax1.errorbar(phase_X_final_csp[1],X_AKS_final_csp[1],yerr=err_X_final_csp[1],marker='s',color='c',markersize=10,linestyle='None',label=r'\textbf{VCSP}')
ax1.errorbar(phase_X_final_csp[2],X_AKS_final_csp[2]-1,yerr=err_X_final_csp[2],marker='<',color='r',markersize=10,linestyle='None',label=r'\textbf{r-1 CSP}')
ax1.errorbar(phase_X_final_csp[3],X_AKS_final_csp[3]-2,yerr=err_X_final_csp[3],marker='^',color='k',markersize=10,linestyle='None',label=r'\textbf{i-2 CSP}')

ax1.errorbar(phase_X_final_kait[0],X_AKS_final_kait[0],yerr=err_X_AKS_final_kait[0],marker='o',markersize=10,mfc='None', mec='b',ecolor='b',linestyle='None',label=r'\textbf{AKS B}')
ax1.errorbar(phase_X_final_kait[1],X_AKS_final_kait[1],yerr=err_X_AKS_final_kait[1],marker='s',markersize=10,mfc='None', mec='c',ecolor='c',linestyle='None',label=r'\textbf{AKS R}')
ax1.errorbar(phase_X_final_kait[2],X_AKS_final_kait[2]-1,yerr=err_X_AKS_final_kait[2],marker='<',markersize=10,mfc='None', mec='r',ecolor='r',linestyle='None',label=r'\textbf{AKS V}')
ax1.errorbar(phase_X_final_kait[3],X_AKS_final_kait[3]-2,yerr=err_X_AKS_final_kait[3],marker='^',markersize=10,mfc='None', mec='k',ecolor='k',linestyle='None',label=r'\textbf{AKS I}')

ax1.errorbar(phase_X_final_kait[0],X_final_kait[0],yerr=err_X_AKS_final_kait[0],marker='X',markersize=10,mfc='None', mec='b',ecolor='b',linestyle='None',label=r'\textbf{Nat B}')
ax1.errorbar(phase_X_final_kait[1],X_final_kait[1],yerr=err_X_AKS_final_kait[1],marker='X',markersize=10,mfc='None', mec='c',ecolor='c',linestyle='None',label=r'\textbf{Nat V}')
ax1.errorbar(phase_X_final_kait[2],X_final_kait[2]-1,yerr=err_X_AKS_final_kait[2],marker='X',markersize=10,mfc='None', mec='r',ecolor='r',linestyle='None',label=r'\textbf{Nat R}')
ax1.errorbar(phase_X_final_kait[3],X_final_kait[3]-2,yerr=err_X_AKS_final_kait[3],marker='X',markersize=10,mfc='None', mec='k',ecolor='k',linestyle='None',label=r'\textbf{Nat I}')



ax1.set_xlabel(r'\textbf{Epoch since explosion}')
ax1.set_ylabel(r'\textbf{mag}')
ax1.minorticks_on()
for tick in ax1.yaxis.get_major_ticks():
	tick.label.set_weight('bold')
for tick in ax1.xaxis.get_major_ticks():
	tick.label.set_weight('bold')
ax1.tick_params(axis='x', pad=10)
ax1.legend(loc=(-0.02, 1.01),title='',markerscale=0.8,ncol=6,prop={'size':10})

max_y=max(max(X_AKS_final_csp[0]),max(X_AKS_final_csp[1]),max(X_AKS_final_csp[2]-1),max(X_AKS_final_csp[3]-2))
min_y=min(min(X_AKS_final_csp[0]),min(X_AKS_final_csp[1]),min(X_AKS_final_csp[2]-1),min(X_AKS_final_csp[3]-2))
ax1.set_ylim([min_y-0.5,max_y+0.5])

max_x=max(max(phase_X_final_csp[0]),max(phase_X_final_csp[1]),max(phase_X_final_csp[2]),max(phase_X_final_csp[3]))
min_x=min(min(phase_X_final_csp[0]),min(phase_X_final_csp[1]),min(phase_X_final_csp[2]),min(phase_X_final_csp[3]))
ax1.set_xlim([min_x-10,max_x+30])
plt.gca().invert_yaxis()
#plt.show()
plt.savefig('Figures/%s.png'%SN_name)



