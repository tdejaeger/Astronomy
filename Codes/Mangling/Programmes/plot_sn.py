n_subplot=num.size(bands)
cols_subplot=2
# Compute Rows required
rows_subplot = n_subplot // cols_subplot 
rows_subplot += n_subplot % cols_subplot
# Create a Position index
pos_subplot = range(1,n_subplot + 1)

# Create main figure

f= plt.figure(1,figsize=(8, 8), facecolor='w', edgecolor='k')
f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.0, hspace=0.0)
for i in range(n_subplot):
	# add every single subplot to the figure with a for loop
	ax = f.add_subplot(rows_subplot,cols_subplot,pos_subplot[i])
	ax.tick_params(axis='both', which='major', labelsize=15,direction='out')
	ax.minorticks_on()
	ax.tick_params(axis='both', which='minor', labelsize=15,direction='out')
	ax.tick_params(which='minor', length=4,width=2, color='k')
	ax.tick_params(which='major',width=2, color='k')

	if pos_subplot[i] % 2 == 0:
		ax.yaxis.tick_right()
	ax.text(0.8, 0.80,'%s'%bands[i],horizontalalignment='center',verticalalignment='center',fontsize=15,transform = ax.transAxes)
	ax.errorbar(MJD['%s'%bands[i]]-MJD['%s'%bands[i]][0],mags['%s'%bands[i]],yerr=emags['%s'%bands[i]],marker='o',color='b',linestyle='None')
	ax.invert_yaxis()

f.suptitle('%s'%sn_name, fontsize=16)
f.text(0.5, 0.01, 'Date since first epoch [days]', ha='center',fontsize=15)
f.text(0.01, 0.5, 'Magnitude [mags]', va='center', rotation='vertical',fontsize=15)
f.savefig('Figures/photometry_%s.png'%sn_name)
show()

