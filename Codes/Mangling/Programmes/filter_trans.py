'''
Get LCOGT filter light-curves

'''
from scipy import interpolate
lambda_eff_tot_LCOGT=[3540,4361,4770,5448,6215,7545,8700] #uBgVriz
#uBgVri
### Filter u ###########
trans_u=np.loadtxt('Filters/u.txt')
lambda_u=np.array(trans_u[:,0]*10)[np.array(trans_u[:,1])>0.01]
s_u=np.array(trans_u[:,1])[np.array(trans_u[:,1])>0.01]
### Filter g ###########
trans_g=np.loadtxt('Filters/g.txt')
lambda_g=np.array(trans_g[:,0]*10)[np.array(trans_g[:,1])>0.01]
s_g=np.array(trans_g[:,1])[np.array(trans_g[:,1])>0.01]
### Filter r ###########
trans_r=np.loadtxt('Filters/r.txt')
lambda_r=np.array(trans_r[:,0]*10)[np.array(trans_r[:,1])>0.01]
s_r=np.array(trans_r[:,1])[np.array(trans_r[:,1])>0.01]
### Fitler i ###########
trans_i=np.loadtxt('Filters/i.txt')
lambda_i=np.array(trans_i[:,0]*10)[np.array(trans_i[:,1])>0.01]
s_i=np.array(trans_i[:,1])[np.array(trans_i[:,1])>0.01]
### Fitler z ###########
trans_z=np.loadtxt('Filters/z.txt')
lambda_z=np.array(trans_z[:,0]*10)[np.array(trans_z[:,1])>0.01]
s_z=np.array(trans_z[:,1])[np.array(trans_z[:,1])>0.01]
### Filter B ###########
trans_B=np.loadtxt('Filters/B.txt')
lambda_B=np.array(trans_B[:,0]*10)[np.array(trans_B[:,1])>0.01]
s_B=np.array(trans_B[:,1])[np.array(trans_B[:,1])>0.01]
### Filter V ###########
trans_V=np.loadtxt('Filters/V.txt')
lambda_V=np.array(trans_V[:,0]*10)[np.array(trans_V[:,1])>0.01]
s_V=np.array(trans_V[:,1])[np.array(trans_V[:,1])>0.01]

filt_B_func=interpolate.interp1d(lambda_B,s_B) #interpolation Filtre
filt_V_func=interpolate.interp1d(lambda_V,s_V) #interpolation Filtre V
filt_u_func=interpolate.interp1d(lambda_u,s_u) #interpolation Filtre u
filt_g_func=interpolate.interp1d(lambda_g,s_g) #interpolation Filtre g
filt_r_func=interpolate.interp1d(lambda_r,s_r) #interpolation Filtre r
filt_i_func=interpolate.interp1d(lambda_i,s_i) #interpolation Filtre i
filt_z_func=interpolate.interp1d(lambda_z,s_z) #interpolation Filtre z

filt_tot_func=[filt_u_func,filt_B_func,filt_g_func,filt_V_func,filt_r_func,filt_i_func,filt_z_func]
lambda_filt_tot=[lambda_u,lambda_B,lambda_g,lambda_V,lambda_r,lambda_i,lambda_z]
bands_LCOGT=['u_2.5m','B','g_2.5m','V','r_2.5m','i_2.5m','z_2.5m']
dict_ZP=dict()
#LCOGT
for i in range(size(bands_LCOGT)):
	dict_ZP["%s"%bands_LCOGT[i]] = [lambda_eff_tot_LCOGT[i],filt_tot_func[i],lambda_filt_tot[i]]

'''
Get ATLAS filter light-curves and ZP

'''
from scipy import interpolate
lambda_eff_tot_ATLAS=[4100,6880] #uBgVriz
#uBgVri
### Filter cyan ###########
trans_c=np.loadtxt('Filters/cyan.dat')
lambda_c=trans_c[:,0]*10
s_c=trans_c[:,1]
### Filter orange ###########
trans_o=np.loadtxt('Filters/orange.dat')
lambda_o=trans_o[:,0]*10
s_o=trans_o[:,1]

filt_c_func=interpolate.interp1d(lambda_c,s_c) #interpolation Filtre c
filt_o_func=interpolate.interp1d(lambda_o,s_o) #interpolation Filtre o

filt_tot_func=[filt_c_func,filt_o_func]
lambda_filt_tot=[lambda_c,lambda_o]
bands_ATLAS=['cyan','orange']
#ATLAS
for i in range(size(bands_ATLAS)):
	dict_ZP["%s"%bands_ATLAS[i]] = [lambda_eff_tot_ATLAS[i],filt_tot_func[i],lambda_filt_tot[i]]
