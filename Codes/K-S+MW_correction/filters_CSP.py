### Function de transmission du filtre g ###########
trans_B_CSP=np.loadtxt('Filters/CSP/B_swope.txt')
lambda_B_CSP=trans_B_CSP[:,0]
s_B_CSP=trans_B_CSP[:,1]
### Function de transmission du filtre r ###########
trans_r_CSP=np.loadtxt('Filters/CSP/r_swope.txt')
lambda_r_CSP=trans_r_CSP[:,0]
s_r_CSP=trans_r_CSP[:,1]
### Function de transmission du filtre i atm ###########
trans_i_CSP=np.loadtxt('Filters/CSP/i_swope.txt')
lambda_i_CSP=trans_i_CSP[:,0]
s_i_CSP=trans_i_CSP[:,1]
### Function de transmission du filtre i atm ###########
trans_V_CSP=np.loadtxt('Filters/CSP/V_swope.txt')
lambda_V_CSP=trans_V_CSP[:,0]
s_V_CSP=trans_V_CSP[:,1]

F_B_func_CSP=interpolate.interp1d(lambda_B_CSP,s_B_CSP) #interpolation Filtre g
F_V_func_CSP=interpolate.interp1d(lambda_V_CSP,s_V_CSP) #interpolation Filtre r
F_i_func_CSP=interpolate.interp1d(lambda_i_CSP,s_i_CSP) #interpolation Filtre i
F_r_func_CSP=interpolate.interp1d(lambda_r_CSP,s_r_CSP) #interpolation Filtre Y

