from PYmodule import *
from PYmodule.l_intg import *

t1 = time.time()

N_Mh = 3 # 3 base halos: 1e11, 1e12, 1e13
Ntr = 10000
eta = 0.3
z = int(6)
tz = t_from_z(z)
alpha = 1.

# LF bins same w/ Matsu18
N_lf = len(bin_cen)

# # initial: lambda_0=0.01, logM0 = 8.
# t_life, d_fit, l_cut, a = 20, .01, 1., 0.1 # f_seed = .01, log_prob= -9.89
# t_life, d_fit, l_cut, a = 25, .01, 1.2, -0.2 # f_seed = .1, log_prob= -36.82
# t_life, d_fit, l_cut, a = 30, .01, 1., -.2 # f_seed = 1., log_prob= -13.88
# # # best:
# # t_life, d_fit, l_cut, a = 19.8, 1.2e-3, 1.1557, -1.8e-01 # f_seed = 1.

# new_nbase initial: lambda_0=0.01, logM0 = 8.
t_life, d_fit, l_cut, a = 30, .01, 1., 0.1 # f_seed = .01, log_prob= -8.14
t_life, d_fit, l_cut, a = 35, .01, 1.2, -0.2 # f_seed = .1, log_prob= -14.76
t_life, d_fit, l_cut, a = 40, .01, .9, -.2 # f_seed = 1., log_prob= -10.33

# # bests
# t_life, d_fit, l_cut, a = 21.9, .1,  0.87,  0.20 # f_seed = .01, log_prob=-3.16

x0 = lambda_0/l_cut
I_toinf = integral_toinf(a,x0)

print('t_life, d_fit, l_cut, a,  f_seed, x0, logM0 = ', 
t_life,', ',d_fit,', ', l_cut,', ',a,', ', f_seed,', ', x0,', ', logM0,', ')

t_life *= Myr


# print('mean Dt:',np.mean((tz-T['t_col']))/Myr)

Chi2_min = 1e10; find_min = False

N_Mh = 3
## --------- Mass Function ---------
dn_MBH = np.zeros(N_mf)
for iM in range(N_Mh):
    for i_bsm in range(Nbsm):
        T = Ts[iM][i_bsm]
        Nt = np.max((tz-T['t_col'])//t_life)
        dP_MBH = np.zeros(N_mf)
        # t_tot = np.zeros(len(T))
        while Nt>=0:
            t_point = tz - Nt*t_life
            T_seed = T[np.logical_and(t_point-t_life<=T['t_col'],T['t_col']<t_point)]
            dt_seed = t_point - T_seed['t_col']
            dP_MBH_prev = dP_MBH.copy()
            # new seeds (using 2d meshgrids)
            if len(T_seed):
                # z_mesh = kernelS_MBHmesh(abin_mf, T_seed['Mstar0'], dt_seed, l_cut)
                z_mesh = kernelS_MBH_M_mesh(abin_mf, T_seed['Mstar0'], dt_seed, 1., l_cut, d_fit)
                z_mesh[z_mesh<x0] = x0
                Ps = integral(a,z_mesh,x0)/I_toinf
                dP_seed = Ps[1:,:] - Ps[:-1,:]
                dP_seed = np.nansum(dP_seed, axis=1)/len(T)
            else:
                dP_seed = 0.
            # prev BHMF
            z_mesh = kernelS_MBH_M_mesh(M_BH, abin_mf, t_life, 1., l_cut, d_fit)
            z_mesh[z_mesh<x0] = x0
            Ps = integral(a,z_mesh,x0)/I_toinf
            dP_MBH = np.nansum( (Ps[:,:-1]-Ps[:,1:])*dP_MBH_prev, axis=1) + dP_seed
            Nt -= 1
        dn_MBH += dP_MBH*n_base[iM]*f_bsm[i_bsm]*f_seed

consv_ratio = np.nansum(dn_MBH)/(np.sum(n_base[:N_Mh]))
print('in Phi: MF conserved fraction=%.10f'%consv_ratio)
# if consv_ratio<.9:
#     print('conserved fraction=%.10f'%consv_ratio)

T = Table(
    [M_BH, dn_MBH/dlog10M, MF(M_BH)],
    names=('M_BH','Phi','W10_MF')
)

MFname = z6datapre+'MFIMF_3n2f'
ascii.write( Table([np.log10(T['M_BH']), T['Phi'], T['W10_MF']],
            names=['M_BH','Phi','W10_MF']),
            MFname,
            formats={'M_BH':'4.2e','Phi':'4.2e','W10_MF':'4.2e'},
            overwrite=True)
# exit(0)
# T  = T[np.logical_and(True,T['M_BH']<2e10)] # select M_BH range
index = np.logical_and(T['M_BH']>1e7, T['M_BH']<1e10)
Chi2_M = np.sum(pow(np.log(T['Phi']/T['W10_MF'])[index], 2))/(np.sum(index)-1)
off_M = np.max(abs(np.log(T['Phi']/T['W10_MF'])[index]))

# 10 N_M in 1e7-1e10 range, plus 12 N_L
index = np.where(np.logical_and(1e7<M_BH,M_BH<1e10))
xs = M_BH[index][::len(index[0])//10]
ys = np.log10( MF(xs)  ) # Willott 2010 30 points as data
y_model = np.log10( (dn_MBH/dlog10M)[index][::len(index[0])//10] )
y_err = pow(np.log10(xs)-8.5,2)/3. + .2 # from 0.2 to 0.95
Chi2_M =  np.sum( pow((ys - y_model)/y_err, 2))

# # --------- Luminosity Function ---------
Phi = np.zeros(N_lf)

T['dn_MBH'] = T['Phi']*dlog10M
z_mesh = kernelS_M1450_mesh(bin_edg, M_BH, l_cut)
z_mesh[z_mesh<x0] = x0
Ps = integral(a,z_mesh,x0)/I_toinf
dPhi_mesh = np.nansum((Ps[:-1,:]-Ps[1:,:])*dn_MBH,axis=1)

Phi = dPhi_mesh/bin_wid

Phi *= 1e9
Phi_DO = Phi/corr_U14D20(bin_cen)
Chi2 = np.nansum(pow( (np.log(Phi_DO) - np.log(Phi_obs))/np.log(Phi_err), 2))/(len(Phi_obs)-1)
off_L = np.nanmax(abs( (np.log(Phi_DO) - np.log(Phi_obs))/np.log(Phi_err)))

T = Table(
    [bin_cen,Phi_obs,Phi_DO,Phi,Chi2*np.ones(N_lf)],
    names=('bin_cen','Phi_obs','Phi_DO','Phi','Chi2')
)

LFname = z6datapre+'LFIMF_3n2f'
ascii.write(T, LFname,
            formats={'bin_cen':'6.2f','Phi_obs':'4.2e','Phi_DO':'4.2e','Phi':'4.2e','Chi2':'4.2e'},
            overwrite=True)

if np.nanmin([Chi2, Chi2_min]) == Chi2:
    find_min = True
    Chi2_min = Chi2
    T_min = T
    LFname_min = LFname

print('Chi2_M=',Chi2_M,'Chi2_L=',Chi2)
print('log_like=',-.5*(Chi2_min*(len(Phi_obs)-1)+Chi2_M))
#, LFname_min,'Chi2_min',Chi2_min, 'Chi2_M',Chi2_M, 'off_L',off_L, 'off_M',off_M)

# print('time=',time.time()-t1)