from PYmodule import *
from PYmodule.l_intg import *

t1 = time.time()
prex = '../data_eq/'

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
t_life, d_fit, l_cut, a = 30, .01, 1., 0.1 # f_seed = .01, log_prob= -9.11
t_life, d_fit, l_cut, a = 35, .01, 1.2, -0.2 # f_seed = .1, log_prob= -16.37
t_life, d_fit, l_cut, a = 40, .01, .9, -.2 # f_seed = 1., log_prob= -11.45

t_life, d_fit, l_cut, a = 21.4, 0., .89, .15 # f_seed = 0.1, M1M0_e
t_life, d_fit, l_cut, a = 21.4, 0.001, .89, .15 # f_seed = 0.1, M1M0_d

# after more abin_mf bins, after calibration w/ direct sampling; 
# easycali initial: (& 3 points)
t_life, logd_fit, l_cut, a = 21.8, -1, .88, .19 # f_seed = 0.01
t_life, logd_fit, l_cut, a = 21.4, -3, .89, .15 # f_seed = 0.1
t_life, logd_fit, l_cut, a = 22.2, -2.98, .99, -.04 # f_seed = 1
# easycali best:
t_life, logd_fit, l_cut, a = 19.9, -1.08, .87, .17 # f_seed = 0.01, easycali
# t_life, logd_fit, l_cut, a = 19.6, -2.96, .87, .12 # f_seed = 0.1


x0 = lambda_0/l_cut
I_toinf = integral_toinf(a,x0)

d_fit = pow(10.,logd_fit)
print('t_life, d_fit, l_cut, a,  f_seed, x0, logM0 = ', 
t_life,', ',d_fit,', ', l_cut,', ',a,', ', f_seed,', ', x0,', ', logM0,', ')

print('mf len:',len(abin_mf))
t_life *= Myr
T = Ts[0][0]
# T['t_col'] = 20*Myr; T['Mstar0'] = 1e3
f_bsm = 1.
n_base = n_base[0]

# print('mean Dt:',np.mean((tz-T['t_col']))/Myr)


## --------- Mass Function ---------
dn_MBH = np.zeros(N_mf)
Nt = np.max((tz-T['t_col'])//t_life)
Nmax = Nt
dP_MBH = np.zeros(N_mf)
# t_tot = np.zeros(len(T))
while Nt>=0:
# for i__ in range(20):
    t_point = tz - Nt*t_life
    T_seed = T[np.logical_and(t_point-t_life<T['t_col'],T['t_col']<=t_point)]
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
        print(np.nansum(dP_seed),'sum of dP_seed')
    else:
        dP_seed = np.zeros(N_mf)
        print('all seeds counted in')
    # prev BHMF
    # print(dP_MBH_prev)
    for iM1 in range(N_mf):
        # l1 = kernelS_MBH_M(M_BH[iM1],         M0s,t_life,1.,l_cut,d_fit)
        x1 = kernelS_MBH_M(M_BH[iM1], M_BH,t_life,1.,l_cut,d_fit)
        dlnldlogM1 = np.log(10.) * M_BH[iM1]/1e8/(x1*l_cut) * .5/(t_life/(0.1*t_Edd)) * (1.e8/M_BH[iM1] + pow(M_BH[iM1]/1.e8,d_fit-1.))
        dPdlnl = pow(x1,a)*np.exp(-x1)/I_toinf
        dPdlnl[x1<x0] = 0
        # print( 'sum of Plambda',np.nansum(klbd/dlnldlogM1),dlog10M )
        dP_MBH[iM1] = np.nansum(dPdlnl*dlnldlogM1*dP_MBH_prev*dlog10M) + dP_seed[iM1]/dlog10M

    print('Nt',Nt,' consv_ratio =',np.nansum(dP_MBH*dlog10M))
    # # renormalization
    # dP_MBH *= 1./np.nansum(dP_MBH*dlog10M)
    Nt -= 1

dn_MBH = dP_MBH*n_base*f_bsm*f_seed*dlog10M

consv_ratio = np.nansum(dP_MBH*dlog10M)
print('len of abin_mf:',len(abin_mf))
print('in Phi_easyeq: MF conserved fraction=%.10f'%consv_ratio)

T = Table(
    [M_BH, dn_MBH/dlog10M, MF(M_BH)],
    names=('M_BH','Phi','W10_MF')
)

MFname = prex+'Phi_easyeqMF'
ascii.write( Table([np.log10(T['M_BH']), T['Phi'], T['W10_MF']],
            names=['M_BH','Phi','W10_MF']),
            MFname,
            formats={'M_BH':'4.2e','Phi':'4.2e','W10_MF':'4.2e'},
            overwrite=True)

lbin = np.arange(-3,3,0.1)
x = np.logspace(lbin[0]-np.log10(l_cut),lbin[-1]-np.log10(l_cut),num=len(lbin)) # for Pana
Pana = integral(a,x,x0)/I_toinf
with open(prex+"Phi_easyeqERDFz6.txt",'w') as f:
    np.savetxt(f, np.array([lbin[:-1],Pana[1:]-Pana[:-1]]).transpose(), fmt='%10.3e')

print('time=',time.time()-t1)