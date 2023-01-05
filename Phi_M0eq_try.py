from PYmodule import *
from PYmodule.l_intg import *
# same as Phi_M0eq, numerically calculate dlnl/dlogM1

t1 = time.time()

N_Mh = 3 # 3 base halos: 1e11, 1e12, 1e13
Ntr = 10000
eta = 0.3
z = int(6)
tz = t_from_z(z)
tz = 300*Myr

d_fit = pow(10., -.5)
l_cut = .1
a = .1
t_life = 15
tz = 500*Myr

x0 = lambda_0/l_cut
I_toinf = integral_toinf(a,x0)

print('t_life, d_fit, l_cut, a,  f_seed, x0, logM0 = ', 
t_life,', ',d_fit,', ', l_cut,', ',a,', ', f_seed,', ', x0,', ', logM0,', ')

t_life *= Myr
T = Ts[0][0]
T['t_col'] = 20*Myr*np.ones(len(T))
T['Mstar0'] = 1000
f_bsm = 1.
n_base = n_base[0]

# print('mean Dt:',np.mean((tz-T['t_col']))/Myr)

i = 0
Chi2_min = 1e10; find_min = False

## --------- Mass Function ---------
dn_MBH = np.zeros(N_mf)
Nt = np.max((tz-T['t_col'])//t_life)
Nmax = Nt
dP_MBH = np.zeros(N_mf)
# t_tot = np.zeros(len(T))
DT = 0
M0s = np.logspace(2,12,num=10000)
dlog10M0 = np.log10(M0s[1]/M0s[0])

while Nt>=0:
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
        DT+=dt_seed
    else:
        dP_seed = np.zeros(N_mf)
    # prev BHMF
    DT+= t_life
    dP_MBH_prev = np.exp(np.interp(np.log(M0s),np.log(M_BH),np.log(dP_MBH))) # M0 grow, consv_ratio=0
    dP_MBH_prev = np.interp(np.log(M0s),np.log(M_BH),dP_MBH)
    # print(dP_MBH_prev)
    eps = 1e-5
    for iM1 in range(N_mf):
        # kernelS_MBH_M(M1, M0, dt, f_duty, l_cut, d_fit, logM_0=logM0):
        l1 = kernelS_MBH_M(M_BH[iM1],         M0s,t_life,1.,l_cut,d_fit)
        l2 = kernelS_MBH_M(M_BH[iM1]*(1.+eps),M0s,t_life,1.,l_cut,d_fit)
        dlnldlogM1 = np.log(l2/l1)/np.log10(1.+eps)
        # dlnldlogM1 = np.log(10.) * M_BH[iM1]/1e8/(l1*l_cut) * .5/(t_life/(0.1*t_Edd)) * (1.e8/M_BH[iM1] + pow(M_BH[iM1]/1.e8,d_fit-1.))

        klbd = dlnldlogM1 * pow(l1,a)*np.exp(-l1)/I_toinf
        # dP = (integral(a,l2,x0)-integral(a,l1,x0))/I_toinf
        # klbd = dP/np.log10(1.+eps)
        # klbd[np.logical_and(l1<x0, l2<x0)] = 0
        klbd[l1<x0] = 0
        # print( 'sum of Plambda',np.nansum(klbd/dlnldlogM1),dlog10M )
        dP_MBH[iM1] = np.nansum(klbd*dP_MBH_prev*dlog10M0) + dP_seed[iM1]/dlog10M

    # print('each cycle: consv_ratio =',np.nansum(dP_MBH*dlog10M))

    Nt -= 1
print('DT=',np.mean((DT-t_life)/Myr))
dn_MBH = dP_MBH*n_base*f_bsm*f_seed
dn_MBH = dP_MBH*n_base*f_bsm*f_seed*dlog10M

consv_ratio = np.nansum(dn_MBH)/(n_base*f_seed)
print('in Phi_easy: MF conserved fraction=%.10f'%consv_ratio)
# if consv_ratio<.9:
#     print('conserved fraction=%.10f'%consv_ratio)

T = Table(
    [M_BH, dn_MBH/dlog10M, MF(M_BH)],
    names=('M_BH','Phi','W10_MF')
)

MFname = z6datapre+'Phi_M0eq_tryMF'
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

lbin = np.arange(-3,3,0.1)
x = np.logspace(lbin[0]-np.log10(l_cut),lbin[-1]-np.log10(l_cut),num=len(lbin)) # for Pana
Pana = integral(a,x,x0)/I_toinf
with open(z6datapre+"Phi_M0eq_tryERDFz6.txt",'w') as f:
    np.savetxt(f, np.array([lbin[:-1],Pana[1:]-Pana[:-1]]).transpose(), fmt='%10.3e')

print('time=',time.time()-t1)
exit(0)

# # --------- Luminosity Function ---------
z_mesh = kernelS_M1450_mesh(bin_edg, M_BH, l_cut)
z_mesh[z_mesh<x0] = x0
Ps = integral(a,z_mesh,x0)/I_toinf
dPhi_mesh = np.nansum((Ps[:-1,:]-Ps[1:,:])*dn_MBH,axis=1)

# using abin_lf: test consv
# consv_ratio = np.nansum(dPhi_mesh)/(n_base*f_seed)
# print('LF consv_ratio:',consv_ratio)

Phi = dPhi_mesh/bin_wid

Phi *= 1e9
Phi_DO = Phi/corr_U14D20(bin_cen)
Chi2 = np.nansum(pow( (np.log(Phi_DO) - np.log(Phi_obs))/np.log(Phi_err), 2))/(len(Phi_obs)-1)
off_L = np.nanmax(abs( (np.log(Phi_DO) - np.log(Phi_obs))/np.log(Phi_err)))

T = Table(
    [bin_cen,Phi_obs,Phi_DO,Phi],
    names=('bin_cen','Phi_obs','Phi_DO','Phi')
)

LFname = z6datapre+'Phi_M0eq_tryLF'
ascii.write(T, LFname,
            formats={'bin_cen':'6.2f','Phi_obs':'4.2e','Phi_DO':'4.2e','Phi':'4.2e','Chi2':'4.2e'},
            overwrite=True)

if np.nanmin([Chi2, Chi2_min]) == Chi2:
    find_min = True
    Chi2_min = Chi2
    T_min = T
    LFname_min = LFname

print('Chi2_M=',Chi2_M,'Chi2_L=',Chi2)
print('log_prob=',-.5*(Chi2_min*(len(Phi_obs)-1)+Chi2_M))
#, LFname_min,'Chi2_min',Chi2_min, 'Chi2_M',Chi2_M, 'off_L',off_L, 'off_M',off_M)

# print('time=',time.time()-t1)