from PYmodule import *
from PYmodule.l_intg import *

t1 = time.time()

N_Mh = 3 # 3 base halos: 1e11, 1e12, 1e13
Ntr = 10000
eta = 0.3
z = int(6)
corr = 'U'
tz = t_from_z(z)

# # M0 best
t_life, logd_fit, l_cut, a = 18.7555167,  -1.2574505,   0.87372563,  0.20389703; f_seed = 0.01
# t_life, logd_fit, l_cut, a = 20.07157851, -2.98140382,  0.89453609,  0.12195823; f_seed = 0.1
# t_life, logd_fit, l_cut, a = 23.12675104, -2.97342483,  0.95753445, -0.06535641; f_seed = 1

t_life = 10

x0 = lambda_0/l_cut
I_toinf = integral_toinf(a,x0)
d_fit = pow(10.,logd_fit)
print('t_life, d_fit, l_cut, a,  f_seed, x0, logM0 = ', 
t_life,', ',d_fit,', ', l_cut,', ',a,', ', f_seed,', ', x0,', ', logM0,', ')

t_life *= Myr
T = Ts[0][0]
print('t_seed: mean={:.1f},median={:.1f}'.format(np.mean(T['t_col']/Myr), np.median(T['t_col'])/Myr))
print('M_seed: mean={:.1e},median={:.1e}'.format(np.mean(T['Mstar0']), np.median(T['Mstar0'])))
print('tz6 = %.1f'%(t_from_z(6)/Myr))
print('mean Dt:',np.mean((tz-T['t_col']))/Myr)

t_col = 80*Myr
Mstar0 = 1e3
T['t_col'] = t_col
T['Mstar0'] = Mstar0
# single M0 mass can grow with local maximum of lambda=a*l_cut to:
print('M1 expect=%.1e'%(Mstar0*np.exp(a*l_cut*(tz-t_col)/(45*Myr))))

f_bsm = 1.
n_base = n_base[0]


## --------- Mass Function ---------
dn_MBH = np.zeros(N_mf)
Nt = np.max((tz-T['t_col'])//t_life)
Nmax = Nt
dP_MBH = np.zeros(N_mf)
# t_tot = np.zeros(len(T))
DT = 0
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
        DT+=dt_seed
    else:
        dP_seed = np.zeros(N_mf)
    # prev BHMF
    DT+= t_life
    z_mesh = kernelS_MBH_M_mesh(abin_mf, M_BH, t_life, 1., l_cut, d_fit)
    z_mesh[z_mesh<x0] = x0
    Ps = integral(a,z_mesh,x0)/I_toinf
    dP_MBH = np.nansum( (Ps[1:,:]-Ps[:-1,:])*dP_MBH_prev, axis=1) + dP_seed
    Nt -= 1
print('DT=',np.mean((DT-t_life)/Myr))
dn_MBH = dP_MBH*n_base*f_bsm*f_seed

consv_ratio = np.nansum(dn_MBH)/(n_base*f_seed)
print('in Phi_M0: MF conserved fraction=%.10f'%consv_ratio)
# if consv_ratio<.9:
#     print('conserved fraction=%.10f'%consv_ratio)

T = Table(
    [M_BH, dn_MBH/dlog10M, MF(M_BH)],
    names=('M_BH','Phi','W10_MF')
)

MFname = z6datapre+'f{:d}Phi_M0{:.1e}MF_z{:d}'.format(int(abs(np.log10(f_seed))),Mstar0,z)
MFname = z6datapre+'f{:d}Phi_M0{:.1e}MFt{:d}'.format(int(abs(np.log10(f_seed))),Mstar0,int(t_life/Myr))
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
if corr=='M':
    Phi_DO = Phi/corr_M14D20(bin_cen)
Chi2 = np.nansum(pow( (np.log(Phi_DO) - np.log(Phi_obs))/np.log(Phi_err), 2))/(len(Phi_obs)-1)
off_L = np.nanmax(abs( (np.log(Phi_DO) - np.log(Phi_obs))/np.log(Phi_err)))

T = Table(
    [bin_cen,Phi_obs,Phi_DO,Phi],
    names=('bin_cen','Phi_obs','Phi_DO','Phi')
)

LFname = z6datapre+'f{:d}Phi_M0{:.1e}LFdot_z{:d}'.format(int(abs(np.log10(f_seed))),Mstar0,z)
LFname = z6datapre+'f{:d}Phi_M0{:.1e}LFdot_t{:d}'.format(int(abs(np.log10(f_seed))),Mstar0,int(t_life/Myr))
ascii.write(T, LFname,
            formats={'bin_cen':'6.2f','Phi_obs':'4.2e','Phi_DO':'4.2e','Phi':'4.2e','Chi2':'4.2e'},
            overwrite=True)


# for spread plot
z_mesh = kernelS_M1450_mesh(abin_lf, M_BH, l_cut)*l_cut
z_mesh[z_mesh<lambda_0] = lambda_0
Ps = integral(a,z_mesh/l_cut,lambda_0/l_cut)/I_toinf
dPhi_mesh = np.nansum((Ps[:-1,:]-Ps[1:,:])*dn_MBH,axis=1)
Phi = dPhi_mesh/dmag
Phi *= 1e9
Phi_DO = Phi/corr_U14D20(M1450)
if corr=='M':
    Phi_DO = Phi/corr_M14D20(bin_cen)

T = Table(
    [M1450,Phi_DO],
    names=('M1450','Phi_DO')
)
LFname = z6datapre+'f{:d}Phi_M0{:.1e}LFt{:d}'.format(int(abs(np.log10(f_seed))),Mstar0,int(t_life/Myr))
ascii.write(T, LFname,
            formats={'M1450':'6.2f','Phi_DO':'4.2e'},
            overwrite=True)

print('time=',time.time()-t1)
