from PYmodule import *
from PYmodule.l_intg import *

t1 = time.time()

N_Mh = 3 # 3 base halos: 1e11, 1e12, 1e13
Ntr = 10000
eta = 0.3
z = int(6)
corr = 'U'
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

# # 300 abin_mf bins, 3 points,  calibration w/ direct sampling; 
# # easycali initial: 
# t_life, logd_fit, l_cut, a = 21.8, -1, .88, .19; f_seed = 0.01
# t_life, logd_fit, l_cut, a = 21.4, -3, .89, .15; f_seed = 0.1
# t_life, logd_fit, l_cut, a = 22.2, -2.98, .99, -.04; f_seed = 1
# easycali best:
t_life, logd_fit, l_cut, a = 19.9, -1.08, .87, .17; f_seed = 0.01
t_life, logd_fit, l_cut, a = 19.6, -2.96, .87, .12; f_seed = 0.1
# t_life, logd_fit, l_cut, a = 26.1, -2.59, .88, -0.05; f_seed = 1

# # 3p init
# logd_fit, l_cut, a = -2, .87, .1; f_seed = 0.01
# logd_fit, l_cut, a = -2.96, .87, .01; f_seed = 0.1
# logd_fit, l_cut, a = -2.59, .88, -0.2; f_seed = 1

# t_life = 50

# # t_life << Dt, 找近似
# t_life, logd_fit, l_cut, a = 1., -3, 1, .1; f_seed = 0.01

# # M0 best
t_life, logd_fit, l_cut, a = 18.7555167,  -1.2574505,   0.87372563,  0.20389703; f_seed = 0.01
# t_life, logd_fit, l_cut, a = 20.07157851, -2.98140382,  0.89453609,  0.12195823; f_seed = 0.1
# t_life, logd_fit, l_cut, a = 23.12675104, -2.97342483,  0.95753445, -0.06535641; f_seed = 1


x0 = lambda_0/l_cut
I_toinf = integral_toinf(a,x0)
# print('I_toinf',I_toinf); exit(0)
d_fit = pow(10.,logd_fit)
print('t_life, d_fit, l_cut, a,  f_seed, x0, logM0 = ', 
t_life,', ',d_fit,', ', l_cut,', ',a,', ', f_seed,', ', x0,', ', logM0,', ')

t_life *= Myr
T = Ts[0][0]
f_bsm = 1.
n_base = n_base[0]
T_stay = T[:-int(f_seed*len(T))]
# add seed M_BHs 
hist_seed, bin_edges = np.histogram(T_stay['Mstar0'],bins=abin_mf,density=False)
dn_seed = hist_seed*n_base*f_bsm/len(T)

print(np.median(T['z_col']))
print('mean t_seed:',np.mean(T['t_col'])/Myr)
print('mean Dt:',np.mean((tz-T['t_col']))/Myr)
print('mean Dt:',np.mean((tz-t_from_z(35)))/Myr)

i = 0
Chi2_min = 1e10; find_min = False

## --------- Mass Function ---------
dn_MBH = np.zeros(N_mf)
Nt = np.max((tz-T['t_col'])//t_life)
Nmax = Nt
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
        dP_seed = np.zeros(N_mf)
    # prev BHMF
    z_mesh = kernelS_MBH_M_mesh(abin_mf, M_BH, t_life, 1., l_cut, d_fit)
    z_mesh[z_mesh<x0] = x0
    Ps = integral(a,z_mesh,x0)/I_toinf
    dP_MBH = np.nansum( (Ps[1:,:]-Ps[:-1,:])*dP_MBH_prev, axis=1) + dP_seed
    # print(np.nansum((Ps[1:,:]-Ps[:-1,:]),axis=0))
    Nt -= 1

dn_MBH = dP_MBH*n_base*f_bsm*f_seed

consv_ratio = np.nansum(dn_MBH)/(n_base*f_seed)
# print('in Phi_easy: MF conserved fraction=%.10f'%consv_ratio)
# if consv_ratio<.9:
#     print('conserved fraction=%.10f'%consv_ratio)

print('consv_ratio',consv_ratio)

PhiM_growth  =   dn_MBH/dlog10M

# # total BHMF
# # average every n values
# n = 5
# M_BH         =    ave_over( M_BH,n); 
# PhiM_growth  =    ave_over( PhiM_growth, n)
# PhiM_tot     =    ave_over( (dn_MBH+dn_seed)/dlog10M, n)
# dn_MBH       =    ave_over( dn_MBH,n); 

# T_tot = Table(
#     [M_BH, PhiM_growth, PhiM_tot, MF(M_BH,6)],
#     names=('M_BH','PhiM_growth','PhiM_tot','PhiM_obs')
# )
# MFname = './totBHMFz{:d}'.format(z)
# ascii.write(T_tot,
#             MFname,
#             formats={'M_BH':'4.2e','y_growth':'4.2e','y_tot':'4.2e','y_obs':'4.2e'},
#             overwrite=True)
# exit(0)

T = Table(
    [M_BH, PhiM_growth, MF(M_BH)],
    names=('M_BH','Phi','W10_MF')
)

MFname = z6datapre+'f{0:d}Phi_easyMF_z{1:d}'.format(int(abs(np.log10(f_seed))),z)
MFname = z6datapre+'f{:d}Phi_easyMF3pt{:d}'.format(int(abs(np.log10(f_seed))),int(t_life/Myr))
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
ys = np.log10( MF(xs)  ) # Willott 2010 as data
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
    [bin_cen,Phi_obs,Phi_DO,Phi,Chi2*np.ones(N_lf)],
    names=('bin_cen','Phi_obs','Phi_DO','Phi','Chi2')
)

LFname = z6datapre+'f{0:d}Phi_easyLFdot_z{1:d}'.format(int(abs(np.log10(f_seed))),z)
LFname = z6datapre+'f{:d}Phi_easyLF3pdot_t{:d}'.format(int(abs(np.log10(f_seed))),int(t_life/Myr))
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
    Phi_DO = Phi/corr_M14D20(M1450)

T = Table(
    [M1450,Phi,Phi/corr_U14D20(M1450),Phi/corr_M14D20(M1450)],
    names=('M1450','Phi','Phi_U','Phi_M')
)
LFname = z6datapre+'f{:d}Phi_easyLF3pt{:d}'.format(int(abs(np.log10(f_seed))),int(t_life/Myr))
# LFname = './f{:d}Phi_easyLF_z{:d}'.format(int(abs(np.log10(f_seed))),z)
ascii.write(T, LFname,
            formats={'M1450':'6.2f','Phi':'4.2e','Phi_U':'4.2e','Phi_M':'4.2e'},
            overwrite=True)

if np.nanmin([Chi2, Chi2_min]) == Chi2:
    find_min = True
    Chi2_min = Chi2
    T_min = T
    LFname_min = LFname

# print('Chi2_M=',Chi2_M,'Chi2_L=',Chi2)
print('ln_like=',-.5*(Chi2_min*(len(Phi_obs)-1)+Chi2_M))
#, LFname_min,'Chi2_min',Chi2_min, 'Chi2_M',Chi2_M, 'off_L',off_L, 'off_M',off_M)

print(z)
print('time=',time.time()-t1)