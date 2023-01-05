from PYmodule import *
from PYmodule.l_intg import *

# pre-determined in __init__.py
# eta = 0.3
# alpha = 1.

z = int(6)
T = Ts[0][0]
f_bsm = 1.
n_base = n_base[0]


def lnlike(theta):
    t_life, l_cut, a = theta
    d_fit = 0
    t_life *= Myr
    x0 = lambda_0/l_cut
    I_toinf = integral_toinf(a,x0)
## --------- Mass Function ---------
    tz = t_from_z(z)
    dn_MBH = np.zeros(N_mf)
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
        # z_mesh = kernelS_MBHmesh(M_BH, abin_mf, t_life, l_cut)
        z_mesh = kernelS_MBH_M_mesh(M_BH, abin_mf, t_life, 1., l_cut, d_fit)
        z_mesh[z_mesh<x0] = x0
        Ps = integral(a,z_mesh,x0)/I_toinf
        dP_MBH = np.nansum( (Ps[:,:-1]-Ps[:,1:])*dP_MBH_prev, axis=1) + dP_seed

        Nt -= 1
    dn_MBH = dP_MBH*n_base*f_bsm*f_seed

    consv_ratio = np.nansum(dn_MBH)/(n_base*f_seed)
    if abs(consv_ratio-1)>.5:
        print('theta: ',theta,'x0',x0)
        print('consv_ratio: ',consv_ratio)
        return -np.inf

    # 10 N_M in 1e7-1e10 range, plus 12 N_L
    index = np.where(np.logical_and(1e7<M_BH,M_BH<1e10))
    xs = M_BH[index][::len(index[0])//10]
    ys = np.log10( MF(xs)  ) # Willott 2010 as data
    y_model = np.log10( (dn_MBH/dlog10M)[index][::len(index[0])//10] )
    y_err = pow(np.log10(xs)-8.5,2)/3. + .2 # from 0.2 to 0.95
    Chi2_M =  np.sum( pow((ys - y_model)/y_err, 2))
    if not np.isfinite(Chi2_M):
        print('theta=',theta)
        print('inf or nan? Chi2_M=',Chi2_M)
        return -np.inf

# # --------- Luminosity Function ---------
    z_mesh = kernelS_M1450_mesh(bin_edg, M_BH, l_cut)
    z_mesh[z_mesh<x0] = x0
    Ps = integral(a,z_mesh,x0)/I_toinf
    dPhi_mesh = np.nansum((Ps[:-1,:]-Ps[1:,:])*dn_MBH,axis=1)
    Phi = dPhi_mesh/bin_wid

    Phi *= 1e9
    Phi_DO = Phi/corr_U14D20(bin_cen)
    # Chi2 = np.nansum(pow( (np.log(Phi_DO) - np.log(Phi_obs))/np.log(Phi_err), 2))/(len(Phi_obs)-1)
    # print('Chi2=%.2e'%Chi2)
    ys = np.log10(Phi_obs)
    y_model = np.log10(Phi_DO)
    y_err = np.log10(Phi_err)
    Chi2_L = np.sum( pow((ys - y_model)/y_err, 2))
    if not np.isfinite(Chi2_L):
        print('theta=',theta)
        print('inf or nan? Chi2_L=',Chi2_L)
        return -np.inf
    return -.5*(Chi2_M + Chi2_L)


# 3prange1: 1e1<t_life<200. and l_cut>0:
# 3prange2: 1e1<t_life<200. and l_cut>0.1:

def lnprior(theta):
    t_life, l_cut, a = theta
    if 1e1<t_life<200. and l_cut>0.1:
        return 0.0 - 0.5*((l_cut-l_mean)/sigma_l)**2 - 0.5*((a-a_mean)/sigma_a)**2
    else:
        return -np.inf

def lnprobab(theta):
    lp = lnprior(theta)
    # print('theta',theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta)
