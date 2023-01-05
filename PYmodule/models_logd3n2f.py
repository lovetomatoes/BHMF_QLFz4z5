from PYmodule import *
from PYmodule.l_intg import *

# pre-determined in __init__.py
# eta = 0.3
# alpha = 1.
# N_Mh = 1 now

z = int(6)

def model(theta, z = int(6), f_seed=f_seed, l_cut= l_cut, a=a):
    t_life, logd_fit, l_cut, a = theta
    d_fit = pow(10.,logd_fit)
    t_life = t_life * Myr # wli: not *=!!!!!! theta changed by this
    I_toinf = integral_toinf(a,lambda_0/l_cut)

## --------- Mass Function ---------
    tz = t_from_z(z)
    dn_MBH = np.zeros(N_mf)
    for iM in range(N_Mh):
        for i_bsm in range(Nbsm):
            T = Ts[iM][i_bsm]
            Nt = np.max((tz-T['t_col'])//t_life)
            dP_MBH = np.zeros(N_mf)
            while Nt>=0:
                t_point = tz - Nt*t_life
                T_seed = T[np.logical_and(t_point-t_life<=T['t_col'],T['t_col']<t_point)]
                dt_seed = t_point - T_seed['t_col']
                dP_MBH_prev = dP_MBH.copy()

                # new seeds (using 2d meshgrids)
                if len(T_seed):
                    # z_mesh = kernelS_MBHmesh(abin_mf, T_seed['Mstar0'], dt_seed, l_cut)
                    z_mesh = kernelS_MBH_M_mesh(abin_mf, T_seed['Mstar0'], dt_seed, 1., l_cut, d_fit)*l_cut
                    z_mesh[z_mesh<lambda_0] = lambda_0
                    Ps = integral(a,z_mesh/l_cut,lambda_0/l_cut)/I_toinf
                    dP_seed = Ps[1:,:] - Ps[:-1,:]
                    dP_seed = np.nansum(dP_seed, axis=1)/len(T)
                else:
                    dP_seed = 0.
                # prev BHMF
                # z_mesh = kernelS_MBHmesh(M_BH, abin_mf, t_life, l_cut)
                z_mesh = kernelS_MBH_M_mesh(M_BH, abin_mf, t_life, 1., l_cut, d_fit)*l_cut
                z_mesh[z_mesh<lambda_0] = lambda_0
                Ps = integral(a,z_mesh/l_cut,lambda_0/l_cut)/I_toinf
                dP_MBH = np.nansum( (Ps[:,:-1]-Ps[:,1:])*dP_MBH_prev, axis=1) + dP_seed

                Nt -= 1
            dn_MBH += dP_MBH*n_base[iM]*f_bsm[i_bsm]*f_seed

    consv_ratio = np.nansum(dn_MBH)/(np.sum(n_base[:N_Mh])*f_seed)
    if abs(consv_ratio-1)>.5:
        print('theta: ',theta,'x0',x0, 'consv_ratio: ',consv_ratio)
        # assert 0

    # 10 N_M in 1e7-1e10 range, plus 12 N_L
    index = np.where(np.logical_and(1e7<M_BH,M_BH<1e10))
    xs = M_BH[index][::len(index[0])//10]
    ys = np.log10( MF(xs)  ) # Willott 2010 as data
    y_model = np.log10( (dn_MBH/dlog10M)[index][::len(index[0])//10] )
    mf_err = pow(np.log10(xs)-8.5,2)/3. + .2 # from 0.2 to 0.95
    Chi2_M =  np.sum( pow((ys - y_model)/mf_err, 2))

# # --------- Luminosity Function ---------
    z_mesh = kernelS_M1450_mesh(bin_edg, M_BH, l_cut)*l_cut
    z_mesh[z_mesh<lambda_0] = lambda_0
    Ps = integral(a,z_mesh/l_cut,lambda_0/l_cut)/I_toinf
    dPhi_mesh = np.nansum((Ps[:-1,:]-Ps[1:,:])*dn_MBH,axis=1)
    Phi = dPhi_mesh/bin_wid
    Phi *= 1e9
    Phi_DO = Phi/corr_U14D20(bin_cen)
    ys = np.log10(Phi_obs)
    y_model = np.log10(Phi_DO)
    lf_err = np.log10(Phi_err)
    Chi2_L = np.sum( pow((ys - y_model)/lf_err, 2))
    # print('in model -.5*(Chi2_L+Chi2_M)',-.5*(Chi2_L+Chi2_M))

    # for spread plot
    z_mesh = kernelS_M1450_mesh(abin_lf, M_BH, l_cut)*l_cut
    z_mesh[z_mesh<lambda_0] = lambda_0
    Ps = integral(a,z_mesh/l_cut,lambda_0/l_cut)/I_toinf
    dPhi_mesh = np.nansum((Ps[:-1,:]-Ps[1:,:])*dn_MBH,axis=1)
    Phi = dPhi_mesh/dmag
    Phi *= 1e9
    Phi_DO = Phi/corr_U14D20(M1450)
    return {'M_BH':M_BH, 'MF':dn_MBH/dlog10M, 'MF_data':MF(M_BH), 'MF_data_err':mf_err,
            'Chi2_M':Chi2_M,
            'M1450_data':bin_cen, 'LF_data':Phi_obs, 'LF_data_err':Phi_err,
            'Chi2_L':Chi2_L,
            'M1450':M1450, 'LF':Phi_DO}