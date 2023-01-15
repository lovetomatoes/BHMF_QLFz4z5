from PYmodule import *
from PYmodule.l_intg import *

# models: from z0=6 (../data/M0r8_f*totBHMFz6) to z

def model(theta, z, f_seed=f_seed, corr='U',l_cut= l_cut, a=a):
    t_life, logd_fit, l_cut, a  = theta
    d_fit = pow(10.,logd_fit)
    t_life = t_life * Myr # wli: not *=!!!!!! theta changed by this
    x0 = lambda_0/l_cut
    I_toinf = integral_toinf(a,x0)


    z0 = 6
    # # ----------   z=6 BHMF   ---------- #
    #  DPL fitting function
    phi_s  = 1.31e-6
    M_star = 6.13e7
    alp    = -1.41
    beta   = -2.58
    Phi_MBH = phi_s/(pow(M_BH/M_star,-(alp+1)) + pow(M_BH/M_star,-(beta+1)))
    # ---  interpolation log-log  --- #
    MF_z6 = ascii.read('../data/M0r8_f2totBHMFz6')
    M_BHa = MF_z6['Mx']; Phi_MBHa = MF_z6['y_tot']
    Phi_MBH = pow(np.e, np.interp(np.log(M_BH), np.log(M_BHa), np.log(Phi_MBHa)))
    # Phi_MBH = np.interp(np.log(M_BH), np.log(M_BHa), Phi_MBHa)
    dn_MBH = dlog10M * Phi_MBH
    n_tot = np.nansum(dn_MBH)
    dP_MBH = dn_MBH/n_tot
    # print(len(dn_MBH), n_tot)


    ## --------- Mass Function ---------
    tz = t_from_z(z)
    # print(t_Edd/ Myr); exit(0)
    # print('t_tot',(tz-t_from_z(z0))/t_tot); exit(0)
    t_tot = tz-t_from_z(z0)
    t = 0
    while t<t_tot:
        # prev BHMF
        dP_MBH_prev = dP_MBH.copy()
        # print('np.nansum(dP_MBH_prev)',np.nansum(dP_MBH_prev))
        z_mesh = kernelS_MBH_M_mesh(abin_mf, M_BH, t_life, 1., l_cut, d_fit)
        z_mesh[z_mesh<x0] = x0
        Ps = integral(a,z_mesh,x0)/I_toinf
        dP_MBH = np.nansum( (Ps[1:,:]-Ps[:-1,:])*dP_MBH_prev, axis=1)
        t += t_life

    dn_MBH = n_tot*dP_MBH * f_seed
    Phi_MBH = dn_MBH/dlog10M

    consv_ratio = np.nansum(dP_MBH)
    # print('consv_ratio',consv_ratio)
    if abs(consv_ratio-1)>.01:
        print('theta: ',theta,'x0',x0, 'MF consv_ratio: ',consv_ratio)
        # assert 0

    # T_ = Table(
    #         [M_BH, Phi_MBH],
    #         names=('M_BH','Phi_MBH')
    #         )
    # MFname = './interpBHMFz{:d}.txt'.format(z)
    # ascii.write(T_,
    #         MFname,
    #         formats={'M_BH':'4.2e','Phi_MBH':'4.2e'},
    #         overwrite=True)

    # # chi2 in BHMF
    # # 10 N_M in 1e7-1e10 range, plus 12 N_L
    # index = np.where(np.logical_and(1e7<M_BH,M_BH<1e10))
    # xs = M_BH[index][::len(index[0])//10]
    # ys = np.log10( MF(xs)  ) # Willott 2010 as data
    # y_model = np.log10( (dn_MBH/dlog10M)[index][::len(index[0])//10] )
    # y_err = pow(np.log10(xs)-8.5,2)/3. + .2 # from 0.2 to 0.95
    # Chi2_M =  np.sum( pow((ys - y_model)/y_err, 2))


    # # ---------  Luminosity Function for spread plot  ---------
    # lf_edg = abin_lf
    # lf_wid = dmag
    # lf_cen = M1450
    # PhiLF_obs = LF_M1450(M1450,z=z)*1e9

    # z_mesh = kernelS_M1450_mesh(lf_edg, M_BH, l_cut)
    # z_mesh[z_mesh<x0] = x0
    # Ps = integral(a,z_mesh,x0)/I_toinf
    # # print(Ps[-1,:]-Ps[0,:])
    # dn_M1450 = np.nansum((Ps[:-1,:]-Ps[1:,:])*dn_MBH,axis=1)
    # Phi = dn_M1450/lf_wid
    # print('LF consv:',np.nansum(dn_M1450)/n_tot) # M1450 range broad enough
    # Phi *= 1e9
    # Phi_DO = Phi/corr_U14D20(lf_cen)
    # if corr=='M':
    #     Phi_DO = Phi/corr_M14D20(lf_cen)
    # T_ = Table(
    #         [lf_cen, Phi_DO, PhiLF_obs],
    #         names=('M1450','Phi_M1450','Phi_obs')
    #         )
    # LFname = './spreadQLFz{:d}.txt'.format(z)
    # ascii.write(T_,
    #         LFname,
    #         formats={'M1450':'4.2e','Phi_M1450':'4.2e','Phi_obs':'4.2e'},
    #         overwrite=True)

    # # --------- binned Luminosity Function ---------
    Chi2_L = 0
    lflabels = {5:['_1','_2'], 4:['']}
    for lflabel in lflabels[z]:
        # print('{:d}{:s}'.format(z,lflabel))
        lf_edg = bin_edg['{:d}{:s}'.format(z,lflabel)]
        lf_wid = bin_wid['{:d}{:s}'.format(z,lflabel)]
        lf_cen = bin_cen['{:d}{:s}'.format(z,lflabel)]
        PhiLF_obs = Phi_obs['{:d}{:s}'.format(z,lflabel)]
        PhiLF_err = Phi_err['{:d}{:s}'.format(z,lflabel)]

        # for spread plot
        z_mesh = kernelS_M1450_mesh(lf_edg, M_BH, l_cut)
        z_mesh[z_mesh<x0] = x0
        Ps = integral(a,z_mesh,x0)/I_toinf
        # print(Ps[-1,:]-Ps[0,:])
        dn_M1450 = np.nansum((Ps[:-1,:]-Ps[1:,:])*dn_MBH,axis=1)
        Phi = dn_M1450/lf_wid
        Phi *= 1e9
        Phi_DO = Phi/corr_U14D20(lf_cen)
        if corr=='M':
            Phi_DO = Phi/corr_M14D20(lf_cen)
        
        # T_ = Table(
        #         [lf_cen, Phi_DO, PhiLF_obs],
        #         names=('M1450','Phi_M1450','Phi_obs')
        #         )
        # LFname = './binnedQLFz{:d}{:s}.txt'.format(z,lflabel)
        # ascii.write(T_,
        #         LFname,
        #         formats={'M1450':'4.2e','Phi_M1450':'4.2e','Phi_obs':'4.2e'},
        #         overwrite=True)

        ys = PhiLF_obs
        y_model = Phi_DO
        y_err = PhiLF_err
        Chi2_L += np.sum( pow((ys - y_model)/y_err, 2))
    # print('Chi2_L=%.4e'%Chi2_L)

    return Chi2_L


# return {'M_BH':M_BH, 'MF':dn_MBH/dlog10M, 'MF_data':MF(M_BH), 
# 'MF_data_err':pow(np.log10(M_BH)-8.5,2)/3. + .2,
# 'Chi2_M':Chi2_M,
# 'M1450_data':bin_cen, 'LF_data':Phi_obs, 'LF_data_err':Phi_err,
# 'Chi2_L':Chi2_L,
# 'M1450':M1450, 'LF':Phi_DO}