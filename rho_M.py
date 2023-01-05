from PYmodule import *
from PYmodule.l_intg import *

t1 = time.time()

# M0 best
t_life, logd_fit, l_cut, a = 18.7555167,  -1.2574505,   0.87372563,  0.20389703; f_seed = 0.01
# t_life, logd_fit, l_cut, a = 20.07157851, -2.98140382,  0.89453609,  0.12195823; f_seed = 0.1
# t_life, logd_fit, l_cut, a = 23.12675104, -2.97342483,  0.95753445, -0.06535641; f_seed = 1

x0 = lambda_0/l_cut
I_toinf = integral_toinf(a,x0)
d_fit = pow(10.,logd_fit)
print('t_life={0:.1f}, d_fit={1:.2f}, l_cut={2:.2f}, a={3:.2f}, f_seed={4:.2f}'.format(t_life,  d_fit, l_cut, a, f_seed))

t_life *= Myr
T = Ts[0][0]
f_bsm = 1.
n_base = n_base[0]

L_ups = [-23,-24,-25,-26,-27]
M_up  = 1e7
fname = '../f%drho_evol.txt'%int(abs(np.log10(f_seed)))
with open(fname,'w') as f:
    f.write('{:>10}'.format('z'))
    f.write('{:>10}'.format('rho_M'+str(int(np.log10(M_up)))))
    for L in L_ups:
        f.write('{:>10}'.format(str(L)))
    f.write('\n')

    for z in np.arange(6,11,1):
        tz = t_from_z(z)
        # print('mean Dt:',np.mean((tz-T['t_col']))/Myr)

        ## --------- Mass Function ---------
        dn_MBH = np.zeros(N_mf)
        Nt = np.max((tz-T['t_col'])//t_life)
        Nmax = Nt
        dP_MBH = np.zeros(N_mf)
        while Nt>=0:
            t_point = tz - Nt*t_life
            T_seed = T[np.logical_and(t_point-t_life<=T['t_col'],T['t_col']<t_point)]
            dt_seed = t_point - T_seed['t_col']
            dP_MBH_prev = dP_MBH.copy()
            # new seeds (using 2d meshgrids)
            if len(T_seed):
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

        rho_M = np.sum(dn_MBH[M_BH>M_up])
        print('rho_1e7= %.1e'%rho_M)
        consv_ratio = np.nansum(dn_MBH)/(n_base*f_seed); print('z={0:d}\tconsv={1:.10e}'.format(z,consv_ratio))

        T_ = Table(
            [M_BH, dn_MBH/dlog10M, MF(M_BH)],
            names=('M_BH','Phi','W10_MF')
        )

        MFname = z6datapre+'f{0:d}Phi_easyMF_z{1:d}'.format(int(abs(np.log10(f_seed))),z)
        ascii.write( Table([np.log10(T_['M_BH']), T_['Phi'], T_['W10_MF']],
                    names=['M_BH','Phi','W10_MF']),
                    MFname,
                    formats={'M_BH':'4.2e','Phi':'4.2e','W10_MF':'4.2e'},
                    overwrite=True)


        # # --------- Luminosity Function ---------
        Phi_obs = LF_M1450(M1450)*1e9
        z_mesh = kernelS_M1450_mesh(abin_lf, M_BH, l_cut)
        z_mesh[z_mesh<x0] = x0
        Ps = integral(a,z_mesh,x0)/I_toinf
        dPhi_mesh = np.nansum((Ps[:-1,:]-Ps[1:,:])*dn_MBH,axis=1)
        dPhi_DO = dPhi_mesh/corr_U14D20(M1450)
        dPhi_DO = dPhi_mesh/corr_M14D20(M1450)
        T_ = Table(
            [M1450,dPhi_DO/dmag,Phi_obs],
            names=('M1450','Phi_DO','Phi_obs')
        )

        LFname = z6datapre+'f{0:d}Phi_easyLF_z{1:d}'.format(int(abs(np.log10(f_seed))),z)
        ascii.write(T_, LFname,
                    formats={'M1450':'6.2f','Phi_DO':'4.2e','Phi_obs':'4.2e'},
                    overwrite=True)
        # # z=6
        # L_low = -24; L_up = -23
        # print('L_up={0:.0f}, n={1:.1e}'.format(L_up,np.sum(dPhi_DO[np.logical_and(M1450>L_low, M1450<L_up)])) ) 
        # L_low = -25; L_up = -24
        # print('L_up={0:.0f}, n={1:.1e}'.format(L_up,np.sum(dPhi_DO[np.logical_and(M1450>L_low, M1450<L_up)])) ) 
        # L_low = -26; L_up = -25
        # print('L_up={0:.0f}, n={1:.1e}'.format(L_up,np.sum(dPhi_DO[np.logical_and(M1450>L_low, M1450<L_up)])) ) 
        # L_low = -27; L_up = -26
        # print('L_up={0:.0f}, n={1:.1e}'.format(L_up,np.sum(dPhi_DO[np.logical_and(M1450>L_low, M1450<L_up)])) ) 
        # L_low = -30; L_up = -27
        # print('L_up={0:.0f}, n={1:.1e}'.format(L_up,np.sum(dPhi_DO[np.logical_and(M1450>L_low, M1450<L_up)])) ) 

        L_ups = [-23,-24,-25,-26,-27]
        ns = []
        for L_up in L_ups:
            L_low = L_up-1 if L_up>-27 else -30
            ns.append(np.sum(dPhi_DO[np.logical_and(M1450>L_low, M1450<L_up)]))
        print(ns)
        np.savetxt(f, [np.concatenate((np.array([z,rho_M]),np.array(ns)))], delimiter='',fmt='%10.3e')

# print('time=',time.time()-t1)