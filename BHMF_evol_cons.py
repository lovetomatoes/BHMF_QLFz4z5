from PYmodule import *
from PYmodule.l_intg import *

# models: from z=6 to z=5; compare with BHMF_evol_cons.f
# gnuplot: p 'BHMFz5.txt' u 1:2 w l, 'fort.10' u 1:3 w l

z0 = 6
# z=6 BHMF 
phi_s  = 1.31e-6
M_star = 6.13e7
alp    = -1.41
beta   = -2.58

# MF bins
abin_mf =  np.logspace(5,11,num=1000) # default endpoint=True
M_BH = abin_mf[:-1]*np.sqrt(abin_mf[1]/abin_mf[0])
bin_left = abin_mf[:-1]; bin_right = abin_mf[1:]
wid_mf = bin_right - bin_left
dlog10M = np.log10(abin_mf[1]/abin_mf[0]) # print('Mbin ratio',abin_mf[1]/abin_mf[0])
N_mf = len(abin_mf)-1

def model(theta, z, f_seed=f_seed, corr='U',l_cut= l_cut, a=a):
    if len(theta)==4:
        t_life, logd_fit, l_cut, a = theta
    elif len(theta)==3:
        t_life, logd_fit, l_cut, a = theta
    else:
        assert 0
    dn_MBH = dlog10M * phi_s/(pow(M_BH/M_star,-(alp+1)) + pow(M_BH/M_star,-(beta+1)))
    n_tot = np.nansum(dn_MBH)
    # print(len(dn_MBH), n_tot)
    dP_MBH = dn_MBH/n_tot
    d_fit = pow(10.,logd_fit)
    t_life = t_life * Myr # wli: not *=!!!!!! theta changed by this
    x0 = lambda_0/l_cut
    I_toinf = integral_toinf(a,x0)

## --------- Mass Function ---------
    t_tot = (1.17e3-935)*Myr
    # print(t_Edd/ Myr); exit(0)
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


    dn_MBH = n_tot*dP_MBH
    Phi_MBH = dn_MBH/dlog10M

    consv_ratio = np.nansum(dP_MBH)
    print('consv_ratio',consv_ratio)
    if abs(consv_ratio-1)>.1:
        print('theta: ',theta,'x0',x0, 'consv_ratio: ',consv_ratio)
        assert 0
    
    T_ = Table(
        [M_BH, Phi_MBH*1e9],
        names=('M_BH','Phi_MBH_Gpc')
    )
    MFname = './BHMFz{:d}.txt'.format(z)
    ascii.write(T_,
                MFname,
                formats={'M_BH':'4.2e','Phi_MBH_Gpc':'4.2e'},
                overwrite=True)

def main():
    t_life, logd_fit, l_cut, a = 10.368000984191895, -2., 0.301, 0.89161005

    x = (t_life, logd_fit, l_cut, a)
    z = 5
    model(x,z)

main()