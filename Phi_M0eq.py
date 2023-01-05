from PYmodule import *
from PYmodule.l_intg import *

t1 = time.time()
prex = '../data_eq/'

d_fit = pow(10., -.5)
# d_fit = 0
l_cut = 1.
a = .1
t_life = 10
tz = 200*Myr

x0 = lambda_0/l_cut
I_toinf = integral_toinf(a,x0)

print('t_life, d_fit, l_cut, a,  f_seed, lambda_0, logM0 = ', 
t_life,', ',d_fit,', ', l_cut,', ',a,', ', f_seed,', ', lambda_0,', ', logM0,', ')
print('mf len:',len(abin_mf))

t_life *= Myr

t_col = 20*Myr
Mstar0 = 1e3

f_bsm = 1.
n_base = n_base[0]

## --------- Mass Function ---------
dn_MBH = np.zeros(N_mf)
Nt = np.max((tz-t_col)//t_life)
Nmax = Nt
dP_MBH = np.zeros(N_mf)
# t_tot = np.zeros(len(T))
DT = 0
# M0s = np.logspace(2,12,num=10000)
# dlog10M0 = np.log10(M0s[1]/M0s[0])

while Nt>=0:
# for i__ in range(3):
    t_point = tz - Nt*t_life
    # new seeds (using 2d meshgrids)
    if t_point-t_life<t_col and t_col <=t_point:
        dt_seed = t_point - t_col
        # kernelS_MBH_M_mesh(M1, M0, dt, f_duty, l_cut, d_fit, logM_0=logM0):
        x1 = kernelS_MBH_M(abin_mf, Mstar0, dt_seed,1.,1.,d_fit)
        x1[x1<x0] = x0

        Ps = integral(a,x1,x0)/I_toinf
        dP_seed = Ps[1:] - Ps[:-1]
        print(np.nansum(dP_seed),'is sum of dP_seed')
        DT+=dt_seed
    else:
        dP_seed = np.zeros(N_mf)
    # prev BHMF
    DT+= t_life
    dP_MBH_prev = dP_MBH.copy()
    dP_MBH = np.zeros(N_mf)

    M1s = M_BH
    for iM0 in range(N_mf): # np.arange(N_mf-1, N_mf-100, -1):
        dP_MBH_prev_ = dP_MBH_prev.copy()
        x1 = kernelS_MBH_M(M1s, M_BH[iM0],t_life,1.,l_cut,d_fit)
        dlnldlogM1 = np.log(10.) * M1s/1e8/(x1*l_cut) * .5/(t_life/(0.1*t_Edd)) * (1.e8/M1s + pow(M1s/1.e8,d_fit-1.))
        dPdlnl = pow(x1,a)*np.exp(-x1)/I_toinf
        dPdlnl[x1<x0] = 0
        dlnldlogM1[x1<x0] = 0
        # print(x1[x1>x0][:10],x0)
        # if np.any(x1>x0):
        #     dPdlnl[x1<x0] = 0
        #     dlnldlogM1[x1<x0] = 0
        #     dP_MBH_prev_[x1<x0] = 0
        #     dlnldlogM1[np.argmax(x1>x0)] = np.log(x1[np.argmax(x1>x0)]/x0)*2./dlog10M
        dP_MBH += dPdlnl*dlnldlogM1*dP_MBH_prev_*dlog10M
        # print('M0 grow,total Pl:',np.nansum(dPdlnl*dlnldlogM1*dlog10M))
        # if iM0==10:
        #     exit(0)
    dP_MBH += dP_seed/dlog10M
    # print('P_tot',P_tot)
    print('Nt: ',Nt,'each cycle: consv_ratio =',np.nansum(dP_MBH*dlog10M))
    # MFname = prex+'Phi_M0eqMFbin{0:d}_{1:d}'.format(len(abin_mf),i__)
    # ascii.write( Table([np.log10(M_BH), dP_MBH/dlog10M],
    #             names=['M_BH','dPdlogM']),
    #             MFname,
    #             formats={'M_BH':'4.2e','dPdlogM':'4.2e'},
    #             overwrite=True)
    Nt -= 1
print('DT=',np.mean((DT-t_life)/Myr))
dn_MBH = dP_MBH*n_base*f_bsm*f_seed*dlog10M

consv_ratio = np.nansum(dP_MBH*dlog10M)
print('len of abin_mf:',len(abin_mf))
print('in Phi_easy: MF conserved fraction=%.10f'%consv_ratio)

T = Table(
    [M_BH, dn_MBH/dlog10M, MF(M_BH)],
    names=('M_BH','Phi','W10_MF')
)

MFname = prex+'Phi_M0eqMF'
ascii.write( Table([np.log10(T['M_BH']), T['Phi'], T['W10_MF']],
            names=['M_BH','Phi','W10_MF']),
            MFname,
            formats={'M_BH':'4.2e','Phi':'4.2e','W10_MF':'4.2e'},
            overwrite=True)

lbin = np.arange(-3,3,0.1)
x = np.logspace(lbin[0]-np.log10(l_cut),lbin[-1]-np.log10(l_cut),num=len(lbin)) # for Pana
Pana = integral(a,x,x0)/I_toinf
with open(prex+"Phi_M0eqERDFz6.txt",'w') as f:
    np.savetxt(f, np.array([lbin[:-1],Pana[1:]-Pana[:-1]]).transpose(), fmt='%10.3e')

print('time=',time.time()-t1)