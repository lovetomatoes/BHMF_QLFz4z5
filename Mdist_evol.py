from PYmodule import *
from PYmodule.l_intg import *
from datetime import datetime
# paras -> BH growth to z_end=6; 
# at z1=7.642, capture BH matching Wang21's BH. (M_BH & lambda)
# output: Mevol, t(z)evol, levol; BH M_BH,lambda,L,M1450 at z=6

time0 = time.time()

z1 = 7.642
t_1 = t_from_z(z1)
z_end = 6
t_end = t_from_z(z_end)

T = Ts[0][0]
Nsite = len(T) # 1e4

f_bsm = 1.
n_base = n_base[0]
N_concatenate = int(1e2) # Nsite * N_concatenate samples
N = Nsite*N_concatenate
M_low = 1e9

# 不同f_seed 不影响BH growth 只是参数不同; 影响BHMF hist normalization
# 用f_seed=0.01 相当于1e2*N_BH sample

# M0r8 bests:
# f1 [20.07157851 -2.98140382  0.89453609  0.12195823] -4.137038049964094
t_life, logd_fit, l_cut, a = 20.07157851, -2.98140382,  0.89453609,  0.12195823; f_seed = 0.1
# f2 [18.7555167  -1.2574505   0.87372563  0.20389703] -3.1054538991409824
t_life, logd_fit, l_cut, a = 18.7555167,  -1.2574505,   0.87372563,  0.20389703; f_seed = 0.01

d_fit = pow(10.,logd_fit)
t_life *= Myr

f_seedlabel = 'f%d'%abs(int(np.log10(f_seed)))
# table stores the cdf of lambda
x0 = lambda_0/l_cut
I_toinf = integral_toinf(a,x0)
x = np.logspace(np.log10(x0),1.2,num=200)
Pa_x = integral(a,x,x0)/I_toinf

Nt = int(np.max((t_end-T['t_col'])//t_life)+1)

# # N_bri: count active bright quasar (L>L_bright) numbers from all sample at time t
# # potentially comparable with luminous quasar density Wang+2019b
# N_bri = [0]; ts = [0]; L_bright = 1e47

# print('Mstar0',T['Mstar0'])
# print("Nsite",Nsite, 't:',t/Myr,'t_end:',t_end/Myr, 'Nt:',Nt)

prex = '../4p/distf{0:d}N{1:d}_'.format(int(abs(np.log10(f_seed))),int(np.log10(N)))+datetime.now().strftime('%m%d%H%M_')
prex = '../temp/distf{0:d}N{1:d}_'.format(int(abs(np.log10(f_seed))),int(np.log10(N)))+datetime.now().strftime('%m%d%H%M_')
Mfname = prex+'Mevol.txt'
tfname = prex+'tevol.txt'
zfname = prex+'zevol.txt'
lfname = prex+'levol.txt'
z6fname= prex+'BHatz6.txt'

with open(Mfname,'w') as fM, open(z6fname, 'w') as fz6, open(lfname, 'w') as fl, open(tfname,'w') as ft, open(zfname,'w') as fz:
    fz6.write('{0:10s}{1:10s}{2:10s}{3:10s}\n'.format('M1','ls','L1','M1450_1'))
    for i_concatenate in range(N_concatenate):
        M0 = T['Mstar0']; t = T['t_col']; 
        M1s = [M0]; ts =  [t/Myr]; l1s = [.01*np.ones(len(M0))]
        zs = [T['z_col']]
        if not i_concatenate%100:
            print('i_concatenate/N_concatenate={:d}/{:d}'.format(i_concatenate,N_concatenate))
        # index store all qualified BHs in i_concatenate
        match_inallcycle = np.zeros(Nsite,dtype=bool)
        cross_inallcycle = np.zeros(Nsite,dtype=bool)
        for i_ in range(Nt):
            if np.any(t + t_life > t_end): 
                dt = t_life * np.ones(Nsite)
                dt[ np.where(t + t_life > t_end) ]  = t_end - t[ np.where(t + t_life > t_end)] 
                t = t + dt
            else:
                t = t + t_life
                dt = t_life
            # #  ---------   generate ls by table    ---------
            # uniform_a = np.random.uniform(size=Nsite)
            # ls = np.zeros(Nsite)
            # for i in range(Nsite):
            #     ls[i] = x[np.argmax(Pa_x>uniform_a[i])]
            # ls = ls*l_cut
            # ---------   a>0; directly use gamma dist   ---------
            ls = np.random.gamma(a, l_cut, Nsite)
            while np.sum(ls<lambda_0):
                ls[ls<lambda_0] = np.random.gamma(a, l_cut, np.sum(ls<lambda_0))

            M1 = M1M0_d(M0,ls,dt,d_fit)
            # M1 = M1M0_e(M0,dt,ls) #试了 MF还是对不上

            M0 = M1
            L1 = L_M(M1,ls)
            M1450_1 = M1450_Lbol(L1)
        # select close to Wang2021: M = 1.6e9, l=0.67
            cross_in1cycle = np.logical_and(t_1 > t-dt, t_1 <= t )
            cross_inallcycle[cross_in1cycle] = True 
            # match and cross
            match_in1cycle = np.logical_and(np.logical_and(np.logical_and(1.5e9<M1,M1<2e9),
                                                           np.logical_and(0.4<ls,ls<0.8)),
                                            cross_in1cycle)
            # match_in1cycle = np.logical_and(np.logical_and(np.logical_and(1.e9<M1,M1<2.e9),
            #                                                np.logical_and(ls,ls)),
            #                                 cross_in1cycle)
            match_inallcycle[match_in1cycle] = True
            # print matched objects
            for _i in range(Nsite):
                if match_in1cycle[_i]:
                    print('M1 match= {:10.3e}, l match = {:10.3e}'.format(M1[_i],ls[_i]))
            M1s.append(M1); ts.append(t/Myr); l1s.append(ls); zs.append(z_tH(t/Myr))
            # N_bri.append(len(np.where(L1>=L_bright)[0]))
        # print(np.allclose(t_end,t))

        M1s=np.array(M1s);l1s=np.array(l1s);ts=np.array(ts);zs=np.array(zs)
        index = np.where(M1>M_low)
        np.savetxt(fM, M1s.transpose()[index], fmt='%10.3e')
        np.savetxt(fl, l1s.transpose()[index], fmt='%10.3e')
        np.savetxt(ft, ts.transpose()[index], fmt='%10.3e')
        np.savetxt(fz, zs.transpose()[index], fmt='%10.3e')
        np.savetxt(fz6, np.array([M1,ls,L1,M1450_1]).transpose(), fmt='%10.3e')

pyname = sys.argv[0]
print(pyname,' time: ',time.time()-time0)
print('np.min(L1)=%.1e'%np.min(L1))
print('np.max(L1)=%.1e'%np.max(L1))
print('np.max(M1)=%.1e'%np.max(M1))

exit(0)
T_z6 = ascii.read(prex+'BHatz6.txt', guess=False, delimiter=' ')
M1, ls, L1, M1450_1 = T_z6['M1'], T_z6['ls'], T_z6['L1'], T_z6['M1450_1']

# MF histogram comparing w/ Phi_easy 
hist, bin_edges = np.histogram(M1,bins=abin_mf,density=False)
hist = hist*n_base*f_bsm*f_seed/N/dlog10M
ascii.write(Table([np.log10(M_BH),hist,MF(M_BH)]), 
z6datapre+f_seedlabel+'Mdist_evolMFz{:d}_t{:d}_d{:d}'.format(int(z1),int(t_life/Myr),int(abs(logd_fit))),
names=['M_BH','hist','W10'],
formats={'M_BH':'10.2e','hist':'10.2e','W10':'10.2e'},
overwrite=True)

# MF histogram comparing w/ Phi_easy 
hist, bin_edges = np.histogram(M1450_1,bins=bin_edg,density=False)
hist = hist*n_base*f_bsm*f_seed/N/bin_wid*1e9
hist_DO = hist/corr_U14D20(bin_cen)
ascii.write(Table([bin_cen,Phi_obs,hist_DO,hist]), 
z6datapre+f_seedlabel+'Mdist_evolLFz{:d}_t{:d}_d{:d}'.format(int(z1),int(t_life/Myr),int(abs(logd_fit))),
names=['M1450','Phi_obs','hist_DO','hist'],
formats={'M1450':'10.2f','hist':'10.2e','hist_DO':'10.2e','Phi_obs':'10.2e'},
overwrite=True)

exit(0)
# write lambda histogram: lbin, hist_tot, hist_L45, hist_L46, Pana; (into 2 files...)
lbin = np.linspace(-2,1.2,num=20)
hist, bin_edges = np.histogram(np.log10(ls),bins=lbin,density=False)
x = np.logspace(lbin[0]-np.log10(l_cut),lbin[-1]-np.log10(l_cut),num=len(lbin)) # for Pana
Pana = integral(a,x,x0)/I_toinf

for L_limit in [1e45, 1e46]:
    index = np.where(L_limit<L1)
    M1_ = M1[index]; L1_ = L1[index]; ls_ = ls[index]
    print('len(ls) after selection:',len(ls_),' / ',Nsite*N_concatenate)
    print('len of all L1',len(L1))
    hist_, bin_edges = np.histogram(np.log10(ls_),bins=lbin,density=False)

    strhist_ = 'hist_L%d'%(int(np.log10(L_limit)))
    # lambda histogram file
    ascii.write(Table([bin_edges[:-1],hist_/len(ls),hist/len(ls),Pana[1:]-Pana[:-1]]), 
    prex+strhist_+'.txt',
    names=['log_l',strhist_,'hist_tot','ana'],
    formats={'log_l':'10.2f',strhist_:'10.2e','hist_tot':'10.2e','ana':'10.2e'},
    overwrite=True)