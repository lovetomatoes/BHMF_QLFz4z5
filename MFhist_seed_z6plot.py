from PYmodule import *
from PYmodule.l_intg import *
# inherited from sk1:BHMF_prev/plotBHMF.py, dist_grow.py;
# plot histogram version of seed BHMF, after grow, Willott2010 curve.
# some files generated from MFhist_seed_z6
t1 = time.time()

z1 = 6
t_end = t_from_z(z1)

T = Ts[0][0]

Nsite = len(T) # N3 <-> 1e7
M0 = T['Mstar0']

f_bsm = 1.
n_base = n_base[0]
N_concatenate = int(1e1)

h_seed, bin_edges = np.histogram(M0,bins=abin_mf,density=False)
Phi_seed = h_seed*n_base/float(Nsite)/dlog10M

t_life, d_fit, l_cut, a = 19.8, 1.2e-3, 1.1557, -1.8e-01 # f_seed = 1.
f_seed = 0.01
t_life, d_fit, l_cut, a = 1.68814272e+01, 7.31406767e-03, 1.02157385e+00, 1.46099557e-01 # f_seed = .01

t_life *= Myr

# # table stores the cdf of lambda 
# I_toinf = integral_toinf(a)
# x = np.logspace(np.log10(x0),1.2,num=200)
# Pa_x = integral(a,x)/I_toinf

# ascii.write(Table([x,Pa_x]),'../Pa.dat',
# names=['x','Pa_x'],
# formats={'x':'10.5e','Pa_x':'10.5e'},
# overwrite=True)

## --------- Mass Function ---------
hist_mf = np.zeros(N_mf)
h_mf = np.zeros(N_mf)

Nt = int(np.max((t_end-T['t_col'])//t_life))

t = T['t_col']
# print(t/Myr, 't_end:', t_end/Myr, 'Nt:',Nt)

# for i_concatenate in range(N_concatenate):
#     for i_ in range(Nt+1):
#         if np.any(t + t_life > t_end): 
#             dt = t_life * np.ones(Nsite)
#             dt[ np.where(t + t_life > t_end) ]  = t_end - t[ np.where(t + t_life > t_end)] 
#             t = t + dt
#         else:
#             t = t + t_life
#             dt = t_life
#         uniform_a = np.random.uniform(size=Nsite)
#         ls = np.zeros(Nsite)
#         for i in range(Nsite):
#             ls[i] = x[np.argmax(Pa_x>uniform_a[i])]
#         ls = ls*l_cut

#         M1 = M1M0_d(M0,ls,dt,d_fit)
#         M0 = M1
#     # ascii.write(Table([t/Myr,np.ones(Nsite)*t_end/Myr]),'../t_evol.dat',
#     # names=['t','t_end'],
#     # formats={'t':'10.5e','t_end':'10.5e'},
#     # overwrite=True)

#     hist_mf, bin_edges = np.histogram(M1,bins=abin_mf,density=False)


#     h_mf += hist_mf

# Phi_mf = h_mf*n_base/float(Nsite*N_concatenate)/dlog10M
T_hist = ascii.read('../MFhist.dat', guess=False, delimiter=' ')
T_Phi = ascii.read('../MFIMF', guess=False, delimiter=' ')
M_BH = T_hist['M_BH']
Phi_seed = T_hist['Phi_seed']
Phi_mf = T_Phi['Phi']

print('time=',time.time()-t1)

plt.figure(figsize=(10,9),dpi=400)
x = (abin_mf[:-1]+abin_mf[1:])/2.
plt.bar(x, Phi_mf,width=wid_mf,color='C'+str(0),alpha=0.5,label='z=6')
plt.bar(x, Phi_seed,width=wid_mf,color='C'+str(1),alpha=0.5,label='seeds')

index = np.logical_and(1e7<x, x<1e10)
plt.plot(x[index],MF(x)[index], '--',c='black',lw=2, label='Willott+ 2010')
plt.tick_params(labelsize=fstick)
plt.xlabel(r'$\mathrm{M_{BH}}$' +' '+ r'$\left(\mathrm{M_\odot}\right)$',fontsize=fslabel)
plt.ylabel(r'$\mathrm{\Phi}$'+' '+r'$\mathrm{\left(Mpc^{-3}dex^{-1}\right)}$',fontsize=fslabel)
plt.xscale('log'); plt.yscale('log')
# plt.title('z='+str(z),fontsize=fslabel)
plt.xlim(1e2,1e10); plt.ylim(1e-10,1e-2)
# plt.grid(True)
plt.legend(fontsize=fslegend,loc='best')
plt.tight_layout()
# plt.savefig(z6figpre+'test.png')
# plt.savefig('../testN{0:.0e}.png'.format(N_concatenate))
plt.savefig('../BHMF.png'.format(N_concatenate))

# ascii.write(Table([x,Phi_seed,Phi_mf,MF(x)]),
# '../MFhist{0:.0e}.dat'.format(N_concatenate),#z6datapre+'MFhist.dat',
# names=['M_BH','Phi_seed','Phi_mf','MF10'],
# formats={'M_BH':'10.2e','Phi_seed':'10.2e','Phi_mf':'10.2e','MF10':'10.2e'},
# overwrite=True)
