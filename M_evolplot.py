from PYmodule import *
from PYmodule.l_intg import *
# analyse BH ERDF from sampling in MBH_evol

N_BH = int(1e5)
prex = z6datapre+'/f{0:.0e}N{1:d}'.format(f_seed,int(np.log10(N_BH)))
prex = '../2022-05-05_14:27:43_'

# z=6 BH mass, λ, L_bol
T_z6 = ascii.read(prex+'BHatz6.dat', guess=False, delimiter=' ')
M1, ls, L1 = T_z6['M1'], T_z6['ls'], T_z6['L1']
SMl1 = ls[M1>1e9] # super-massive BHs >1e9

T_evol = ascii.read(prex+'tN_evol.dat', guess=False, delimiter=' ')
zs, ts, N_bri = T_evol['zs'],T_evol['ts'],T_evol['N_bri']

T_M = np.loadtxt(prex+'Mevol.dat') # each line a BH growth track; all end with M>1e8
plt.figure(figsize=(10,10),dpi=400)
# all BHs>1e8
for i in range(len(T_M)):
    plt.plot(zs,T_M[i], c='grey', alpha=0.5)

M1 = T_M.transpose()[-1]
# print('M1',M1,T_M.shape,'len zs',len(zs))
# i = np.argmax(M1); print('max M1 = %.1e'%M1[i])
# plt.plot(zs,T_M[i], '--',c='black',lw=2)

SMM1 = M1[M1>1e9]
print('number of >1e9 MBH:',sum(M1>1e9))
for i in range(len(SMM1)):
    plt.plot(zs,T_M[M1>1e9][i],c='C%d'%i)
    plt.plot(zs,M1M0_e(SMM1[i],-(ts[-1]-ts)*Myr,SMl1[i]),'--',c='C%d'%i)
    plt.plot(zs,M1M0_e(SMM1[i],-(ts[-1]-ts)*Myr,1),'.',c='C%d'%i)
    plt.scatter(6,SMM1[i],marker='D',s=30,c='C%d'%i)

plt.yscale('log')
plt.xlim(np.max(zs),np.min(zs))
plt.xlim(30,5.8)
plt.ylim(1e2,1e10)
plt.grid(True)
plt.xlabel(r'$\mathrm{z}$',fontsize=fslabel)
plt.ylabel(r'$\mathrm{M_{BH}~(M_\odot)}$',fontsize=fslabel)
# plt.legend(loc='best',fontsize=fslabel)
plt.xticks(fontsize=fstick);plt.yticks(fontsize=fstick)
plt.savefig(prex+'Mevol.png', dpi=300, bbox_inches='tight')

exit(0)

# plot z=6 λ dist.
L_limit = 1e45
strhist_ = 'hist_L%d'%(int(np.log10(L_limit)))
lbin = np.linspace(-2,1.2,num=20)
T_hist = ascii.read(prex+strhist_+'.dat', guess=False, delimiter=' ')
left_of_lbin, hist_, hist_tot, ana = T_hist['log_l'],T_hist[strhist_],T_hist['hist_tot'],T_hist['ana']

print('sum(hist_tot)={0:.2e},sum(hist_):{1:.2e}'.format(np.sum(hist_tot),np.sum(hist_)))

plt.figure(figsize=(10,9))
x = (lbin[:-1]+lbin[1:])/2.
wid = lbin[1:]-lbin[:-1]
plt.bar(x, hist_tot,width=wid,color='C'+str(0),alpha=0.5,label='total')
plt.bar(x, 10*hist_,   width=wid,color='C'+str(1),alpha=0.5,label=r'$L>10^{45}$ erg/s')
# plt.bar(x, hist_,   width=wid,color='C'+str(1),alpha=0.5,label=r'$L>10^{46}$ erg/s')
# plt.yscale('log'); plt.ylim(bottom=1e-5)

# plot z=6 λ dist.
L_limit = 1e46
strhist_ = 'hist_L%d'%(int(np.log10(L_limit)))
lbin = np.linspace(-2,1.2,num=20)
T_hist = ascii.read(prex+strhist_+'.dat', guess=False, delimiter=' ')
left_of_lbin, hist_, hist_tot, ana = T_hist['log_l'],T_hist[strhist_],T_hist['hist_tot'],T_hist['ana']
plt.bar(x, 100*hist_,   width=wid,color='C'+str(2),alpha=0.5,label=r'$L>10^{46}$ erg/s')

# index = np.logical_and(1e7<x, x<1e10)
plt.plot(x, ana, '--',c='black',lw=2, label=r'$\mathrm{P(\lambda)}$')
plt.tick_params(labelsize=fstick)
plt.xlabel(r'$\log\lambda$',fontsize=fslabel)
# plt.ylabel(r'$\mathrm{\Phi}$'+' '+r'$\mathrm{\left(Mpc^{-3}dex^{-1}\right)}$',fontsize=fslabel)
# plt.xlim(1e2,1e10); plt.ylim(1e-10,1e-2)
plt.legend(fontsize=fslegend,loc='best')
plt.savefig(prex+'l_hist.png', dpi=300, bbox_inches='tight')