from PYmodule import *
# seed BHMF for all Mh, v_bsm; migrated from src_prev/Mdist_seed.py
# no big diff: all seed MF compared w/ Mh=1e11,vbsm=0 seed MF (generated from src_prev/Phi_seed.py)
# current plot directly, better write hist file later
log10Ms = [11,12,13]
typenames = ['H'+r'$_2$', 'H-H'+r'$_2$', 'H-H']
matplotlib.rcParams['hatch.linewidth'] = .5 
from matplotlib.ticker import LogLocator

print('z:35, t_Hubble: ', t_from_z(35)/Myr)
print('z:20, t_Hubble: ', t_from_z(25)/Myr)
print('t_edd: ', t_Edd/Myr)

abin_mf =  np.logspace(2,12,num=100) # default endpoint=True
M_BH = abin_mf[:-1]*np.sqrt(abin_mf[1]/abin_mf[0])
bin_left = abin_mf[:-1]; bin_right = abin_mf[1:]
wid_mf = bin_right - bin_left
dlog10M = np.log10(abin_mf[1]/abin_mf[0]) # print('Mbin ratio',abin_mf[1]/abin_mf[0])
N_mf = len(abin_mf)-1

# print(len(Ts[0])); exit(0)
# Ts shape: [3][2]
Ts_iso = [] # n_tell = 1.e4, T_tell = 4000K
Ts_H2 = []
Ts_isofail = []
Ts_isoOK = []
pres = ['../data/1e11','../data/1e12','../data/1e13']
figprefix = '../'

dz = 5
zs = np.arange(5,51,dz)
# [ 5 10 15 20 25 30 35 40 45 50]
dz = -10
zs = np.arange(50,5,dz)
# [10 20 30 40 50]

Nz = len(zs)-1
NM = 3
Nbsm = 2

# # check eta=0.6, 0.3 difference
# for iM in range(NM):
#     for i_bsm in range(Nbsm):
#         for eta in [.6,.3]:
#             T = Ts[iM][i_bsm]
#             for j in range(len(T)):
#                 T['Mstar0'][j] = Mdot2M(T['Mdot'][j]*eta)
#             print('eta={0:.1f} iM={1:d} i_bsm={2:d} ;'.format(eta,iM,i_bsm), 'min, max of Mstar0: {0:.1e}, {1:.1e} Msun'.format(np.min(T['Mstar0']),np.max(T['Mstar0'])))
# # fiducially eta=0.3

Tzs = [ [[0 for i in range(NM)] for j in range(Nbsm)] for k in range(Nz)]

whole = 0
for iM in range(NM):
    for i_bsm in range(Nbsm):
        T = Ts[iM][i_bsm]
        for iz in range(Nz):
            Tzs[iz][i_bsm][iM] = T[np.logical_and(zs[iz+1]<=T['z_col'], T['z_col']<zs[iz])]
            # whole += len(Tzs[iz][i_bsm][iM])
        #  #------------- print check selected table correct ----------------
        #     if Tzs[iz][i_bsm][iM]:
        #         print(zs[iz],' to ',zs[iz+1],
        #         np.max(Tzs[iz][i_bsm][iM]['z_col']),np.min(Tzs[iz][i_bsm][iM]['z_col']),
        #         'ibsm',i_bsm,np.mean(Tzs[iz][i_bsm][iM]['i_bsm']),
        #         'iM',iM,log10Ms[iM] )
        # if iM==0 and i_bsm==0:
        #     print(np.max(T['Mstar0'])); exit(0)
# print(whole); exit(0)

T_tell = 8000
x = (bin_left+bin_right)/2.


# ----------- seperate by fbsm ---------------
for i_bsm in range(Nbsm):
    h_0 = np.zeros(N_mf)
    h_1 = np.zeros(N_mf)
    h_2 = np.zeros(N_mf)
    for iM in range(NM):
        T = Ts[iM][i_bsm]
        # print(np.min(T['z_col']))
        T_H2 = T[T['Tg_max']<=T_tell] 
        T_isofail = T[np.logical_and(T['Tg_max']>T_tell, T['iso_col']==0)]
        T_isoOK =  T[np.logical_and(T['Tg_max']>T_tell, T['iso_col']==1)]
        # if T_isoOK:
        #     print(np.min(T_isoOK['Mstar0']), np.max(T_isoOK['Mstar0']))
        hist0, bin_edges = np.histogram(T_H2['Mstar0'],bins=abin_mf,density=False)
        hist1, bin_edges = np.histogram(T_isofail['Mstar0'],bins=abin_mf,density=False)
        hist2, bin_edges = np.histogram(T_isoOK['Mstar0'],bins=abin_mf,density=False)

        h_0 += hist0*n_base[iM]; h_1 += hist1*n_base[iM]; h_2 += hist2*n_base[iM]
    h_0 *= f_bsm[i_bsm]; h_2 *= f_bsm[i_bsm]; h_2 *= f_bsm[i_bsm]
    
    plt.figure(figsize=(10,8),dpi=400)
    plt.bar(x,h_0/1.e4,width=wid_mf,color='C'+str(0),alpha=0.5,label=typenames[0])
    plt.bar(x,h_1/1.e4,width=wid_mf,bottom=h_0/1e4,color='C'+str(1),alpha=0.5,label=typenames[1])
    plt.bar(x,h_2/1.e4,width=wid_mf,bottom=(h_0+h_1)/1e4,color='C'+str(2),alpha=0.5,label=typenames[2])
    plt.tick_params(labelsize=fstick)
    plt.xlabel(r'$\mathrm{M_{\bullet}}$'+r' $(M_{\odot})$',fontsize=fslabel)
    plt.ylabel(r'$\mathrm{dn/d\logM~[Mpc^{-3}dex^{-1}]}$',fontsize=fslabel)
    plt.xscale('log'); plt.yscale('log')
    plt.title(r'$\mathrm{v_{bsm}}=$'+str(i_bsm)+r' $\sigma$',fontsize=fslabel)
    plt.xlim(1e2,1e6); plt.ylim(1.e-12,1e-2)
    plt.tight_layout()
    plt.savefig(figprefix+'eta3bsm'+str(i_bsm)+'.png')

# ----------- seperate by z_col ---------------
# ----------- + cumulative over z ---------------
fig, axs = plt.subplots(2, 2, sharex=True, sharey=True)
h_allz = np.zeros(N_mf)
for iz in range(Nz):
    h = np.zeros(N_mf)
    for i_bsm in range(Nbsm):
        hh = np.zeros(N_mf)
        for iM in range(NM):
            T = Tzs[iz][i_bsm][iM]
            hist, bin_edges = np.histogram(T['Mstar0'],bins=abin_mf,density=False)
            hh = hh + hist*n_base[iM]/1e4
        h += hh*f_bsm[i_bsm]
    h_allz += h
    # try colors: darkturquoise, hatched region not good
    # axs[iz//2,iz%2].bar(x,h_allz/dlog10M,width=wid_mf,color='darkturquoise',alpha=.6,edgecolor=None,label=typenames[0])
    # axs[iz//2,iz%2].bar(x,h/dlog10M, bottom=(h_allz-h)/dlog10M,width=wid_mf,edgecolor=None, hatch="//",alpha=0.0,label=typenames[0], linewidth=.5)
    axs[iz//2,iz%2].bar(x,h/dlog10M, bottom=(h_allz-h)/dlog10M,width=wid_mf,edgecolor='k', color='darkorchid',alpha=0.6,label=typenames[0], linewidth=.5)
    axs[iz//2,iz%2].bar(x,(h_allz-h)/dlog10M,width=wid_mf,color='lightseagreen',edgecolor='k',alpha=0.6,label=typenames[0],linewidth=.5)
    # order matters!
    # ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(10))
    # ax.set_yscale('log')
    axs[iz//2,iz%2].set_xscale('log'); axs[iz//2,iz%2].set_yscale('log')
    axs[iz//2,iz%2].set_xlim(1e2,1e6); axs[iz//2,iz%2].set_ylim(1.e-8,1e-2)
    axs[iz//2,iz%2].text(3e4,1e-3,'z='+str(int(zs[iz]))+r'$\to$'+str(int(zs[iz+1])),fontsize=10)


    # axs[iz//2,iz%2].tick_params(which='minor', length=4, color='r')
    locmajx = LogLocator(base=10,numticks=100) 
    locminx = LogLocator(base=10,subs=np.arange(2, 10) * .1,numticks=100) # subs=(0.2,0.4,0.6,0.8)
    locmajy = LogLocator(base=100,numticks=100) 
    locminy = LogLocator(base=10,subs=np.arange(2, 10) * .1,numticks=100) # subs=(0.2,0.4,0.6,0.8)
    axs[iz//2,iz%2].xaxis.set_major_locator(locmajx)
    axs[iz//2,iz%2].xaxis.set_minor_locator(locminx)
    axs[iz//2,iz%2].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    axs[iz//2,iz%2].yaxis.set_major_locator(locmajy)
    axs[iz//2,iz%2].yaxis.set_minor_locator(locminy)
    axs[iz//2,iz%2].yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())


# axs[0,0].set_ylabel(r'$\mathrm{\Phi_{M_\bullet}^{seed}~(Mpc^{-3}dex^{-1})}$',fontsize=11)
# axs[1,0].set_ylabel(r'$\mathrm{\Phi_{M_\bullet}^{seed}~(Mpc^{-3}dex^{-1})}$',fontsize=11)
# axs[1,0].set_xlabel(r'$\mathrm{{M_\bullet}~(M_\odot)}$',fontsize=11)
# axs[1,1].set_xlabel(r'$\mathrm{{M_\bullet}~(M_\odot)}$',fontsize=11)

fig.savefig(figpre+'seed.pdf',dpi=300,bbox_inches='tight')
exit(0)

# ----------- all z_col sample ---------------
h_0 = np.zeros(N_mf)
h_1 = np.zeros(N_mf)
h_2 = np.zeros(N_mf)
for i_bsm in range(Nbsm):
    hh_0 = np.zeros(N_mf)
    hh_1 = np.zeros(N_mf)
    hh_2 = np.zeros(N_mf)
    for iM in range(NM):
        T = Ts[iM][i_bsm]

        T_H2 = T[T['Tg_max']<=T_tell] 
        T_isofail = T[np.logical_and(T['Tg_max']>T_tell, T['iso_col']==0)]
        T_isoOK =  T[np.logical_and(T['Tg_max']>T_tell, T['iso_col']==1)]

        hist0, bin_edges = np.histogram(T_H2['Mstar0'],bins=abin_mf,density=False)
        hist1, bin_edges = np.histogram(T_isofail['Mstar0'],bins=abin_mf,density=False)
        hist2, bin_edges = np.histogram(T_isoOK['Mstar0'],bins=abin_mf,density=False)
        hh_0 = hh_0 + hist0*n_base[iM]/1e4
        hh_1 = hh_1 + hist1*n_base[iM]/1e4
        hh_2 = hh_2 + hist2*n_base[iM]/1e4
    h_0 = h_0 + hh_0*f_bsm[i_bsm]
    h_1 = h_1 + hh_1*f_bsm[i_bsm]
    h_2 = h_2 + hh_2*f_bsm[i_bsm]
plt.figure(figsize=(10,8),dpi=400)
plt.bar(x,h_0,width=wid_mf,color='C'+str(0),alpha=0.5,label=typenames[0])
plt.bar(x,h_1,width=wid_mf,bottom=h_0,color='C'+str(1),alpha=0.5,label=typenames[1])
plt.bar(x,h_2,width=wid_mf,bottom=(h_0+h_1),color='C'+str(2),alpha=0.5,label=typenames[2])
plt.tick_params(labelsize=fstick)
s=r'$v_{bsm}=$'+str(i_bsm)+r'$\sigma$'
plt.xlabel(r'$\mathrm{M_{\bullet}}$'+r' $(M_{\odot})$',fontsize=fslabel)
plt.ylabel(r'$\mathrm{dn/d\logM~[Mpc^{-3}dex^{-1}]}$',fontsize=fslabel)
plt.xscale('log'); plt.yscale('log')
# plt.title(r'$\eta=$'+'0.'+str(int(10*eta)),fontsize=fslabel)
plt.xlim(1e2,1e6); plt.ylim(1.e-8,1e-2)
plt.legend(fontsize=fslegend,loc='best')
plt.tight_layout()
plt.savefig(figprefix+'eta'+str(int(10*eta))+'allz.png')
# all z,Mh,vbsm seed BHMF
ascii.write(Table([x,(h_0+h_1+h_2)/dlog10M]),'../IMFhist.dat',names=['M','Phi'],
formats={'M':'10.2e','Phi':'10.2e'},overwrite=True)

# ----------- seperate by Mbase ---------------
h_allM = np.zeros(N_mf)
for iM in range(NM):
    h = np.zeros(N_mf)
    for i_bsm in range(Nbsm):
        T = Ts[iM][i_bsm]
        hist, bin_edges = np.histogram(T['Mstar0'],bins=abin_mf,density=False)
        h += hist * f_bsm[i_bsm]
    # h *= 1.e-4
    h *= n_base[iM]/1e4
    h_allM += h

    plt.figure(figsize=(10,8),dpi=400)
    plt.bar(x,h_allM/dlog10M,width=wid_mf,color='C0',alpha=0.5,edgecolor=None,label=typenames[0])
    plt.bar(x,h/dlog10M, bottom=(h_allM-h)/dlog10M,width=wid_mf,edgecolor=None, hatch="//",alpha=0.5,label=typenames[0])
    # plt.bar(x,h_allz/dlog10M,width=wid_mf,color='C0',alpha=0.5,label=typenames[0])
    # plt.bar(x,h/dlog10M, bottom=(h_allz-h)/dlog10M,width=wid_mf,edgecolor='black',hatch="//",alpha=0.5,label=typenames[0])
    plt.tick_params(labelsize=fstick)
    plt.xlabel(r'$\mathrm{M_{\bullet}}$'+r' $(M_{\odot})$',fontsize=fslabel)
    plt.ylabel(r'$\mathrm{dn/d\logM~[Mpc^{-3}dex^{-1}]}$',fontsize=fslabel)
    # plt.ylabel('fraction',fontsize=fslabel)
    plt.xscale('log'); plt.ylim(0,3)
    plt.yscale('log'); plt.ylim(1.e-7,1e-2)
    plt.title(r'$\mathrm{M_h=1e}$'+str(int(log10Ms[iM])),fontsize=fslabel)
    plt.xlim(1e2,1e6)
    plt.tight_layout()
    plt.savefig(figprefix+'M'+str(int(log10Ms[iM]))+'.png')