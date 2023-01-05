from PYmodule import *
from scipy.optimize import curve_fit

t1 = time.time()

f_seed = 0.01
f_seedlabel = 'f%d'%abs(int(np.log10(f_seed)))
prex = '../4p/M0r8_' + f_seedlabel + 'lin' # on sk1
# prex = '../M0r8_' + f_seedlabel + 'lin' # local

Mmin = 1e7
outMprex = prex +'Mmin%d'%(int(np.log10(Mmin)))
Lmax = -26
outLprex = prex +'Lmax%d'%abs(int(Lmax))


figm, axm = plt.subplots(figsize=(10, 10))

rhoM = []; rhoM_min = []; rhoM_max = []
nL = []; nL_min = []; nL_max = []
# labels = {'NIRCam_deep','NIRCam_med','Roman_deep','Roman_wide','Euclid_deep','Euclid_wide'}
styles = ['-','--','-.','.']
colors = ['k','red','orange','orange']
labels = ['LSST_deep','Roman_wide','Euclid_wide','Roman_deep','Euclid_deep']
N_detect = {'LSST_deep':[], 'Roman_wide':[], 'Euclid_wide':[]}

zs = np.arange(6,11)
for z in zs:
    figl, axl = plt.subplots(figsize=(10, 10))
    fMname = prex+'BHMFz{:d}'.format(z)
    fLname = prex+'QLFz{:d}'.format(z)
    TM = ascii.read(fMname, guess=False, delimiter=' ')
    TL = ascii.read(fLname, guess=False, delimiter=' ')

#   plot best and spread for BHMF/QLF at z
    y_best, med, spread = TM['y_best'], TM['med'], TM['spread']
    axm.fill_between(TM['Mx'],med-spread,med+spread,color='C%d'%(z%10),alpha=0.3,label=r'_$1\sigma$ Posterior Spread')
    axm.plot(TM['Mx'], y_best,linewidth=lw, color='C%d'%(z%10), label='z=%d'%z)
    y_best, med, spread = TL['y_best'], TL['med'], TL['spread']
    y_noobs = y_best * corr_U14D20(M1450)
    axl.fill_between(TL['M1450'],med-spread,med+spread,color='C%d'%(z%10),alpha=0.3,label=r'_$1\sigma$ Posterior Spread')
    axl.plot(TL['M1450'], y_best, linewidth=lw,color='C%d'%(z%10), label='z=%d'%z)
    axl.plot(TL['M1450'], y_noobs, '--',linewidth=lw,color='C%d'%(z%10), label='z=%d'%z)
    # QLF add detection limit
    # area same, Roman_deep deeper than Euclid_deep
    for ilabel in range(3):
        label = labels[ilabel]
        deep_mag = M_absolute(Depth[label],z)
        Vc_cover = Vc(Area[label],z,1)/1e9
        detection = [np.linspace(-32,deep_mag,num=50),np.linspace(1./Vc_cover,1e10,num=50)]
        # horizontal
        axl.plot(detection[0],1./Vc_cover*np.ones(50),styles[0],linewidth=1,color=colors[ilabel])#'C%d'%(z%10))
        # vertical
        axl.plot(deep_mag*np.ones(50),detection[1],styles[0],linewidth=1,color=colors[ilabel])#'C%d'%(z%10))
        if ilabel==0 and z==6:
            print(label+' depth mag={:.1f}, Vc={:.1e}, 1/Vc={:.1e}'.format(deep_mag,Vc_cover,1./Vc_cover))
            print('faintest abundance: ',np.max(y_best[TL['M1450']<deep_mag]))
            print('integrate QLF: ',np.nansum(y_best[TL['M1450']<deep_mag]*dmag))
            print('Vc_cover: ',Vc_cover)
        N_detect[label].append(np.nansum(y_best[TL['M1450']<deep_mag]*dmag*Vc_cover))

    index = np.where(TM['Mx']>Mmin)
    Mx = TM['Mx'][index]
    My_best,My_med,My_std = TM['y_best'][index],TM['med'][index],TM['spread'][index]
    dlog10Mx = np.log10(TM['Mx'][1]/TM['Mx'][0])
    rhoM.append( np.nansum(Mx*My_best*dlog10Mx) )
    # neglect My_min<0
    My_min = My_med-My_std
    rhoM_min.append( np.nansum((Mx*My_min*dlog10Mx)[My_min>0]) )
    My_max = My_med+My_std
    rhoM_max.append( np.nansum((Mx*My_max*dlog10Mx)) )

    index = np.where(TL['M1450']<Lmax)
    Ly_best,Ly_med,Ly_std = TL['y_best'][index],TL['med'][index],TL['spread'][index]
    nL.append( np.nansum(Ly_best*dmag) )
    Ly_min = (Ly_med-Ly_std)[(Ly_med-Ly_std)>0]
    nL_min.append( np.nansum(Ly_min*dmag) )
    Ly_max = (Ly_med+Ly_std)
    nL_max.append( np.nansum(Ly_max*dmag) )

    axl.set_xlim(np.max(TL['M1450']),np.min(TL['M1450']))
    axl.set_xlim(-18,-30)
    axl.set_ylim(1e-3,1e3)
    # axl.grid(1)
    axl.set_yscale('log')
    # axl.set_xlabel(r'$\mathrm{M_{1450}}$',fontsize=fslabel)
    # axl.set_ylabel(r'$\mathrm{\Phi~(Gpc^{-3}mag^{-1})}$',fontsize=fslabel)
    axl.tick_params(labelsize=fstick)
    # plt.tight_layout()
    figl.savefig(prex+'LF_spreadz%d_.png'%z,dpi=300,bbox_inches='tight')

axm.set_xlim(1e7,1e10); axm.set_xscale('log')
axm.set_ylim(1e-10,1e-4); axm.set_yscale('log')
axm.legend(fontsize=fslegend)
# axm.set_xlabel(r'$\mathrm{M_{BH}~(M_\odot)}$',fontsize=fslabel)
# axm.set_ylabel(r'$\mathrm{\Phi_M~(Mpc^{-3}dex^{-1})}$',fontsize=fslabel)
axm.tick_params(labelsize=fstick)
figm.savefig(prex+'MF_spread.png',dpi=300,bbox_inches='tight')

# minimum integrant, may be 0
print('rhoM_min',rhoM_min)
print('nL_min',nL_min)

print('N_detect:',N_detect)
# ----------------------------------------------------------------
#                     fitting rho v.s. z
# ----------------------------------------------------------------
#  fitting function, rho=rho_0*10^(k(z-6))
def func(z, logrho0, k):
    return logrho0 + k*(z-6.)

logrhoM,logrhoM_min,logrhoM_max =  np.log10(rhoM), np.log10(rhoM_min), np.log10(rhoM_max)
lognL,lognL_min,lognL_max =  np.log10(nL), np.log10(nL_min), np.log10(nL_max)

# logrho_min = -inf: logrho_err=2.*(logrhoM_max[i]-logrhoM[i])
logrhoM_err=[logrhoM_max[i]-logrhoM_min[i] if not np.isinf(logrhoM_min[i]) else 2.*(logrhoM_max[i]-logrhoM[i])for i in range(len(zs))]
lognL_err=[lognL_max[i]-lognL_min[i] if not np.isinf(lognL_min[i]) else 2.*(lognL_max[i]-lognL[i])for i in range(len(zs))]

print('logrhoM ',*logrhoM,sep=",")
print('logrhoM_err ',*logrhoM_err,sep=",")
print('lognL ',*lognL,sep=",")
print('lognL_err ',*lognL_err,sep=",")
# above as input for fit_rhok.ipynb

# fit M
popt, pcov = curve_fit(func, zs, logrhoM, sigma=logrhoM_err)
logrho0, k = popt
perr = np.sqrt(np.diag(pcov)) # error of fitted para, but better use MCMC
print(popt,perr)

# fit L
popt, pcov = curve_fit(func, zs, lognL, sigma=lognL_err)
logn0, kl = popt
perr = np.sqrt(np.diag(pcov)) # error of fitted para, but better use MCMC
print(popt,perr)

# or use k, kl from fit_rhok.ipynb calculation
# ----------------------------------------------------------------

# plot
# figl, axl = plt.subplots(num='LF',figsize=(10, 10))
figl = plt.figure(num='LF',figsize=(10, 10))
axl = figl.add_subplot(1, 1, 1)
z_ = np.arange(5.6,7.5,.1)
dash_k = 0.39*pow(10., -.78*(z_-6.7))
axl.plot(z_,dash_k, '--',c='black')
axl.errorbar(6.,1.33,yerr=0.33,fmt='D',elinewidth=2,capsize=4,label='Jiang+ 2016',
            c='red',ms=10,mfc='none')
axl.errorbar(6.7,0.39,yerr=0.11,fmt='s',elinewidth=2,capsize=4,label='Wang+ 2019',
            c='red',ms=10,mfc='none')

for i in range(len(zs)):
    z = zs[i]
    axl.errorbar(z,nL[i],
        yerr=[[nL[i]-nL_min[i]],[nL_max[i]-nL[i]]],
        fmt='o',ms=10,elinewidth=2,capsize=4,
        color='C%d'%(z%10),label='z=%d'%z,mfc='none')
# axl.plot(zs,nL,c='grey')
axl.plot(zs, pow(10., logn0+kl*(zs-6)),'--',c='grey')
# axl.text(6,np.median(nL),'log n0: {:.2f}, kL: {:.2f}'.format(logn0,kl))
print('log n0: {:.2f}, kL: {:.2f}'.format(logn0,kl))
axl.set_xlim(5.5,10.5)
# axl.set_ylim(1e-10,1e-4); 
axl.set_yscale('log')
axl.legend(fontsize=fslegend)
# axl.set_xlabel(r'$\mathrm{M_{BH}~(M_\odot)}$',fontsize=fslabel)
# axl.set_ylabel(r'$\mathrm{\Phi_M~(Mpc^{-3}dex^{-1})}$',fontsize=fslabel)
axl.tick_params(labelsize=fslabel)
figl.savefig(outLprex+'nL.png',dpi=300,bbox_inches='tight')


# figm, axm = plt.subplots(num='MF',figsize=(10, 10))
figm = plt.figure(num='MF',figsize=(10, 10))
axm = figm.add_subplot(1, 1, 1)
for i in range(len(zs)):
    z = zs[i]
    axm.errorbar(z, rhoM[i],
        yerr=[[rhoM[i]-rhoM_min[i]],[rhoM_max[i]-rhoM[i]]],
        fmt='o',elinewidth=2,capsize=4,
        color='C%d'%(z%10),label='z=%d'%z)
# axm.plot(zs,rhoM,c='grey')
axm.plot(zs, pow(10., logrho0+k*(zs-6)),'--',c='grey')
# axm.text(6,np.median(rhoM),'log rho0: {:.2f}, kM: {:.2f}'.format(logrho0,k))
print('Mmin={:.1e}, log rho0: {:.2f}, kM: {:.2f}'.format(Mmin,logrho0,k))
axm.set_xlim(5.5,10.5)
# axm.set_ylim(1e-10,1e-4); 
axm.set_yscale('log')
axm.legend(fontsize=fslegend)
# axm.set_xlabel(r'$\mathrm{M_{BH}~(M_\odot)}$',fontsize=fslabel)
# axm.set_ylabel(r'$\mathrm{\Phi_M~(Mpc^{-3}dex^{-1})}$',fontsize=fslabel)
axm.tick_params(labelsize=fslabel)
figm.savefig(outMprex+'rhoM.png',dpi=300,bbox_inches='tight')