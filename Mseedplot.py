from PYmodule import *

t1 = time.time()

N_Mh = 3 # 3 base halos: 1e11, 1e12, 1e13
Ntr = 10000
eta = 0.3
z = int(6)
tz = t_from_z(z)
alpha = 1.



T = Ts[0][0]
f_bsm = 1.
n_base = n_base[0]

Mseed = T['Mstar0']
tseed = T['t_col']/Myr
Tvir = Tv(T['Mh_col'],T['z_col'])

print('mean Dt:',np.mean((tz-T['t_col']))/Myr)
print('Tvir of collapse min, max = {:.1e}, {:.1e}'.format(np.min(Tvir),np.max(Tvir)))


figm = plt.figure(num='MF',figsize=(10, 10))
axm = figm.add_subplot(1, 1, 1)

cbar = axm.scatter(tseed,np.log10(Mseed),
    # color=my_cmap(np.log10(np.min(np.log10(Tvir)))+rescale(np.log10(Tvir))), 
    color=my_cmap(rescale(np.log10(Tvir))), 
    # color=my_cmap( np.log10(Tvir) ), 
    alpha = 0.1,label='_Myr')
# axm.plot(M_BH, TM['W10_MF'],c='C0')
figm.colorbar(cbar,ticks=[0,1,2,3,4,5])
# axm.set_xlim(1e7,1e10); axm.set_xscale('log')
# axm.set_ylim(1e-10,1e-4); axm.set_yscale('log')
axm.legend(fontsize=fslegend)
axm.tick_params(labelsize=fstick)
figm.savefig('../Mseed_t.png')