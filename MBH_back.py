from PYmodule import *
# Phiz6 paras: f_0, t_life, lbd ~Schechter(l_cut,a)

# quasar J1205 in Onoue 2019
M1 = 2.2e9
z1 = 6.8
l1 = .16 #由于δ, M1并非以λ_Edd=0.16 Eddington trace back
tz = t_from_z(z1)

# print( z_tH(tz/Myr) ); exit(0)
z0 = 30
Dt = t_from_z(z1) - t_from_z(z0)

nt = 100
ts = np.linspace(0.01, 1, nt)*Dt

t_life, d_fit = 100 ,  0.25
t_life *= Myr

Nt = int(Dt//t_life)
ls = np.append([l1],np.random.gamma(shape=a, scale=l_cut, size=Nt))
# ls = 2.*np.ones(Nt) # Eddington accretion when δ<<1
ls = np.array([ls[int(ts[i]//t_life)] for i in range(nt)])
M0s = [M1]

print(Nt, l_cut, a, ls)
print('ls[0]={0:e}\nM0={1:e}'.format(ls[0], M0M1(M1, ls[0], t_life, d_fit)))

for i in range(1,nt):
    dt = ts[i] - ts[i-1]
    M0 = M0M1(M1, ls[i], dt, d_fit)
    M0s.append(M0)
    M1 = M0


prex = '../MBH_back'


plt.figure(figsize=(10,10),dpi=400)
plt.plot(ts/Myr,M0s)
plt.xscale('log'); plt.yscale('log')
# plt.xlim(-21,-29)
# plt.ylim(1e-2,1e2)
plt.grid(True)
plt.xlabel(r'$\mathrm{t_{b}}$',fontsize=fslabel)
plt.ylabel(r'$\mathrm{M_{BH}}$',fontsize=fslabel)
plt.legend(loc='best',fontsize=fslabel)
plt.savefig(prex+'_t.png')

zs = z_tH((tz-ts)/Myr)
plt.figure(figsize=(10,10),dpi=400)
plt.scatter(zs,M0s)
plt.yscale('log')
plt.xlim(30,6)
plt.ylim(1e2,1e10)
plt.grid(True)
plt.xlabel(r'$\mathrm{t_{b}}$',fontsize=fslabel)
plt.ylabel(r'$\mathrm{M_{BH}}$',fontsize=fslabel)
plt.legend(loc='best',fontsize=fslabel)
plt.savefig(prex+'_z.png')