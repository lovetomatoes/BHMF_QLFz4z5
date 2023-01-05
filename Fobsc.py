from PYmodule import *
from PYmodule.l_intg import *
from scipy.stats import norm, uniform

for m in np.arange(-29,-21,1):
    print('%.0e'%Lbol_M1450(m))

Lbols = np.logspace(44,50,num=100)
M1450s= M1450_Lbol(Lbols)
Lxs = Lbols/K_AVE20(Lbols)
print('faintest: ',Lbols[0],Lxs[0],M1450s[0],
'brightest: ',Lbols[-1],Lxs[-1],M1450s[-1])

f_U = 1 - 1./corr_U14D20(M1450s)
f_M = 0.56+1/np.pi*np.arctan((43.89-np.log10(Lxs))/0.46)
f_M[Lxs<1e43] = f_M[np.argmax(Lxs>1e43)]

def phi(Lx):
    phi_min, phi_max, beta, a1 = .2, .84, .24, .48
    phi_4375_0 = .43
    phi_4375_z = phi_4375_0*pow(1+2.,a1)
    print('phi_4375_z=%.2f'%phi_4375_z)
    if isinstance(Lx,float):
        return min( phi_max, max(phi_4375_z - beta*(np.log10(Lx)-43.75), phi_min))
    else:
        return np.array([min( phi_max, max(phi_4375_z - beta*(np.log10(Lx[i])-43.75), phi_min)) for i in range(len(Lx))])

fig, ax = plt.subplots(figsize=(10,10),dpi=400)
c_M = corr_M14D20(M1450s)
c_U = corr_U14D20(M1450s)
ax.plot(Lbols,c_M)
ax.plot(Lbols,c_U)
ax.set_xscale('log')
ax.set_yticks(np.arange(0,5,.5))
ax.set_xlim(np.min(Lbols),np.max(Lbols))
ax.grid(True)
ax2 = ax.twiny()
ax2.plot(M1450s, c_M); ax2.plot(M1450s, c_U)
ax2.plot(M1450s, c_M/c_U)
ax2.set_xlim(np.max(M1450s),np.min(M1450s))
ax2.set_xlabel('M1450')
fig.savefig('../corr_M1450.png')


fig, ax = plt.subplots(figsize=(10,10),dpi=400)
ax.plot(Lbols,f_M)
ax.plot(Lbols,f_U)
ax.set_xscale('log')
ax.set_xlim(np.min(Lbols),np.max(Lbols))
ax2 = ax.twiny()
ax2.plot(M1450s, f_M); ax2.plot(M1450s, f_U)
ax2.plot(M1450s, f_M/f_U)
ax2.set_xlim(np.max(M1450s),np.min(M1450s))
ax2.set_xlabel('M1450')
fig.savefig('../fobsc_M1450.png')

fig, ax = plt.subplots(figsize=(10,10),dpi=400)
ax.plot(np.log10(Lxs),f_M)
ax.plot(np.log10(Lxs),f_U)
ax.plot(np.log10(Lxs),phi(Lxs))
ax.plot(np.log10(Lxs),f_M/f_U)
print(np.min(f_M/f_U))
ax.set_xlim(np.min(np.log10(Lxs)),46)
ax.set_xlabel('logLxs')
ax.set_ylim(0,1)
ax.grid(True)
fig.savefig('../fobsc_Lx.png')

# see difference of U14, M14; correct Matsu18
fig, ax = plt.subplots(figsize=(10,10),dpi=400)
ax.scatter(bin_cen,Phi_obs)
ax.plot(M1450,1e9*LF_M1450(M1450),label='obs fit')
ax.plot(M1450,1e9*LF_M1450(M1450)*corr_M14D20(M1450),label='M14 corr')
ax.plot(M1450,1e9*LF_M1450(M1450)*corr_U14D20(M1450),label='U14 corr')
ax.set_xlim(-21,-30)
ax.set_yscale('log')
ax.grid(True)
ax.legend(loc='best',fontsize=25)
fig.savefig('../LF_fUfM.png')

exit(0)

ratio = np.array([(f_M/f_U)[i] if (f_M/f_U)[i]>1 else (f_U/f_M)[i]for i in range(len(f_M))])
print(np.max(ratio))
print('-20',f_U[np.argmax(M1450s<-20)])
print('-24',f_U[np.argmax(M1450s<-24)])
print('-28',f_U[np.argmax(M1450s<-28)])

fig, ax = plt.subplots(figsize=(10,10),dpi=400)
ax.plot(x,1./(1-f_M))
ax.plot(x,1./(1-f_U))
ax.set_xscale('log')
ax.set_xlim(np.min(Lbols),np.max(Lbols))
ax2 = ax.twiny()
ax2.plot(M1450s,1./(1-f_M)); ax2.plot(M1450s, 1./(1-f_U))
ax2.set_xlim(np.max(M1450s),np.min(M1450s))
ax2.set_xlabel('M1450')
plt.savefig('../corr.png')

fig, ax = plt.subplots(constrained_layout=True)
x = np.arange(0, 360, 1)
y = np.sin(2 * x * np.pi / 180)
ax.plot(x, y)
ax.set_xlabel('angle [degrees]')
ax.set_ylabel('signal')
ax.set_title('Sine wave')
def deg2rad(x):
    return x * np.pi / 180
def rad2deg(x):
    return x * 180 / np.pi
secax = ax.secondary_xaxis('top', functions=(deg2rad, rad2deg))
secax.set_xlabel('angle [rad]')
plt.savefig('../angle.png')



# closed = np.all(np.isclose(lhs,rhs,rtol=1e-2))
# closed_ = np.all(np.isclose(lhs,nume,rtol=1e-2))
# print('closed=',closed,closed_)

fig, ax = plt.subplots(1, 1, dpi=400)
x = np.linspace(uniform.ppf(0.01),
                uniform.ppf(0.99), 100)
y = uniform.pdf(x)
ax.plot(x, y,
       'r-', lw=5, alpha=0.6, label='uniform pdf')
rv = uniform()
ax.plot(x, rv.pdf(x), 'k-', lw=2, label='frozen pdf')
r = uni
x = np.arange(-21,-29,-1)
y = corr_U14D20(x)
print(x,y)
plt.figure(figsize=(10,10),dpi=400)
z = 4 + 2.5*np.tanh(.5*(x+21))
plt.plot(x,y)
plt.plot(x,z)
plt.ylim(1,7)
plt.savefig('../corr.png')
exit(0)