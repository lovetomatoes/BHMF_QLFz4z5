# fit paras: t_life, l_cut, a
# Gauss prior for parameter l_cut and a

from PYmodule import *
from PYmodule.MLF3p import *
from PYmodule.models import *
from emcee import EnsembleSampler as EnsembleSampler
import corner
import os
os.environ["OMP_NUM_THREADS"] = "1"
from schwimmbad import MPIPool

t1 = time.time()

# initial guess
# t_life, l_cut, a = 80, 0.5, 0. # x0=0.01, f_seed = 1
# t_life, l_cut, a = 30, 1.5, 0. # x0=0.001, f_seed = 1

# lambda_0=0.01
t_life, l_cut, a = 30, 1., -.2 # f_seed = 1.

initial = np.array([t_life,l_cut,a])

ndim = len(initial)
nwalkers = 100
nsteps = 10000
rball = 1e-4

prex='../3p/MLF3prange2_Gs.3_l{0:.2f}_a{1:.2f}_sl{2:.2f}_sa{3:.2f}_fseed_{4:.1e}_ns_{5:.0e}'.format(l_mean,a_mean,sigma_l,sigma_a,f_seed,nsteps)
# LFbin, LFcur, MF1e8 

fname =prex+'.h5'

# nsteps = 10000
# prex += '_xu'

with MPIPool() as pool:
    if not pool.is_master():
        pool.wait()
        sys.exit(0)

    np.random.seed(42)

    backend = emcee.backends.HDFBackend(fname)

    # --initialize-- clear  output and reset
    backend.reset(nwalkers, ndim)
    p0 = [initial + rball*np.random.randn(ndim) for i in range(nwalkers)]
    sampler = EnsembleSampler(nwalkers, ndim, lnprobab, pool=pool, backend=backend)
    sampler.run_mcmc(p0, nsteps, progress=True)

    # # --resume-- 
    # sampler = EnsembleSampler(nwalkers, ndim, lnprobab, pool=pool, backend=backend)
    # sampler.run_mcmc(None, nsteps, progress=True)


fig, axes = plt.subplots(ndim+1, figsize=(10, 7), sharex=True)
samples = sampler.get_chain()
# print('samples',samples)
probs = sampler.get_log_prob()
# print('probs',probs)

labels = ['t_life', 'l_cut', 'a', 'prob']
for i in range(ndim):
    ax = axes[i]
    ax.plot(samples[:, :, i], "k", alpha=0.3)
    ax.set_xlim(0, len(samples))
    ax.set_ylabel(labels[i])
    ax.yaxis.set_label_coords(-0.1, 0.5)

i += 1
ax = axes[i]
ax.plot(probs[:, :], "k", alpha=0.3)
ax.set_xlim(0, len(samples))
ax.set_ylabel(labels[i])
ax.yaxis.set_label_coords(-0.1, 0.5)
axes[-1].set_xlabel("step number")

plt.savefig(prex+'_chain.png')


samples = sampler.flatchain
probs = sampler.flatlnprobability
print('len of samples:', len(samples))
theta_max = samples[np.argmax(probs)]
print('initial paras: t_life={0:.1e}, l_cut={1:.1f}, a={2:.1f}, x0{3:.0e}, prob{4:.1e}'.format(t_life,l_cut,a,x0,probs[0]))
print('Gaussian scatter sigma_l,sigma_a:',sigma_l,sigma_a)
print('best paras:',labels,theta_max,np.max(probs))

all_samples = np.concatenate(
    (samples, probs[:, None]), axis=1
)
fig = corner.corner(all_samples,show_titles=True,labels=labels,plot_datapoints=True,quantiles=[0.16, 0.5, 0.84])
plt.savefig(prex+'_corner.png')

# pyname = sys.argv[0][:-3] # current .py file name
print('nwalkers={0:d}, nsteps={1:.0e}, rball={2:.0e}'.format(int(nwalkers),int(nsteps),rball))

print(
    "Mean acceptance fraction: {0:.3f}".format(
        np.mean(sampler.acceptance_fraction))
)

tau = sampler.get_autocorr_time(tol=1)
print('Autocorrelation timescale: ',tau)

print('running time: {0:.1f} hrs'.format((time.time()-t1)/3600))

ndraw = 100
fig, axes = plt.subplots(1,2, figsize=(12, 6),dpi=400)
ax = axes[0]; curve_name = 'MF'
best_model = model(theta_max)
xs = best_model['M_BH']
y_data = best_model[curve_name+'_data']
y_best = best_model[curve_name]
ax.plot(xs, y_data, label='data')
draw = np.floor(np.random.uniform(0,len(samples),size=ndraw)).astype(int)
thetas = samples[draw]
for i in thetas:
    mod = model(i)[curve_name]
    ax.plot(xs, mod, c='grey',label='_',alpha=.2)
ax.plot(xs, y_best, c='C1', label='Highest Likelihood Model')
ax.set_xlim(1e7,1e10); ax.set_xscale('log')
ax.set_ylim(1e-10,1e-4); ax.set_yscale('log')
ax.legend()
ax = axes[1]; curve_name = 'LF'
xs = best_model['M1450']
y_data = best_model[curve_name+'_data']
y_best = best_model[curve_name]
ax.scatter(xs, y_data, label='_')
for i in thetas:
    mod = model(i)[curve_name]
    ax.plot(xs, mod, c='grey',label='_',alpha=.2)
ax.plot(xs, y_best, c='C1', label='_')
ax.text(-26,10,
r'$t_{life}=$'+'{0:.1e}Myr\n'.format(theta_max[0])+r'$\lambda_{cut}=$'+'{0:.2f}\n'.format(theta_max[1])\
+r'$\alpha=$'+'{0:.1e}\n'.format(theta_max[2])
)
ax.set_xlim(-22,-29)
ax.set_ylim(1e-2,1e2)
ax.set_yscale('log')
ax.legend()
plt.savefig(prex+'_spread.png')