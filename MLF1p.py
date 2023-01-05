from PYmodule import *
from PYmodule.MLF1p import *
from PYmodule.models import *
from emcee import EnsembleSampler as EnsembleSampler
import corner
import os
os.environ["OMP_NUM_THREADS"] = "1"
from schwimmbad import MPIPool


# initial guess
t_life = 80.
initial = np.array([t_life])

ndim = len(initial)
nwalkers = 100
nsteps = 5000
rball = 1e-4

prex='../1p/d_{0:d}_l{1:.1f}_a{2:.1f}'.format(d_fit,l_cut,a)
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

labels = ['t_life', 'prob']
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
print('initial paras: t_life, d_fit, logM0, l_cut, a, prob',labels,t_life,d_fit,logM0,l_cut,a,probs[0])
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
print(tau)

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
ax.text(-26,3e1,r'$t_{life}=$'+'{0:.1e}Myr\n'.format(theta_max[0]))
ax.set_xlim(-22,-29)
ax.set_ylim(1e-2,1e2)
ax.set_yscale('log')
ax.legend()
plt.savefig(prex+'_spread.png')