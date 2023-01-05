from PYmodule import *
from PYmodule.models_logd import *


f_seedlabel = 'f%d'%abs(int(np.log10(f_seed)))
f_seedlegend = r'$\mathrm{f_{seed}}=$'+str(f_seed)
fname = datapre+'/M0r8_'+f_seedlabel+'.h5'
prex = '../' + f_seedlabel

labels = [r'$\mathrm{\tau}$', r'$\log \delta$', r'$\lambda_0$', r'$\alpha$']

print('fname=',fname)
labels = [r'$\tau$', r'$\log \delta$', r'$\lambda_0$', r'$\alpha$']

reader = emcee.backends.HDFBackend(fname)
ndim = len(labels)

# tau = reader.get_autocorr_time(tol=1)
# tau = np.max(tau); print('max tau:',tau)
# Nburnin = int(3*tau)
# Nthin = int(tau/2)

# v.s. Nburnin=500, Nthin=1 no big difference
Nburnin = 0
Nthin = 1

samples = reader.get_chain(discard=Nburnin)
probs = reader.get_log_prob(discard=Nburnin)
print('len of samples:', len(samples))

samples = reader.get_chain(discard=Nburnin, thin=Nthin, flat=True)
probs = reader.get_log_prob(discard=Nburnin, thin=Nthin, flat=True)
print('len of extracted samples:', len(samples))
theta_max = samples[np.argmax(probs)]
print('best paras:',theta_max,np.max(probs))

# fig = corner.corner(samples,show_titles=True,title_kwargs={'fontsize':15},
# label_kwargs={'fontsize':20},max_n_ticks=4,labels=labels,plot_datapoints=True,
# quantiles=[0.16, 0.5, 0.84])
# axes = np.array(fig.axes).reshape((ndim, ndim))
# for i in range(ndim):
#     for j in range(ndim):
#         axes[i][j].tick_params(labelsize=12)
# plt.savefig(prex+'_corner.png',dpi=400,rasterized=True)
# exit(0)


ndraw = 60

best_model = model(theta_max)
draw = np.floor(np.random.uniform(0,len(samples),size=ndraw)).astype(int)
thetas = samples[draw]
model_thetas = [model(theta_i) for theta_i in thetas]

curve_name = 'MF'
xs = best_model['M_BH'][::int(N_mf/100)] # ~100 points
y_data = best_model[curve_name+'_data'][::int(N_mf/100)]
y_logdata = np.log10(y_data)
y_best = best_model[curve_name][::int(N_mf/100)]
y_err = best_model[curve_name+'_data_err'][::int(N_mf/100)]

mod_list = []
for i in range(ndraw):
    mod = model_thetas[i][curve_name][::int(N_mf/100)]
    mod_list.append(mod)
    # ax.plot(xs, mod, c='grey',label='_',alpha=.3)
spread = np.std(mod_list,axis=0)
med_model = np.median(mod_list,axis=0)

fMname = f_seedlabel+'ndraw%d'%ndraw+curve_name
ascii.write(Table([xs,y_data,y_best,y_err,med_model,spread]),
            fMname,
            names=['xs','y_data','y_best','y_err','med_model','spread'],
            formats={'xs':'10.3e','y_data':'10.3e','y_best':'10.3e','y_err':'10.3e','med_model':'10.3e','spread':'10.3e'},
            overwrite=True)

curve_name = 'LF'
xs = best_model['M1450']
Nappend = len(best_model['M1450'])-len(best_model['M1450_data'])
x_data = np.append(best_model['M1450_data'],np.ones(Nappend))
y_data = np.append(best_model[curve_name+'_data'],np.ones(Nappend))
y_data_err = np.append(best_model[curve_name+'_data_err'],np.ones(Nappend))

y_best = best_model[curve_name]
# spread
mod_list = []
for i in range(ndraw):
    mod = model_thetas[i][curve_name]
    mod_list.append(mod)
spread = np.std(mod_list,axis=0)
med_model = np.median(mod_list,axis=0)

fLname = f_seedlabel+'ndraw%d'%ndraw+curve_name
ascii.write(Table([xs,x_data,y_data,y_data_err,y_best,med_model,spread]),
            fLname,
            names=['xs','x_data','y_data','y_data_err','y_best','med_model','spread'],
            formats={'xs':'10.3e','x_data':'10.3e','y_data':'10.3e','y_data_err':'10.3e','y_best':'10.3e','med_model':'10.3e','spread':'10.3e'},
            overwrite=True)