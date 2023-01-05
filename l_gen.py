from PYmodule import *
from PYmodule.l_intg import *

t1 = time.time()

l_cut = .5
a = 0.5
x0 = lambda_0/l_cut
I_toinf = integral_toinf(a,x0)

x = np.logspace(np.log10(x0),1.2,num=200)
Pa_x = integral(a,x,x0)/I_toinf
# print(x,Pa_x)

N_BH = int(1e5)

uniform_a = np.random.uniform(size=N_BH)
xx,yy = np.meshgrid(uniform_a,Pa_x)
ls = x[(yy>xx).argmax(axis=0)] * l_cut

print('np.min(ls)',np.min(ls))

abin = np.log10(x)
hist, bin_edges = np.histogram(np.log10(ls/l_cut),bins=abin,density=False)

plt.figure(figsize=(10,8),dpi=400)
plt.scatter( bin_edges[:-1],hist/len(ls))
print(np.sum(hist)/len(ls))
print(np.sum((Pa_x[1:]-Pa_x[:-1])))
plt.plot(np.log10(x[:-1]),(Pa_x[1:]-Pa_x[:-1]),c='C1')
plt.yscale('log')
plt.ylim(1e-8, 1e-1)
plt.savefig('../Plambda_res{0:d}Nsamp{1:d}.png'.format(len(x),int(np.log10(N_BH))))
print('time: ',time.time()-t1)

# time to save memory --> OMG also save time!

t1 = time.time()
uniform_a = np.random.uniform(size=N_BH)
ls = np.zeros(N_BH)
for i in range(N_BH):
    ls[i] = x[np.argmax(Pa_x>uniform_a[i])] # argmax: first True of a bool array

abin = np.log10(x)
hist, bin_edges = np.histogram(np.log10(ls),bins=abin,density=False)

plt.figure(figsize=(10,8),dpi=400)
plt.scatter( bin_edges[:-1],hist/len(ls))
print(np.sum(hist)/len(ls))
print(np.sum((Pa_x[1:]-Pa_x[:-1])))
plt.plot(np.log10(x[:-1]),(Pa_x[1:]-Pa_x[:-1]),c='C1')
plt.yscale('log')
plt.ylim(1e-8, 1e-1)
plt.savefig('../Plambda_indv_res{0:d}Nsamp{1:d}.png'.format(len(x),int(np.log10(N_BH))))
print('time: ',time.time()-t1)