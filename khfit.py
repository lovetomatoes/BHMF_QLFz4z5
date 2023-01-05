from PYmodule import *

t1 = time.time()


ftrprex = '/Users/wli/treefiles1e11/'
ftrprex = '/Users/wli/Desktop//Projects/treefiles/'
ftrprex = '/Users/wli/treefiles1e13/'

Ntree = 10000

ks = []
for i in range(Ntree):
    ftree = ftrprex+'tree_%dmer'%i
    T = ascii.read(ftree, guess=False, delimiter=' ')
    # print(T.info)
    z0 = T['z'][0]; M0 =  T['Mh_Ms'][0]
    ks.append( (np.log10(M0)-11.)/(z0-6.) ) 

print(np.mean(ks))
print(np.median(ks))
print(np.std(ks))

# 1e11: -0.16310373802687905; -0.1613976260263008; 0.01890051942654577
# 1e12: -0.14851385449635732; -0.14698968274194465; 0.016083836537797068
# 1e13: -0.1475742026265627; -0.1458100691025007; 0.01616766306666514