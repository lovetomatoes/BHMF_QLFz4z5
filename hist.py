from PYmodule import *
# analyse BH ERDF from sampling in MBH_evol

# M1450=-26, Lbol=1e47; M1450=-22, Lbol=2e45
L_limit = 1e46

lbin = np.linspace(-2,1.2,num=20)
hist_tot = np.zeros(len(lbin)-1)
Nfiles = int(10)
for ii in range(1):
    Tz6 = ascii.read('../BHatz6.dat', guess=False, delimiter=' ')
    # print(Tz6.info)
    Mmin = np.min(Tz6['M1']); Mmax = np.max(Tz6['M1'])
    # print('min max of M1 {0:.1e}, {1:.1e}'.format(Mmin,Mmax))
    lmin = np.min(Tz6['ls']); lmax = np.max(Tz6['ls'])
    print('min max of ls {0:.1e}, {1:.1e}'.format(lmin,lmax))

    index = np.logical_and(Mmin<=Tz6['M1'],Tz6['M1']<=Mmax) # all
    M1 = Tz6['M1'][index]; L1 = Tz6['L1'][index]; ls = Tz6['ls'][index]

    index = np.where(L_limit<L1)
    M1_ = M1[index]; L1_ = L1[index]; ls_ = ls[index]
    hist, bin_edges = np.histogram(np.log10(ls_),bins=lbin,density=False)
    hist_tot += hist

    # plt.figure(figsize=(10,8),dpi=400)
    # plt.scatter( bin_edges[:-1],hist/len(ls))
    # # plt.yscale('log')
    # plt.ylim(.1/len(ls), 1)
    # plt.savefig('../hist_l_Llim%d.png'%ii)

# histogram file
ascii.write(Table([bin_edges[:-1],hist_tot/Nfiles/len(ls)]),'../hist_tot%d.dat'%(int(np.log10(L_limit))),
names=['log_l','hist'],formats={'log_l':'10.2f','hist':'10.2e'},overwrite=True)