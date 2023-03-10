# Gauss prior for parameter l_cut and a;

from PYmodule import *
from PYmodule.models import *
# python mpi
import sys 
sys.path.append('/gpfs/share/home/1801110214')
from mpi4py import MPI


# mpiexec -n 2 python -Wi MLF_paraSurvey.py

t1 = time.time()


logds = np.logspace(-2.5,-.5,num=10)
ls = np.arange(0.1,1.4,0.1)
aas = np.arange(-0.3,0.31,.1)
ts = np.arange(1,101,5)

# compare with try1.py
ts = [10]
logds = [-2.,-3.,-4.,-2.,-3.,-4.,-2.,-3.,-4.,-4.]
ls = [1.]
aas = [0.2]


N1=len(ts); N2=len(logds); N3=len(ls); N4=len(aas)
N_tot=N1*N2*N3*N4

# N_tot=5720; Ncore=220; Ntot//Ncore=26; Ntot%Ncore=0
# print(Ntot//Ncore,Ntot%Ncore); exit(0)

# def probab_anytau(para,like):
#     lp = lnprior_anytau(para)
#     prob = like + lp if np.isfinite(lp) else -np.inf
#     return prob

def chi2_prior(para):
    t_life, logd_fit, l_cut, a = para
    # t: uniform, no bound; fitting used 1e1<t_life<200.
    # logd: 
    # return -2*ln(P_prior)
    sigma_a = 1
    a_mean = 0
    sigma_l = .3
    l_mean = .5
   
    if -4<=logd_fit<=-.3 and l_cut>0.1:
        # return 0.0 - 0.5*((l_cut-l_mean)/sigma_l)**2 - 0.5*((a-a_mean)/sigma_a)**2 #- 0.5*((logd_fit+3.)/.1)**2
        if logd_fit< -3.:
            return ((l_cut-l_mean)/sigma_l)**2 + ((a-a_mean)/sigma_a)**2 + ((logd_fit+3.)/.1)**2
        else:
            return ((l_cut-l_mean)/sigma_l)**2 + ((a-a_mean)/sigma_a)**2
    else:
        return np.inf

def main():
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    
    N_each = N_tot//size

#     print('i1,i2,i3,i4={:d},{:d},{:d},{:d}'.format(i0,i1,i2,i3))
# while i<N_tot:
#     print('i1,i2,i3,i4={:d}'.format(i))
#     i += 1

    b = []
    for i in range(rank*N_each,(rank+1)*N_each):
        i1 = i//N4//N3//N2%N1
        i2 = i//N4//N3%N2
        i3 = i//N4%N3
        i4 = i%N4
        t_life = ts[i1]
        logd_fit, l_cut, a = logds[i2], ls[i3], aas[i4]
        para = np.array([t_life,logd_fit,l_cut,a])
        # like = lnlike(para)
        # prob = probab_anytau(para,like)
        # b.append( np.array([t_life,logd_fit,l_cut,a,-2.*like,-2.*prob]) )
        chi2_z4 = model(para,z=4)
        chi2_z5 = model(para,z=5)
        b.append( np.array([t_life,logd_fit,l_cut,a,chi2_z4,chi2_z5]) )

        # print('rank={:d}, i1,i2,i3,i4={:d},{:d},{:d},{:d}'.format(rank,i1,i2,i3,i4))


    b=np.array(b)
    sendbuf = b
    recvbuf = None

    if (rank==0):
        recvbuf = np.empty([N_each*size, len(para)+2])
    comm.Gather(sendbuf, recvbuf, root=0)

    if (rank==0):
        for i in range(N_each*size,N_tot):
            i1 = i//N4//N3//N2%N1
            i2 = i//N4//N3%N2
            i3 = i//N4%N3
            i4 = i%N4
            t_life = ts[i1]
            logd_fit, l_cut, a = logds[i2], ls[i3], aas[i4]
            para = np.array([t_life,logd_fit,l_cut,a])
            # like = lnlike(para)
            # prob = probab_anytau(para,like)
            # recvbuf = np.concatenate((recvbuf,np.array([[t_life,logd_fit,l_cut,a,-2.*like,-2.*prob]])),axis=0)
            chi2_z4 = model(para,z=4)
            chi2_z5 = model(para,z=5)
            recvbuf = np.concatenate((recvbuf,np.array([[t_life,logd_fit,l_cut,a,chi2_z4,chi2_z5]])),axis=0)
        np.savetxt('./chi2z4z5_Ntot{:d}Ncore{:d}.txt'.format(N_tot,size), recvbuf, delimiter=',', fmt='%10.3e')
        print('time=',time.time()-t1)
    

main()