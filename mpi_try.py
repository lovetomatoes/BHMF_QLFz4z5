from mpi4py import MPI
import numpy as np
# mpiexec -n 4 python -Wi mpi_try.py 
# test mpi; in main: 1st column rank=0 takes the remaining largest t_life;

##  -------------   send & receive  --------------------
def try_():
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    sendbuf = np.zeros((1,3), dtype='i') + rank
    recvbuf = None

    print(rank,sendbuf) 
    if rank == 0:
        recvbuf = 100*np.ones([size, 3], dtype='i')

    comm.Gather(sendbuf, recvbuf, root=0)
    if rank == 0:
        print(recvbuf)
        recvbuf
        for i in range(size):
            assert np.allclose(recvbuf[i,:], i)
try_()


##  -------------   para survey   --------------------

ts = np.arange(5)
logds = np.arange(1)
ls = np.arange(2)
aas = np.arange(3)

N1=len(ts); N2=len(logds); N3=len(ls); N4=len(aas)
N_tot=N1*N2*N3*N4

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
        # b.append( np.array([t_life, logd_fit, l_cut, a, -2.*lnlike(para), -2.*lnprobab(para)]) )
        b.append( np.array([t_life, logd_fit, l_cut, a, 0., 1.]) )

        print('rank={:d}, i1,i2,i3,i4={:d},{:d},{:d},{:d}'.format(rank,i1,i2,i3,i4))


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
            # b.append( np.array([t_life, logd_fit, l_cut, a, -2.*lnlike(para), -2.*lnprobab(para)]) )
            recvbuf = np.concatenate((recvbuf,np.array([[t_life, logd_fit, l_cut, a, 0., 1.]])),axis=0)
        np.savetxt('test.out', recvbuf, delimiter=',', fmt='%10.3e')
    

main()