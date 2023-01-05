from PYmodule import *

def integral(a,x,x0):
    if a>0:
        return gamma(a)*(gammainc(a,x)-gammainc(a,x0))
    elif a ==0:
        a= 1e-8
        return integral(a,x,x0)
    else:
        return 1./a * (integral(a+1,x,x0) + pow(x,a)*np.exp(-x)-pow(x0,a)*np.exp(-x0))

def integral_toinf(a,x0):
    if a>0:
        return gamma(a)*gammaincc(a,x0)
    elif a==0:
        a = 1e-8
        return integral_toinf(a,x0)
    else:
        return 1./a * (integral_toinf(a+1,x0)-pow(x0,a)*np.exp(-x0))

dlogx = .0001
logx_min = np.log10(x0); logx_max = 2.+dlogx

def P_left2d(a,l2d):
    # integration of dP ~ x^a exp(-x) dlogx; normalized by integral over (-2,1)
    x = pow(10., np.arange(logx_min,logx_max,dlogx))
    m,n = l2d.shape
    l = l2d.reshape(-1)
    xx, ll = np.meshgrid(x,l)
    zz = (xx<ll)*(logx_min<=np.log10(ll)) * pow(xx,a)*np.exp(-xx) * dlogx
    # print(xx.shape)
    return np.sum(zz,axis=1).reshape(m,n)


def P_tot(a):
    logls = np.arange(logx_min,logx_max,dlogx)
    ls = pow(10., logls)
    return np.sum(pow(ls,a) * np.exp(-ls) * dlogx)

def P_left(a,l):
    # integration of dP ~ x^a exp(-x) dlogx; normalized by integral over (-2,1)
    x = pow(10., np.arange(logx_min,logx_max,dlogx))
    xx, ll = np.meshgrid(x, l)
    # zz = (xx<ll)*np.logical_and(logx_min<=np.log10(ll),np.log10(ll)<=logx_max) * pow(xx,a)*np.exp(-xx) * dlogx
    zz = (xx<ll)*(logx_min<=np.log10(ll)) * pow(xx,a)*np.exp(-xx) * dlogx
    return np.sum(zz,axis=1)

def P_left_norm(a,l):
    return P_left(a,l)/P_tot(a)

""" 
def P_tot(a):
    logx = logx_min
    P_tot = 0
    while logx<1:
        x = pow(10., logx)
        P_tot += pow(x,a)*np.exp(-x) * dlogx
        logx += dlogx
    # print('P_tot',P_tot)
    return P_tot

def P_left(a,l):
    logl = np.log10(l)
    if logx_min <= logl <= logx_max:
        logx = logx_min
        P = 0
        while logx<logl:
            x = pow(10., logx)
            P += pow(x,a)*np.exp(-x) * dlogx
            logx += dlogx
        return P
    else:
        return np.nan
 """