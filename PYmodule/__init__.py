import numpy as np
from ctypes import * # c 类型库
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.table import Table, vstack
import os
import sys
from matplotlib.ticker import LogLocator,LinearLocator,MultipleLocator,AutoLocator,FixedLocator
from matplotlib.legend_handler import HandlerLine2D, HandlerTuple

from scipy.stats import *
import time
import math
import numpy.ma as ma
from scipy import special
from scipy.special import gamma, gammainc, gammaincc
from emcee import EnsembleSampler
import corner
import emcee

z4figpre = '../z4/figs/'
z4datapre = '../z4/data/'
z5figpre = '../z5/figs/'
z5datapre = '../z5/data/'
z6figpre = '../z6/figs/'
z6datapre = '../z6/data/'
datapre = '../data/'
figpre = '../figs/'

d_fit = 0.
logM0 = 8.
l_cut = 1. # l_cut=2., l_cut' = l_cut/2; M=M_cut=1e7 grow as Eddington
a = -.7

lambda_0 = 0.01 # starting point of lambda; x0 now changable
x0 = 0.01 # integration of l/l_cut, starting point; previous context
t_life = 50
t_life = 1000

f_seed = .01
corr = 'U'

l_mean, a_mean = 0.6, 0.
sigma_l, sigma_a = .4, .3

typenames = ['H'+r'$_2$', 'H-H'+r'$_2$', 'H-H']
lfnames = {'4':'Akiyama_18','5':'Niida_20','6':'Matsuoka_18'}

global G, h0, H0, Omega_m0, Omega_L0, m_H, mu, Ms, pi, km, pc, Myr, alpha_T
G, c, k_B, m_H = 6.67408e-8, 2.9979245e10, 1.38064852e-16, 1.66053904e-24
pi = 3.141593
mu = 1.2
Ms = 2.e33
Lsun = 3.828e33
pc = 3.e18
Mpc = 1.e6*pc
km = 1.e5
yr = 365*24*3600
Myr = 1.e6*(365*24*3600)
Omega_m0 = 0.307
Omega_L0 = 1 - Omega_m0
h0 = .677
H0 = h0*100*km/Mpc


from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmo = FlatLambdaCDM(H0=h0*100 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=Omega_m0)

def M_absolute(m_apparent,z):
    dl = cosmo.luminosity_distance(z)
    return m_apparent+2.5*np.log10((10/(dl.value*(10**6)))**2*(1+z))

def m_apparent(M_absolute,z):
    dl = cosmo.luminosity_distance(z)
    return M_absolute-2.5*np.log10((10/(dl.value*(10**6)))**2*(1+z))

# comoving volumes
def Vc(A_deg2,z,dz):
    dp = cosmo.comoving_distance(z+dz/2.)
    dm = cosmo.comoving_distance(z-dz/2.)
    return 4*np.pi/3. * A_deg2/41253.* ( pow(dp.value,3) - pow(dm.value,3) )

Area = {'NIRCam_deep':0.013,'NIRCam_med':0.053,'Roman_deep':40,'Roman wide':2000,'Euclid_deep':40,'Euclid wide':15000}
Depth = {'NIRCam_deep':30.6,'NIRCam_med':29.7,'Roman_deep':29,'Roman wide':27,'Euclid_deep':26,'Euclid wide':23}
Area['LSST_deep'] = 18000
Depth['LSST_deep'] = 27.5

t_Edd = 1./(4*pi*G/.4/c)
fbol_1450 = 4.4

log10Ms = [11,12,13]
# n_base = [1.63,1.09e-01,4.02e-03,3.87e-05,1.07e-08]
n_base = [4.02e-03,3.87e-05,1.07e-08]
n_base = [9.87e-4, 6.15e-6, 8.92e-10] # from dNdlnM.ipynb & Jdist.cpp
Mhpres = [datapre+'1e11',datapre+'1e12',datapre+'1e13']

f_bsm = [.6,.4]
Nbsm = 2 # control bsm case to use

W37 = 1e44

alpha_T = 2.324e4

fstick = 20
fstxt = 20
fslabel = 23
fstitle = 20
fslegend = 20
lw = 2

# color map
my_cmap = plt.get_cmap("viridis")
rescale = lambda y: (y - np.min(y)) / (np.max(y) - np.min(y))

# -----------------  binned data in Gpc^-3 mag^-1     ---------------------
# -----------------  '6': z=6 Matsuoka 2018; Table 4     ---------------------
# -----------------  '5': z=5 Niida 2020; Table 3 excluding m_i<23.1    ---------------------
# -----------------  '4': z=4 Akiyama 2018;     ---------------------    
bin_cen = {'6':np.array([-29,-27.5,-26.75,-26.25,-25.75,-25.25,-24.75,-24.25,-23.75,-23.25,-22.75,-22]),
           '5':np.append(np.arange(-28.63,-26.,.25),np.arange(-27.07,-23.5,.5)),
           '4':np.arange(-28.875, -21.625, .25)
           }
bin_wid = {'6':np.array([2, 1, .5, .5, .5, .5, .5, .5, .5, .5, .5, 1]),
           '5':np.append(0.25*np.ones(len(np.arange(-28.63,-26.,.25))),
                         0.5*np.ones(len(np.arange(-27.07,-23.5,.5)))),
           '4':.25*np.ones(len(bin_cen['4']))
           }
bin_edg = {'6':np.append(bin_cen['6'] - bin_wid['6']/2., bin_cen['6'][-1]+bin_wid['6'][-1]/2.),
           '5':np.append(bin_cen['5'] - bin_wid['5']/2., bin_cen['5'][-1]+bin_wid['5'][-1]/2.),
           '4':np.append(bin_cen['4'] - bin_wid['4']/2., bin_cen['4'][-1]+bin_wid['4'][-1]/2.)
           }
Phi_obs = {'6':np.array([.0079, .242, .58, .9, 1.33, 4.6, 7.,  6.6, 8.3, 10.9, 23., 16.2]),
           '5':10.*np.array([.018,.018,.0092,.055,.083,.212,.277,.499,.628,.772,.639,
                             1.2,.6,1.8,2.98,7.78,5.39,10.7,12.5]),
           '4':10**np.array([-9.710,-9.534,-9.710,-9.057,-8.781,-8.534,-8.192,-7.966,-7.853,-7.642,
                             -7.552,-7.387,-7.077,-7.189,-6.899,-6.756,-6.745,-6.577,-6.487,-6.493,
                             -6.405,-6.341,-6.369,-6.297,-6.382,-6.298,-6.219,-6.088,-6.253]) * 1e9
           }
Phi_err = {'6':np.array([.0079, .061, .17, .32, .6, 1.2, 1.7, 2., 2.6, 3.6, 8.1, 16.2]),
           '5':10.*np.array([.013,.013,.0092,.023,.028,.044,.051,.068,.076,.087,.171,
                             1.58,1.39,1.75,2.01,2.81,2.46,3.2,3.4]),
           '4':np.array([0.138,0.169,0.138,0.292,0.402,0.534,0.791,1.026,1.169,1.490,1.657,2.545,
                         29.581,17.287,23.032,27.088,27.422,33.384,37.385,37.343,42.783,47.960,47.366,
                         51.538,46.855,52.208,59.367,71.892,80.733])
           }     


# !!!!!!!!!!!! eta_max should be used when 1/eta \simeq (1-eta)/eta breaks
eta_max = .5; eta_min = 0.057 # 0.057

# def M1M0(M0,dt,f_duty,mu_fit,eta8,delta):
#     M1 = 1e8*pow(mu_fit*f_duty*delta*dt/(eta8*t_Edd)+pow(M0/1e8,delta),1./delta)

#     ## eta mimicking Ueda14 empirical formula: eta = eta8*(M/M8)**delta
#     eta = eta8*pow(M0/1e8,delta)
#     ## M1: mass after growth following t^(1/delta) power
#     M1 = 1e8*pow(mu_fit*f_duty*delta*dt/(eta8*t_Edd)+pow(M0/1e8,delta),1./delta)

#     ## if eta > maximum -> use Eddington accretion -- M(t) \propto M0*exp(...)
#     eta = ma.masked_greater(eta, eta_max)
#     M1[eta.mask] = M0[eta.mask]*np.exp(mu_fit*f_duty*dt/(eta_max*t_Edd))
#     # i = ma.argmax(eta)
#     # M1[eta.mask] = M0[eta.mask]/M0[i]*M1[i]
#     eta[eta.mask] = eta_max

#     ## if eta < minimum -> use Eddington accretion -- M(t) \propto M0*exp(...)
#     eta = ma.masked_less(eta, eta_min)
#     i = ma.argmin(eta)
#     M1[eta.mask] = M0[eta.mask]/M0[i]*M1[i]
#     # the following exponential formula not continuous
#     # !!!!!!!! M1[eta.mask] = M0[eta.mask]*np.exp(mu_fit*f_duty*dt/(eta_min*t_Edd))    
#     return M1

def linear(xs,ys,x):
    i = np.argmax(xs>x)
    if i==0:
        i = 1
    t = (x-xs[i-1])/(xs[i]-xs[i-1])
    return ys[i-1]*(1-t)+ys[i]*t
    # return np.nansum([ys[i]*(1-t), ys[i+1]*t])

def M1M0_e(M0,dt,l):
    M1 = M0*np.exp(l*dt/(.1*t_Edd))
    # print('efolding:',l*dt/(.1*t_Edd))
    return M1

eta_0 = 0.1
M_cut = 10**logM0
def M1M0_d(M0, l, dt, delta, M_cut = M_cut):
    # dM_tilt/dt_tilt = l M_tilt/(eta_0*.5*(1+M_tilt**delta)); M_tilt=M/M_cut, t_tilt=t/t_Edd
    M0 = M0/M_cut
    M1 = M0
    # print('M0={0:e}\n'.format(M0))
    niter = 20
    for i in range(niter):
        fM1 = np.log(M1/M0) + (pow(M1, delta)-pow(M0, delta))/delta - l*dt/(.5*eta_0*t_Edd)
        dfdM1 = 1./M1 + pow(M1, delta-1)
        dM1 = - fM1/dfdM1
        M1 = M1 + dM1
        # print('i={0:d}, fM0={1:e}, dfdM={2:e}, M0={3:e}\n'.format(i,fM0,dfdM0,M0))

        if np.all(abs(dM1) < .01*M1):
            break
    if i==niter-1:
        print('not converge...')
    # print('final match?',np.log(M1/M0) + (pow(M1, delta)-pow(M0, delta))/delta - l*dt/eta_0/t_Edd )
    M1 *= M_cut
    return M1

# ---------------- kernel*: kernel of calculating P(lbd) integral ----------------
def kernelS_MBHmesh(M1, M0, dt, l_cut):
    # eta = 0.1; (1-eta)/eta
    xx,yy = np.meshgrid(M0, M1)
    lbd = np.log(yy/xx) / ( 9.*dt/t_Edd )
    return lbd/l_cut


# exponential growth
def kernel_MBH1(Mgrow_ratio, dt, f_duty, mu, sigma_dex):
    lbd = np.log(Mgrow_ratio)/( f_duty*dt/(0.1*t_Edd) )
    return np.log(lbd/mu) / (sigma_dex*np.log(10.)*math.sqrt(2.))

# power law growth + exp extrapolation
def kernel_MBH2(M1, M0, dt, f_duty, mu, sigma_dex, eta8, delta):
    lbd = ( pow(M1/1e8, delta) - pow(M0/1e8, delta) )/(f_duty*delta*dt)*(eta8*t_Edd)
    eta = eta8*pow(M0/1e8,delta)
    # if eta > maximum -> use Eddington accretion -- M(t) \propto M0*exp(...)
    eta = ma.masked_greater(eta, eta_max)
    # M1[eta.mask] = M0[eta.mask]*np.exp(mu*f_duty*dt/(eta_max*t_Edd))
    i = ma.argmax(eta)
    lbd[eta.mask] = lbd[i] * np.log(M1/M0[eta.mask]) / np.log(M1/M0[i])
    # lbd[eta.mask] = np.log(M1/M0[eta.mask])/ (f_duty*dt/(eta_max*t_Edd))
    eta[eta.mask] = eta_max

    # if eta < minimum -> use Eddington accretion -- M(t) \propto M0
    eta = ma.masked_less(eta, eta_min)
    i = ma.argmin(eta)
    # lbd[i] \propto log(M1/M0[i]) from Eddington accretion
    lbd[eta.mask] = lbd[i] * np.log(M1/M0[eta.mask]) / np.log(M1/M0[i])
    return np.log(lbd/mu) / (sigma_dex*np.log(10.)*math.sqrt(2.))

# piece-wise lbd; 2 exp + 1 pow
def kernel_MBH3(M1, M0, dt, f_duty, mu, sigma_dex, eta8, delta):
    lbd = ( pow(M1/1e8, delta) - pow(M0/1e8, delta) )/(f_duty*delta*dt)*(eta8*t_Edd)
    eta = eta8*pow(M0/1e8,delta)
    # if eta > maximum -> use Eddington accretion -- M(t) \propto M0*exp(...)
    eta = ma.masked_greater(eta, eta_max)
    # M1[eta.mask] = M0[eta.mask]*np.exp(mu*f_duty*dt/(eta_max*t_Edd))
    lbd[eta.mask] = np.log(M1/M0[eta.mask])/ (f_duty*dt/(eta_max*t_Edd))
    eta[eta.mask] = eta_max
    # if eta < minimum -> use Eddington accretion -- M(t) \propto M0
    eta = ma.masked_less(eta, eta_min)
    lbd[eta.mask] = np.log(M1/M0[eta.mask])/ (f_duty*dt/(eta_min*t_Edd))
    return np.log(lbd/mu) / (sigma_dex*np.log(10.)*math.sqrt(2.))

# lambda from Schechter function 
def kernelS_MBH(Mgrow_ratio, dt, f_duty, l_cut):
    lbd = np.log(Mgrow_ratio)/( f_duty*dt/(0.1*t_Edd) )
    return lbd/l_cut

# def kernelS_MBH_M(M1, M0, dt, f_duty, l_cut, d_fit, logM_0=logM0):
#     M_cut = pow(10., logM_0)
#     if d_fit:
#         lbd = .5*(np.log(M1/M0) + (pow(M1/M_cut,d_fit)-pow(M0/M_cut,d_fit))/d_fit) / ( f_duty*dt/(0.1*t_Edd) )
#     else:
#         lbd = np.log(M1/M0)  / ( f_duty*dt/(0.1*t_Edd) )
#     return lbd/l_cut

def kernelS_MBH_M(M1, M0, dt, f_duty, l_cut, d_fit, logM_0=logM0):
    M_cut = 1e8
    lbd = .5*(np.log(M1/M0) + (pow(M1/M_cut,d_fit)-pow(M0/M_cut,d_fit))/d_fit) / ( f_duty*dt/(0.1*t_Edd) )
    return lbd/l_cut

def kernelS_MBH_M_mesh(M1, M0, dt, f_duty, l_cut, d_fit, logM_0=logM0):
    M_cut = pow(10., logM_0)
    xx,yy = np.meshgrid(M0, M1)
    if d_fit:
        lbd = .5*(np.log(yy/xx) + (pow(yy/M_cut,d_fit)-pow(xx/M_cut,d_fit))/d_fit) / ( f_duty*dt/(0.1*t_Edd) )
    else:
        lbd = np.log(yy/xx)  / ( f_duty*dt/(0.1*t_Edd) )
    return lbd/l_cut


def kernelS_M1450(M1450, MBH, l_cut):
    lbd = Lbol_M1450(M1450)/(1.25e38*MBH)
    return lbd/l_cut

def kernelS_M1450_mesh(M1450, MBH, l_cut):
    xx,yy = np.meshgrid(MBH,M1450)
    zz = Lbol_M1450(yy)/(1.25e38*xx)
    return zz/l_cut

def kernel_M1450(M1450, MBH, mu, sigma_dex):
    lbd = Lbol_M1450(M1450)/(1.25e38*MBH)
    return np.log(lbd/mu) / (sigma_dex*np.log(10.)*math.sqrt(2.))

def ave_over(arr,n):
    # -1: value inferred from the length of the array and remaining dimensions.
    return np.mean(arr[:(len(arr)//n)*n].reshape(-1,n), axis=1)

# Willot+ 2010
def MF(M,z=6):
    alpha = -1.03
    Phi_star = 1.23e-8
    M_star = 2.24e9
    if z==6:
        return Phi_star*pow(M/M_star,alpha)*np.exp(-M/M_star)
    else:
        M_star *= pow(10,.5*(6-z))
        return Phi_star*pow(M/M_star,alpha)*np.exp(-M/M_star)

def L_M(M,Edd_ratio): # L_bol from M_BH in Msun
    return 1.25e38*Edd_ratio*M
def M_L(L,Edd_ratio): # M_BH in Msun from L_bol
    return L/(1.25e38*Edd_ratio)
def E_ML(M,L):
    return L/(1.25e38*M)

def Mdot2M(Mdot):
    eta = 1
    beta = 2.775e-6*(1.5)**.5
    Mdot1 = 0.04
    Mdot2 = 0.1
    if Mdot<Mdot1:
        M = eta*Mdot/beta
    elif Mdot>Mdot2:
        M = (0.83*np.log10(Mdot)+2.48)*1.e5
    else:
        M1 = eta*Mdot1/beta
        M2 = (0.83*np.log10(Mdot2)+2.48)*1.e5
        t = (np.log(Mdot)-np.log(Mdot1))/(np.log(Mdot2)-np.log(Mdot1))
        M = np.exp( t*np.log(M2) + (1-t)*np.log(M1) )
    return M

def LF(l): # dn/dlogL in Mpc^-3 dex^-1
    Phi_M_star = 1.14e-8
    M_star = -25.13
    alpha  = -1.5; beta = -2.81
    Phi_L_star = Phi_M_star * 2.5
    L_star = pow(10,-.4*(M_star-34.1)) * 3e18/1450 *1e7 * fbol_1450
    L_1 = pow(10,-.4*(-27.2-34.1)) * 3e18/1450 *1e7 * fbol_1450 
    L_2 = pow(10,-.4*(-20.7-34.1)) * 3e18/1450 *1e7 * fbol_1450 
    # print('break L',L_star/W37, 'Phi_L_star', Phi_L_star)
    t = (np.log10(l) - np.log10(L_1)) / (np.log10(L_2) - np.log10(L_1))
    return Phi_L_star/( pow(l/L_star,-(alpha+1)) + pow(l/L_star,-(beta+1)) ) * (2*(1-t)+3*t) 

def LF_Gal(x,z,mod='Schechter'): # dn/dmag in Gpc^-3 mag^-1
    if z==7:
        if mod == 'Schechter':     
            M=-20.49
            phi = 10**(-3.14)
            a=-1.88
            return np.log(10)/2.5*phi*(10**(0.4*(M-x)))**(a+1.0)*np.exp(-10**(0.4*(M-x)))*1e9
        if mod == 'DPL':
            Phi = 10**(-3.05)
            A=-1.89
            B=-3.81
            MM=-20.12
            return np.log(10)/2.5*Phi/(10**(0.4*(A+1)*(x-MM)) + 10**(0.4*(B+1)*(x-MM))) * 1e9
    if z==8:
        if mod == 'Schechter':
            M=-20.12
            phi = 10**(-3.35)
            a=-2.02
            return np.log(10)/2.5*phi*(10**(0.4*(M-x)))**(a+1.0)*np.exp(-10**(0.4*(M-x)))*1e9
        if mod == 'DPL':
            Phi = 3.3e-4
            A=-2.04
            B=-4.26
            MM=-20.02
            return np.log(10)/2.5*Phi/(10**(0.4*(A+1)*(x-MM)) + 10**(0.4*(B+1)*(x-MM))) * 1e9
    if z==9:
        if mod == 'Schechter':
            M = -21.29
            Phi = 10**(-4.88)
            a = -2.35
            return np.log(10)/2.5*Phi*(10**(0.4*(M-x)))**(a+1.0)*np.exp(-10**(0.4*(M-x)))*1e9
        if mod == 'DPL':
            Phi = 10**(-3.7)
            A=-2.1
            B=-3.34
            M=-19.62
            return np.log(10)/2.5*Phi/(10**(0.4*(A+1)*(x-M)) + 10**(0.4*(B+1)*(x-M))) * 1e9
    if z==10:
        if mod == 'Schechter':     
            M=-22.81
            phi = 10**(-6.56)
            a=-2.35
            return np.log(10)/2.5*phi*(10**(0.4*(M-x)))**(a+1.0)*np.exp(-10**(0.4*(M-x)))*1e9
        if mod == 'DPL':
            Phi = 10**(-4.44)
            A=-2.1
            B=-2.74
            MM=-19.60
            return np.log(10)/2.5*Phi/(10**(0.4*(A+1)*(x-MM)) + 10**(0.4*(B+1)*(x-MM))) * 1e9


def LF_M1450(M,z=6,W10=False): # dn/dmag in Mpc^-3 mag^-1
    if z==6: 
        # Matsuoka 2018
        Phi_M_star = 1.09e-8
        M_star = -24.9
        alpha  = -1.23; beta = -2.73
        # Willot 2010 CFHQS + SDSS
        if W10:
            Phi_M_star = 1.14e-8
            M_star = -25.13
            alpha  = -1.5; beta = -2.81
    elif z==5:
        # McGreer 2018 DPL; 
        Phi_M_star = pow(10., -8.97+0.47)
        M_star = -27.47
        alpha  = -1.97; beta = -4.
        # refit by Matsuoka 2018 (beta & M_star); me: (alpha & Phi_M_star)
        Phi_M_star = 3.8e-8
        M_star = -25.6
        alpha  = -1.23; beta = -3.
        # Niida 2020 Table 4
        Phi_M_star = 1.01e-7
        M_star = -25.05
        alpha  = -1.22; beta = -2.9
    elif z==4: # Akiyama 2018
        Phi_M_star = 2.66e-7
        M_star = -25.36
        alpha  = -1.3; beta = -3.11
    else:
        print("wrong redshift")
    return Phi_M_star/( pow(10., 0.4*(alpha+1)*(M-M_star)) + pow(10., 0.4*(beta+1)*(M-M_star)) ) #* (2*(1-t)+3*t) 

def M1450_Lbol(L):
    return 34.1-2.5*np.log10(L/(fbol_1450*3e18/1450*1e7))
def Lbol_M1450(M):
    return pow(10., -0.4*(M-34.1)) * (fbol_1450*3e18/1450*1e7)

# X-ray bolometric correction; Hopkins+07 & Duras+20
def K_AVE07(Lbol):
    return 10.83*pow(Lbol/(1e10*Lsun),0.28)+6.08*pow(Lbol/(1e10*Lsun),-0.02)
def K_AVE20(Lbol):
    a = 10.96
    b = 11.93
    c = 17.79
    return a*( 1 + pow(np.log10(Lbol/Lsun)/b,c) )

def KX_AVE20(Lx):
    a = 15.33
    b = 11.48
    c = 16.2
    return a*( 1 + pow(np.log10(Lx/Lsun)/b,c) )

# obscured fraction = Type II AGN fraction
def f_obsc_U14(logLx,z): # Ueda 14; 22< log NH < 24 fraction; as a func of Lx
    a1 = .48
    phi4375_0 = .43
    phi4375_z = phi4375_0*(1+z)**a1 if z<=2. else phi4375_0*(1+2.)**a1
    phimax = .84
    phimin = .2
    beta = .24
    phi = min( phimax, max(phi4375_z - beta*(logLx-43.75), phimin))
    f_obsc_sum = phi # sum over 22< log NH < 24 range
    print(0.72,phi4375_z)
    return f_obsc_sum

# constant obscured fraction; motivated by Vito+ 2018
f_obsc_const = .8

# correction factor including Compton thick AGNs; different fbol_Xray
def corr_U14H07(M1450): # Ueda+14 & Shankar+09
    L_bol = Lbol_M1450(M1450)
    f_bol = K_AVE07(L_bol)
    Lx = L_bol/f_bol
    eta = 1.7
    a1 = .48
    phi4375_0 = .43
    phi4375_z = phi4375_0*(1+2.)**a1
    phimax = .84
    phimin = .2
    beta = .24
    phi = min( phimax, max(phi4375_z - beta*(np.log10(Lx)-43.75), phimin))
    f_obsc_sum = phi # sum over 22< log NH < 24 range
    f_CTK = phi/2.
    return 1./(1-f_obsc_sum)
def corr_U14D20(M1450): # Ueda 14
    L_bol = Lbol_M1450(M1450)
    f_bol = K_AVE20(L_bol)
    Lx = L_bol/f_bol
    eta = 1.7
    a1 = .48
    phi4375_0 = .43
    phi4375_z = phi4375_0*(1+2.)**a1
    phimax = .84 #(1+eta)/(3+eta)
    phimin = .2
    beta = .24
    if isinstance(M1450,float):
        phi = min( phimax, max(phi4375_z - beta*(np.log10(Lx)-43.75), phimin))
    else:
        phi = np.zeros(len(M1450))
        for i in range(len(M1450)):
            phi[i] = min( phimax, max(phi4375_z - beta*(np.log10(Lx[i])-43.75), phimin))
    f_obsc_sum = phi # sum over 22< log NH < 24 range
    f_CTK = phi/2.
    return 1./(1-f_obsc_sum)

def corr_M14D20(M1450): # Merloni 14
    L_bol = Lbol_M1450(M1450)
    f_bol = K_AVE20(L_bol)
    Lx = L_bol/f_bol
    f_M = 0.56+1/np.pi*np.arctan((43.89-np.log10(Lx))/0.46)
    return 1./(1-f_M)

def LF_M1450_CO(M,z): # dn/dmag in Mpc^-3 mag^-1
    # Matsuoka 2018
    return LF_M1450(M,z)/(1-f_obsc_const)

def LF_M1450_DO(M,z): # dn/dmag in Mpc^-3 mag^-1
    # Matsuoka 2018
    return LF_M1450(M,z)*corr_U14D20(M)

def t_freefall(nH):
    C = np.sqrt( 32*G*(mu*m_H)/ (3*pi) )
    return 1./C/np.sqrt(nH)

def t_from_z(z): # age of universe at redshift z: tH = 2/(3Hz)
    return 2./(3*H0*np.sqrt(Omega_m0)) * pow(1+z, -1.5)

def z_tH(t): # t in Myr
    return pow( 2/(3*H0*np.sqrt(Omega_m0)*t*Myr), 2./3.)-1

def Tv(Mh,z):
    return alpha_T * (Mh/1.e8)**(2./3.) *  (1+z)/11.

def Mh_Tv(Tv,z):
    return 1.e8*(Tv/alpha_T/(1+z)*11.)**1.5

def Omega_mz(z):
    return Omega_m0*(1+z)**3 /(Omega_m0*(1+z)**3 + Omega_L0)

def Hz(z):
    return H0*np.sqrt( Omega_m0*(1+z)**3 + Omega_L0 ) 
 
def RHO_crit(z):
    return 3*pow(H0,2)/(8*pi*G)*(1+z)**3*Omega_m0/Omega_mz(z) 

class HALO:
    def __init__(self,M,z0):
        self.Mh = M
        self.z = z0
        self.c = 18*pow(self.Mh/(1.e11*Ms), -0.13)/(1+self.z) #concentration parameter c from Dekel & Birnboim 2006 Eq(22)
        c, z = self.c, self.z
        self.d = Omega_mz(z) - 1 
        d = self.d
        self.Delta_crit = 18.0*pi*pi + 82*d - 39*d*d  # Delta_crit ~ 200, overdensity
        Delta_crit = self.Delta_crit

        self.delta0 = self.Delta_crit/3.*pow(c,3)/(-c/(1+c) + np.log(1+c)) # characteristic overdensity parameter 
        delta0 = self.delta0
        self.rho_crit = RHO_crit(z)  # mean density of DM at z
        self.rho_c = self.rho_crit * delta0 

        self.Rvir = pow( self.Mh/(4./3*pi*Delta_crit*self.rho_crit),1./3. ) 
        self.Rs = self.Rvir/self.c 
        self.Vc = np.sqrt(G*self.Mh/self.Rvir) 

        self.t_dyn = self.Rvir/self.Vc 
        self.Tvir = G*self.Mh*(mu*m_H)/(2.*k_B*self.Rvir)
        self.gc = 2*c/(np.log(1+c) - c/(1+c)) 
        self.alpha = self.Tvir/self.Mh**(2./3)

    def Rho_r(self, r):
        rho_crit, delta0, Rvir = self.rho_crit, self.delta0, self.Rvir
        c, x = self.c, r/Rvir
        return rho_crit*delta0/( c*x * (1+c*x)**2 )

    # x = r/Rvir  c = Rvir/Rs
    def F_NFW(self,x):
        c = self.c
        return -c*x/(1+c*x) + np.log(1+c*x)
        
    def M_enc(self,r):
        rho_crit, delta0, Rs, Rvir = self.rho_crit, self.delta0, self.Rs, self.Rvir
        M_r = 4*pi*rho_crit*delta0*pow(Rs,3)*self.F_NFW(r/Rvir)
        return M_r

    def Phi(self, r):
        # lim r -> 0
        #return -4*pi*G*rho_crit*delta0*Rs*Rs 
        rho_crit, delta0, Rs = self.rho_crit,  self.delta0, self.Rs
        return -4*pi*G*rho_crit*delta0*(Rs**3)/r*np.log(1+r/Rs) 

# seeding files for Mh=1e11, 12, 13 Msun halos; v_bsm = 0, 1 sigma
Ts = [] # bsm=0,1 two files
Ntr = 10000
N_Mh = 3
eta = 0.3
for iM in range(N_Mh):
    TM = []
    for i_bsm in range(2):
        T=ascii.read(Mhpres[iM]+'Jcol_'+str(i_bsm)+'.txt', guess=False,delimiter=' ') #  None has np.where(T['z_col']==-1)
        T['Mdot'] = (k_B*T['Tg_loi']*T['f_loi']/(mu*m_H))**1.5/G/(Ms/yr)
        T['Mstar0'] = np.zeros(len(T))
        for i in range(len(T)):
            T['Mstar0'][i] = Mdot2M(T['Mdot'][i]*eta)
        T['t_col'] = t_from_z(T['z_col'])
        TM.append(T)
        # print(np.min(T['z_col']),np.max(T['z_col']),np.min(T['Mstar0']),np.max(T['Mstar0']))
    Ts.append(TM)

# MF bins
abin_mf =  np.logspace(2,12,num=800) # default endpoint=True
M_BH = abin_mf[:-1]*np.sqrt(abin_mf[1]/abin_mf[0])
bin_left = abin_mf[:-1]; bin_right = abin_mf[1:]
wid_mf = bin_right - bin_left
dlog10M = np.log10(abin_mf[1]/abin_mf[0]) # print('Mbin ratio',abin_mf[1]/abin_mf[0])
N_mf = len(abin_mf)-1

M0s = np.logspace(2,12,num=10000)
dlog10M0 = np.log10(M0s[1]/M0s[0])
eps = 1e-5

# LF bins
abin_lf = np.linspace(-32,-15,num=171)
dmag = abin_lf[1]-abin_lf[0]
L_left = abin_lf[:-1]; L_right = abin_lf[1:]
M1450  = (L_left+L_right)/2.
N_lf = len(M1450)

# bin_edg = abin_lf
# bin_wid = dmag
# bin_cen = M1450
# Phi_obs = LF_M1450(M1450)*1e9
# Phi_err = 10.**.5

# LF bins same w/ Matsu18
z = int(6)
bin_edg = bin_edg[str(z)]
bin_wid = bin_wid[str(z)]
bin_cen = bin_cen[str(z)]
Phi_obs = Phi_obs[str(z)]
Phi_err = Phi_err[str(z)]
