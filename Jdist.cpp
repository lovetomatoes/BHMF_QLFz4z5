#include <stdio.h>
#include <string.h>
#include <cmath>
#include <iostream>
#include <iomanip>

#include <fstream>
#include <algorithm>
#include <sstream>
#include <cstdlib>
#include <time.h>

using namespace std;

//   g++ Jdist.cpp -o J.out && ./J.out 
//  "cj"

static int Nk;
double pi = 3.14159265358979, G=6.67e-8, c=2.9979245e10,
        yr=3.155e7, h_p=6.63e-27, Ms=2.e33, km = 1.e5, mu = 1.22,
        pc=3.0856775e18, kpc = 1.e3*pc, Mpc=1.e6*pc;

double xk_a[2000], Pk_a[2000];
static double h, Om, OL, Ob, sigma_8, Tcmb0, z_eq, rho_0, rho_m,
    z1, Om_z, OL_z, D_z, d_c, d, ratio_d, nu, nu_p, tff_6, M_ac,
    // rho_0=1.9e-29*(h*h),
    // rho_m=8.534251287809294e10*Ms*(h*h)/pow(Mpc,3),
    eps = 1.e-5, S = 1./3.;


/* ------------------------- Global Variables ------------------------ */
/* set_cosmo() -- sets labeled cosmological parameter models;
    the routine sets up all of the scalar quantites needed; 
	computation of the fitting formula.  
    The input parameters are: 
	1) Om -- Density of CDM, baryons, and massive neutrinos @ z=0,
				in units of the critical density. 
	2) Ob -- Density of baryons @ z=0, in units of critical. 
	3) OL -- Cosmological constant 
	4) h       -- Hubble constant, in units of 100 km/s/Mpc 
    5) d_c = 1.686 : critical density fluctuation 
    6) D_z growth factor (Carroll et al. 1992) delta(z) = delta(0) * D(z) / D(0), 
            where the normalization is set to D(0)=1 
        Iliev+03 defined it as  delta(z)=delta(0)/D(z)
        in 2nd paragraph in Sec 3.1.
    7) nu = d_c/D_z : d_c value extropolated to z=0
        In the PS formalism, nu is defined by nu = d_c / D_z / sigma_M,
        but Iliev et al. (2003) define nu = d_c/D_z
        Since the definition of nu is correct, 
        the calculation except g_fac is correct
        as long as D_z is defined as in Carroll+92.
*/

double W(double x){
    return 3.0/pow(x,3)*(sin(x)-x*cos(x));
}
double Gauss(double x, double y){
    return 1./sqrt(2.*pi*y) * exp(-x*x/2./y);
}
// x: delta_c, y: sigma(M)
double ff(double x, double y){
    return 1./sqrt(2.0*pi)*(x/pow(y,1.5))*exp(-x*x/y/2.0);
}
double Ez(double z){
    return sqrt(Om*pow(z1,3)+OL);
}
// g(z) in Carroll 1992; Galaxy formation & evolution Eq.4.76
double gg(double z){ // Eisenstein & Hu 1998 ApJ Eq(A4)
    double z1 = 1.+z, Om_z, OL_z; // must!!! since gg(z)/gg(0) uses 2 redshifts
    // double z_eq = 2.5e4*Om*hubble*hubble/pow(Tcmb0/2.7,4);
    Om_z = Om*pow(z1,3) / (Om*pow(z1,3)+OL);
    OL_z = 1.0 - Om_z;
    return 2.5 * Om_z / ( pow(Om_z,4./7.) - OL_z + (1.0 + Om_z/2.0)*(1.0 + OL_z/70.0));
}
void growth(double z,double& D_z){
    D_z = gg(z)/gg(0.0)/(1.0+z);
}

void linear(double* xa, double* ya, int m, double x, double& y){
   int ms;
   int i,j;
   double y1, y2, t;
   //if( x<xa[0] or x>xa[m-1]) printf("xa_min=%3.2e, x=%3.2e, xa_max=%3.2e\n", xa[0],x,xa[m-1]);
   for (i=0;i<m;i++){
      if (x-xa[i]<=0.) {
         ms = i;
         break;
      }
      else ms = m-1; //xa[m-1],xa[m-2]
   }
   if (ms==0) ms=1; //xa[1],xa[0]
   y1 = ya[ms-1];
   y2 = ya[ms];
   t=(x-xa[ms-1])/(xa[ms]-xa[ms-1]);
   y=(1.-t)*y1+t*y2;
}

void set_cosmo(double z, int i_cosmo=0){
    //****************** cosmo paras ******************
    // i_cosmo=0 : Planck 2015 EH98 Carroll 92
    // i_cosmo=1 : Planck 2020 CAMB Carroll 92
    if (i_cosmo==0) {
        h=0.677;
        Om=0.307;
        Ob = 0.0486;
        Tcmb0=2.725;
        z_eq = 2.5e4*Om*h*h/pow(Tcmb0/2.7,4);
        sigma_8 = 0.8159;
    }
    if (i_cosmo==1) {
        h=0.6732;
        Om=0.3158;
        Ob = 0.156*Om;
        sigma_8 = 0.8120;
    }
    OL=1.0-Om;
    rho_0 = 3*pow(100.*h*km/Mpc,2)/(8.*pi*G);
    rho_m = Om *rho_0;

    //**************** f(z) ***************
    z1 = 1. + z;
    Om_z = Om*pow(z1,3) / (Om*pow(z1,3)+OL);
    OL_z = 1.0 - Om_z;
    d_c = 1.686;
    d = Om_z - 1.;
    ratio_d = (Om/Om_z)*(18.*pi*pi +82.*d -39.*d*d)/(18.*pi*pi);
    tff_6 = 83./pow( z1/11., 1.5); // free-fall timescale in Myr
    M_ac = pow( 1.e4/1.98e4 / (mu/.6) * pow(ratio_d,-S) / (z1/10) ,1.5) * 1.e8/h*Ms;

    //**************** fluctuation growth ***************
    growth(z, D_z);
    nu = d_c/D_z;
    nu_p = nu * sqrt(0.707);
    //**************** store Power Spectrum array ***************
    string line, k_str, Pk_str, name;
    // i_cosmo==0
    name = "../data/hmf_Pk_Planck15_EH98.dat";
    Nk = 567;
    if (i_cosmo==1) name = "../data/hmf_Pk_Planck20_CAMB.dat";
    ifstream inFile(name); if (!inFile) cout<<"read error\n";
    int j = 0;
    getline(inFile,line); //skipping 1st row in case it is column name line
    while (getline(inFile, line)){
        istringstream isolated(line);
        isolated>>k_str>>Pk_str;
        xk_a[j] = stod(k_str)/Mpc; Pk_a[j] = stod(Pk_str);
        j++;
    }
    inFile.close();

    double krat, dk, dk_vol, Power, xk, sigma_M, sigma_M_0, xx, radius;
    radius = 8*Mpc/h;
    xk           =    1e-8/Mpc;
    krat         =    1.001;
    sigma_M      =    0.;
    sigma_M_0    =    1e-99;
    while (abs(sigma_M/sigma_M_0-1.)>1.0e-12){
        dk = xk*(krat-1.);
        xk += dk;
        dk_vol = xk * xk * dk / (2.*pi*pi);
        xx = xk * radius;
        linear(xk_a,Pk_a,Nk,xk,Power);
        Power=Power/pow(xk,3)*2.*pi*pi;
    // mass variance 1
        sigma_M_0 = sigma_M;
        sigma_M   += Power * W(xx) * W(xx) * dk_vol;
    }
    // // normalized to sigma_8; no need already consistent
    // for (j=0;j<Nk;j++) Pk_a[j] *= sigma_8/sigma_M;
    // // consistent with cosmo para sigma8= 0.8159
    printf("sigma8 for icosmo=0 is :%10.4f\n", sqrt(sigma_M));

}

// return Sigma_M = sigma(M)^2
void SIGMA_M(double Mass, double& sigma_M){
    double radius, xk, krat, dk, dk_vol, xx, Power, sigma_M_0;

    radius = pow(3.*Mass/4./pi/rho_m, S);
    xk           = 0.001/Mpc;
    krat         = 1.001;
    sigma_M_0    = 1.0e-99;
    sigma_M      = 0.;

    while (abs(sigma_M/sigma_M_0-1.)>1.0e-12){
        dk = xk*(krat-1.);
        xk += dk;
        dk_vol = xk * xk * dk / (2.*pi*pi);
        xx = xk * radius;
        linear(xk_a,Pk_a,Nk,xk,Power);
        Power=Power/pow(xk,3)*2.*pi*pi;
    // mass variance 1
        sigma_M_0 = sigma_M;
        sigma_M   += Power * W(xx) * W(xx) * dk_vol;
    }
// sigma_M^2 -> sigma_M
    sigma_M  =  sqrt(sigma_M);
}

void dN_ST_dM(double Mass, double& n_ST){
    double xk, krat, dk, dk_vol, Power;
    double radius, radius_p, xx, xx_p, sigma_M_0, sigma_M, sigma_Mp_0, sigma_Mp;

    radius    = pow(3.*Mass/4./pi/rho_m, S);
    radius_p  = pow(3.*Mass*(1.+eps)/4./pi/rho_m, S);

    xk            = 0.001/Mpc;
    krat          = 1.001;
    sigma_M_0     = 1.0e-99;
    sigma_M       = 0.;
    sigma_Mp_0    = 1.0e-99;
    sigma_Mp      = 0.;

    while (abs(sigma_M/sigma_M_0-1.)>1.0e-12 or abs(sigma_Mp/sigma_Mp_0-1.)>1.0e-12){
        dk = xk*(krat-1.);
        xk += dk;
        dk_vol = xk * xk * dk / (2.*pi*pi);
        xx = xk * radius;
        xx_p = xk * radius_p;
        linear(xk_a,Pk_a,Nk,xk,Power);
        Power=Power/pow(xk,3)*2.*pi*pi;
    // mass variance M & Mp
        sigma_M_0   =    sigma_M;
        sigma_M     +=   Power * W(xx) * W(xx) * dk_vol;
        sigma_Mp_0  =    sigma_Mp;
        sigma_Mp    +=   Power * W(xx_p) * W(xx_p) * dk_vol;
    }
// sigma_M^2 -> sigma_M
    sigma_M = sqrt(sigma_M);
    sigma_Mp = sqrt(sigma_Mp);
    n_ST = -sqrt(2./pi) * (rho_m/Mass)
                    * nu_p / pow(sigma_M, 2)
                    * (sigma_Mp-sigma_M)/(Mass*eps)              // dSigma/dM
                    * exp(-pow(nu_p/sigma_M,2)/2.)
                    * 0.3222 * ( 1.0 + pow(nu_p/sigma_M, -0.6) );
                    // / pow(1./Mpc,3);                                    // dndm in [1/Mpc]^3 g^-1
                    // * Mass_2 / pow(h/Mpc,3);                         // dndlnm in [h/Mpc]^3
                    // * log(10.) * Mass_2 / pow(1./Mpc,3);                         // dndlog10m in [1/Mpc]^3
}

// fraction of universe @ z contained in mass > Mass
double F_Mz(double Mass, double z){
    double sigma, d_c=1.69, D_z, nu, F;
    set_cosmo(z);
    SIGMA_M(Mass, sigma);
    growth(z, D_z);
    nu = d_c/sigma/D_z;
    printf("nu = %5.4e\n", nu);
    F = 0.4* (1.+0.4/pow(nu,0.4))*erfc(0.85*nu/sqrt(2.));
    return F;
}

int main(){
    clock_t t0 = clock();
    fstream file;
    file<<setiosflags(ios::scientific)<<setprecision(5);
    string infile, outfile;
// ----------------------------------- Sigma_M ---------------------------------------
    // double z = 0, Mass = 2e9*Ms, dN_dlogM_ST, n_ST;
    // set_cosmo(z);
    // growth(z,D_z);
    // double Mratio = pow(10.,.01);
    // Mass = 1e8*Ms;
    // double s1;
    // file.open("./sigmaM.txt", ios::out | ios::trunc);
    // file<<setw(15)<<"Mass"<<setw(15)<<"sigma"<<setw(15)<<"nu"<<endl;
    // while (Mass<1.e14*Ms){
    //     SIGMA_M(Mass,s1);
    //     nu = d_c/s1/D_z;
    //     // dN_dlogM_ST = n_ST * Mass / pow(1./Mpc,3);             // dndlnm in [1/Mpc]^3
    //     file<<setw(15)<<Mass/Ms<<setw(15)<<s1<<setw(15)<<nu<<endl;
    //     Mass *= Mratio;
    // }
    // file.close();
// ----------------------------------- D_z ---------------------------------------
    // z = 0.;
    // file.open("./growth_factor.txt", ios::out | ios::trunc);
    // file<<setw(15)<<"Mass"<<setw(15)<<"sigma"<<setw(15)<<"nu"<<endl;
    // while (z<30){
    //     set_cosmo(z);
    //     growth(z,D_z);
    //     file<<setw(15)<<z<<setw(15)<<D_z<<endl;
    //     z += .1;
    // }
    // file.close();

// ----------------------------------- DNDM_ST ---------------------------------------
    //********** use file hmf_Pk.dat to compare w/ hmf ***************//
    double n_ST;
    double sigma_M_1, sigma_M_2, d_sigma_M_2, fsigma;
    double z = 17;
    set_cosmo(z,1);
    double nu_p = nu * sqrt(0.707);
    // 只看Mass_2
    double Mass_2 = 1e6*Ms; 
    int Nrat = 1000;
    double qrat = pow(1e7, 1./double(Nrat));
    double dlog10M = log10(qrat), n_tot=0.;
    file.open("./dnSTdM_z"+to_string(int(z))+"cpp.txt", ios::out | ios::trunc);
    file<<setw(15)<<"M"<<setw(15)<<"n_ST"<<endl;
    for (int i=0; i<Nrat; i++){
        dN_ST_dM(Mass_2,n_ST);
        n_ST = n_ST  * log(10.) *Mass_2 / pow(1./Mpc,3);         // dndlogm in [1/Mpc]^3
        n_tot += n_ST*dlog10M;
        file<<setw(15)<<Mass_2/Ms<<setw(15)<<n_ST<<endl;
        Mass_2 *= qrat;
    }
    file.close();
    printf("M2=%.1e, n_tot=%.3e, dlog10M=%.1e\n",Mass_2/Ms,n_tot,dlog10M);
    clock_t t1 = clock();
    printf("dt=%.2e\n", (double)(t1-t0)/CLOCKS_PER_SEC);
}