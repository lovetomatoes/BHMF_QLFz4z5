#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// grid resolution and iteration number
#define NM 600
#define Niter 40000

// set constant prams.
double const_pc = 3e18; //[cm]
double const_msun = 2e33; //[g]
double const_G = 6.67e-8; //[g^-1 cm^3 s^-2]
double const_pi = 3.1415926;
double const_rsun = 7e10; //[cm]
double const_Stefan = 5.670373e-5; //[erg/s/cm^2/K^4]
double const_c = 3e10; //[cm/s]
double const_kappa = 0.4; //[cm^2/g]
double const_yr = 3.15e7; //[s]

void cal_phi(double *phi, double *mass){

    double dm;
    double tot = 0.0;
    int i;

    for (i=0; i < NM-1; i++){
        dm = mass[i+1]-mass[i];
        tot += (phi[i+1]+phi[i])*dm*0.5;
    }

    for (i=0; i < NM; i++){
        phi[i] /= tot;
        if (phi[i] < 1e-10) phi[i] = 1e-10;
    }

}

double cal_mdot(double mass, double lmdot01, double alpha, double lmb){
    double mdot, mdot01, mb;

    mdot01  = pow(10.0,lmdot01);
    mb      = pow(10.0,lmb);

    mdot    = mdot01*pow(mass/0.1,2.0)*pow(1.0+mass/mb,alpha);

    return mdot;
}

double nrand(){
        static int sw=0;
        static double r1,r2,s;

        if (sw==0){
                sw=1;
                do {
                        r1=2.0*drand48()-1.0;
                        r2=2.0*drand48()-1.0;
                        s=r1*r1+r2*r2;
                } while (s>1.0 || s==0.0);
                        s=sqrt(-2.0*log(s)/s);
                        return r1*s;
        }
        else {
                sw=0;
                return r2*s;
        }
}

int main(int argc, char *argv[]){ //argv[1] = flg, argv[2] = mmin_chem, argv[3] = mmax_chem, argv[4] = gnm_chem

    int  i, j, k;
    int flg;
    double gnm = -2.35, mmin = 0.05, mmax = 1000.0, mmin_i = 0.1, mmax_i = 50.0; //[Msun]
    double gnm_chem, mmin_chem, mmax_chem; //[Msun]
    double dm;
    double phi[NM] = {0.0}, phi_sal[NM] = {0.0}, phi_chem[NM] = {0.0}, mass[NM] = {0.0};
    double time, dt, dt0;
    double tmax = 5.0e6; //[yr]
    int iter_max = 100000000;
    double mdot0, mdotm; //[Msun/yr]
    
    // Free parameters in MCMC cal.
    double lmdot01, alpha, lmb, lhd; 
    double lmdot01_p, alpha_p, lmb_p, lhd_p;
    double lmdot01_max, alpha_max, lmb_max, lhd_max;
    double s_lmdot01 = 1e-2, s_alpha = 1e-2, s_lmb = 1e-2;
    double x,y,z,xx;

    if (argc < 2) flg = 0;
    else flg = atoi(argv[1]);

    if (argc < 3) mmin_chem = 1.0;
    else mmin_chem = atof(argv[2]);

    if (argc < 4) mmax_chem = 130.0;
    else mmax_chem = atof(argv[3]);

    if (argc < 5) gnm_chem = -0.5;
    else gnm_chem = atof(argv[4]);

    printf("IMF No. %d, mmin = %3.1e, mmax = %3.1e, gnm = %3.1e\n",flg+1,mmin_chem,mmax_chem,gnm_chem);

    FILE *outputfile1, *outputfile2;
    if (flg == 0){
        outputfile1 = fopen("param_mcmc_1.txt", "w");
        outputfile2 = fopen("imf_mcmc_1.txt", "w");
    }else if (flg == 1){
        outputfile1 = fopen("param_mcmc_2.txt", "w");
        outputfile2 = fopen("imf_mcmc_2.txt", "w");
    }else if (flg == 2){
        outputfile1 = fopen("param_mcmc_3.txt", "w");
        outputfile2 = fopen("imf_mcmc_3.txt", "w");
    }

    // intialize mass and IMF
    for (i = 0; i < NM; i++) {
        mass[i] = mmin*pow(10.0,log10(mmax/mmin)*i/(NM-1));
        if (mass[i] < mmin_i){
            phi_sal[i] = 0.0;
        }else if (mass[i] >= mmin_i && mass[i] <= mmax_i){
            phi_sal[i] = pow(mass[i],gnm);
        }else{
            phi_sal[i] = 0.0;
        }
        if (mass[i] < mmin_chem){
            phi_chem[i] = 0.0;
        }else if (mass[i] >= mmin_chem && mass[i] <= mmax_chem){
            phi_chem[i] = pow(mass[i],gnm_chem);
        }else{
            phi_chem[i] = 0.0;
        }
    }
    cal_phi(phi_sal,mass);
    cal_phi(phi_chem,mass);

    // main calculation
    k = 0;
    lhd_p = -1e100;
    lhd_max = -1e100;
    lmdot01_p = -7.7;
    alpha_p = -2.0;
    lmb_p = 1.0;

    // start MCMC procedure
    for (k = 0; k < Niter; k++){

        x = nrand();
        y = nrand();
        z = nrand();
        lmdot01  = x*s_lmdot01+lmdot01_p;
        alpha    = y*s_alpha+alpha_p;
        lmb      = z*s_lmb+lmb_p;

        for (i = 0; i < NM; i++) phi[i] = phi_sal[i];

        j = 0;
        time = 0.0;
        while (time < tmax && j < iter_max){

            // set time step
            dt = 1e20;
            for (i = 1; i < NM; i++) {
                dm = mass[i]-mass[i-1];
                mdot0 = cal_mdot(mass[i],lmdot01,alpha,lmb);
                dt0 = 0.5*dm/fabs(mdot0);
                if (dt0 < dt){
                    dt = dt0;
                }
            }
            time += dt;
            j++;

            // calculate IMF
            for (i = 1; i < NM; i++) {
                dm = mass[i]-mass[i-1];
                mdot0 = cal_mdot(mass[i],lmdot01,alpha,lmb);
                mdotm = cal_mdot(mass[i-1],lmdot01,alpha,lmb);
                phi[i] -= dt/dm*mdot0*(phi[i]-phi[i-1]*mdotm/mdot0); //Upwind difference scheme
            }

            // normalize IMF
            cal_phi(phi,mass);

        }

        // evaluate the likelihood value
        lhd = 0.0;
        for (i = 0; i < NM; i++) {
            if (mass[i] > mmin_chem && mass[i] < 0.7*mmax_chem){
                lhd += -(log10(phi[i])-log10(phi_chem[i]))*(log10(phi[i])-log10(phi_chem[i]))/(2.0*5.0*5.0);
            }
            else if (mass[i] > 1.3*mmax_chem){
                lhd += -(log10(phi[i])-log10(phi_chem[i]))*(log10(phi[i])-log10(phi_chem[i]))/(2.0*5.0*5.0);
            }
        }

        xx = drand48();
        if (lhd > lhd_p){
            lmdot01_p = lmdot01;
            alpha_p = alpha;
            lmb_p = lmb;
            lhd_p = lhd;
        } else {
            if (xx < exp(lhd-lhd_p)){
                lmdot01_p = lmdot01;
                alpha_p = alpha;
                lmb_p = lmb;
                lhd_p = lhd;
            }
        }

        if (k%10 == 0) printf("%d, %e, %e, %e, %e\n",k, lhd_p, pow(10.0,lmdot01_p), alpha_p, pow(10.0,lmb_p));

        if (k > Niter/2){
            fprintf(outputfile1,"%e %e %e %e\n",lhd_p, pow(10.0,lmdot01_p), alpha_p, pow(10.0,lmb_p));
            for (i = 0; i < NM; i++) {
                fprintf(outputfile2," %e ", phi[i]);
            }
            fprintf(outputfile2,"\n");
        }

    }

    return 0;

}
