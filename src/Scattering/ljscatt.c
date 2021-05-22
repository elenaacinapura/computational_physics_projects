#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <print_routines.h>

#include "util.h"

int main(){

    /* system parameters */
    double m_D = 2.014; // uma
    double m_Ar = 39.948; // uma
    double s = 3.075; // armstrong
    double eps = 54.4; // kelvin

    /* conversion, do not edit */
    double kb = 8.617333262e-2; // meV/K    
    double mu = m_D * m_Ar / (m_D + m_Ar); // uma
    double uma = 931.49410242e9; // meV/c^2
    double s_ev = (s / 1.97327)*1e-6; // ħ*c/meV
    eps *= kb; // meV 

    /* numerical parameter xi */
    double xi = 1/(2*mu*uma*s_ev*s_ev*eps);
    printf("---------------------------------");
    printf("\n\n Value of the adimensional parameter xi = %lf\n\n",xi);
    printf("---------------------------------\n");

    /* energy in [E_start,E_end] meV */
    double E_start = 0.5 / eps;
    double E_end = 5.0 / eps;
    double dE = 0.01 / eps;
    
    /* other parameters */
    int l_max = 25;
    double L = 20.0;
    double dr = 0.001;

    int dim = (int)(L / dr);
    double u[dim];
    double r[dim];

    /* initial condition */
    r[0] = 0.5;
    r[1] = r[0] + dr;
    u[0] = exp(-sqrt(4.0 / (xi*25.0) ) * pow(r[0], -5));
    u[1] = exp(-sqrt(4.0 / (xi*25.0) ) * pow(r[1], -5));


    int idx1 = dim-8;
    int idx2 = dim-20;

    FILE* f_sigma;
    f_sigma = fopen("sigma_tot.csv","w");
    
    double E = E_start;
    while( E<E_end ) {
        double k = sqrt(E/xi);
        
        double sigma = 0.0;
        for(int l=0; l < l_max; l++){
            Param p = {xi,E,l};
            solve_numerov(r,u,dim,dr,F,&p);
            double delta = phase_shift(k,l,r[idx1],r[idx2],u[idx1],u[idx2]);
            sigma += (2*l+1) * pow(sin(delta),2);
        }
        sigma *= 4*M_PI/(k*k);

        fprint_double(f_sigma,E*eps);
        fprint_double(f_sigma,sigma*s*s);
        fprintf(f_sigma,"\n");

        E += dE;
    }

    fclose(f_sigma);

}