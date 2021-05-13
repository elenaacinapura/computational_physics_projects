#include <stdio.h>
#include <assert.h>
#include <math.h>

#include <differential_eq/runge_kutta.h>
#include <print_routines.h>

#include "util.h"

int main(){

    /* physical parameters */
    double G = 6.67e-11; // N * m^2 * kg^-2
    double rho = 1.622e5; // kg * m^-3
    double R = 7e8; // m
    double M = 2e30; // kg

    /* parameters, now I use reduced units */
    double dxi = 1e-3;

    Vec xi,theta,eta;
    Param p;

    double n = 3.0;
    p.n = n;

    /* initial conditions */
    set_initial_condition(&xi,&theta,&eta);  

    /* solve equation */
    int i = 1;
    do{
        p.eta = eta.v[i-1];
        p.theta = theta.v[i-1];
        theta.v[i] = RK4(theta.v[i-1],xi.v[i-1],dxi,F_theta,&p);
        eta.v[i] = RK4(eta.v[i-1],xi.v[i-1],dxi,F_eta,&p);
        xi.v[i] = xi.v[i-1] + dxi;
        increment_dimension(&xi,&theta,&eta);
        i++;
    }while( theta.v[i-1]>EPS );

    double xi_0 = xi.v[i-1];
    printf("xi_0 = %lf\n",xi_0);

    /* calculating the integral */
    double rho_star = rho * R*R*R / M;
    printf("rho = %lf\n",rho_star);

    double integ = integral(xi_0,dxi,&theta,&p);
    printf("integral = %lf\n",integ);
    
    double result = pow(xi_0,3) / (4.0 * M_PI * integ);
    printf("confronta = %lf\n",result);



}