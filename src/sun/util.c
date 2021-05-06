#include <stdio.h>
#include <assert.h>
#include <math.h>

#include <print_routines.h>

#include "util.h"

double F_theta(double useless, double xi, void *p){
    assert(xi > 0.0);
    Param *par = (Param *)p;
    double eta = par->eta;
    return -eta/(xi*xi);
}

double F_eta(double useless, double xi, void *p){
    assert(xi > 0.0);
    Param *par = (Param *)p;
    double n = par->n;
    double theta = par->theta;
    return xi*xi*pow(theta,n);
}

void set_initial_condition(Vec *xi, Vec *theta, Vec *eta){
    xi->v[0] = START;
    theta->v[0] = 1.0;
    eta->v[0] = 0.0;
    xi->dim = 1;
    theta->dim = 1;
    eta->dim = 1;
}

void increment_dimension(Vec *xi, Vec *theta, Vec *eta){
    (xi->dim)++;
    (theta->dim)++;
    (eta->dim)++;
}

double integral(double xi_0, double dxi, Vec *theta, Param *p){
    int dim = theta->dim;
    dim -= 1;
    double n = p->n;
    
    double res = 0.0;
    for(int i=0;i<dim;i++){
        double xi = dxi*i;
        res += xi * xi * pow(theta->v[i],n) * dxi; 
    }
    return res;
}