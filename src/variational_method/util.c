#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gnuplot_i.h>
#include <integrate_notypesafe.h>
#include <print_routines.h>

#include "util.h"

double cutoff = 50.0;
int density = 500.0;
double alpha_max = 0.5;
double dalpha = 0.01;

double integrand(double x,void *p){
    Param_E *par = p;
    double xi = par->xi;
    double alpha = par->alpha;
    
    double res = -sqrt(2/ (M_PI*alpha*alpha) ) * exp(-2*x*x/ (alpha*alpha) );
    res *= xi*( (4*x*x/pow(alpha,4)) - (2/ (alpha*alpha) ) ) + pow(cosh(x),-4);
    return res;

}

double integrate_trap (double f (double, void *), double low, double high, int density, void *param) {
    int n = (int)ceil((high - low) * density);  //number of points
    double dx = (high - low) / (n + 1);         // increment in integration variable
    double k = 0.0;

	//trapezoidal rule 
    for(int i=1;i<n;i++){
        k += f(low+i*dx,param);
    }
    k = (dx/2)*(f(low,param)+f(high,param)+2*k);
    return k;

/*	//modified
	for (int i = 1; i < n; i++) { 
		double x1 = low + dx*(i);
		double x2 = x1 + dx;             
		I += f((x1 + x2)/2.0, param);
	}
	return I *= dx;
*/

}

void fill_alpha(double alpha[],int dim){
    alpha[0] = 5*dalpha;
    for(int i=1;i<dim;i++){
        alpha[i] = alpha[i-1] + dalpha;
    }
}









