#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gnuplot_i.h>
#include <integrate_notypesafe.h>
#include <print_routines.h>

#include "util.h"

int main(){

    /* valori di csi */
    double xi[] = {0.2,0.05,0.01};
    
    /* valori di alpha */
    int dim = (int)rint(alpha_max/dalpha);
    double alpha[dim];
    fill_alpha(alpha,dim);


    /*  */
    double cutoff = 50.0;
    int density = 500;
    Param_E par;
    double E;
    
    /* runno */
    FILE *f_E;
    f_E = fopen("E.csv","w");

    for(int i=0;i<dim;i++){
        
        fprint_double(f_E,alpha[i]);
        
        /* primo valore di xi */
        par.xi = xi[0];
        par.alpha = alpha[i];
        E = integrate_trap(integrand,-cutoff,cutoff,density,&par);
        fprint_double(f_E,E);
        
        /* secondo valore di xi */
        par.xi = xi[1];
        par.alpha = alpha[i];
        E = integrate_trap(integrand,-cutoff,cutoff,density,&par);
        fprint_double(f_E,E);
        
        /* primo valore di xi */
        par.xi = xi[2];
        par.alpha = alpha[i];
        E = integrate_trap(integrand,-cutoff,cutoff,density,&par);
        fprint_double(f_E,E);

        fprintf(f_E,"\n");

    }
    fclose(f_E);

    /* Plot */
    gnuplot_ctrl *h;
    h = gnuplot_init();
	gnuplot_cmd(h, "load \'plot_E.gp\'");
	system("eog E.png");



}