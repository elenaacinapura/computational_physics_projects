#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gnuplot_i.h>
#include <integrate_notypesafe.h>
#include <print_routines.h>

#include "util.h"

int main(){
    
    /* vettore di alpha */
    int dim = (int)rint(alpha_max/dalpha);
    double alpha[dim];
    fill_alpha(alpha,dim);
    
    /* runno */
    double xi[] = {0.1,0.025,0.005};
    Param_E par;
    double E;
    FILE *f_E;
    f_E = fopen("E.csv","w");

    for(int i=0;i<dim;i++){
        
        fprint_double(f_E,alpha[i]);
        
        par.alpha = alpha[i];
        
        /* primo valore di xi */
        par.xi = xi[0];
        E = calculate_E(par);
        fprint_double(f_E,E);
        
        /* secondo valore di xi */
        par.xi = xi[1];
        E = calculate_E(par);
        fprint_double(f_E,E);
        
        /* primo valore di xi */
        par.xi = xi[2];
        E = calculate_E(par);
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