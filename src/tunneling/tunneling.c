#include "util.h"
#include <gnuplot_i.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <assert.h>
#include <stdlib.h>

#include <sqrt_babylon.h>

int main () {
    printf("babylon sqrt (2) = %lf\n", sqrt_babylon(2.0));
    /* INITIAL CONDITIONS */
    L = 10.0;
    E = 0.5;
    dx = 1e-3;

    /* initialize X */
    fill_X();

    /* initialize phi */
    phi.dim = 2;
    double k = 2.0*E;
    phi.v[0] = cexp(-I*k*L);
    phi.v[1] = cexp(-I*k*(L+dx));

    double V0 = 1.0;
    double a = 1.0;

    FILE * outputfile;
    outputfile = fopen("tunneling_phi.csv", "w");

    int idx = 1;
    while (idx < X.dim-1) {
        numerov_step(idx, V_square, &(Param_V_square) {V0, a});
        fprint_double(outputfile, X.v[idx]);
        fprint_double(outputfile, V_square(X.v[idx], &(Param_V_square) {V0, a}));
        fprint_double(outputfile, creal(phi.v[idx]));
        fprintf(outputfile, "\n");
        idx++;
    }

    double T = calculate_T(X.dim-1, X.dim-5);

    printf("|T|^2 = %lf\n", square_cabs(T));

    gnuplot_ctrl * h ;      /* un po' come FILE *f */
    h = gnuplot_init() ;

    gnuplot_set_term_png(h, "phi.png");
    gnuplot_cmd(h, "set xlabel \'x\'");

    gnuplot_cmd(h, "plot \'tunneling_phi.csv\' using 1:2 title \'V(x)\' w lines, \'tunneling_phi.csv\' using 1:3 title \'Re(phi(x))\' w lines");

    system("eog phi.png"); /* open png from terminal */

    fclose(outputfile);
    
}