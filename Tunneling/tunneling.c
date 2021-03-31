#include "util.h"
#include "../libraries/gnuplot_i/gnuplot_i.h"
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <assert.h>

int main () {
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
    outputfile = fopen("Output/phi.csv", "w");

    int idx = 1;
    while (idx < X.dim-1) {
        numerov_step(idx, V_square, &(Param_V_square) {V0, a});
        fprint_double(outputfile, X.v[idx]);
        fprint_double(outputfile, V_square(X.v[idx], &(Param_V_square) {V0, a}));
        fprint_double(outputfile, phi.v[idx]);
        fprintf(outputfile, "\n");
        idx++;
    }

    double T = calculate_T(X.dim-1, X.dim-5);

    printf("T = %lf\n", T);

    gnuplot_ctrl * h ;      /* un po' come FILE *f */
    h = gnuplot_init() ;

    gnuplot_set_term_png(h, "phi.png");

    gnuplot_cmd(h, "plot \'Output/phi.csv\' using 1:2 title \'V(x)\' w lines, \'Output/phi.csv\' using 1:3 title \'phi(x)\' w lines");

    fclose(outputfile);
    
}