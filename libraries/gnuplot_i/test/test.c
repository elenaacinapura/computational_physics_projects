#include "../gnuplot_i.h"
#include <string.h>


int main () {
    double x [] = {1.0, 2.0, 3.0, 4.0};
    double y [] = {2.0, 4.0, 6.0, 4.0};

    gnuplot_ctrl * h ;      /* un po' come FILE *f */
    h = gnuplot_init() ;    /* inizializza */

    gnuplot_set_xlabel(h, "x");
    gnuplot_set_ylabel(h, "y");

    // gnuplot_cmd(h, "set term pdfcairo");
    // gnuplot_cmd(h, "set output \'figure.pdf\'");
    gnuplot_set_term_png(h, "figura_nuova.png");

    gnuplot_plot_xy(h, x, y, 4, "Il mio plot!");

    gnuplot_close(h);
}