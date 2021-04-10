#include "util.h"

#include <gnuplot_i.h>
#include <integrate_notypesafe.h>
#include <print_routines.h>
#include <zeros_newton.h>

#include <assert.h>
#include <math.h>
#include <setjmp.h>
#include <stdbool.h>
#include <stdio.h>

int main () {
    Params p;
    double xi[] = {0.05, 0.01, 0.005};

    bool stop = 0;
    for (int i = 0; i < 3; i++) {
        p.xi = xi[i];
        printf("---------------------------------------------\nxi = %.03lf\n", xi[i]);
        for (int n = 0; n < 10; n++) {
            if (stop) {
                stop = 0;
                break;
            }
            p.n = n;
            if (setjmp(JumpBuffer)) {
                printf("End of bound states reached for n = %d\n\n", n);
                stop = 1;
                continue;
            }
            double E = zero_newton(f, -0.0001, &p);
            printf("E = %lf\n", E);
        }
    }

}