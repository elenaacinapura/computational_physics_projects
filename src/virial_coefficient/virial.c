#include <math.h>
#include <assert.h>
#include <print_routines.h>
#include <numerical_methods/integrate_notypesafe.h>

typedef struct Param {
    double T;
} Param;

double integrand (double r, void *param) {
    Param *p = (Param *) param;
    double T = p->T;

    double V = 4.0 * (pow(r, -12) - pow (r, -6));

    return -2 * M_PI * r*r * (exp(-V/T) - 1.0);
}

int main () {
    Param p;
    double T_start = 1.0;
    double T_end = 50.0;
    double dT = 0.1;

    FILE *file;
    file = fopen("virial.csv", "w");
    fprintf(file, "T\tB\n");

    double T = T_start;
    while (T <= T_end + 0.1* T_end) {
        p.T = T;
        double B = integrate(integrand, 0.0, 50.0, 100, &p);
        fprint_double(file, T);
        fprint_double_newline(file, B);
        T += dT;
    }
    fclose(file);

}