#include <math.h>
#include <assert.h>
#include <numerical_methods/integrate_notypesafe.h>
#include <numerical_methods/zero_bisection.h>
#include <stdio.h>

/*=============== PARAMETER STRUCT ============*/
typedef struct Param {
    double T;
} Param;
/*=============== FUNCTION HEADERS ============*/
double integrand (double r, void *param);
double B (double T, void *param);
/*==================== MAIN ===================*/
int main () {
    Param useless_struct;
    double boyle_T = zero_bisection(B, 0.5, 5.0, &useless_struct);
    printf("================================================\n");
    printf("BOYLE TEMPERATURE\n");
    printf("================================================\n");
    printf("I found the following Boyle Temperature (in reduced units):\nT = %.3lf\n", boyle_T);
    printf("================================================\n");

}
/*=============== FUNCTION BODIES ============*/
double integrand (double r, void *param){
    Param *p = (Param *) param;
    double T = p->T;

    double V = 4.0 * (pow(r, -12) - pow (r, -6));

    return -2 * M_PI * r*r * (exp(-V/T) - 1.0);
}
double B (double T, void *p) {
    Param par;
    par.T = T;
    return integrate(integrand, 0.0, 50.0, 100, &par);
}