#include <assert.h>
#include <complex.h>
#include <differential_eq/numerov.h>
#include <math.h>
#include <numerical_methods/zero_bisection.h>
#include <print_routines.h>

#include "../routines.h"

int main() {
	/*=================== COMMON VARIABLES =================*/
	double dx = 1e-3;
	double a, xi;
	double L, x0, A, B;

	double E_start = -1.0;
	double dE = 0.0005;
	double E_end = 0.0;

	FILE *file;

	/*=================== Welcome ==================*/
	printf("=================================================\n");
	printf("BOUND STATES WITH NUMEROV ALGORITHM\n");
	printf("=================================================\n");

	/***********************************************
		Cosh potential
	***********************************************/
	/*============ Assign parameters ============*/
	xi = 0.05;
	a = 2.0 / xi;
	
	L = 5.0;
	x0 = 0.3;
	

	/***********************************************
		Lennard Jones Potential
	***********************************************/
	/*============ Assign parameters ============*/
	double kB = 1.38e-23;
	double eps = 35.7 * kB;
	double s = 2.79e-10;
	double hbar = 1.05e-34;
	double m = 10 * 1.66e-27;

	
}