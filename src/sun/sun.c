#include <assert.h>
#include <differential_eq/runge_kutta.h>
#include <math.h>
#include <print_routines.h>
#include <stdio.h>

#include "util.h"

int main() {
	/*========================= WELCOME ========================*/
	printf("=================================================\n");
	printf("A MODEL FOR THE SUN\n");
	printf("=================================================\n");
	printf("Calculating...\n\n");
	/*========================= PARAMETERS ========================*/
	double G = 6.67e-11;   // N * m^2 * kg^-2
	double rho = 1.622e5;  // kg * m^-3
	double R = 7e8;		   // m
	double M = 2e30;	   // kg

	double rho_star = rho * R * R * R / M;

	double dxi = 1e-3;

	Vec xi, theta, eta;
	Param p;

	/*========================= CYCLE ON n ========================*/
	double n_min = 3.0;
	double dn = 0.0001;
	double n = n_min;
	double rho_estimated, xi_0;
	do {
		n += dn;
		p.n = n;
		/* initial conditions */
		set_initial_condition(&xi, &theta, &eta);
		/* solve equation */
		int i = 1;
		do {
			p.eta = eta.v[i - 1];
			p.theta = theta.v[i - 1];
			theta.v[i] = RK4(theta.v[i - 1], xi.v[i - 1], dxi, F_theta, &p);
			eta.v[i] = RK4(eta.v[i - 1], xi.v[i - 1], dxi, F_eta, &p);
			xi.v[i] = xi.v[i - 1] + dxi;
			increment_dimension(&xi, &theta, &eta);
			i++;
		} while (theta.v[i - 1] > EPS);
		/* set xi_0 */
		xi_0 = xi.v[i - 1];
		/* calculating the integral */
		double integ = integral(xi_0, dxi, &xi, &theta, &p);
		rho_estimated = pow(xi_0, 3) / (4.0 * M_PI * integ);
	} while (fabs(rho_estimated - rho_star) > 0.01);
	double alpha = R / xi_0;
	double K = alpha * alpha * 4.0 * M_PI * G / ((n + 1.0) * pow(rho, 1.0 / n - 1.0));
	double P = K * pow(rho, (n + 1.0) / n);
	printf("Calculations ended successfully!\n");
	printf("=================================================\n");
	printf("Result of the calculations: \n\t n = %.3lf\n\t P = %lf billion bar\n", n, P / 1e14);
	printf("=================================================\n");
}
