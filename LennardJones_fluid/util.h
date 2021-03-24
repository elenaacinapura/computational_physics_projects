#ifndef __UTIL_H__
#define __UTIL_H__

#define N 108   // numero particelle

double x [N][3];
double v [N][3];
double a [N][3];

extern double L, rho, cutoff, T, dt, t, duration;

void print_mat (double m[N][3]);

void print_double (double d, char filepath []);

double r_polari (double x, double y, double z);

double lj_u (double r);

double lj_part (double r);

void calculate_acc ();

void verlet_step ();

double calculate_U ();

double calculate_K ();

double calculate_T ();

void rescale_velocities ();

#endif