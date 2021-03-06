#ifndef __UTIL_H__
#define __UTIL_H__

#include <stdio.h>
#define N 108   // numero particelle
#define S 200
// #define EPS 0.0001

extern double x[N][3], v[N][3], a[N][3];
extern double dx[N][N][3];
extern double g[S];

extern double L,rho, cutoff, T, dt, t, duration;

void print_mat (double m[N][3]);

void print_double (double d, FILE *f);

double min_double (double a, double b);

double r_polari (double x, double y, double z);

double lj_u (double r);

double lj_part (double r);

void calculate_acc ();

void calculate_distance();

void calculate_forces();

void verlet_step ();

double calculate_U ();

double calculate_K ();

double calculate_T ();

double calculate_P();

void rescale_velocities ();

#endif