#ifndef __LJSCATT_H__
#define __LJSCATT_H__

/* TYPEDEF */
typedef struct Param{
    double xi;
    double E;
    int l;
}Param;

/* FUNCTIONS */
double lj(double r);
double F(double r, void *par);
void solve_numerov(double r[], double u[], int dim, double dr, double F(double, void*), void *p);
double phase_shift(double k, int l, double r1, double r2, double u1, double u2);


#endif