#include <assert.h>
#include <math.h>
#include <print_routines.h>

/*========================== PARAMETERS ==========================*/
const double R0 = 2.0;
const double S0 = 1.0 - 1e-5;
const double I0 = 1e-5;
const double r0 = 0.0;
const double dt = 1e-3;
const double T = 40;

/*======================= FUNCTION HEADERS =======================*/
double f_S(double S, double I, double R);
double f_I(double S, double I, double R);
double f_R(double S, double I, double R);
void RK4_step(double *t, double *S, double *I, double *R);

/*============================= MAIN =============================*/
int main() {
    double S, I, R;
    S = S0;
    I = I0;
    R = r0;

    FILE *file;
    file = fopen("sir.csv", "w");
    fprintf(file, "t\tS\tI\tR\n");

    double t = 0.0;
    while(t < T + 0.1*dt) {
        RK4_step(&t, &S, &I, &R);

        /* Print to file */
        fprint_double(file, t);
        fprint_double(file, S);
        fprint_double(file, I);
        fprint_double_newline(file, R);
    }
    fclose(file);
}

/*========================== FUNCTIONS ===========================*/
double f_S(double S, double I, double R) {
	return -R0 * I * S;
}
double f_I(double S, double I, double R) {
	return R0 * I * S - I;
}
double f_R(double S, double I, double R) {
	return I;
}
void RK4_step(double *t, double *S, double *I, double *R) {
    double s = *S;
    double i = *I;
    double r = *R;
	double S_new, I_new, R_new;
	double K1, K2, K3, K4;
	/*======= S ======*/
	K1 = f_S(s, i, r);
	K2 = f_S(s + 0.5 * dt * K1, i, r);
	K3 = f_S(s + 0.5 * dt * K2, i, r);
	K4 = f_S(s + dt * K3, i, r);
	S_new = s + dt / 6.0 * (K1 + 2.0 * K2 + 2.0 * K3 + K4);
	/*======= I ======*/
	K1 = f_I(s, i, r);
	K2 = f_I(s, i + 0.5 * dt * K1, r);
	K3 = f_I(s, i + 0.5 * dt * K2, r);
	K4 = f_I(s, i + dt * K3, r);
	I_new = i + dt / 6.0 * (K1 + 2.0 * K2 + 2.0 * K3 + K4);
    /*======= R ======*/
    K1 = f_R(s, i, r);
	K2 = f_R(s, i, r + 0.5 * dt * K1);
	K3 = f_R(s, i, r + 0.5 * dt * K2);
	K4 = f_R(s, i, r + dt * K3);
	R_new = r + dt / 6.0 * (K1 + 2.0 * K2 + 2.0 * K3 + K4);

    *S = S_new;
    *I = I_new;
    *R = R_new;
    *t += dt;
}