#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <assert.h>
#include <gsl/gsl_fft_complex.h>
#include <print_routines.h>

/* --------------------------- */
#define EPS 1e-6

/* --------------------------- */
// typedef struct Param_pot{
//     double v,a;
// }Param_pot;

/* --------------------------- */
//const double T = 250.0;

/* --------------------------- */
void read_ground(double x[], double psi[], int N);
double potential(double x);
// complex double potential(const double x, Param_pot *p);
// void evolution_step(complex double psi[], complex double rho[], complex double eta[]);
// void calculate_norms(double x[], complex double psi[], double *n1, double *n2);

/************************** MAIN *********************************/
int main(){

    /* parameters */
    double xi = 0.015;
    double L = 2.0 * 4.096;

    /* ground state energy */
    double E_0 = -0.018028;
    int N = 2 * 4096;
    double dx, x[N], psi[N];
    read_ground(x,psi,N);
    dx = x[1] - x[0];

    /* potential */
    double V[N];
    for(int i=0;i<N;i++){
        V[i] = potential(x[i]);
    }

    /* split operator method */
    double T = 1000.0;
    double dt = 0.001;
    complex double K, rho[N], eta[N];
    
    for(int i=0;i<N;i++){

        V[i] = potential(x[i]);
        rho[i] = cexp(-I * V[i] * dt);
        K = I * 2 * M_PI / L * ( i <= N/2 ? i : i - N );
        eta[i] = cexp( I * xi * K * K * dt );

    }

    
    /* test */
    FILE *fp;
    fp = fopen("ground.csv","w");
    for(int i=0;i<N;i++){
        fprint_double(fp,x[i]);
        fprint_double(fp,V[i]);
        fprint_double_newline(fp,psi[i]);
    }
    fclose(fp);


}
/*****************************************************************/
void read_ground(double x[], double psi[], int N){
    
    /* read from file the ground state wave function */
    assert(N = 2 * 4096);
    int N_prov = N / 2;
    double x_prov[N_prov], psi_prov[N_prov];
    FILE *f_psi;
    f_psi = fopen("eigenfunction.csv","r");
    assert(f_psi != NULL);
    for(int i=0;i<N_prov;i++){
        fscanf(f_psi,"%lf%lf",&x_prov[i],&psi_prov[i]);
    }
    fclose(f_psi);

    /* initialize vectors */
    double dx = x_prov[1] - x_prov[0];
    for(int i=0;i<N;i++){
        if(i<N_prov){
            x[i] = x_prov[i];
            psi[i] = psi_prov[i];
        }else{
            x[i] = x[i-1] + dx;
            psi[i] = 0.0;
        }
    }

    /* check normalization */
    double norm = 0.0;
    for(int i=0;i<N;i++){
        norm += psi[i] * psi[i] * dx;
    }
    assert(fabs(norm - 1.0) < EPS);

}

double potential(double x){
    if(x < 0){
        return - x * x * (x + 1);
    }
    return x * x * (x - 1);
}

// complex double potential(const double x, Param_pot *p){
//     return 0.0;    
// }

// void evolution_step(complex double psi[], complex double rho[], complex double eta[]){

//     /* potential phase shift */
//     for(int i=0;i<N;i++){
//         psi[i] *= rho[i];
//     }

//     /* kinetic energy phase shift */
//     gsl_fft_complex_radix2_forward((double *)psi,1,N);

//     for(int i=0;i<N;i++){
//         psi[i] *= eta[i];
//     }

//     gsl_fft_complex_radix2_inverse((double *)psi,1,N);
// }

// void calculate_norms(double x[], complex double psi[], double *n1, double *n2){
//     int cnt = 0;

//     *n1 = 0.0;
//     while(x[cnt] < 0.0){
//         *n1 += pow(cabs(psi[cnt]),2);
//         cnt++;
//     }
//     *n2 = 0.0;
//     while(cnt < N){
//         *n2 += pow(cabs(psi[cnt]),2);
//}