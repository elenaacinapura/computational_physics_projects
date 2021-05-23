#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <assert.h>
#include <gsl/gsl_fft_complex.h>
#include <print_routines.h>
#include <linear_algebra/blas_wrappers.h>

/* --------------------------- */
#define EPS 1e-6
#define DIM_MAX (int)1e6
#define TYPE 0 // 1 for animation
#define POTENTIALTYPE 0 // 0 for Tunneling(1) and 1 for Tunneling(2)

/* --------------------------- */
void read_ground(double x[], double psi[], int N);
double potential(double x);
void evolution_step(complex double psi[], complex double rho[], complex double eta[], int N);
void run_for_animation(double T, double dt, double x[], double V[], complex double Psi[], 
                complex double rho[], complex double eta[], int N);
double integrate_for_p(complex double Psi[], double x[], int N);
void fourier_analysis(FILE *file, complex double f[], double L, int N, double *m1, double *m2);
double calculate_norm(double x[], complex double Psi[], int N);
void points_of_maximum(double x[], double f[], int N, double *m1, double *m2);

/************************** MAIN *********************************/
int main(){

    /* parameters */
    double xi = 0.015;
    double L = 2.0 * 4.096;

    /* ground state energy */
    int N = 2 * 4096;
    double x[N], psi[N], V[N];
    read_ground(x,psi,N);

    /* split operator method */
    double T = 500.0;
    double dt = 1.0;
    complex double K, rho[N], eta[N];
    
    for(int i=0;i<N;i++){

        V[i] = potential(x[i]);
        rho[i] = cexp(-I * V[i] * dt);
        K = I * 2 * M_PI / L * ( i <= N/2 ? i : i - N );
        eta[i] = cexp( I * xi * K * K * dt );

    }

    /* complexification of psi */
    complex double Psi[N];
    for(int i=0;i<N;i++){
        Psi[i] = psi[i];
    }

    /******************** EVOLUTION *******************/
    if(!TYPE){
        
        T = pow(2,12);
        dt = 1.0;

        int tau = (int)(T/dt);
        assert(tau % 2 == 0);
        complex double p_time[tau];
        double prob;
        
        FILE *f_p;
        f_p = fopen("probability.csv","w");
        assert(f_p != NULL);
        
        double t = 0.0;
        int cnt = 0;
        do{

            prob = integrate_for_p(Psi,x,N);
            p_time[cnt] = prob;

            fprint_double(f_p,t);
            fprint_double_newline(f_p,prob);
            
            evolution_step(Psi,rho,eta,N);
            t += dt;
            cnt++;

        }while(t < T);
        fclose(f_p);
    
        /* spectral analysis */ 
        double m1,m2;
        FILE *f_spectrum;
        f_spectrum = fopen("spectrum.csv","w");
        fourier_analysis(f_spectrum,p_time,T,tau,&m1,&m2);
        fclose(f_spectrum);

        /* print period */
        printf("\nT = %lf\n",2*M_PI/m1);
    
    }

    /********************* ANIMATION *************************/
    if(TYPE){
        run_for_animation(T,dt,x,V,Psi,rho,eta,N);
    }
    

}
/*****************************************************************/
void read_ground(double x[], double psi[], int N){
    
    /* read from file the ground state wave function */
    assert(N == 2 * 4096);
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
    
    if(x > 0 && !POTENTIALTYPE){
       
        return  x * x * (x - 1.0);
    
    }else if(x > 0 && POTENTIALTYPE){
    
        return  2.0 * x * x * (x - 1.0);
    
    }
    return -x * x * (x + 1.0);
}

void evolution_step(complex double psi[], complex double rho[], complex double eta[], int N){

    /* potential phase shift */
    for(int i=0;i<N;i++){
        psi[i] *= rho[i];
    }

    /* kinetic energy phase shift */
    gsl_fft_complex_radix2_forward((double *)psi,1,N);

    for(int i=0;i<N;i++){
        psi[i] *= eta[i];
    }

    gsl_fft_complex_radix2_inverse((double *)psi,1,N);
}

void run_for_animation(double T, double dt, double x[], double V[], complex double Psi[], 
                complex double rho[], complex double eta[], int N){

    /* evolution, for the plot */
    const int stp = (int)ceil(T/dt);

    int density = (int)stp / 100;

    double data[N][stp/density+2];
    for(int n=0;n<N;n++){
        data[n][0] = x[n];
        data[n][1] = V[n];
    }

    int cnt = 0, col = 0;
    do{
        if(cnt % density == 0){
            for(int n=0;n<N;n++){
                data[n][col+2] = pow(cabs(Psi[n]),2);
            }
            col++;
        }
        evolution_step(Psi,rho,eta,N);
        cnt++;
    }while(cnt<stp);
    assert(cnt == stp);
    assert(col == stp/density);

    /* print */
    FILE *file;
    file = fopen("animation.csv","w");

    for(int i=0;i<N;i++){
        for(int j=0;j<stp/density+2;j++){
            fprint_double(file,data[i][j]);
        }
        fprintf(file,"\n");
    }
    fclose(file);

}

double integrate_for_p(complex double Psi[], double x[], int N){
    
    int cnt = 0;
    double dx = x[1] - x[0];
    double res = 0.0;
    while(x[cnt] < 0.0){
        res += 0.0;
        cnt++;
    }
    while(cnt < N){
        res += pow(cabs(Psi[cnt]),2) * dx;
        cnt++;
    }

    return res;
}

void fourier_analysis(FILE *file, complex double f[], double L, int N, double *m1, double *m2){

    gsl_fft_complex_radix2_forward((double *)f,1,N);

    double k_ord[N];
    for(int n=0;n<N;n++){
        k_ord[n] = 2 * M_PI / L * n - M_PI * N / L;
    }

    double f_abs_ord[N];
    for(int n=0;n<=N/2-2;n++){
        f_abs_ord[n] = cabs(f[N/2+1+n]);
    }
    for(int n=N/2-2+1;n<N;n++){
        f_abs_ord[n] = cabs(f[n-N/2+1]);
    }

    for(int n=0;n<N;n++){
        fprint_double(file,k_ord[n]);
        fprint_double_newline(file,f_abs_ord[n]);
    }

    points_of_maximum(k_ord,f_abs_ord,N,m1,m2);

}

double calculate_norm(double x[], complex double Psi[], int N){
    double dx = x[1] - x[0];
    double norm = 0.0;
    for(int i=0;i<N;i++){
        norm += pow(cabs(Psi[i]),2) * dx;
    }
    return norm;
}

void points_of_maximum(double x[], double f[], int N, double *m1, double *m2){

    int cnt = 0;
    while(x[cnt] < EPS){
        cnt++;
    }

    /* first point of maximum */
    double point = 0.0;
    int pos = 0;
    double max1 = 0.0;
    for(int i=cnt;i<N;i++){
        if(f[i] >= max1){
            max1 = f[i];
            point = x[i];
            pos = i;
        }
    }
    *m1 = point;

    /* second point of maximum */
    *m2 = 0.0;

}