// The Lagrange polynomial interpolation
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

void m_lag(FILE *fp);

int main(void)
{
    FILE *fp;
    fp = fopen("my_lag.txt", "rw");
    m_lag(fp);

    return 0;
}

void m_lag(FILE *fp){
    // assignment the order of polynomial and independent variables
    int m, n;
    // initialize
    fscanf(fp, "%d%d", &m, &n);
    // assignemnt independent variables
    double *u;
    u = (double *)malloc(n*sizeof(double));
    for(int i = 0; i < n; i++){
        fscanf(fp, "%lf", &u[i]);
    }
    // assignment samples
    struct Sam{
        double *x;
        double *y;
    } samples;
    samples.x = (double *)malloc((m+1)*sizeof(double));
    samples.y = (double *)malloc((m+1)*sizeof(double));
    for(int i = 0; i <= m; i++){
        fscanf(fp, "%lf", &samples.x[i]);
    }
    for(int i = 0; i <= m; i++){
        fscanf(fp, "%lf", &samples.y[i]);
    }
    // calculate
    double ans;
    int b;
    // n independent variables
    for(int i = 0; i < n; i++,b=0,ans=0){
        // wether u == xi
        for(int j = 0; j <= m; j++){
            if(u[i]==samples.x[j]){
                b = j;
            }
        }
        if(b){
            printf("P(%lf) = %lf\n", samples.x[b], samples.y[b]);
            break;
        }
        // get ans
        double w = 1, dw = 1;
        for(int j = 0; j <= m; j++, w = 1, dw = 1){
            // assignment w, dw
            for(int k = 0; k <= m; k++){
                w *= u[i] - samples.x[k];
                if(k!=j){
                    dw *= samples.x[j] - samples.x[k];
                }
                // printf("%lf %lf\n", w, dw);
            }
            ans += samples.y[j]*w/((u[i]-samples.x[j])*dw);
        }
        printf("P(%lf) = %lf\n", u[i], ans);
    }
    free(u);
    free(samples.x);
    free(samples.y);
}