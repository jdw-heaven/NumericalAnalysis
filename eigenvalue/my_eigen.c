//  Jacobi's way to get eigenvalues and eigenvector
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define E 10e-12

void m_eigen(FILE *fp);
int m_jd(int n, double *A);
void m_max(int n, double *A, int *pp, int *pq);

int main(void)
{
    FILE *fp;
    fp = fopen("my_eigen.txt", "r");
    m_eigen(fp);

    return 0;
}
// judge whether to loop out the circulation (the sum of nondiagonal element's square less than the given error)
int m_jd(int n, double *A){
    int b = 1;
    double sum = 0;
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            if(i != j){
                sum += A[i*n+j]*A[i*n+j];
            }
        }
    }
    if(sum < E){
        b = 0;
    }
    return b;
}
// sezrch for the maximum value
void m_max(int n, double *A, int *pp, int *pq){
    *pp = 1;*pq = 0;
    double max = fabs(A[1*n+0]);
    for(int i = 1; i < n; i++){
        for(int j = 0; j < i; j++){
            if( fabs(A[i*n+j]) > max ){
                max = fabs(A[i*n+j]);
                *pp = i;*pq = j;
            }
        }
    }
}
// you can find the way in 张韵华 数值计算方法与算法 P177-P178,but he has some problem, I'll mark them below
void m_eigen(FILE *fp){
    int n;
    fscanf(fp, "%d", &n);
    double *A, *I;
    A = (double *)malloc(n*n*sizeof(double));
    I = (double *)malloc(n*n*sizeof(double));
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            fscanf(fp, "%lf", &A[i*n+j]);
            //printf("%4.3lf ", A[i*n+j]);
        }//printf("\n");
    }
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            fscanf(fp, "%lf", &I[i*n+j]);
            //printf("%4.3lf ", I[i*n+j]);
        }//printf("\n");
    }
    int p = 0, q = 0;
    int *pp, *pq;
    pp = &p; pq = &q;
    double s = 0, t = 0, c = 0, d = 0, t1 = 0, t2 = 0, x = 0, y = 0;
    do{
        m_max(n, A, pp, pq);
        //printf("p = %d, q = %d\n", p, q);
        s = (A[q*n+q]-A[p*n+p])/(2*A[p*n+q]);
        if( s==0 ){
            t = 1;
        }else{
            t1 = -s-sqrt(s*s+1); t2 = -s+sqrt(s*s+1);
            if( fabs(t1) > fabs(t2) ){
                t = t1;
            }else{
                t = t2;
            }
        }
        c = 1/sqrt(1+t*t); d = t/sqrt(1+t*t);
        for(int i = 0; i < n; i++){
            if( i!=p && i!=q ){
                x = c*A[p*n+i]-d*A[q*n+i]; y = c*A[q*n+i]+d*A[p*n+i];
                A[i*n+p] = x; A[p*n+i] = x;
                A[i*n+q] = y; A[q*n+i] = y;
            }
            A[p*n+p] = A[p*n+p]-t*A[p*n+q];
            A[q*n+q] = A[q*n+q]+t*A[p*n+q];
            A[p*n+q] = 0; // the book don't have this
            A[q*n+p] = 0; // the book don't have this
        }
        for(int i = 0; i < n; i++){
            x = I[i*n+p]; y = I[i*n+q];// ****tempory variable****
            I[i*n+p] = c*x-d*y;
            I[i*n+q] = d*x+c*y;

        }
    }while( m_jd(n, A) );
    printf("the eigenvalues are:\n");
    for(int i = 0; i < n; i++){
        printf("%.3lf  ", A[i*n+i]);
    }printf("\n");
    printf("the eigenvectors are:\n");
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            printf("%.3lf  ", I[i*n+j]);
        }printf("\n");
    }

    free(A);
    free(I);
}