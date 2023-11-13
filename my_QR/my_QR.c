/* QR fractorization_HouseHolder-algorithm I'll give A's QR fractorization matrix --- QT and R, then test that QT is orthometric and 
R is Superior-triangular*/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

void m_QR(FILE *fp);

void main(void)
{
    FILE *fp;
    fp = fopen("my_QR.txt", "r");

    m_QR(fp);
}

void m_QR(FILE *fp){
    int n, np;
    //the order of matrix
    fscanf(fp, "%d", &n);
    //define matrix A and QT which initialize as an unit matrix
    double *A, *QT0, *QT1, *I, *H;
    A = (double *)malloc(n*n*sizeof(double));
    QT0 = (double *)malloc(n*n*sizeof(double));
    QT1 = (double *)malloc(n*n*sizeof(double));
    I = (double *)malloc(n*n*sizeof(double));
    H = (double *)malloc(n*n*sizeof(double));
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            fscanf(fp, "%lf", &A[i*n+j]);
            if(i==j)QT0[i*n+j] = 1;else QT0[i*n+j] = 0;
            if(i==j)QT0[i*n+j] = 1;else QT0[i*n+j] = 0;
            if(i==j)I[i*n+j] = 1;else I[i*n+j] = 0;
            if(i==j)H[i*n+j] = 1;else H[i*n+j] = 0;
        }
    }
    //check
    /*printf("A is:\n");
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            printf("%lf ", A[i*n+j]);
        }printf("\n");
    }printf("\n");
    printf("QT is:\n");
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            printf("%lf ", QT[i*n+j]);
        }printf("\n");
    }printf("\n");*/
    //calculate vector v and parameter b, then give HA
    double *v, *w, *wp;
    double y, t, b, a, v0;
    for(int i = 0; i < n; i++){
        v = (double *)malloc((n-i)*sizeof(double));
        w = (double *)malloc((n-i)*sizeof(double));
        wp = (double *)malloc((n-i)*sizeof(double));
        //find the max value and normalize
        y = fabs(A[i*n+i]);
        for(int j = i; j < n; j++){
            if( y < fabs(A[j*n+i]) )y = fabs(A[j*n+i]);
        }
        //check
        //printf("%lf\n", y);
        t = 0;
        for(int j = i+1; j < n; j++){
            t += A[j*n+i]*A[j*n+i];
            v[j-i] = A[j*n+i];
        }
        //check
        //printf("%lf\n", t);
        //for(int j = 0; j < n-i; j++)printf("%lf ", v[j]);printf("\n");
        if( t==0 )b = 0;
        else{
            a = sqrt( A[i*n+i]*A[i*n+i]+t );
            if( A[i*n+i] <= 0 )v[0] = A[i*n+i]-a; else v[0] = -t/(A[i*n+i]+a);
            b = 2*v[0]*v[0]/(t+v[0]*v[0]);
            v0 = v[0];
            for(int j = 0; j < n-i; j++)v[j] = v[j]/v0;
        }
        //check
        //printf("%lf\n", b);
        //for(int j = 0; j < n-i; j++)printf("%lf ", v[j]);printf("\n");
        //calculate w = b*AT*v
        for(int j = i; j < n; j++){
            w[j-i] = 0;
            for(int k = i; k < n; k++){
                    w[j-i] += A[k*n+j]*v[k-i];
            }
            w[j-i] = b*w[j-i];
        }
        for(int j = i; j < n; j++){
            wp[j-i] = 0;
            for(int k = i; k < n; k++){
                    wp[j-i] += I[k*n+j]*v[k-i];
            }
            wp[j-i] = b*wp[j-i];
        }
        for(int j = i; j < n; j++){
            for(int k = i; k < n; k++){
                A[j*n+k] = A[j*n+k] - v[j-i]*w[k-i];
                H[j*n+k] = I[j*n+k] - v[j-i]*wp[k-i];                
            }
        }
        for(int j = 0; j < i; j++){
            for(int k = 0; k < n; k++){
                if(j==k)H[j*n+k] = 1;else H[j*n+k] = 0;
                if(j==k)H[k*n+j] = 1;else H[k*n+j] = 0;
            }
        }
        //check
        /*for(int j = 0; j < n; j++){
            for(int k = 0; k < n;k++){
                printf("%lf ", H[j*n+k]);
            }printf("\n");
        }printf("\n");*/

        for(int j = 0; j < n; j++){
            for(int k = 0; k < n; k++){
                QT1[j*n+k] = 0;
                for(int l = 0; l < n; l++){
                    QT1[j*n+k] += H[j*n+l]*QT0[l*n+k];
                }
            }
        }

        for(int j = 0; j < n; j++){
            for(int k = 0; k < n; k++){
                QT0[j*n+k] = QT1[j*n+k];
            }
        }


        free(v);
        free(w);
        free(wp);
    }


    printf("A is:\n");
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            printf("%lf ", A[i*n+j]);
        }printf("\n");
    }printf("\n");

    printf("QT is:\n");
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            printf("%lf ", QT0[i*n+j]);
        }printf("\n");
    }printf("\n");
    printf("QTQ is(theorically an unit matrix):\n");
    double x;
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            x = 0;
            for(int k = 0; k < n; k++){
                x += QT0[j*n+k]*QT0[i*n+k];
            }
            printf("%lf ", x);
        }printf("\n");
    }printf("\n");

    
    free(A);
    free(QT0);
    free(QT1);
    free(I);
    free(H);
}