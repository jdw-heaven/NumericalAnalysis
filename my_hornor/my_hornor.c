// Hornor algorithm to solve polynomial
#include<stdio.h>
#include<math.h>
#include<stdlib.h>

int main(void)
{
    int n, m;
    double *a;
    double *x;
    double *ans;
    FILE *fp;
    fp = fopen("my_hornor.txt", "r");
    //get the order of polynomial and x
    fscanf(fp, "%d", &n);
    fscanf(fp, "%d", &m);
    //initialize the co and x
    a = (double *)malloc((n+1)*sizeof(double));
    x = (double *)malloc(m*sizeof(double));
    ans = (double *)malloc(m*sizeof(double));
    for(int i = 0; i < n+1; i++){
        fscanf(fp, "%lf", &a[i]);
    }
    for(int i = 0; i < m; i++){
        fscanf(fp, "%lf", &x[i]);
    }
    for(int i = 0; i < m; i++){
        ans[i] = x[i]*a[n];
        for(int j = n-1; j > 0; j--){
            ans[i] = x[i]*(a[j]+ans[i]);
        }
        ans[i] += a[0];
        printf("when x = %lf, ans = %lf\n", x[i], ans[i]);
    }

    free(a);
    free(x);
    free(ans);
    return 0;
}