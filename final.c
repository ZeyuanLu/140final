#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#define  EULER 2.718281
int n;
double exp(double ex){
    return pow(EULER, ex);
}
void Get_args(
      char*    argv[]        /* in  */,
      int*     n) 
      {

   n = strtol(argv[1], NULL, 10);
}

int choice(){
    srand(time(NULL));
    return floor(2*(double)rand()/ (double)RAND_MAX);
}
void print_mat(double*** mat, int n){

    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            printf("%d ",mat[i][j]);
        }
        printf("\n");
    }
}

void vector_make(double*** mat,int n)
{
     // make vector 
    (*mat)=malloc(n*sizeof(double*));
    for(int i=0;i<n;i++)
    {  
        (*mat)[i]=malloc(sizeof(double)*n);
        for(int j=0;j<k;j++)
        {
            (*mat)[i][j]=choice();
        }
    }
}
int neighbor(double*** mat, int ki, int kj, int n){
    int adj=0;
    if(ki>0)
        adj+=mat[ki-1][kj];
    if(ki<n-1)
        adj+=mat[ki+1][kj];
    if(kj>0)
        adj+=mat[ki][kj-1];
    if(kj<n-1)
        adj+=mat[ki][kj+1];
    return adj;
}
void markov_ising_anneal(double*** mat,int n, int thread_count){
    int ki,kj;
    int E=0;
    double B=0.01;
    int h;
    int del_E;
    double gamma;
  # pragma omp parallel for num_threads(thread_count)\
   private(i) schedule(static,2) 
   for (i = 0; i < thread_count; i++) {
    srand(time(NULL));
    ki=floor(n*(double)rand()/ (double)RAND_MAX);
    kj=floor(n*(double)rand()/ (double)RAND_MAX);
    h=neighbor(mat,ki, kj, n);
    del_E=2*h*mat[ki][kj];
    gamma=exp(-B*del_E);
    srandom(time(NULL));
    if(((double)rand()/ (double)RAND_MAX)<gamma){
        # pragma omp critical
            mat[ki][kj]=-1*mat[ki][kj];
            E+=del_E;
    }
    B=1.00005*B;
   
   }
}