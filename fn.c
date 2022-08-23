#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#define  EULER 2.718281
int n;
//int** adj;

double exp(double ex){
    return pow(EULER, ex);
}

void Get_args(
      char*    argv[]        /* in  */,
      int*     n,
      int*     times      /* out */,
      int*     nthreads      /* out */
      
      )
    {
        *n = strtol(argv[1], NULL, 10);
        *times = strtol(argv[2], NULL, 10);
        *nthreads = strtol(argv[3], NULL, 10);
    }

int choice(){
    //srand(time(NULL));
    int a=-1;
    if(floor(2*(double)rand()/ (double)RAND_MAX)>=1)
    a=1;
    return a;
}

void vector_make(int*** mat,int n)
{
     // make vector 
    (*mat)=malloc(n*sizeof(int*));
    for(int i=0;i<n;i++)
    {  
        (*mat)[i]=malloc(sizeof(int)*n);
        for(int j=0;j<n;j++)
        {   
            //printf("making\n");
            (*mat)[i][j]=choice();
        }
    }
}
int neighbor(int** mat, int ki, int kj, int n){
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
int adj(int** mat,int ki,int kj,int n){
    int tar;
    int** adj;
    adj=malloc(sizeof(int*)*n);
    for (int i = 0; i < n; i++)
    {
        adj[i]=malloc(sizeof(int)*n);
        for (int j = 0; j < n; j++)
        {
            adj[i][j]=neighbor(mat,i,j,n);
        }
    }
    tar=adj[ki][kj];

    for (int i = 0; i < n; i++)
    {
        free(adj[i]);
    }
    free(adj);
    return tar;
}
void markov_ising_anneal(int*** mat,int n, int times, int thread){
    int ki,kj;
    int E=0;
    double B=0.01;
    int h;
    int del_E;
    double gamma;
    int i;
  
   # pragma omp parallel for num_threads(thread)\
    private(i) schedule(static,2) 
   for (i = 0; i < times; i++) {
    //srand(time(NULL));
    ki=floor(n*(double)rand()/ (double)RAND_MAX);
    kj=floor(n*(double)rand()/ (double)RAND_MAX);
    h=adj((*mat),ki, kj, n);
    del_E=2*h*(*mat)[ki][kj];
    gamma=exp(-B*del_E);
    srandom(i*time(NULL));
        //# pragma omp atomic
    if(((double)rand()/ (double)RAND_MAX)<gamma){
            
                (*mat)[ki][kj]=-1*(*mat)[ki][kj];
                //# pragma omp atomic
                E+=del_E;
    }
        //# pragma omp atomic
    B=1.00005*B;
    }
   
}


void Usage(char prog_name[] /* in */) {
   fprintf(stderr, "usage: %s ", prog_name); 
   fprintf(stderr, "<n_side> <times> <omp_count>\n");
   exit(0);
}  /* Usage */



void print_mat(int** mat, int n,char* filename){
    FILE* pgmimg;
    pgmimg = fopen(filename, "wb");
    // Writing Magic Number to the File
    fprintf(pgmimg, "P2\n"); 
    // Writing Width and Height
    fprintf(pgmimg, "%d %d\n", n, n); 
  
    int temp;
    // Writing the maximum gray value
    fprintf(pgmimg, "255\n"); 
    // int count = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            temp = mat[i][j];
            // Writing the gray values in the 2D array to the file
            fprintf(pgmimg, "%d ", temp);
        }
        fprintf(pgmimg, "\n");
    }
    fclose(pgmimg);  
}


int main(int argc, char* argv[])
{
    // size of the matrix
    int n;
    // number of threads
    int thread_count;
    int times;
    int** mat;
    struct timespec begin,end;
    double elapsed;
    srand(time(0));
    if (argc != 4) Usage(argv[0]); 
    Get_args(argv, &n,&times,&thread_count);
    // thread_count=omp_get_max_threads();
    //printf("1\n");
    vector_make(&mat,n);
    printf("matrix before simulation\n");
    //print_mat(&mat,n);
    print_mat(mat,n,"before_omp.pgm");
    clock_gettime(CLOCK_MONOTONIC,&begin);
    markov_ising_anneal(&mat,n,times,thread_count);
    clock_gettime(CLOCK_MONOTONIC,&end);
    printf("matrix after simulation\n");
    //print_mat(&mat,n);
    elapsed=(end.tv_nsec-begin.tv_nsec)/1e9+(end.tv_sec-begin.tv_sec);
    print_mat(mat,n,"after_omp.pgm");
    printf("Elapsed time: %f\n",elapsed);
    for (int i=0; i<n;i++){
        free(mat[n]);
    }
    free(mat);
    return 0;
}