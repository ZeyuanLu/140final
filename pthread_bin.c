#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#define  EULER 2.718281

int** mat;
int n,load;
int times;
double E=0;
double B=0.01;
pthread_mutex_t* mutexes;
struct timespec begin, end;
double elapsed;

double exp(double ex){
    return pow(EULER, ex);
}

int min(int a, int b) {
    return a < b ? a : b;
}

void Get_args(
      char*    argv[]        /* in  */,
      int*     n,
     int*    times, 
      int*     nthreads
)
    {
        *n = strtol(argv[1], NULL, 10);
        *times = strtol(argv[2], NULL, 10);
        *nthreads = strtol(argv[3], NULL, 10);
    }

int choice(){
    // srand(time(NULL));
    // printf("%f\n",floor(2*(double)rand()/ (double)RAND_MAX));
    return floor(2*(double)rand()/ (double)RAND_MAX);
}

void vector_make(int*** mat,int n)
{
     // make vector 
    (*mat)=malloc(n*sizeof(double*));
    for(int i=0;i<n;i++)
    {  
        (*mat)[i]=malloc(sizeof(double)*n);
        for(int j=0;j<n;j++)
        {
            (*mat)[i][j]=choice();
            if((*mat)[i][j]==0)
            {
                (*mat)[i][j]=-1;
            }
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

void* markov_ising_anneal(void* rank){
    int ki,kj;
    int h;
    double del_E;
    double gamma;
    long thread_num=(long)rank;
    long i;
    // printf("check\n");
    // printf("%d\n",load);
   for (i = thread_num*load; i < min((thread_num+1)*load,times); i++) {
        ki=floor(n*(double)rand()/ (double)RAND_MAX);
        kj=floor(n*(double)rand()/ (double)RAND_MAX);
        // printf("%d %d\n",ki,kj);
        pthread_mutex_lock(&mutexes[ki]);
        h=adj(mat,ki, kj, n);
        del_E=2*h*mat[ki][kj];
        gamma=exp(-B*del_E);
        // printf("%d %d %d %d %f\n",ki,kj,h,del_E,gamma);
        srandom(thread_num*time(NULL));
        // printf("thread:%ld,%f\n",thread_num,((double)rand()/ (double)RAND_MAX));
        if(((double)rand()/ (double)RAND_MAX)<gamma){
                // printf("%ld\n",thread_num);
                mat[ki][kj]=-1*mat[ki][kj];
                E+=del_E;
        }
        B=1.00005*B;
        pthread_mutex_unlock(&mutexes[ki]);
        
    }
    
   return NULL;
}

void Usage(char prog_name[] /* in */) {
   fprintf(stderr, "usage: %s ", prog_name); 
   fprintf(stderr, "<n_side> <times> <threads_count>\n");
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
    // number of threads
    int thread_count;

    long thread;
    pthread_t* thread_handles;
    if (argc != 4) Usage(argv[0]); 
    Get_args(argv, &n,&times,&thread_count);
    thread_handles=malloc(thread_count*sizeof(pthread_t));

    load=times/thread_count+1;
    srand(time(0));
    vector_make(&mat,n);
    printf("matrix before simulation\n");
    print_mat(mat,n,"before_pthread_bin.pgm");

    mutexes = (pthread_mutex_t*)malloc(sizeof(pthread_mutex_t)*n);
    for(int i = 0; i < n; i++) {
        pthread_mutex_init(&mutexes[i],NULL);
    }
    clock_gettime(CLOCK_MONOTONIC, &begin);
    for (thread=0;thread<thread_count;thread++){
        printf("thread:%ld\n",thread);
        pthread_create(&thread_handles[thread], NULL, markov_ising_anneal, (void*)thread);
    }
    for (thread=0;thread<thread_count;thread++){
        pthread_join(thread_handles[thread], NULL);
    }
    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed = (end.tv_nsec - begin.tv_nsec)/1e9 + (end.tv_sec - begin.tv_sec); 

    for(int i = 0; i < n; i++) {
      pthread_mutex_destroy(&mutexes[i]);
    }
    // markov_ising_anneal(mat,n,thread_count);
    printf("matrix after simulation\n");
    // -1 as white, 1 as black
    print_mat(mat,n,"after_pthread_bin.pgm");
    printf("Elapsed time: %f\n", elapsed);
    for (int i = 0; i < n; i++) {
        free(mat[i]);
    }
    free(mat);
    return 0;
}