//#include "readCSVFile.h"
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<limits.h>
#include<pthread.h>

#define K 100
#define M 10
//#define 

float * data = NULL;
int * centre = NULL;
float * medoids = NULL;
int * labels = NULL;

char *filename = "submatrix.txt";

struct worker{
    float* start;
    int len;
    float *medoid;
    //int n_simples;
    int n_features;
    int n_cluters;
    double *mycost;
};


void showData(float *d, int r, int c);
int search(int *p, int x, int n);
float euclidean_dist(float* a, float* b, int n_features);
int getdata(int k);
void* get_sub_cost(void* param);
double getCosts(float *X, float *medoid, int n_simples, int n_features, int n_cluters, int n_jobs);
void init(float *X, int n_simples, int n_features, int n_cluters);
void kmedoids(float *X, int n_simples, int n_features, int n_cluters, int max_iter);

double MINC = 0;

int main()
{
    
    double  duration, d2; 
    int lineNum = getdata(K); 
    // printf("%d\n", lineNum);
    // return 0;
    clock_t start = clock(), end;
    /**/
    //showData(data, lineNum, K);
    
    kmedoids(data, lineNum, 100, K, 20);
    
    
    printf("====medoids====\n");
    for(int i=0;i<K;i++) printf("%d ", centre[i]);
    printf("\n\n\n====label====\n");
    for(int i=0;i<lineNum;i++) printf("%d ", labels[i]);
    printf("%f \n",MINC);
    end = clock();
    duration = (double)(end - start) / CLOCKS_PER_SEC;    
    printf( "\nruns %f seconds\n", duration );
    
    
    free(data);
    return 0;
}


void showData(float *d, int r, int c)
{
    for(int i=0;i<r;i++)
    {
        printf("no.%d\n", i+1);
        for(int j=0;j<c;j++)
            printf("%f ", (d + i * c)[j]);
        printf("\n\n\n");
    }
}

int search(int *p, int x, int n)
{
    for(int i=0;i<n;i++)
        if(p[i]==x) return 1;
    
    return 0;
}

float euclidean_dist(float* a, float* b, int n_features)
{
    double result = 0, tmp = 0;
    for(int i=0;i<n_features;i++)
    {
        tmp = (a[i] - b[i]);
        result += tmp*tmp;
    }
    return (float)sqrt(result);
}


int getdata(int k)
{
    char temp[5000];
    FILE *f = NULL;
    int lineNum = 0;
    if(NULL == (f = fopen(filename, "r")))
    {
        printf("cannot open %s", filename);
        exit(0);
    }
    
    
    while(fgets(temp, 2048, f) != NULL) lineNum++;

    /* */
    if(NULL == (data = (float *) malloc(sizeof(float) * lineNum * k)))
    {
        printf("Allocated memory failure!\n");
        exit(0);
    }
   
    fseek(f, 0, SEEK_SET);  //
    char *p_num; //= (char *)malloc(sizeof(char) * 50);//free(p_num) 
    int curline=0;
    while(fgets(temp, 2048, f) != NULL)
    {
        //printf("%s", temp);
        /**/
        //char *line = (char *) malloc(5000);
        //memcpy(line, temp, sizeof(char) * strlen(temp));
        p_num = strtok(temp, ",");
        int n=0;
        while(p_num != NULL)
        {
             
            p_num = strtok(NULL, ",");
            if(p_num!=NULL)
                (data + curline * K) [n++] = atof(p_num); 
        }
        curline++;
        //free(line);  
    }/**/
 
    free(p_num);
    
    
    if(fclose(f))
        printf("%scanot closed!\n", filename);
    else
        printf("%s closed!\n", filename);
    
    return lineNum;
}

void* get_sub_cost(void* param)
{
    struct worker *p = (struct worker*)param;
    float *X = p->start;
    int n_simples = p->len;
    int n_features = p->n_features;
    float *medoid = p->medoid;
    int n_cluters = p->n_cluters;
    double cost = 0, temp = 0, mindis;
    for(int m=0;m<n_simples;m++)
    {
        mindis = euclidean_dist(X + m*n_features, medoid, n_features);
        for(int i=1;i<n_cluters;i++)
        {
             temp = euclidean_dist(X + m*n_features, medoid + i*n_features, n_features);
             if(temp<mindis) mindis = temp;
        }
        cost += mindis;
    }
    *(p->mycost) = cost;
    //return (void*)cost;
}

double getCosts(float *X, float *medoid, int n_simples, int n_features, int n_cluters, int thread_num)
{
    double cost = 0, temp = 0, mindis;
    int jobs = n_simples/thread_num;
    struct worker *myworker = (struct worker *)malloc(sizeof(struct worker)*thread_num);
    double *worker_cost = (double*)malloc(sizeof(float)*thread_num);

    pthread_t *threads = (pthread_t*)malloc(sizeof(pthread_t)*thread_num);

    for(int i=0;i<thread_num;i++)
    {
        if(i<thread_num-1||thread_num==1)
        {
            myworker[i].len = jobs;

        }
        else{
            myworker[i].len = n_simples - jobs*(i);
        }
        myworker[i].start = (X + i*jobs*n_features);
        myworker[i].medoid = medoid;
        myworker[i].n_cluters = n_cluters;
        myworker[i].n_features = n_features;
        myworker[i].mycost = worker_cost + i;
        //myworker[i].

        //pthread_create(&threads[i], 0, &get_sub_cost, (void*)(myworker + i));
        pthread_create(&threads[i], 0, &get_sub_cost, (void*)(myworker + i));
        
    }
    for(int i=0;i<thread_num;i++)
    {
        pthread_join(threads[i], 0);
    }

    for(int i=0;i<thread_num;i++)
    {
        cost += worker_cost[i];
    }
    return cost;
    /*
    for(int m=0;m<n_simples;m++)
    {
        mindis = euclidean_dist(X + m*n_features, medoid, n_features);
        for(int i=1;i<n_cluters;i++)
        {
             temp = euclidean_dist(X + m*n_features, medoid + i*n_features, n_features);
             if(temp<mindis) mindis = temp;
        }
        cost += mindis;
    }
    return cost;
    */
}

void init(float *X, int n_simples, int n_features, int n_cluters)
{
    centre = (int *) malloc(sizeof(int) * n_cluters);
    medoids = (float *) malloc(sizeof(float)*n_cluters*n_features);
    labels = (int *) malloc(sizeof(int)*n_simples);


    double temp = 0, mindis;
    int index = 0;
    /*随机生成一个数*/
    //srand((unsigned int)time(NULL));  //设置随机种子
    int r = rand()%n_simples + 1;
    r = 6;
    *centre = r;
    /*随机生成一个数*/
    memcpy(medoids, X + r*n_features, sizeof(float)*n_features);


    //kmeans++
    for(int k=1;k<n_cluters;k++)
    {   
        int Index=0;
        double Maxdis=-1;
        for(int m=0;m<n_simples;m++)
        {      
            index = -1;
            mindis = INT_MAX;
            for(int i=0;i<k;i++)
            {
                if(search(centre, m, k)) continue;
                temp = euclidean_dist(X + m*n_features, medoids + i*n_features, n_features);
                if(temp<mindis)
                {
                     mindis = temp;
                     index = m;
                }
            }
            if(Maxdis<mindis&&index!=-1)
            {
                Maxdis = mindis;
                Index = index;
            }
            
        }

        //
        centre[k] = Index;
        memcpy(medoids + k*n_features, X + Index*n_features, sizeof(float)*n_features);

    }
 

}


void kmedoids(float *X, int n_simples, int n_features, int n_cluters, int max_iter)
{
    init(X, n_simples, n_features, n_cluters);
    
    for(int i=0;i<K;i++) printf("%d ", centre[i]);

    //exit(0);


    //PAM start
    
    int convergence = 0;
    float *temp = (float *) malloc(sizeof(float)*n_cluters*n_features);
    int iter = 0;
    for(;iter<max_iter&&!convergence;iter++)
    {
        convergence = 1;
        for(int p=0;p<n_cluters;p++)
        {
            double mincosts = getCosts(X, medoids, n_simples, n_features, n_cluters, 4);
            memcpy(temp, medoids, sizeof(float)*n_cluters*n_features);
            for(int m=0;m<n_simples;m++)
            {
                if(search(centre, m, n_cluters)) continue;
                memcpy(temp + p*n_features, X + m*n_features, sizeof(float)*n_features);
                double cost = getCosts(X, temp, n_simples, n_features, n_cluters, 4);
                if(cost<mincosts)
                {
                    mincosts = cost;
                    MINC = cost;
                    memcpy(medoids + p*n_features, X + m*n_features, sizeof(float)*n_features);
                    centre[p] = m;
                    convergence = 0;
                }
                else{
                    memcpy(temp + p*n_features, medoids + p*n_features, sizeof(float)*n_features);
                }
            }
        }
    }

    for(int i=0;i<n_simples;i++)
    {
        double mindist = INT_MAX, dist;
        for(int j=0;j<n_cluters;j++)
        {
            dist = euclidean_dist(X + i*n_features, medoids + j*n_features, n_features);
            if(dist<mindist)
            {
                mindist = dist;
                labels[i] = j;
            }
        }
    }

    printf("\nkmedoids run %d iter.\n\n", iter);

}


/*
CLARA (Clustering LARge Applications，大型应用中的聚类方法)(Kaufmann and Rousseeuw in 1990)
*/

void scara()
{


}
