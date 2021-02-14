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
int * labels2 = NULL;
int * indexs = NULL;
int * sub_indexs = NULL;

char *filename = "100wsubmatrix.txt";


struct worker{
    float* data;
    int len;
    int n_simples;
    int n_features;
    int n_cluters;
    int *myindexs;
    int *mycentre;
    int bias;
    int *min_i;
    int *min_h;
    double *min_cijh;
    
};


void showData(float *d, int r, int c);
int search(int *p, int x, int n);
float euclidean_dist(float* a, float* b, int n_features);
void savedata(char *filename, int n_simples, int n_cluters);
int getdata(int k);
double getCosts(float *X, float *medoid, int n_simples, int n_features, int n_cluters);
void init(float *X, int n_simples, int n_features, int n_cluters);
void allocate(float *X, int n_simples, int n_features, int n_cluters);
void kmedoids(float *X, int n_simples, int n_features, int n_cluters, int max_iter);
void pam(float *X, int n_simples, int n_features, int n_cluters, int max_iter, int *myindexs);
void *sub_pam(void *param);
void pam_pal(float *X, int n_simples, int n_features, int n_cluters, int max_iter, int *myindexs, int thread_num);
void clara_init(float *X, int n_simples, int sub_n_simples, int n_features, int n_cluters, int iter);
void clara_allocate(float *X, int sub_n_simples, int n_features, int n_cluters);
void clara_pam(float *X, int n_simples, int n_features, int n_cluters, int max_iter, int *myindexs);
void clara_pam_pal(float *X, int n_simples, int n_features, int n_cluters, int max_iter, int *myindexs, int thread_num);
void clara(float *X, int n_simples, int sub_n_simples, int n_features, int n_cluters, int max_iter, int pam_num);
double allocate_atall(float *X, int n_simples, int n_features, int n_cluters);

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
    
    //kmedoids(data, lineNum, 100, K, 500);
    
    clara(data, lineNum, 2000, 100, K, 300, 8);
    
    //printf("====medoids====\n");
    //for(int i=0;i<K;i++) printf("%d ", centre[i]);
    //printf("\n\n\n====label====\n");
    //for(int i=0;i<lineNum;i++) printf("%d ", labels[i]);
    printf("\ncost=%f \n",getCosts(data, medoids, lineNum, 100, K));
    end = clock();
    duration = (double)(end - start) / CLOCKS_PER_SEC;    
    printf( "\nruns %f seconds\n", duration );
    
    savedata("results.txt", lineNum, K);
    
    /*??*/
    free(data);
    free(centre);
    free(medoids);
    free(labels);
    free(labels2);
    free(indexs);
    free(sub_indexs);
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

//find x in p[0-n]
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


void savedata(char *filename, int n_simples, int n_cluters)
{
    char temp[5000];
    FILE *f = NULL;
    if(NULL == (f = fopen(filename, "a")))
    {
        printf("cannot open %s", filename);
        exit(0);
    }
    char buf[100];
    fputs("====medoids====\n", f);
    for(int i=0;i<n_cluters;i++)
    {
        sprintf(buf, "%d", centre[i]);
        fputs(buf, f);
        fputs(" ", f);
    }

    fputs("\n\n====label====\n", f);
    for(int i=0;i<n_simples;i++)
    {
        sprintf(buf, "%d", labels[i]);
        fputs(buf, f);
        fputs(" ", f);
    }
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

double getCosts(float *X, float *medoid, int n_simples, int n_features, int n_cluters)
{
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
    return cost;
}

void init(float *X, int n_simples, int n_features, int n_cluters)
{
    centre = (int *) malloc(sizeof(int) * n_cluters);
    medoids = (float *) malloc(sizeof(float)*n_cluters*n_features);
    labels = (int *) malloc(sizeof(int)*n_simples);
    labels2 = (int *) malloc(sizeof(int)*n_simples);
    indexs = (int *) malloc(sizeof(int)*(n_simples - n_cluters));

    double temp = 0, mindis;
    int index = 0;
    /*随机生成一个数*/
    srand((unsigned int)time(NULL));  //设置随机种子
    int r = rand()%n_simples + 1;
    //r = 6;
    *centre = r;
    /*随机生成一个数*/
    memcpy(medoids, X + r*n_features, sizeof(float)*n_features);

    /*中心点完全随机

    for(int k=1;k<n_cluters;k++)
    {
        int Index=0;
        int x = rand()%n_simples;
        while(search(centre, x, k)) x = rand()%n_simples;
        
        Index = x;
        centre[k] = Index;
        memcpy(medoids + k*n_features, X + Index*n_features, sizeof(float)*n_features);
    }*/
    /**/
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
    /**/
 

}


void kmedoids(float *X, int n_simples, int n_features, int n_cluters, int max_iter)
{
    init(X, n_simples, n_features, n_cluters);

    for(int i=0;i<K;i++) printf("%d ", centre[i]);


    //pam(X, n_simples, n_features, n_cluters, max_iter, indexs);
    pam_pal(X, n_simples, n_features, n_cluters, max_iter, indexs, 2);
    

}

void allocate(float *X, int n_simples, int n_features, int n_cluters)
{
    for(int i=0;i<n_simples;i++)
    {
        double mindist = INT_MAX, dist;
        int flag=-1;
        for(int j=0;j<n_cluters;j++)
        {
            dist = euclidean_dist(X + i*n_features, medoids + j*n_features, n_features);
            if(dist<mindist)
            {
                mindist = dist;
                labels2[i] = labels[i];
                labels[i] = j;
                flag++;
            }
            else{
                if(flag<=0)
                {
                    labels2[i] = j;
                }
            }
        }
    }
    int n=0;
    for(int i=0;i<n_simples;i++)
    {
        if(!search(centre, i, n_cluters)) indexs[n++]=i;
    }

}


void pam(float *X, int n_simples, int n_features, int n_cluters, int max_iter, int *myindexs)
{
    allocate(X, n_simples, n_features, n_cluters);
    int min_h = 0, min_i = 0, iter = 0;
    double min_cijh;
    double C_ijh = 0;
    for(;iter<max_iter;iter++)
    {
        min_cijh = INT_MAX;
        for(int i=0;i<n_cluters;i++)
        {
            // o_i -> o_h  find minumize o_h
            for(int h=0;h<n_simples-n_cluters;h++)
            {
                C_ijh = 0;
                /*swap o_i and o_h*/
                int t = myindexs[h];
                myindexs[h] = centre[i];
                centre[i] = t;
                /*swap o_i and o_h*/
                for(int j=0;j<n_simples-n_cluters;j++)
                {
                    double d_jh = euclidean_dist(X + myindexs[j]*n_features, X + centre[i]*n_features, n_features);                  
                   
                    if(labels[myindexs[j]] == i)
                    {
                        double d_jj2 = euclidean_dist(X + myindexs[j]*n_features, X + centre[labels2[myindexs[j]]]*n_features, n_features);
                        double d_ji = euclidean_dist(X + myindexs[j]*n_features, X + myindexs[h]*n_features, n_features);
                         //case 1: o_j belong to o_i && d_jh>=d_jj2
                        if (d_jh>d_jj2||d_jh==d_jj2)
                        {
                            C_ijh += (d_jj2 - d_ji);
                        }
                         //case 2: o_j belong to o_i && d_jh<d_jj2
                        else{
                            C_ijh += (d_jh - d_ji);
                        }
                    }
                    else 
                    {
                        //case 3: o_j not belong to o_i && d_jh>=d_jj2
                        double d_jj2 = euclidean_dist(X + myindexs[j]*n_features, X + centre[labels[myindexs[j]]]*n_features, n_features);
                        if(d_jh>d_jj2||d_jh==d_jj2)
                        {
                            ;
                        }
                        //case 4: o_j not belong to o_i && d_jh<d_jj2
                        else{
                            
                            C_ijh += (d_jh - d_jj2);
                        }
                    }
                }
                /*swap o_i and o_h*/
                t = myindexs[h];
                myindexs[h] = centre[i];
                centre[i] = t;
                /*swap o_i and o_h*/
                
                // double d = euclidean_dist(X + centre[labels[myindexs[h]]]*n_features, X + myindexs[h]*n_features, n_features);
                // C_ijh -= d;

                if(C_ijh < min_cijh)
                {
                    min_i = i;
                    min_h = h;
                    min_cijh = C_ijh;
                }
            }
        }
        if(min_cijh > 0 || min_cijh == 0) break;
        
        /*swap minumize o_i and o_h*/
        int t = myindexs[min_h];
        myindexs[min_h] = centre[min_i];
        centre[min_i] = t;

        memcpy(medoids + min_i*n_features, X + t*n_features, sizeof(float)*n_features);
        allocate(X, n_simples, n_features, n_cluters);

        /*swap minumize o_i and o_h*/

    }

    printf("\nPAM has %d s iters!\n", iter);

}


void *sub_pam(void *param)
{
    struct worker *p = (struct worker*)param;
    float *X = p->data;
    int jobsnum = p->len;
    int n_simples = p->n_simples;
    int n_features = p->n_features;
    //float *medoid = p->medoid;
    int n_cluters = p->n_cluters;
    int *myindexs = p->myindexs;
    int *mycentre = p->mycentre;
    int bias = p->bias;
    double cost = 0, temp = 0, mindis;

    double min_cijh = INT_MAX;
    double C_ijh = 0;


    int min_h = 0, min_i = 0;
    for(int i=bias;i<(jobsnum+bias);i++)
    {
        // o_i -> o_h  find minumize o_h
        for(int h=0;h<n_simples-n_cluters;h++)
        {
            C_ijh = 0;
            
            /*swap o_i and o_h*/
            int t = myindexs[h];
            myindexs[h] = mycentre[i];
            mycentre[i] = t;
            /*swap o_i and o_h*/
            for(int j=0;j<n_simples-n_cluters;j++)
            {
                double d_jh = euclidean_dist(X + myindexs[j]*n_features, X + mycentre[i]*n_features, n_features);                  
                
                if(labels[myindexs[j]] == i)
                {
                    double d_jj2 = euclidean_dist(X + myindexs[j]*n_features, X + mycentre[labels2[myindexs[j]]]*n_features, n_features);
                    double d_ji = euclidean_dist(X + myindexs[j]*n_features, X + myindexs[h]*n_features, n_features);
                    //case 1: o_j belong to o_i && d_jh>=d_jj2
                    if (d_jh>d_jj2||d_jh==d_jj2)
                    {
                        C_ijh += (d_jj2 - d_ji);
                    }
                    //case 2: o_j belong to o_i && d_jh<d_jj2
                    else{
                        C_ijh += (d_jh - d_ji);
                    }
                }
                else 
                {
                    //case 3: o_j not belong to o_i && d_jh>=d_jj2
                    double d_jj2 = euclidean_dist(X + myindexs[j]*n_features, X + mycentre[labels[myindexs[j]]]*n_features, n_features);
                    if(d_jh>d_jj2||d_jh==d_jj2)
                    {
                        ;
                    }
                    //case 4: o_j not belong to o_i && d_jh<d_jj2
                    else{
                        
                        C_ijh += (d_jh - d_jj2);
                    }

                }

            }
            /*swap o_i and o_h*/
            t = myindexs[h];
            myindexs[h] = mycentre[i];
            mycentre[i] = t;
            /*swap o_i and o_h*/
            
            //double d = euclidean_dist(X + mycentre[labels[myindexs[h]]]*n_features, X + myindexs[h]*n_features, n_features);
            //C_ijh -= d;

            if(C_ijh < min_cijh)
            {
                min_i = i;
                min_h = h;
                min_cijh = C_ijh;
            }
        }
    }
    *(p->min_i) = min_i;
    *(p->min_h) = min_h;
    *(p->min_cijh) = min_cijh;

}

void pam_pal(float *X, int n_simples, int n_features, int n_cluters, int max_iter, int *myindexs, int thread_num)
{
    allocate(X, n_simples, n_features, n_cluters);

    int jobs = n_cluters/thread_num;
    struct worker *myworker = (struct worker *)malloc(sizeof(struct worker)*thread_num);
    double *worker_cijh = (double*)malloc(sizeof(double)*thread_num);
    int *worker_min_i = (int*)malloc(sizeof(int)*thread_num);
    int *worker_min_h = (int*)malloc(sizeof(int)*thread_num);

    int *worker_indexs = (int*)malloc(sizeof(int)*(n_simples-n_cluters)*thread_num);
    int *worker_centre = (int*)malloc(sizeof(int)*(n_cluters)*thread_num);

    pthread_t *threads = (pthread_t*)malloc(sizeof(pthread_t)*thread_num);

    int min_h = 0, min_i = 0, iter = 0;
    double min_cijh;
    double C_ijh = 0;
    for(;iter<max_iter;iter++)
    {
        min_cijh = INT_MAX;
        
        for(int i=0;i<thread_num;i++)
        {
            if(i<thread_num-1||thread_num==1)
            {
                myworker[i].len = jobs;
            }
            else{
                myworker[i].len = n_cluters - jobs*(i);
            }
            myworker[i].data = X;
            //myworker[i].medoid = medoid;
            myworker[i].n_simples = n_simples;
            myworker[i].n_cluters = n_cluters;
            myworker[i].n_features = n_features;
            myworker[i].min_cijh = worker_cijh + i;
            myworker[i].min_i = worker_min_i + i;
            myworker[i].min_h = worker_min_h + i;
            myworker[i].bias = jobs*(i);
            
            memcpy(worker_indexs + i*(n_simples-n_cluters), myindexs, sizeof(int)*(n_simples-n_cluters));
            myworker[i].myindexs = worker_indexs + i*(n_simples-n_cluters);
            
            memcpy(worker_centre + i*(n_cluters), centre, sizeof(int)*(n_cluters));
            myworker[i].mycentre = worker_centre + i*(n_cluters);
            
            
            //myworker[i].

            //pthread_create(&threads[i], 0, &get_sub_cost, (void*)(myworker + i));
            pthread_create(&threads[i], 0, &sub_pam, (void*)(myworker + i));
            
        }

        for(int i=0;i<thread_num;i++)
        {
            pthread_join(threads[i], 0);
        }
            
        min_cijh = *(worker_cijh);
        min_i = *(worker_min_i);
        min_h = *(worker_min_h);
        for(int p=1;p<thread_num;p++)
        {
            if(*(worker_cijh+p)<min_cijh)
            {
                min_cijh = *(worker_cijh+p);
                min_i = *(worker_min_i+p);
                min_h = *(worker_min_h+p);
            }
        }

        if(min_cijh > 0 || min_cijh == 0) break;
        
        /*swap minumize o_i and o_h*/
        int t = myindexs[min_h];
        myindexs[min_h] = centre[min_i];
        centre[min_i] = t;

        memcpy(medoids + min_i*n_features, X + t*n_features, sizeof(float)*n_features);
        allocate(X, n_simples, n_features, n_cluters);

        /*swap minumize o_i and o_h*/

    }
    free(myworker);
    free(worker_cijh);
    free(worker_min_i);
    free(worker_min_h);
    free(worker_indexs);
    free(worker_centre);
    free(threads);

    printf("\nPAM has %d s iters!\n", iter);

}


void clara_init(float *X, int n_simples, int sub_n_simples, int n_features, int n_cluters, int iter)
{
    if(!iter)
    {
        centre = (int *) malloc(sizeof(int) * n_cluters);
        medoids = (float *) malloc(sizeof(float)*n_cluters*n_features);
        labels = (int *) malloc(sizeof(int)*n_simples);
        labels2 = (int *) malloc(sizeof(int)*n_simples);
        indexs = (int *) malloc(sizeof(int)*(sub_n_simples - n_cluters));
        sub_indexs = (int *) malloc(sizeof(int)*(sub_n_simples));

    }

    double temp = 0, mindis;
    int index = 0;

    int *myindexs = sub_indexs;
    
    //如果用clara样本，可以这样
    /**/
    for(int i=0;i<sub_n_simples;i++)
    {
        int r = rand()%n_simples;
        while(search(myindexs, r, i)) r = rand()%n_simples;
        myindexs[i] = r;

    }
    
    /*
    //如果用全部样本，可以这样
    for(int i=0;i<sub_n_simples;i++)
    {
        myindexs[i] = i;
    }*/


    /*随机生成一个数*/
    srand((unsigned int)time(NULL));  //设置随机种子
    int r = rand()%sub_n_simples;
    printf("\n----r=%d----\n", myindexs[r]);
    //r = 6;
    *centre = myindexs[r];
    /*随机生成一个数*/
    memcpy(medoids, X + myindexs[r]*n_features, sizeof(float)*n_features);


    /*中心点完全随机

    for(int k=1;k<n_cluters;k++)
    {
        int Index=0;
        int x = rand()%sub_n_simples;
        while(search(centre, x, k)) x = rand()%n_simples;
        
        Index = myindexs[x];
        centre[k] = Index;
        memcpy(medoids + k*n_features, X + Index*n_features, sizeof(float)*n_features);
    }*/

    /*
    //kmeans++*/
    for(int k=1;k<n_cluters;k++)
    {   
        int Index=0;
        double Maxdis=-1;
        for(int m=0;m<sub_n_simples;m++)
        {      
            index = -1;
            mindis = INT_MAX;
            for(int i=0;i<k;i++)
            {
                if(search(centre, myindexs[m], k)) continue;
                temp = euclidean_dist(X + myindexs[m]*n_features, medoids + i*n_features, n_features);
                if(temp<mindis)
                {
                     mindis = temp;
                     index = myindexs[m];
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
    
    //free(myindexs);

}

void clara_allocate(float *X, int sub_n_simples, int n_features, int n_cluters)
{
    for(int i=0;i<sub_n_simples;i++)
    {
        double mindist = INT_MAX, dist;
        int flag=-1;
        for(int j=0;j<n_cluters;j++)
        {
            dist = euclidean_dist(X + sub_indexs[i]*n_features, medoids + j*n_features, n_features);
            if(dist<mindist)
            {
                mindist = dist;
                labels2[sub_indexs[i]] = labels[sub_indexs[i]];
                labels[sub_indexs[i]] = j;
                flag++;
            }
            else{
                if(flag<=0)
                {
                    labels2[sub_indexs[i]] = j;
                }
            }
        }
    }
    int n=0;
    for(int i=0;i<sub_n_simples;i++)
    {
        if(!search(centre, sub_indexs[i], n_cluters)) indexs[n++]=sub_indexs[i];
    }

}


void clara_pam(float *X, int n_simples, int n_features, int n_cluters, int max_iter, int *myindexs)
{
    clara_allocate(X, n_simples, n_features, n_cluters);
    int min_h = 0, min_i = 0, iter = 0;
    double min_cijh;
    double C_ijh = 0;
    for(;iter<max_iter;iter++)
    {
        min_cijh = INT_MAX;
        for(int i=0;i<n_cluters;i++)
        {
            // o_i -> o_h  find minumize o_h
            for(int h=0;h<n_simples-n_cluters;h++)
            {
                C_ijh = 0;
                /*swap o_i and o_h*/
                int t = myindexs[h];
                myindexs[h] = centre[i];
                centre[i] = t;
                /*swap o_i and o_h*/
                for(int j=0;j<n_simples-n_cluters;j++)
                {
                    double d_jh = euclidean_dist(X + myindexs[j]*n_features, X + centre[i]*n_features, n_features);                  
                   
                    if(labels[myindexs[j]] == i)
                    {
                        double d_jj2 = euclidean_dist(X + myindexs[j]*n_features, X + centre[labels2[myindexs[j]]]*n_features, n_features);
                        double d_ji = euclidean_dist(X + myindexs[j]*n_features, X + myindexs[h]*n_features, n_features);
                         //case 1: o_j belong to o_i && d_jh>=d_jj2
                        if (d_jh>d_jj2||d_jh==d_jj2)
                        {
                            C_ijh += (d_jj2 - d_ji);
                        }
                         //case 2: o_j belong to o_i && d_jh<d_jj2
                        else{
                            C_ijh += (d_jh - d_ji);
                        }
                    }
                    else 
                    {
                        //case 3: o_j not belong to o_i && d_jh>=d_jj2
                        double d_jj2 = euclidean_dist(X + myindexs[j]*n_features, X + centre[labels[myindexs[j]]]*n_features, n_features);
                        if(d_jh>d_jj2||d_jh==d_jj2)
                        {
                            ;
                        }
                        //case 4: o_j not belong to o_i && d_jh<d_jj2
                        else{
                            
                            C_ijh += (d_jh - d_jj2);
                        }
                    }
                }
                /*swap o_i and o_h*/
                t = myindexs[h];
                myindexs[h] = centre[i];
                centre[i] = t;
                /*swap o_i and o_h*/
                
                // double d = euclidean_dist(X + centre[labels[myindexs[h]]]*n_features, X + myindexs[h]*n_features, n_features);
                // C_ijh -= d;

                if(C_ijh < min_cijh)
                {
                    min_i = i;
                    min_h = h;
                    min_cijh = C_ijh;
                }
            }
        }
        if(min_cijh > 0 || min_cijh == 0) break;
        
        /*swap minumize o_i and o_h*/
        int t = myindexs[min_h];
        myindexs[min_h] = centre[min_i];
        centre[min_i] = t;

        memcpy(medoids + min_i*n_features, X + t*n_features, sizeof(float)*n_features);
        clara_allocate(X, n_simples, n_features, n_cluters);

        /*swap minumize o_i and o_h*/

    }

    printf("\nPAM has %d s iters!\n", iter);

}

void clara_pam_pal(float *X, int n_simples, int n_features, int n_cluters, int max_iter, int *myindexs, int thread_num)
{
    clara_allocate(X, n_simples, n_features, n_cluters);

    int jobs = n_cluters/thread_num;
    
    struct worker *myworker = (struct worker *)malloc(sizeof(struct worker)*thread_num);
    double *worker_cijh = (double*)malloc(sizeof(double)*thread_num);
    int *worker_min_i = (int*)malloc(sizeof(int)*thread_num);
    int *worker_min_h = (int*)malloc(sizeof(int)*thread_num);

    int *worker_indexs = (int*)malloc(sizeof(int)*(n_simples-n_cluters)*thread_num);
    int *worker_centre = (int*)malloc(sizeof(int)*(n_cluters)*thread_num);

    pthread_t *threads = (pthread_t*)malloc(sizeof(pthread_t)*thread_num);

    int min_h = 0, min_i = 0, iter = 0;
    double min_cijh;
    double C_ijh = 0;
    for(;iter<max_iter;iter++)
    {
        min_cijh = INT_MAX;
        
        for(int i=0;i<thread_num;i++)
        {
            if(i<thread_num-1||thread_num==1)
            {
                myworker[i].len = jobs;
            }
            else{
                myworker[i].len = n_cluters - jobs*(i);
            }
            myworker[i].data = X;
            //myworker[i].medoid = medoid;
            myworker[i].n_simples = n_simples;
            myworker[i].n_cluters = n_cluters;
            myworker[i].n_features = n_features;
            myworker[i].min_cijh = worker_cijh + i;
            myworker[i].min_i = worker_min_i + i;
            myworker[i].min_h = worker_min_h + i;
            myworker[i].bias = jobs*(i);
            
            memcpy(worker_indexs + i*(n_simples-n_cluters), myindexs, sizeof(int)*(n_simples-n_cluters));
            myworker[i].myindexs = worker_indexs + i*(n_simples-n_cluters);
            
            memcpy(worker_centre + i*(n_cluters), centre, sizeof(int)*(n_cluters));
            myworker[i].mycentre = worker_centre + i*(n_cluters);
            
            
            //myworker[i].

            //pthread_create(&threads[i], 0, &get_sub_cost, (void*)(myworker + i));
            pthread_create(&threads[i], 0, &sub_pam, (void*)(myworker + i));
            
        }

        for(int i=0;i<thread_num;i++)
        {
            pthread_join(threads[i], 0);
        }
            
        min_cijh = *(worker_cijh);
        min_i = *(worker_min_i);
        min_h = *(worker_min_h);
        for(int p=1;p<thread_num;p++)
        {
            if(*(worker_cijh+p)<min_cijh)
            {
                min_cijh = *(worker_cijh+p);
                min_i = *(worker_min_i+p);
                min_h = *(worker_min_h+p);
            }
        }

        if(min_cijh > 0 || min_cijh == 0) break;
        
        /*swap minumize o_i and o_h*/
        int t = myindexs[min_h];
        myindexs[min_h] = centre[min_i];
        centre[min_i] = t;

        memcpy(medoids + min_i*n_features, X + t*n_features, sizeof(float)*n_features);
        clara_allocate(X, n_simples, n_features, n_cluters);

        /*swap minumize o_i and o_h*/

    }

    /* 至今用clara并行处理的时候不明白为什么释放内存会出错 */
    free(myworker);
    free(worker_cijh);
    free(worker_min_i);
    free(worker_min_h);
    free(worker_indexs);
    free(worker_centre);
    free(threads);

    printf("\nPAM has %d s iters!\n", iter);

}
/*
CLARA (Clustering LARge Applications，大型应用中的聚类方法)(Kaufmann and Rousseeuw in 1990)
*/
void clara(float *X, int n_simples, int sub_n_simples, int n_features, int n_cluters, int max_iter, int pam_num)
{
    int *mincentre = (int *) malloc(sizeof(int) * n_cluters);
    float *minmedoids = (float *) malloc(sizeof(float)*n_cluters*n_features);
    int *minlabels = (int *) malloc(sizeof(int)*n_simples);
    
    double mincost = INT_MAX, curcost;

    for(int i=0;i<pam_num;i++)
    {
        clara_init(X, n_simples, sub_n_simples, n_features, n_cluters, i);


        //clara_pam(X, sub_n_simples, n_features, n_cluters, max_iter, indexs);
        clara_pam_pal(X, sub_n_simples, n_features, n_cluters, max_iter, indexs, 25);
        curcost = allocate_atall(X, n_simples, n_features, n_cluters);
        
        //curcost = getCosts(data, medoids, n_simples, 100, K);
        if(curcost<mincost)
        {
            mincost = curcost;
            memcpy(mincentre, centre, sizeof(int)*n_cluters);
            memcpy(minmedoids, medoids, sizeof(float)*n_cluters*n_features);
            memcpy(minlabels, labels, sizeof(int)*n_simples);
        }
        printf("\niter %d cost=%f \n", i+1, curcost);
    }
    memcpy(centre, mincentre, sizeof(int)*n_cluters);
    memcpy(medoids, minmedoids, sizeof(float)*n_cluters*n_features);
    memcpy(labels, minlabels, sizeof(int)*n_simples);

    free(mincentre);
    free(minmedoids);
    free(minlabels);
}


double allocate_atall(float *X, int n_simples, int n_features, int n_cluters)
{
    double cost = 0;
    for(int i=0;i<n_simples;i++)
    {
        double mindist = INT_MAX, dist;
        int flag=-1;
        for(int j=0;j<n_cluters;j++)
        {
            dist = euclidean_dist(X + i*n_features, medoids + j*n_features, n_features);
            if(dist<mindist)
            {
                mindist = dist;
                labels2[i] = labels[i];
                labels[i] = j;
                flag++;
            }
            else{
                if(flag<=0)
                {
                    labels2[i] = j;
                }
            }
        }
        cost += mindist;
    }
    

}
