#include <stdio.h>
#include <sched.h>
#include <sys/sysinfo.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <unistd.h>
#include <pthread.h>
#include <time.h>
#include <cmath>
#include <cstring>

#define EPS 1e-15


struct Args{	
	double* a = nullptr;
	double* x = nullptr;
	int n = 0;
	int m = 0;
	int p = 0;
	int k = 0;
	int r = 0;
    int s = 0; 
    const char* filename = nullptr;

	double cpu_time = 0.0;
	int status = 0;
	
};


double get_full_time();
double get_cpu_time();
double func(int s, int n, int i, int j);
void init_a(double* a, int n, int s);
int read(double* a, int n, const char* name);
void print(double* a, int l, int n, int r);
double matrix_norm(double* a, int n);
double norm1(double* a, double* x, double* id, int n);
double norm2(double* a, double* x, double* id, int n);
void multiply(double* a, double* b,double* c, int n1, int m1, int n2, int m2);
void* thread_func(void* ptr);


template<class T>
void reduce_sum(int total_threads, T* a=nullptr, int n=0){
	static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
	static pthread_cond_t cond_in = PTHREAD_COND_INITIALIZER;
	static pthread_cond_t cond_out = PTHREAD_COND_INITIALIZER;
	static int t_in = 0;
	static int t_out = 0;
	static T* res = nullptr;
	int i;

	pthread_mutex_lock(&mutex);
	if(res == nullptr) res = a;
	else{
		for(i = 0; i < n; i++){
			res[i] += a[i];
		}
	}
	t_in++;
	if(t_in >= total_threads){
		t_out = 0;
		pthread_cond_broadcast(&cond_in);
	}
	else{
		while(t_in < total_threads) pthread_cond_wait(&cond_in, &mutex);
	}
	if(res != a){
		for(i = 0; i < n; i++){
			a[i] = res[i];
		}
	}
	t_out++;
	if(t_out >= total_threads){
		t_in = 0;
		res = nullptr;
		pthread_cond_broadcast(&cond_out);
	}
	else{
		while(t_out < total_threads) pthread_cond_wait(&cond_out, &mutex);
	}
	pthread_mutex_unlock(&mutex);
}
