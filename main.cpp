#include"func.h"

int main(int argc, char** argv){
	double *a, *x;
	const char* filename = nullptr;
	int n, m, p, r, s;
  int task = 20;
  double full_time, norm_time, r1, r2;

  if(argc != 7 && argc != 6)
    {
      printf ("Usage: %s n m p r s [file]\n", argv[0]);
      return -1;
    }
  if(sscanf(argv[1], "%d", &n) != 1 || sscanf(argv[2], "%d", &m) != 1 || sscanf(argv[3], "%d", &p) != 1
           || sscanf(argv[4], "%d", &r) != 1 || sscanf(argv[5], "%d", &s) != 1)
    {
      printf ("Incorrect input args\n");
      return -1;
    }

  if(argc == 7)
    {
      filename = argv[6];
      if(s != 0)
        {
          return -1;
        }
    }

  a = new double[n*n];
  x = new double[n*n];
  Args* arg = new Args[p];
  pthread_t* tid = new pthread_t[p];
  for(int k = 0; k < p; k++)
  {
    arg[k].a = a;
    arg[k].x = x;
    arg[k].k = k;
    arg[k].n = n;
    arg[k].m = m;
    arg[k].p = p;
    arg[k].r = r;
    arg[k].s = s;
    arg[k].filename = filename;
  }
  full_time = get_full_time();
  for(int k = 1; k < p; k++){
    if(pthread_create(tid + k, 0, thread_func, arg+k)){
      printf("Can not create thread number: %d\n",k);
      abort();
    }
  }
  thread_func(arg);
  for(int k = 1; k < p; k++){
    pthread_join(tid[k], 0);
  }
  full_time = get_full_time() - full_time;


  for (int k = 0; k < p; k++){
    if(arg[k].status < 0){
      printf ("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d P = %d\n",
        argv[0], task, -1.0, -1.0, 0.0, 0.0, s, n, m, p);
      delete[] a;
      delete[] x;
      delete[] arg;
      delete[] tid;
      return -1;
    }
  }

  for(int k = 0; k < p; k++){
    printf("Thread %d = %.2f\n", k, arg[k].cpu_time);
  }

  if(argc == 6) init_a(a, n, s);
  else read(a, n, filename);

  double* id = new double[n*n];
  norm_time = clock();
  r1 = norm1(a, x, id, n);
  r2 = norm2(a, x, id, n);
  norm_time = (clock() - norm_time) / CLOCKS_PER_SEC;
  delete[] id;

  printf("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d P = %d\n",
    argv[0], task, r1, r2, full_time, norm_time, s, n, m, p);
  delete[] a;
  delete[] x;
  delete[] arg;
  delete[] tid;
	return 0;
}