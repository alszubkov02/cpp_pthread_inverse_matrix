#include"func.h"

void get_block (double* a, double* b, int n, int m, int row, int col)
{
  int k = n / m;
  int l = n - k*m;
  int h = (row < k ? m : l);
  int w = (col < k ? m : l);
  double* block = a + m * (row * n + col);
  for(int i = 0; i < h; i ++)
    {
      for(int j = 0; j < w; j ++)
        {
          b[i*w+j] = block[i*n+j];
        }
    }
}

void put_block (double* a, double* b, int n, int m, int row, int col)
{
  int k = n / m;
  int l = n - k*m;
  int h = (row < k ? m : l);
  int w = (col < k ? m : l);
  double* block = a + m * (row * n + col);
  for(int i = 0; i < h; i ++)
    {
      for(int j = 0; j < w; j ++)
        {
          block[i * n + j] = b[i * w + j];
        }
    }
}

void init_a_thread(double* a, int n, int m, int s, int k, int p) {
    for (int j = k*m; j < n; j += p*m) {
        int h = (j + m < n) ? m : n - j;
        for (int i = 0; i < n; ++i) {
            for (int v = 0; v < h; ++v) {
                a[i*n + j + v] = func(s, n, i, j + v);
            }   
        }
    }
    reduce_sum<int>(p);
}

void init_x_thread(double* x, int n, int m, int k, int p) {
    for (int j = k*m; j < n; j += p*m) {
        int h = (j + m < n) ? m : n - j;
        for (int i = 0; i < n; ++i) {
            for (int v = 0; v < h; ++v) {
                x[i*n + j + v] = (i == j + v ? 1.0 : 0.0);
            }   
        }
    }
    reduce_sum<int>(p);
}


void copy_matrix(double *a, double *b, int n, int m){
  for (int i = 0; i < n; ++i){
    for (int j = 0; j < m; ++j){
      a[i*m+j] = b[i*m+j];
    }
  }
}



void swap_rows (double* a, int n, int i, int j)
{
  for(int k = 0; k < n; k++)
    {
      std::swap(a[i * n + k], a[j * n + k]);
    }
}

int jordan_reverse (double* a, double* x, int n, double norm)
{
	double lead, temp;
  int lead_j;

  for(int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      x[i*n+j] = (i == j) ? 1.0 : 0.0;
    }
  }

  for (int k = 0; k < n; k++) 
    {
      lead = a[k*n+k];
      lead_j = k;

      for(int j = k; j < n; j++)
        {
          if(fabs(a[j*n+k]) > fabs(lead))
            {
              lead = a[j*n+k];
              lead_j = j;
            }
        }
      if(fabs(lead) < EPS*norm)
        {
          return -1;
        }

      swap_rows(a,n,k,lead_j);
      swap_rows(x,n,k,lead_j);

      temp = a[k*n+k];
      for (int j = k; j < n; j++)
        {
          a[k*n+j] /= temp;
        }
      for (int j = 0; j < n; j++)
        {
          x[k*n+j] /= temp;
        }

      for (int i = 0;i < n; i++)
        {
          if(i == k) continue;
          temp = a[i*n+k];
          for(int j = k; j < n; j++) a[i*n+j] -= temp * a[k*n+j];
          for(int j = 0; j < n; j++) x[i*n+j] -= temp * x[k*n+j];
        }
    }

  return 0;
}


void swap_rows_for_blocks (double* a, int n, int m, int i, int j, double* A, double* B)
{
  int k = n / m;
  for(int s = 0; s < k; s++)
    {
      get_block(a, A, n, m, i, s);
      get_block(a, B, n, m, j, s);
      put_block(a, B, n, m, i, s);
      put_block(a, A, n, m, j, s);
    }
}

void swap_columns_for_blocks (double* a, int n, int m, int i, int j, double* A, double* B)
{
  int k = n / m;
  for(int s = 0; s < k; s++)
    {
      get_block(a, A, n, m, s, i);
      get_block(a, B, n, m, s, j);
      put_block(a, B, n, m, s, i);
      put_block(a, A, n, m, s, j);
    }
}



int solve(double* a, double* x, int n, int m, int p, int k, double norm, double* A, double* B, double* C, double* copy_row1, double* copy_row2)
{
    int d = n / m;
    int l = n - d * m;
    int h = l ? d + 1 : d;

    for (int t = 0; t < d; ++t) {
        if(t%p != k){
            copy_matrix(copy_row1, a+t*m*n, m, n);
            copy_matrix(copy_row2, x+t*m*n, m, n);
        }
        reduce_sum<int>(p);

        if(t%p == k) get_block(a, A, n, m, t, t);
        else get_block(copy_row1, A, n, m, 0, t);
        int res = jordan_reverse(A, B, m, norm);
        if (res == -1) {
            return -1;
        }

        reduce_sum<int>(p);

        for(int u = t+1; u < d; u++){
            if(t%p == k){
                get_block(a, A, n, m, t, u);
                multiply(B, A, C, m, m, m, m);
                put_block(a, C, n, m, t, u);
            }
            else{
                get_block(copy_row1, A, n, m, 0, u);
                multiply(B, A, C, m, m, m, m);
                put_block(copy_row1, C, n, m, 0, u);
            }
        }
        if(l > 0){
            if(t%p == k){
                get_block(a, A, n, m, t, d);
                multiply(B, A, C, m, m, m, l);
                put_block(a, C, n, m, t, d);
            }
            else{
                get_block(copy_row1, A, n, m, 0, d);
                multiply(B, A, C, m, m, m, l);
                put_block(copy_row1, C, n, m, 0, d);
            }
        }
        for(int u = 0; u < d; u++){
            if(t%p == k){
                get_block(x, A, n, m, t, u);
                multiply(B, A, C, m, m, m, m);
                put_block(x, C, n, m, t, u);
            }
            else{
                get_block(copy_row2, A, n, m, 0, u);
                multiply(B, A, C, m, m, m, m);
                put_block(copy_row2, C, n, m, 0, u);
            }
        }

        if(l > 0){
            if(t%p == k){
                get_block(x, A, n, m, t, d);
                multiply(B, A, C, m, m, m, l);
                put_block(x, C, n, m, t, d);
            }
            else{
                get_block(copy_row2, A, n, m, 0, d);
                multiply(B, A, C, m, m, m, l);
                put_block(copy_row2, C, n, m, 0, d);
            }
        }

        
        reduce_sum<int>(p);
        for(int q = k; q < h; q+=p){
            if(q == t) continue;
            get_block(a, A, n, m, q, t);

            int mult_rows = q < d ? m : l;

            for(int v = t+1; v < h; v++){
                if(t%p == k) get_block(a, B, n, m, t, v);
                else  get_block(copy_row1, B, n, m, 0, v);
                int mult_cols = v < d ? m : l;
                multiply(A, B, C, mult_rows, m, m, mult_cols);
                get_block(a, B, n, m, q, v);
                for (int I = 0; I < mult_rows; ++I) {
                    for (int J = 0; J < mult_cols; ++J) {
                        B[mult_cols*I + J] -= C[mult_cols*I + J];           
                    }
                }
                put_block(a, B, n, m, q, v);
            }
            for(int v = 0; v < h; v++){
                if(t%p == k) get_block(x, B, n, m, t, v);
                else  get_block(copy_row2, B, n, m, 0, v);
                int mult_cols = v < d ? m : l;
                multiply(A, B, C, mult_rows, m, m, mult_cols);
                get_block(x, B, n, m, q, v);
                for (int I = 0; I < mult_rows; ++I) {
                    for (int J = 0; J < mult_cols; ++J) {
                        B[mult_cols*I + J] -= C[mult_cols*I + J];           
                    }
                }
                put_block(x, B, n, m, q, v);
            }
        }
        reduce_sum<int>(p);
    }
    if(l != 0){
        int res = 0;
        if(d%p == k){
            get_block(a, A, n, m, d, d);
            res = jordan_reverse(A, B, l, norm);
            if (res == 0) {
                for(int i = 0; i < d; i++){
                  get_block(x, A, n, m, d, i);
                  multiply(B, A, C, l, l, l, m);
                  put_block(x, C, n, m, d, i);
                }
                get_block(x, A, n, m, d, d);
                multiply(B, A, C, l, l, l, l);
                put_block(x, C, n, m, d, d);                
            }
        }

        reduce_sum(p, &res, 1);
        if (res < 0) { 
            return -1;
        }
        if(d%p != k){
            copy_matrix(copy_row2, x+d*m*n, l, n);
        }
        reduce_sum<int>(p);
        for (int q = k; q < d; q += p) { 
            get_block(a, A, n, m, q, d);
            for(int v = 0; v < h; v++){
                if(d%p == k) get_block(x, B, n, m, d, v);
                else  get_block(copy_row2, B, n, m, 0, v);
                int mult_cols = v < d ? m : l;
                multiply(A, B, C, m, l, l, mult_cols);
                get_block(x, B, n, m, q, v);
                for (int I = 0; I < m; ++I) {
                    for (int J = 0; J < mult_cols; ++J) {
                        B[mult_cols*I + J] -= C[mult_cols*I + J];           
                    }
                }
                put_block(x, B, n, m, q, v);
            }
        }
    }
    reduce_sum<int>(p);

    return 0;

}














void* thread_func(void* ptr){
    Args* arg = (Args*)ptr;
    double* a = arg->a;
    double* x = arg->x;

    int n = arg->n;
    int m = arg->m;
    int p = arg->p;
    int k = arg->k;
    int r = arg->r;
    int s = arg->s;
    const char* file_name = arg->filename;

    cpu_set_t cpu;
    CPU_ZERO(&cpu);
    int n_cpus = get_nprocs();
    int cpu_id = n_cpus-1-(k%n_cpus);
    CPU_SET(cpu_id, &cpu);
    pthread_t tid = pthread_self();
    pthread_setaffinity_np(tid, sizeof(cpu), &cpu);

    for (int i = k*m; i < n; i += p*m) {
        int h = (i + m < n) ? m : n - i;
        memset(a + i*h, 0, n*h*sizeof(double));
        memset(x + i*h, 0, n*h*sizeof(double));
    }

    reduce_sum<int>(p);

    double* A = new double[m*m];
    double* B = new double[m*m];
    double* C = new double[m*m];
    double* copy_row1 = new double[m*n];
    double* copy_row2 = new double[m*n];

    if(file_name != nullptr){
        int res = 0;
        if(k == 0){
          res = read(a, n, file_name);
        }
        reduce_sum(p, &res, 1);
        if(res < 0){
            arg->status = res;
            delete[] A;
            delete[] B;
            delete[] C;
            delete[] copy_row1;
            delete[] copy_row2;
            return nullptr;
        }
    }
    else init_a_thread(a, n, m, s, k, p);
    init_x_thread(x, n, m, k, p);

    double norm = 0;

    if(k == 0){
        print(a, n, n, r);

        norm = matrix_norm(a, n);
    }

    reduce_sum(p, &norm, 1);

    double cpu_time = get_cpu_time();
    int res = solve(a, x, n, m, p, k, norm, A, B, C, copy_row1, copy_row2);
    cpu_time = get_cpu_time() - cpu_time;

    reduce_sum(p,&res,1);

    arg->status = res;
    arg->cpu_time = cpu_time;
    if (res < 0) {
      delete[] A;
      delete[] B;
      delete[] C;
      delete[] copy_row1;
      delete[] copy_row2;
      return nullptr;
    }
    if(k == 0){
      print(x, n, n, r);
    }
    delete[] A;
    delete[] B;
    delete[] C;
    delete[] copy_row1;
    delete[] copy_row2;
    reduce_sum<int>(p);
    return nullptr;

}