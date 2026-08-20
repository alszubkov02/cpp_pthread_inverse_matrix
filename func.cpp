#include "func.h"


double get_full_time()
{
    struct timeval buf;
    gettimeofday(&buf, 0);
    return buf.tv_sec + buf.tv_usec / 1.e6;
}

double get_cpu_time()
{
    struct rusage buf;
    getrusage(RUSAGE_THREAD, &buf);
    return buf.ru_utime.tv_sec + buf.ru_utime.tv_usec /1.e6;
}

double func(int s, int n, int i, int j)
{
  if(s == 1)
    {
      return n - std::max (i,j);
    }
  else if(s == 2)
    {
      return std::max (i,j)+1;
    }
  else if(s == 3)
    {
      return abs (i-j);
    }
  else
    {
      return 1.0 / ((i + j + 1)*1.0);
    }
}

void init_a (double* a, int n, int s)
{
  for(int i = 0; i < n; i++)
    {
      for(int j = 0; j < n; j++)
        {
          a[i*n + j] = func (s, n, i, j);
        }
    }
}

int read (double* a, int n, const char* name)
{
    FILE* f;
    f = fopen (name, "r");
    if(!f)
      {
        printf ("Can not read file\n");
        return -1;
      }
    for(int i = 0; i < n; i++)
      {
        for(int j = 0; j < n; j++)
          {
            if(fscanf (f, "%lf", &a[i*n+j]) != 1)
              {
                fclose (f);
                return -1;
            }
        }
    }
    fclose (f);
    return 1;
}

void print (double* a, int n, int m, int r)
{
  int n_ = (n < r) ? n : r;
  int m_ = (m < r) ? m : r;

  for(int i = 0; i < n_; i++)
    {
      for(int j = 0; j < m_; j++)
        {
          printf ("%10.3e ", a[i*n+j]);
        }
      printf ("\n");
    }
  printf ("\n");
}


double matrix_norm (double* a, int n)
{
  double norm = 0.0, s;
  for(int j = 0; j < n; j++)
    {
      s = 0.0;
      for(int i = 0; i < n; i++)
        {
          s += fabs (a[i*n+j]);
        }
      norm = std::max (s, norm);
    }
  return norm;
}

double norm1 (double* a, double* x, double* id, int n)
{
	if(n > 11000)
    {
      return 0.0;
    }
	multiply (a, x, id, n, n, n, n);
	for(int i = 0; i < n; i++)
    {
		  id[i*n+i] -= 1.0;
	  }

	double norm = matrix_norm (id,n);
	return norm;
}

double norm2 (double* a, double* x, double* id, int n)
{
	if(n > 11000)
    {
      return 0.0;
    }
	multiply (x, a, id, n, n, n, n);
	for(int i = 0; i < n; i++)
    {
		  id[i*n+i] -= 1.0;
	  }
	double norm = matrix_norm (id,n);
	return norm;
}

void multiply (double* a, double* b,double* c, int n1, int m1, int n2, int m2)
{
  for(int i = 0; i < n1; i++){
    for(int j = 0; j < m2; j++){
      c[i*m2+j] = 0.0;
    }
  }
  int v3 = n1%3;
  int h3 = m2%3;
  int i = 0, j = 0, k = 0;
  double s, tmp;
  double s00, s01, s02,
          s10, s11, s12,
          s20, s21, s22;
  for(i = 0; i < v3; i++)
    {
      for(j = 0; j < h3; j++)
        {
          s = 0;
          for(k = 0; k < n2; k++)
            {
              s += a[i*m1+k]*b[k*m2+j];
            }
          c[i*m2+j] = s;
        }
      for(; j < m2; j += 3)
        {
          s00 = 0;
          s01 = 0;
          s02 = 0;
          for(k = 0; k < m1; k++)
            {
              tmp = a[i*m1+k];
              s00 += tmp*b[k*m2+j];
              s01 += tmp*b[k*m2+j+1];
              s02 += tmp*b[k*m2+j+2];
            }
          c[i*m2+j] += s00;
          c[i*m2+j+1] += s01;
          c[i*m2+j+2] += s02;
        }
    }
  for(; i < n1; i += 3)
    {
      for(j = 0; j < h3; j++)
        {
          s00 = 0;
          s10 = 0;
          s20 = 0;
          for(k = 0; k < m1; k++)
            {
              tmp = b[k * m2 + j];
              s00 += a[i * m1 + k]*tmp;
              s10 += a[(i+1) * m1 + k]*tmp;
              s20 += a[(i+2) * m1 + k]*tmp;
            }
          c[i * m2 + j] += s00;
          c[(i+1) * m2 + j] += s10;
          c[(i+2) * m2 + j] += s20;
        }
      for(; j < m2; j += 3)
        {
          s00 = 0;
          s01 = 0;
          s02 = 0;
          s10 = 0;
          s11 = 0;
          s12 = 0;
          s20 = 0;
          s21 = 0;
          s22 = 0;
          for(k = 0; k < m1; k++)
            {
              s00 += a[i*m1+k]*b[k*m2+j];
              s01 += a[i*m1+k]*b[k*m2+j+1];
              s02 += a[i*m1+k]*b[k*m2+j+2];
              s10 += a[(i+1)*m1+k]*b[k*m2+j];
              s11 += a[(i+1)*m1+k]*b[k*m2+j+1];
              s12 += a[(i+1)*m1+k]*b[k*m2+j+2];
              s20 += a[(i+2)*m1+k]*b[k*m2+j];
              s21 += a[(i+2)*m1+k]*b[k*m2+j+1];
              s22 += a[(i+2)*m1+k]*b[k*m2+j+2];
            }
          c[i*m2+j] += s00;
          c[i*m2+j+1] += s01;
          c[i*m2+j+2] += s02;
          c[(i+1)*m2+j] += s10;
          c[(i+1)*m2+j+1] += s11;
          c[(i+1)*m2+j+2] += s12;
          c[(i+2)*m2+j] += s20;
          c[(i+2)*m2+j+1] += s21;
          c[(i+2)*m2+j+2] += s22;
        }
    }
}