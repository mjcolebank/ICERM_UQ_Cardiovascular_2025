/***************************************************************************/
/*                                                                         */
/*  Program: tools.C                                                       */
/*  Version: 2.0                                                           */
/*  By: Mette Olufsen, Math-Tech                                           */
/*  Date: 14. Jan. 1997                                                    */
/*                                                                         */
/*  This module includes auxiliary functions such as error-handling,       */
/*                      a version of the Numerical Recipes in C routine    */
/*  for solving a nonlinear set of equations using Newton Raphson's met-   */
/*  hod, and finally one for solving a one-dimensional nonlinear equation. */
/*  Since these functions are taken directly from Numerical Recipes they   */
/*  will not be commented further.                                         */
/*                                                                         */
/*  The program tools uses a number of vector and matrix handling functions*/
/*  that are available in the library nrutil, therefore I have included    */
/*  nrutil.h. Apart from that only standard c-libraries stdio.h, stdlib.h, */
/*  and math.h are used.                                                   */
/*                                                                         */
/***************************************************************************/


#include "tools.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>

using namespace std;

// Standard parameters.
const double EPS = 1.0e-20;

// The general error handling of this program prints out an error message
// and terminates the program.
void error (const char *s, const char *u)
{
  fprintf (stdout, "ERROR (%s): %s\n", s, u);
  //exit (0);
}

// The remaining functions are from the Numerical Recipes in C library.
// They are all used with Newton Raphson's method when solving a set of
// non-linear equations, and when solving a one-dimensional nonlinear equation,
// also using Newton Raphson's method.

void lubksb (double *a[], int n, int indx[], double b[])
{
  int ii = 0;

  for (int i=0; i<n; i++)
  {
    int ip = indx[i];
    double sum = b[ip-1];
    b[ip-1] = b[i];
    if (ii != 0) {
      for (int j=ii-1; j<i; j++)
      {
        sum = sum - a[i][j]*b[j];
      }
    }
    else if (sum) ii=i+1;
    b[i] = sum;
  }
  for (int i=n-1; i>=0; i--)
  {
    double sum = b[i];
    for (int j=i+1; j<n; j++)
    {
      sum = sum - a[i][j]*b[j];
    }
    b[i] = sum / a[i][i];
  }
}


void ludcmp (double *a[], int n, int indx[], double *d)
{
  int imax = -1;
  double vv[n];
  *d   = 1.0;
  for (int i=0; i<n; i++)
  {
    double big = 0.0;
    for (int j=0; j<n; j++)
    {
	big = max(big, fabs(a[i][j]));
    }
    if (big == 0.0) {
      const char *error_text = "Singular matrix in routine LUDCMP";
      fprintf (stdout, "Numerical Recipes run-time error...\n");
      fprintf (stdout, "%s\n",error_text);
      fprintf (stdout, "...now exiting to system...\n");
      exit (1);
    }
    vv[i] = 1.0/big;
  }

  for (int j=0; j<n; j++)
  {
    for (int i=0; i<j; i++)
    {
      double sum = a[i][j];
      for (int k=0; k<i; k++)
      {
        sum = sum - a[i][k]*a[k][j];
      }
      a[i][j] = sum;
    }

    double big = 0.0;
    for (int i=j; i<n; i++)
    {
      double sum = a[i][j];
      for (int k=0; k<j; k++)
      {
	sum = sum - a[i][k]*a[k][j];
      }
      a[i][j] = sum;

      double dum = vv[i]*fabs(sum);
      if (dum >= big) {
	big  = dum;
	imax = i;
      }
    }

    if (j != imax) {
      for (int k=0; k<n; k++)
      {
	swap(a[imax][k], a[j][k]);
      }
      *d       = -(*d);
      vv[imax] = vv[j];
    }
    indx[j] = imax+1;
    if (a[j][j] == 0.0) a[j][j] = EPS;
    if (j+1 != n) {
      double dum = 1.0 / (a[j][j]);
      for (int i=j+1; i<n; i++)
      {
        a[i][j] = a[i][j]*dum;
      }
    }
  }
}


// Takes one step with Newton Raphson's method. Assumes to find zeros
// for a multi-dimensional problem.
int zero (double* x, int n, double tolx, double tolf,
          double fvec[], double *fjac[])
{
  int indx[n];
  double p[n];
  double errf = 0.0;
  for (int i=0; i<n; i++)
  {
    errf = errf + fabs(fvec[i]);
  }
  if (errf <= tolf)
  {
    return (1);
  }

  for (int i=0; i<n; i++) {
    p[i] = - fvec[i];
  }

  double d; // unused return value
  ludcmp (fjac, n, indx, &d);
  lubksb (fjac, n, indx, p);
  double errx = 0.0;
  for (int i=0; i<n; i++)
  {
    errx = errx + fabs(p[i]);
    x[i] = x[i] + p[i];
  }
  if (errx <= tolx)
  {
    return (1);
  }
  return (2);
}


// Takes one step with Newton Raphson's method. Assumes to find zeros for
// a one-dimensional problem.
bool zero_1d (double *x, double f, double df, double tolx)
{
  double dx  = f/df;
  *x  = *x-dx;
  //if (fabs(dx) < tolx) return(1); else   // Original statement.
  return (fabs(dx) < tolx && fabs(f) < tolx);
}

