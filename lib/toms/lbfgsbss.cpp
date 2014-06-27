/*
Copyright (C) 1998, 1999, 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                       http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#include <math.h>
#include <time.h>
#include "toms.h"
#include "f2c.h"
#include "f2clib.h"

double lbfgsb_SS(const VOF &fnc, long n, long np, double *x, double
                 *l, double *u, long *nbd, const SetLbfgsb *set)
{
  long m;
  long iprint;
  int iter;
  double factr;
  double pgtol;
  double small;
  double step;
  if (set)
  {
    m = set->m;
    iprint = set->iprint;
    iter = set->iter;
    factr = set->factr;
    pgtol = set->pgtol;
    small = set->small;
    step = set->step;
  }
  else
  {
    static SetLbfgsb set;
    m = set.m;
    iprint = set.iprint;
    iter = set.iter;
    factr = set.factr;
    pgtol = set.pgtol;
    small = set.small;
    step = set.step;
  }
  char task[60];
  char csave[60];
  long lsave[4];
  long isave[44];
  double dsave[29];
  double f;
  vec_double g(n);
  vec_double fn(np);
  vec_double fn1(np);
  vec_double fn2(np);
  vec_double wa((2*m + 4)*n + 11*m*m + 8*m);
  vector<long> iwa(3*n); 
  
  bool central = false;
  size_t i, j;
  double oldx, dx;
  s_copy(task, "START", 60L, 5L);
  while (1)
  {
    setulb_(&n, &m, x, l, u, nbd, &f, &*g.begin(), &factr, &pgtol, &*wa.begin(), 
        &*iwa.begin(), task, &iprint, csave, lsave, isave, dsave);
    if (s_cmp(task, "FG", 2L, 2L) == 0)
    {
      fnc(&n, x, &*fn.begin(), &np);
      f = 0.;
      for (i = 0; i < np; ++i)
      {
        f += fn[i]*fn[i];
      }
      for (j = 0; j < n; ++j)
      {
        oldx = x[j];
        dx = step*(x[j] + step);
        x[j] += dx;
        fnc(&n, x, &*fn1.begin(), &np);
        g[j] = 0.;
        if (central)
        {
          x[j] = oldx - dx;
          fnc(&n, x, &*fn2.begin(), &np);
          for (i = 0; i < np; ++i)
          {
            g[j] += fn[i]*(fn1[i] - fn2[i])/dx/2.0;
          }
        }
        else
        {
          for (i = 0; i < np; ++i)
          {
            g[j] += fn[i]*(fn1[i] - fn[i])/dx;
          }
        }
        g[j] *= 2.;
        x[j] = oldx;
        if (!central && fabs(g[j]) < small)
          central = true;
      }
    }
    else if (s_cmp(task, "NEW_X", 5L, 5L) == 0)
    {
      if (isave[34] > iter)
      s_copy(task, "STOP - iter", 60L, 11L);
      continue;
    }
    else
      break;
  }
  return f;
}

#ifndef CLK_TCK
#include <unistd.h>
long CLK_TCK= sysconf(_SC_CLK_TCK);
#endif

extern "C" 
{
  double dpmeps_(void)
  {
     return DBL_EPSILON;
  }
  int timer_(double *ttime)
  {
      *ttime = (double)clock()/double(CLK_TCK);
      return 0;
  } 
}
