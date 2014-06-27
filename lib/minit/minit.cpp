/****************************** MINIT ******************************************
 * This file contains procedures to solve a Linear Programming
 * problem of n variables and m constraints, the last p of which
 * are equality constraints by the dual simplex method.
 *
 * The code was originally in Algol 60 and was published in Collected
 * Algorithms from CACM (alg #333) by Rudolfo C. Salazar and Subrata
 * K. Sen in 1968 under the title MINIT (MINimum ITerations) algorithm
 * for Linear Programming.
 * It was directly translated into C by Badri Lokanathan, Dept. of EE,
 * University of Rochester, with no modification to code structure. 
 *
 * The problem statement is
 * Maximise z = cX
 *
 * subject to
 * AX <= b
 * X >=0
 *
 * c is a (1*n) row vector
 * A is a (m*n) matrix
 * b is a (m*1) column vector.
 * e is a matrix with with (m+1) rows and lcol columns (lcol = m+n-p+1)
 *   and forms the initial tableau of the algorithm.
 * td is read into the procedure and should be a very small number,
 *   say 10**-8
 * 
 * The condition of optimality is the non-negativity of
 * e[1,j] for j = 1 ... lcol-1 and of e[1,lcol] = 2 ... m+1.
 * If the e[i,j] values are greater than or equal to -td they are
 * considered to be non-negative. The value of td should reflect the 
 * magnitude of the co-efficient matrix.
 *
 ******************************************************************************/

/*
Changes to minit by Evgenii Rudnyi:

  1) the program is made reentrent - converted to class
  2) float converted to double
  3) the internal format of matirices A and E is changed
  4) the memory allocation is changed to C++
  
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                       http://www.chem.msu.su/~rudnyi/

The modifications are a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#include <float.h>
#include <iomanip>
#include "minit.h"

/************************ CONSTANT DEFINITIONS ********************************/
#define LARGENUMBER 1e6
#define VLARGENUMBER 1e8
#define VSMALLNUMBER 1e-8
//#define LARGENUMBER 1./(DBL_EPSILON*10)
//#define VLARGENUMBER 1./(DBL_EPSILON*10000)
//#define VSMALLNUMBER DBL_EPSILON*10000

double minit::td=VSMALLNUMBER;   /* This is application specific. */

void minit::minit_space_handler(size_t N, size_t M, size_t P)
{
  /* Initialize globals. */
  m = M;
  n = N;
  p = P;
  lcol = m + n - p + 1;
  int MS = m + 1;
  int LS = m + n - p + 1;
  
  jmax.resize(MS);
  ind1.resize(MS);
  chk.resize(MS);
  delmax.resize(MS);
  imin.resize(LS);
  ind.resize(LS);
  thmin.resize(LS);
  e.resize(MS, LS);

  /* Finally, initialize all arrays to 0. */
  for(int i = 0; i <= m; i++)
  {
    jmax[i] = ind1[i] = chk[i] = 0;
    delmax[i] = 0.0;
    for(int j = 0; j < lcol; j++)
      e[i][j] = 0.0;
  }

  for(int j = 0; j < lcol; j++)
  {
    imin[j] = ind[j] = 0;
    thmin[j] = 0.0;
  }
}

int minit::minit_solve(Cmatrix &A, double *b, double *c, double *x, 
    double *w, double &z)
{
  register int i, j;

  /* Generate initial tableau. */
  for(j = 0; j < n; j++)
    e[0][j] = -c[j];

  for(i = 0; i < m; i++)
  {
    for(j = 0; j < n; j++)
      e[i+1][j] = A[i][j];
  }

  for(j = 0; j < m; j++)
    e[j+1][lcol-1] = b[j];

  for(i = 0; i < m - p; i++)
    e[1+i][n+i] = 1.0;

  /* Now begins the actual algorithm. */
  for(i = 1; i < m+1; i++)
    chk[i] = -1;  /* Indexing problem; Algol
         * version treates 0 as out
         * of bounds, in C we prefer -1.
         */
  if(!p)  /* There are no equality constraints */
    goto RCS;
  else
  {
  if(phase1())
    return(1);
  }
    
  RCS: L = 1; k = 1;

  if(debug)
    tab();  /* Print the tableau. */

  for(j = 0; j < lcol - 1; j++)
  {
    if(e[0][j] < -td)
      ind[(L++)-1] = j; /* ind[L] keeps track of the
             * columns in which e[0][j]
             * is negative.
             */
  }

  for(i = 1; i < m + 1; i++)
  {
    if(e[i][lcol-1] < -td)
      ind1[(k++)-1] = i; /* ind1[k] keeps track of the
              * rows in which e[i][lcol]
              * is negative.
              */
  }

  if(L == 1)
  {
    if(k == 1) /* Results */
      goto RESULTS;
    else
    {
      if(k == 2)
      {
        for(j = 0; j < lcol - 1; j++)
        {
          if(e[ind1[0]][j] < 0)
            goto R;
        }
        if (debug)
          *debug << "the bad row is " << ind1[0] << endl;
        /* Primal problem has no feasible solutions. */
        return(2);
      }
      else
        goto R;
    }
  }
  else
  {
    if(L == 2)
    {
      if(k == 1)
      {
        for(i = 1; i < m + 1; i++)
        {
          if(e[i][ind[0]] > 0)
            goto C;
        }

        /* Primal problem is unbounded. */
        return(3);
      }
      else
        goto S;
    }

    if(k == 1)
      goto C;
    else
      goto S;
  }

  R:
  prophi();
  if(rowtrans(imax,jm))
    return(1);
  goto RCS;

  C:
  progamma();
  if(rowtrans(im,jmin))
    return(1);
  goto RCS;

  S:
  progamma();
  prophi();
  if(gmin == LARGENUMBER)
  {
    if(rowtrans(imax,jm))
      return(1);
    goto RCS;
  }

  if(phimax == - LARGENUMBER)
  {
    if(rowtrans(im,jmin))
      return(1);
    goto RCS;
  }

  if(fabs(phimax) > fabs(gmin))
  {
    if(rowtrans(imax,jm))
      return(1);
  }
  else
  {
    if(rowtrans(im,jmin))
      return(1);
  }
  goto RCS;

  RESULTS:
  /* Output results here. */
  z = e[0][lcol-1];
  for(i = 0; i < n; i++)
    x[i] = 0.0;

  for(j = 0; j < m; j++)
    w[j] = 0.0;

  for(i = 1; i < m + 1; i++)
  {
    if(chk[i] >= n)
      chk[i] = -1;

    if(chk[i] > -1)
      x[chk[i]] = e[i][lcol-1];
  }

  for(j = n; j < lcol - 1; j++)
    w[j-n] = e[0][j];
  
  return(0);
}

int minit::rowtrans(int p, int q)
{
  register int i, j;
  double dummy;

  if(p == -1 || q == -1) /* No solution. */
    return(1);

  dummy = e[p][q];

  if(debug)
    *debug << "--\nPivot element is e[" << p << "][" << q << "] = " 
      << dummy << endl;

  for(j = 0; j < lcol; j++)
    e[p][j] /= dummy; 

  for(i = 0; i < m + 1; i++)
  {
    if(i != p)
    {
      if(e[i][q] != 0.0)
      {
        dummy = e[i][q];
        for(j = 0; j < lcol; j++)
          e[i][j] -= e[p][j] * dummy;
      }
    }
  }

  chk[p] = q;
  return(0);
} /* rowtrans */

int minit::progamma()
{
  double theta, gamma;
  register int i, L1;

  gmin = LARGENUMBER; /* gmin is set to a large no. for
         * initialization purposes.
         */
  jmin = -1;    /* Array indices in C start from 0 */

  for(L1 = 0; L1 < L - 1; L1++)
  {
    imin[ind[L1]] = 0;
    thmin[ind[L1]] = LARGENUMBER;
    for(i = 1; i < m + 1; i++)
    {
      if(e[i][ind[L1]] > td && e[i][lcol-1] >= -td)
      {
        theta = e[i][lcol-1] / e[i][ind[L1]];
        if(theta < thmin[ind[L1]])
        {
          thmin[ind[L1]] = theta;
          imin[ind[L1]] = i;
        }
      }
    }

    if(thmin[ind[L1]] == LARGENUMBER)
      gamma = VLARGENUMBER;
    else
      gamma = thmin[ind[L1]] * e[0][ind[L1]];

    if(gamma < gmin)
    {
      gmin = gamma;
      jmin = ind[L1];
    }
  }
  if(jmin > -1)
    im = imin[jmin];
} /* progamma */

int minit::prophi()
{
  double delta, phi;
  register int j, k1;

  phimax = - LARGENUMBER; /* phimax is set to a small no. for
         * initialization purposes.
         */
  imax = -1;    /* Array indices in C start from 0 */
  
  for(k1 = 0; k1 < k - 1; k1++)
  {
    jmax[ind1[k1]] = 0;
    delmax[ind1[k1]] = - LARGENUMBER;
    for(j = 0; j < lcol - 1; j++)
    {
      if(e[ind1[k1]][j] < -td && e[0][j] >= -td)
      {
        delta = e[0][j] / e[ind1[k1]][j];
        if(delta > delmax[ind1[k1]])
        {
          delmax[ind1[k1]] = delta;
          jmax[ind1[k1]] = j;
        }
      }
    }

    if(delmax[ind1[k1]] == - LARGENUMBER)
      phi = - VLARGENUMBER;
    else
      phi = delmax[ind1[k1]] * e[ind1[k1]][lcol-1];
    
    if(phi > phimax)
    {
      phimax = phi;
      imax = ind1[k1];
    }
  }
  if(imax > -1)
    jm = jmax[imax];
} /* prophi */

int minit::phase1()
{
  double theta, gamma;
  register int i, j, r;
  /* Fix suggested by Holmgren, Obradovic, Kolm. */
  register int im1, jmin1, first;

  im1 = jmin1 = -1;
  /* Fix suggested by Messham to allow negative coeffs. in
   * equality constraints.
   */
  for(i = m - p + 1; i < m + 1; i++)
  {
    if(e[i][lcol - 1] < 0)
    {
      for(j = 0; j < lcol; j++)
        e[i][j] = -e[i][j];
    }
  }

  for(r = 0; r < p; r++)
  {
    gmin = LARGENUMBER; /* gmin is set to a large no. for
           * initialization purposes.
           */
    L = 1;
    jmin = -1;
    first = 1;
    for(j = 0; j < n; j++)
    {
      thmin[j] = LARGENUMBER;
      /* Fix suggested by Kolm and Dahlstrand */
      /* if(e[0,j] < 0) */
      if(e[0][j] < -td)
        ind[(L++)-1] = j;
    }
  L1: if(L == 1)
    {
      for(j = 0; j < n; j++)
        ind[j] = j;

      L = n + 1;
    }

    for(k = 0; k < L - 1; k++)
    {
      for(i = m - p + 1; i < m + 1; i++)
        if(chk[i] == -1)
        {
          /* Fix suggested by Kolm
           * and Dahlstrand
           */
          /* if(e[i][ind[k]] > 0.0) */
          if(e[i][ind[k]] > td)
          {
            if((theta = e[i][lcol-1] /
              e[i][ind[k]]) < thmin[ind[k]])
            {
              thmin[ind[k]] = theta;
              imin[ind[k]] = i;
            }
          }
        }

      /* Fix suggested by Obradovic overrides
       * fixes suggested by Kolm and Dahstrand
       * as well as Messham.
       */
      if(thmin[ind[k]] < LARGENUMBER)
      {
        gamma = thmin[ind[k]] * e[0][ind[k]];
        if(gamma < gmin)
        {
          gmin = gamma;
          jmin = ind[k];
        }
      }
    }
    if(jmin == -1)
    {
      if(first)
      {
        first = 0;
        L = 1;
        goto L1;
      }
      else
        im = -1;
    }
    else
      im = imin[jmin];
    
    if(im == im1 && jmin == jmin1)
    {
      L = 1;
      goto L1;
    }

    if(debug)
    {
      tab();  /* Print the tableau. */
      *debug << "im is " << im << " and jmin is " << jmin << endl;
    }

    if(rowtrans(im,jmin))
      return(1);
    
    im1 = im;
    jmin1 = jmin;
  }
  return(0);
} /* phase1 */

void minit::tab()
{
  if (!debug)
    return;
  register int i, j;

  *debug << endl << "-------------------------- TABLEAU ----------------------" 
    << endl;

  for(i = 0; i < m+1; i++)
  {
    for(j = 0; j < lcol; j++)
      *debug << setw(12) << e[i][j];
    *debug << endl;
  }
}
