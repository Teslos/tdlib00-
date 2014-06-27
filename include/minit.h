#ifndef __MINIT_H
#define __MINIT_H

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

/************************ CONSTANT DEFINITIONS ********************************/
#include <general.h>
#include <fstream>

class Cmatrix : public vec_double
{
  size_t m;
  size_t n;

public:
  Cmatrix() : m(0), n(0) {}
  Cmatrix(size_t m_, size_t n_) : vec_double(m_*n_), m(m_), n(n_) {}

  void resize(size_t m_, size_t n_)
    {m = m_; n = n_; vec_double::resize(m_*n_);}
  void reserve(size_t m_, size_t n_)
    {vec_double::reserve(m_*n_);}
  void clear()
    {m = n = 0; vec_double::clear();}

  size_t NRows() const {return m;}
  size_t NCols() const {return n;}
  
  double* operator[](size_t i)
  {
    return &*begin() + i*n;
  }
  const double* operator[](size_t i) const
  {
    return &*begin() + i*n;
  }

  double& operator()(size_t i, size_t j)
  {
    return *(begin() + i*n + j);
  }
  const double& operator()(size_t i, size_t j) const
  {
    return *(begin() + i*n + j);
  }

};

class minit
{
  int k, m, n, p, lcol, L, im, jmin, jm, imax;
  static double td;
  double gmin, phimax;
/************************** FILE-LIMITED ARRAYS *******************************/
  vec_int imin, jmax, ind, ind1, chk;
  Cmatrix e;
  vec_double thmin, delmax;

/************************ minit_space_handler **********************************
 * This routine performs space maintenance and initializes global arrays.     
 * It looks at n, m and p, which are static in this file, to look for
 * previously allocated space.
 ******************************************************************************/
  void minit_space_handler(size_t N, size_t M, size_t P);
/********************************* MINIT ***************************************
 * This is the main procedure. It is invoked with the various matrices,
 * after space for other arrays has been generated elsewhere (in zx3lp.)
 * It returns
 * 0 if everything went OK and x, w, z as the solutions.
 * 1 if no solution existed.
 * 2 if primal had no feasible solutions.
 * 3 if primal had unbounded solutions.
 ******************************************************************************/
  int minit_solve(Cmatrix &A, double* b, double* c, double* x, 
      double* w, double &z);
  
/****************************** rowtrans ***************************************
 * Performs the usual tableau transformations in a linear programming
 * problem, (p,q) being the pivotal element.
 * Returns the following error codes:
 * 0: Everything was OK.
 * 1: No solution.
 ******************************************************************************/
  int rowtrans(int p, int q);

/****************************** progamma ***************************************
 * Performs calculations over columns to determine the pivot element.
 ******************************************************************************/
  int progamma();
/****************************** prophi *****************************************
 * Performs calculations over rows to determine the pivot element.
 ******************************************************************************/
  int prophi();
/****************************** phase1 *****************************************
 * Applied only to equality constraints if any.
 ******************************************************************************/
  int phase1();

/****************************** tab *****************************************
 * The following procedure is for debugging. It simply prints the
 * current tableau.
 ******************************************************************************/
 void tab();

public:
  ofstream* debug;
  
  minit() : debug(0) {}
    
/********************************* ZX3LP ***************************************
 * This is a user-specific front end to minit.
 * In its present form, it has been written for compatibility with ifplan. 
 * Note that zx3lp was originally a fortran subroutine. The weird return
 * values and flags (such as IER) are inherited from there.
 * Variable IA and arrays RW, IW are not necessary.
 ******************************************************************************/
  void zx3lp(Cmatrix &A, vec_double &B, vec_double &C, size_t N, size_t M1, 
      size_t M2, double &S, vec_double &PSOL, vec_double &DSOL, int &IER)
  {
    zx3lp(A, &*B.begin(), &*C.begin(), N, M1, M2, S, &*PSOL.begin(), 
        &*DSOL.begin(), IER);
  }
  
  void zx3lp(Cmatrix &A, double* B, double* C, size_t N, size_t M1, 
      size_t M2, double &S, double* PSOL, double* DSOL, int &IER)
  {
    int err;
    minit_space_handler(N, M1+M2, M2);
    err = minit_solve(A, B, C, PSOL, DSOL, S);
    switch(err)
    {
      case 1:
        IER = 134;
        break;
      case 2:
        IER = 133;
        break;
      case 3:
        IER = 131;
        break;
      default:
        IER = 0;
    }
  }
};

#endif

