/*
Copyright (C) 1998, 1999, 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                       http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#ifndef __TOMS_H
#define __TOMS_H

#include <string>
#include <callback.hpp>
#include <float.h>
#include <general.h>

typedef CBFunctor1wRet<double, double> D2D;
typedef CBFunctor2wRet<double*, long, double> DL2D;
typedef CBFunctor4<long*, double*, double*, long*> VOF;

// from lapack

extern "C"
{
  long idamax_(long *n, double *x, long *incx);
  int dcopy_(long *n, double *x, long *incx, double *y, long *incy);
  double dnrm2_(long *n, double *x, long *incx);
  int daxpy_(long *n, double *da, double *dx, long *incx, double *dy, 
      long *incy);
  double ddot_(long *n, double *dx, long *incx, double *dy, long *incy);
  int dscal_(long *n, double *da, double *dx, long *incx);
}

// from enclofx.cpp

class err_rroot : public gError
{
public:
  err_rroot() : gError("rroot: No change in sign of f(x) was found") {}
};

int rroot_(const D2D &f, long *neps, double *eps, double *a, double *b,
           double *root);

inline double rroot(const D2D &f, double a, double b, double eps = 0.)
{
  double root;
  long neps = 1000;
  if (!eps)
    eps = DBL_EPSILON;
  rroot_(f, &neps, &eps, &a, &b, &root);
  return root;
}

// end enclofx

// from brent1.cpp

class err_brent1 : public gError
{
public:
  int ier;
  err_brent1(const string &msg, int ier) : gError(msg), ier(ier) {}
};

int brent1_(const VOF &fcn, long *n, double *x, double *fvec, double
            *tol, long *info, double *wa, long *lwa);

inline int brent1(const VOF &fcn, long n, double *x, double *fvec,
                  double tol = 0)
{
  if (!tol)
    tol = DBL_EPSILON*100.;
  long info;
  long lwa = n*(n + 3);
  vec_double wa(lwa);
  brent1_(fcn, &n, x, fvec, &tol, &info, &*wa.begin(), &lwa);
  if (info < 1 || info > 3)
    throw err_brent1("brent1: no convergence", info);
  return info;
}

// end brent1.cpp

// from chabis.cpp

class err_chabis : public gError
{
public:
  int ier;
  err_chabis(const string &msg, int ier) : gError(msg), ier(ier) {}
};

int chabis(const DL2D &fnc, long n, double *x, double *h,
           double *vas, double delta = 0);

// end chabis.cpp

// lbfgsb

struct SetLbfgsb
{
  int m;
  int iprint;
  int iter;
  double factr;
  double pgtol;
  double small;
  double step;
  SetLbfgsb()
  {
    m = 10;
    iprint = -1;
    iter = 100;
    factr = 1e5;
    pgtol = DBL_EPSILON*1e8;
    small = 0.1;
    step = 100.*sqrt(DBL_EPSILON);
  }
};

double lbfgsb_SS(const VOF &fnc, long n, long np, double *x, double
                 *l, double *u, long *nbd, const SetLbfgsb *set = 0);

extern "C" 
{
  int setulb_(long *n, long *m, double *x, double *l, double *u, long
            *nbd, double *f, double *g, double *factr, double *pgtol,
            double *wa, long *iwa, char *task, long *iprint, char
            *csave, long *lsave, long *isave, double *dsave);

  int timer_(double *ttime);
  double dpmeps_(void);
}

// end lbfgsb

extern "C"
{
  void tsdumj_(double *x, double *aha, int *maxn, int *m, int *n); 
  int tsnesv_(int *maxm, // LEADING DIMENSION OF AJA, ANLS, AND FV
      int *maxn, 
      int *maxp, // LEADING DIMENSION OF S AND SHAT
      double *xc, // XC(N) - INITIAL ESTIMATE OF SOLUTION
      int *m, // DIMENSIONS OF PROBLEM
      int *n,
      double *typex, // TYPX(N) - TYPICAL SIZE FOR EACH COMPONENT OF X
      double *typef, // TYPF(M) - TYPICAL SIZE FOR EACH COMPONENT OF F
      int *itnlim, // MAXIMUM NUMBER OF ALLOWABLE ITERATIONS
      
      int *jacflg, // JACFLG = 0 IF ANALYTIC JACOBIAN NOT SUPPLIED
      double *gradtl, // tolerances for gradient, step and function
      double *steptl, 
      double *ftol,
      int *method, // METHOD = 1 : TENSOR METHOD IS USED
      int *global, // GLOBAL = 0 : LINE SEARCH, GLOBAL = 1 : 2-DIMENSIONAL TRUST REGION

      double *stepmx, // MAXIMUM ALLOWABLE STEP SIZE
      double *dlt, // TRUST REGION RADIUS
      int *ipr, // DEVICE TO WHICH TO SEND OUTPUT
      double *x, // X(SQRN) - ESTIMATE TO A ROOT OF FCN ( USED BY UNCMIN)
      double *typxu, // TYPXU(SQRN) - TYPICAL SIZE FOR EACH COMPONENT OF X (USED BY UNCMIN)
      double *xpls, // XPLS(SQRN) - LOCAL MINIMUM OF OPTIMIZATION FUNCTION FCN USED BY UNCMIN
      double *gpls, // GPLS(SQRN) - GRADIENT AT SOLUTION XPLS (USED BY UNCMIN)
      double *a, // A(MAXP,SQRN) - WORKSPACE FOR HESSIAN (OR ESTIMATE) (USED BY UNCMIN)
      double *wrk, // WRK(MAXP,SQRN) - WORKSPACE (USED BY UNCMIN)
      double *dfn, // DFN(M) - DIAGONAL SCALING MATRIX FOR F

      double *wrk1, // WRK1(M) - workspaces
      double *wrk2, // WRK2(M)
      double *wrk3, // WRK3(M)
      double *wrk4, // WRK4(M)
      double *wrk5, //WRK5(M)
      double *fq, // FQ(M)
      double *fqq, // FQQ(M)
      double *fc, // FC(M) - FUNCTION VALUE AT CURRENT ITERATE
      double *fhat, //FHAT(M) - WORKSPACE

      double *anls, // ANLS(MAXM,SQRN) - TENSOR TERM MATRIX
      double *fv, // FV(MAXM,SQRN) - WORKSPACE USED TO STORE PAST FUNCTION VALUES
      double *aja, // AJA(MAXM,N) - JACOBIAN MATRIX
      double *dxn, // DXN(N) - DIAGONAL SCALING MATRIX FOR X
      double *dn, // DN(N) - STANDARD STEP
      double *dt, // DT(N) - TENSOR STEP
      double *df, // DF(N) - workspaces
      double *d, // D(N)
      double *gbar, // GBAR(N)
      double *dbar, // DBAR(N)
      double *dbarp, // DBARP(N)
      
      double *s, // S(MAXN,SQRN) - MATRIX OF PREVIOUS DIRECTIONS
      double *shat, // SHAT(MAXN,SQRN) - MATRIX OF PAST LINEARLY INDEPENDENT DIRECTIONS
      int *curpos,    // CURPOS(N) - pivots
      int *pivot,     // PIVOT(N)
      int *pbar,      // PBAR(N)
      double *epsm, // MACHINE PRECISION
      int *sqrn, // MAXIMUM COLUMN DIMENSION OF ANLS, S, AND SHAT
      void (*fvec)(double *x, double *f, int *m, int *n), 
      
      void (*jac)(double *x, double *aha, int *maxn, int *m, int *n), 
      int *msg, // MESSAGE TO INHIBIT CERTAIN AUTOMATIC CHECKS + OUTPUT
      double *xp, // XP(N)
      double *fp, // FP(M)
      double *gp, // GP(N)
      int *termcd); 
}
/*
       SUBROUTINE TSNESV(MAXM,MAXN,MAXP,XC,M,N,TYPX,TYPF,ITNLIM,
     +                  JACFLG,GRADTL,STEPTL,FTOL,METHOD,GLOBAL,
     +                  STEPMX,DLT,IPR,X,TYPXU,XPLS,GPLS,A,WRK,DFN,
     +                  WRK1,WRK2,WRK3,WRK4,WRK5,FQ,FQQ,FC,FHAT,
     +                  ANLS,FV,AJA,DXN,DN,DT,DF,D,GBAR,DBAR,DBARP,
     +                  S,SHAT,CURPOS,PIVOT,PBAR,EPSM,SQRN,FVEC,
     +                  JAC,MSG,XP,FP,GP,TERMCD)

        INTEGER MAXM,MAXN,MAXP,M,N,SQRN,TERMCD
        INTEGER ITNLIM,JACFLG,METHOD,GLOBAL,MSG,IPR
        DOUBLE PRECISION GRADTL,STEPTL,FTOL,STEPMX,DLT,FPLS,EPSM
        EXTERNAL FVEC,JAC

C**********************************************************************
C THIS IS THE DRIVER FOR NONLINEAR EQUATIONS/NONLINEAR LEAST SQUARES
C PROBLEMS.
C**********************************************************************
C
C       INPUT PARAMETERS :
C       -----------------
C
C       MAXM   : LEADING DIMENSION OF AJA, ANLS, AND FV
C       MAXN   : LEADING DIMENSION OF S AND SHAT
C       XC     : INITIAL ESTIMATE OF SOLUTION
C       M,N    : DIMENSIONS OF PROBLEM
C       TYPX   : TYPICAL SIZE FOR EACH COMPONENT OF X
C       TYPF   : TYPICAL SIZE FOR EACH COMPONENT OF F
C       ITNLIM : MAXIMUM NUMBER OF ALLOWABLE ITERATIONS
C       JACFLG : JACOBIAN FLAG WITH THE FOLLOWING MEANINGS:
C                JACFLG = 1 IF ANALYTIC JACOBIAN SUPPLIED
C                JACFLG = 0 IF ANALYTIC JACOBIAN NOT SUPPLIED
C       GRADTL : TOLERANCE AT WHICH GRADIENT IS CONSIDERED CLOSE ENOUGH
C                TO ZERO TO TERMINATE ALGORITHM
C       STEPTL : TOLERANCE AT WHICH SUCCESSIVE ITERATES ARE CONSIDERED
C                CLOSE ENOUGH TO TERMINATE ALGORITHM
C       FTOL : TOLERANCE AT WHICH FUNCTION VALUE IS CONSIDERED CLOSE
C                ENOUGH TO ZERO
C       METHOD : METHOD TO USE
C                METHOD = 0 : STANDARD METHOD IS USED
C                METHOD = 1 : TENSOR METHOD IS USED
C       GLOBAL : GLOBAL STRATEGY TO USE
C                GLOBAL = 0 : LINE SEARCH
C                GLOBAL = 1 : 2-DIMENSIONAL TRUST REGION
C       STEPMX : MAXIMUM ALLOWABLE STEP SIZE
C       DLT    : TRUST REGION RADIUS
C       IPR    : DEVICE TO WHICH TO SEND OUTPUT
C       X      : ESTIMATE TO A ROOT OF FCN ( USED BY UNCMIN)
C       TYPXU  : TYPICAL SIZE FOR EACH COMPONENT OF X (USED BY UNCMIN)
C       XPLS   : LOCAL MINIMUM OF OPTIMIZATION FUNCTION FCN USED BY
C                UNCMIN
C       GPLS   : GRADIENT AT SOLUTION XPLS (USED BY UNCMIN)
C       A      : WORKSPACE FOR HESSIAN (OR ESTIMATE) (USED BY UNCMIN)
C       WRK    : WORKSPACE (USED BY UNCMIN)
C       WRK1,WRK2,WRK3,WRK4,WRK5,FQ,FQQ:  WORKSPACE
C       FC     : FUNCTION VALUE AT CURRENT ITERATE
C       FHAT   : WORKSPACE
C       DFN    : DIAGONAL SCALING MATRIX FOR F
C       ANLS   : TENSOR TERM MATRIX
C       FV     : WORKSPACE USED TO STORE PAST FUNCTION VALUES
C       AJA    : JACOBIAN MATRIX
C       DN     : STANDARD STEP
C       DT     : TENSOR STEP
C       DF,D,GBAR,DBAR,DBARP : WORKSPACE
C       DXN    : DIAGONAL SCALING MATRIX FOR X
C       S      : MATRIX OF PREVIOUS DIRECTIONS
C       SHAT   : MATRIX OF PAST LINEARLY INDEPENDENT DIRECTIONS
C       CURPOS,PIVOT,PBAR : PIVOT VECTORS
C       SQRN   : MAXIMUM COLUMN DIMENSION OF ANLS, S, AND SHAT
C       EPSM   : MACHINE PRECISION
C       FVEC   : NAME OF SUBROUTINE TO EVALUATE FUNCTION
C       JAC    : (OPTIONAL) NAME OF SUBROUTINE TO EVALUATE JACOBIAN.
C                MUST BE DECLARED EXTERNAL IN CALLING ROUTINE
C
C
C       INPUT-OUTPUT PARAMETERS :
C       ------------------------
C
C       MSG : MESSAGE TO INHIBIT CERTAIN AUTOMATIC CHECKS + OUTPUT
C
C
C       OUTPUT PARAMETERS :
C       -----------------
C
C       XP : SOLUTION TO THE SYSTEM OF NONLINEAR EQUATIONS
C       FP : FUNCTION VALUE AT THE SOLUTION
C       GP : GRADIENT AT THE SOLUTION
C       TERMCD : TERMINATION CODE
*/

#endif
