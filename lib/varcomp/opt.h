/*
Copyright (C) 1994 - 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#ifndef _OPT_H___
#define _OPT_H___

#include "sumsqr.h"
#include <toms.h>
#include <species.h>

class RefOptimizer
{
  RefOptimizer(const RefOptimizer &old)
  {
    throw gError("RefOptimizer: coping is not allowed");
  }
  RefOptimizer& operator=(const RefOptimizer &old)
  {
    throw gError("RefOptimizer: coping is not allowed");
  }
public:
  RefOptimizer() {}
  virtual ~RefOptimizer() {}
  virtual void set(size_t n, size_t m, const SGML &e) = 0;
  virtual void solve(double *x_, double *l, double *u) = 0;
  virtual void print(ostream &out) = 0;
  virtual matrix& jac() = 0;
};

class optimizer
{
  RefOptimizer *ptr;
public:
  static map<string, SGML> opt;
  static void (*InitSGML[])();
  
  optimizer();
  static string algorithm;
  static void ReadSGML(const SGML &e);
  static void WriteSGML(ostream &out, size_t shift = 0);
  
  void init(size_t n, size_t m);
  void solve(double *x, double *l, double *u) {ptr->solve(x, l, u);}
  void print(ostream &out) {ptr->print(out);}
  matrix& jac() {return ptr->jac();}
};

class tensolve : public RefOptimizer
{
  int maxm;
  int maxn;
  int m; // DIMENSIONS OF PROBLEM
  int n;
  vec_double typex; // TYPX(N) - TYPICAL SIZE FOR EACH COMPONENT OF X
  vec_double typef; // TYPF(M) - TYPICAL SIZE FOR EACH COMPONENT OF F
  int itnlim; // MAXIMUM NUMBER OF ALLOWABLE ITERATIONS
  
  int jacflg; // JACFLG = 0 IF ANALYTIC JACOBIAN NOT SUPPLIED
  double gradtl; // tolerances for gradient, step and function
  double steptl; 
  double ftol;
  int method; // METHOD = 1 : TENSOR METHOD IS USED
  int global; // GLOBAL = 0 : LINE SEARCH, 
              //GLOBAL = 1 : 2-DIMENSIONAL TRUST REGION
  double stepmx; // MAXIMUM ALLOWABLE STEP SIZE
  double dlt; // TRUST REGION RADIUS
  double dltin;
  int ipr; // DEVICE TO WHICH TO SEND OUTPUT
  vec_double x; // X(SQRN) - ESTIMATE TO A ROOT OF FCN ( USED BY UNCMIN)
  vec_double typxu; // TYPXU(SQRN) - TYPICAL SIZE FOR EACH COMPONENT OF X (USED BY UNCMIN)
  vec_double xpls; // XPLS(SQRN) - LOCAL MINIMUM OF OPTIMIZATION FUNCTION FCN USED BY UNCMIN
  vec_double gpls; // GPLS(SQRN) - GRADIENT AT SOLUTION XPLS (USED BY UNCMIN)
  vec_double a; // A(MAXP,SQRN) - WORKSPACE FOR HESSIAN (OR ESTIMATE) (USED BY UNCMIN)
  vec_double wrk; // WRK(MAXP,SQRN) - WORKSPACE (USED BY UNCMIN)
  vec_double dfn; // DFN(M) - DIAGONAL SCALING MATRIX FOR F
  vec_double wrk1; // WRK1(M) - workspaces
  vec_double wrk2; // WRK2(M)
  vec_double wrk3; // WRK3(M)
  vec_double wrk4; // WRK4(M)
  vec_double wrk5; //WRK5(M)
  vec_double fq; // FQ(M)
  vec_double fqq; // FQQ(M)
  vec_double fc; // FC(M) - FUNCTION VALUE AT CURRENT ITERATE
  vec_double fhat; //FHAT(M) - WORKSPACE
  vec_double anls; // ANLS(MAXM,SQRN) - TENSOR TERM MATRIX
  vec_double fv; // FV(MAXM,SQRN) - WORKSPACE USED TO STORE PAST FUNCTION VALUES
  matrix aja; // AJA(MAXM,N) - JACOBIAN MATRIX
  vec_double dxn; // DXN(N) - DIAGONAL SCALING MATRIX FOR X
  vec_double dn; // DN(N) - STANDARD STEP
  vec_double dt; // DT(N) - TENSOR STEP
  vec_double df; // DF(N) - workspaces
  vec_double d; // D(N)
  vec_double gbar; // GBAR(N)
  vec_double dbar; // DBAR(N)
  vec_double dbarp; // DBARP(N)
  vec_double s; // S(MAXN,SQRN) - MATRIX OF PREVIOUS DIRECTIONS
  vec_double shat; // SHAT(MAXN,SQRN) - MATRIX OF PAST LINEARLY INDEPENDENT DIRECTIONS
  vec_int curpos;     // CURPOS(N) - pivots
  vec_int pivot;    // PIVOT(N)
  vec_int pbar;       // PBAR(N)
  double epsm; // MACHINE PRECISION
  int sqrn; // MAXIMUM COLUMN DIMENSION OF ANLS, S, AND SHAT
  int msg; // MESSAGE TO INHIBIT CERTAIN AUTOMATIC CHECKS + OUTPUT
  vec_double xp; // XP(N)
  vec_double fp; // FP(M)
  vec_double gp; // GP(N)
  int termcd;
  
  static double ssbest;
  static vec_double xbest;

  static void ssn(double *x, double *f, int *m, int *n);
public:
  static void AddSGML();
  
  tensolve() : m(0), n(0), sqrn(0) {}
  
  virtual void set(size_t n_, size_t m_, const SGML &e);
  virtual void solve(double *x_, double *l, double *u);
  virtual void print(ostream &out);
  virtual matrix& jac()
  {
    return aja;
  }
};

class LBFGSB : public RefOptimizer
{
  int m;
  int n;
  int ier;
  matrix xjac;
  double ssq;
  vec_double f;
  vec_double work;
  vector<long> nbd;
  char task[60];
  char csave[60];
  long lsave[4];
  long isave[44];
  double dsave[29];
  vec_double grad;
  vec_double xold;
  vec_double fn2;
  vector<long> iwa;
  long opt_m;
  double factor;
  double pgtol;
  long iprint;
  int iter;
  double dsig;
  double step;
  double small;

public:
  static void AddSGML();
  
  LBFGSB() : m(0), n(0) {}
  
  virtual void set(size_t n_, size_t m_, const SGML &e);
  virtual void solve(double *x, double *l, double *u);
  virtual void print(ostream &out);
  virtual matrix& jac()
  {
    return xjac;
  }
};

class Hooke : public RefOptimizer
{
  int m;
  int n;
  matrix xjac;
  vec_double z;
  vec_double fv;
  vec_double delta;
  vec_double xbefore;
  vec_double newx;
  vec_double endpt;
  double rho;
  double eps;
  int iter;
  int ii;
  
  double best_nearby(vec_double &delta, vec_double &point, double prevbest, 
      int nvars)
  {
    return best_nearby(&*delta.begin(), &*point.begin(), prevbest, nvars);
  }
  double best_nearby(double *delta, double *point, double prevbest, int nvars);
  int hooke(int nvars, double *startpt, double *endpt, double rho, 
    double epsilon, int itermax);

public:
  double f(double *x, int n)
  {
    ss(x, m, n, &*fv.begin());
    double sum = 0.;
    for (int i = 0; i < m; ++i)
      sum += fv[i]*fv[i];
    return sum;
  }
  double f(vec_double &xv, int n)
  {
    return f(&*xv.begin(), n);
  }
  static void AddSGML();
  
  Hooke() : m(0), n(0) {}
  
  virtual void set(size_t n_, size_t m_, const SGML &e);
  virtual void solve(double *x, double *l, double *u)
  {
    ii = hooke(n, x, &*endpt.begin(), rho, eps, iter);
    for (int i = 0; i < n; ++i)
      x[i] = endpt[i];
  }
  virtual void print(ostream &out);
  virtual matrix& jac()
  {
    return xjac;
  }
};

#ifdef INCLUDE_ZXSSQ
#include <../imsl_c/imsl_c.h>

class ZXSSQ : public RefOptimizer
{
  int m;
  int n;
  int ier;
  matrix xjac;
  double ssq;
  vec_double f;
  vec_double work;
  vec_double xjtj;
  int nsig;
  double eps;
  int delta;
  int maxfn;
  int iopt;
  int ixjac;
  int infer;
  double parm[4];
public:
  static void AddSGML();
  
  ZXSSQ() : m(0), n(0) {}
  
  virtual void set(size_t n_, size_t m_, const SGML &e);
  virtual void solve(double *x, double *l, double *u)
  {
    zxssq(ss, m, n, nsig, eps, delta, maxfn, iopt, parm, x, ssq, 
        &*f.begin(), &*xjac.begin(), ixjac, &*xjtj.begin(), &*work.begin(), 
        infer, ier);
  }
  virtual void print(ostream &out);
  virtual matrix& jac()
  {
    return xjac;
  }
  
};
#endif

#endif
