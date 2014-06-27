/*
Copyright (C) 1998, 1999, 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                       http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#ifndef __PHASE_EQ_H
#define __PHASE_EQ_H

#include <iomanip>
#include <algorithm>
#include <toms.h>
#include <td_algo.h>

void GuessByGrid(const VOF &fnc, long n, long np, double *x, const double *h,
                   double *fv, size_t N = 5);

typedef vector<StateX> vec_statex;
typedef vec_statex::iterator vec_statex_i;
typedef vec_statex::const_iterator vec_statex_ci;

class PhaseEquilibrium : public RefAlgorithm
{
  enum StatusOfVariableBase {unknown, constraint, HardConstraint, dependent, 
    functional, value};
  static const char *StatusOfVariableName[];
  typedef Enum<StatusOfVariableBase, StatusOfVariableName> StatusOfVariable;

  vec_phase sys;
  vec_string nm_sys;
  StateTp Tp;
  vec_statex vx;

  vec_string nm_st;
  vec_string nm_st2;
  vector<StatusOfVariable> fix;
  vec_double_2ptr st;

  mutable bool IsXSet;
  mutable bool solved;
  mutable enum {NoEst, PartEst, FullEst} IsEst;
  mutable bool ThrowException;
  mutable bool SaveSolution;

  matrix eq;
  vec_int IsBasis;
  mutable vec_double fv;
  mutable vec_double mu;
  mutable vec_double mu_basis;
  mutable double fmin;
  mutable vec_double prop;
  mutable long Neq; // = eq.NCols() + prop.size()

  mutable vec_double x_;
  mutable long Nx; // Nx = xest.size();
  mutable vec_double xest;
  mutable vec_double h;
  mutable vec_double l;
  mutable vec_double u;
  mutable vector <long> nbd;
  mutable vec_int iun;

  mutable vec_double ph;
  mutable vec_double pl;
  mutable vec_double pu;

  vec_convert vc;
  vec_int i_vc;
  vec_int i_dmf;

//for output
  vec_double prop_in;
  vec_double st_in;
  
  enum source {Property, Partial, AllState, Residual};
  struct descriptor
  {
    string name;
    function f;
    index i;
    int ix;
    int ip;
    size_t N;
    source from;
  };
  mutable vector<descriptor> vd;

  int SearchState(const string& str) const
  {
    for (size_t j = 0; j < nm_st.size(); ++j)
      if (nm_st[j] == str || nm_st2[j] == str)
        return j;
    return -1;
  }

  void fun(long *n, double *x, double *f, long *np) const;
  void by_grid() const
  {
    if (debug)
      *debug << "search of the etismates by grid" << endl;
    x_ = l;
    int dim;
    if (x_.size() == 1)
      dim = 15;
    else if (x_.size() == 2)
      dim = 5;
    else
      dim = 3;
    GuessByGrid(makeFunctor((VOF*)0, *this, &PhaseEquilibrium::fun), 
        Nx, Neq, &*x_.begin(), &*h.begin(), &*fv.begin(), dim);
  }
  void pby_grid() const
  {
    if (debug)
      *debug << "search of the partial etismates by grid" << endl;
    copy(pl.begin(), pl.end(), x_.begin());
    int dim;
    if (pl.size() == 1)
      dim = 15;
    else if (pl.size() == 2)
      dim = 5;
    else
      dim = 3;
    GuessByGrid(makeFunctor((VOF*)0, *this, &PhaseEquilibrium::fun),
        pl.size(), Neq, &*x_.begin(), &*ph.begin(), &*fv.begin(), dim);
  }
  void lbfgsb() const
  {
    if (debug)
      *debug << "search of the solution by lbfgsb" << endl;
    fmin = lbfgsb_SS(makeFunctor((VOF*)0, *this, &PhaseEquilibrium::fun), 
        Nx, Neq, &*x_.begin(), &*l.begin(), &*u.begin(), &*nbd.begin(), &set);
    fmin = sqrt(fmin);
  }
  void plbfgsb() const
  {
    if (debug)
      *debug << "search of the partial solution by lbfgsb" << endl;
    fmin = lbfgsb_SS(makeFunctor((VOF*)0, *this, &PhaseEquilibrium::fun),
        ph.size(), Neq, &*x_.begin(), &*pl.begin(), &*pu.begin(), 
        &*nbd.begin(), &set);
    fmin = sqrt(fmin);
  }
  void solver() const;

  void SetX() const;

  void SetDependent() const;
  void SetPhases();
  void SetEq();
  void SetVar(const vec_double &lo, const vec_double &up, 
      const vec_convert &vctmp);

public:
  static double eps;
  static double penalty;
  static double Tmin;
  static double Tmax;
  static SetLbfgsb GlSet;
  SetLbfgsb set;

  PhaseEquilibrium() 
  {
    set = GlSet;
    Nx = Neq = 0;
    IsXSet = false;
    solved = false;
  }
  PhaseEquilibrium(const PhaseEquilibrium &old) : RefAlgorithm(old)
  {
    operator=(old);
  }
  PhaseEquilibrium& operator=(const PhaseEquilibrium &old);

  virtual void clear();
  virtual PhaseEquilibrium* clone() const 
    {return new PhaseEquilibrium(*this);}

  virtual void read(const SGML &e);
  virtual void WriteAttributes(ostream &out, size_t shift = 0) const
  {
    RefAlgorithm::WriteAttributes(out, shift);
    out << endl << PutTab(shift) << "SaveSolution=" << SaveSolution;
    out << endl << PutTab(shift) << "ThrowException=" << ThrowException;
  }
  virtual void WriteBody(ostream &out, size_t shift = 0) const;

  virtual int GetInputIndex(const string& s) const;
  virtual void SetVariable(int i, const double &x) const;
  virtual const string& name(int i) const;
  virtual int GetOutputIndex(const string& s, size_t &N) const;
  virtual void OutNames(int i, vec_string &vs) const;
  void OutProperty(int i, double *x) const;
  virtual void out(int i, double *x) const;

  virtual void SetSolved(bool t) const
  {
    solved = t;
  }
};

#endif
