/*
Copyright (C) 1998, 1999, 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                       http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#ifndef __INTVARS_H
#define __INTVARS_H

#include <phase.h>
#include <fstream>

class RefInternalVariables : public RefPhase
{
protected:
  bool same;
  
  string FileName;
  bool debug_;
  mutable ofstream *debug;

// vint contains all the information about internal variables, 
// so the real dimension of the problem might be less than vint.size()
  mutable vec_double vint;
  vec_string vsint;
  
  virtual double Zex(function f, const StateTp &Tp, const StateX &x) const = 0;
//add to res
  virtual void dZexdx(function f, const StateTp &Tp, const StateX &x, 
      vec_double& res) const = 0;

//the next two functions take const vec_double &y from vint 
//do not forget to set vint (for example by calling yeq)
  virtual double Zin(function f, const StateTp &Tp, const StateX &x) const = 0;

//add to res
  virtual void dZindx(function f, const StateTp &Tp, const StateX &x, 
      vec_double &res) const = 0;
  
//result in vint - do not forget to put to vint all the information
  virtual const vec_double& yeq(const StateTp &Tp, const StateX &x) const = 0;

  virtual void CopyInternalVariables(const vec_double &y) const
  {
    if (vint.size() != y.size())
      throw gError("RefInternalVariables: dimension of y != vint");
    vint = y;
  }

public:
  static double step;
  
  RefInternalVariables() : same(false), debug(0) {}
  RefInternalVariables(size_t n, size_t m) : same(false), debug(0) {clear(n, m);}
  ~RefInternalVariables() {DeleteDebug();}
  
// n - number of components, m - dimension of vint (all the information)
  void virtual clear(size_t n = 0, size_t m = 0, size_t l = 0)
  {
    RefPhase::clear(n);
    vint.resize(m);
    vsint.resize(m);
  }

  void DeleteDebug();
  void SetDebug(const SGML &e);
  virtual void WriteAttributes(ostream &out, size_t shift = 0) const
  {
    out << endl << PutTab(shift) << "debug=" << debug_;
    if (!FileName.empty())
      out << endl << PutTab(shift) << "DebugFile=" << FileName;
  }
  
  void SetIfEqTheSame(bool t) {same = t;}

  virtual size_t NIntVars() const {return vint.size();}
  virtual const vec_string& NamesIntVars() const {return vsint;}
  virtual const vec_double& IntVarsEq(const StateTp &Tp, const StateX &x) const 
  {
    return yeq(Tp, x);
  }

  virtual double Z(function f, index i, const StateTp &Tp, 
      const StateX &x) const
  {
    return Z_y(f, i, Tp, x, yeq(Tp, x));
  }
  virtual const vec_double& dZdx(function f, index i, const StateTp &Tp, 
      const StateX &x) const
  {
    return dZdx_y(f, i, Tp, x, yeq(Tp, x));
  }
};

class OneInternalVariable : public RefInternalVariables
{
protected:
//the next functions take const vec_double &y from vint 
//do not forget to set vint (for example by calling yeq)
  virtual double dEqdy(const StateTp &Tp, const StateX &x) const 
    {throw gError("OneInternalVariable: dEqdy is not defined");}
  virtual double dEqdT(const StateTp &Tp, const StateX &x) const
    {throw gError("OneInternalVariable: dEqdT is not defined");}
  virtual double dEqdp(const StateTp &Tp, const StateX &x) const
    {throw gError("OneInternalVariable: dEqdp is not defined");}
//rewrite res
  virtual void dEqdx(const StateTp &Tp, const StateX &x, vec_double &res) const
    {throw gError("OneInternalVariable: dEqdx is not defined");}
  virtual double dydT(const StateTp &Tp, const StateX &x) const
    {return -dEqdT(Tp, x)/dEqdy(Tp, x);}
  virtual double dydp(const StateTp &Tp, const StateX &x) const
    {return -dEqdp(Tp, x)/dEqdy(Tp, x);}
//rewrite res
  virtual void dydx(const StateTp &Tp, const StateX &x, vec_double &res) const
  {
    dEqdx(Tp, x, res);
    double d = dEqdy(Tp, x);
    for (size_t i = 0; i < size(); ++i)
      res[i] = - res[i]/d;
  }

  virtual double dZindy(function f, const StateTp &Tp, 
      const StateX &x) const = 0;

public:
  OneInternalVariable(size_t n = 0, size_t m = 0) 
    : RefInternalVariables(n, m) {}
  
// const vec_double y is considered to be at equilibrium at Tp and x
// it must be obtained by yeq(Tp, x)
  virtual double Z_y(function f, index i, const StateTp &Tp, 
      const StateX &x, const vec_double &y) const;
  virtual const vec_double& dZdx_y(function f, index i, const StateTp &Tp, 
      const StateX &x, const vec_double &y) const;
};

class InternalVariables : public RefInternalVariables
{
protected:
  mutable vec_double e;
  mutable vec_int ipiv;
  mutable matrix dedy;
  mutable matrix res;

//the next functions take const vec_double &y from vint 
//do not forget to set vint (for example by calling yeq)

//result in dedy
  virtual matrix& dEqdy(const StateTp &Tp, const StateX &x) const
    {throw gError("InternalVariables: dEqdy is not defined");}
//result in res
  virtual vec_double&  dEqdT(const StateTp &Tp, const StateX &x) const
    {throw gError("InternalVariables: dEqdy is not defined");}
//result in res
  virtual vec_double&  dEqdp(const StateTp &Tp, const StateX &x) const
    {throw gError("InternalVariables: dEqdy is not defined");}
//result in res
  virtual matrix& dEqdx(const StateTp &Tp, const StateX &x) const
    {throw gError("InternalVariables: dEqdy is not defined");}
  
//result in res
  virtual const vec_double& dydT(const StateTp &Tp, const StateX &x) const;
  virtual const vec_double& dydp(const StateTp &Tp, const StateX &x) const;
  virtual const matrix& dydx(const StateTp &Tp, const StateX &x) const;
  
//result in e
  virtual const vec_double& dZindy(function f, const StateTp &Tp, 
      const StateX &x) const = 0;
  
public:
  InternalVariables(size_t n = 0, size_t m = 0, size_t l = 0) 
    : RefInternalVariables(n, m), e(l), ipiv(l), dedy(l, l), res(l, n)  {}
  
// n - number of components, m - dimension of vint (all the information)
// l - is the number of internal variables to solve
  void clear(size_t n = 0, size_t m = 0, size_t l = 0)
  {
    RefInternalVariables::clear(n, m, l);
    e.resize(l);
    ipiv.resize(l);
    dedy.resize(l, l);
    res.resize(l, n);
  }

// const vec_double y is considered to be at equilibrium at Tp and x
// it must be obtained by yeq(Tp, x)
  virtual double Z_y(function f, index i, const StateTp &Tp, 
      const StateX &x, const vec_double &y) const;
  virtual const vec_double& dZdx_y(function f, index i, const StateTp &Tp, 
      const StateX &x, const vec_double &y) const;
};

#endif
