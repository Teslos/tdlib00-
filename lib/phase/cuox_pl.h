/*
Copyright (C) 1998, 1999, 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                       http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#ifndef __CUOX_PL_H
#define __CUOX_PL_H

#include <intvars.h>

class CuOx_plane : public OneInternalVariable
{
  size_t na_;
  size_t nb_;
  vec_func_Tp f;
  bool old_model;
  FuncTpx ref;
  vec_formula vf;
  vec_string vp;

  virtual double Zex(function f, const StateTp &Tp, const StateX &x) const
    {return ref.Z(f, Tp, x);}
//add to res
  virtual void dZexdx(function f, const StateTp &Tp, const StateX &x, 
      vec_double& res) const
    {ref.dZdx(f, Tp, x, res);}
  
//result in vint
  virtual const vec_double& yeq(const StateTp &Tp, const StateX &x) const
  {
    vint[0] = x_eq(Tp, x[1]);
    return vint;
  }
  
  const double& y() const {return vint[0];}
  
  virtual double Zin(function f, const StateTp &Tp, const StateX &x) const
    {return Zmix(f, Tp, x[1], y());}

//add to res
  virtual void dZindx(function f, const StateTp &Tp, const StateX &x, 
      vec_double &res) const
    {res[1] += dZdz(f, Tp, x[1], y());}
  
  virtual double dydT(const StateTp &Tp, const StateX &x) const
    {return dxdT(Tp, x[1], y());}
  virtual double dydp(const StateTp &Tp, const StateX &x) const
    {return dxdp(Tp, x[1], y());}
//rewrite res
  virtual void dydx(const StateTp &Tp, const StateX &x, vec_double &res) const
  {
    res[0] = 0.;
    res[1] = dxdz(Tp, x[1], y());
  }

  virtual double dZindy(function f, const StateTp &Tp, const StateX &x) const
    {return dZdx_(f, Tp, x[1], y());}

public:
  CuOx_plane() : OneInternalVariable(2, 1), na_(0), nb_(0), f(2), 
    old_model(true), vf(2), vp(2)
  {
    vp[0] = "tetragonal";
    vp[1] = "orthorombic";
    SetIfEqTheSame(true);
    vsint[0] = "x";
  }
  ~CuOx_plane() {}

  virtual void clear(size_t n = 0, size_t m = 0, size_t l = 0)
  {
    na_ = 0;
    nb_ = 0;
    f.clear();
    f.push_back(func_Tp());
    f.push_back(func_Tp());
    ref.clear();
    RefInternalVariables::clear(2, 1);
    vsint[0] = "x";
  }
  virtual CuOx_plane* clone() const {return new CuOx_plane(*this);}

  virtual void read(const SGML &e);
  virtual void WriteBody(ostream& out, size_t shift = 0) const;

  size_t na() const
    {return na_;}
  size_t nb() const
    {return nb_;}
  const func_Tp& fi(size_t i) const
    {return f[i];}
  const func_Tp& ai(size_t i) const
    {return f[2 + i];}
  const func_Tp& bi(size_t i) const
    {return f[2 + na_ + i];}

  double Zmix(function f, const StateTp &Tp,
              const double &z, const double &x) const;
  double dZdx_(function f, const StateTp &Tp,
               const double &z, const double &x) const;
  double dZdz(function f, const StateTp &Tp,
              const double &z, const double &x) const;
  double x_eq(const StateTp &Tp, const double &z) const;
  double dxdz(const StateTp &Tp, const double &z, const double &x) const;
  double dxdT(const StateTp &Tp, const double &z, const double &x) const;
  double dxdp(const StateTp &Tp, const double &z, const double &x) const;

  virtual const vec_formula& comps() const
    {return vf;}
  
  virtual size_t NIntPhases() const {return 2;}
  virtual const vec_string& NamesIntPhases() const
    {return vp;}
  virtual bool IsStable(size_t i, const StateTp &Tp,
      const StateX &x) const;
  virtual const string& StableIntPhase(const StateTp &Tp, 
      const StateX &x) const
  {
    if (x_eq(Tp, x[1]) == 0.)
      return vp[0];
    else
      return vp[1];
  }
  virtual double CheckBoundary(size_t i, size_t j, const StateTp &Tp, 
      const StateX &x) const;
};

#endif
