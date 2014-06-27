/*
Copyright (C) 1998, 1999, 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                       http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#ifndef __CUOX_OR_H
#define __CUOX_OR_H

#include <phase.h>
#include <func_tpx.h>

// x = z/2

class CuOx_ordered_plane : public RefSimple
{
  size_t na_;
  vec_func_Tp f;
  FuncTpx ref;
  size_t mul_ent;

public:
  CuOx_ordered_plane() : RefSimple(2), na_(0), f(2), mul_ent(2) {}
  ~CuOx_ordered_plane() {}

  virtual void clear(size_t n = 0, size_t m = 0, size_t l = 0)
  {
    na_ = 0;
    f.clear();
    f.push_back(func_Tp());
    f.push_back(func_Tp());
    RefSimple::clear(2);
  }
  virtual CuOx_ordered_plane* clone() const 
    {return new CuOx_ordered_plane(*this);}

  virtual void read(const SGML &e);
  virtual void WriteBody(ostream &out, size_t shift = 0) const;

  size_t na() const
    {return na_;}
  const func_Tp& fi(size_t i) const
    {return f[i];}
  const func_Tp& ai(size_t i) const
    {return f[2 + i];}

  double dZdz(function f, const StateTp &Tp, const double &z) const;

  double Zmix(function f, const StateTp &Tp, const double &z) const;

  virtual double Z(function f, index i, const StateTp &Tp,
                   const StateX &x) const
  {
    switch (i)
    {
      case ::full:
        return ref.Z(f, Tp, x) + Zmix(f, Tp, x[1]);
      case ::ref:
        return ref.Z(f, Tp, x);
      case ::mix:
        return Zmix(f, Tp, x[1]);
      default:
        return 0.;
    }
  }
  virtual const vec_double& dZdx(function f, index i, const StateTp &Tp,
                              const StateX &x) const;
  virtual const vec_double& z(function f, index i, const StateTp &Tp,
                              const StateX &x) const;
};

#endif
