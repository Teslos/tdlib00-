/*
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#ifndef __FUNC_TPX_IMP_H
#define __FUNC_TPX_IMP_H

#include <func_tpx.h>

class NullFuncTpx : public RefFuncTpx
{
public:
  NullFuncTpx() {}

  virtual NullFuncTpx* clone() const
  {
    return new NullFuncTpx(*this);
  }

  virtual void read(const SGML &e, const vec_formula *f = 0) {}
  virtual void WriteBody(ostream &out, const vec_formula *f = 0, 
      size_t shift = 0) const {}

  virtual double Z(function f, const StateTp &Tp, const StateX &x) const
    {return 0.;}
  virtual void dZdx(function f, const StateTp &Tp,
                    const StateX &x, vec_double &res) const {}
  virtual void d2Zdx2(function f, const StateTp &Tp,
                    const StateX &x, SymMatrix &res) const {}
};

class IdealMixing : public RefFuncTpx
{
  vec_coef vc;
public:
  IdealMixing() 
  {
    SetIndex(::ideal);
  }

  virtual void clear()
  {
    vc.clear();
    RefFuncTpx::clear();
  }
  virtual IdealMixing* clone() const
  {
    return new IdealMixing(*this);
  }

  virtual void read(const SGML &e, const vec_formula *f = 0);
  virtual void WriteBody(ostream &out, const vec_formula *f = 0, 
      size_t shift = 0) const;

  virtual void set(const vec_int &vec, const func_Tp &f = func_Tp("Cp_zero"), 
      size_t power = 0)
  {
    SetVars(vec);
    vc.resize(vec.size());
    fill(vc.begin(), vc.end(), coef(1.));
    for (size_t i = 0; i < size(); ++i)
      vars[i] = i;
  }

  const vec_coef& GetCoef() const
  {
    return vc;
  }

  virtual double Z(function f, const StateTp &Tp, const StateX &x) const;
  virtual void dZdx(function f, const StateTp &Tp,
                    const StateX &x, vec_double &res) const;
  virtual void d2Zdx2(function f, const StateTp &Tp,
                    const StateX &x, SymMatrix &res) const;
};

class Reference : public RefFuncTpx
{
  vec_species sps;
public:
  Reference() 
  {
    SetIndex(::ref);
  }

  virtual void clear()
  {
    sps.clear();
    RefFuncTpx::clear();
  }
  virtual Reference* clone() const
  {
    return new Reference(*this);
  }

  virtual void read(const SGML &e, const vec_formula *f = 0);
  virtual void WriteBody(ostream &out, const vec_formula *f = 0, 
      size_t shift = 0) const;

  virtual void set(const vec_int &vec, const func_Tp &f = func_Tp("Cp_zero"), 
      size_t power = 0)
  {
    SetVars(vec);
    sps.resize(vec.size());
  }
  const vec_species& GetSpecies() const
  {
    return sps;
  }
  virtual double Z(function f, const StateTp &Tp, const StateX &x) const;
  virtual void dZdx(function f, const StateTp &Tp, const StateX &x, 
      vec_double &res) const;
  virtual void d2Zdx2(function f, const StateTp &Tp, const StateX &x, 
      SymMatrix &res) const {}
  void z(function f, const StateTp &Tp, vec_double &res) const;
};

#endif

