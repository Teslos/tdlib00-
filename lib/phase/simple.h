/*
Copyright (C) 1998, 1999, 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                       http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#ifndef __SIMPLE_H
#define __SIMPLE_H

#include <func_tpx_imp.h>
#include <phase.h>

class index_x
{
  size_t m;     // number of interactions
  size_t n;     // number of components
  size_t nm;    // maximum number of different sets
  vec_int in;   // vector of indeces in[m]
  size_t ilin;  // number by order 0, ..., nm - 1

public:

  index_x(size_t m, size_t n) : m(m), n(n), in(m)
  {
    if (!(m > 1 && n > 1 && m <= n)) throw gError("index_x - wrong");
    init();
  }
  int operator[](size_t i) const
    {return in[i];}
  void init();
  void operator++();
  operator vec_int&() {return in;}
  operator int() const {return nm-ilin;}
};

class SimpleSolution : public RefSimple
{
  vec_FuncTpx p;
  mutable SymMatrix mat;

  void insert(const FuncTpx &y);
public:

  SimpleSolution() {}
  virtual ~SimpleSolution() {}

  void set(size_t nc, size_t m_max = 0);

  virtual void clear(size_t n = 0, size_t m = 0, size_t l = 0)
  {
    mat.clear();
    p.clear(); 
    RefSimple::clear();
  }
  virtual SimpleSolution* clone() const {return new SimpleSolution(*this);}

  virtual void read(const SGML &e);
  virtual void WriteBody(ostream &out, size_t shift = 0) const;

  const IdealMixing& FindIdealMixing();
 
  virtual double Z(function f, index i, const StateTp &Tp, 
      const StateX &x) const;
  virtual const vec_double& dZdx(function f, index i, const StateTp &Tp, 
      const StateX &x) const;
  const SymMatrix& d2Zdx2(function f, index i, const StateTp &Tp, 
      const StateX &x) const;
};

#endif
