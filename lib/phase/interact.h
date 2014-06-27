/*
Copyright (C) 1998, 1999, 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                       http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#ifndef __INTERACT_H
#define __INTERACT_H

#include <func_tpx.h>
#include <iterator>

class interaction : public RefFuncTpx
{
protected:
  size_t m;       // m = vars.size();
  size_t p;       // p = f.size();
  vec_func_Tp f;  // f.size() is maximum number of polynomial members
  formalism frm;
  mutable vec_double dSdv;
  mutable vec_double dSdx;
  mutable SymMatrix d2Sdv2;
  mutable SymMatrix d2Sdx2;
  mutable vec_double dvdx_;
  mutable vec_double d2vdx2_;
  mutable vec_double vec;
  mutable double SumX;
  mutable double Sum;

// these functions are inlined in interact.cpp
  bool check(const StateX &x) const;
  void SetVec(const StateX &x) const;
  void Setdvdx(const StateX &x) const; // after SetVec
  const double& dvdx(size_t i, size_t j) const;
  void Setd2vdx2(const StateX &x) const; // after Setdvdx
  const double& d2vdx2(size_t i, size_t j, size_t k) const;

public:
  interaction() {SetIndex(::excess); m = p = 0; frm = ::Muggianu;}
  explicit interaction(size_t n)
    {SetIndex(::excess); clear(n);}

  virtual void clear(size_t i = 0)
  {
    RefFuncTpx::clear(i);
    p = 0;
    f.clear();
    dSdv.resize(i);
    dSdx.resize(i);
    d2Sdv2.resize(i);
    d2Sdx2.resize(i);
    dvdx_.resize(2*i);
    d2vdx2_.resize(3*i);
    vec.resize(i);
    frm = Muggianu;
    m = i;
  }

  void ReadFormalism(const SGML& e)
  {
    string f = e.FindString("formalism");
// compatibility    
    if (f == "all_comps")
      f = "NoChange";
    if (!f.empty())
      frm = f;
  }
  virtual void WriteAttributes(ostream &out, const vec_formula *f = 0,
      size_t shift = 0) const 
  {
    out << " formalism=" << frm;
  }
  
  virtual formalism SetFormalism(formalism frm_) 
  {
    formalism old = frm;
    frm = frm_;
    return old;
  }
};

class RedlichKister : public interaction
{
  formalism old;
  
// these functions are inlined in interact.cpp
  void SetSum(function fun, const StateTp &Tp, const StateX &x) const;
  void SetdSdx(function fun, const StateTp &Tp, const StateX &x) const;
  void Setd2Sdx2(function fun, const StateTp &Tp, const StateX &x) const;

public:
  RedlichKister() : interaction(2) {old = frm;}

  virtual RedlichKister* clone() const
  {
    return new RedlichKister(*this);
  }

  virtual void clear()
  {
    interaction::clear(2);
    old = frm;
  }
  virtual void set(const vec_int &vec, const func_Tp &fun, size_t power = 0)
  {
    clear();
    if (vec.size() == 2)
    {
      vars = vec;
      p = power + 1;
      f.insert(f.begin(), p, fun);
    }
    else
      throw gError("Redlich_Kister: m > 2");
  }

  virtual void read(const SGML &e, const vec_formula *f = 0);
  virtual void WriteBody(ostream &out, const vec_formula *f = 0, 
      size_t shift = 0) const;
  virtual void WriteAttributes(ostream &out, const vec_formula *f = 0,
      size_t shift = 0) const 
  {
    out << " formalism=" << old;
  }

  virtual double Z(function f, const StateTp &Tp, const StateX &x) const;
  virtual void dZdx(function f, const StateTp &Tp,
                    const StateX &x, vec_double &res) const;
  virtual void d2Zdx2(function f, const StateTp &Tp,
                    const StateX &x, SymMatrix &res) const;
};

class HochArpshofen : public interaction
{
  vec_int pwr;

// these functions are inlined in interact.cpp
  void SetSum(function fun, const StateTp &Tp, const StateX &x) const;
  void SetdSdx(function fun, const StateTp &Tp, const StateX &x) const;
  void Setd2Sdx2(function fun, const StateTp &Tp, const StateX &x) const;

public:
  HochArpshofen() : interaction(2) {}
  ~HochArpshofen() {}

  virtual HochArpshofen* clone() const
  {
    return new HochArpshofen(*this);
  }

  virtual void clear()
  {
    interaction::clear(2);
    pwr.clear();
  }
  virtual void set(const vec_int &vec, const func_Tp &fun, size_t power = 0)
  {
    clear();
    if (vec.size() == 2)
    {
      vars = vec;
      p = power + 1;
      f.insert(f.end(), p, fun);
      for (size_t i = 0; i < f.size(); ++i)
        pwr.push_back(i + 2);
    }
    else
      throw gError("Hoch_Arpshofen: m > 2");
  }

  virtual void read(const SGML &e, const vec_formula *f = 0);
  virtual void WriteBody(ostream &out, const vec_formula *f = 0, 
      size_t shift = 0) const;

  virtual double Z(function f, const StateTp &Tp, const StateX &x) const;
  virtual void dZdx(function f, const StateTp &Tp,
                    const StateX &x, vec_double &res) const;
  virtual void d2Zdx2(function f, const StateTp &Tp,
                    const StateX &x, SymMatrix &res) const;
};

class Polynomial : public interaction
{
  vec_int pwr;
  size_t power;

// these functions are inlined in interact.cpp
  void SetSum(function fun, const StateTp &Tp, const StateX &x) const;
  void SetdSdx(function fun, const StateTp &Tp, const StateX &x) const;
  void Setd2Sdx2(function fun, const StateTp &Tp, const StateX &x) const;

public:
  class iterator : public std::iterator<forward_iterator_tag, func_Tp>
  {
  protected:
    vec_int indx;

  public:
    friend class Polynomial;

    iterator() {}
    bool operator==(const iterator &y) const
      {return indx == y.indx;}
    vec_int& operator*() {return indx;}
    iterator& operator++();
    iterator operator++(int);
    size_t operator[](size_t i) const {return indx[i];}
    vec_int& index() {return indx;}
  };
  iterator begin();
  iterator end();

  Polynomial() : power(0) {}
  virtual void clear(size_t m)
  {
    power = 0;
    pwr.clear();
    interaction::clear(m);
  }
  virtual Polynomial* clone() const
  {
    return new Polynomial(*this);
  }
  virtual void set(const vec_int &vct, const func_Tp &fun, size_t power_ = 0)
  {
    clear(vct.size());
    if (m > 0)
    {
      vars = vct;
      power = power_;
      if (m == 1) power = 0;
      p = bin_coef(m - 1, power + m - 1);
      f.insert(f.begin(), p, fun);
      iterator i = begin();
      while (true)
      {
        pwr.insert(pwr.end(), (*i).begin(), (*i).end());
        if (i == end())
          break;
        ++i;
      }
    }
  }

  virtual void read(const SGML &e, const vec_formula *f = 0);
  virtual void WriteBody(ostream &out, const vec_formula *f = 0, 
      size_t shift = 0) const;

  virtual double Z(function f, const StateTp &Tp, const StateX &x) const;
  virtual void dZdx(function f, const StateTp &Tp,
                    const StateX &x, vec_double &res) const;
  virtual void d2Zdx2(function f, const StateTp &Tp,
                    const StateX &x, SymMatrix &res) const;
};

#endif

