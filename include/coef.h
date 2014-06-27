/*
Copyright (C) 1998, 1999, 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                       http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#ifndef __COEF_H
#define __COEF_H

#include <calc.h>

struct coefficient
{
  double x;
  double *ptr;
  double lower;
  double upper;
  double scale;
  bool unknown;
  mutable bool print;
  void clear()
  {
    x = 0.;
    ptr = &x;
    lower = -HUGE_VAL;
    upper = HUGE_VAL;
    scale = 1.;
    unknown = false;
    print = true;
  }
  coefficient()
  {
    clear();
  }
  coefficient(const coefficient &x_)
  {
    operator=(x_);
  }
  coefficient& operator=(const coefficient &x_)
  {
    if (this != &x_)
    {
      x = x_.x;
      ptr = &x;
      lower = x_.lower;
      upper = x_.upper;
      scale = x_.scale;
      unknown = x_.unknown;
      print = x_.print;
    }
    return *this;
  }
};

class UnknownCoefs
{
  map<string, coefficient, less<string> > coefs;
  mutable vec_double x;
  mutable vec_double l;
  mutable vec_double u;
  mutable vec_double sc;

public:
  typedef map<string, coefficient, less<string> >::iterator iterator; 
  typedef map<string, coefficient, less<string> >::const_iterator const_iterator;

  iterator begin() {return coefs.begin();}
  const_iterator begin() const {return coefs.begin();}
  iterator end() {return coefs.end();}
  const_iterator end() const {return coefs.end();}
  
  UnknownCoefs() {}
  bool add(const string &s, const coefficient &x)
  {
    pair<iterator, bool> i;
    i = coefs.insert(pair<const string, coefficient>(s, x));
    return i.second;
  }
  bool defined(const string &s) const
  {
    if (coefs.find(s) != coefs.end())
      return true;
    else
      return false;
  }
  coefficient& find(const string &s)
  {
    iterator i = coefs.find(s);
    if (i != coefs.end())
    {
      return (*i).second;
    }
    else
      throw gError(string("UnknownCoefs: coef ").append(s).
          append(" is not defined"));
  }
  const coefficient& find(const string &s) const
  {
    const_iterator i = coefs.find(s);
    if (i != coefs.end())
    {
      return (*i).second;
    }
    else
      throw gError(string("UnknownCoefs: coef ").append(s).
          append(" is not defined"));
  }
  void Unknowns(vec_double &x, vec_double &l, vec_double &u, vec_double &sc,
      vec_string &vs);
  void SetUnknowns(double *x);
};

struct ComputedCoefficient : public calculator
{
  UnknownCoefs *u;
  string s;
  double x;
  double *ptr;
  mutable double print;
  ComputedCoefficient() : u(0), x(0.), ptr(&x), print(true) {}
  ComputedCoefficient(const ComputedCoefficient &x)
  {
    operator=(x);
  }
  ComputedCoefficient& operator=(const ComputedCoefficient &x_)
  {
    if (this != &x_)
    {
      s = x_.s;
      x = x_.x;
      print = x_.print;
      ptr = &x;
    }
    return *this;
  }
  virtual void AnalizeId(const string &id)
  {
    coefficient& c = u->find(id);
    curr_tok.key = DOUBLE_PTR_PTR;
    curr_tok.x_ptr_ptr = &c.ptr;
  }
  void precompile(UnknownCoefs *u_)
  {
    u = u_;
    FromString(s);
  }
  void est()
  {
    x = calculator::est();
  }
};

class ComputedCoefs
{
  map<string, ComputedCoefficient, less<string> > coefs;
  
public:
  typedef map<string, ComputedCoefficient, less<string> >::iterator iterator;
  typedef map<string, ComputedCoefficient, less<string> >::const_iterator const_iterator;

  iterator begin() {return coefs.begin();}
  const_iterator begin() const {return coefs.begin();}
  iterator end() {return coefs.end();}
  const_iterator end() const {return coefs.end();}

  bool add(const string &s, const ComputedCoefficient &x)
  {
    pair<iterator, bool> i;
    i = coefs.insert(pair<const string, ComputedCoefficient>(s, x));
    return i.second;
  }
  bool defined(const string &s) const
  {
    if (coefs.find(s) != coefs.end())
      return true;
    else
      return false;
  }
  ComputedCoefficient& find(const string &s)
  {
    iterator i = coefs.find(s);
    if (i != coefs.end())
    {
      return (*i).second;
    }
    else
      throw gError(string("ComputedCoefs: coef ").append(s).append(" is not defined"));
  }
  const ComputedCoefficient& find(const string &s) const
  {
    const_iterator i = coefs.find(s);
    if (i != coefs.end())
    {
      return (*i).second;
    }
    else
      throw gError(string("ComputedCoefs: coef ").append(s).append(" is not defined"));
  }
  void precompile(UnknownCoefs *u)
  {
    for (iterator i = coefs.begin(); i != coefs.end(); ++i)
      (*i).second.precompile(u);
  }
  void est()
  {
    for (iterator i = coefs.begin(); i != coefs.end(); ++i)
      (*i).second.est();
  }
};

class CalculatorWithCoefs;

class coef
{
  friend class CalculatorWithCoefs;
  static UnknownCoefs *u;
  static ComputedCoefs *c;
  static bool compiled;
  static void init()
  {
    if (!u)
    {
      u = new UnknownCoefs;
      c = new ComputedCoefs;
    }
  }

  class Destruct;
  friend class Destruct;
  class Destruct
  {
  public:
    ~Destruct()
    {
      if (u)
      {
        delete u;
        delete c;
        u = 0;
        c = 0;
      }
    }
  };
  static Destruct clean;

  double x_;
  double **ptr;
  bool computed;
  string id;

  static void WriteUnknown(ostream &out, const string &id);
  static void WriteComputed(ostream &out, const string &id);
  
public:
  static void WriteAll(ostream &out);
  static void WriteUnknowns(ostream &out);
  static void SetPrinted();
  static void SetNonPrinted();
  static void Unknowns(vec_double &x, vec_double &l, vec_double &u_, 
      vec_double &sc, vec_string &vs)
  {
    if (u)
      u->Unknowns(x, l, u_, sc, vs);
    else
      x.clear();
  }
  static void SetUnknowns(double *x)
  {
    if (u)
      u->SetUnknowns(x);
  }
  static void EstimateComputed();
  static void InitialValues(istream &in, const string &str = string());
  
  coef(const double& x = 0.) : x_(x), ptr(0), computed(false)
  {
    init();
  }
  istream& read(istream &in)
  {
    SGML e(in);
    read(e);
    return in;
  }
  void clear()
  {
    ptr = 0;
    x_ = 0.;
    id.erase();
    computed = false;
  }
  void read(const SGML &e);
  ostream& write(ostream &out) const;
  const double& x() const
  {
    if (ptr)
      return **ptr;
    else
      return x_;
  }
};

OVERLOAD_STREAMS(coef)

typedef vector<coef> vec_coef;
typedef vec_coef::iterator vec_coef_i;
typedef vec_coef::const_iterator vec_coef_ci;

class CalculatorWithCoefs : public calculator
{
  vec_coef v;
  string AnalizeSGML(const SGML &e)
  {
    coef c;
    c.read(e);
    if (!c.ptr)
      throw gError("CalculatorWithCoefs: anonymous coefs are not allowed");
    v.push_back(c);
    curr_tok.key = DOUBLE_PTR_PTR;
    curr_tok.x_ptr_ptr = c.ptr;
    return "xxx";
  }

public:
  CalculatorWithCoefs() {}
  CalculatorWithCoefs(const CalculatorWithCoefs &x)
  {
    FromString(x.body());
  }
  CalculatorWithCoefs& operator=(const CalculatorWithCoefs &x)
  {
    if (this != &x)
    {
      clear();
      FromString(x.body());
    }
    return *this;
  }
  virtual void clear()
  {
    v.clear();
    calculator::clear();
  } 
  
  ostream& write(ostream &out, size_t shift = 0) const;
  string body() const;
};

OVERLOAD_STREAMS(CalculatorWithCoefs)

#endif

