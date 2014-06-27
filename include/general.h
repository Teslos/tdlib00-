/*
Copyright (C) 1998, 1999, 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                       http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#ifndef __GENERAL_H
#define __GENERAL_H

#include <sstream>
#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <string>
#include <cmath>
#include <cstdlib>

using namespace std;

typedef vector<double> vec_double;
typedef vec_double::iterator vec_double_i;
typedef vec_double::const_iterator vec_double_ci;
typedef vector<double*> vec_double_ptr;
typedef vec_double_ptr::iterator vec_double_ptr_i;
typedef vec_double_ptr::const_iterator vec_double_ptr_ci;
typedef vector<double**> vec_double_2ptr;
typedef vec_double_2ptr::iterator vec_double_2ptr_i;
typedef vec_double_2ptr::const_iterator vec_double_2ptr_ci;
typedef vector<string> vec_string;
typedef vec_string::iterator vec_string_i;
typedef vec_string::const_iterator vec_string_ci;
typedef vector<int> vec_int;
typedef vec_int::iterator vec_int_i;
typedef vec_int::const_iterator vec_int_ci;
typedef map<string, string, less<string> > map_string_string;
typedef map_string_string::iterator map_string_string_i;
typedef map_string_string::const_iterator map_string_string_ci;

#define index index_

void WriteDebug(string &s); //from parser

//the idea of this class is taken from gnussl
class gError
{
public:
  string message;         // the error message

  gError() {WriteDebug(message);}
  gError(const char *msg) : message(msg) {WriteDebug(message);}
  gError(string msg) : message(msg) {WriteDebug(message);}
};

class CanNotCompute : public gError
{
public:
  double fmin;
  CanNotCompute(double f)
    : gError("The residual was not computed"), fmin(f) {}
};

template <class T, const char **name> class Enum
{
  size_t x;
public:
  Enum() : x(0) {}
  Enum(T xx) : x(xx) {}
  Enum(const string &s) {operator=(s);}
  operator T() const {return T(x);}
  Enum& operator=(const string &s)
  {
    size_t i = 0;
    while (name[i])
    {
      if (name[i][0] && s == name[i])
      {
        x = i;
        return *this;
      }
      ++i;
    }
    throw gError(string("Enum: can not convert string - ").append(s));
  }
  friend ostream& operator<<(ostream &out, Enum xx)
  {
    return out << name[xx.x];
  }
};

enum Function {G, H, S, Cp, V, dVdT, dVdp};
extern const char *FunctionName[];
typedef Enum<Function, FunctionName> function;

enum Index {excess = 1, ideal = 2, ref = 4, mix = 3, ref_ideal = 6, full = 7};
extern const char *IndexName[];
typedef Enum<Index, IndexName> index;

#define OVERLOAD_STREAMS(T)                     \
inline istream& operator>>(istream &in, T &to)  \
{                                               \
  in >> ws;                                     \
  return to.read(in);                           \
}                                               \
                                                            \
inline ostream& operator<<(ostream &out, const T &old)      \
{                                                           \
  return old.write(out);                                    \
}

// was not able to make template for above
// conflict with ostream::operator<<(char *)

class CheckInf
{
  double *xr;
public:
  CheckInf(double &xr_) : xr(&xr_) {}
  CheckInf(const double& xr_) : xr(&(const_cast<double&>(xr_))) {}

  istream& read(istream &in) const;
  ostream& write(ostream &out) const;
};

inline istream& operator>>(istream &in, const CheckInf &to)
{
  in >> ws;
  return to.read(in);
}
                                                            
inline ostream& operator<<(ostream &out, const CheckInf &old)
{
  return old.write(out);
}

class PutTab
{
public:
  size_t n;
  PutTab(size_t n) : n(n) {}
};

inline ostream& operator<<(ostream &out, const PutTab &old)
{
  for (size_t i = 0; i < old.n; ++i)
    out << "  ";
  return out;
}

template <class T>
string ObjToString(const T &x, int precision = 7)
{
  ostringstream out;
	out.precision(precision);
  out << x;
  return out.str();
}

struct limits
{
  double lower;
  double upper;
  limits() {}
  limits(double x, double y) {lower = x; upper = y;}
};

class global
{
public:
  static double R;     // R = 8.31441 by default
  static double T;     // T = 1000. by default
  static double p;     // p = 1. by default
  static double step1;   // 10.*sqrt(DBL_EPSILON) by default
  static double step2;   // 10.*pow(DBL_EPSILON, 1./3.) by default
  static bool neg_log;
  static double eps;
};

class StateTp
{
  double *Tptr;
  double *pptr;
  double T_;
  double p_;
public: 
  StateTp(const double &T = global::T, const double &p = global::p)
  {
    SetOwn();
    T_ = T;
    p_ = p;
  }
  StateTp(const StateTp &x)
  {
    SetOwn();
    T_ = x.T();
    p_ = x.p();
  }
  StateTp& operator=(const StateTp &x)
  {
    if (this != &x)
    {
      SetOwn();
      T_ = x.T();
      p_ = x.p();
    }
    return *this;
  }
  double& T()
  {
    return *Tptr;
  }
  const double& T() const
  {
    return *Tptr;
  }
  double& p()
  {
    return *pptr;
  }
  const double& p() const
  {
    return *pptr;
  }
  double*& ptrT()
  {
    return Tptr;
  }
  double* ptrT() const
  {
    return Tptr;
  }
  double*& ptrp()
  {
    return pptr;
  }
  double* ptrp() const
  {
    return pptr;
  }
  void SetOwn()
  {
    Tptr = &T_;
    pptr = &p_;
  }
  void SetOwnT(const double &T)
  {
    Tptr = &T_;
    T_ = T;
  }
  void SetOwnp(const double &p)
  {
    pptr = &p_;
    p_ = p;
  }
  void SetT(double &T)
  {
    T = *Tptr;
    Tptr = &T;
  }
  void Setp(double &p)
  {
    p = *pptr;
    pptr = &p;
  }
};

class StateX
{
  vec_double x_;
  vec_double_ptr xptr;
public:
  StateX() {}
  explicit StateX(size_t n) : x_(n), xptr(n)
  {
    SetOwn();
  }
  StateX(const StateX &x) : x_(x.size()), xptr(x.size())
  {
    SetOwn();
    x_ = x.x_;
  }
  StateX& operator=(const StateX &x)
  {
    if (this != &x)
    {
      if (x_.size() != x.size())
      {
        x_.resize(x.size());
        xptr.resize(x.size());
      }
      SetOwn();
      x_ = x.x_;
    }
    return *this;
  }
  double& operator[](size_t i)
  {
    return *(xptr[i]);
  }
  const double& operator[](size_t i) const
  {
    return *(xptr[i]);
  }
  double*& ptr(size_t i)
  {
    return xptr[i];
  }
  double* ptr(size_t i) const
  {
    return xptr[i];
  }
  size_t size() const
  {
    return x_.size();
  }
  void resize(size_t n)
  {
    x_.resize(n);
    xptr.resize(n);
    SetOwn();
  }
  void SetOwn()
  {
    for (size_t i = 0; i < x_.size() ; ++i)
      xptr[i] = &x_[i];
  }
  void SetX(size_t i, double &x)
  {
    x = *(xptr[i]);
    xptr[i] = &x;
  }
  void SetOwnX(size_t i, const double &x)
  {
    xptr[i] = &x_[i];
    x_[i] = x;
  }
};

template <class F>
inline double FirstDerivative(F f, const double &x, const double step = global::step1)
{
  double stepx = step*(x + step);
  double f1 = f(x - stepx);
  double f2 = f(x + stepx);
  return (f2 - f1)/stepx/2.;
}

template <class F>
inline double SecondDerivative(F f, const double &x, const double step = global::step2)
{
  double stepx = step*(x + step);
  double f1 = f(x - stepx);
  double f2 = f(x);
  double f3 = f(x + stepx);
  return ((f3 - f2) - (f2 - f1))/stepx/stepx;
}

template <class F>
inline double MixedDerivative(F f, const double &x, const double &y,
                      const double step = global::step2)
{
  double stepx = step*(x + step);
  double stepy = step*(y + step);
  double f1 = f(x - stepx, y - stepy);
  double f2 = f(x + stepx, y - stepy);
  double f3 = f(x - stepx, y + stepy);
  double f4 = f(x + stepx, y + stepy);
  return ((f4 - f3) - (f2 - f1))/stepx/stepy/4.;
}

class SGML
{
  static size_t dummy;
  static string dummystring;
public:
  
  string name;
  string body;
  map_string_string attr;

  SGML() {}
  SGML(istream &in)
  {
    read(in);
  }
  istream& read(istream &in, size_t &lf = dummy);
  ostream& write(ostream &out, size_t shift = 0) const;

  void compare(const string &s) const
  {
    if (!(name == s))
      throw gError(string("SGML: expected element ").append(s));
  }
  bool defined(const string &s) const
  {
    if (attr.find(s) != attr.end())
      return true;
    else
      return false;
  }
  const string& FindString(const string &x, const string &y = dummystring) const
  {
    map_string_string_ci i;
    if ((i = attr.find(x)) != attr.end())
      return (*i).second;
    else
      return y;
  }
  double FindDouble(const string &x, const double &y = 0.) const
  {
    map_string_string_ci i;
    if ((i = attr.find(x)) != attr.end())
      return atof((*i).second.c_str());
    else
      return y;
  }
  int FindInt(const string &x, int y = 0) const
  {
    map_string_string_ci i;
    if ((i = attr.find(x)) != attr.end())
      return atoi((*i).second.c_str());
    else
      return y;
  }
};

OVERLOAD_STREAMS(SGML)

class parser
{
  struct descriptor
  {
    string from;
    string from_id;
    string last;
    string last_id;
    size_t lf;
  };
  static list<descriptor> *debug;
  
  istream *in;
  bool del;
  list<descriptor>::iterator idebug;

  parser(const parser& old)
    {throw gError("parser: parser can not be copied");}
  parser& operator=(const parser& old)
    {throw gError("parser: parser can not be copied");}
  void init(const string &from, const string &from_id);
public:
  friend void WriteDebug(string &s);
  
  string token;
  string space;

  parser(istream &in_, const string &file = string()) : in(&in_), del(false) 
  {
    init("stream", file);
  }
  parser(const SGML &e)
  {
    init(e.name, e.FindString("id"));
    in = new istringstream(e.body);
    del = true;
  }
  ~parser();

  operator void*()
  {
    return *in;
  }

  string& GetToken();
  void compare(const string &s);
  void ReadCompare(const string &s)
  {
    GetToken();
    compare(s);
  }
// returns in->peek() 
  int SkipWS();
  bool eof()
  {
    return SkipWS() == EOF;
  }
  bool SkipChar(char t)
  {
    SkipWS();
    if (in->peek() == t)
    {
      in->ignore();
      return true;
    }
    return false;
  }
  bool SkipUntil(char t);
  size_t npos()
  {
    return in->tellg();
  }
  
  SGML& GetSGML(SGML &e);
  double GetDouble(const double &y = 0.)
  {
    double d = y;
    SkipWS();
    *in >> CheckInf(d);
    if (in->fail())
      in->clear();
    return d;
  }
  int GetInt(int y = 0)
  {
    int i = y;
    SkipWS();
    *in >> i;
    if (in->fail())
      in->clear();
    return i;
  }
  size_t lf() const
  {
    return (*idebug).lf;
  }
};

/* 
// old approximation
inline double log_z(const double &x)
{
  if (global::neg_log && x < global::eps)
    return log(global::eps) - 1. + x/global::eps;
  else
    return log(x);
}

//dlog_x/dx
inline double dlog_z(const double &x)
{
  if (global::neg_log && x < global::eps)
    return 1./global::eps;
  else
    return 1./x;
}
*/

inline double log_z(const double &x)
{
  if (global::neg_log && x < global::eps)
  {
    double div = x/global::eps;
    return log(global::eps) - 1.5 + div*(2. - 0.5*div);
  }
  else
    return log(x);
}

//dlog_x/dx
inline double dlog_z(const double &x)
{
  if (global::neg_log && x < global::eps)
    return (2. - x/global::eps)/global::eps;
  else
    return 1./x;
}

#endif

