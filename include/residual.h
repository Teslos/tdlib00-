/*
Copyright (C) 1999 - 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#ifndef __RESIDUAL_H
#define __RESIDUAL_H

#include <td_algo.h>

class residual
{
  string id_;
  
  string yname;
  string xname;
  
  mutable int iy; //the ordinal number of the column for the y
  double sc; //the scale of x: x = x/scale

  compute comp;
  
public:
  residual() {clear();}

  void clear()
  {
    iy = -1;
    sc = 1.;
  }
  
  istream& read(istream& in)
  {
    parser p(in);
    SGML el;
    p.GetSGML(el);
    read(el);
    return in;
  }
  void read(const SGML &el);
  ostream& write(ostream& out, size_t shift = 0) const;

  const string& id() const
  {
    return id_;
  }
  void SetInput(const vec_string &vs) const
  {
    iy = SearchString(vs, yname);
    if (iy == -1)
      throw gError(string("residual: no yname - ").append(yname));
    comp.SetInput(vs);
  }
  void SetOnceInput(const vec_string &vs) const
  {
    comp.SetOnceInput(vs);
  }
  void SetOnce(double *reg) const
  {
    comp.SetOnce(reg);
  }
  double f(double *reg) const
  {
    return comp.f(reg);
  }
  double res(double *reg) const;
  double operator()(double *reg) const
  {
    return res(reg);
  }
  size_t NOfX(const vec_string &vs) const
  {
    return SearchString(vs, xname);
  }
  const string& NameOfX() const
  {
    return xname;
  }
  size_t NOfY() const
  {
    return iy;
  }
  const string& NameOfY() const
  {
    return yname;
  }
  double ScaleOfX() const
  {
    return sc;
  }
};

OVERLOAD_STREAMS(residual)

typedef map<string, residual> map_residual;
typedef map_residual::iterator map_residual_i;
typedef map_residual::const_iterator map_residual_ci;

extern map_residual res;
residual& FindResidual(const string &str);

#endif
