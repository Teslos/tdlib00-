/*
Copyright (C) 1998, 1999, 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                       http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#include "interact.h"
#include <algorithm>

void RegisterInteraction()
{
  FuncTpx::RegisterType("Redlich_Kister", new RedlichKister);
  FuncTpx::RegisterType("RedlichKister", new RedlichKister);
  FuncTpx::RegisterType("Hoch_Arpshofen", new HochArpshofen);
  FuncTpx::RegisterType("HochArpshofen", new HochArpshofen);
  FuncTpx::RegisterType("Borelius_interaction", new Polynomial);
  FuncTpx::RegisterType("Polynomial", new Polynomial);
}

inline bool interaction::check(const StateX &x) const
{
  if (!p)
    return false;
/* the problem with Kohler - it is necessary to divide by SumX,
   on the other hand the derivative d2Sdx1dx2 from Ax1x2 is
   not zero
  SumX = 0.;
  for (size_t i = 0; i < m; ++i)
    SumX += x[vars[i]];
  if (!SumX)
    return false;
*/    
  return true;
}

inline void interaction::SetVec(const StateX &x) const
{
  switch (frm)
  {
    case ::NoChange:
    {
      for (size_t i = 0; i < m; ++i)
        vec[i] = x[vars[i]];
      return;
    }
    case ::Kohler:
    {
      SumX = 0.;
      for (size_t i = 0; i < m; ++i)
        SumX += x[vars[i]];
      if (!SumX)
        SumX = 1.; //to think it over
      for (size_t i = 0; i < m; ++i)
        vec[i] = x[vars[i]]/SumX;
      return;
    }
    case ::Muggianu:
    {
      SumX = 1.;
      for (size_t i = 0; i < m; ++i)
        SumX -= x[vars[i]];
      SumX /= double(m);
      for (size_t i = 0; i < m; ++i)
        vec[i] = x[vars[i]] + SumX;
      return;
    }
  }
}

inline void interaction::Setdvdx(const StateX &x) const // after SetVec
{
  switch (frm)
  {
    case ::NoChange:
      dvdx_[0] = 0.;
      dvdx_[1] = 1.;
      return;
    case ::Kohler:
    {
      for (size_t i = 0; i < m; ++i)
      {
        dvdx_[2*i] = -vec[i]/SumX;
        dvdx_[2*i + 1] = dvdx_[2*i] + 1/SumX;
      }
      return;
    }
    case ::Muggianu:
    {
      dvdx_[0] = -1./double(m);
      dvdx_[1] = dvdx_[0] + 1.;
      return;
    }
  }
}

inline const double& interaction::dvdx(size_t i, size_t j) const
// dvi/dxj
{
  switch (frm)
  {
    case ::Kohler:      // fi = xi/Sum
      return (i == j) ? dvdx_[2*i + 1] : dvdx_[2*i];
    case ::NoChange:   // fi = xi
    case ::Muggianu:    // fi = xi + (1 - Sum)/m
      return (i == j) ? dvdx_[1] : dvdx_[0];
    default:
      throw gError(string("interaction : not known formalism"));
  }
}

inline void interaction::Setd2vdx2(const StateX &x) const // after Setdvdx
{
  switch (frm)
  {
    case ::NoChange:
    case ::Muggianu:
      d2vdx2_[0] = 0.;
      return;
    case ::Kohler:
    {
      for (size_t i = 0; i < m; ++i)
      {
        d2vdx2_[3*i] = -2.*dvdx_[2*i]/SumX;
        d2vdx2_[3*i + 1] = d2vdx2_[3*i] - 1/SumX/SumX;
        d2vdx2_[3*i + 2] = d2vdx2_[3*i + 1] - 1/SumX/SumX;
      }
      return;
    }
  }
}

inline const double& interaction::d2vdx2(size_t i, size_t j, size_t k) const
// d2vi/dxj/dxk
{
  switch (frm)
  {
    case ::Kohler:
      if (i == j && i == k)
        return d2vdx2_[i*3 + 2];
      else if (i == j || i == k)
        return d2vdx2_[i*3 + 1];
      else
        return d2vdx2_[i*3];
    case ::NoChange: 
    case ::Muggianu:
      return d2vdx2_[0];
    default:
      throw gError(string("interaction : not known formalism"));
  }
}

inline void RedlichKister::SetSum(function fun, const StateTp &Tp, 
    const StateX &x) const
{
  Sum = f[0].Z(fun, Tp);
  if (p > 1)
  {
    SetVec(x);
    for (int i = 1; i < p; ++i)
      Sum += f[i].Z(fun, Tp)*pow(vec[0] - vec[1], i);
  }
}

inline void RedlichKister::SetdSdx(function fun, const StateTp &Tp, 
    const StateX &x) const
{
  dSdx[0] = dSdx[1] = dSdv[0] = dSdv[1] = 0.;
  if (p > 1)
  {
    dSdv[0] = f[1].Z(fun, Tp);
    for (int i = 2; i < p; ++i)
      dSdv[0] += double(i)*f[i].Z(fun, Tp)*pow(vec[0] - vec[1], i - 1);
    dSdv[1] = -dSdv[0];
    switch (frm)
    {
      case ::NoChange:
        dSdx[0] = dSdv[0];
        dSdx[1] = dSdv[1];
        return;
      case ::Kohler:
      case ::Muggianu:
        Setdvdx(x);
        dSdx[0] = dSdv[0]*(dvdx(0, 0) - dvdx(1, 0));
        dSdx[1] = dSdv[0]*(dvdx(0, 1) - dvdx(1, 1));
        return;
    }
  }
}

inline void RedlichKister::Setd2Sdx2(function fun, const StateTp &Tp, 
    const StateX &x) const
{
  d2Sdv2(0, 0) = d2Sdv2(0, 1) = d2Sdv2(1, 1) = 0.;
  d2Sdx2(0, 0) = d2Sdx2(0, 1) = d2Sdx2(1, 1) = 0.;
  if (p > 2)
  {
    d2Sdv2(0, 0) = 2.*f[2].Z(fun, Tp);
    for (int i = 3; i < p; ++i)
      d2Sdv2(0, 0) += double(i*(i - 1))
        *f[i].Z(fun, Tp)*pow(vec[0] - vec[1], i - 2);
    d2Sdv2(1, 1) = d2Sdv2(0, 0);
    d2Sdv2(0, 1) = -d2Sdv2(0, 0);
  }
  if (p > 1)
  {
    Setd2vdx2(x);
    switch (frm)
    {
      case ::NoChange:
        d2Sdx2(0, 0) = d2Sdv2(0, 0);
        d2Sdx2(0, 1) = d2Sdv2(0, 1);
        d2Sdx2(1, 1) = d2Sdv2(1, 1);
        return;
      case ::Kohler:
        d2Sdx2(0, 0) = dSdv[0]*d2vdx2(0, 0, 0) + dSdv[1]*d2vdx2(1, 0, 0);
        d2Sdx2(0, 1) = dSdv[0]*d2vdx2(0, 0, 1) + dSdv[1]*d2vdx2(1, 0, 1);
        d2Sdx2(1, 1) = dSdv[0]*d2vdx2(0, 1, 1) + dSdv[1]*d2vdx2(1, 1, 1);
      case ::Muggianu:
        d2Sdx2(0, 0) += d2Sdv2(0, 0)
          *(pow(dvdx(0, 0),2) - 2.*dvdx(0, 0)*dvdx(1, 0) + pow(dvdx(1, 0),2));
        d2Sdx2(1, 1) += d2Sdv2(0, 0)
          *(pow(dvdx(0, 1),2) - 2.*dvdx(0, 1)*dvdx(1, 1) + pow(dvdx(1, 1),2));
        d2Sdx2(0, 1) += d2Sdv2(0, 0)
          *(dvdx(0, 0)*dvdx(0, 1) - dvdx(0, 0)*dvdx(1, 1) 
              - dvdx(0, 1)*dvdx(1, 0) + dvdx(1, 0)*dvdx(1, 1));
      return;
    }
  }
}

void RedlichKister::read(const SGML &e, const vec_formula *fm)
{
  clear();
  ReadFormalism(e);
  old = frm;
  if (frm == Muggianu)
    frm = NoChange;
  parser par(e);
  par.SkipChar('+');
  vars[0] = FindFormula(par, fm);
  par.SkipChar('*');
  vars[1] = FindFormula(par, fm);
  par.SkipChar('*');
  par.SkipUntil('<');
  if (vars[0] == vars[1])
    throw gError("RedlichKister - wrong input");
  SGML el;
  int i;
  do
  {
    par.GetSGML(el);
    el.compare("func_x");
    parser par2(el);
    par2.SkipChar('+');
    par2.SkipChar('1');
    par2.SkipChar('*');
    if (par2.eof())
      i = 0;
    else
    {
      if (par2.SkipUntil('^'))
      {
        par2.SkipChar('^');
        i = par2.GetInt();
        if (i < 0)
          throw gError("RedlichKister - wrong input");
      }
      else
        i = 1;
    }
    if (i + 1 > f.size())
      f.insert(f.end(), i - f.size() + 1, func_Tp());
    f[i].read(par.GetSGML(el));
  }
  while (!par.SkipChar(')') && !par.eof());
  p = f.size();
}

void RedlichKister::WriteBody(ostream &out, const vec_formula *fm, 
      size_t shift) const
{
  size_t i;
  out << PutTab(shift) << '+';
  for (i = 0; i != vars.size(); ++i)
  {
    out << "x(";
    WriteFormula(out, i, fm);
    out << ")*";
  }
  out << '(' << endl;
  char nm = 'v';
  if (frm == NoChange)
    nm = 'x';
  for (i = 0; i < f.size(); ++i)
  {
    out << PutTab(shift + 1) << "<func_x> +";
    if (i != 0)
    {
      out << '(' << nm << '(';
      WriteFormula(out, 0, fm);
      out << ")-" << nm << '(';
      WriteFormula(out, 1, fm);
      out << "))";
      if (i == 1)
        out << '*';
      else
        out << '^' << i << '*';
    }
    out << " </func_x> " << endl;
    f[i].write(out, shift + 1);
  }
  out << PutTab(shift) << ")" << endl;
}

double RedlichKister::Z(function fun, const StateTp &Tp,
                        const StateX &x) const
{
  if (!check(x))
    return 0.;
  SetSum(fun, Tp, x);
  return x[vars[0]]*x[vars[1]]*Sum;
}

void RedlichKister::dZdx(function fun, const StateTp &Tp,
                         const StateX &x, vec_double& res) const
{
  if (!check(x))
    return;
  double P = x[vars[0]]*x[vars[1]];
  SetSum(fun, Tp, x); //SetVec
  SetdSdx(fun, Tp, x); // Setdvdx
  res[vars[0]] += x[vars[1]]*Sum + P*dSdx[0];
  res[vars[1]] += x[vars[0]]*Sum + P*dSdx[1];
}

void RedlichKister::d2Zdx2(function fun, const StateTp &Tp,
                           const StateX &x, SymMatrix& res) const
{
  if (!check(x))
    return;
  double P = x[vars[0]]*x[vars[1]];
  SetSum(fun, Tp, x); // SetVec
  SetdSdx(fun, Tp, x); // Setdvdx 
  Setd2Sdx2(fun, Tp, x); // Setdvdx
  res(vars[0], vars[0]) += 2.*x[vars[1]]*dSdx[0] + P*d2Sdx2(0, 0);
  res(vars[0], vars[1]) += Sum + x[vars[1]]*dSdx[1] 
    + x[vars[0]]*dSdx[0] + P*d2Sdx2(0, 1);
  res(vars[1], vars[1]) += 2.*x[vars[0]]*dSdx[1] + P*d2Sdx2(1, 1);
}

inline void HochArpshofen::SetSum(function fun, const StateTp &Tp, 
    const StateX &x) const
{
  Sum = 0;
  for (size_t i = 0; i < p; ++i)
    if (pwr[i] > 0.)
      Sum += f[i].Z(fun, Tp)*(vec[0] - pow(vec[0], pwr[i]));
    else
      Sum += f[i].Z(fun, Tp)*(vec[1] - pow(vec[1], -pwr[i]));
}

inline void HochArpshofen::SetdSdx(function fun, const StateTp &Tp, 
    const StateX &x) const
{
  dSdx[0] = dSdx[1] = dSdv[0] = dSdv[1] = 0.;
  for (size_t i = 0; i < p; ++i)
    if (pwr[i] > 0.)
      dSdv[0] += f[i].Z(fun, Tp)*(1. - pwr[i]*pow(vec[0], pwr[i] - 1));
    else
      dSdv[1] += f[i].Z(fun, Tp)*(1. + pwr[i]*pow(vec[1], -pwr[i] - 1));
  switch (frm)
  {
    case ::NoChange:
      dSdx[0] = dSdv[0];
      dSdx[1] = dSdv[1];
      return;
    case ::Kohler:
    case ::Muggianu:
      Setdvdx(x);
      dSdx[0] = dSdv[0]*dvdx(0, 0) + dSdv[1]*dvdx(1, 0);
      dSdx[1] = dSdv[0]*dvdx(0, 1) + dSdv[1]*dvdx(1, 1);
      return;
  }
}

inline void HochArpshofen::Setd2Sdx2(function fun, const StateTp &Tp, 
    const StateX &x) const
{
  d2Sdv2(0, 0) = d2Sdv2(0, 1) = d2Sdv2(1, 1) = 0.;
  d2Sdx2(0, 0) = d2Sdx2(0, 1) = d2Sdx2(1, 1) = 0.;
  for (size_t i = 0; i < p; ++i)
    if (pwr[i] > 0.)
      d2Sdv2(0, 0) += f[i].Z(fun, Tp)
        *(-pwr[i]*(pwr[i] - 1)*pow(vec[0], pwr[i] - 2));
    else
      d2Sdv2(1, 1) += f[i].Z(fun, Tp)
        *(+pwr[i]*(-pwr[i] - 1)*pow(vec[1], -pwr[i] - 2));
  Setd2vdx2(x);
  switch (frm)
  {
    case ::NoChange:
      d2Sdx2(0, 0) = d2Sdv2(0, 0);
      d2Sdx2(0, 1) = d2Sdv2(0, 1);
      d2Sdx2(1, 1) = d2Sdv2(1, 1);
      return;
    case ::Kohler:
      d2Sdx2(0, 0) = dSdv[0]*d2vdx2(0, 0, 0) + dSdv[1]*d2vdx2(1, 0, 0);
      d2Sdx2(0, 1) = dSdv[0]*d2vdx2(0, 0, 1) + dSdv[1]*d2vdx2(1, 0, 1);
      d2Sdx2(1, 1) = dSdv[0]*d2vdx2(0, 1, 1) + dSdv[1]*d2vdx2(1, 1, 1);
    case ::Muggianu:
      d2Sdx2(0, 0) += d2Sdv2(0, 0)*pow(dvdx(0, 0),2) 
        + d2Sdv2(1, 1)*pow(dvdx(1, 0),2);
      d2Sdx2(1, 1) += d2Sdv2(0, 0)*pow(dvdx(0, 1),2) 
        + d2Sdv2(1, 1)*pow(dvdx(1, 1),2);
      d2Sdx2(0, 1) += d2Sdv2(0, 0)*dvdx(0, 0)*dvdx(0, 1) 
        + d2Sdv2(1, 1)*dvdx(1, 0)*dvdx(1, 1);
      return;
  }
}

void HochArpshofen::read(const SGML &e, const vec_formula *fm)
{
  clear();
  ReadFormalism(e);
  parser par(e);
  par.SkipChar('+');
  vars[0] = FindFormula(par, fm);
  par.SkipChar('*');
  vars[1] = FindFormula(par, fm);
  par.SkipChar('*');
  par.SkipUntil('<');
  sort(vars.begin(), vars.end());
  if (vars[0] == vars[1])
    throw gError("HochArpshofen - wrong input");
  SGML el;
  func_Tp ftp;
  do
  {
    par.GetSGML(el);
    el.compare("func_x");
    int sign, num, i = 0;
    parser par2(el);
    par2.SkipChar('+');
    par2.SkipChar('(');
    num = FindFormula(par2, fm);
    if (num == vars[0])
      sign = 1;
    else if (num == vars[1])
      sign = -1;
    else
      throw gError("HochArpshofen - wrong input");
    par2.SkipUntil('^');
    par2.SkipChar('^');
    i = par2.GetInt();
    if (i < 2)
      throw gError("HochArpshofen - wrong input");
    ftp.read(par.GetSGML(el));
    size_t j;
    for (j = 0; j < pwr.size(); ++j)
      if (i < pwr[j])
        break;
    pwr.insert(pwr.begin() + j, sign*i);
    f.insert(f.begin() + j, ftp);
  }
  while (!par.SkipChar(')') && !par.eof());
  p = f.size();
}

void HochArpshofen::WriteBody(ostream &out, const vec_formula *fm, 
      size_t shift) const
{
  size_t i;
  out << PutTab(shift) << '+';
  out << "x(";
  WriteFormula(out, 0, fm);
  out << ")*x(";
  WriteFormula(out, 1, fm);
  out << ")/v(";
  WriteFormula(out, 0, fm);
  out << ")/v(";
  WriteFormula(out, 1, fm);
  out << ")*(" << endl;
  for (i = 0; i < f.size(); ++i)
  {
    out << PutTab(shift + 1) << "<func_x> +";
    if (pwr[i] > 0.)
    {
      out << "(v(";
      WriteFormula(out, 0, fm);
      out << ")-v(";
      WriteFormula(out, 0, fm);
      out << ")^" << pwr[i] << ")*";
    }
    else
    {
      out << "(v(";
      WriteFormula(out, 1, fm);
      out << ")-v(";
      WriteFormula(out, 1, fm);
      out << ")^" << -pwr[i] << ")*";
    }
    out << " </func_x> " << endl;
    f[i].write(out, shift + 1);
  }
  out << PutTab(shift) << ")" << endl;
}

double HochArpshofen::Z(function fun, const StateTp &Tp,
                         const StateX &x) const
{
  if (!check(x))
    return 0.;
  SetVec(x);
  SetSum(fun, Tp, x);
  if (frm != ::NoChange && vec[0] != 0. && vec[1] != 0.)
    Sum = Sum/vec[0]/vec[1]*x[vars[0]]*x[vars[1]];
  return Sum;
}

void HochArpshofen::dZdx(function fun, const StateTp &Tp,
                          const StateX &x, vec_double& res) const
{
  if (!check(x))
    return;
  SetVec(x);
  SetdSdx(fun, Tp, x);
  if (frm == ::NoChange)
  {
    res[vars[0]] += dSdx[0];
    res[vars[1]] += dSdx[1];
  }
  else
  {
    double P = x[vars[0]]*x[vars[1]];
    double P1 = vec[0]*vec[1];
    if (P1 == 0.)
      return;
    SetSum(fun, Tp, x);
    res[vars[0]] += P/P1*dSdx[0] + Sum*(x[vars[1]]/P1 
        - P/P1/vec[0]*dvdx(0, 0) - P/P1/vec[1]*dvdx(1, 0));
    res[vars[1]] += P/P1*dSdx[1] + Sum*(x[vars[0]]/P1 
        - P/P1/vec[0]*dvdx(0, 1) - P/P1/vec[1]*dvdx(1, 1));
  }
}

void HochArpshofen::d2Zdx2(function fun, const StateTp &Tp,
                           const StateX &x, SymMatrix& res) const
{
  if (!check(x))
    return;
  SetVec(x);
  SetdSdx(fun, Tp, x);
  Setd2Sdx2(fun, Tp, x); 
  if (frm == ::NoChange)
  {
    res(vars[0],vars[0]) += d2Sdx2(0, 0);
    res(vars[0],vars[1]) += d2Sdx2(0, 1);
    res(vars[1],vars[1]) += d2Sdx2(1, 1);
  }
  else
  {
    double P = x[vars[0]]*x[vars[1]];
    double P1 = vec[0]*vec[1];
    if (P1 == 0.)
      return;
    SetSum(fun, Tp, x);
    double S1 = x[vars[1]]/P1 - P/P1/vec[0]*dvdx(0, 0) - P/P1/vec[1]*dvdx(1, 0);
    double S2 = x[vars[0]]/P1 - P/P1/vec[0]*dvdx(0, 1) - P/P1/vec[1]*dvdx(1, 1);
//    res[vars[0]] += P/P1*dSdx[0] + Sum*S1; first derivatives
    res(vars[0],vars[0]) += P/P1*d2Sdx2(0, 0) + 2.*dSdx[0]*S1 + Sum*(
// first member
        - x[vars[1]]/P1/vec[0]*dvdx(0, 0) - x[vars[1]]/P1/vec[1]*dvdx(1, 0)
// second member
        - x[vars[1]]/P1/vec[0]*dvdx(0, 0) - P/P1/vec[0]*d2vdx2(0, 0, 0)
        + 2.*P/P1/vec[0]/vec[0]*dvdx(0, 0)*dvdx(0, 0)
        + P/P1/vec[0]/vec[1]*dvdx(0, 0)*dvdx(1, 0)
// third member
        - x[vars[1]]/P1/vec[1]*dvdx(1, 0) - P/P1/vec[1]*d2vdx2(1, 0, 0)
        + P/P1/vec[1]/vec[0]*dvdx(1, 0)*dvdx(0, 0)
        + 2.*P/P1/vec[1]/vec[1]*dvdx(1, 0)*dvdx(1, 0)
        );
    res(vars[0],vars[1]) += P/P1*d2Sdx2(0, 1) + dSdx[0]*S2 + dSdx[1]*S1 + Sum*(
// first member
        1/P1 - x[vars[1]]/P1/vec[0]*dvdx(0, 1) - x[vars[1]]/P1/vec[1]*dvdx(1, 1)
// second member
        - x[vars[0]]/P1/vec[0]*dvdx(0, 0) - P/P1/vec[0]*d2vdx2(0, 0, 1)
        + 2.*P/P1/vec[0]/vec[0]*dvdx(0, 0)*dvdx(0, 1)
        + P/P1/vec[0]/vec[1]*dvdx(0, 0)*dvdx(1, 1)
// third member
        - x[vars[0]]/P1/vec[1]*dvdx(1, 0) - P/P1/vec[1]*d2vdx2(1, 0, 1)
        + P/P1/vec[1]/vec[0]*dvdx(1, 0)*dvdx(0, 1)
        + 2.*P/P1/vec[1]/vec[1]*dvdx(1, 0)*dvdx(1, 1)
        );
//    res[vars[1]] += P/P1*dSdx[1] + Sum*S2; first derivatives
    res(vars[1],vars[1]) += P/P1*d2Sdx2(1, 1) + 2.*dSdx[1]*S2 + Sum*(
// first member
        - x[vars[0]]/P1/vec[0]*dvdx(0, 1) - x[vars[0]]/P1/vec[1]*dvdx(1, 1)
// second member
        - x[vars[0]]/P1/vec[0]*dvdx(0, 1) - P/P1/vec[0]*d2vdx2(0, 1, 1)
        + 2.*P/P1/vec[0]/vec[0]*dvdx(0, 1)*dvdx(0, 1)
        + P/P1/vec[0]/vec[1]*dvdx(0, 1)*dvdx(1, 1)
// third member
        - x[vars[0]]/P1/vec[1]*dvdx(1, 1) - P/P1/vec[1]*d2vdx2(1, 1, 1)
        + P/P1/vec[1]/vec[0]*dvdx(1, 1)*dvdx(0, 1)
        + 2.*P/P1/vec[1]/vec[1]*dvdx(1, 1)*dvdx(1, 1)
        );
  }
}

inline void Polynomial::SetSum(function fun, const StateTp &Tp, 
    const StateX &x) const
{
  SetVec(x);
  Sum = 0.;
  size_t ipwr = 0;
  for (vec_func_Tp_ci ip = f.begin(); ip != f.end(); ++ip)
  {
    double P1 = (*ip).Z(fun, Tp);
    for (size_t i = 0; i < m; ++i, ++ipwr)
      if (pwr[ipwr])
        P1 *= pow(vec[i], pwr[ipwr]);
    Sum += P1;
  }
}

inline void Polynomial::SetdSdx(function fun, const StateTp &Tp, 
    const StateX &x) const
{
  fill(dSdv.begin(), dSdv.end(), 0.);   
  size_t i;
  size_t ipwr = 0;
  for (vec_func_Tp_ci ip = f.begin(); ip != f.end(); ++ip)
  {
    // used as a temporary
    fill(dSdx.begin(), dSdx.end(), (*ip).Z(fun, Tp));   
    double mul1, mul2;
    for (i = 0; i < m; ++i, ++ipwr)
    {
      if (pwr[ipwr] == 0)
      {
        mul1 = 0.;
        mul2 = 1.;
      }
      else if (pwr[ipwr] == 1)
      {
        mul1 = 1.;
        mul2 = vec[i];
      }
      else if (pwr[ipwr] == 2)
      {
        mul1 = 2.*vec[i];
        mul2 = vec[i]*vec[i];
      }
      else
      {
        mul1 = pow(vec[i], pwr[ipwr] - 1);
        mul2 = mul1*vec[i];
        mul1 *= pwr[ipwr];
      }
      for (size_t j = 0; j < m; ++j)
        if (j != i)
          dSdx[j] *= mul2;
        else
          dSdx[j] *= mul1;
    }
    for (i = 0; i < m; ++i)
      dSdv[i] += dSdx[i];
  }
  switch (frm)
  {
    case ::NoChange:
      for (i = 0; i < m; ++i)
        dSdx[i] = dSdv[i];
      return;
    case ::Kohler:
    case ::Muggianu:
      Setdvdx(x);
      fill(dSdx.begin(), dSdx.end(), 0.);   
      for (i = 0; i < m; ++i)
        for (size_t j = 0; j < m; ++j)
          dSdx[i] += dSdv[j]*dvdx(j, i);
      return;
  }
}

inline void Polynomial::Setd2Sdx2(function fun, const StateTp &Tp, 
    const StateX &x) const
{
  fill(d2Sdv2.begin(), d2Sdv2.end(), 0.);
  size_t i;
  size_t ipwr = 0;
  for (vec_func_Tp_ci ip = f.begin(); ip != f.end(); ++ip)
  {
    // used as a temporary
    fill(d2Sdx2.begin(), d2Sdx2.end(), (*ip).Z(fun, Tp));
    double mul1, mul2, mul3;
    for (i = 0; i < m; ++i, ++ipwr)
    {
      if (pwr[ipwr] == 0)
      {
        mul1 = 0.;
        mul2 = 0.;
        mul3 = 1.;
      }
      else if (pwr[ipwr] == 1)
      {
        mul1 = 0.;
        mul2 = 1.;
        mul3 = vec[i];
      }
      else if (pwr[ipwr] == 2)
      {
        mul1 = 2.;
        mul2 = 2.*vec[i];
        mul3 = vec[i]*vec[i];
      }
      else
      {
        mul1 = pow(vec[i], pwr[ipwr] - 2);
        mul2 = mul1*vec[i];
        mul3 = mul2*vec[i];
        mul2 *= pwr[ipwr];
        mul1 *= pwr[ipwr]*(pwr[ipwr] - 1);
      }
      for (size_t j = 0; j < m; ++j)
        for (size_t k = 0; k <= j; ++k)
          if (i == j && i == k)
            d2Sdx2(j, k) *= mul1;
          else if (i == j || i == k)
            d2Sdx2(j, k) *= mul2;
          else
            d2Sdx2(j, k) *= mul3;
    }
    for (i = 0; i < d2Sdx2.size(); ++i)
      d2Sdv2[i] += d2Sdx2[i];
  }
  Setd2vdx2(x);
  fill(d2Sdx2.begin(), d2Sdx2.end(), 0.);
  switch (frm)
  {
    case ::NoChange:
      for (i = 0; i < d2Sdx2.size(); ++i)
        d2Sdx2[i] = d2Sdv2[i];
      return;
    case ::Kohler:
      for (i = 0; i < m; ++i)
        for (size_t j = 0; j <= i; ++j)
          for (size_t k = 0; k < m; ++k)
            d2Sdx2(i, j) += dSdv[k]*d2vdx2(k, i, j);
    case ::Muggianu:
      for (i = 0; i < m; ++i)
        for (size_t j = 0; j <= i; ++j)
          for (size_t k = 0; k < m; ++k)
            for (size_t l = 0; l < m; ++l)
              d2Sdx2(i, j) += d2Sdv2(k, l)*dvdx(k, i)*dvdx(l, j);
      return;
  }
}

Polynomial::iterator Polynomial::begin()
{
  iterator t;
  if (m)
  {
    t.indx.push_back(power);
    t.indx.insert(t.indx.end(), m - 1, 0);
  }
  return t;
}

Polynomial::iterator Polynomial::end()
{
  iterator t;
  if (m)
  {
    t.indx.insert(t.indx.end(), m - 1, 0);
    t.indx.push_back(power);
  }
  return t;
}

Polynomial::iterator& Polynomial::iterator::operator++()
{
  for (size_t j = 0; j < indx.size() - 1; ++j)
    if (indx[j])
    {
      if (j == 0)
      {
        --indx[0];
        ++indx[1];
      }
      else
      {
        --indx[j];
        ++indx[j+1];
        indx[0] = indx[j];
        indx[j] = 0;
      }
      return *this;
    }
  return *this;
}

Polynomial::iterator Polynomial::iterator::operator++(int)
{
  iterator tmp = *this;
  ++*this;
  return tmp;
}

void Polynomial::read(const SGML &e, const vec_formula *fm)
{
  istringstream in(e.body);
  parser par(in);
  par.SkipChar('+');
  string str;
  char ch;
  while (in.get(ch), in && in.peek() != '<')
    str.append(1, ch);
  str.erase(str.find_last_of('('));
  size_t i, ix, pi;
  istringstream in2(str);
  parser par2(in2);
  while (!par2.eof())
  {
    ix = FindFormula(par2, fm);
    vars.push_back(ix);
    par2.SkipChar('*');
  }
  sort(vars.begin(), vars.end());
  if (adjacent_find(vars.begin(), vars.end()) != vars.end())
    throw gError("Polynomial - wrong input of mf");
  clear(vars.size());
  vec_int q(vars.size());
  SGML el;
  func_Tp ftp;
  do
  {
    par.GetSGML(el);
    el.compare("func_x");
    parser par2(el);
    par2.SkipChar('+');
    par2.SkipChar('1');
    par2.SkipChar('*');
    fill(q.begin(), q.end(), 0);
    while (!par2.eof())
    {
      ix = FindFormula(par2, fm);
      if (par2.SkipChar('^'))
        pi = par2.GetInt();
      else
        pi = 1;
      for (i = 0; i != vars.size(); ++i)
        if (ix == vars[i]) 
        {
          q[i] = pi;
          break;
        }
      par2.SkipChar('*');
    }
    ftp.read(par.GetSGML(el));
    size_t j;
    for (j = 0; j < pwr.size(); j += m)
    {
      if (lexicographical_compare(pwr.begin() + j, pwr.begin() + j + m,
          q.begin(), q.end()))
        break;
    }
    pwr.insert(pwr.begin() + j, q.begin(), q.end());
    f.insert(f.begin() + j/m, ftp);
  }
  while (!par.SkipChar(')') && !par.eof());
  p = f.size();
  ReadFormalism(e);
}

void Polynomial::WriteBody(ostream &out, const vec_formula *fm, 
    size_t shift) const
{
  if (m == 0) return;
  size_t i;
  out << PutTab(shift) << '+';
  for (i = 0; i != m; ++i)
  {
    out << "x(";
    WriteFormula(out, i, fm);
    out << ")*";
  }
  out << '(' << endl;
  size_t ipwr = 0;
  for (vec_func_Tp_ci j = f.begin(); j != f.end(); ++j)
  {
    out << PutTab(shift + 1) << "<func_x> +";
    for (i = 0; i != m; ++i, ++ipwr)
      if (pwr[ipwr] == 1)
      {
        out << "v(";
        WriteFormula(out, i, fm);
        out << ")*";
      }
      else if (pwr[ipwr] > 1)
      {
        out << "v(";
        WriteFormula(out, i, fm);
        out << ")^" << pwr[ipwr] << '*';
      }
    out << " </func_x> " << endl;
    (*j).write(out, shift + 1);
  }
  out << PutTab(shift) << ")" << endl;
}

double Polynomial::Z(function fun, const StateTp &Tp, const StateX &x) const
{
  if (!check(x))
    return 0.;
  double P = 1.;
  for (size_t i = 0; i < m; ++i)
    P *= x[vars[i]];
  SetSum(fun, Tp, x); //SetVec
  return P*Sum;
}

void Polynomial::dZdx(function fun, const StateTp &Tp, const StateX &x, 
    vec_double& res) const
{
  if (!check(x))
    return;
  size_t i;
  double P = 1.;
  for (i = 0; i < m; ++i)
    P *= x[vars[i]];
  SetSum(fun, Tp, x);
  SetdSdx(fun, Tp, x);
  for (i = 0; i < m; ++i)
  {
    double Pi = 1.;
    for (size_t j = 0; j < m; ++j)
      if (j != i)
        Pi *= x[vars[j]];
    res[vars[i]] += Pi*Sum + P*dSdx[i];
  }
}

void Polynomial::d2Zdx2(function fun, const StateTp &Tp, const StateX &x, 
    SymMatrix& res) const
{
  if (!check(x))
    return;
  size_t i;
  double P = 1.;
  for (i = 0; i < m; ++i)
    P *= x[vars[i]];
  SetSum(fun, Tp, x);
  SetdSdx(fun, Tp, x);
  Setd2Sdx2(fun, Tp, x);
  vec_double dPdx(m);
  for (i = 0; i < m; ++i)
  {
    dPdx[i] = 1.;
    for (size_t k = 0; k < m; ++k)
      if (k != i)
        dPdx[i] *= x[vars[k]];
  }
  for (i = 0; i < m; ++i)
    for (size_t j = 0; j <= i; ++j)
    {
      double P2ij;
      if (i == j)
        P2ij = 0.;
      else
      {
        P2ij = 1.;
        for (size_t k = 0; k < m; ++k)
          if (k != i && k != j )
            P2ij *= x[vars[k]];
      }
      res(vars[i], vars[j]) += P2ij*Sum + dPdx[i]*dSdx[j] + dPdx[j]*dSdx[i] + 
          P*d2Sdx2(i, j);
    }
}

