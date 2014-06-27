/*
Copyright (C) 1998, 1999, 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                       http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#include <iomanip>
#include <toms.h>
#include "cuox_pl.h"

void RegisterCuOx_plane()
{
  phase::RegisterType("CuOx_plane", new CuOx_plane);
}

void CuOx_plane::read(const SGML &e_)
{
  clear();
  size_t i, k;
  SGML e = e_;
  ReadComponents(vf, e, ref);
  parser p(e);
  SGML el;
  MolecularFormula_ci i0 = vf[0].mf().begin();
  MolecularFormula_ci i1 = vf[1].mf().begin();
  elem O("O");
  while (i1 != vf[1].mf().end())
  {
    if ((*i1).first == O)
    {
      if ((*i0).first == O)
      {
        if ((*i1++).second - (*i0++).second != 1)
          throw gError("CuOx_plane: components are incompatible");
      }
      else
      {
        if ((*i1++).second != 1)
          throw gError("CuOx_plane: components are incompatible");
      }
    }
    else
    {
      if (*i1++ != *i0++)
        throw gError("CuOx_plane: components are incompatible");
    }
  }
  if (p.GetToken() == "Gid")
    old_model = true;
  else if (p.token == "R*T*")
  {
    p.SkipChar('(');
    if (p.SkipChar('('))
      old_model = false;
    else
      old_model = true;
  }
  else
    throw gError("CuOx_plane: expected ideal Gibbs energy");
  p.SkipUntil('<');
  while (!p.eof())
  {
    p.GetSGML(el);
    el.compare("func_x");
    parser p2(el);
    p2.SkipChar('+');
    p2.SkipChar('1');
    p2.SkipChar('*');
    p2.GetToken();
    if (p2.token == "z*")
    {
      if (!p2.SkipUntil('^'))
        i = 1;
      else
      {
        p2.SkipChar('^');
        k = p2.GetInt();
        if (k == 0)
          throw gError("CuOx_plane: read: wrong token");
        i = k + 1;
        if (k > na_)
        {
          f.insert(f.begin() + 2 + na_, k - na_, func_Tp());
          na_ = k;
        }
      }
    }
    else if (p2.token == "(")
    {
      if (!p2.SkipUntil('z'))
      {
        i = 2 + na_;
        if (!nb_)
        {
          f.insert(f.begin() + 2 + na_, func_Tp());
          nb_++;
        }
      }
      else 
      {
        p2.SkipUntil('^');
        p2.SkipChar('^');
        k = p2.GetInt();
        if (k == 0)
          throw gError("CuOx_plane: read: wrong token");
        i = k + 1 + na_;
        if (k > nb_)
        {
          f.insert(f.begin() + 2 + na_ + nb_, k - nb_, func_Tp());
          nb_ = k;
        }
      }
    }
    else if (p2.eof())
    {
      i = 0; 
    }
    else 
      throw gError("CuOx_plane: read: wrong token");
    f[i].read(p.GetSGML(el));
  }
  SetDebug(e);
}

void CuOx_plane::WriteBody(ostream& out, size_t shift) const
{
  out << PutTab(shift) << "<components> " << vf[0] << " " << vf[1]
    << " </components>" << endl;
  ref.write(out, &vf, shift);
  out << PutTab(shift) << "R*T*(";
  if (old_model)
  {
    out << "z*log(z)+(1-z)*log(1-z)" << endl;
    out << PutTab(shift) << "+";
  }
  out << "(c+x)*log(c+x)+(c-x)*log(c-x)" << endl;
  out << PutTab(shift) << "+(1-c+x)*log(1-c+x)+(1-c-x)*log(1-c-x))" << endl;
  out << PutTab(shift) << "<func_x> + </func_x> " << endl;
  fi(0).write(out, shift);
  out << PutTab(shift) << "<func_x> +z* </func_x> " << endl;
  fi(1).write(out, shift);
  size_t i;
  for (i = 0; i < na_; ++i)
  {
    out << PutTab(shift) << "<func_x> +z*(1-z)^" << i + 1
        << "* </func_x> " << endl;
    ai(i).write(out, shift);
  }
  if (nb_)
  {
    out << PutTab(shift) << "<func_x> +(c^2-x^2)* </func_x> " << endl;
    bi(0).write(out, shift);
    for (i = 1; i < nb_; ++i)
    {
      out << PutTab(shift) << "<func_x> +(c^2-x^2)*(1-z)^" << i
          << "* </func_x> " << endl;
      bi(i).write(out, shift);
    }
  }
}

double CuOx_plane::Zmix(function f, const StateTp &Tp, const double &z, 
    const double &x) const
{
  double sum = 0.;
  double c = z/2.;
  double lg;
  if (f == ::G || f == ::S)
  {
    if (old_model)
    {
      lg = log_z(z);
      if (lg != -HUGE_VAL)
        sum += z*lg;
      lg = log_z(1. - z);
      if (lg != -HUGE_VAL)
        sum += (1. - z)*lg;
    }
    lg = log_z(c + x);
    if (lg != -HUGE_VAL)
      sum += (c + x)*lg;
    lg = log_z(c - x);
    if (lg != -HUGE_VAL)
      sum += (c - x)*lg;
    lg = log_z(1. - c + x);
    if (lg != -HUGE_VAL)
      sum += (1. - c + x)*lg;
    lg = log_z(1. - c - x);
    if (lg != -HUGE_VAL)
      sum += (1. - c - x)*lg;
    if (f == ::G)
      sum *= global::R*Tp.T();
    else
      sum *= -global::R;
  }
  int i;
  sum += fi(0).Z(f, Tp);
  sum += z*fi(1).Z(f, Tp);
  if (na_)
  {
    sum += z*(1. - z)*ai(0).Z(f, Tp);
    for (i = 1; i < na_; ++i)
      sum += z*(1. - z)*ai(i).Z(f, Tp)*pow(1. - z, i);
  }
  if (nb_)
  {
    sum += (c*c - x*x)*bi(0).Z(f, Tp);
    for (i = 1; i < nb_; ++i)
      sum += (c*c - x*x)*bi(i).Z(f, Tp)*pow(1. - z, i);
  }
  return sum;
}

double CuOx_plane::dZdx_(function f, const StateTp &Tp, const double &z, 
    const double &x) const
{
  double sum = 0.;
  double c = z/2.;
  if (f == ::G || f == ::S)
  {
    double lg;
    lg = log_z(c + x);
    if (lg != -HUGE_VAL)
      sum += lg;
    else
      goto M_HUGE;
    lg = log_z(c - x);
    if (lg != -HUGE_VAL)
      sum -= lg;
    else
      goto P_HUGE;
    lg = log_z(1. - c + x);
    sum += lg;
    lg = log_z(1. - c - x);
    if (lg != -HUGE_VAL)
      sum -= lg;
    else
      goto P_HUGE;
    goto GOOD;
  M_HUGE:
    if (f == ::S)
      return HUGE_VAL;
    else
      return -HUGE_VAL;
  P_HUGE:
    if (f == ::S)
      return -HUGE_VAL;
    else
      return +HUGE_VAL;
  GOOD:
    if (f == ::G)
      sum *= global::R*Tp.T();
    else
      sum *= -global::R;
  }
  if (nb_)
  {
    sum -= 2.*x*bi(0).Z(f, Tp);
    for (int i = 1; i < nb_; ++i)
      sum -= 2.*x*bi(i).Z(f, Tp)*pow(1. - z, i);
  }
  return sum;
}

double CuOx_plane::dZdz(function f, const StateTp &Tp, const double &z, 
    const double &x) const
{
  double c = z/2.;
  double sum = 0.;
  if (f == ::G || f == ::S)
  {
    double lg;
    if (old_model)
    {
      lg = log_z(z);
      if (lg != -HUGE_VAL)
        sum += lg;
      else
        goto M_HUGE;
      lg = log_z(1. - z);
      if (lg != -HUGE_VAL)
        sum -= lg;
      else
        goto P_HUGE;
    }
    lg = 0.5*log_z(c + x);
    if (lg != -HUGE_VAL)
      sum += lg;
    else
      goto M_HUGE;
    lg = 0.5*log_z(c - x);
    if (lg != -HUGE_VAL)
      sum += lg;
    else
      goto M_HUGE;
    lg = 0.5*log_z(1. - c + x);
    sum -= lg;
    lg = 0.5*log_z(1. - c - x);
    if (lg != -HUGE_VAL)
      sum -= lg;
    else
      goto P_HUGE;
    goto GOOD;
  M_HUGE:
    if (f == ::S)
      return HUGE_VAL;
    else
      return -HUGE_VAL;
  P_HUGE:
    if (f == ::S)
      return -HUGE_VAL;
    else
      return +HUGE_VAL;
  GOOD:
    if (f == ::G)
      sum *= global::R*Tp.T();
    else
      sum *= -global::R;
  }
  sum += fi(1).Z(f, Tp);
  if (na_)
    sum += (1. - 2.*z)*ai(0).Z(f, Tp);
  if (na_ > 1)
  {
    sum += ai(1).Z(f, Tp)*((1. - 2*z)*(1. - z) - z*(1. - z));
    for (int i = 2; i < na_; ++i)
      sum += ai(i).Z(f, Tp)*
        ((1. - 2*z)*pow(1. - z, i) - z*(1. - z)*i*pow(1. - z, i - 1));
  }
  if (nb_)
    sum += c*bi(0).Z(f, Tp);
  if (nb_ > 1)
  {
    sum += bi(1).Z(f, Tp)*(c*(1. - z) - (c*c - x*x));
    for (int i = 2; i < nb_; ++i)
      sum += bi(i).Z(f, Tp)*
        (c*pow(1. - z, i) - (c*c - x*x)*i*pow(1. - z, i - 1));
  }
  return sum;
}

class solve_x
{
public:
  double T;
  double z;
  double sum;
  ofstream *debug;
  solve_x(const double &T_, const double &z_, const double &sum_, 
      ofstream *debug_)
    : T(T_), z(z_), sum(sum_), debug(debug_) {}
  double f(const double &x)
  {
    double c = z/2.;
    double res = global::R*T*log((c + x)*(1. - c + x)/(c - x)/(1. - c - x)) 
      - 2.*x*sum;
    if (debug)
      *debug << setw(12) << x << setw(12) << res << endl;
    return res;
  }
};

double CuOx_plane::x_eq(const StateTp &Tp, const double &z) const
{
  if (debug)
  {
    *debug << "T=" << Tp.T() << " p=" << Tp.p() << " z=" << z << endl;
    *debug << setw(12) << "x" << setw(12) << "f" << endl;
  }
  int i;
  double c = z/2.;
  double eps = sqrt(DBL_EPSILON);
  double x;
  if (c < eps || c > 1. - eps)
  {
    if (debug)
      *debug << "x is zero because of c = " << c << endl << endl;
    return 0.;
  }
  double sum = 0.;
  if (nb_)
  {
    sum += bi(0).G(Tp);
    for (i = 1; i < nb_; ++i)
      sum += bi(i).G(Tp)*pow(1. - z, i);
  }
  if (global::R*Tp.T()/c/(1. - c) - sum >= 0.)
  {
    x = 0.;
    if (debug)
      *debug << "x is zero because the derivative is positive" << endl 
        << endl;
  }
  else
  {
    try
    {
      double a = c*eps;
      double b;
      if (c <= 0.5)
        b = c*(1. - eps);
      else
        b = (1. - c)*(1. - eps);
      solve_x sol(Tp.T(), z, sum, debug);
      x = rroot(makeFunctor((D2D*)0, sol, &solve_x::f), a, b);
    }
    catch (err_rroot)
    {
      if (Zmix(::G, Tp, z, 0.) < Zmix(::G, Tp, z, c))
        x = 0.;
      else
        x = c;
    }
    if (debug)
      *debug << "final x is " << x << endl << endl;
  }
  return x;
}

/*
  DEPEND sum, z;
  DEPEND c, z;
  OFF EXP;
  ON GCD;
  ON EZGCD;
  f := R*T*log((c + x)*(1 - c + x)/(c - x)/(1 - c - x)) - 2*x*sum;
  %  c := z/2;
  dfdx := df(f, x);
  dfdz := df(f, z);
  LET df(c, z) = 1/2;
  dfdz := dfdz;
  dxdz := -dfdz/dfdx;
  ON FORT;
  dxdz;
      ANS=(((2.*C-1.)*R*T-2.*(C+X-1.)*(C+X)*(C-X-1.)*(C-X)*
     . DF(SUM,Z))*X)/(2.*((C**2-C+X**2)*R*T+(C+X-1.)*(C+X)*
     . (C-X-1.)*(C-X)*SUM))
*/

double CuOx_plane::dxdz(const StateTp &Tp, const double &z, 
    const double &x) const
{
  int i;
  double c = z/2.;
  double sum = 0.;
  if (z == 0.)
    return 0.;
  if (nb_)
  {
    sum = bi(0).G(Tp);
    for (i = 1; i < nb_; ++i)
      sum += bi(i).G(Tp)*pow(1. - z, i);
  }
  double dsumdz = 0.;
  if (nb_ > 1)
  {
    dsumdz = -bi(1).G(Tp);
    for (i = 2; i < nb_; ++i)
      dsumdz += -i*bi(i).G(Tp)*pow(1. - z, i - 1);
  }
  return
     (((2.*c-1.)*global::R*Tp.T()-2.*(c+x-1.)*(c+x)*(c-x-1.)*(c-x)*
     dsumdz)*x)/(2.*((pow(c, 2)-c+pow(x, 2))*global::R*Tp.T()+(c+x-1.)*(c+x)*
     (c-x-1.)*(c-x)*sum));
}

/*
  DEPEND sum, T;
  OFF EXP;
  ON GCD;
  ON EZGCD;
  f := R*T*log((c + x)*(1 - c + x)/(c - x)/(1 - c - x)) - 2*x*sum;
  dfdx := df(f, x);
  dfdT := df(f, T);
  dxdT := -dfdT/dfdx;
  ON FORT;
  dxdT;
      ANS=-((2.*DF(SUM,T)*X-LOG(((C+X)*(C-X-1.))/((C+X-1.)*
     . (C-X)))*R)*(C+X-1.)*(C+X)*(C-X-1.)*(C-X))/(2.*((C**2
     . -C+X**2)*R*T+(C+X-1.)*(C+X)*(C-X-1.)*(C-X)*SUM))
*/

double CuOx_plane::dxdT(const StateTp &Tp, const double &z, 
    const double &x) const
{
  int i;
  double c = z/2.;
  if (z == 0)
    return 0;
  double sum = 0.;
  double dsumdT = 0.;
  if (nb_)
  {
    sum = bi(0).G(Tp);
    dsumdT = -bi(0).Z(::S, Tp);
    for (i = 1; i < nb_; ++i)
    {
      sum += bi(i).G(Tp)*pow(1. - z, i);
      dsumdT -= bi(i).Z(::S, Tp)*pow(1. - z, i);
    }
  }
//  if (x == c/2.) then if would be good to do something
//  let us hope that this will not happen
  return
      -((2.*dsumdT*x-log(((c+x)*(c-x-1.))/((c+x-1.)*
      (c-x)))*global::R)*(c+x-1.)*(c+x)*(c-x-1.)*(c-x))/(2.*((pow(c, 2)
      -c+pow(x, 2))*global::R*Tp.T()+(c+x-1.)*(c+x)*(c-x-1.)*(c-x)*sum));
}

/*
  DEPEND sum, p;
  OFF EXP;
  ON GCD;
  ON EZGCD;
  f := R*T*log((c + x)*(1 - c + x)/(c - x)/(1 - c - x)) - 2*x*sum;
  dfdx := df(f, x);
  dfdp := df(f, p);
  dxdp := -dfdp/dfdx;
  ON FORT;
  dxdp;

      ANS=-((C+X-1.)*(C+X)*(C-X-1.)*(C-X)*DF(SUM,P)*X)/((C
     . **2-C+X**2)*R*T+(C+X-1.)*(C+X)*(C-X-1.)*(C-X)*SUM)

*/

double CuOx_plane::dxdp(const StateTp &Tp, const double &z, 
    const double &x) const
{
  int i;
  double c = z/2.;
  if (z == 0)
    return 0;
  double sum = 0.;
  double dsumdp = 0.;
  if (nb_)
  {
    sum = bi(0).G(Tp);
    dsumdp = bi(0).Z(::V, Tp);
    for (i = 1; i < nb_; ++i)
    {
      sum += bi(i).G(Tp)*pow(1. - z, i);
      dsumdp += bi(i).Z(::V, Tp)*pow(1. - z, i);
    }
  }
//  if (x == c/2.) then if would be good to do something
//  let us hope that this will not happen
  return
      -((c+x-1.)*(c+x)*(c-x-1.)*(c-x)*dsumdp*x)/((pow(c, 2.)
      -c+pow(x, 2.))*global::R*Tp.T()+(c+x-1.)*(c+x)*(c-x-1.)*(c-x)*sum);
}

bool CuOx_plane::IsStable(size_t i, const StateTp &Tp,
    const StateX &x) const
{
  double xeq = x_eq(Tp, x[1]);
  switch (i)
  {
    case 0:
      if (xeq)
        return false;
      else
        return true;
    case 1:
      if (xeq)
        return true;
      else
        return false;
    default:
      throw gError("CuOx_plane: phase not defined");
  }
}

double CuOx_plane::CheckBoundary(size_t i, size_t j, const StateTp &Tp, 
      const StateX &x) const
{
  double c = x[1]/2.;
  double s1 = 0.;
  if (nb_)
  {
    s1 = bi(0).G(Tp);
    for (int i = 1; i < nb_; ++i)
      s1 += bi(i).G(Tp)*pow((1. - x[1]), i);
  }
  return global::R*Tp.T() - c*(1. - c)*s1;
}

