/*
Copyright (C) 1998, 1999, 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                       http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#include "cuox_or.h"

void RegisterCuOx_ordered_plane()
{
  phase::RegisterType("CuOx_ordered_plane", new CuOx_ordered_plane);
}

void CuOx_ordered_plane::read(const SGML &e_)
{
  clear();
  size_t i, k;
  SGML e = e_;
  ReadComponents(vf, e, ref);
  parser p(e);
  SGML el;
  if (ref.size() != 2)
  {
    p.GetSGML(el);
    el.compare("Reference");
    ref.read(el);
  }
  if (p.GetToken() == "Gid")
    mul_ent = 2;
  else if (p.token == "2*R*T*")
  {
    mul_ent = 2;
  }
  else if (p.token == "R*T*")
  {
    mul_ent = 1;
  }
  else
    throw gError("CuOx_ordered_plane: expected ideal Gibbs energy");
  p.SkipUntil('<');
  while (!p.eof())
  {
    p.GetSGML(el);
    el.compare("func_x");
    parser p2(el);
    if (p2.GetToken() == "+1*" || p2.token == "+")
      i = 0;
    else if (p2.token == "+z*")
    {
      if (!p2.SkipUntil('^'))
        i = 1;
      else
      {
        p2.SkipChar('^');
        k = p2.GetInt();
        if (k == 0)
          throw gError("CuOx_ordered_plane: read: wrong token");
        i = k + 1;
        if (k > na_)
        {
          f.insert(f.begin() + 2 + na_, k - na_, func_Tp());
          na_ = k;
        }
      }
    }
    else
      throw gError("CuOx_ordered_plane: read: wrong token");
    f[i].read(p.GetSGML(el));
  }
}

void CuOx_ordered_plane::WriteBody(ostream& out, size_t shift) const
{
  WriteComponents(out, shift);
  ref.write(out, &vf, shift);
  out << PutTab(shift);
  if (mul_ent != 1)
  {
    out << mul_ent << "*";
  }
  out << "R*T*(z*log(z)+(1-z)*log(1-z))" << endl;
  out << PutTab(shift) << "<func_x> + </func_x> " << endl;
  fi(0).write(out, shift + 1);
  out << PutTab(shift) << "<func_x> +z* </func_x> " << endl;
  fi(1).write(out, shift);
  size_t i;
  for (i = 0; i < na_; ++i)
  {
    out << PutTab(shift) << "<func_x> +z*(1-z)^" << i + 1
        << "* </func_x> " << endl;
    ai(i).write(out, shift);
  }
}

double CuOx_ordered_plane::Zmix(function f, const StateTp &Tp,
                                const double &z) const
{
  double sum = 0.;
  double lg;
  if (f == ::G || f == ::S)
  {
    lg = log_z(z);
    if (lg != -HUGE_VAL)
      sum += z*lg;
    lg = log_z(1. - z);
    if (lg != -HUGE_VAL)
      sum += (1. - z)*lg;

    if (f == ::G)
      sum *= global::R*Tp.T();
    else
      sum *= -global::R;
    sum *= mul_ent;
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
  return sum;
}

double CuOx_ordered_plane::dZdz(function f, const StateTp &Tp,
                                const double &z) const
{
  double sum = 0.;
  if (f == ::G || f == ::S)
  {
    double lg;
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
    sum *= mul_ent;
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
  return sum;
}

const vec_double& CuOx_ordered_plane::dZdx(function f, index i, 
    const StateTp &Tp, const StateX &x) const
{
  vz[0] = vz[1] = 0.;
  if (i & ::ref)
    ref.dZdx(f, Tp, x, vz);
  if ((i & ::mix) == ::mix)
    vz[1] += dZdz(f, Tp, x[1]);
  return vz;
}

const vec_double& CuOx_ordered_plane::z(function f, index i, const StateTp &Tp,
                                const StateX &x) const
{
  vz[0] = vz[1] = 0.;
  if (i & ::ref)
    ref.dZdx(f, Tp, x, vz);
  if ((i & ::mix) == ::mix)
  {
    double Zm = Zmix(f, Tp, x[1]);
    double dZx = dZdz(f, Tp, x[1]);
    vz[0] += Zm - x[1]*dZx;
    vz[1] += Zm + (1. - x[1])*dZx;
  }
  return vz;
}

