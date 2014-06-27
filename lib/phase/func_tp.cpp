/*
Copyright (C) 1998, 1999, 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                       http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#include "func_tp_imp.h"
#include <typeinfo>
#include <set>
#include <algorithm>

static map_string_string *type;

Ref_func_Tp *func_Tp::PTR_NULL;
func_Tp::map_types *func_Tp::types = 0;
func_Tp::map_types *func_Tp::obj = 0;
func_Tp::Destruct func_Tp::clean;

void registerTp()
{
  func_Tp::RegisterType("compound_Tp", new compound_Tp);
  func_Tp::RegisterType("complex_Tp", new complex_Tp);
  func_Tp::RegisterType("calc_Tp", new calc_Tp);
}

void (*func_Tp::init_ptr[])() = {registerTp, 0};

istream& Ref_func_Tp::read(istream& in)
{
  func_Tp::init();
  SGML e(in);
  e.compare((*type)[typeid(*this).name()]);
  read(e);
  return in;
}

ostream& Ref_func_Tp::write(ostream& out, size_t shift) const
{
  func_Tp::init();
  out << PutTab(shift) << "<" << (*type)[typeid(*this).name()] 
    << " class=func_Tp";
  WriteAttributes(out, shift + 1);
  out << ">" << endl;
  WriteBody(out, shift + 1);
  out << PutTab(shift) << "</" << (*type)[typeid(*this).name()] << "> " << endl;
  return out;
}

struct GT_T
{
  const Ref_func_Tp &f;
  StateTp Tp;
  GT_T(const Ref_func_Tp &f_, const double &p) : f(f_) {Tp.p() = p;}
  double operator()(const double &T)
  {
    Tp.T() = T;
    return f.G(Tp)/T;
  }
};

struct G_T
{
  const Ref_func_Tp &f;
  StateTp Tp;
  G_T(const Ref_func_Tp &f_, const double &p) : f(f_) {Tp.p() = p;}
  double operator()(const double &T)
  {
    Tp.T() = T;
    return f.G(Tp);
  }
};

struct G_p
{
  const Ref_func_Tp &f;
  StateTp Tp;
  G_p(const Ref_func_Tp &f_, const double &T) : f(f_) {Tp.T() = T;}
  double operator()(const double &p)
  {
    Tp.p() = p;
    return f.G(Tp);
  }
};

struct G_Tp
{
  const Ref_func_Tp &f;
  StateTp Tp;
  G_Tp(const Ref_func_Tp &f_) : f(f_) {}
  double operator()(const double &T, const double &p)
  {
    Tp.T() = T;
    Tp.p() = p;
    return f.G(Tp);
  }
};

double Ref_func_Tp::num_H(const StateTp &Tp) const
{
  GT_T num(*this, Tp.p());
  return -Tp.T()*Tp.T()*FirstDerivative(num, Tp.T());
}

double Ref_func_Tp::num_S(const StateTp &Tp) const
{
  G_T num(*this, Tp.p());
  return -FirstDerivative(num, Tp.T());
}

double Ref_func_Tp::num_Cp(const StateTp &Tp) const
{
  G_T num(*this, Tp.p());
  return -Tp.T()*SecondDerivative(num, Tp.T());
}


double Ref_func_Tp::num_V(const StateTp &Tp) const
{
  G_p num(*this, Tp.T());
  return FirstDerivative(num, Tp.p());
}

double Ref_func_Tp::num_dVdT(const StateTp &Tp) const
{
  G_Tp num(*this);
  return MixedDerivative(num, Tp.T(), Tp.p());
}

double Ref_func_Tp::num_dVdp(const StateTp &Tp) const
{
  G_p num(*this, Tp.T());
  return SecondDerivative(num, Tp.p());
}

func_Tp::Destruct::~Destruct()
{
  if (func_Tp::types)
  {
    delete func_Tp::types;
    delete func_Tp::obj;
    func_Tp::types = 0;
    func_Tp::obj = 0;
    delete type;
  }
}

void func_Tp::init()
{
  if (!types)
  {
    type = new map_string_string;
    PTR_NULL = new null_Tp;
    types = new map_types;
    obj = new map_types;
    RegisterType("null_Tp", PTR_NULL);
    RegisterType("Cp_zero", new Cp_zero);
    RegisterType("Cp_const", new Cp_const);
    RegisterType("Cp_BB2", new Cp_BB2);
    RegisterType("Cp_BB2_Tref", new Cp_BB2_Tref);
    RegisterType("Cp_BB4", new Cp_BB4);
    RegisterType("IVT_Tp", new IVT_Tp);
    RegisterType("SGTE_Tp", new SGTE_Tp);
    RegisterType("ideal_gas", new ideal_gas);
    RegisterType("V_const", new V_const);
    RegisterType("alpha_const", new alpha_const);
    RegisterType("alpha_kappa_const", new alpha_kappa_const);
    RegisterType("alpha_kappa_const2", new alpha_kappa_const2);
    for (size_t i = 0; init_ptr[i]; ++i)
      (*init_ptr[i])();
  }
}

void func_Tp::RegisterType(const string& name, Ref_func_Tp* ptr)
{
  init();
  (*types)[name] = func_Tp(ptr);
  (*type)[typeid(*ptr).name()] = name;
}

void func_Tp::read(const SGML &e)
{
  string id;
  if (!(id = e.FindString("IDREF")).empty())
  {
    find(id);
  }
  else
  {
    create(e.name);
    cow();
    ptr->read(e);
    if (!(ptr->id = e.FindString("id")).empty())
    {
      pair<map_types_i, bool> i;
      i = obj->insert(map_types::value_type(ptr->id, *this));
      if (!i.second)
      {
        cout << "func_Tp: " << ptr->id << " is already defined - "
          << "set to anonymous" << endl;
        ptr->id.erase();
      }
    }
  }
}

ostream& func_Tp::write(ostream& out, size_t shift) const
{
  if (ptr->id.empty() || ptr->print)
  {
    out << PutTab(shift) << "<" << (*type)[typeid(*ptr).name()] 
      << " class=func_Tp";
    if (!ptr->id.empty())
    {
      out << " id=" << ptr->id;
    }
    ptr->WriteAttributes(out, shift + 1);
    out << ">" << endl;
    ptr->WriteBody(out, shift + 1);
    out << PutTab(shift) << "</" << (*type)[typeid(*ptr).name()] << "> " << endl;
    ptr->print = false;
  }
  else
  {
    out << PutTab(shift) << "<func_Tp class=func_Tp IDREF="
      << ptr->id << "></func_Tp>"
      << endl;
  }
  return out;
}

void null_Tp::read(const SGML &e)
{
  parser p(e);
  SGML el;
  if (!label_[0].empty())
    p.GetToken();
  for (size_t i = 0; i != coef_.size(); ++i)
  {
    coef_[i].read(p.GetSGML(el));
    if (!label_[i + 1].empty())
      p.SkipUntil('<');
  }
}

void null_Tp::WriteBody(ostream& out, size_t shift) const
{
  out << PutTab(shift) << label_[0] << endl;
  for (size_t i = 0; i != coef_.size(); ++i)
    out << PutTab(shift) << coef_[i] << " " << label_[i + 1] << endl;
}

double Cp_zero::Z(function f, const StateTp &Tp) const
{
  switch (f)
  {
    case ::G:
      return x(0) + Tp.T()*x(1);
    case ::H:
      return x(0);
    case ::S:
      return -x(1);
    case ::Cp:
    case ::V:
    case ::dVdT:
    case ::dVdp:
      return 0.;
    default:
      throw gError("Cp_zero: unknown code for the function");
  }
}

double Cp_const::Z(function f, const StateTp &Tp) const
{
  switch (f)
  {
    case ::G:
      return x(0) + Tp.T()*x(1) + Tp.T()*log(Tp.T())*x(2);
    case ::H:
      return x(0) - Tp.T()*x(2);
    case ::S:
      return -(x(1) + (log(Tp.T())+1.)*x(2));
    case ::Cp:
      return -x(2);
    case ::V:
    case ::dVdT:
    case ::dVdp:
      return 0.;
    default:
      throw gError("Cp_const: unknown code for the function");
  }
}

double Cp_BB2::Z(function f, const StateTp &Tp) const
{
  switch (f)
  {
    case ::G:
      return x(0) + Tp.T()*x(1) + Tp.T()*log(Tp.T())*x(2)
                 + sqrt(Tp.T())*x(3);
    case ::H:
      return x(0) - Tp.T()*x(2) + 0.5*sqrt(Tp.T())*x(3);
    case ::S:
      return -(x(1) + (log(Tp.T())+1.)*x(2) + 0.5/sqrt(Tp.T())*x(3));
    case ::Cp:
      return -x(2) + 0.25/sqrt(Tp.T())*x(3);
    case ::V:
    case ::dVdT:
    case ::dVdp:
      return 0.;
    default:
      throw gError("Cp_BB2: unknown code for the function");
  }
}

void Cp_BB2_Tref::read(const SGML &e)
{
  To = e.FindDouble("To", 1000.);
  null_Tp::read(e);
}

void Cp_BB2_Tref::WriteAttributes(ostream& out, size_t shift) const
{
  out << endl << PutTab(shift) << "To=" << To;
}

double Cp_BB2_Tref::Z(function f, const StateTp &Tp) const
{
  switch (f)
  {
    case ::G:
      return x(0) - Tp.T()*x(1) + (Tp.T() - To - Tp.T()*log(Tp.T()/To))*x(2)
                 + (sqrt(Tp.T()) - sqrt(To))*4.*x(3);
    case ::H:
      return x(0) + (Tp.T() - To)*x(2) + 2.*(sqrt(Tp.T()) - sqrt(To))*x(3);
    case ::S:
      return x(1) + log(Tp.T()/To)*x(2) - 2.*(1./sqrt(Tp.T()) - 1./sqrt(To))*x(3);
    case ::Cp:
      return x(2) + 1./sqrt(Tp.T())*x(3);
    case ::V:
    case ::dVdT:
    case ::dVdp:
      return 0.;
    default:
      throw gError("Cp_BB2_Tref: unknown code for the function");
  }
}

double Cp_BB4::Z(function f, const StateTp &Tp) const
{
  switch (f)
  {
    case ::G:
// g := x0 + x1*Te + x2*Te*log(Te) + x3*sqrt(Te) + x4/Te + x5/Te^2;
      return x(0) + Tp.T()*x(1) + Tp.T()*log(Tp.T())*x(2) + sqrt(Tp.T())*x(3) +
            x(4)/Tp.T() + x(5)*pow(Tp.T(), -2);
    case ::H:
//h := g + Te*s;
//             1          -1                     -1          -2
//H :=  - ( - ---*SQRT(TE)  *TE*X3 + TE*X2 - 2*TE  *X4 - 3*TE  *X5 - X0)
//             2
// -Te^2*df(g/Te,Te);
//        1          -1                     -1          -2
// - ( - ---*SQRT(TE)  *TE*X3 + TE*X2 - 2*TE  *X4 - 3*TE  *X5 - X0)
//        2
      return x(0) - Tp.T()*x(2) + 0.5*sqrt(Tp.T())*x(3) +
            2.*x(4)/Tp.T() + 3.*x(5)*pow(Tp.T(), -2);
    case ::S:
//s := -df(g, te);
//          1          -1                     -2          -3
//S :=  - (---*SQRT(TE)  *X3 + LOG(TE)*X2 - TE  *X4 - 2*TE  *X5 + X1 +  X2)
//          2
      return -(x(1) + (log(Tp.T())+1.)*x(2) + 0.5/sqrt(Tp.T())*x(3)
            -x(4)*pow(Tp.T(), -2) - 2.*x(5)*pow(Tp.T(), -3));
    case ::Cp:
//cp := df(h, Te);
//              1          -1          -2          -3
//CP :=  - ( - ---*SQRT(TE)  *X3 + 2*TE  *X4 + 6*TE  *X5 + X2)
//              4
//Te*df(s, Te);
//        1          -1          -2          -3
// - ( - ---*SQRT(TE)  *X3 + 2*TE  *X4 + 6*TE  *X5 + X2)
//        4
      return -x(2) + 0.25/sqrt(Tp.T())*x(3)
            - 2.0*x(4)*pow(Tp.T(), -2) - 6.*x(5)*pow(Tp.T(), -3);
    case ::V:
    case ::dVdT:
    case ::dVdp:
      return 0.;
    default:
      throw gError("Cp_BB4: unknown code for the function");
  }
}

double IVT_Tp::Z(function f, const StateTp &Tp) const
{
  switch (f)
  {
    case ::G:
// g1 := x0 + x1*Te + x2*Te*log(Te) + x3*Te^2 + x4*Te^3 + x5*Te^4 + x6/Te;
      return x(0) + x(1)*Tp.T() + x(2)*Tp.T()*log(Tp.T())
                 + x(3)*pow(Tp.T(), 2) + x(4)*pow(Tp.T(), 3)
                 + x(5)*pow(Tp.T(), 4) + x(6)/Tp.T();
    case ::H:
//h1 := g1 + Te*s1;
//              4          3        2                  -1
//H1 :=  - (3*TE *X5 + 2*TE *X4 + TE *X3 + TE*X2 - 2*TE  *X6 - X0)
//
//-Te^2*df(g1/Te,Te);
//        4          3        2                  -1
// - (3*TE *X5 + 2*TE *X4 + TE *X3 + TE*X2 - 2*TE  *X6 - X0)
      return x(0) - x(2)*Tp.T() - x(3)*pow(Tp.T(), 2)
                 - x(4)*2.*pow(Tp.T(), 3)
                 - x(5)*3.*pow(Tp.T(), 4)
                 + x(6)*2./Tp.T();
    case ::S:
//s1 := -df(g1, te);
//                           3          2                  -2
//S1 :=  - (LOG(TE)*X2 + 4*TE *X5 + 3*TE *X4 + 2*TE*X3 - TE  *X6 + X1 +
//
//           X2)
      return -(x(1) + x(2)*(log(Tp.T())+1.)
                   + x(3)*2.*Tp.T()
                   + x(4)*3.*pow(Tp.T(), 2)
                   + x(5)*4.*pow(Tp.T(), 3)
                   - x(6)/pow(Tp.T(), 2));
    case ::Cp:
//cp1 := df(h1, Te);
//                3          2                    -2
//CP1 :=  - (12*TE *X5 + 6*TE *X4 + 2*TE*X3 + 2*TE  *X6 + X2)
//
//Te*df(s1, Te);
//         3          2                    -2
// - (12*TE *X5 + 6*TE *X4 + 2*TE*X3 + 2*TE  *X6 + X2)
      return -(x(2) + x(3)*2.*Tp.T()
                   + x(4)*6.*pow(Tp.T(), 2)
                   + x(5)*12.*pow(Tp.T(), 3)
                   + x(6)*2./pow(Tp.T(), 2));
    case ::V:
    case ::dVdT:
    case ::dVdp:
      return 0.;
    default:
      throw gError("IVT_Tp: unknown code for the function");
  }
}

void SGTE_Tp::read(const SGML &e)
{
  Bo = e.FindDouble("Bo");
  if (Bo != 0.)
  {
    pm = e.FindDouble("pm");
    if (pm == 0.)
      throw gError("SGTE_Tp: pm is not defined");
    D = 518./1125. + 11692./15975.*(1./pm - 1.);
    Tc = e.FindDouble("Tc", -1.);
    if (Tc < 0.)
      throw gError("SGTE_Tp: Tc is not defined");
    label(8) = "/T^9+Gmag)";
  }
  else
    label(8) = "/T^9)";
  null_Tp::read(e);
}

void SGTE_Tp::WriteAttributes(ostream& out, size_t shift) const
{
  if (Bo != 0.)
    out << endl << PutTab(shift) << " Bo=" << Bo 
      << endl << PutTab(shift) << " Tc=" << Tc 
      << endl << PutTab(shift) << " pm=" << pm;
}

double SGTE_Tp::Z(function f, const StateTp &Tp) const
{
  double mag;
  double tau;
  switch (f)
  {
    case ::G:
//g := x0 + x1*Te + x2*Te*log(Te) + x3*Te^2 + x4*Te^3 + x5*Te^7 + x6/Te +
//     x7/Te^9;
      if (Bo != 0.)
      {
        tau = Tp.T()/Tc;
        mag = global::R*Tp.T()*log(1. + Bo);
        if (tau <= 1.)
          mag *= 1. - (79./tau/140./pm + 474./497.*(1./pm - 1.)
                      *(pow(tau, 3)/6. + pow(tau, 9)/135.
                      + pow(tau, 15)/600.))/D;
        else
          mag *= -(pow(tau, -5)/10. + pow(tau, -15)/315.
                   + pow(tau, -25)/1500.)/D;
      }
      else
        mag = 0.;
      return x(0) + x(1)*Tp.T() + x(2)*Tp.T()*log(Tp.T())
                 + x(3)*pow(Tp.T(), 2) + x(4)*pow(Tp.T(), 3)
                 + x(5)*pow(Tp.T(), 7) + x(6)/Tp.T()
                 + x(7)/pow(Tp.T(), 9) + mag;
    case ::H:
//h := g + Te*s;
//             7          3        2                  -1           -9
//H :=  - (6*TE *X5 + 2*TE *X4 + TE *X3 + TE*X2 - 2*TE  *X6 - 10*TE  *
//
//         X7 - X0)
//
//-Te^2*df(g/Te,Te);
//        7          3        2                  -1           -9
// - (6*TE *X5 + 2*TE *X4 + TE *X3 + TE*X2 - 2*TE  *X6 - 10*TE  *X7 -
//
//    X0)
      if (Bo != 0.)
      {
        tau = Tp.T()/Tc;
        mag = global::R*Tp.T()*log(1. + Bo);
        if (tau <= 1.)
          mag *= (-79./tau/140./pm + 474./497.*(1./pm - 1.)
                  *(pow(tau, 3)/2. + pow(tau, 9)/15.
                  + pow(tau, 15)/40.))/D;
        else
          mag *= -(pow(tau, -5)/2. + pow(tau, -15)/21.
                  + pow(tau, -25)/60.)/D;
      }
      else
        mag = 0.;
      return x(0) - x(2)*Tp.T() - x(3)*pow(Tp.T(), 2)
                 - x(4)*2.*pow(Tp.T(), 3)
                 - x(5)*6.*pow(Tp.T(), 7)
                 + x(6)*2./Tp.T()
                 + x(7)*10./pow(Tp.T(), 9) + mag;
    case ::S:
//s := -df(g, te);
//                          6          2                  -2
//S :=  - (LOG(TE)*X2 + 7*TE *X5 + 3*TE *X4 + 2*TE*X3 - TE  *X6 - 9*
//
//           -10
//         TE   *X7 + X1 + X2)
      if (Bo != 0.)
      {
        tau = Tp.T()/Tc;
        mag = -global::R*log(1. + Bo);
        if (tau <= 1.)
          mag *= 1. - (474./497.*(1./pm - 1.)
                      *(2.*pow(tau, 3)/3. + 2.*pow(tau, 9)/27.
                      + 2.*pow(tau, 15)/75.))/D;
        else
          mag *= (2.*pow(tau, -5)/5. + 2.*pow(tau, -15)/45.
                  + 2.*pow(tau, -25)/125.)/D;
      }
      else
        mag = 0.;
      return -(x(1) + x(2)*(log(Tp.T())+1.)
                   + x(3)*2.*Tp.T()
                   + x(4)*3.*pow(Tp.T(), 2)
                   + x(5)*7.*pow(Tp.T(), 6)
                   - x(6)/pow(Tp.T(), 2)
                   - x(7)*9./pow(Tp.T(), 10)) + mag;
    case ::Cp:
//cp := df(h, Te);
//               6          2                    -2           -10
//CP :=  - (42*TE *X5 + 6*TE *X4 + 2*TE*X3 + 2*TE  *X6 + 90*TE   *X7 +
//
//          X2)
//
//Te*df(s, Te);
//         6          2                    -2           -10
// - (42*TE *X5 + 6*TE *X4 + 2*TE*X3 + 2*TE  *X6 + 90*TE   *X7 + X2)
      if (Bo != 0.)
      {
        tau = Tp.T()/Tc;
        mag = global::R*log(1. + Bo);
        if (tau <= 1.)
          mag *= (474./497.*(1./pm - 1.)
                  *(2.*pow(tau, 3) + 1.*pow(tau, 9)/3.
                      + 2.*pow(tau, 15)/5.))/D;
        else
          mag *= (2.*pow(tau, -5) + 2.*pow(tau, -15)/3.
                 + 2.*pow(tau, -25)/5.)/D;
      }
      else
        mag = 0.;
      return -(x(2) + x(3)*2.*Tp.T()
                   + x(4)*6.*pow(Tp.T(), 2)
                   + x(5)*42.*pow(Tp.T(), 6)
                   + x(6)*2./pow(Tp.T(), 2)
                   + x(7)*90./pow(Tp.T(), 10)) + mag;
    case ::V:
    case ::dVdT:
    case ::dVdp:
      return 0.;
    default:
      throw gError("SGTE_Tp: unknown code for the function");
  }
}

double ideal_gas::Z(function f, const StateTp &Tp) const
{
  switch (f)
  {
    case ::G:
      return global::R*Tp.T()*log(Tp.p());
    case ::H:
      return 0.;
    case ::S:
      return -global::R*log(Tp.p());
    case ::V:
      return global::R*Tp.T()/Tp.p();
    case ::Cp:
      return 0.;
    case ::dVdT:
      return global::R/Tp.p();
    case ::dVdp:
      return -global::R*Tp.T()/Tp.p()/Tp.p();
    default:
      throw gError("ideal_gas: unknown code for the function");
  }
}

double V_const::Z(function f, const StateTp &Tp) const
{
  switch (f)
  {
    case ::G:
    case ::H:
      return Tp.p()*x(0);
    case ::S:
      return 0;
    case ::V:
      return x(0);
    case ::Cp:
    case ::dVdT:
    case ::dVdp:
      return 0.;
    default:
      throw gError("V_const: unknown code for the function");
  }
}

double alpha_const::Z(function f, const StateTp &Tp) const
{
  switch (f)
  {
    case ::G:
      return x(0)*exp(x(1)*Tp.T())*(Tp.p() - 1.);
    case ::H:
      return x(0)*exp(x(1)*Tp.T())*(Tp.p() - 1.)*(1. - x(1)*Tp.T());
    case ::S:
      return -x(0)*x(1)*exp(x(1)*Tp.T())*(Tp.p() - 1.);
    case ::Cp:
      return -x(0)*exp(x(1)*Tp.T())*(Tp.p() - 1.)*x(1)*x(1)*Tp.T();
    case ::V:
      return x(0)*exp(x(1)*Tp.T());
    case ::dVdT:
      return x(1)*x(0)*exp(x(1)*Tp.T());
    case ::dVdp:
      return 0.;
    default:
      throw gError("alpha_const: unknown code for the function");
  }
}

double alpha_kappa_const::Z(function f, const StateTp &Tp) const
{
  switch (f)
  {
    case ::G:
      return x(0)*exp(x(1)*Tp.T())*(1. - exp(x(2)*(1. - Tp.p())));
    case ::H:
      return x(0)*exp(x(1)*Tp.T())*(1. - exp(x(2)*(1. - Tp.p())))*(1. - x(1)*Tp.T());
    case ::S:
      return -x(0)*x(1)*exp(x(1)*Tp.T())*(1. - exp(x(2)*(1. - Tp.p())));
    case ::Cp:
      return -x(0)*exp(x(1)*Tp.T())*(1. - exp(x(2)*(1. - Tp.p())))*x(1)*x(1)*Tp.T();
    case ::V:
      return x(0)*x(2)*exp(x(1)*Tp.T())*exp(x(2)*(1. - Tp.p()));
    case ::dVdT:
      return x(1)*x(0)*x(2)*exp(x(1)*Tp.T())*exp(x(2)*(1. - Tp.p()));
    case ::dVdp:
      return -x(2)*x(0)*x(2)*exp(x(1)*Tp.T())*exp(x(2)*(1. - Tp.p()));
    default:
      throw gError("alpha_kappa_const: unknown code for the function");
  }
}

double alpha_kappa_const2::Z(function f, const StateTp &Tp) const
{
  switch (f)
  {
    case ::G:
      if (x(2))
        return x(0)/x(2)*exp(x(1)*Tp.T())*(1. - exp(x(2)*(1. - Tp.p())));
      else
        return x(0)*exp(x(1)*Tp.T())*(Tp.p() - 1.);
    case ::H:
      if (x(2))
        return x(0)/x(2)*exp(x(1)*Tp.T())*(1. - exp(x(2)*(1. - Tp.p())))*(1. - x(1)*Tp.T());
      else
        return x(0)*exp(x(1)*Tp.T())*(Tp.p() - 1.)*(1. - x(1)*Tp.T());
    case ::S:
      if (x(2))
        return -x(0)/x(2)*x(1)*exp(x(1)*Tp.T())*(1. - exp(x(2)*(1. - Tp.p())));
      else
        return -x(0)*x(1)*exp(x(1)*Tp.T())*(Tp.p() - 1.);
    case ::Cp:
      if (x(2))
        return -x(0)/x(2)*exp(x(1)*Tp.T())*(1. - exp(x(2)*(1. - Tp.p())))*x(1)*x(1)*Tp.T();
      else
        return -x(0)*exp(x(1)*Tp.T())*(Tp.p() - 1.)*x(1)*x(1)*Tp.T();
    case ::V:
      return x(0)*exp(x(1)*Tp.T())*exp(x(2)*(1. - Tp.p()));
    case ::dVdT:
      return x(1)*x(0)*exp(x(1)*Tp.T())*exp(x(2)*(1. - Tp.p()));
    case ::dVdp:
      return -x(2)*x(0)*exp(x(1)*Tp.T())*exp(x(2)*(1. - Tp.p()));
    default:
      throw gError("alpha_kappa_const2: unknown code for the function");
  }
}

struct lim
{
  func_Tp f;
  double Ti, Tf, pi, pf;
  lim() {Ti = 0.; Tf = HUGE_VAL; pi = 0.; pf = HUGE_VAL;}
  friend bool operator<(const lim &x, const lim &y)
  {
    if (x.Ti < y.Ti && x.Tf <= y.Ti)
      return true;
    else if (x.Ti == y.Ti && x.Tf == y.Tf)
    {
      if (x.pi < y.pi && x.pf <= y.pi)
        return true;
      else if (x.pi >= y.pf)
        return false;
      else
        throw gError("compound_Tp: crossing intervals");
    }
    else if (x.Ti >= y.Tf)
      return false;
    else
      throw gError("compound_Tp: crossing intervals");
  }
};

void compound_Tp::read(const SGML &e)
{
  clear();
  parser p(e);
  SGML el;
  set<lim, less<lim> > t;
  while (!p.eof())
  {
    lim x;
    x.f.read(p.GetSGML(el));
    if (!dynamic_cast<const simple_Tp*>(x.f.pointer()))
      throw gError("compound_Tp: not derived from simple_Tp");
    p.GetSGML(el);
    el.compare("limit_Tp");
    parser current(el);
    for (size_t j = 0; j < 2; ++j)
    {
      if (current.GetToken() == "T")
      {
        current.SkipChar('=');
        current.SkipChar('{');
        x.Ti = current.GetDouble();
        current.SkipChar(',');
        x.Tf = current.GetDouble();
        current.SkipChar('}');
      }
      else if (current.token == "p")
      {
        current.SkipChar('=');
        current.SkipChar('{');
        x.pi = current.GetDouble();
        current.SkipChar(',');
        x.pf = current.GetDouble();
        current.SkipChar('}');
      }
    }
    t.insert(t.end(), x);
  }
  set<lim, less<lim> >::iterator i;
  size_t iT = 1, ip = 1;
  for (i = t.begin(); i != t.end(); ++i)
  {
    if (i == t.begin())
    {
      limT_[0] = (*i).Ti;
      limT_[1] = (*i).Tf;
      limp_[0] = (*i).pi;
      limp_[1] = (*i).pf;
    }
    else if (iT == 1)
    {
      if (limT_[0] == (*i).Ti && limT_[1] == (*i).Tf)
      {
        ++ip;
        if (limp_[ip - 1] == (*i).pi)
          limp_.push_back((*i).pf);
        else
          throw gError("compound_Tp: wrong intervals");
      }
      else
      {
        ++iT;
        ip = 1;
        if (limT_[iT - 1] == (*i).Ti &&
            limp_[ip - 1] == (*i).pi &&
            limp_[ip] == (*i).pf)
          limT_.push_back((*i).Tf);
        else
          throw gError("compound_Tp: wrong intervals");
      }
    }
    else
    {
      if (++ip < limp_.size())
      {
        if (!(limp_[ip - 1] == (*i).pi && limp_[ip] == (*i).pf))
          throw gError("compound_Tp: wrong intervals");
      }
      else
      {
        ++iT;
        ip = 1;
        if (limT_[iT - 1] == (*i).Ti &&
            limp_[ip - 1] == (*i).pi &&
            limp_[ip] == (*i).pf)
          limT_.push_back((*i).Tf);
        else
          throw gError("compound_Tp: wrong intervals");
      }
    }
  }
  if ((limT_.size() - 1)*(limp_.size() - 1) == t.size())
    for (i = t.begin(); i != t.end(); ++i)
    {
      if (i == t.begin())
        fptr_[0] = (*i).f;
      else
        fptr_.push_back((*i).f);
    }
  else
    throw gError("compound_Tp: wrong intervals");
}

void compound_Tp::WriteBody(ostream& out, size_t shift) const
{
  size_t iT = 1, ip = 1;
  for (vec_func_Tp_ci i = fptr_.begin(); i != fptr_.end(); ++i)
  {
    (*i).write(out, shift);
    out << PutTab(shift) << "<limit_Tp> ";
    out << "T = {" << limT_[iT - 1] << ", " << CheckInf(limT_[iT]) << "} ";
    out << "p = {" << limp_[ip - 1] << ", " << CheckInf(limp_[ip]) << "} ";
    out << "</limit_Tp>" << endl;
    if (++ip == limp_.size())
    {
      ++iT;
      ip = 1;
    }
  }
  if (iT != limT_.size() && ip != 1)
    throw gError("compound_Tp: corrupt structure");
}

const func_Tp& compound_Tp::find(const StateTp &Tp) const
{
  size_t iT =
    lower_bound(limT_.begin(), limT_.end(), Tp.T()) - limT_.begin();
  if (iT == 0 && Tp.T() == limT_.front()) ++iT;
  if (iT == 0 || iT == limT_.size())
    throw gError(string("compound_Tp: T out of range = ").append(ObjToString(Tp.T(), 6)));
  size_t ip =
    lower_bound(limp_.begin(), limp_.end(), Tp.p()) - limp_.begin();
  if (ip == 0 && Tp.p() == limp_.front()) ++ip;
  if (ip == 0 || ip == limp_.size())
    throw gError(string("compound_Tp: p out of range = ").append(ObjToString(Tp.T(), 6)));
  return fptr_[(limp_.size() - 1)*(iT - 1) + (ip - 1)];
}

void calc_Tp::AnalizeId(const string &id)
{
  if (id == "T")
  {
    curr_tok.key = DOUBLE_PTR_PTR;
    curr_tok.x_ptr_ptr = &ptrT;
  }
  else if (id == "p")
  {
    curr_tok.key = DOUBLE_PTR_PTR;
    curr_tok.x_ptr_ptr = &ptrp;
  }
  else
   throw gError("ID_UNDEF");
}

double calc_Tp::Z(function f, const StateTp &Tp) const
{
  switch (f)
  {
    case ::G:
      ptrT = Tp.ptrT();
      ptrp = Tp.ptrp();
      return est();
    case ::H:
      return num_H(Tp);
    case ::S:
      return num_S(Tp);
    case ::Cp:
      return num_Cp(Tp);
    case ::V:
      return num_V(Tp);
    case ::dVdT:
      return num_dVdT(Tp);
    case ::dVdp:
      return num_dVdp(Tp);
    default:
      throw gError("calc_Tp: unknown code for the function");
  }
}


