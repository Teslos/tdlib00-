/*
Copyright (C) 1998, 1999, 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                       http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#ifndef __FUNC_TP_IMP_H
#define __FUNC_TP_IMP_H

#include <func_tp.h>

class null_Tp : public simple_Tp
{
protected:
  vec_coef coef_;
  vec_string label_;

  const double& x(size_t i) const {return coef_[i].x();}
  string& label(size_t i) {return label_[i];}

public:

  null_Tp() : label_(1) {label_[0] = "(0)";}
  explicit null_Tp(size_t n) : coef_(n), label_(n + 1) {}

  virtual null_Tp* clone() const {return new null_Tp(*this);}

  virtual void read(const SGML &e);
  virtual void WriteBody(ostream &out, size_t shift = 0) const;

  virtual double Z(function f, const StateTp& Tp) const {return 0;}
};

class Cp_zero : public null_Tp
{
public:
  Cp_zero() : null_Tp(2)
  {
    label(0) = "(";
//    name(0) = "a";
    label(1) = "+";
//    name(1) = "b";
    label(2) = "*T)";
  }
  virtual Cp_zero* clone() const {return new Cp_zero(*this);}
  virtual double Z(function f, const StateTp &Tp) const;
};

class Cp_const : public null_Tp
{
public:
  Cp_const() : null_Tp(3)
  {
    label(0) = "(";
//    name(0) = "a";
    label(1) = "+";
//    name(1) = "b";
    label(2) = "*T+";
//    name(2) = "c";
    label(3) = "*T*log(T))";
  }
  virtual Cp_const* clone() const {return new Cp_const(*this);}
  virtual double Z(function f, const StateTp &Tp) const;
};

class Cp_BB2 : public null_Tp
{
public:
  Cp_BB2() : null_Tp(4)
  {
    label(0) = "(";
//    name(0) = "a";
    label(1) = "+";
//    name(1) = "b";
    label(2) = "*T+";
//    name(2) = "c";
    label(3) = "*T*log(T)+";
//    name(3) = "d";
    label(4) = "*sqrt(T))";
  }
  virtual Cp_BB2* clone() const {return new Cp_BB2(*this);}
  virtual double Z(function f, const StateTp &Tp) const;
};

class Cp_BB2_Tref : public null_Tp
{
  double To;
public:
  Cp_BB2_Tref() : null_Tp(4)
  {
    label(0) = "(";
//    name(0) = "Ho";
    label(1) = "-";
//    name(1) = "So";
    label(2) = "*T+";
//    name(2) = "a";
    label(3) = "*(T-To-T*log(T/To))+4*";
//    name(3) = "b";
    label(4) = "*(sqrt(T)-sqrt(To)))";
    To = 1000.;
  }
  virtual Cp_BB2_Tref* clone() const {return new Cp_BB2_Tref(*this);}
  virtual void read(const SGML& e);
  virtual void WriteAttributes(ostream& out, size_t shift = 0) const;
  virtual double Z(function f, const StateTp &Tp) const;
};

class Cp_BB4 : public null_Tp
{
public:
  Cp_BB4() : null_Tp(6)
  {
    label(0) = "(";
//    name(0) = "a";
    label(1) = "+";
//    name(1) = "b";
    label(2) = "*T+";
//    name(2) = "c";
    label(3) = "*T*log(T)+";
//    name(3) = "d";
    label(4) = "*sqrt(T)+";
//    name(4) = "e";
    label(5) = "/T+";
//    name(5) = "f";
    label(6) = "/T^2)";
  }
  virtual Cp_BB4* clone() const {return new Cp_BB4(*this);}
  virtual double Z(function f, const StateTp &Tp) const;
};

class IVT_Tp : public null_Tp
{
public:
  IVT_Tp() : null_Tp(7)
  {
    label(0) = "(";
//    name(0) = "a";
    label(1) = "+";
//    name(1) = "b";
    label(2) = "*T+";
//    name(2) = "c";
    label(3) = "*T*log(T)+";
//    name(3) = "d";
    label(4) = "*T^2+";
//    name(4) = "e";
    label(5) = "*T^3+";
//    name(5) = "f";
    label(6) = "*T^4+";
//    name(6) = "g";
    label(7) = "/T)";
  }
  virtual IVT_Tp* clone() const {return new IVT_Tp(*this);}
  virtual double Z(function f, const StateTp &Tp) const;
};

class SGTE_Tp : public null_Tp
{
  double Bo, Tc, pm, D;
public:
  SGTE_Tp() : null_Tp(8)
  {
    label(0) = "(";
//    name(0) = "a";
    label(1) = "+";
//    name(1) = "b";
    label(2) = "*T+";
//    name(2) = "c";
    label(3) = "*T*log(T)+";
//    name(3) = "d";
    label(4) = "*T^2+";
//    name(4) = "e";
    label(5) = "*T^3+";
//    name(5) = "f";
    label(6) = "*T^7+";
//    name(6) = "g";
    label(7) = "/T+";
//    name(7) = "h";
    label(8) = "/T^9)";
    Bo = 0.;
  }
  virtual SGTE_Tp* clone() const {return new SGTE_Tp(*this);}
  virtual void read(const SGML &e);
  virtual void WriteAttributes(ostream& out, size_t shift = 0) const;
  virtual double Z(function f, const StateTp &Tp) const;
};

class ideal_gas : public null_Tp
{
public:
  ideal_gas() : null_Tp(0)
  {
    label(0) = "(R*T*log(p))";
  }
  virtual ideal_gas* clone() const {return new ideal_gas(*this);}
  virtual double Z(function f, const StateTp &Tp) const;
};

class V_const : public null_Tp
{
public:
  V_const() : null_Tp(1)
  {
    label(0) = "(";
//    name(0) = "V";
    label(1) = "*(p/po))";
  }
  virtual V_const* clone() const {return new V_const(*this);}
  virtual double Z(function f, const StateTp &Tp) const;
};

class alpha_const : public null_Tp
{
public:
  alpha_const() : null_Tp(2)
  {
    label(0) = "(";
//    name(0) = "Vo";
    label(1) = "*exp(";
//    name(1) = "alpha";
    label(2) = "*T)*(p/po-1)";
  }
  virtual alpha_const* clone() const {return new alpha_const(*this);}
  virtual double Z(function f, const StateTp &Tp) const;
};

class alpha_kappa_const : public null_Tp
{
public:
  alpha_kappa_const() : null_Tp(3)
  {
    label(0) = "(";
//    name(0) = "Vb";
    label(1) = "*exp(";
//    name(1) = "alpha";
    label(2) = "*T)*(1.-exp(";
//    name(2) = "kappa";
    label(3) = "*(1-p/po)))";
  }
  virtual alpha_kappa_const* clone() const {return new alpha_kappa_const(*this);}
  virtual double Z(function f, const StateTp &Tp) const;
};

class alpha_kappa_const2 : public null_Tp
{
public:
  alpha_kappa_const2() : null_Tp(3)
  {
    label(0) = "(";
//    name(0) = "V";
    label(1) = "/kappa*exp(";
//    name(1) = "alpha";
    label(2) = "*T)*(1.-exp(";
//    name(2) = "kappa";
    label(3) = "*(1-p/po)))";
  }
  virtual alpha_kappa_const2* clone() const {return new alpha_kappa_const2(*this);}
  virtual double Z(function f, const StateTp &Tp) const;
};

class compound_Tp : public Ref_func_Tp
{
  vec_double limT_;
  vec_double limp_;
  vec_func_Tp fptr_;

public:
  compound_Tp() : limT_(2), limp_(2), fptr_(1)
    {limT_[0] = 0.; limT_[1] = HUGE_VAL;
     limp_[0] = 0.; limp_[1] = HUGE_VAL;}

  virtual void clear()
  {
    vec_double t(2);
    t[0] = 0.;
    t[1] = HUGE_VAL;
    limT_ = t;
    limp_ = t;
    fptr_ = vec_func_Tp(1);
  }
  virtual compound_Tp* clone() const {return new compound_Tp(*this);}

  virtual void read(const SGML &e);
  virtual void WriteBody(ostream& out, size_t shift = 0) const;

  virtual limits limitT() const
    {return limits(limT_.front(), limT_.back());}
  virtual limits limitp() const
    {return limits(limp_.front(), limp_.back());}
  virtual bool IsInside(const StateTp &Tp) const
    {return Tp.T() >= limT_.front() &&
            Tp.T() <= limT_.back() &&
            Tp.p() >= limp_.front() &&
            Tp.p() <= limp_.back() ;}

  const func_Tp& find(const StateTp &Tp) const;
  virtual double Z(function f, const StateTp &Tp) const
    {return find(Tp).Z(f, Tp);}
};

class complex_Tp : public Ref_func_Tp
{
  func_Tp a;
  func_Tp b;

public:
  complex_Tp() {}

  virtual complex_Tp* clone() const {return new complex_Tp(*this);}

  virtual void read(const SGML &e)
  {
    istringstream in(e.body);
    char t;
    in >> t >> a >> t >> b;
  }
  virtual void WriteBody(ostream& out, size_t shift = 0) const
  {
    out << PutTab(shift) << "(" << endl;
    a.write(out, shift);
    out << PutTab(shift) << "+" << endl;
    b.write(out, shift);
    out << PutTab(shift) << ")" << endl;
  }

  virtual limits limitT() const
  {
    return limits(max(a.limitT().lower, b.limitT().lower),
                  min(a.limitT().upper, b.limitT().upper));
  }
  virtual limits limitp() const
  {
    return limits(max(a.limitp().lower, b.limitp().lower),
                  min(a.limitp().upper, b.limitp().upper));
  }
  virtual bool IsInside(const StateTp &Tp) const
  {
    return a.IsInside(Tp) && b.IsInside(Tp);
  }

  virtual double Z(function f, const StateTp &Tp) const
    {return a.Z(f, Tp) + b.Z(f, Tp);}
};

class calc_Tp : public simple_Tp, public CalculatorWithCoefs
{
  mutable double* ptrT;
  mutable double* ptrp;
  virtual void AnalizeId(const string &id);
public:
  calc_Tp() {}

  virtual void clear() {CalculatorWithCoefs::clear();}
  virtual calc_Tp* clone() const {return new calc_Tp(*this);}

  virtual void read(const SGML &e)
  {
    clear();
    FromString(e.body);
  }

  virtual void WriteBody(ostream& out, size_t shift = 0) const
  {
    CalculatorWithCoefs::write(out, shift);
  }

  virtual double Z(function f, const StateTp &Tp) const;
};

#endif
