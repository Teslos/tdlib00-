/*
Copyright (C) 1998, 1999, 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                       http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#include <algorithm>
#include <typeinfo>
#include "phase.h"
#include "func_tpx_imp.h"

static map_string_string *type;

RefPhase *phase::PTR_NULL;
phase::map_types *phase::types = 0;
phase::map_types *phase::obj = 0;
phase::Destruct phase::clean;

vec_double RefPhase::dempty;
vec_string RefPhase::sempty;

bool phase::debug = false;

void RegisterSimple();
void RegisterCuOx_plane();
void RegisterCuOx_ordered_plane();
void RegisterApBq();
void RegisterAssociatedSolution();

void (*phase::init_ptr[])() = {
  RegisterSimple, 
  RegisterCuOx_plane, 
  RegisterCuOx_ordered_plane, 
//  RegisterApBq, 
  RegisterAssociatedSolution, 
  0};

istream& RefPhase::read(istream& in)
{
  phase::init();
  SGML e(in);
  e.compare((*type)[typeid(*this).name()]);
  ReadID(e);
  read(e);
  return in;
}

ostream& RefPhase::write(ostream& out, size_t shift) const
{
  phase::init();
  out << PutTab(shift) << "<" << (*type)[typeid(*this).name()] 
    << " class=phase";
  WriteAttributes(out, shift + 1);
  out << ">" << endl;
  WriteBody(out, shift + 1);
  out << PutTab(shift) << "</" << (*type)[typeid(*this).name()] << "> " << endl;
  return out;
}

void RefPhase::ReadID(SGML &e)
{
// for compatibility with the old format
  parser p(e);
  SGML el;
  p.GetSGML(el);
  if (el.name == "name")
  {
    e.body.erase(0, p.npos());
    parser p2(el);
    e.attr["id"] = p2.GetToken();
  }
}

void RefPhase::ReadComponents(vec_formula &vf, SGML &e, FuncTpx &ref)
{
  parser p(e);
  SGML el;
  p.GetSGML(el);
  ref.clear();
  if (el.name == "components")
  {
    e.body.erase(0, p.npos());
    if (el.body.find("<") != string::npos)
    {
// compatibility  
      el.name = "Reference";
      ref.read(el);
      if (!vf.size())
        vf.resize(ref.size());
      if (ref.size() != vf.size())
        throw 
          gError("RefPhase: number of components in ref is not equal to vf");
      else
        for (size_t i = 0; i < ref.size(); ++i)
          vf[i] = dynamic_cast<const Reference*>(ref.pointer())->
            GetSpecies()[i].GetFormula();
      }
    else
    {
// new format
      istringstream in(el.body);
      if (vf.size() > 0)
        for (size_t i = 0; i < vf.size(); ++i)
          in >> vf[i];
      else
      {
        formula f;
        while (in >> ws, in.peek() != EOF)
        {
          in >> f;
          if (f.cf().empty())
            break;
          vf.push_back(f);
        }
      }
      parser p(e);
      p.GetSGML(el);
      if (el.name == "Reference")
      {
        ref.read(el, &vf);
        e.body.erase(0, p.npos());
      }
    }
    vz.resize(vf.size());
  }
}

struct GT_T2
{
  const RefPhase &f;
  index i;
  StateTp Tp;
  const StateX &x;
  GT_T2(const RefPhase &f_, index i_, const double &p, const StateX &x_) 
    : f(f_), i(i_), x(x_) {Tp.p() = p;}
  double operator()(const double &T)
  {
    Tp.T() = T;
    return f.Z(::G, i, Tp, x)/T;
  }
};

struct G_Tp2
{
  const RefPhase &f;
  index i;
  StateTp Tp;
  const StateX &x;
  G_Tp2(const RefPhase &f_, index i_, const StateX &x_) 
    : f(f_), i(i_), x(x_) {}
  double operator()(const double &T, const double &p)
  {
    Tp.T() = T;
    Tp.p() = p;
    return f.Z(::G, i, Tp, x);
  }
};

double RefPhase::num_H(index i, const StateTp &Tp, const StateX &x) const
{
  GT_T2 num(*this, i, Tp.p(), x);
  return -Tp.T()*Tp.T()*FirstDerivative(num, Tp.T());
}

double RefPhase::num_S(index i, const StateTp &Tp, const StateX &x) const
{
  Z_T num(::G, *this, i, Tp.p(), x);
  return -FirstDerivative(num, Tp.T());
}

double RefPhase::num_Cp(index i, const StateTp &Tp, const StateX &x) const
{
  Z_T num(::G, *this, i, Tp.p(), x);
  return -Tp.T()*SecondDerivative(num, Tp.T());
}

double RefPhase::num_V(index i, const StateTp &Tp, const StateX &x) const
{
  Z_p num(::G, *this, i, Tp.T(), x);
  return FirstDerivative(num, Tp.p());
}

double RefPhase::num_dVdT(index i, const StateTp &Tp, const StateX &x) const
{
  G_Tp2 num(*this, i, x);
  return MixedDerivative(num, Tp.T(), Tp.p());
}

double RefPhase::num_dVdp(index i, const StateTp &Tp, const StateX &x) const
{
  Z_p num(::G, *this, i, Tp.T(), x);
  return SecondDerivative(num, Tp.p());
}

struct GTA_T
{
  const RefPhase &f;
  index i;
  StateTp Tp;
  GTA_T(const RefPhase &f_, index i_, const double &p) 
    : f(f_), i(i_) {Tp.p() = p;}
  double operator()(const StateX &x, const double &T)
  {
    Tp.T() = T;
    return f.Z(::G, i, Tp, x)/T;
  }
};

struct GA_T
{
  const RefPhase &f;
  index i;
  StateTp Tp;
  GA_T(const RefPhase &f_, index i_, const double &p) 
    : f(f_), i(i_) {Tp.p() = p;}
  double operator()(const StateX &x, const double &T)
  {
    Tp.T() = T;
    return f.Z(::G, i, Tp, x);
  }
};

struct GA_p
{
  const RefPhase &f;
  index i;
  StateTp Tp;
  GA_p(const RefPhase &f_, index i_, const double &T) 
    : f(f_), i(i_) {Tp.T() = T;}
  double operator()(const StateX &x, const double &p)
  {
    Tp.p() = p;
    return f.Z(::G, i, Tp, x);
  }
};

const vec_double& RefPhase::num_dGdx(index i, const StateTp &Tp,
                                     const StateX &x) const
{
  ZA_x num(::G, *this, i, Tp);
  ArrayFirstDerivative(num, x, vz);
  return vz;
}

const vec_double& RefPhase::num_dHdx(index i, const StateTp &Tp,
                                     const StateX &x) const
{
  GTA_T num(*this, i, Tp.p());
  ArrayMixedDerivative(num, x, Tp.T(), vz, -Tp.T()*Tp.T());
  return vz;
}

const vec_double& RefPhase::num_dSdx(index i, const StateTp &Tp,
                                     const StateX &x) const
{
  GA_T num(*this, i, Tp.p());
  ArrayMixedDerivative(num, x, Tp.T(), vz, -1.);
  return vz;
}

const vec_double& RefPhase::num_dVdx(index i, const StateTp &Tp,
                                     const StateX &x) const
{
  GA_p num(*this, i, Tp.T());
  ArrayMixedDerivative(num, x, Tp.p(), vz);
  return vz;
}

const vec_double& RefPhase::dZdx(function f, index i, const StateTp &Tp,
                               const StateX &x) const
{
  if (size() == 0)
    return vz;
  if (size() == 1)
  {
    vz[0] = 0.;
    return vz;
  }
  fill(vz.begin(), vz.end(), 0.);
  switch (f) // result in vz
  {
    case ::G:
      num_dGdx(i, Tp, x);
      break;
    case ::H:
      num_dHdx(i, Tp, x);
      break;
    case ::S:
      num_dSdx(i, Tp, x);
      break;
    case ::V:
      num_dVdx(i, Tp, x);
      break;
    default:
      throw gError("RefPhase: function is not supported");
  }
  return vz;
}

const vec_double& RefPhase::z(function f, index i, const StateTp &Tp,
                               const StateX &x) const
{
  if (size() == 0)
    return vz;
  if (size() == 1)
  {
    vz[0] = Z(f, i, Tp, x);
    return vz;
  }
  dZdx(f, i, Tp, x); // result in vz
  size_t j;
  double Sum = Z(f, i, Tp, x);
  for (j = 0; j != size(); ++j)
    Sum -= vz[j]*x[j];
  for (j = 0; j != size(); ++j)
    vz[j] += Sum;
  return vz;
}

const vec_double& RefPhase::z_y(function f, index i, const StateTp &Tp, 
    const StateX &x, const vec_double &y) const
{
  if (size() == 0)
    return vz;
  if (size() == 1)
  {
    vz[0] = Z_y(f, i, Tp, x, y);
    return vz;
  }
  dZdx_y(f, i, Tp, x, y); // result in vz
  size_t j;
  double Sum = Z_y(f, i, Tp, x, y);
  for (j = 0; j != size(); ++j)
    Sum -= vz[j]*x[j];
  for (j = 0; j != size(); ++j)
    vz[j] += Sum;
  return vz;
}

phase::Destruct::~Destruct()
{
  if (phase::types)
  {
    delete phase::types;
    delete phase::obj;
    phase::types = 0;
    phase::obj = 0;
    delete type;
  }
}

void phase::init()
{
  if (!types)
  {
    type = new map_string_string;
    PTR_NULL = new NullPhase;
    types = new map_types;
    obj = new map_types;
    RegisterType("NullPhase", PTR_NULL);
    RegisterType("point_phase", new PointPhase);
    RegisterType("PointPhase", new PointPhase);
    RegisterType("NumericalDerivatives", new NumericalDerivatives);
    for (size_t i = 0; init_ptr[i]; ++i)
      (*init_ptr[i])();
  }
}

void phase::RegisterType(const string& name, RefPhase* ptr)
{
  init();
  (*types)[name] = phase(ptr);
  (*type)[typeid(*ptr).name()] = name;
}

void phase::find(const string &name)
{
  map_types_i i = obj->find(name);
  if (i == obj->end())
  {
    species sps;
    try
    {
      sps.find(name);
      create("PointPhase");
      cow();
      dynamic_cast<PointPhase*>(ptr)->assign(sps);
    }
    catch (gError)
    {
      throw gError(string("phase: id ").append(name).append
          (" is not defined"));
    }
  }
  else
    *this = (*i).second;
}

void phase::read(const SGML &e_)
{
  string id;
  if (!(id = e_.FindString("IDREF")).empty())
  {
    find(id);
  }
  else
  {
    create(e_.name);
    cow();
    SGML e = e_;
    RefPhase::ReadID(e);
    try
    {
      ptr->read(e);
    }
    catch (ObjFromProxy &obj)
    {
      *this = phase(obj.obj);
    }
    if (!(ptr->id = e.FindString("id")).empty())
    {
      pair<map_types_i, bool> i;
      i = obj->insert(map_types::value_type(ptr->id, *this));
      if (!i.second)
      {
        cout << "phase: " << ptr->id << " is already defined - "
          << "set to anonymous" << endl;
        ptr->id.erase();
      }
    }
  }
}

ostream& phase::write(ostream& out, size_t shift) const
{
  if (ptr->id.empty() || ptr->print)
  {
    out << PutTab(shift) << "<" << (*type)[typeid(*ptr).name()] 
      << " class=phase";
    if (!ptr->id.empty())
    {
      out << " id=" << ptr->id;
    }
    ptr->WriteAttributes(out, shift + 1);
    out << ">" << endl;
    ptr->WriteBody(out, shift + 1);
    out << PutTab(shift) << "</" << (*type)[typeid(*ptr).name()] << "> " 
      << endl;
    ptr->print = false;
  }
  else
  {
    out << PutTab(shift) << "<phase class=phase IDREF="
      << ptr->id << "></phase>"
      << endl;
  }
  return out;
}

void RefSimple::WriteComponents(ostream &out, size_t shift) const
{
  out << PutTab(shift) << "<components> ";
  for (vec_formula_ci i = vf.begin(); i != vf.end(); ++i)
    out << (*i) << " ";
  out << "</components>" << endl;
}


void PointPhase::read(const SGML &e)
{
  parser p(e);
  SGML el;
  p.GetSGML(el);
  if (el.name == "components")
  {
    parser p2(el);
    SGML elem;
    p2.GetSGML(elem);
    sps.read(elem);
  }
  else
    sps.read(el);
  vf[0] = sps.GetFormula();
}

void PointPhase::WriteBody(ostream &out, size_t shift) const
{
  sps.write(out, shift);
}

double NumericalDerivatives::Z(function f, index i, const StateTp &Tp, 
    const StateX &x) const
{
  switch (f)
  {
    case ::G:
      return ph.Z(f, i, Tp, x);
    case ::H:
      return num_H(i, Tp, x);
    case ::S:
      return num_S(i, Tp, x);
    case ::Cp:
      return num_Cp(i, Tp, x);
    case ::V:
      return num_V(i, Tp, x);
    case ::dVdT:
      return num_dVdT(i, Tp, x);
    case ::dVdp:
      return num_dVdp(i, Tp, x);
    default:
      throw gError("NumericalDerivatives: function is not supported");
  }
}

