/*
Copyright (C) 1998, 1999, 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                       http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#ifndef __PHASE_H
#define __PHASE_H

#include <func_tpx.h>

template <class F>
inline void ArrayFirstDerivative(F f, const StateX &x, vec_double &res,
    const double step = global::step1)
{
  double oldx;
  StateX &xx = const_cast<StateX &>(x);
  for (size_t i = 0; i < x.size(); ++i)
  {
    oldx = xx[i];
    double dx = step*(xx[i] + step);
    xx[i] += dx;
    double f2 = f(xx);
    xx[i] = oldx - dx;
    double f1 = f(xx);
    res[i] = (f2 - f1)/dx/2.0;
    xx[i] = oldx;
  }
}

template <class F>
inline void ArrayMixedDerivative(F f, const StateX &x, const double &y, 
    vec_double &res, const double &mul = 1., const double &step = global::step2)
{
  double oldx;
  StateX &xx = const_cast<StateX &>(x);
  double dy = step*(y + step);
  for (size_t i = 0; i < x.size(); ++i)
  {
    oldx = xx[i];
    double dx = step*(xx[i] + step);
    xx[i] -= dx;
    double f1 = f(x, y - dy);
    double f2 = f(x, y + dy);
    xx[i] = oldx + dx;
    double f3 = f(x, y - dy);
    double f4 = f(x, y + dy);
    res[i] = mul*((f4 - f3) - (f2 - f1))/dy/dx/4.0;
    xx[i] = oldx;
  }
}

class RefPhase
{
  friend class phase;
  size_t Refs;

  void AddReference()
    {++Refs;}
  size_t References()
    {return Refs;}
  size_t RemoveReference()
    {return --Refs;}

  string id;
  bool print;
  
  static vec_double dempty;
  static vec_string sempty;
protected:
  mutable vec_double vz;

public:
  RefPhase() : Refs(1), print(true) {}
  RefPhase(const RefPhase &old)
    : Refs(1), id(old.id), print(old.print), vz(old.vz.size()) {}
  explicit RefPhase(size_t n) : Refs(1) {clear(n);}
  virtual ~RefPhase() {}

  virtual void clear(size_t n = 0, size_t m = 0, size_t l = 0)
  {
    vz.resize(n);
    print = true;
  }

  virtual RefPhase* clone() const = 0;

  istream& read(istream& in);
  ostream& write(ostream& out, size_t shift = 0) const;
  virtual void read(const SGML &e) = 0;
  virtual void WriteAttributes(ostream &out, size_t shift = 0) const {}
  virtual void WriteBody(ostream &out, size_t shift = 0) const = 0;
  void static ReadID(SGML &e);
  void ReadComponents(vec_formula &vf, SGML &e, FuncTpx &ref);

  virtual limits limitT() {return limits(0., HUGE_VAL);}
  virtual limits limitp() {return limits(0., HUGE_VAL);}
  virtual bool IsInside(const StateTp &Tp)
    {return Tp.T() >= 0. && Tp.p() >= 0.;}

  virtual const vec_formula& comps() const = 0;
  size_t size() const {return comps().size();}

  virtual double Z(function f, index i, const StateTp &Tp, 
      const StateX &x) const = 0;
  double G(const StateTp &Tp, const StateX &x) const
    {return Z(::G, ::full, Tp, x);}

  double num_H(index i, const StateTp &Tp, const StateX &x) const;
  double num_S(index i, const StateTp &Tp, const StateX &x) const;
  double num_Cp(index i, const StateTp &Tp, const StateX &x) const;
  double num_V(index i, const StateTp &Tp, const StateX &x) const;
  double num_dVdT(index i, const StateTp &Tp, const StateX &x) const;
  double num_dVdp(index i, const StateTp &Tp, const StateX &x) const;

  virtual const vec_double& dZdx(function f, index i, const StateTp &Tp, 
      const StateX &x) const;
  virtual const vec_double& z(function f, index i, const StateTp &Tp,
      const StateX &x) const;
  const vec_double& mu(const StateTp &Tp, const StateX &x) const
    {return z(::G, ::full, Tp, x);}

  const vec_double& num_dGdx(index i, const StateTp &Tp, const StateX &x) const;
  const vec_double& num_dHdx(index i, const StateTp &Tp, const StateX &x) const;
  const vec_double& num_dSdx(index i, const StateTp &Tp, const StateX &x) const;
  const vec_double& num_dVdx(index i, const StateTp &Tp, const StateX &x) const;

  virtual size_t NIntVars() const {return 0;}
  virtual const vec_string& NamesIntVars() const {return sempty;}
  virtual const vec_double& IntVarsEq(const StateTp &Tp, const StateX &x) const 
    {return dempty;}
  virtual const vec_double& IntVarsFromGuess(const StateTp &Tp, const StateX &x,
      const vec_double& y) const 
    {return IntVarsEq(Tp, x);}
  virtual double Z_y(function f, index i, const StateTp &Tp,
      const StateX &x, const vec_double &y) const {return Z(f, i, Tp, x);}
  virtual const vec_double& dZdx_y(function f, index i, const StateTp &Tp, 
      const StateX &x, const vec_double &y) const {return dZdx(f, i, Tp, x);}
  virtual const vec_double& z_y(function f, index i, const StateTp &Tp,
      const StateX &x, const vec_double &y) const;

  virtual size_t NIntPhases() const {return 0;}
  virtual const vec_string& NamesIntPhases() const {return sempty;}
  virtual bool IsStable(size_t i, const StateTp &Tp, 
      const StateX &x) const {return true;}
  virtual const string& StableIntPhase(const StateTp &Tp, 
      const StateX &x) const {return id;}
  virtual double CheckBoundary(size_t i, size_t j, const StateTp &Tp, 
      const StateX &x) const {return 0.;};
};

OVERLOAD_STREAMS(RefPhase)

class ObjFromProxy : public gError
{
public:
  RefPhase *obj;
  ObjFromProxy(RefPhase *obj_) : gError("ObjFromProxy"), obj(obj_) {}
};

struct Z_T
{
  const RefPhase &f;
  index i;
  StateTp Tp;
  const StateX &x;
  function fun;
  Z_T(function fun_, const RefPhase &f_, index i_, const double &p, 
      const StateX &x_) 
    : f(f_), i(i_), x(x_), fun(fun_) {Tp.p() = p;}
  double operator()(const double &T)
  {
    Tp.T() = T;
    return f.Z(fun, i, Tp, x);
  }
};

struct Z_p
{
  const RefPhase &f;
  index i;
  StateTp Tp;
  const StateX &x;
  function fun;
  Z_p(function fun_, const RefPhase &f_, index i_, const double &T, 
      const StateX &x_) 
    : f(f_), i(i_), x(x_), fun(fun_) {Tp.T() = T;}
  double operator()(const double &p)
  {
    Tp.p() = p;
    return f.Z(fun, i, Tp, x);
  }
};

struct ZA_x
{
  const RefPhase &f;
  index i;
  const StateTp &Tp;
  function fun;
  ZA_x(function fun_, const RefPhase &f_, index i_, const StateTp &Tp_)
    : f(f_), i(i_), Tp(Tp_), fun(fun_) {}
  double operator()(const StateX &x)
  {
    return f.Z(fun, i, Tp, x);
  }
};

class phase
{
  typedef map<string, phase, less<string> > map_types;
  typedef map_types::iterator map_types_i;
  typedef map_types::const_iterator map_types_ci;

  class Destruct;
  friend class phase::Destruct;
  class Destruct
  {
  public:
    ~Destruct();
  };
  static Destruct clean;

  static RefPhase *PTR_NULL;
  static map_types *types;
  static void (*init_ptr[])();
  static map_types *obj;

  friend istream& RefPhase::read(istream& in);
  friend ostream& RefPhase::write(ostream& out, size_t shift = 0) const;
  static void init();

  RefPhase *ptr;

  void cow()
  {
    if (ptr->References() > 1)
    {
      ptr->RemoveReference();
      ptr = ptr->clone();
    }
  }
  phase(RefPhase* ptr) : ptr(ptr) {}

public:
  static bool debug;
  static void RegisterType(const string& name, RefPhase* ptr);
  static void WriteAll(ostream &out)
  {
    SetNonPrinted();
    if (obj)
      for (map_types_i i = obj->begin(); i != obj->end(); ++i)
        (*i).second.write(out);
  }
  static void SetPrinted()
  {
    if (obj)
      for (map_types_i i = obj->begin(); i != obj->end(); ++i)
        (*i).second.ptr->print = false;
  }
  static void SetNonPrinted()
  {
    if (obj)
      for (map_types_i i = obj->begin(); i != obj->end(); ++i)
        (*i).second.ptr->print = true;
  }

  phase()
  {
    init();
    ptr = PTR_NULL;
    ptr->AddReference();
  }
  phase(const phase& old)
  {
    ptr = old.ptr;
    ptr->AddReference();
  }
  phase(const string &name)
  {
    init();
    map_types_ci i = types->find(name);
    if (i == types->end())
      throw gError(string("phase: unknown type ") + name);
    else
    {
      ptr = (*i).second.ptr;
      ptr->AddReference();
    }
  }
  ~phase()
  {
    if (ptr->RemoveReference() == 0)
      delete ptr;
  }
  phase& operator=(const phase &old)
  {
    if (this != &old)
    {
      if (ptr->RemoveReference() == 0)
        delete ptr;
      ptr = old.ptr;
      ptr->AddReference();
    }
    return *this;
  }
  void create(const string &name)
  {
    map_types_ci i = types->find(name);
    if (i == types->end())
      throw gError(string("phase: unknown type ") + name);
    else
      *this = (*i).second;
  }
  void find(const string &name);
  void clear()
  {
    cow();
    ptr->clear();
  }
  const RefPhase* pointer() const
    {return ptr;}

  istream& read(istream& in)
  {
    SGML e(in);
    read(e);
    return in;
  }
  void read(const SGML &e);
  ostream& write(ostream& out, size_t shift = 0) const;

  limits limitT() const
    {return ptr->limitT();}
  limits limitp() const
    {return ptr->limitp();}
  bool IsInside(const StateTp &Tp) const
    {return ptr->IsInside(Tp);}

  const string& id() const
    {return ptr->id;}
  const vec_formula& comps() const
    {return ptr->comps();}
  size_t size() const
    {return ptr->size();}
  
  double Z(function f, index i, const StateTp &Tp, const StateX &x) const
    {return ptr->Z(f, i, Tp, x);}
  double Z(function f, const StateTp &Tp, const StateX &x) const
    {return ptr->Z(f, ::full, Tp, x);}
  double G(const StateTp &Tp, const StateX &x) const
    {return ptr->Z(::G, ::full, Tp, x);}
  
  const vec_double& z(function f, index i, const StateTp &Tp,
                      const StateX &x) const
    {return ptr->z(f, i, Tp, x);}
  const vec_double& z(function f, const StateTp &Tp, const StateX &x) const
    {return ptr->z(f, ::full, Tp, x);}
  const vec_double& mu(const StateTp &Tp, const StateX &x) const
    {return ptr->z(::G, ::full, Tp, x);}

  size_t NIntVars() const 
    {return ptr->NIntVars();}
  const vec_string& NamesIntVars() const 
    {return ptr->NamesIntVars();}
  const vec_double& IntVarsEq(const StateTp &Tp, const StateX &x) const 
    {return ptr->IntVarsEq(Tp, x);}
  const vec_double& IntVarsEq(const StateTp &Tp, const StateX &x, 
      const vec_double& y) const 
    {return ptr->IntVarsFromGuess(Tp, x, y);}
  
  double Z(function f, index i, const StateTp &Tp,
      const StateX &x, const vec_double &y) const 
    {return ptr->Z_y(f, i, Tp, x, y);}
  const vec_double& dZdx(function f, index i, const StateTp &Tp, 
      const StateX &x, const vec_double &y) const 
    {return ptr->dZdx_y(f, i, Tp, x, y);}
  const vec_double& z(function f, index i, const StateTp &Tp,
      const StateX &x, const vec_double &y) const
    {return ptr->z_y(f, i, Tp, x, y);}

  size_t NIntPhases() const 
    {return ptr->NIntPhases();}
  const vec_string& NamesIntPhases() const 
    {return ptr->NamesIntPhases();}
  bool IsStable(size_t i, const StateTp &Tp, const StateX &x) const 
    {return ptr->IsStable(i, Tp, x);}
  const string& StableIntPhase(const StateTp &Tp, const StateX &x) const 
    {return ptr->StableIntPhase(Tp, x);}
  double CheckBoundary(size_t i, size_t j, const StateTp &Tp, 
      const StateX &x) const
    {return ptr->CheckBoundary(i, j, Tp, x);}
};

OVERLOAD_STREAMS(phase)

typedef vector<phase> vec_phase;
typedef vec_phase::iterator vec_phase_i;
typedef vec_phase::const_iterator vec_phase_ci;

class RefSimple : public RefPhase
{
protected:
  vec_formula vf;

public:
  RefSimple() {}
  explicit RefSimple(size_t n) {clear(n);}

  virtual void clear(size_t n = 0, size_t m = 0, size_t l = 0)
  {
    vf.resize(n);
    RefPhase::clear(n);
  }

  void WriteComponents(ostream &out, size_t shift = 0) const;

  virtual const vec_formula& comps() const
    {return vf;}
};

class NullPhase : public RefSimple
{
public:
  NullPhase() {}

  virtual NullPhase* clone() const {return new NullPhase(*this);}

  virtual void read(const SGML &e) {}
  virtual void WriteBody(ostream &out, size_t shift = 0) const {}

  virtual double Z(function f, index i, const StateTp &Tp,
                   const StateX &x) const
    {return 0.;}
};

class PointPhase : public RefSimple
{
  species sps;
public:
  PointPhase() : RefSimple(1) {}

  virtual PointPhase* clone() const {return new PointPhase(*this);}

  virtual void clear(size_t n = 0, size_t m = 0, size_t l = 0)
    {RefSimple::clear(1);}

  virtual void read(const SGML &e);
  virtual void WriteBody(ostream &out, size_t shift = 0) const;

  void assign(const species &s)
  {
    sps = s;
    vf[0] = sps.GetFormula();
  }

  virtual double Z(function f, index i, const StateTp &Tp,
                   const StateX &x) const
    {return sps.Z(f, i, Tp);}
};

class NumericalDerivatives : public RefPhase
{
  phase ph;
public:
  virtual NumericalDerivatives* clone() const 
    {return new NumericalDerivatives(*this);}

  virtual void read(const SGML &e)
  {
    parser p(e);
    SGML el;
    ph.read(p.GetSGML(el));
    clear(ph.size());
  }
  
  virtual void WriteBody(ostream &out, size_t shift = 0) const
  {
    ph.write(out, shift);
  }
 
  virtual const vec_formula& comps() const {return ph.comps();}

  virtual double Z(function f, index i, const StateTp &Tp, 
      const StateX &x) const;
  virtual const vec_double& dZdx(function f, index i, const StateTp &Tp, 
      const StateX &x) const {return RefPhase::dZdx(f, i, Tp, x);}
};

#endif
