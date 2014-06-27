/*
Copyright (C) 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#ifndef __FUNC_TPX_H
#define __FUNC_TPX_H

#include <species.h>

size_t bin_coef(size_t m, size_t n);

class SymMatrix : public vec_double
{
  size_t n;

public:
  SymMatrix() : n(0) {}
  explicit SymMatrix(size_t n_) : vec_double(n_*(n_ + 1)/2), n(n_) {}

  void resize(size_t n_)
    {n = n_; vec_double::resize(n*(n + 1)/2);}
  void reserve(size_t n_)
    {vec_double::reserve(n_*(n_ + 1)/2);}
  void clear()
    {n = 0; vec_double::clear();}

  double& operator()(size_t i, size_t j)
  {
    return operator[]((i < j) ? i + (j*j + j)/2 : j + (i*i + i)/2);
  }
  const double& operator()(size_t i, size_t j) const
  {
    return operator[]((i < j) ? i + (j*j + j)/2 : j + (i*i + i)/2);
  }
};

enum Formalism {NoChange, Kohler, Muggianu};
extern const char *FormalismName[];
typedef Enum<Formalism, FormalismName> formalism;

class FuncTpx;

class RefFuncTpx
{
private:
  friend class FuncTpx;
  size_t Refs;

  void AddReference()
    {++Refs;}
  size_t References()
    {return Refs;}
  size_t RemoveReference()
    {return --Refs;}

  string id;
  bool print;

  index idx;
  
protected:
  vec_int vars;
  void SetIndex(index i) {idx = i;}
  void SetVars(const vec_int& v) {vars = v;}

  void WriteFormula(ostream& out, size_t i, const vec_formula *f) const
  {
    if (f)
      out << (*f)[vars[i]];
    else
      out << vars[i] + 1;
  }
  static size_t FindFormula(const string  &id, const vec_formula *f)
  {
    int num = -2;
    if (f)
      for (size_t i = 0; i < (*f).size(); ++i)
			{
        if (id == (*f)[i].cf())
          num = i;
			}
    if (num == -2)
      throw gError(string("func_Tpx: no formula ").append(id));
    return num;
  }
  static size_t FindFormula(parser &p, const vec_formula *f)
  {
    int num = -2;
    if (!p.SkipChar('x'))
      p.SkipChar('v');
    p.SkipChar('(');
    num = p.GetInt(num) - 1;
    if (num < 0) 
    {
      p.GetToken();
      num = FindFormula(p.token, f);
    }
    p.SkipChar(')');
    return num;
  }
  
public:
  RefFuncTpx() : Refs(1), print(true), idx((Index)0) {}
  RefFuncTpx(const RefFuncTpx &old)
    : Refs(1), print(old.print), idx(old.idx), vars(old.vars){}
  explicit RefFuncTpx(size_t n) : Refs(1), print(true), idx((Index)0), 
    vars(n) {}
  virtual ~RefFuncTpx() {}

  virtual void clear(size_t i = 0)
  {
    vars.resize(i);
    print = true;
  }
  virtual RefFuncTpx* clone() const = 0;

  istream& read(istream& in);
  ostream& write(ostream& out, size_t shift = 0) const;
  virtual void read(const SGML &e, const vec_formula *f = 0) = 0;
  virtual void WriteAttributes(ostream &out, const vec_formula *f = 0, 
      size_t shift = 0) const {} 
  virtual void WriteBody(ostream &out, const vec_formula *f = 0, 
      size_t shift = 0) const = 0;

  virtual limits limitT() {return limits(0., HUGE_VAL);}
  virtual limits limitp() {return limits(0., HUGE_VAL);}
  virtual bool IsInside(const StateTp &Tp)
    {return Tp.T() >= 0. && Tp.p() >= 0.;}

  virtual formalism SetFormalism(formalism frm) {return ::NoChange;}
  virtual void set(const vec_int &vec, const func_Tp &f = func_Tp("Cp_zero"), 
      size_t power = 0) {}
  
  index GetIndex() const {return idx;}
  size_t size() const {return vars.size();}
  const vec_int& GetVars() const {return vars;}

  virtual double Z(function f, const StateTp &Tp, const StateX &x) const = 0;

// ATTENTION!! these functions should add the result to res
// Do not overwrite res
  virtual void dZdx(function f, const StateTp &Tp,
                    const StateX &x, vec_double &res) const = 0;
  virtual void d2Zdx2(function f, const StateTp &Tp,
                    const StateX &x, SymMatrix &res) const = 0;
};

OVERLOAD_STREAMS(RefFuncTpx)

class FuncTpx
{
  typedef map<string, FuncTpx, less<string> > map_types;
  typedef map_types::iterator map_types_i;
  typedef map_types::const_iterator map_types_ci;

  class Destruct;
  friend class FuncTpx::Destruct;
  class Destruct
  {
  public:
    ~Destruct();
  };
  static Destruct clean;

  static RefFuncTpx *PTR_NULL;
  static map_types *types;
  static void (*init_ptr[])();
  static map_types *obj;

  friend istream& RefFuncTpx::read(istream& in);
  friend ostream& RefFuncTpx::write(ostream& out, size_t shift = 0) const;
  static void init();

  RefFuncTpx *ptr;

  void cow()
  {
    if (ptr->References() > 1)
    {
      ptr->RemoveReference();
      ptr = ptr->clone();
    }
  }
  FuncTpx(RefFuncTpx* ptr) : ptr(ptr) {}

public:
  static void RegisterType(const string& name, RefFuncTpx* ptr);
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

  FuncTpx()
  {
    init();
    ptr = PTR_NULL;
    ptr->AddReference();
  }
  FuncTpx(const FuncTpx& old)
  {
    ptr = old.ptr;
    ptr->AddReference();
  }
  FuncTpx(const string &name)
  {
    init();
    map_types_ci i = types->find(name);
    if (i == types->end())
      throw gError(string("FuncTpx: unknown type ") + name);
    else
    {
      ptr = (*i).second.ptr;
      ptr->AddReference();
    }
  }
  ~FuncTpx()
  {
    if (ptr->RemoveReference() == 0)
      delete ptr;
  }
  FuncTpx& operator=(const FuncTpx &old)
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
      throw gError(string("FuncTpx: unknown type ") + name);
    else
      *this = (*i).second;
  }
  void find(const string &name)
  {
    map_types_i i = obj->find(name);
    if (i == obj->end())
      throw gError(string("FuncTpx: id ").append(ptr->id).append
          (" is not defined"));
    else
      *this = (*i).second;
  }
  void clear()
  {
    cow();
    ptr->clear();
  }
  const RefFuncTpx* pointer() const
    {return ptr;}

  istream& read(istream& in)
  {
    SGML e(in);
    read(e);
    return in;
  }
  void read(const SGML &e, const vec_formula *f = 0);
  ostream& write(ostream& out, const vec_formula *f = 0, 
      size_t shift = 0) const;

  limits limitT() const
    {return ptr->limitT();}
  limits limitp() const
    {return ptr->limitp();}
  bool IsInside(const StateTp &Tp) const
    {return ptr->IsInside(Tp);}

  formalism SetFormalism(formalism frm) 
  {
    cow();
    return ptr->SetFormalism(frm);
  }
  void set(const vec_int &vec, const func_Tp &f = func_Tp("Cp_zero"), 
      size_t power = 0)
  {
    cow();
    return ptr->set(vec, f, power);
  }

  index GetIndex() const {return ptr->idx;}
  size_t size() const {return ptr->vars.size();}
  const vec_int& GetVars() const {return ptr->vars;}

  double Z(function f, const StateTp &Tp, const StateX &x) const
  {
    return ptr->Z(f, Tp, x);
  }
  void dZdx(function f, const StateTp &Tp,
            const StateX &x, vec_double &res) const
  {
    ptr->dZdx(f, Tp, x, res);
  }
  void d2Zdx2(function f, const StateTp &Tp,
            const StateX &x, SymMatrix &res) const
  {
    ptr->d2Zdx2(f, Tp, x, res);
  }
};

OVERLOAD_STREAMS(FuncTpx)
  
typedef vector<FuncTpx> vec_FuncTpx;
typedef vec_FuncTpx::iterator vec_FuncTpx_i;
typedef vec_FuncTpx::const_iterator vec_FuncTpx_ci;

#endif

