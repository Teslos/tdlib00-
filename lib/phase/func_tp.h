/*
Copyright (C) 1998, 1999, 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                       http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#ifndef __FUNC_TP_H
#define __FUNC_TP_H

#include <coef.h>

class func_Tp;

class Ref_func_Tp
{
  friend class func_Tp;
  size_t Refs;
  string id;
  bool print;

  void AddReference()
    {++Refs;}
  size_t References()
    {return Refs;}
  size_t RemoveReference()
    {return --Refs;}
public:

  Ref_func_Tp() : Refs(1), print(true) {}
  Ref_func_Tp(const Ref_func_Tp &old) : Refs(1), print(true) {}
  virtual ~Ref_func_Tp() {}
  Ref_func_Tp& operator=(const Ref_func_Tp &old)
  {
    if (this != &old)
    {
      Refs = 1;
      print = old.print;
    }
    return *this;
  }

  virtual void clear() {}
  virtual Ref_func_Tp* clone() const = 0;

  istream& read(istream& in);
  ostream& write(ostream& out, size_t shift = 0) const;
  virtual void read(const SGML &e) = 0;
  virtual void WriteAttributes(ostream &out, size_t shift = 0) const {}
  virtual void WriteBody(ostream &out, size_t shift = 0) const = 0;
  
  virtual limits limitT() const
    {return limits(0., HUGE_VAL);}
  virtual limits limitp() const
    {return limits(0., HUGE_VAL);}
  virtual bool IsInside(const StateTp &Tp) const
    {return Tp.T() >= 0. && Tp.p() >= 0.;}

  virtual double Z(function f, const StateTp &Tp) const = 0;
  double G(const StateTp &Tp) const
    {return Z(::G, Tp);}
  double num_H(const StateTp &Tp) const;
  double num_S(const StateTp &Tp) const;
  double num_Cp(const StateTp &Tp) const;
  double num_V(const StateTp &Tp) const;
  double num_dVdT(const StateTp &Tp) const;
  double num_dVdp(const StateTp &Tp) const;
};

OVERLOAD_STREAMS(Ref_func_Tp)

class func_Tp
{
  typedef map<string, func_Tp, less<string> > map_types;
  typedef map_types::iterator map_types_i;
  typedef map_types::const_iterator map_types_ci;

  class Destruct;
  friend class func_Tp::Destruct;
  class Destruct
  {
  public:
    ~Destruct();
  };
  static Destruct clean;
  
  static Ref_func_Tp *PTR_NULL;
  static map_types *types;
  static void (*init_ptr[])();
  static map_types *obj;

  friend istream& Ref_func_Tp::read(istream& in);
  friend ostream& Ref_func_Tp::write(ostream& out, size_t shift = 0) const;
  static void init();

  Ref_func_Tp *ptr;

  void cow()
  {
    if (ptr->References() > 1)
    {
      ptr->RemoveReference();
      ptr = ptr->clone();
    }
  }
  func_Tp(Ref_func_Tp *ptr) : ptr(ptr) {}

public: 
  static void RegisterType(const string& name, Ref_func_Tp* ptr);
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

  func_Tp()
  {
    init();
    ptr = PTR_NULL;
    ptr->AddReference();
  }
  func_Tp(const func_Tp& old)
  {
    ptr = old.ptr;
    ptr->AddReference();
  }
  func_Tp(const string &name)
  {
    init();
    map_types_ci i = types->find(name);
    if (i == types->end())
      throw gError(string("func_Tp: unknown type ") + name);
    else
    {
      ptr = (*i).second.ptr;
      ptr->AddReference();
    }
  }
  ~func_Tp()
  {
    if (ptr->RemoveReference() == 0)
      delete ptr;
  }
  func_Tp& operator=(const func_Tp& old)
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
  void create(const string& name)
  {
    map_types_ci i = types->find(name);
    if (i == types->end())
      throw gError(string("func_Tp: unknown type ") + name);
    else
      *this = (*i).second;
  }
  void find(const string& name)
  {
    map_types_i i = obj->find(name);
    if (i == obj->end())
      throw gError(string("func_Tp: id ").append(ptr->id).append
          (" is not defined"));
    else
      *this = (*i).second;
  }
  void clear()
  {
    cow();
    ptr->clear();
  }
  const Ref_func_Tp* pointer() const
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

  double Z(function f, const StateTp &Tp) const
    {return ptr->Z(f, Tp);}
  double G(const StateTp &Tp) const
    {return ptr->Z(::G, Tp);}
};

OVERLOAD_STREAMS(func_Tp)

typedef vector<func_Tp> vec_func_Tp;
typedef vec_func_Tp::iterator vec_func_Tp_i;
typedef vec_func_Tp::const_iterator vec_func_Tp_ci;

class simple_Tp: public Ref_func_Tp
{
public:
  simple_Tp() {}
};

#endif
