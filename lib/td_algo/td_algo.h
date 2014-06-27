/*
Copyright (C) 1998, 1999, 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                       http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#ifndef __TD_ALGO_H
#define __TD_ALGO_H

#include <phase.h>
#include <fstream>

inline int SearchString(const vec_string &vs, const string &str)
{
  for (size_t i = 0; i < vs.size(); ++i)
    if (vs[i] == str)
      return i;
  return -1;
}

inline int SearchFormula(const vec_formula &vs, const string &str)
{
  for (size_t i = 0; i < vs.size(); ++i)
    if (vs[i].cf() == str)
      return i;
  return -1;
}

class convert : public CalculatorWithCoefs
{
  struct var
  {
    double *ptr;
    int i;
    var() : i(0) {}
  };
  typedef map<string, var, less<string> > varmap;
  typedef varmap::iterator varmap_i;
  typedef varmap::const_iterator varmap_ci;
  mutable varmap v;
    
  string name_;
  virtual void AnalizeId(const string &id);

public:
  convert() {}
  convert(const convert &old) : v(old.v), name_(old.name_)
  {
    FromString(old.body());
  }
  convert& operator=(const convert &old)
  {
    if (this != &old)
    {
      v = old.v;
      name_ = old.name_;
      CalculatorWithCoefs::clear();
      FromString(old.body());
    }
    return *this;
  }
  
  istream& read(istream& in)
  {
    SGML e(in);
    read(e);
    return in;
  }
  void read(const SGML &e)
  {
    e.compare("convert");
    name_ = e.FindString("name", "");
    v.clear();
    CalculatorWithCoefs::clear();
    FromString(e.body);
  }
  ostream& write(ostream& out, size_t shift = 0) const;

  const string& name() const
  {
    return name_;
  }
  void SetID(const vec_string &vs) const;
  void SetPtr(const double *vec) const
  {
    for (varmap_i i = v.begin(); i != v.end(); ++i)
      (*i).second.ptr = const_cast<double*>(vec) + (*i).second.i;
  }
  void SetPtr(const vec_double &vec) const
  {
    SetPtr(&*vec.begin());
  }
  double operator()(const double *vec) const
  {
    SetPtr(vec);
    return est();
  }
  double operator()(const vec_double &vec) const
  {
    SetPtr(vec);
    return est();
  }
  double operator()(const double *vec, const vec_string &vs) const
  {
    SetID(vs);
    return operator()(vec);
  }

  void SetPtr(const vec_double_ptr &vec) const
  {
    for (varmap_i i = v.begin(); i != v.end(); ++i)
      (*i).second.ptr = const_cast<double*>(vec[(*i).second.i]);
  }
  double operator()() const
  {
    return est();
  }
  double operator()(const vec_double_ptr &vec, const vec_string &vs) const
  {
    SetID(vs);
    SetPtr(vec);
    return est();
  }
  
  void SetPtr(const vec_double_2ptr &vec) const
  {
    for (varmap_i i = v.begin(); i != v.end(); ++i)
      (*i).second.ptr = const_cast<double*>(*vec[(*i).second.i]);
  }
};

OVERLOAD_STREAMS(convert)

typedef vector<convert> vec_convert;
typedef vec_convert::iterator vec_convert_i;
typedef vec_convert::const_iterator vec_convert_ci;

class RefAlgorithm
{
  friend class algorithm;
  size_t Refs;

  void AddReference()
    {++Refs;}
  size_t References()
    {return Refs;}
  size_t RemoveReference()
    {return --Refs;}

  string id;
  bool print;

protected:
  string FileName;
  ofstream *debug;
  bool debug_;

public:
  RefAlgorithm() : Refs(1), print(true), debug(0), debug_(0) {}
  RefAlgorithm(const RefAlgorithm &old) 
    : Refs(1), id(old.id), print(old.print), debug(0), debug_(0) {}
  virtual ~RefAlgorithm() {DeleteDebug();}

  virtual void clear() = 0;
  virtual RefAlgorithm* clone() const = 0;

  void DeleteDebug();
  void SetDebug(const SGML &e);

  istream& read(istream& in);
  ostream& write(ostream& out, size_t shift = 0) const;
  virtual void read(const SGML &e) = 0;
  virtual void WriteAttributes(ostream &out, size_t shift = 0) const
  {
    out << endl << PutTab(shift) << "debug=" << debug_;
    if (!FileName.empty())
      out << endl << PutTab(shift) << "DebugFile=" << FileName;
  }
  virtual void WriteBody(ostream &out, size_t shift = 0) const = 0;

  virtual int GetInputIndex(const string& s) const = 0;
  virtual void SetVariable(int i, const double &x) const = 0;
  virtual const string& name(int i) const = 0;
  virtual int GetOutputIndex(const string& s, size_t &N) const = 0;
  virtual void OutNames(int i, vec_string &vs) const = 0;
  virtual void out(int i, double *x) const = 0;

  virtual void SetSolved(bool t) const {}
};

OVERLOAD_STREAMS(RefAlgorithm)

class algorithm
{
  typedef map<string, algorithm, less<string> > map_types;
  typedef map_types::iterator map_types_i;
  typedef map_types::const_iterator map_types_ci;

  class Destruct;
  friend class algorithm::Destruct;
  class Destruct
  {
  public:
    ~Destruct();
  };
  static Destruct clean;

  static RefAlgorithm *PTR_NULL;
  static map_types *types;
  static void (*init_ptr[])();
  static map_types *obj;

  friend istream& RefAlgorithm::read(istream& in);
  friend ostream& RefAlgorithm::write(ostream& out, size_t shift = 0) const;
  static void init();

  RefAlgorithm *ptr;

  void cow()
  {
    if (ptr->References() > 1)
    {
      ptr->RemoveReference();
      ptr = ptr->clone();
    }
  }
  algorithm(RefAlgorithm* ptr) : ptr(ptr) {}

public:
  static bool debug;
  
  static void RegisterType(const string& name, RefAlgorithm* ptr);
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
  static void SetSolved(bool t)
  {
    if (obj)
      for (map_types_i i = obj->begin(); i != obj->end(); ++i)
        (*i).second.ptr->SetSolved(t);
  }

  algorithm()
  {
    init();
    ptr = PTR_NULL;
    ptr->AddReference();
  }
  algorithm(const algorithm& old)
  {
    ptr = old.ptr;
    ptr->AddReference();
  }
  algorithm(const string &name)
  {
    init();
    map_types_ci i = types->find(name);
    if (i == types->end())
      throw gError(string("algorithm: unknown type ") + name);
    else
    {
      ptr = (*i).second.ptr;
      ptr->AddReference();
    }
  }
  ~algorithm()
  {
    if (ptr->RemoveReference() == 0)
      delete ptr;
  }
  algorithm& operator=(const algorithm &old)
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
      throw gError(string("algorithm: unknown type ") + name);
    else
      *this = (*i).second;
  }
  void find(const string &name);
  void clear()
  {
    cow();
    ptr->clear();
  }
  const RefAlgorithm* pointer() const {return ptr;}

  istream& read(istream& in)
  {
    SGML e(in);
    read(e);
    return in;
  }
  void read(const SGML &e);
  ostream& write(ostream& out, size_t shift = 0) const;
  
  const string& id() {return ptr->id;}
  
  int GetInputIndex(const string& s) const {return ptr->GetInputIndex(s);}
  void SetVariable(int i, const double &x) const {ptr->SetVariable(i, x);}
  const string& name(int i) const {return ptr->name(i);}
  int GetOutputIndex(const string& s, size_t &N) const
    {return ptr->GetOutputIndex(s, N);}
  void OutNames(int i, vec_string &vs) const {return ptr->OutNames(i, vs);}
  void out(int i, double *x) const {return ptr->out(i, x);}
};

OVERLOAD_STREAMS(algorithm)

typedef vector<algorithm> vec_algorithm;
typedef vec_algorithm::iterator vec_algorithm_i;
typedef vec_algorithm::const_iterator vec_algorithm_ci;

class SimpleAlgorithm : public RefAlgorithm
{
protected:	
  mutable vec_double x;
  vec_string vs;

public:

  virtual void clear()
  {
    x.clear();
    vs.clear();
  }

  virtual void read(const SGML &e);
  virtual void WriteBody(ostream &out, size_t shift = 0) const;

  virtual int GetInputIndex(const string& s) const;
  virtual void SetVariable(int i, const double &x_) const
  {
    x[i] = x_;
  }
  virtual const string& name(int i) const;
  virtual int GetOutputIndex(const string& s, size_t &N) const
  {
    N = 1;
    return -1;
  }
  virtual void OutNames(int i, vec_string &vs_) const
  {
    vs_.push_back("f");
  }
};

class InputValue
{
  enum {Name, NumberOf, Value, Convert} source;
  int ind;
  union
  {
    mutable size_t iy;
    double x;
  };
  string iname;
  convert expr;
  
public:
  void read(const SGML &e, const algorithm &alg);
  void write(ostream &out, size_t shift, const algorithm &alg, 
      bool once = false) const;
  
  void set(const vec_string& vs) const;
  void set(const algorithm &alg, double *reg) const;
  void set(const algorithm &alg, double *reg, const vec_string &vs) const
  {
    set(vs);
    set(alg, reg);
  }
  friend bool operator<(const InputValue &x, const InputValue &y)
  {
    return x.ind < y.ind;
  }
};

class OutputValue
{
  int ind;
  size_t N;
  string name_;
  bool noname;
  
public:
  void read(const SGML &e, const algorithm &alg);
  void write(ostream &out, size_t shift, const algorithm &alg) const;
  void out(const algorithm &alg, double *x) const
  {
    alg.out(ind, x);
  }
  size_t size() const
  {
    return N;
  }
  const string& name() const
  {
    return name_;
  }
  void names(vec_string& vs, const algorithm &alg) const;
};

class compute
{
  algorithm alg;

  typedef set<InputValue, less<InputValue> > set_input;
  typedef set_input::iterator set_input_i;
  typedef set_input::const_iterator set_input_ci;

  typedef vector<OutputValue> vec_output;
  typedef vec_output::iterator vec_output_i;
  typedef vec_output::const_iterator vec_output_ci;

  mutable vec_double v1;
  mutable vec_double v2;

  vec_string nm_out;
  vec_string nm_con;
  set_input once;
  set_input in;
  vec_output out;
  vec_convert con;

  void SetIn(double *reg) const
  {
    for (set_input_ci i = in.begin(); i != in.end(); ++i)
      (*i).set(alg, reg);
  }
  void SetV1() const
  {
    double *j = &*v1.begin();
    for (vec_output_ci i = out.begin(); i != out.end(); ++i)
    {
      (*i).out(alg, j);
      j += (*i).size();
    }
  }
  void SetV2() const
  {
    for (size_t i = 0; i != con.size(); ++i)
      v2[i] = con[i]();
  }

  void SetOutConvert()
  {
    for (vec_convert_ci i = con.begin(); i != con.end(); ++i)
    {
      (*i).SetID(nm_out);
      (*i).SetPtr(v1);
    }
  }
public:
  compute()
  {
    out.reserve(4);
    con.reserve(4);
  }
  compute(const compute &x)
  {
    operator=(x);
  }
  compute& operator=(const compute &x)
  {
    if (this != &x)
    {
      alg = x.alg;
      v1 = x.v1;
      v2 = x.v2;
      nm_out = x.nm_out;
      nm_con = x.nm_con;
      once = x.once;
      in = x.in;
      out = x.out;
      con = x.con;
      SetOutConvert();
    }
    return *this;
  }

  void clear()
  {
    v1.clear();
    v2.clear();
    nm_out.clear();
    nm_con.clear();
    once.clear();
    in.clear();
    out.clear();
    con.clear();
  }
  istream& read(istream& in)
  {
    parser p(in);
    SGML e;
    read(p.GetSGML(e));
    return in;
  }
  void read(const SGML &e);
  ostream& write(ostream& out, size_t shift = 0) const;

  void SetInput(const vec_string &vs) const
  {
    for (set_input_i i = in.begin(); i != in.end(); ++i)
      (*i).set(vs);
  }
  void SetOnce(const vec_string &vs = vec_string(), double *reg = 0) const
  {
    for (set_input_i i = once.begin(); i != once.end(); ++i)
      (*i).set(alg, reg, vs);
  }
  void SetOnceInput(const vec_string &vs) const
  {
    for (set_input_i i = once.begin(); i != once.end(); ++i)
      (*i).set(vs);
  }
  void SetOnce(double *reg) const
  {
    for (set_input_i i = once.begin(); i != once.end(); ++i)
      (*i).set(alg, reg);
  }
  const double& f(double *reg = 0) const;
  const vec_double& vec(double *reg = 0) const;
  size_t size() const
  {
    if (con.empty())
      return v1.size();
    else
      return v2.size();
  }
  const vec_string& names() const
  {
    if (con.empty())
      return nm_out;
    else
      return nm_con;
  }
};

OVERLOAD_STREAMS(compute)

typedef vector<compute> vec_compute;
typedef vec_compute::iterator vec_compute_i;
typedef vec_compute::const_iterator vec_compute_ci;

#endif
