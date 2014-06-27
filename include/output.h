/*
Copyright (C) 1998, 1999, 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                       http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#ifndef _OUTPUT_H
#define _OUTPUT_H

#include "td_algo.h"

class RefOutput
{
  friend class output;
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
  mutable vec_string nm;
  mutable vec_double vals;
  mutable size_t nc;
  
public:
  RefOutput() : Refs(1), print(true) {vals.reserve(200);}
  RefOutput(const RefOutput &old) 
    : Refs(1), id(old.id), print(old.print), nm(old.nm), vals(old.vals), 
      nc(old.nc) {vals.reserve(200);}
  virtual ~RefOutput() {}

  virtual void clear()
  {
    nc = 0;
    nm.clear();
    vals.clear();
  }
  virtual RefOutput* clone() const = 0;

  istream& read(istream& in);
  ostream& write(ostream& out, size_t shift = 0) const;
  virtual void read(const SGML &e) = 0;
  virtual void WriteAttributes(ostream &out, size_t shift = 0) const {}
  virtual void WriteBody(ostream &out, size_t shift = 0) const = 0;
  
  virtual bool IsLine() const {return true;}
  virtual bool ChangeXY() const {return false;}
  virtual void est() const = 0;
  
  size_t NRows() const {return vals.size()/nc;}
  size_t NCols() const {return nc;}
  const string& name(size_t i) const {return nm[i];}
  const double& value(size_t i, size_t j) const {return vals[i*nc + j];}
};

OVERLOAD_STREAMS(RefOutput)

class output
{
  typedef map<string, output, less<string> > map_types;
  typedef map_types::iterator map_types_i;
  typedef map_types::const_iterator map_types_ci;

  class Destruct;
  friend class output::Destruct;
  class Destruct
  {
  public:
    ~Destruct();
  };
  static Destruct clean;

  static RefOutput *PTR_NULL;
  static map_types *types;
  static void (*init_ptr[])();
  static map_types *obj;

  friend istream& RefOutput::read(istream& in);
  friend ostream& RefOutput::write(ostream& out, size_t shift = 0) const;
  static void init();

  RefOutput *ptr;

  void cow()
  {
    if (ptr->References() > 1)
    {
      ptr->RemoveReference();
      ptr = ptr->clone();
    }
  }
  output(RefOutput* ptr) : ptr(ptr) {}

public:
  static void RegisterType(const string& name, RefOutput* ptr);
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

  output()
  {
    init();
    ptr = PTR_NULL;
    ptr->AddReference();
  }
  output(const output& old)
  {
    ptr = old.ptr;
    ptr->AddReference();
  }
  output(const string &name)
  {
    init();
    map_types_ci i = types->find(name);
    if (i == types->end())
      throw gError(string("output: unknown type ") + name);
    else
    {
      ptr = (*i).second.ptr;
      ptr->AddReference();
    }
  }
  ~output()
  {
    if (ptr->RemoveReference() == 0)
      delete ptr;
  }
  output& operator=(const output &old)
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
      throw gError(string("output: unknown type ") + name);
    else
      *this = (*i).second;
  }
  void find(const string &name);
  void clear()
  {
    cow();
    ptr->clear();
  }
  const RefOutput* pointer() const {return ptr;}

  istream& read(istream& in)
  {
    SGML e(in);
    read(e);
    return in;
  }
  void read(const SGML &e);
  ostream& write(ostream& out, size_t shift = 0) const;
  
  const string& id() {return ptr->id;}
  
  bool IsLine() const {return ptr->IsLine();}
  bool ChangeXY() const {return ptr->ChangeXY();}
  void est() const {ptr->est();}
  
  size_t NRows() const {return ptr->NRows();}
  size_t NCols() const {return ptr->NCols();}
  const string& name(size_t i) const {return ptr->nm[i];}
  const double& value(size_t i, size_t j) const {return ptr->value(i, j);}
};

OVERLOAD_STREAMS(output)

typedef vector<output> vec_output;
typedef vec_output::iterator vec_output_i;
typedef vec_output::const_iterator vec_output_ci;

class NullOutput : public RefOutput
{
public:
  virtual NullOutput* clone() const {return new NullOutput(*this);}

  virtual void read(const SGML &e) {}
  virtual void WriteBody(ostream &out, size_t shift = 0) const {}
  
  virtual void est() const {}
};

class comparison
{
  enum Operation {LE, EQ, GE};
  static const char *OperationName[];
  typedef Enum<Operation, OperationName> operation;
  
  operation op;
  
public:
  comparison() {}
  comparison(const string &s) {operator=(s);}
  comparison& operator=(const string &s)
  {
    op = s;
    return *this;
  }
  friend ostream& operator<<(ostream &out, const comparison& x)
  {
    return out << x.op;
  }
  
  bool operator()(const double& x, const double &y) const
  {
    switch (op)
    {
      case LE: 
        return x <= y;
      case EQ:
//        if (isnan(x) && isnan(y))
        if (x == HUGE_VAL && y == HUGE_VAL)
          return true;
        else
          return x == y;
      case GE: 
        return x >= y;
    }
  }
};

typedef vector<comparison> vec_comparison;
typedef vec_comparison::iterator vec_comparison_i;
typedef vec_comparison::const_iterator vec_comparison_ci;

class CycleOutput : public RefOutput
{
protected:  
// start
  vec_string nm_cycle;
  mutable vec_double vds;
  compute start;
// finish
  vec_int vif;
  mutable vec_double vdf;
  vec_comparison vof;
  compute finish;
  
  size_t iout;
  size_t jout;
  comparison comp;
// step
  vec_int vic;
  vec_convert vc;

  mutable int abort;
  mutable vec_double vd;

public:
  void clear()
  {
    RefOutput::clear();
    nm_cycle.clear();
    vds.clear();
    start.clear();
    vif.clear();
    vdf.clear();
    vof.clear();
    finish.clear();
    vic.clear();
    vc.clear();
    vd.clear();
  }
  
  void ReadStart(const SGML &e);
  void ReadFinish(const SGML &e);
  void ReadStep(const SGML &e);
  void WriteCycle(ostream &out, size_t shift) const;
  
  virtual void ComputeCurrentPoint(double *reg, double*) const = 0;
  virtual void SetOnce() const = 0;
  
  virtual void est() const;
};

class ComputeOutput : public CycleOutput
{
  vec_compute vec;
  bool NoCycle;
  bool change;
public:
  ComputeOutput() {vec.reserve(6);}
  virtual void clear()
  {
    CycleOutput::clear();
    vec.clear();
  }
  virtual ComputeOutput* clone() const {return new ComputeOutput(*this);}

  virtual void read(const SGML &e);
  virtual void WriteAttributes(ostream &out, size_t shift = 0) const
  {
    if (change)
      out << endl << PutTab(shift) << "ChangeXY";
  }
  virtual void WriteBody(ostream &out, size_t shift = 0) const;
  
  virtual bool ChangeXY() const {return change;}
  virtual void est() const;
  virtual void ComputeCurrentPoint(double *reg, double *x) const;
  virtual void SetOnce() const
  {
    for (vec_compute_ci i = vec.begin(); i != vec.end(); ++i)
      (*i).SetOnce(&*vds.begin());
  }
};

class OutputFile
{
  string ext;
  vec_output vo;
  string format;

public:
  static int width;

  OutputFile() : format("file") 
  {
    vo.reserve(6); 
    ext = "axm";
  }

  void clear()
  {
    ext = "axm";
    vo.clear();
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
  
  void out(const string &basename) const;
};

OVERLOAD_STREAMS(OutputFile)

typedef vector<OutputFile> vec_of;
typedef vec_of::iterator vec_of_i;
typedef vec_of::const_iterator vec_of_ci;

#endif
