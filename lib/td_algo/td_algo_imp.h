/*
Copyright (C) 1998, 1999, 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                       http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#ifndef __TD_ALGO_IMP_H
#define __TD_ALGO_IMP_H

#include <td_algo.h>

class PassThrough : public RefAlgorithm
{
  mutable vec_double x;
  mutable vec_string vs;

public:

  virtual void clear()
  {
    x.clear();
    vs.clear();
  }
  virtual PassThrough* clone() const {return new PassThrough(*this);}

  virtual void read(const SGML &e) {}
  virtual void WriteBody(ostream &out, size_t shift = 0) const {}

  virtual int GetInputIndex(const string& s) const;
  virtual void SetVariable(int i, const double &x_) const
  {
    x[i] = x_;
  }
  virtual const string& name(int i) const
  {
    return vs[i];
  }
  virtual int GetOutputIndex(const string& s, size_t &N) const
  {
    N = 1;
    return GetInputIndex(s);
  }
  virtual void OutNames(int i, vec_string &vs_) const
  {
    vs_.push_back(vs[i]);
  }
  virtual void out(int i, double *x_) const
  {
    *x_ = x[i];
  }
};

class PhaseProperty : public RefAlgorithm
{
  phase ph;
  StateTp Tp;
  StateX x;
  mutable vec_double st;
  mutable vec_double vint;
  vec_string nm_st;
  vec_string nm_st2;
  size_t dmf;

  mutable bool eval;
  mutable bool eval_vint;

  enum source {Property, Partial, Internal, Stability, AllState};
  struct descriptor
  {
    string name;
    function f;
    index i;
    int ix;
    size_t N;
    source from;
  };
  mutable vector<descriptor> vd;

  int SearchState(const string& str) const
  {
    for (size_t j = 0; j < nm_st.size(); ++j)
      if (nm_st[j] == str || nm_st2[j] == str)
        return j;
    return -1;
  }
public:
  virtual void clear();
  virtual PhaseProperty* clone() const {return new PhaseProperty(*this);}

  virtual void read(const SGML &e);
  virtual void WriteAttributes(ostream &out, size_t shift = 0) const;
  virtual void WriteBody(ostream &out, size_t shift = 0) const;

  void assign(const phase &ph_)
  {
    ph = ph_;
    clear();
  }
  virtual int GetInputIndex(const string& s) const;
  virtual void SetVariable(int i, const double &x_) const
  {
    st[-(i + 1)] = x_;
    if (i < -2)
      eval = false;
    eval_vint = false;
  }
  virtual const string& name(int i) const;
  virtual int GetOutputIndex(const string& s, size_t &N) const;
  virtual void OutNames(int i, vec_string &vs) const;
  virtual void out(int i, double *x) const;
  virtual void SetSolved(bool t) const
  {
    eval_vint = t;
  }
};

class reaction : public RefAlgorithm
{
  vec_compute reac;
  mutable vec_string vs;
  mutable vec_double reg;

  mutable bool set;
  size_t OutN;

  void SetOnce(const vec_string &vs = vec_string(), double *reg = 0) const
  {
    for (vec_compute_ci i = reac.begin(); i != reac.end(); ++i)
      (*i).SetOnce(vs, reg);
  }
  void SetInput(const vec_string &vs) const
  {
    for (vec_compute_ci i = reac.begin(); i != reac.end(); ++i)
      (*i).SetInput(vs);
  }

public:
  reaction() : OutN(0) {reac.reserve(6);}
  
  void clear()
  {
    reac.clear();
    vs.clear();
    reg.clear();
    set = false;
    OutN = 0;
  }
  
  virtual reaction* clone() const {return new reaction(*this);}

  virtual void read(const SGML &e);
  virtual void WriteBody(ostream &out, size_t shift = 0) const;

  virtual int GetInputIndex(const string& s) const;
  
  virtual void SetVariable(int i, const double &x) const 
  {
    reg[i] = x;
  }
  virtual const string& name(int i) const;
  virtual int GetOutputIndex(const string& s, size_t &N) const
  {
    N = OutN;
    return -1;
  }
  virtual void OutNames(int i, vec_string &vs_) const
  {
    for (vec_compute_ci i = reac.begin(); i != reac.end(); ++i)
      vs_.insert(vs_.end(), (*i).names().begin(), (*i).names().end());
  }
  virtual void out(int i, double *x) const;
};

#endif
