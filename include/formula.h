/*
Copyright (C) 1998, 1999, 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                       http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#ifndef __FORMULA_H
#define __FORMULA_H

#include <elem.h>

typedef map <elem, double, greater<elem> > map_elem_double;
typedef map_elem_double::iterator MolecularFormula_i;
typedef map_elem_double::const_iterator MolecularFormula_ci;

class MolecularFormula : public map_elem_double
{
mutable vec_double fv;
public:
  MolecularFormula() {}
  
  ostream& write(ostream& out) const;
  size_t nel() const
    {return size();}
  double noa(elem el) const;
  double mass() const;
  const vec_double& FormulaVector() const;
};

class formula
{
  typedef set<MolecularFormula> set_mf;
  typedef set_mf::iterator set_mf_i;
  typedef set_mf::const_iterator set_mf_ci;
  static set_mf *smf;
  static const MolecularFormula *PTR_NULL;
  static void init()
  {
    if (!smf)
    {
      smf = new set_mf;
      PTR_NULL = &*(smf->insert(smf->begin(), MolecularFormula()));
    }
  }
  class Destruct;
  friend class formula::Destruct;
  class Destruct
  {
  public:
    ~Destruct()
    {
      if (smf)
      {
        delete formula::smf;
        formula::smf = 0;
        formula::PTR_NULL = 0;
      }
    }
  };
  static Destruct clean;

  const MolecularFormula *first;
  string second;

public:
  formula()
  {
    init();
    first = PTR_NULL;
  }
  formula(const char *in)
  {
    init(); 
    first = PTR_NULL;
    istringstream inp(in); 
    read(inp);
  }
  ~formula()
    {}
  istream& read(istream& in);
  ostream& write(ostream& out) const
    {return out << second;}
  ostream& write_mf(ostream& out) const
    {return first->write(out);}
  bool empty() const
    {return first->empty();}
  void clear()
    {first = PTR_NULL; second.erase();}

  const MolecularFormula& mf() const
    {return *first;}
  const string& cf() const
    {return second;}

  operator const string&() const
    {return second;}
    
  size_t nel() const
    {return first->size();}
  double noa(elem el) const
    {return first->noa(el);}
  double mass() const
    {return first->mass();}
  const vec_double& FormulaVector() const
    {return first->FormulaVector();}
};

inline bool operator==(const formula &x, const formula &y)
{
  return x.cf() == y.cf();
}

inline bool operator!=(const formula &x, const formula &y)
{
  return x.cf() != y.cf();
}

inline bool operator<(const formula &x, const formula &y)
{
  return x.mf() < y.mf() || 
         (!(y.mf() < x.mf()) && x.cf() < y.cf());
}

OVERLOAD_STREAMS(formula)

typedef vector<formula> vec_formula;
typedef vec_formula::iterator vec_formula_i;
typedef vec_formula::const_iterator vec_formula_ci;

#endif
