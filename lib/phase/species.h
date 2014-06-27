/*
Copyright (C) 1998, 1999, 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                       http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#ifndef __SPECIES_H
#define __SPECIES_H

#include <func_tp.h>
#include <formula.h>

class Ref_species;

class species
{
private:
  typedef map<string, species, less<string> > map_species;
  typedef map_species::iterator map_species_i;
  typedef map_species::const_iterator map_species_ci;

  class Destruct;
  friend class species::Destruct;
  class Destruct
  {
  public:
    ~Destruct();
  };
  static Destruct clean;

  Ref_species *ptr;

  static Ref_species *PTR_NULL;
  static map_species *s;

  static void init();                         //inline
  void cow();                                 //inline

public:
  static void SetPrinted();                       //inline
  static void SetNonPrinted();                     //inline
  static void WriteAll(ostream& out);  //inline

  species();                                    // inline
  species(const species& old);                  // inline
  species(const string &name);                  // inline
  ~species();                                   //inline
  species& operator=(const species &old);       // inline

  void create(const string &name);            //inline
  void find(const string &id) {create (id);}

  void clear();                               // inline

  istream& read(istream& in)
  {
    SGML e(in);
    read(e);
    return in;
  }
  void read(const SGML &e);
  ostream& write(ostream& out, size_t shift = 0) const;

  limits limitT() const;                     // inline
  limits limitp() const;                     // inline
  bool IsInside(const StateTp &Tp) const;   // inline

  const string& id() const;                   // inline
  const MolecularFormula& mf() const;               // inline
  const formula& GetFormula() const;          // inline
  const string& cf() const;                   // inline
  const vector<species>& ref_plane() const;       // inline
  const vec_double& coefs() const;            // inline

  double Z(function f, index i, const StateTp &Tp) const;
  double Z(function f, const StateTp &Tp) const
    {return Z(f, ::full, Tp);}
  double G(const StateTp &Tp) const
    {return Z(::G, ::full, Tp);}
};

OVERLOAD_STREAMS(species)

typedef vector<species> vec_species;
typedef vec_species::iterator vec_species_i;
typedef vec_species::const_iterator vec_species_ci;

class Ref_species
{
  friend class species;

  size_t Refs;
  void AddReference()
    {++Refs;}
  size_t References()
    {return Refs;}
  size_t RemoveReference()
    {return --Refs;}

  func_Tp f;
  formula fml;
  vec_species r;
  vec_double c;
  string id;
  bool p;

  Ref_species() : Refs(1), p(true) {}
  Ref_species(const Ref_species &old) :
    Refs(1), f(old.f), fml(old.fml), r(old.r), c(old.c), id(old.id), p(old.p) {}
  void clear()
  {
    p = true;
    r.clear();
    c.clear();
    id.erase();
  }
};

inline species::Destruct::~Destruct()
{
  if (species::s)
  {
    delete species::s;
    delete species::PTR_NULL;
    species::s = 0;
    species::PTR_NULL = 0;
  }
}

inline void species::init()
{
  if (!s)
  {
    s = new map_species;
    PTR_NULL = new Ref_species;
  }
}

inline void species::cow()
{
  if (ptr->References() > 1)
  {
    ptr->RemoveReference();
    ptr = new Ref_species(*ptr);
  }
}

inline void species::SetNonPrinted()
{
  if (s)
    for (map_species_i i = s->begin(); i != s->end(); ++i)
      (*i).second.ptr->p = true;
}

inline void species::SetPrinted()
{
  if (s)
    for (map_species_i i = s->begin(); i != s->end(); ++i)
      (*i).second.ptr->p = false;
}

inline void species::WriteAll(ostream& out)
{
  SetNonPrinted();
  if (s)
    for (map_species_i i = s->begin(); i != s->end(); ++i)
      (*i).second.write(out);
}

inline species::species()
{
  init();
  ptr = PTR_NULL;
  ptr->AddReference();
}

inline species::species(const species& old)
{
  init();
  ptr = old.ptr;
  ptr->AddReference();
}

inline species::species(const string &id)
{
  init();
  map_species_ci i = s->find(id);
  if (i == s->end())
    throw gError(string("species: unknown id ") + id);
  else
  {
    ptr = (*i).second.ptr;
    ptr->AddReference();
  }
}

inline species::~species()
{
  if (ptr->RemoveReference() == 0)
    delete ptr;
}

inline species& species::operator=(const species &old)
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

inline void species::create(const string &name)
{
  map_species_ci i = s->find(name);
  if (i == s->end())
    throw gError(string("spesies: unknown type ") + name);
  else
    *this = (*i).second;
}

inline void species::clear()
{
  if (ptr->References() > 1)
  {
    ptr->RemoveReference();
    ptr = new Ref_species;
  }
  else
    ptr->clear();
}

inline limits species::limitT() const
{
  return ptr->f.limitT();
}

inline limits species::limitp() const
{
  return ptr->f.limitp();
}

inline bool species::IsInside(const StateTp &Tp) const
{
  return ptr->f.IsInside(Tp);
}


inline const string& species::id() const
{
  return ptr->id;
}

inline const MolecularFormula& species::mf() const
{
  return ptr->fml.mf();
}

inline const string& species::cf() const
{
  return ptr->fml.cf();
}

inline const formula& species::GetFormula() const
{
  return ptr->fml;
}

inline const vector<species>& species::ref_plane() const
{
  return ptr->r;
}

inline const vec_double& species::coefs() const
{
  return ptr->c;
}

class GetNuError : public gError
{
public:
  GetNuError(const string &str) : gError(str) {}
};

class matrix : public vec_double
{
  size_t m;
  size_t n;

public:
  matrix() : m(0), n(0) {}
  matrix(size_t m_, size_t n_) : vec_double(m_*n_), m(m_), n(n_) {}

  void resize(size_t m_, size_t n_)
    {m = m_; n = n_; vec_double::resize(m_*n_);}
  void reserve(size_t m_, size_t n_)
    {vec_double::reserve(m_*n_);}
  void clear()
    {m = n = 0; vec_double::clear();}

  size_t NRows() const {return m;}
  size_t NCols() const {return n;}
  
  double& operator()(size_t i, size_t j)
  {
    return operator[](i + j*m);
  }
  const double& operator()(size_t i, size_t j) const
  {
    return operator[](i + j*m);
  }
};

class GetNu
{
  const set_elem &els;
  vec_species bs;
  vec_formula fbs;
  bool from_species;
  mutable vec_double fm;
  mutable vec_int ipiv;
  mutable vec_int indx;
  mutable size_t Nel;
  mutable size_t Nc;
  mutable int info;
  mutable vec_double b;
  void SetFormulaMatrix() const;
public:
  static double tol; //DBL_EPSILON*1000.
  
  GetNu() : els(elem::elements()), Nel(elem::elements().size()), Nc(0) {}
  GetNu(const vec_species &r) : els(elem::elements())
    {reset(r);}
  GetNu(const vec_formula &r) : els(elem::elements())
    {reset(r);}
  const vec_species& basis() const {return bs;}
  void reset(const vec_species &r);
  const vec_int& reset(const vec_formula &r);
  size_t NOfComps() const {return Nc;}
  void nu(const species &s, double *n) const  {nu(s.cf(), n);}
  void nu(const formula &frm, double *n) const;
  void FormulaMatrix(const vec_formula &vf, matrix &A);
};

#endif

