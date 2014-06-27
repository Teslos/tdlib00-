/*
Copyright (C) 1998, 1999, 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                       http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#ifndef __ELEM_H
#define __ELEM_H

#include <general.h>
#include <set>

class elem
{
public:
  struct elem_property
  {
    char name[2];
    double mass;

    elem_property()
      {name[0] = name[1] = 0; mass = 0.;}
    istream& read(istream& in)
    {
      if (!in)
        return in;
      in >> name[0];
      if (!isspace(in.peek()))
        in.get(name[1]);
      if (in && !isspace(in.peek()))
        throw gError("Too long name for the element");
      in >> mass;
      return in;
    }
    ostream& write(ostream& out) const
    {
      if (name[0])
        out << name[0];
      else
        return out;
      if (name[1])
        out << name[1];
      return out << '\t' << mass;
    }
  };
  typedef vector<elem_property> vec_el;
  typedef vec_el::iterator vec_el_i;
  typedef vec_el::const_iterator vec_el_ci;

private:
  class Destruct;
  friend class elem::Destruct;
  class Destruct
  {
  public:
    ~Destruct()
    {
      if (elem::v)
      {
        delete elem::v;
        elem::v = 0;
        if (elem::s)
        {
          delete elem::s;
          elem::s = 0;
        }
      }
    }
  };
  static Destruct clean;

  short int el;

  static vec_el *v;
  static set<elem, greater<elem> > *s;
  static void init();

public:
  static void ReadElements(istream& in);
  static void ReadElements(const SGML &e);
  static ostream& WriteElements(ostream& out);
  static const set<elem, greater<elem> >& elements()
  {
    if (!s)
      s = new set<elem, greater<elem> >;
      return *s;
  }

  elem()
    {init(); el = -1;}
  elem(const char *name, bool reg = true)
    {init(); istringstream in(name); read(in, reg);}
  elem(const elem& old)
    {el = old.el;}
  elem(const pair<const elem, double> &y)
    {el = y.first.el;}
  elem& operator=(elem old)
    {el = old.el; return *this;}
  elem& operator=(const pair<const elem, double> &y)
    {return *this = y.first;}

  istream& read(istream& in, bool reg = true);
  ostream& write(ostream& out) const
  {
    if (el == -1) return out;
    if ((*v)[el].name[0])
      out << (*v)[el].name[0];
    else
      return out;
    if ((*v)[el].name[1])
      out << (*v)[el].name[1];
    return out;
  }

  operator int() const
    {return el + 1;}
  void clear()
    {el = -1;}
  double mass() const
    {return el == -1 ? 0 : (*v)[el].mass;}
};

OVERLOAD_STREAMS(elem)

inline void destroy(elem*) {}
inline void destroy(pair<const elem, double>*) {}

inline bool operator<(const elem &x, const pair<const elem, double> &y)
{
  return x < y.first;
}

typedef set<elem, greater<elem> > set_elem;
typedef set_elem::iterator set_elem_i;
typedef set_elem::const_iterator set_elem_ci;

#endif
