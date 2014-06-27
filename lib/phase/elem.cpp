/*
Copyright (C) 1998, 1999, 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                       http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#include "elem.h"

elem::Destruct elem::clean;
elem::vec_el *elem::v = 0;
set_elem *elem::s = 0;

static const char* LIST_ELEMENTS = "\
e  0.0005485   \
O  15.9994     \
H  1.0079      \
D  2.014       \
T  3.016       \
F  18.9984     \
Cl 35.453      \
Br 79.904      \
I  126.904     \
At 210         \
He 4.0026      \
Ne 20.179      \
Ar 39.948      \
Kr 83.8        \
Xe 131.3       \
Rn 222         \
S  32.06       \
Se 78.9        \
Te 127.6       \
Po 209         \
N  14.0067     \
P  30.9738     \
As 74.9216     \
Sb 121.75      \
Bi 208.98      \
C  12.011      \
Si 28.0855     \
Ge 72.59       \
Sn 118.69      \
Pb 207.2       \
B  10.81       \
Al 26.9815     \
Ga 69.72       \
In 114.82      \
Tl 204.37      \
Zn 65.38       \
Cd 112.41      \
Hg 200.59      \
Cu 63.546      \
Ag 107.868     \
Au 196.967     \
Fe 55.847      \
Co 58.9332     \
Ni 58.7        \
Ru 101.07      \
Rh 102.906     \
Pd 106.4       \
Os 190.2       \
Ir 192.22      \
Pt 195.09      \
Mn 54.938      \
Tc 97          \
Re 186.207     \
Cr 51.996      \
Mo 95.94       \
W  183.85      \
V  50.9414     \
Nb 92.9064     \
Ta 180.948     \
Ti 47.9        \
Zr 91.22       \
Hf 178.49      \
Sc 44.9559     \
Y  88.9059     \
La 138.906     \
Ce 140.12      \
Pr 140.908     \
Nd 144.24      \
Pm 145         \
Sm 150.4       \
Eu 151.96      \
Gd 157.25      \
Tb 158.925     \
Dy 162.5       \
Ho 164.93      \
Er 167.26      \
Tm 168.934     \
Yb 173.04      \
Lu 174.97      \
Ac 227         \
Th 232.038     \
Pa 231.036     \
U  238.029     \
Np 237.048     \
Pu 244         \
Am 243         \
Cm 247         \
Bk 247         \
Cf 251         \
Es 254         \
Fm 257         \
Md 258         \
No 259         \
Be 9.01218     \
Mg 24.305      \
Ca 40.08       \
Sr 87.62       \
Ba 137.33      \
Ra 226.025     \
Li 6.941       \
Na 22.9898     \
K  39.0983     \
Rb 85.4678     \
Cs 132.905     \
Fr 223 ";

void elem::init()
{
  if (!v)
  {
    v = new vec_el;
    v->reserve(110);
    istringstream in2(LIST_ELEMENTS);
    while (in2)
    {
      elem_property t;
      t.read(in2);
      if (t.name[0])
        v->push_back(t);
    }
    if (!s)
      s = new set_elem;
  }
}

void elem::ReadElements(istream &in)
{
  SGML e(in);
  ReadElements(e);
}

void elem::ReadElements(const SGML &e)
{
  e.compare("elements");
  if (v)
    throw gError("elem: the list of elements is already initialized");
  else
  {
    v = new vec_el;
    if (!s)
      s = new set_elem;
  }
  v->reserve(110);
  istringstream in2(e.body);
  while (in2)
  {
    elem_property t;
    t.read(in2);
    if (t.name[0])
      v->push_back(t);
  }
}

ostream& elem::WriteElements(ostream& out)
{
  if (v)
  {
    out << "<elements>" << endl;
    for (vec_el_ci i = v->begin(); i != v->end(); ++i)
    {
      out << '\t';
      i->write(out);
      out << endl;
    }
    out << "</elements>" << endl;
  }
  return out;
}

istream& elem::read(istream& in, bool reg)
{
  el = -1;
  if (!in) return in;
  char name = 0;
  in.get(name);
  for (vec_el_ci i = v->begin(); i != v->end(); ++i)
  {
    if (name == i->name[0])
    {
      if (i->name[1])
      {
        if (i->name[1] == in.peek())
        {
          el = i - v->begin();
          in.ignore();
          break;
        }
      }
      else
        el = i - v->begin();
    }
  }
  if (el == -1)
    in.putback(name);
  else
  {
    if (reg)
      s->insert(*this);
  }
  return in;
}

