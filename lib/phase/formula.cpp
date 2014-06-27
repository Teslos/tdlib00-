/*
Copyright (C) 1998, 1999, 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                       http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#include "formula.h"

formula::set_mf *formula::smf = 0;
const MolecularFormula *formula::PTR_NULL = 0;
formula::Destruct formula::clean;

ostream& MolecularFormula::write(ostream& out) const
{
  for (MolecularFormula_ci i = begin(); i != end(); ++i)
  {
    elem e("e", false);
    if ((*i).first == e)
    {
      if ((*i).second == 1)
      {
        out << '-';
        continue;
      }
      else if ((*i).second == -1)
      {
        out << '+';
        continue;
      }
      else
      {
        out << '(' << (*i).first << (*i).second << ')';
        continue;
      }
    }
    out << (*i).first;
    if ((*i).second != 1.)
      out << (*i).second;
  }
  return out;
}

double MolecularFormula::noa(elem el) const
{
  MolecularFormula_ci i = find(el);
  if (i == end())
    return 0;
  else
    return (*i).second;
}

double MolecularFormula::mass() const
{
  double mass = 0.;
  for (MolecularFormula_ci i = begin(); i != end(); ++i)
    mass += (*i).second*(*i).first.mass();
  return mass;
}

const vec_double& MolecularFormula::FormulaVector() const
{
  const set_elem &els = elem::elements();
  if (fv.size() != els.size())
  {
    fv.resize(els.size());
    size_t i = 0;
    for (set_elem_ci is = els.begin(); is != els.end(); ++i, ++is)
      fv[i] = noa(*is);
  }
  return fv;
}

istream& formula::read(istream& in)
{
  clear();
  if (!in) return in;
  elem t;
  elem e("e", false);
  MolecularFormula mf;
  double coef;
  ostringstream tmp;
  while (true)
  {
    t.read(in);
    if (!t)
    {
      if (in.peek() == '(')
      {
        in.ignore();
        formula t;
        t.read(in);
        if (in.peek() == ')')
        {
          in.ignore();
          if (isspace(in.peek()))
            coef = 1.;
          else
          {
            char sign = in.peek();
            in >> coef;
            if (in.fail())
            {
              in.clear();
              coef = 1.;
              if (sign == '+' || sign == '-')
                in.putback(sign);
            }
          }
          if (coef)
          {
            tmp << '(' << t.second << ')';
            if (coef != 1.)
              tmp << coef;
            for (MolecularFormula_ci i = t.first->begin(); i != t.first->end(); ++i)
              mf[(*i).first] += (*i).second*coef;
          }
        }
        else
        {
          clear();
          throw gError("formula: no right )");
        }
        continue;
      }
      if (in.peek() == '+')
      {
        tmp << '+';
        in.ignore();
        if (e)
          --mf[e];
      }
      if (in.peek() == '-')
      {
        tmp << '-';
        in.ignore();
        if (e)
          ++mf[e];
      }
      if (in.peek() == '_')
      {
        parser p(in);
        tmp << p.GetToken();
      }
      second = tmp.str();
      first = &*(smf->insert(smf->begin(), mf));
      return in;
    }
    else
    {
      if (isspace(in.peek()))
        coef = 1.;
      else
      {
        char sign = in.peek();
        in >> coef;
        if (in.fail())
        {
          in.clear();
          coef = 1.;
          if (sign == '+' || sign == '-')
            in.putback(sign);
        }
      }
      if (coef)
      {
        tmp << t;
        if (coef != 1.)
          tmp << coef;
        mf[t] += coef;
      }
    }
  }
}

