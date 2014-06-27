/*
Copyright (C) 1998, 1999, 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                       http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#include "coef.h"
#include <iterator>

UnknownCoefs *coef::u = 0;
ComputedCoefs *coef::c = 0;
coef::Destruct coef::clean;
bool coef::compiled = false;

void UnknownCoefs::Unknowns(vec_double &x, vec_double &l, vec_double &u, 
    vec_double &sc, vec_string &vs)
{
  x.clear();
  l.clear();
  u.clear();
  sc.clear();
  vs.clear();
  for (iterator i = coefs.begin(); i != coefs.end(); ++i)
    if ((*i).second.unknown)
    {
      x.push_back((*i).second.x/(*i).second.scale);
      u.push_back((*i).second.upper);
      l.push_back((*i).second.lower);
      sc.push_back((*i).second.scale);
      vs.push_back((*i).first);
    }
}

void UnknownCoefs::SetUnknowns(double *x)
{
  size_t n = 0;
  for (iterator i = coefs.begin(); i != coefs.end(); ++i)
    if ((*i).second.unknown)
      (*i).second.x = x[n++]*(*i).second.scale;
}

void coef::WriteUnknown(ostream &out, const string &id)
{
  if (!u)
    return;
  coefficient &cf = u->find(id);
  if (cf.print)
  {
    cf.print = false;
    out << "<coef id=" << id 
        << " unknown=" << cf.unknown;
    if (cf.lower != -HUGE_VAL)
      out << " lower=" << cf.lower;
    if (cf.upper != HUGE_VAL)
      out << " upper=" << cf.upper;
    if (cf.scale != 1.)
      out << " scale=" << cf.scale;
    out << "> " << *cf.ptr << " </coef>"; 
  }
  else
  {
    out << "<coef IDREF=" << id << "></coef>"; 
  }
}

void coef::WriteComputed(ostream &out, const string &id)
{
  ComputedCoefficient &cf = c->find(id);
  if (cf.print)
  {
    cf.print = false;
    out << "<coef id=" << id 
        << " computed> "
        << cf << " </coef>"; 
  }
  else
  {
    out << "<coef IDREF=" << id << " computed></coef>"; 
  }
}

void coef::WriteAll(ostream &out)
{
  SetNonPrinted();
  if (u)
  {
    for (UnknownCoefs::iterator i = u->begin(); i != u->end(); ++i)
    {
      WriteUnknown(out, (*i).first);
      out << endl;
    }
    for (ComputedCoefs::iterator i = c->begin(); i != c->end(); ++i)
    {
      WriteComputed(out, (*i).first);
      out << endl;
    }
  }
}

void coef::WriteUnknowns(ostream &out)
{
  if (u)
  {
    for (UnknownCoefs::iterator i = u->begin(); i != u->end(); ++i)
    {
      const coefficient &cf = (*i).second;
      out << "<coef id=" << (*i).first 
          << " unknown=" << cf.unknown;
      if (cf.lower != -HUGE_VAL)
        out << " lower=" << cf.lower;
      if (cf.upper != HUGE_VAL)
        out << " upper=" << cf.upper;
      if (cf.scale != 1.)
        out << " scale=" << cf.scale;
      out << "> " << *cf.ptr << " </coef>"; 
      out << endl;
    }
  }
}

void coef::SetPrinted()
{
  if (u)
  {
    for (UnknownCoefs::iterator i = u->begin(); i != u->end(); ++i)
      (*i).second.print = false;
    for (ComputedCoefs::iterator i = c->begin(); i != c->end(); ++i)
      (*i).second.print = false;
  }
}

void coef::SetNonPrinted()
{
  if (u)
  {
    for (UnknownCoefs::iterator i = u->begin(); i != u->end(); ++i)
      (*i).second.print = true;
    for (ComputedCoefs::iterator i = c->begin(); i != c->end(); ++i)
      (*i).second.print = true;
  }
}

void coef::EstimateComputed()
{
  if (c)
  {
    if (!compiled)
    {
      c->precompile(u);
      compiled = true;
    }
    c->est();
  }
}

void coef::InitialValues(istream &in, const string &str)
{
  calculator calc;
  string id;
  parser p(in, str);
  SGML el;
  while (!p.eof())
  {
    p.GetSGML(el);
    if (el.name.empty())
      break;
    if (el.defined("IDREF") 
        || el.defined("computed")
        || (id = el.FindString("id")).empty())
      continue;
    coefficient &cf = u->find(id);
    el.body = RemoveSpaces(el.body);
    if (!el.body.empty())
    {
      calc.FromString(el.body);
      cf.x = calc.est();
    }
    cf.unknown = el.FindInt("unknown", cf.unknown);
    cf.upper = el.FindDouble("upper", cf.upper);
    cf.lower = el.FindDouble("lower", cf.lower);
    cf.scale = el.FindDouble("scale", cf.scale);
  }
}

void coef::read(const SGML &e)
{
  clear();
  if (e.name != "v")
    e.compare("coef");
  calculator calc;
  coefficient cf;
  ComputedCoefficient cocf;
  if ((id = e.FindString("IDREF")).empty())
  {
    if ((id = e.FindString("id")).empty() && 
        (id = e.FindString("name")).empty())
    {
      calc.FromString(e.body);
      x_ = calc.est();
    }
    else if (e.defined("computed"))
    {
      computed = true;
      cocf.s = e.body;
      if (!c->add(id, cocf))
        cout << "coef: " << id << " is already defined - set to previous" 
          << endl;
      else
        compiled = false;
      ptr = &c->find(id).ptr;
    }
    else
    {
      calc.FromString(e.body);
      cf.x = calc.est();
      cf.unknown = e.FindInt("unknown", cf.unknown);
      cf.upper = e.FindDouble("upper", cf.upper);
      cf.lower = e.FindDouble("lower", cf.lower);
      cf.scale = e.FindDouble("scale", cf.scale);
      if (!u->add(id, cf))
      {
        cout << "coef: " << id << " is already defined - set to anonymous" 
          << endl;
        x_ = cf.x;
      }
      else
        ptr = &u->find(id).ptr;
    }
  }
  else
  {
    if (e.defined("computed"))
    {
      computed = true;
      ptr = &c->find(id).ptr;
    }
    else
      ptr = &u->find(id).ptr;
  }
}

ostream& coef::write(ostream &out) const
{
  if (!ptr)
    out << "<coef> " << x_ << " </coef>";
  else if (computed)
    WriteComputed(out, id);
  else
    WriteUnknown(out, id);
  return out;
}

ostream& CalculatorWithCoefs::write(ostream &out, size_t shift) const
{
  size_t start = 0;
  size_t end;
  for (vec_coef_ci i = v.begin(); i != v.end(); ++i)
  {
    end = text.find("xxx", start);
    if (end == string::npos)
      throw gError("CalculatorWithCoefs: corrupted structure");
    out << PutTab(shift);
    copy(text.begin() + start, text.begin() + end, 
				ostream_iterator<char>(out));
    out << (*i) << endl;
    start = end + 3;
  }
  if (text.begin() + start < text.end())
  {
    out << PutTab(shift);
    copy(text.begin() + start, text.end(), ostream_iterator<char>(out));
    out << endl;
  }
  return out;
}

string CalculatorWithCoefs::body() const
{
  ostringstream out;
  size_t start = 0;
  size_t end;
  for (vec_coef_ci i = v.begin(); i != v.end(); ++i)
  {
    end = text.find("xxx", start);
    if (end == string::npos)
      throw gError("CalculatorWithCoefs: corrupted structure");
    copy(text.begin() + start, text.begin() + end, 
      ostream_iterator<char>(out));
    out << "<coef IDREF=" << (*i).id;
    if ((*i).computed)
      out << " computed";
    out << "></coef>";
    start = end + 3;
  }
  copy(text.begin() + start, text.end(), ostream_iterator<char>(out));
  return out.str();
}

