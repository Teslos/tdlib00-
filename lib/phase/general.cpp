/*
Copyright (C) 1998, 1999, 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                       http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#include <math.h>
#include <float.h>
#include <algorithm>
#include "general.h"

const char *FunctionName[] = {"G", "H", "S", "Cp", "V", "dVdT", "dVdp", 0};
const char *IndexName[] = {"", "excess", "ideal", "mix", "ref", "", "ref_ideal", 
  "full", 0};
double global::R = 8.31441;
double global::T = 1000.;
double global::p = 1;
double global::step1 = 10*sqrt(DBL_EPSILON);
double global::step2 = 10*pow(DBL_EPSILON, 1./3.);
bool global::neg_log = true;
double global::eps = 1e-7;

size_t SGML::dummy;
string SGML::dummystring;

list<parser::descriptor> *parser::debug = 0;

istream& CheckInf::read(istream &in) const
{
  in >> ws;
  int sign = 1;
  if (in.peek() == '-')
  {
    sign = -1;
    in.ignore();
  }
  else if (in.peek() == '+')
    in.ignore();
  char t1, t2;
  in.get(t1);
  in.get(t2);
  if (toupper(t1) == 'I' &&
      toupper(t2) == 'N' &&
      toupper(in.peek()) == 'F')
  {
    in.ignore();
    if (sign == 1)
      *xr = HUGE_VAL;
    else
      *xr = -HUGE_VAL;
  }
  else
  {
    in.putback(t2);
    in.putback(t1);
    in >> *xr;
    if (!in.fail())
      *xr *= sign;
  }
  return in;
}

ostream& CheckInf::write(ostream &out) const
{
  if (*xr == HUGE_VAL) out << "inf";
  else if (*xr == -HUGE_VAL) out << "-inf";
  else out << *xr;
  return out;
}

istream& SGML::read(istream &in, size_t &lf)
{
  name.erase();
  body.erase();
  attr.clear();
  parser p(in);
  if (p.SkipWS() != '<')
  {
    lf += p.lf();
    return in;
  }
  in.ignore();
  while ((name = p.GetToken()) == "!--")
// treating comments
  {
    while (in && (p.GetToken() != "--" || in.peek() != '>'));
    in.ignore();
    if (in >> ws, in.peek() != '<')
    {
      name.erase();
      return in;
    }
    in.ignore();
  }
  string key, value;
  p.GetToken();
  while (in && (p.token != ">"))
  {
    if (p.token == "=")
      throw gError("SGML: no_key");
    key = p.token;
    p.GetToken();
    if (p.token == "=")
    {
      value = p.GetToken();
      if (value == "=" || value == ">")
        throw gError("SGML: no_value");
      p.GetToken();
    }
    else
      value = "";
    attr[key] = value;
  }
  string finish = '/' + name;
  size_t count = 1;
  string spc;
  for (;;)
  {
    while (in && (p.GetToken() != "<"))
    {
      body.append(p.space).append(p.token);
    }
    spc = p.space;
    if (p.GetToken() == finish)
    {
      if (!--count)
      {
        p.SkipUntil('>');
        p.ReadCompare(">");
        break;
      }
    }
    else if (!in)
    {
      throw gError(string("SGML: no end ").append(finish));
    }
    else if (p.token == name)
      ++count;
    body.append(spc).append("<");
    body.append(p.space).append(p.token);
  }
  lf += p.lf();
  return in;
}

ostream& SGML::write(ostream &out, size_t shift) const
{
  out << PutTab(shift) << "<" << name;
  for (map_string_string_ci i = attr.begin(); i != attr.end(); ++i)
  {
    out << endl << PutTab(shift + 1) << (*i).first;
    if (!(*i).second.empty())
      out << "=" << (*i).second;
  }
  out << ">" << endl;
  if (!body.empty())
    out << PutTab(shift + 1) << body << endl;
  out << PutTab(shift) << "</" << name << ">" << endl;
  return out;
}

void WriteDebug(string &s)
{
  if (!parser::debug)
    return;
  ostringstream out;
  out << endl;
  for (list<parser::descriptor>::iterator i = parser::debug->begin(); 
      i != parser::debug->end(); ++i)
  {
    out << "from: " << (*i).from << "(" << (*i).from_id << ")"
      << " last: " << (*i).last << "(" << (*i).last_id << ")"
      << " lf: " << (*i).lf << endl;
  }
  s += out.str();
}

void parser::init(const string &from, const string &from_id)
{
  if (!debug)
  {
    debug = new list<descriptor>;
  }
  idebug = debug->insert(debug->end(), descriptor());
  (*idebug).from = from;
  (*idebug).from_id = from_id;
  (*idebug).lf = 0;
}

parser::~parser()
{
  if (del)
    delete in;
  debug->erase(idebug);
  if (debug->empty())
  {
    delete debug;
    debug = 0;
  }
}

string& parser::GetToken()
{
  char ch = 0;
  space.erase();
  while (in->get(ch), *in && isspace(ch)) 
  {
    if (ch == '\n')
      (*idebug).lf++;
    space.append(1, ch);
  }
  token.assign(1, ch);
  if (ch == '<' 
      || ch == '>' 
      || ch == '=' 
      || ch == '{' 
      || ch == '}'
      || ch == '('
      || ch == ','
      || ch == ')'
      )
    return token;
  while (in->get(ch), *in 
         && !isspace(ch) 
         && ch != '<' 
         && ch != '>' 
         && ch != '=' 
         && ch != '{' 
         && ch != '}'
         && ch != '('
         && ch != ','
         && ch != ')'
         )
    token.append(1, ch);
  in->putback(ch);
  return token;
}

void parser::compare(const string &old)
{
  string t = old;
  size_t pos;
  t.erase(0, t.find_first_not_of(" "));
  while (1)
  {
    pos = t.find(' ');
    if (t.compare(0, pos, token) != 0)
      throw gError(string("parser: expected ").append(t, 0, pos));
    t.erase(0, pos);
    t.erase(0, t.find_first_not_of(" "));
    if (t.empty())
      break;
    GetToken();
  }
}

int parser::SkipWS()
{
  char ch = 0;
  while (in->get(ch), *in && isspace(ch)) 
  {
    if (ch == '\n')
      (*idebug).lf++;
  }
  in->putback(ch);
  return in->peek();
}

bool parser::SkipUntil(char t)
{
  while (in->peek() != EOF && in->peek() != t)
  {
    if (in->peek() == '\n')
      (*idebug).lf++;
    in->ignore();
  }
  if (in->peek() == t)
    return true;
  else
    return false;
}

SGML& parser::GetSGML(SGML &e)
{
  e.read(*in, (*idebug).lf);
  (*idebug).last = e.name;
  (*idebug).last_id = e.FindString("id");
  return e;
}

