/*
Copyright (C) 1994 - 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                       		         http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#include <iomanip>
#include <math.h>
#include "data.h"
#include <td_algo.h>

const char WARN2[] =
      "you try to exceed the number of colomns allowed in this set";
const char WARN3[] =
      "the number of colomns is not enough in this set - it is ignored";

char get_token(istream &in, string &str)
{
  char ch;
  str.erase();
  in >> ws;
  while (ch = in.get(), ch != ',' &&
                        ch != ';' &&
                        !isspace(ch) &&
                        ch != EOF)
    str += ch;
	return ch;
}

char skip_until(istream &in)
{
	char ch;
  while (ch = in.get(), ch != ',' &&
                        ch != ';' &&
                        ch != EOF);
  return ch;
}


void data::init()
{
  s = 0;
  nreg = 0;
  Nt = 0;
  N = 0;
  P = 0.;
  xav = 0.;
	nx = 0;
	sc = 1.;
}

void data::clear()
{
	if (s)
	  delete [] s;
	init();
	nm.clear();
}

void data::copy(const data &old)
{
  clear();
  if (!old) 
		return;
  nreg = old.nreg;
  Nt = old.Nt;
  N = old.N;
  P = old.P;
  xav = old.xav;
  nm = old.nm;
	sc = old.sc;
	nx = old.nx;
  int len = Nt*(sizeof(double)*nreg + sizeof(double));
  s = new char[len];
  if (s)
    memcpy(s, old.s, len);
  else
    clear();
}

ostream& operator<<(ostream &out, data &old)
{
  if (!old) 
		return out;
  size_t i, j;
  out << ' ';
  for (j = 0; j < old.nreg; j++)
    out << setw(15) << old.nm[j].c_str();
  out << ',' << endl;
  for (i = 0; i < old.Nt; i++)
  {
    if (old.atr(i))
      out << ' ';
    else
      out << '*';
    for (j = 0; j < old.nreg; j++)
      out << setw(15) << old.val(i, j);
    if (i < old.Nt - 1)
      out << ',' << endl;
    else
      out << ';' << endl;
  }
  return out;
}

istream& data::read(istream &in, size_t np, const string &xname, 
		const double &sc_)
{
  char ch;
	string str;
  clear();
	in >> ws;
	if (ch = in.peek(), ch == ','
			|| ch == ';'
			|| ch == EOF)
		return in;
  do
  {
		ch = get_token(in, str);
		if (str.empty())
			break;
		nm.push_back(str);
    nreg++;
    if(np && nreg == np)
    {
      warning(WARN2);
      cout << "allowed number is " << np
           << ", current number is " << nreg
           << ", left characters are ";
      while (ch = in.get(),
             ch != ',' &&
             ch != ';' &&
             ch != EOF)
        cout << ch;
      cout << endl;
      break;
    }
  }
	while (isspace(ch));
	size_t nx_ = SearchString(nm, xname);
  if (!nreg)
  {
    if (ch == ',')
      while (ch = in.get(), ch != ';' &&
                            ch != EOF);
    return in;
  }
  if (np && nreg < np)
  {
      warning(WARN3);
      if (ch == ',')
        while (ch = in.get(), ch != ';' &&
                              ch != EOF);
      clear();
      return in;
  }
  if (ch != ',')
  {
    clear();
    return in;
  }
  bool incl;
  double work;
  size_t j;
	double imit_bool;
  do
  {
    incl = true;
    in >> ws;
    if (in.peek() == '*')
    {
      incl = false;
      in.ignore();
    }
    for (j = 0; j < nreg; j++)
    {
      work = 0;
      in >> work;
      in.clear();
      Buffer.write(&work, sizeof(work));
    }
		memcpy(&imit_bool, &incl, sizeof(bool));
    Buffer.write(&imit_bool, sizeof(imit_bool));
		ch = skip_until(in);
    Nt++;
  }
  while (ch == ',');
  int len = Nt*(sizeof(double)*nreg + sizeof(double));
  s = new char[len];
  if (s && len == Buffer.busy())
  {
    Buffer.read(s);
    nx = nx_;
    sc = sc_;
    set_av();
  }
  else
  {
    clear();
    Buffer.clear();
  }
  return in;
}

void data::set(size_t Nt, size_t nreg)
{
  clear();
  s = new char[Nt*(sizeof(double)*nreg + sizeof(double))];
  if (s)
  {
		nm.resize(nreg);
    data::Nt = Nt;
    data::nreg = nreg;
		nx = 1;
		sc = 1.;
    size_t i, j;
    for (i = 0; i < Nt; i++)
    {
      atr(i) = true;
      for (j = 0; j < nreg; j++)
        val(i, j) = 0.;
    }
  }
}

void data::set_av()
{
  size_t i;
  N = 0;
  for (i = 0; i < Nt; i++)
    if (atr(i)) N++;
  xav = 0.;
  P = 0.;
  if (nreg <= 1 || nx >= nreg)
    return;
	if (nreg > 1)
  {
    for (i = 0; i < Nt; i++)
		{
      if (atr(i)) 
				xav += val(i, nx);
		}
		if (N)
    {
      xav /= N;
      for (i = 0; i < Nt; i++)
        if (atr(i)) 
					P += pow((val(i, nx) - xav), 2);
      P /= sc*sc;
    }
  }
}

