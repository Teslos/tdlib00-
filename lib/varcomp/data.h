/*
Copyright (C) 1994 - 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                       		         http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

// sizeof(bool) is changed to sizeof(double) to make alignment for 64-bit
// SPARC processor

#ifndef _DATA_H
#define _DATA_H

#include <iostream>
#include "common.h"
#include <general.h>

char get_token(istream &in, string &str);
char skip_until(istream &in);

class data
{
  char* s;

protected:
  vec_string nm;
  size_t nreg;
  size_t Nt;
  size_t N;
  double P;
  double xav;

  size_t nx;
  double sc;

public:
	data() {init();}
	data(const data& old) {init(); copy(old);}
	data& operator=(const data &old)
	{
		if (this != &old)
			copy(old);
		return *this;
	}
  data(size_t Nt, size_t nreg)
	{
		init();
		set(Nt, nreg);
	}
	~data() {clear();}

  void init();
  void clear();
  void copy(const data &old);

	operator int() const {return s ? 1 : 0;} 
  friend ostream& operator<<(ostream &out, data &old);
  friend istream& operator>>(istream &in, data &to) {return to.read(in);}

  istream& read(istream &in, size_t np = 0, const string &xname = string(),
                const double &sc = 1.);

  void set(size_t Nt, size_t nreg);
  int Ntot() const {return Nt;}
  int Nvar() const {return nreg;}
  int Ns() const {return N;}
  double Ps() const {return P;}
  double xs() const {return xav;}
	double& scale() {return sc;}
	size_t& NOfX() {return nx;}
  double& operator()(size_t i, size_t j)
	{
		return val(i, j);
	}
  double& val(size_t i, size_t j)
	{
		return *((double*)(s + i*(sizeof(double)*nreg + sizeof(double)) 
				+ sizeof(double)*j));
	}
  double* operator()(size_t i)
	{
		return (double*)(s + i*(sizeof(double)*nreg + sizeof(double)));
	}
  bool& atr(size_t i)
	{
		return *(bool*)(s + (i + 1)*sizeof(double)*nreg + i*sizeof(double));
	}
  string& name(size_t j)
	{
		return nm[j];
	}
	const vec_string& names() const
	{
		return nm;
	}
  void set_av();
};

#endif
