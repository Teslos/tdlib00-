/*
Copyright (C) 1994 - 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#ifndef _POST_AN_H
#define _POST_AN_H

#include <iterator>
#include "sumsqr.h"
#include <output.h>

extern vec_of vec;

void ReadGlobals(const SGML &e);
ostream& WriteGlobals(ostream &out, size_t shift = 0);

void ReadModelFile(istream &in, const string &str, vec_of &vec);
void WriteModelFile(const vec_of &vec);

void errors1(ostream &out, bool all = false, bool points = false);
void errors2(ostream &out, bool all = false);

class SeriesOutput : public RefOutput
{
  bool all;
  bool points;
  vec_string vs;
  bool change;
  
public:
  virtual SeriesOutput* clone() const {return new SeriesOutput(*this);}

  virtual void read(const SGML &e);
  virtual void WriteAttributes(ostream &out, size_t shift = 0) const 
  {
    out << endl << PutTab(shift) << "AllPoints=" << all
      << endl << PutTab(shift) << "Residuals=" << !points;
    if (change)
      out << endl << PutTab(shift) << "ChangeXY";
  }
  virtual void WriteBody(ostream &out, size_t shift = 0) const 
  {
    if (vs.size())
    {
      out << PutTab(shift);
      copy(vs.begin(), vs.end(), ostream_iterator<string>(out, " "));
      out << endl;
    }
  } 
  virtual bool IsLine() const {return false;}
  virtual bool ChangeXY() const {return change;}
  virtual void est() const;
};

class ResidualOutput : public CycleOutput
{
  residual res;
  vec_string vs;
  vec_double vx;
  int ix;
public:
  virtual ResidualOutput* clone() const {return new ResidualOutput(*this);}

  virtual void read(const SGML &e);
  virtual void WriteAttributes(ostream &out, size_t shift = 0) const
  {
    if (nm.size())
      out << endl << PutTab(shift) << "InName=" << nm[0] <<
        endl << PutTab(shift) << "OutName=" << nm[1];
  }
  virtual void WriteBody(ostream &out, size_t shift = 0) const;
  virtual void SetOnce() const
  {
    res.SetOnce(const_cast<double*>(&*vx.begin()));
  }
  virtual void ComputeCurrentPoint(double *reg, double *x) const
  {
    *x++ = reg[ix];
    *x = res.f(reg);
  }
};

class SpinodalOutput : public RefOutput
{
  phase p;
  mutable double Tstep;
  double xstep;
  bool up;
  
public:
  virtual SpinodalOutput* clone() const {return new SpinodalOutput(*this);}

  virtual void read(const SGML &e)
  {
    p.find(e.FindString("phase"));
    Tstep = e.FindDouble("StepT", 50);
    xstep = e.FindDouble("StepX", 0.01);
    up = true;
    if (e.FindString("direction") == "down")
      up = false;
  }
  virtual void WriteAttributes(ostream &out, size_t shift = 0) const 
  {
    out << endl << PutTab(shift) << "phase=" << p.id()
      << endl << PutTab(shift) << "StepT=" << Tstep
      << endl << PutTab(shift) << "StepX=" << xstep;
    out << endl << PutTab(shift) << "direction=";
    if (up)
      out << "up";
    else
      out << "down";
  }
  virtual void WriteBody(ostream &out, size_t shift = 0) const {} 
  virtual void est() const;
};

#endif
