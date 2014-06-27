/*
Copyright (C) 1994 - 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#include <../td_algo/phase_eq.h>
#include <iomanip>
#include "post_an.h"

vec_of vec;

void RegisterPostAnalysis()
{
  output::RegisterType("SeriesOutput", new SeriesOutput);
  output::RegisterType("ResidualOutput", new ResidualOutput);
  output::RegisterType("SpinodalOutput", new SpinodalOutput);
}

void errors1(ostream &out, bool all, bool points)
{
  for (iser = ser.begin(); iser != ser.end(); ++iser)
  {
    series& t = *iser;
    out << setw(wdth - t.ID.length()) << "x" << t.ID.c_str() << " ";
    out << setw(wdth) << t.ID.c_str() << " ";
  }
  out << endl;
  int end = 0;
  for (int i = 0; !end; ++i)
  {
    end = 1;
    ostringstream of1;
    for (iser = ser.begin(); iser != ser.end(); ++iser)
    {
      double* reg;
      series& t = *iser;
      if (i < t.Ntot())
      {
        end = 0;
        if (t.atr(i) || all)
        {
          reg = t(i);
          if (t.NOfX() < t.Nvar())
            of1 << setw(wdth) << reg[t.NOfX()] << " ";
          else
            of1 << setw(wdth) << 0. << " ";
          if (points)
            of1 << setw(wdth) << reg[t.f.NOfY()] << " ";
          else
            of1 << setw(wdth) << t.f(t(i))/sqrt(t.sr2) << " ";
        }
        else
        {
          of1 << setw(wdth) << "miss" << " ";
          of1 << setw(wdth) << "miss" << " ";
        }
      }
      else
      {
        of1 << setw(wdth) << "miss" << " ";
        of1 << setw(wdth) << "miss" << " ";
      }
    }
    if (end == 0)
    {
      out << of1.str() << endl;
    }
  }
}

void errors2(ostream &out, bool all)
{
  out << "ID err_a err_b" << endl;
  for (iser = ser.begin(); iser != ser.end(); ++iser)
  {
    series& t = *iser;
    if (!t.hide || all)
    {
      out << setw(wdth) << t.ID;
      out << setw(wdth) << t.suma/t.Ns()/sqrt(t.sr2);
      if (t.Ps())
        out << setw(wdth) << t.sumb/t.Ps()*sqrt(t.Ps()/t.Ns())/sqrt(t.sr2);
      else
        out << setw(wdth) << 0.;
      out << endl;
    }
  }
}

void SeriesOutput::read(const SGML &e)
{
  if (e.defined("ChangeXY"))
    change = true;
  else
    change = false;
  all = e.FindInt("AllPoints", 0);
  points = e.FindInt("Residuals", 0);
  points = !points;
  vs.clear();
  parser p(e);
  while (!p.eof())
  {
    vs.push_back(p.GetToken());
  }
}

void SeriesOutput::est() const
{
  nm.clear();
  if (vs.empty())
    for (iser = ser.begin(); iser != ser.end(); ++iser)
    {
      series& t = *iser;
      nm.push_back(string("x") + t.ID);
      nm.push_back(t.ID);
    }
  else
    for (size_t i = 0; i < vs.size(); ++i)
    {
      nm.push_back(string("x") + vs[i]);
      nm.push_back(vs[i]);
    }
  nc = nm.size();
  int end = 0;
  for (int i = 0; !end; ++i)
  {
    end = 1;
    vals.insert(vals.end(), nc, 0.);
    double *x = &*vals.end() - nc;
    size_t icur = 0;
    bool was_found = false;
    do
    {
      for (iser = ser.begin(); iser != ser.end(); ++iser)
      {
        double* reg;
        series& t = *iser;
        if (vs.empty() || t.ID == vs[icur])
        {
          was_found = true;
          if (i < t.Ntot())
          {
            end = 0;
            if (t.atr(i) || all)
            {
              reg = t(i);
              if (t.NOfX() < t.Nvar())
                *x++ = reg[t.NOfX()];
              else
                *x++ = 0.;
              if (points)
                *x++ = reg[t.f.NOfY()];
              else
                *x++ = t.f(t(i))/sqrt(t.sr2);
            }
            else
            {
              *x++ = HUGE_VAL;
              *x++ = HUGE_VAL;
            }
          }
          else
          {
            *x++ = HUGE_VAL;
            *x++ = HUGE_VAL;
          }
        }
      }
      if (vs.size() && !was_found)
      {
// series not fount - skipped
        *x++ = HUGE_VAL;
        *x++ = HUGE_VAL;
      }
    }
    while (vs.size() && ++icur < vs.size());
  }
  vals.erase(vals.end() - nc, vals.end());
}

void ResidualOutput::read(const SGML &e)
{
  clear();
  vs.clear();
  vx.clear();
  SGML start, finish, step;
  parser p(e);
  SGML el;
  p.GetSGML(start);
  start.compare("start");
  p.GetSGML(finish);
  finish.compare("finish");
  p.GetSGML(step);
  step.compare("step");
  res = FindResidual(p.GetToken());
  nc = 2;
  nm.push_back(e.FindString("InName", "x"));
  nm.push_back(e.FindString("OutName", "f"));
  string str;
  while (!p.eof())
  {
    p.GetSGML(el);
    el.compare("var");
    str = el.FindString("name");
    if (!str.empty())
    {
      vs.push_back(str);
      vx.push_back(el.FindDouble("value"));
    }
  }
  ReadStart(start);
  ReadFinish(finish);
  ReadStep(step);
  ix = SearchString(nm_cycle, nm[0]);
  if (ix == -1)
    throw gError(string("ResidualOutput: no name - ").append(nm[0]));
  nm_cycle.push_back(res.NameOfY());
  vds.push_back(0.);
  vd.push_back(0.);
  res.SetOnceInput(vs);
  res.SetInput(nm_cycle);
}

void ResidualOutput::WriteBody(ostream &out, size_t shift) const
{
  WriteCycle(out, shift);
  out << PutTab(shift) << res.id() << endl;
  for (size_t i = 0; i < vs.size(); ++i)
    out << "<var name=" << vs[i] << " value=" << vx[i] << "></var>" << endl;
}

void SpinodalOutput::est() const
{
  nm.clear();
  vals.clear();
  nc = 0;
  if (p.size() != 2)
    return;
  double Tstart;
  double Tfinish;
  StateTp Tp;
  StateX x(2);
  if (up)
  {
    Tstart = PhaseEquilibrium::Tmin;
    Tfinish = PhaseEquilibrium::Tmax;
    if (Tstep < 0)
      Tstep = -Tstep;
  }
  else
  {
    Tstart = PhaseEquilibrium::Tmax;
    Tfinish = PhaseEquilibrium::Tmin;
    if (Tstep > 0)
      Tstep = -Tstep;
  }
  double old, next;
  bool notgap = true;
  old = -HUGE_VAL;
  size_t NOfGaps = 0;
  nm.push_back("T");
  vals.push_back(Tstart);
  Tp.T() = Tstart;
  for (double x2 = xstep/2.; x2 < 1.; x2 += xstep)
  {
    x[1] = x2;
    x[0] = 1. - x[1];
    const vec_double &mu = p.mu(Tp, x);
    next = mu[1] - mu[0];
    if (next < old)
    {
      if (notgap)
      {
        ++NOfGaps;
        vals.push_back(x2 - xstep/2.);
        notgap = false;
      }
      if (x2 > 1. - xstep/2. - xstep/4.)
        vals.push_back(x2);
    }
    else
    {
      if (!notgap)
      {
        vals.push_back(x2 - xstep/2.);
        notgap = true;
      }
    }
    old = next;
  }
  for (size_t i = 0; i < NOfGaps; ++i)
  {
    nm.push_back(string("x2_1"));
    nm.push_back(string("x2_2"));
  }
  nc = nm.size();
  for (double T = Tstart + Tstep; 
      up ? T < Tfinish + Tstep/10. : T > Tfinish + Tstep/10.
      ; T += Tstep)
  {
    Tp.T() = T;
    vals.insert(vals.end(), nc, 0.);
    double *xcur = &*vals.end() - nc;
    double *xprev = &*vals.end() - 2*nc;
    *xcur++ = T;
    size_t n = 0;
    for (size_t i = 0; i < NOfGaps; ++i)
    {
      old = -HUGE_VAL;
      bool wasgap = false;
      for (double x2 = xprev[1 + i*2] - xstep/2.; 
          x2 <= xprev[2 + i*2] + xstep/2. + xstep/4.; 
          x2 += xstep)
      {
        x[1] = x2;
        x[0] = 1. - x[1];
        const vec_double &mu = p.mu(Tp, x);
        next = mu[1] - mu[0];
        if (next < old)
        {
          *xcur++ = x2 - xstep/2.;
          wasgap = true;
          break;
        }
        old = next;
      }
      if (wasgap)
      {
        ++n;
        for (double x2 = xprev[2 + i*2] + xstep/2.; 
            x2 >= xprev[1 + i*2] - xstep/2. - xstep/4.; 
            x2 -= xstep)
        {
          x[1] = x2;
          x[0] = 1. - x[1];
          const vec_double &mu = p.mu(Tp, x);
          next = mu[1] - mu[0];
          if (next > old)
          {
            *xcur++ = x2 + xstep/2.;
            break;
          }
          old = next;
        }
      }
      else
      {
        *xcur++ = HUGE_VAL;
        *xcur++ = HUGE_VAL;
      }
    }
    if (!n)
      break;
  }
}

void post_analysis()
{
  if (write_model_file)
  {
    WriteModelFile(vec);
  }
  if (!file_name.empty())
  {
    ofstream out((file_name + ".par").c_str());
    out.precision(16);
    coef::WriteUnknowns(out);
  }
  if (write_e2 && !file_name.empty() && !ser.empty())
  {
    ofstream out((file_name + ".ae2").c_str());
    out.precision(6);
    out.setf(ios::showpoint);
  // print eb vs ea
    errors2(out, false);
    out.close();
  }
  for (vec_of_ci i = vec.begin(); i != vec.end(); ++i)
    (*i).out(file_name);
}

