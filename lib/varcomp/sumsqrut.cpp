/*
Copyright (C) 1994 - 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#include <iomanip>
#include <math.h>
#include "sumsqr.h"

static int itmp;

void series::init()
{
  vs.clear();
  vx.clear();
  hide = false;
  fl_sr2 = own;
  fl_ga = fixed;
  fl_gb = fixed;
  isr2 = 0;
  igat = 0;
  igbt = 0;
  adjust = 0.;
  sr2 = 1.;
  sr2old = 1.;
  ga = 0.;
  gaold = 0.;
  gb = 0.;
  gbold = 0.;
}

void anal_series()
{
  int i;
  iser = ser.begin();
  while(iser != ser.end())
  {
    series& t = *iser;
    if (!t.Ns())
      t.hide = true;
    else if (!t.hide)
    {
      if (!t.Ps())
      {
        t.fl_gb = series::fixed;
        t.gb = 0.;
      }
      if (t.fl_sr2 == series::same)
        Nst[t.isr2]++;
      if (t.fl_ga == series::same)
        Ngat[t.igat]++;
      if (t.fl_gb == series::same)
        Ngbt[t.igbt]++;
    }
    iser++;
  }
  for (i = 0; i < nvar; i++)
  {
    if (Nst[i] == 1)
    {
      Nst[i] = 0;
      warning("number of the same sr2[i] is one - set to own");
      iser = ser.begin();
      while(iser != ser.end())
      {
        series& t = *iser;
        if (!t.hide && t.isr2 == i && t.fl_sr2 == series::same)
        {
          t.fl_sr2 = series::own;
          break;
        }
        iser++;
      }
    }
    if (Ngat[i] == 1)
    {
      Ngat[i] = 0;
      warning("number of the same ga[i] is one - set to own");
      iser = ser.begin();
      while(iser != ser.end())
      {
        series& t = *iser;
        if (!t.hide && t.igat == i && t.fl_ga == series::same)
        {
          t.fl_ga = series::own;
          break;
        }
        iser++;
      }
    }
    if (Ngbt[i] == 1)
    {
      Ngbt[i] = 0;
      warning("number of the same gb[i] is one - set to own");
      iser = ser.begin();
      while(iser != ser.end())
      {
        series& t = *iser;
        if (!t.hide && t.igbt == i && t.fl_gb == series::same)
        {
          t.fl_gb = series::own;
          break;
        }
        iser++;
      }
    }
  }
}

void set_sums_1(series& t)
{
  double f;
  int k;
  t.sumi = 0.;
  t.suma = 0.;
  t.sumb = 0.;
  size_t nx = t.NOfX();
  t.f.SetOnce(&*t.vx.begin());
  for (k = 0; k < t.Ntot(); k++)
  {
    if (t.atr(k))
    {
      f = t.f(t(k));
      t.sumi += f*f;
      t.suma += f;
      if (t.Ps())
        t.sumb += f*(t(k, nx) - t.xs());
    }
  }
  t.sumb /= t.f.ScaleOfX();
  if (t.adjust)
  {
    if (t.Ps())
      t.sumi += t.adjust*(t.Ns()-2);
    else if (t.Ns() > 0)
      t.sumi += t.adjust*(t.Ns()-1);
  }
}

void set_sums(series& t, double* f, int& i)
{
  int is = i;
  double da, db;
  int k;
  t.sumi = 0.;
  t.suma = 0.;
  t.sumb = 0.;
  size_t nx = t.NOfX();
  t.f.SetOnce(&*t.vx.begin());
  for (k = 0; k < t.Ntot(); k++)
  {
    if (t.atr(k))
    {
      f[i] = t.f(t(k));
      t.sumi += f[i]*f[i];
      t.suma += f[i];
      if (t.Ps())
        t.sumb += f[i]*(t(k, nx) - t.xs());
      i++;
    }
  }
  t.sumb /= t.f.ScaleOfX();
  da = (1. - 1./sqrt(1. + t.Ns()*t.ga))/t.Ns();
  if (t.Ps())
    db = (1. - 1./sqrt(1. + t.Ps()*t.gb))/t.Ps();
  else
    db = 0.;
  t.SS = 0.;
  for (k = 0; k < t.Ntot(); k++)
  {
    if (t.atr(k))
    {
      f[is] -= da*t.suma;
      if (t.Ps())
        f[is] -= db*(t(k, nx) - t.xs())/t.f.ScaleOfX()*t.sumb;
      f[is] /= sqrt(t.sr2);
      t.SS += f[is]*f[is];
      is++;
    }
  }
  if (t.adjust)
  {
    if (t.Ps())
    {
      t.sumi += t.adjust*(t.Ns()-2);
      t.SS += t.adjust*(t.Ns()-2);
    }
    else if (t.Ns() > 0)
    {
      t.sumi += t.adjust*(t.Ns()-1);
      t.SS += t.adjust*(t.Ns()-1);
    }
  }
}

static double SSold = HUGE_VAL;

void ss(double* x, int m, int n, double* f)
{
  int i = 0;
  algorithm::SetSolved(false);
  coef::SetUnknowns(x);
  coef::EstimateComputed();
  iser = ser.begin();
  while(iser != ser.end())
  {
    series& t = *iser;
    if (!t.hide)
      set_sums(t, f, i);
    iser++;
  }
  iss++;
  if (print == all_iter)
  {
    *out_file << "  SS  function eval. " << setw(3) << iss << ", ";
    print_ss_par_1(x, n);
  }
  double SS = 0.;
  for (size_t i = 0; i < m; ++i)
    SS += f[i]*f[i];
  if (!file_name.empty() && SS < SSold)
  {
    SSold = SS;
    ofstream out((file_name + ".par").c_str());
    out.precision(16);
    coef::WriteUnknowns(out);
  }
  cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << setw(4) << iss 
    << setw(12) << SS;
  cout.flush();
}

void eval_L(double& L, double& SS)
{
  L = 0.;
  SS = 0.;
  iser = ser.begin();
  while(iser != ser.end())
  {
    series& t = *iser;
    if (!t.hide)
    {
      L -= t.Ns()*log(t.sr2) + log(1. + t.Ns()*t.ga) + log(1. + t.Ps()*t.gb);
      SS += (t.sumi - t.ga/(1. + t.Ns()*t.ga)*t.suma*t.suma
                    - t.gb/(1. + t.Ps()*t.gb)*t.sumb*t.sumb)/t.sr2;
    }
    iser++;
  }
  L -= SS;
}

double gat_solve(double x)
{
  if (x < 0.) return 0.;
  double f = 0.;
  int i = itmp;
  iser = ser.begin();
  while (iser != ser.end())
  {
    series& t = *iser;
    if (t.fl_ga == series::same && !t.hide && t.igat == i)
    {
      if (t.fl_sr2 == series::same)
        f += t.suma*t.suma/sr2t[t.isr2]/pow((1. + t.Ns()*x), 2.)
                      - t.Ns()/(1. + t.Ns()*x);
      else
        f += t.suma*t.suma/t.sr2/pow((1. + t.Ns()*x), 2.)
                      - t.Ns()/(1. + t.Ns()*x);
    }
    iser++;
  }
  if (print == all_iter)
    *out_file << "gat_solve " << i << setw(wdth) << x
       << setw(wdth) << f << endl;
  return f;
}

double gbt_solve(double x)
{
  if (x < 0.) return 0.;
  double f = 0.;
  int i = itmp;
  iser = ser.begin();
  while (iser != ser.end())
  {
    series& t = *iser;
    if (t.fl_gb == series::same && !t.hide && t.igbt == i)
    {
      if (t.fl_sr2 == series::same)
        f += t.sumb*t.sumb/sr2t[t.isr2]/pow((1. + t.Ps()*x), 2.)
                      - t.Ps()/(1. + t.Ps()*x);
      else
        f += t.sumb*t.sumb/t.sr2/pow((1. + t.Ps()*x), 2.)
                      - t.Ps()/(1. + t.Ps()*x);
    }
    iser++;
  }
  if (print == all_iter)
    *out_file << "gbt_solve " << i << setw(wdth) << x
       << setw(wdth) << f << endl;
  return f;
}

void eval_var_comp()
{
  ostream &of = *out_file;
  double L = 1.;
  double Lold = 0.;
  double SS;
  double epsL = 1e-7;
  double a, b;
  int i;

  int iter = 0;

  while (fabs((L - Lold)/L) > epsL && iter < 50)
  {
    iter++;
    Lold = L;
    for (i = 0; i < nvar; i++)
    {
      sr2t[i] = 0.;
      Nst[i] = 0;
    }
    iser = ser.begin();
    while (iser != ser.end())
    {
      series& t = *iser;
      if (!t.hide)
      {
        switch (t.fl_sr2)
        {
          case series::fixed:
            break;
          case series::same:
            sr2t[t.isr2] += t.sumi - t.ga/(1. + t.Ns()*t.ga)*t.suma*t.suma
                            - t.gb/(1. + t.Ps()*t.gb)*t.sumb*t.sumb;
            Nst[t.isr2] += t.Ns();
            break;
          case series::own:
            t.sr2 = (t.sumi - t.ga/(1. + t.Ns()*t.ga)*t.suma*t.suma
                            - t.gb/(1. + t.Ps()*t.gb)*t.sumb*t.sumb)/t.Ns();
            if (t.sr2 <=0)
              throw gError(string("eval_var_comp: sr2 is zero for ") + t.ID);
            break;
          default:
            warning("check fl._sr2");
        }
      }
      iser++;
    }
    for (i = 0; i < nvar; i++)
      if (Nst[i])
      {
        if (sr2t[i] <= 0.)
        {
          throw gError(string("eval_var_comp: sr2t[")
                       + ObjToString(i) + "] is zero");
        }
        sr2t[i] /= Nst[i];
      }

    for (i = 0; i < nvar; i++)
    {
      itmp = i;
      if (Ngat[i] > 1)
      {
        a = 0.;
        b = 1000.;
        try
        {
          gat[i] = rroot(makeFunctor((D2D*)0, gat_solve), a, b);
        }
        catch (err_rroot)
        {
          gat[i] = 0.;
        }
      }

      if (Ngbt[i] > 1)
      {
        a = 0.;
        b = 10000.;
        try
        {
          gbt[i] = rroot(makeFunctor((D2D*)0, gbt_solve), a, b);
        }
        catch (err_rroot)
        {
          gbt[i] = 0.;
        }
      }
    }

    iser = ser.begin();
    while (iser != ser.end())
    {
      series& t = *iser;
      if (!t.hide)
      {
        if (t.fl_sr2 == series::same)
          t.sr2 = sr2t[t.isr2];
        if (t.Ps())
        {
          switch (t.fl_gb)
          {
            case series::fixed:
              break;
            case series::same:
              t.gb = gbt[t.igbt];
              break;
            case series::own:
              t.gb = (t.sumb*t.sumb/t.Ps()/t.sr2 - 1.)/t.Ps();
              break;
            default:
              warning("check fl_gb");
          }
        }
        switch (t.fl_ga)
        {
          case series::fixed:
            break;
          case series::same:
            t.ga = gat[t.igat];
            break;
          case series::own:
            t.ga = (t.suma*t.suma/t.Ns()/t.sr2 - 1.)/t.Ns();
            break;
          default:
            warning("check fl_ga");
        }
        if (t.ga < 0.) t.ga = 0;
        if (t.gb < 0.) t.gb = 0;
      }
      iser++;
    }
    eval_L(L, SS);
    if (print == all_iter)
    {
      of << "var_comp - iteration " << setw(4) << iter << endl;
      of << "value of L is " << setw(15) << L
         << "  and of SS is " << setw(15)
         << SS << endl;
    }
    cout << "\b\b\b\b" << setw(4) << iter;
  }
  if (print != results)
  {
    of << "New variance components after " << setw(4) << iter
       << "  iterations" << endl;
    if ( print != all_iter)
      of << "value of L is " << setw(15) << L
         << "  and of SS is " << setw(15)
         << SS << endl;
    print_all_points();
  }
  else
  {
    of << "new variance components after " << iter << " iteration" << endl;
    of << "value of L is " << setw(15) << L
       << "  and of SS is " << setw(15)
       << SS << endl;
  }
}

