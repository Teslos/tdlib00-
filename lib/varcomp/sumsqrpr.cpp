/*
Copyright (C) 1994 - 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#include <iomanip>
#include <math.h>
#include "sumsqr.h"

void print_ss_par_1(double *x, size_t n)
{
  ostream &of = *out_file;
  double ssq = 0.;
  iser = ser.begin();
  while (iser != ser.end())
  {
    if (!(*iser).hide)
      ssq += (*iser).SS;
    iser++;
  }
  of << "SSQ " << setw(wdth) << ssq << endl;
  of << "scaled unknows" << endl;
  for (size_t j = 0; j < n; j++)
    of << setw(wdth) << xnm[j].c_str() << setw(wdth) << x[j] << endl;
  of << endl;
}

void print_variances(ostream &of, bool print_always_gb)
{
  iser = ser.begin();
  while(iser != ser.end())
  {
    series& t = *iser;
    if (t.hide && (print_series == all_series))
    {
      of << '*';
      set_sums_1(t);
      if (t.Ns() > 1)
      {
        if (t.Ps())
        {
          if (t.Ns() == 2)
          {
            t.sr2 = t.sumi/2.;
            t.gb = 0.;
            t.ga = 0.;
          }
          else
          {
            t.sr2 = (t.sumi - t.suma*t.suma/t.Ns() - t.sumb*t.sumb/t.Ps())
                               /(t.Ns() - 2);
            t.gb = (t.sumb*t.sumb/t.Ps()/t.sr2 - 1)/t.Ps();
            t.ga = (t.suma*t.suma/t.Ns()/t.sr2 - 1)/t.Ns();
          }
        }
        else
        {
          t.sr2 = (t.sumi - t.suma*t.suma/t.Ns())/(t.Ns() - 1);
          t.gb = 0.;
          t.ga = (t.suma*t.suma/t.Ns()/t.sr2 - 1)/t.Ns();
        }
      }
      else
      {
        t.sr2 = t.sumi;
        t.ga = 0.;
        t.gb = 0.;
      }
      if (t.ga < 0.) t.ga = 0.;
      if (t.gb < 0.) t.gb = 0.;
    }
    else if (!t.hide)
      of << ' ';

    if (!t.hide || (print_series == all_series))
    {
      of << setw(6) << t.ID.c_str() << setw(wdth+2) << t.f.id().c_str() << ", ";
      switch(t.fl_sr2)
      {
        case series::fixed:
          of << "*  ";
          break;
        case series::same:
          of << "#";
          if (t.isr2 == 0)
            of << "  ";
          else if (t.isr2 < 10)
            of << t.isr2 << " ";
          else
            of << t.isr2;
          break;
        case series::own:
          of << "%  ";
          break;
      }
      of << setw(wdth-2) << sqrt(t.sr2);
      if (t.adjust)
        of << setw(wdth-2) << sqrt(t.adjust) << ", ";
      else
        of << ", ";
      switch(t.fl_ga)
      {
        case series::fixed:
          of << "*  ";
          break;
        case series::same:
          of << "#";
          if (t.igat == 0)
            of << "  ";
          else if (t.igat < 10)
            of << t.igat << " ";
          else
            of << t.igat;
          break;
        case series::own:
          of << "%  ";
          break;
      }
      of << setw(wdth-2) << sqrt(t.ga) << ", ";
      if (t.Ps() || print_always_gb)
      {
        switch(t.fl_gb)
        {
          case series::fixed:
            of << "*  ";
            break;
          case series::same:
            of << "#";
            if (t.igbt == 0)
              of << "  ";
            else if (t.igbt < 10)
              of << t.igbt << " ";
            else
              of << t.igbt;
            break;
          case series::own:
            of << "%  ";
            break;
        }
        of << setw(wdth-2) << sqrt(t.gb);
      }
      else
        of << "       n/a  ";
      of << ';' << endl;
    }
    iser++;
  }
  of << endl;
}

void print_all_points()
{
  ostream &of = *out_file;
  double tmp;
  of << " ID           eq            sr             sga            sgb" << endl;
  print_variances(of);
  of << endl << "Results " << endl;
  of << "     ID         av_dev        sri      err_a        err_b      err1_b" << endl;
  iser = ser.begin();
  while(iser != ser.end())
  {
    series& t = *iser;
    if (t.hide && (print_series == all_series))
      of << '*';
    else if (!t.hide)
      of << ' ';

    if (!t.hide || (print_series == all_series))
    {
      of << setw(wdth) << t.ID.c_str();
      of << setw(wdth) << sqrt(t.sumi/t.Ns());
      if (t.Ns() == 1)
        of << setw(wdth) << 0.;
      else if (t.Ps())
        {
          if (t.Ns() == 2)
            of << setw(wdth) << 0.;
          else
          {
          tmp = (t.sumi - t.suma*t.suma/t.Ns() - t.sumb*t.sumb/t.Ps())/t.Ns();
          if (tmp > 0.)
            of << setw(wdth) << sqrt(tmp);
          else
            of << setw(wdth) << 0.;
          }
        }
      else
      {
        tmp = (t.sumi - t.suma*t.suma/t.Ns())/t.Ns();
        if (tmp > 0.)
          of << setw(wdth) << sqrt(tmp);
        else
          of << setw(wdth) << 0.;
      }
      of << setw(wdth) << t.suma/t.Ns();
      if (t.Ps())
        of << setw(wdth) << t.sumb/t.Ps()
           << setw(wdth) << t.sumb/t.Ps()*sqrt(t.Ps()/t.Ns()) << endl;
      else
        of << "       n/a         n/a " << endl;
    }
    iser++;
  }
  of << endl;

  of << "     ID       Ntot     Ns      NOfX      xav       Ps       ScaleOfX" 
    << endl;
  iser = ser.begin();
  while(iser != ser.end())
  {
    series& t = *iser;
    if (t.hide && (print_series == all_series))
      of << '*';
    else if (!t.hide)
      of << ' ';

    if (!t.hide || (print_series == all_series))
    {
      of << setw(wdth) << t.ID.c_str();
      of << setw(5) << t.Ntot() << setw(5) << t.Ns();
      of << setw(wdth) << t.NOfX() + 1
         << setw(wdth) << t.xs()
         << setw(wdth) << sqrt(t.Ps())
         << setw(wdth) << t.scale() << endl;
    }
    iser++;
  }
  of << endl;
  if (print_series == no_series) return;

  iser = ser.begin();
  while(iser != ser.end())
  {
    int i, j;
    double f;
    series& t = *iser;
    if (!t.hide || (print_series == all_series))
    {
      of << endl;
      for (i = 0; i < t.Nvar()*wdth + 2; i++)
        of << ' ';
      of << " err_full    err_a";
      if (t.Ps())
        of << "   err_ab" << endl;
      else
        of << endl;
    }
    if (t.hide && (print_series == all_series))
      of << '*';
    else if (!t.hide)
      of << ' ';

    if (!t.hide || (print_series == all_series))
    {
      of << setw(wdth) << t.ID.c_str() << ',' << endl << t.f.id().c_str() 
        << ',' << endl;
      of << ' ';
      for (j = 0; j < t.Nvar(); j++)
        of << setw(wdth) << t.name(j).c_str();
      of << ',' << endl;
      size_t k = t.NOfX();
      for (i = 0; i < t.Ntot(); i++)
      {
        if (t.atr(i))
          of << ' ';
        else
          of << '*';
        for (j = 0; j < t.Nvar(); j++)
          of << setw(wdth) << t.val(i, j);
        f = t.f(t(i));
        of << setw(wdth) << f;
        if (t.Ns())
        {
          of << setw(wdth) << f - t.suma/t.Ns();
          if (t.Ps())
            of << setw(wdth)
            << f - t.suma/t.Ns() - t.sumb/t.Ps()
               *(t.val(i, k) - t.xs())/t.scale();
        }
        if (i < t.Ntot() - 1)
          of << ',' << endl;
        else
          of << ';' << endl;
      }
    }
    iser++;
  }
}
