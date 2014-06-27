/*
Copyright (C) 1994 - 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

//#include <new.h>
#include <iomanip>
#include <math.h>
#include <toms.h>
#include <lapack.h>
#include "opt.h"

/*
void (*old_handler)();

void my_handler()
{
  warning("memory is out - you are going to have problems");
  set_new_handler(old_handler);
}
*/

extern const char* HELP;

int main(int argc, char *argv[])
{
  try
  {
//    old_handler = set_new_handler(my_handler);
    optimizer opt;

    for (int ii = 0; ii < nvar; ++ii)
    {
      sr2t[ii] = 1.;
      Nst[ii] = 0;
      gat[ii] = 0.;
      Ngat[ii] = 0;
      gbt[ii] = 0.;
      Ngbt[ii] = 0;
    }
    
    init_before_reading();
    int arg = argc - 1;
    process_ini("assess", arg, false);  
    if (!open_streams(argc, argv, true, arg))
      return 1;
    if (!arg)
    {
      cout << HELP;
      return 0;
    }
    if (!create_output())
      return 1;
    init_after_reading();
    ostream &of = *out_file;
    if (ser.empty())
    {
      post_analysis();
      return 0;
    }
    anal_series();

    size_t n; //number of unknowns
    int m; //number of points
    size_t i, j;

  // parameters for main loop
    int iter;
    double eps_L = 1e-5;
    double eps_v = 1e-5;
    double difL;
    double difv;

  // other stuff
    double SS;
    double L = 0.;
    double Lold = 0.;

    difL = eps_L + 1.;
    difv = eps_v + 1.;
    vec_double x;
    vec_double u;
    vec_double l;
    vec_double sc;
    coef::Unknowns(x, l, u, sc, xnm);
    n = x.size();
    if (n)
    {
      of << endl << endl << "-------------Start of calculation-------------"
         << endl;
      of << "initial values of parameters and bounds (lower - upper)" << endl;
      for (j = 0; j < n; j++)
      {
        of << setw(wdth) << xnm[j].c_str() 
          << setw(wdth) << x[j]*sc[j]
          << setw(wdth) << l[j] 
          << setw(wdth) << u[j] 
          << endl;
      }
      of << endl;
    }

  // calculations
    iter = 1;
    if (fixmodel)
      n = 0;
    m = 0;
    iser = ser.begin();
    while (iser != ser.end())
    {
      if (!(*iser).hide)
        m += (*iser).Ns();
      iser++;
    }
    if (m && n)
    {
      vec_double fv(m);
      opt.init(n, m);
      vec_double xjtj(n*(n+1)/2);
      if (start == variances)
      {
        iser = ser.begin();
        while(iser != ser.end())
        {
          series& t = *iser;
          if (!t.hide)
          set_sums_1(t);
          iser++;
        }
        cout << endl << "   var_comp       ";
        eval_var_comp();
      }

      while (iter <= max_iter && (difL > eps_L || difv > eps_v))
      {
        iss = 0;
        cout << endl;
        cout << "iteration " << setw(4) << iter << "    SS                    ";
        if (print != results)
          of << endl << "***** iteration" << setw(4) << iter
             << "   *******" << endl;
        opt.solve(&*x.begin(), &*l.begin(), &*u.begin());
        ss(&*x.begin(), m, n, &*fv.begin());
        if (print != results)
        {
          opt.print(of);
          print_ss_par_1(&*x.begin(), n);
          eval_L(L, SS);
          of << "value of L is " << setw(15)
             << L;
          of << "  and of SS is " << setw(15)
             << SS << endl << endl;
        }
        else
        {
          of << "new parameters after " << iss << " iterations" << endl;
          eval_L(L, SS);
          of << "value of L is " << setw(15) << L
             << "  and of SS is " << setw(15)
             << SS << endl;
        }

        cout << "   var_comp       ";
        eval_var_comp();
        eval_L(L, SS);
        difL = fabs((L - Lold)/L);
        Lold = L;
        difv = 0.;
        iser = ser.begin();
        while(iser != ser.end())
        {
          series& t = *iser;
          if (!t.hide)
          {
            if (t.fl_sr2 != series::fixed)
              difv = max(difv, fabs((t.sr2 - t.sr2old)/t.sr2));
            if (t.fl_ga != series::fixed)
            {
              if (t.ga)
                difv = max(difv, fabs((t.ga - t.gaold)/t.ga));
              else
                difv = max(difv, fabs(t.ga - t.gaold));
            }
            if (t.fl_gb != series::fixed)
            {
              if (t.gb)
                difv = max(difv, fabs((t.gb - t.gbold)/t.gb));
              else
                difv = max(difv, fabs(t.gb - t.gbold));
            }
            t.sr2old = t.sr2;
            t.gaold = t.ga;
            t.gbold = t.gb;
          }
          iser++;
        }
        iter++;
      }
      of << endl << "Process ended after " << setw(4) << iter - 1
         << " iterations" << endl;
      opt.print(of);
      const matrix &jac = opt.jac();
      size_t k = 0;
      for (i = 0; i < n; ++i)
        for (j = 0; j <=i; ++j)
        {
          xjtj[k] = 0;
          for (size_t l = 0; l < m; ++l)
            xjtj[k] += jac(l, i)*jac(l, j);
          ++k;
        }
      of << endl;
      of << "difL " << setw(wdth) << difL;
      of << ", difv " << setw(wdth) << difv << endl;
      of << "value of L is " << setw(15) << L;
      of << "  and of SS is " << setw(15) << SS
         << endl << endl;
      cout << endl;
  // estimate confidence limits
      int ier;
      dpptrin(n, &*xjtj.begin(), ier);
      if (ier)
      {
        warning("can't invert (X'*X) - skipped confidence intervals");
        of << "Final values of parameters" << endl;
        for (j = 0; j < n; j++)
          of << setw(wdth) << xnm[j].c_str() << setw(wdth) << x[j]*sc[j] 
            << endl;
        of << endl;
      }
      else
      {
        of << endl << "Parameter +/- standard deviation" << endl;
        i = 1;   //Fortran i
        for (j = 0; j < n; j++)
        {
          double sx = sqrt(xjtj[i*(i+1)/2 - 1]);
          of << setw(wdth) << xnm[j].c_str() 
            << setw(wdth) << x[j]*sc[j] << " +/-"
            << setw(wdth) << sx 
            << endl;
          i++;
        }
        of << endl << "Correlation matrix for free parameters" << endl;
        of << setw(wdth) << 1 << endl;
        for (i = 2; i <= n; i++)      // Fortran i
        {
          for (j = 1 ; j <= i - 1; j++)  // Fortran j
            of << setw(wdth)
               << xjtj[i*(i-1)/2+j-1]/sqrt(xjtj[i*(i+1)/2-1]*xjtj[j*(j+1)/2-1]);
            of << setw(wdth) << 1 << endl;
        }
      }
    }
    else
    {
      if (m)
      {
        iser = ser.begin();
        while(iser != ser.end())
        {
          series& t = *iser;
          if (!t.hide)
            set_sums_1(t);
          iser++;
        }
        cout << "   var_comp       ";
        eval_var_comp();
        cout << endl;
        eval_L(L, SS);
      }
    }

  // print all the points
    if (print == results)
      if (m)
        print_all_points();
    post_analysis();
  // end calculations
  }
  catch (gError &t)
  {
    cout << endl << t.message << endl;
    return 1;
  }
  catch (...)
  {
    cout << endl << "Unknown exception in main" << endl;
    return 2;
  }
  return 0;
}

