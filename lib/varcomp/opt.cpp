/*
Copyright (C) 1994 - 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#include <iomanip>
#include "opt.h"

void (*optimizer::InitSGML[])() = {tensolve::AddSGML, 
  LBFGSB::AddSGML,
  Hooke::AddSGML,
#ifdef INCLUDE_ZXSSQ
  ZXSSQ::AddSGML, 
#endif  
  0};

string optimizer::algorithm = "tensolve";
map<string, SGML> optimizer::opt;

double tensolve::ssbest;
vec_double tensolve::xbest;

optimizer::optimizer()
{
  if (opt.empty())
    for (size_t i = 0; InitSGML[i]; ++i)
      (*InitSGML[i])();
}

void optimizer::ReadSGML(const SGML &e)
{
  algorithm = e.FindString("algorithm", "tensolve");
  max_iter = e.FindInt("NBigIterations", max_iter);
  SGML el;
  parser p(e);
  while (!p.eof())
  {
    p.GetSGML(el);
    if (el.name == algorithm)
    {
      SGML &o = opt[algorithm];
      for (map_string_string_ci i = el.attr.begin(); i != el.attr.end(); ++i)
      {
        map_string_string_i j = o.attr.find((*i).first);
        if (j != o.attr.end())
          (*j).second = (*i).second;
      }
    }
  }
}

void optimizer::WriteSGML(ostream &out, size_t shift)
{
  out << PutTab(shift) << "<optimizer algorithm=" << algorithm << endl;
  out << PutTab(shift + 1) << "NBigIterations=" << max_iter << ">" 
    << endl;
  for (map<string, SGML>::iterator i = opt.begin(); i != opt.end(); ++i)
    (*i).second.write(out, shift + 1);
  out << PutTab(shift) << "</optimizer>" << endl;
}

void optimizer::init(size_t n, size_t m)
{
  if (algorithm == "LBFGSB")
    ptr = new LBFGSB;
  else if (algorithm == "Hooke")
    ptr = new Hooke;
#ifdef INCLUDE_ZXSSQ
  else if (algorithm == "ZXSSQ")
    ptr = new ZXSSQ;
#endif
  else
    ptr = new tensolve;
  ptr->set(n, m, opt[algorithm]);
}

void tensolve::AddSGML()
{
  static const char *str = "<tensolve NIterations=15  MachineEps=1e-10 \
    GradTol=1e-5 StepTol=1e-5 FunctionTol=1e-7 method=Newton \
    strategy=LineSearch MaxStep=10000 TrustRegionRadius=-1 \
    ></tensolve>";
  istringstream in(str);
  SGML e;
  e.read(in);
  optimizer::opt["tensolve"] = e;
}

//int i_dnnt(double *); //from g2c.a

void tensolve::set(size_t n_, size_t m_, const SGML &e)
{
  termcd = 0;
  n = n_;
  m = m_;
  maxm = m + n + 2;
  maxn = n + 2;
  double tmp = sqrt(double(n));
  sqrn = i_dnnt(&tmp);
  typex.resize(n);
  fill(typex.begin(), typex.end(), 1.);
  typef.resize(m);
  fill(typef.begin(), typef.end(), 1.);
  itnlim = e.FindInt("NIterations");
  epsm = e.FindDouble("MachineEps");
  gradtl = e.FindDouble("GradTol");
  steptl = e.FindDouble("StepTol");
  ftol = e.FindDouble("FunctionTol");
  if (e.FindString("method") == "Newton")
    method = 0;
  else
    method = 1;
  if (e.FindString("strategy") == "LineSearch")
    global = 0;
  else
    global = 1;
  stepmx = e.FindDouble("MaxStep");
  dltin = e.FindDouble("TrustRegionRadius");
  jacflg = 0;
/*  
  itnlim = 150;
  epsm = DBL_EPSILON;
  gradtl = pow(epsm, 1./3);
  ftol = steptl = pow(epsm, 2./3.);
  method = 1;
  global = 0;
  stepmx = 1000.;
  dlt = -1.;
*/  
  ipr = 6;
  msg = 8;
  x.resize(sqrn);
  typxu.resize(sqrn);
  xpls.resize(sqrn);
  gpls.resize(sqrn);
  a.resize(sqrn*sqrn);
  wrk.resize(sqrn*sqrn);
  dfn.resize(maxm);
  wrk1.resize(maxm);
  wrk2.resize(maxm);  
  wrk3.resize(maxm);  
  wrk4.resize(maxm);  
  wrk5.resize(maxm);
  fq.resize(maxm);
  fqq.resize(maxm);
  fc.resize(maxm);
  fhat.resize(maxm);
  anls.resize(maxm*sqrn);
  fv.resize(maxm*sqrn);
  aja.resize(maxm,n);
  dxn.resize(maxn);
  dn.resize(maxn);
  dt.resize(maxn);
  df.resize(maxn);
  d.resize(maxn);
  gbar.resize(maxn);
  dbar.resize(maxn);
  dbarp.resize(maxn);
  s.resize(maxn*sqrn);
  shat.resize(maxn*sqrn);
  curpos.resize(maxn);
  pivot.resize(maxn);
  pbar.resize(maxn);
  fp.resize(m);
  gp.resize(n);
  xp.resize(n);
  xbest.resize(n);
}

void tensolve::ssn(double *x, double *f, int *m, int *n)
{
  ss(x, *m, *n, f);
  double sum = 0.;
  for (int i = 0; i < *m; ++i)
    sum += f[i]*f[i];
  if (sum < ssbest)
  {
    ssbest = sum;
    for (int i = 0; i < *n; ++i)
      xbest[i] = x[i];
  }
}

void tensolve::solve(double *x_, double *l, double *u)
{
  ssbest = HUGE_VAL;
  termcd = 0;
  dlt = dltin;
  fill(dfn.begin(), dfn.end(), 1.);
  fill(dxn.begin(), dxn.end(), 1.);
  tsnesv_(&maxm, &maxn, &sqrn, x_, &m, &n, &*typex.begin(), &*typef.begin(), 
      &itnlim,
      &jacflg, &gradtl, &steptl, &ftol, &method, &global, 
      &stepmx, &dlt, &ipr, &*x.begin(), &*typxu.begin(), &*xpls.begin(), 
      &*gpls.begin(), &*a.begin(), &*wrk.begin(), &*dfn.begin(), 
      &*wrk1.begin(), &*wrk2.begin(), &*wrk3.begin(), &*wrk4.begin(), &*wrk5.begin(),
      &*fq.begin(), &*fqq.begin(), &*fc.begin(), &*fhat.begin(), 
      &*anls.begin(), &*fv.begin(), &*aja.begin(), &*dxn.begin(), &*dn.begin(), 
      &*dt.begin(), &*df.begin(), &*d.begin(), &*gbar.begin(), &*dbar.begin(), 
      &*dbarp.begin(), 
      &*s.begin(), &*shat.begin(), &*curpos.begin(), &*pivot.begin(), &*pbar.begin(),
      &epsm, &sqrn, &ssn, &tsdumj_, &msg, &*xp.begin(), &*fp.begin(), 
      &*gp.begin(), &termcd);
  for (int i = 0; i < n; ++i)
    x_[i] = xbest[i];
//    x_[i] = xp[i];
// somehow tensolve might return not the best solution
// I do not know what is wrong
// this is a quick but not that good fix with the use of static members
}

void tensolve::print(ostream& out)
{
  out << "tensolve report" << endl;
  switch (termcd)
  {
    case 0:
      out << "illegal input" << endl;
      break;
    case 1:
      out << "function tolerance reached" << endl;
      break;
    case 2:
      out << "gradient tolerance reached" << endl;
      break;
    case 3:
      out << "Successive iterates within step tolerance" << endl;
      break;
    case 4:
      out << "Last global step failed to locate a point lower than XP" << endl;
      break;
    case 5:
      out << "Iteration limit exceeded" << endl;
      break;
    case 6:
      out << "Five consecutive steps of length STEPMX have been taken" << endl;
      break;
    default:
      out << "Unknown code" << endl;
      break;
  }
  out << "gradient is" << endl;
  for (size_t j = 0; j < gp.size(); j++)
    out << setw(wdth) << xnm[j].c_str() << setw(wdth) << gp[j] << endl;
  out << endl;
}

void LBFGSB::AddSGML()
{
  static const char *str = "<LBFGSB NIterations=100 memory=10 small=0.1 \
    StepTol=1e-5 factor=1e7 pgtol=1e-8 debug=0 step=1e-5 \
    ></LBFGSB>";
  istringstream in(str);
  SGML e;
  e.read(in);
  optimizer::opt["LBFGSB"] = e;
}

void LBFGSB::set(size_t n_, size_t m_, const SGML &e)
{
  n = n_;
  m = m_;
  xjac.resize(m, n);
  f.resize(m);
  opt_m = e.FindInt("memory");
  work.resize((2*opt_m + 4)*n + 11*opt_m*opt_m + 8*opt_m);
  nbd.resize(n);
  xold.resize(n);
  grad.resize(n);
  fn2.resize(m);
  iwa.resize(3*n);
  iter = e.FindInt("NIterations");
  factor = e.FindDouble("factor");
  dsig = e.FindDouble("StepTol");
  pgtol = e.FindDouble("pgtol");
  if (e.FindInt("debug"))
    iprint = 1;
  else
    iprint = -1;
  step = e.FindDouble("step");
  small = e.FindDouble("small");
}

void LBFGSB::solve(double *x, double *l, double *u)
{
  for (int j = 0; j < n; ++j)
  {
    xold[j] = x[j];
    nbd[j] = 0;
    if (l[j] != -HUGE_VAL)
      nbd[j] += 1;
    if (u[j] != HUGE_VAL)
      nbd[j] += 2;
    if (nbd[j] == 2)
      nbd[j] = 3;
    else if (nbd[j] == 3)
      nbd[j] = 2;
  }
  bool central = false;
  size_t i, j;
  double oldx, dx;
  s_copy(task, "START", 60L, 5L);
  long nn = n;
  while (1)
  {
    setulb_(&nn, &opt_m, x, l, u, &*nbd.begin(), &ssq, &*grad.begin(), 
        &factor, &pgtol, &*work.begin(), &*iwa.begin(), task, 
        &iprint, csave, lsave, isave, dsave);
    if (s_cmp(task, "FG", 2L, 2L) == 0)
    {
      ss(x, m, n, &*f.begin());
      ssq = 0.;
      for (i = 0; i < m; ++i)
      {
        ssq += f[i]*f[i];
      }
      for (j = 0; j < n; ++j)
      {
        oldx = x[j];
        dx = step*(x[j] + step);
        x[j] += dx;
        ss(x, m, n, &*xjac.begin() + m*j);
        grad[j] = 0.;
        if (central)
        {
          x[j] = oldx - dx;
          ss(x, m, n, &*fn2.begin());
          for (i = 0; i < m; ++i)
          {
            xjac(i, j) = (xjac(i, j) - fn2[i])/(dx*2.);
            grad[j] += 2.*f[i]*xjac(i, j);
          }
        }
        else
        {
          for (i = 0; i < m; ++i)
          {
            xjac(i, j) = (xjac(i, j) - f[i])/dx;
            grad[j] += 2.*f[i]*xjac(i, j);
          }
        }
        x[j] = oldx;
        if (!central && fabs(grad[j]) < small)
          central = true;
      }
    }
    else if (s_cmp(task, "NEW_X", 5L, 5L) == 0)
    {
      if (iss > iter)
        s_copy(task, "STOP - iter", 60L, 11L);
      bool no_diff = true;
      double dif;
      for (size_t i = 0; i < n; ++i)
      {
        dif = x[i] - xold[i];
        if (x[i] != 0)
          dif /= x[i];
        if (fabs(dif) > dsig)
          no_diff = false;
        xold[i] = x[i];
      }
      if (no_diff)
        s_copy(task, "STOP - diff", 60L, 11L);
      continue;
    }
    else
      break;
  }
}

void LBFGSB::print(ostream &of)
{
  of << "LBFGSB report" << endl;
  for (size_t i = 0; i < 60; ++i)
    of << task[i];
  of << endl;
  of << "the number of current iteration is " << isave[29] << endl;
  of << "the total number of function and gradient evaluations is " << 
    isave[33] << endl;
  of << "the number of free variables is " << isave[37] << endl;
  of << "the number of active constraints is " << isave[38] << endl;
  of << "the accumulated time spent on searching for Cauchy points is " 
    << dsave[6] << endl;
  of << "the accumulated time spent on subspace minimization is " 
    << dsave[7] << endl;
  of << "the accumulated time spent on line search is " << dsave[8] << endl;
  of << "the infinity norm of the projected gradient is " << dsave[12] 
    << endl;
  of << "gradient is" << endl;
  for (size_t j = 0; j < grad.size(); j++)
    of << setw(wdth) << xnm[j].c_str() << setw(wdth) << grad[j] << endl;
  of << endl;
}

void Hooke::AddSGML()
{
  static const char *str = "<Hooke NIterations=30 rho=0.5 \
    StepTol=0.1 \
    ></Hooke>";
  istringstream in(str);
  SGML e;
  e.read(in);
  optimizer::opt["Hooke"] = e;
}

void Hooke::set(size_t n_, size_t m_, const SGML &e)
{
  n = n_;
  m = m_;
  iter = e.FindInt("NIterations");
  eps = e.FindDouble("StepTol");
  rho = e.FindDouble("rho");
  fv.resize(m);
  xjac.resize(m, n);
  z.resize(n);
  delta.resize(n);
  xbefore.resize(n);
  newx.resize(n);
  endpt.resize(n);
}

void Hooke::print(ostream &of)
{
  of << "The search has been done by Hooke in " << ii << " iterations" << endl;
}

#ifdef INCLUDE_ZXSSQ

void ZXSSQ::AddSGML()
{
  static const char *str = "<ZXSSQ NIterations=100 MachineSignificantDigits=10 \
    SignificantDigits=5 FunctionTol=1e-7 \
    ></ZXSSQ>";
  istringstream in(str);
  SGML e;
  e.read(in);
  optimizer::opt["ZXSSQ"] = e;
}

void ZXSSQ::set(size_t n_, size_t m_, const SGML &e)
{
  n = n_;
  m = m_;
  int level = 0, levold;
  uerset(level,levold);
  istream *in;
  ostream *out = &off;
  ugetio(3, in, out);
  ixjac = m;
  maxfn = e.FindInt("NIterations");
  imsl::zxssq_machine_eps = e.FindInt("MachineSignificantDigits");
  nsig = e.FindInt("SignificantDigits");
  eps = e.FindDouble("FunctionTol");
  delta = 0;
  iopt = 1;
  work.resize(5*n + 2*m + (n+1)*n/2);
  f.resize(m);
  xjtj.resize(n*(n+1)/2);
  xjac.resize(m, n);
}

void ZXSSQ::print(ostream &of)
{
  of << "ZXSSQ  reports:" << endl;
  if (ier == 38)
    of << "jacobian is zero - the solution x is a stationary point"
           << endl;
  switch (infer)
  {
    case 0:
      of << "convergence failed - " << endl;
      switch (ier)
      {
        case 129:
          of << "\ta singularity was detected \
                 in the jacobian and recovery failed" << endl;
          break;
        case 130:
          of << "\tat least one of m, n, iopt, parm(1), or parm(2) \n\
                  was specified incorrectly" << endl;
          break;
        case 131:
          of << "\tthe marquardt parameter exceeded parm(3)" << endl;
          break;
        case 132:
          of << "\tafter a successful recovery from a singular jacobian, \n\
                 the vector x has cycled back to the first \
                 singularity" << endl;
          break;
        case 133:
          of << "\tmaxfn was exceeded" << endl;
          break;
        default:
          of << "\tit is not clear what has happened" << endl;
      }
      break;
    case 1:
      of << "criterion NSIG was satisfied" << endl;
      break;
    case 2:
      of << "criterion EPS was satisfied" << endl;
      break;
    case 3:
      of << "criterions NSIG and EPS were satisfied" << endl;
      break;
    case 4:
      of << "criterion DELTA was satisfied" << endl;
      break;
    case 5:
      of << "criterions NSIG and DELTA were satisfied" << endl;
      break;
    case 6:
      of << "criterions EPS and DELTA were satisfied" << endl;
      break;
    case 7:
      of << "criterions NSIG, EPS and DELTA were satisfied" << endl;
      break;
    default:
      of << "it is not clear what has happened" << endl;
  }
  of << "Gragient            " << setw(wdth)
     << work[0] << endl;
  of << "Function evaluations" << setw(wdth)
     << work[1] << endl;
  of << "Significant digits  " << setw(wdth)
     << work[2] << endl;
  of << "Marquardt parameter " << setw(wdth)
     << work[3] << endl;
  of << "Number iterations   " << setw(wdth)
     << work[4] << endl;
}
#endif

