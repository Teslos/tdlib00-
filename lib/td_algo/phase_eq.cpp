/*
Copyright (C) 1998, 1999, 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                       http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#include "phase_eq.h"
#include <iterator>

const char *PhaseEquilibrium::StatusOfVariableName[] = {"unknown", "constraint","HardConstraint", "dependent", "functional", "value", 0};

double PhaseEquilibrium::Tmin = 300.;
double PhaseEquilibrium::Tmax = 2000.;
double PhaseEquilibrium::eps = 1e-2;
double PhaseEquilibrium::penalty = 0.1;
SetLbfgsb PhaseEquilibrium::GlSet;

void RegisterPhaseEquilibrium()
{
  PhaseEquilibrium::GlSet.iter = 300;
  algorithm::RegisterType("PhaseEquilibrium", new PhaseEquilibrium);
}

void GuessByGrid(const VOF &fnc, long n, long np, double *x, const double *h,
                   double *fv, size_t N)
{
  double SSmin = 1e300, SS;
  vec_double xbest(n);
  vec_double xini(n);
  vec_double step(n);
  vec_double xfin(n);
  int i;
  for (i = 0; i < n; ++i)
  {
    step[i] = h[i]/N;
    xfin[i] = x[i] + h[i] + step[i]/2.;
    xbest[i] = xini[i] = x[i];
  }
  while (x[n - 1] < xfin[n - 1])
  {
    fnc(&n, x, fv, &np);
    SS = 0.;
    for (i = 0; i < n; ++i)
      SS = fv[i]*fv[i];
    if (SS < SSmin)
    {
      SSmin = SS;
      copy(x, x + n, xbest.begin());
    }
    for (i = 0; i < n; ++i)
    {
      x[i] += step[i];
        if (x[i] < xfin[i])
          break;
        else
        {
          if (i < n - 1)
            x[i] = xini[i];
        }
    }
  }
  copy(xbest.begin(), xbest.end(), x);
  fnc(&n, x, fv, &np);
}

void PhaseEquilibrium::fun(long *n, double *x, double *f, long *np) const
{
  if (debug)
  {
    *debug << "x";
    for (int i = 0; i < *n; ++i)
      *debug << setw(12) << x[i];
  }
  SetDependent();
  size_t ii = 0;
  if (eq.NCols())
  {
    vec_double_i ib = mu_basis.begin();
    vec_double_i im = mu.begin();
    vec_int_ci it = IsBasis.begin();
    for (size_t i = 0; i < sys.size(); ++i)
    {
      const vec_double &mui = sys[i].mu(Tp, vx[i]);
      for (size_t k = 0; k < mui.size(); ++k)
      {
        if (*it++)
          *ib++ = mui[k];
        else
          *im++ = mui[k];
      }
    }
    for (ii = 0; ii < eq.NCols(); ++ii)
    {
      f[ii] = mu[ii];
      for (size_t j = 0; j < mu_basis.size(); ++j)
        f[ii] -= eq(j, ii)*mu_basis[j];
    }
  }
  for (size_t j = 0; j < prop.size(); ++j)
  {
    double x;
    OutProperty(j, &x);
    f[ii++] = prop[j] - x;
  }
  if (debug)
  {
    *debug << "   f";
    for (int i = 0; i < *np; ++i)
      *debug << setw(12) << f[i];
    *debug << endl;
  }
}

void PhaseEquilibrium::solver() const
{
  if (!solved)
  {
    if (!IsXSet)
      SetX();
    if (debug)
    {
      *debug << "PhaseEquilibrium: ";
      for (size_t i = 0; i < sys.size(); ++i)
        *debug << sys[i].id() << " ";
      *debug << endl;
    }
    if (Neq == 0)
    {
      fmin = 0.;
      SetDependent();
      if (debug)
        *debug << "Neq is zero" << endl;
    }
    else if (Nx == 0)
    {
      fun(&Nx, &*x_.begin(), &*fv.begin(), &Neq);
      fmin = 0.;
      for (size_t i = 0; i < fv.size(); ++i)
        fmin += fv[i]*fv[i];
      fmin = sqrt(fmin);
      if (debug)
        *debug << "Nx is zero, fmin is " << fmin << endl;
    }
    else
      switch (IsEst)
      {
        case FullEst:
        {
          x_ = xest;
          fun(&Nx, &*x_.begin(), &*fv.begin(), &Neq);
          double f = 0;
          for (size_t i = 0; i < fv.size(); ++i)
            f += fv[i]*fv[i];
          f = sqrt(f);
          lbfgsb();
          if (fmin > eps)
          {
            fmin = f;
            x_ = xest;
          }
          else
            fun(&Nx, &*x_.begin(), &*fv.begin(), &Neq);
          if (debug)
            *debug << "full estimate, final fmin is " << fmin << endl;
          break;
        }
        case NoEst:
        {
          by_grid();
          lbfgsb();
          fun(&Nx, &*x_.begin(), &*fv.begin(), &Neq);
          if (debug)
            *debug << "there were no estimates, final fmin is " << fmin << endl;
          break;
        }
        case PartEst:
        {
          size_t j = 0;
          for (int i = 0; i < Nx; ++i)
            if (xest[i] == -1.)
              *st[iun[i]] = &x_[j++];
            else
              *st[iun[i]] = &xest[i];
          for (size_t i = 0; i < vc.size(); ++i)
            vc[i].SetPtr(st);
          pby_grid();
          plbfgsb();
          double f = fmin;
          j = pl.size() - 1;
          for (int i = Nx - 1; i >= 0; --i)
          {
            *st[iun[i]] = &x_[i];
            if (xest[i] == -1.)
              x_[i] = x_[j--];
            else
              x_[i] = xest[i];
          }
          for (size_t i = 0; i < vc.size(); ++i)
            vc[i].SetPtr(st);
          vec_double xtmp(x_);
          fun(&Nx, &*x_.begin(), &*fv.begin(), &Neq);
          if (debug)
            *debug << "now trying to find full solution " << fmin << endl;
          lbfgsb();
          if (fmin > eps)
          {
            fmin = f;
            x_ = xtmp;
          }
          else
            fun(&Nx, &*x_.begin(), &*fv.begin(), &Neq);
          if (debug)
            *debug << "partial estimate, final fmin is " << fmin << endl;
          break;
        }
      }
    solved = true;
    if (debug)
    {
      for (size_t i = 0; i < fix.size(); ++i)
        *debug << nm_st[i] << " " << fix[i] << " " << **st[i] 
          << endl;
      for (size_t i = 0; i < prop.size(); ++i)
      {
        double x;
        OutProperty(i, &x);
        *debug << name(i) << " " << x << endl;
      }
    }
    if (SaveSolution && fmin < eps)
      xest = x_;
  }
  return;
}

void PhaseEquilibrium::SetX() const
{
  int n = 0;
  for (int j = 0; j < Nx; ++j)
  {
    if (xest[j] == -1.)
      ++n;
    else if (xest[j] < l[j] || xest[j] > u[j])
    {
      xest[j] = -1.;
      ++n;
    }
  }
  if (n == Nx)
    IsEst = NoEst;
  else if (n == 0)
    IsEst = FullEst;
  else
  {
    ph.resize(n);
    pl.resize(n);
    pu.resize(n);
    size_t j = 0;
    for (int i = 0; i < Nx; ++i)
    {
      if (xest[i] == -1.)
      {
        ph[j] = h[i];
        pl[j] = l[i];
        pu[j] = u[i];
        j++;
      }
    }
    IsEst = PartEst;
  }
  IsXSet = true;
  if (debug)
  {
    for (size_t i = 0; i < sys.size(); ++i)
      *debug << sys[i].id() << " ";
    *debug << ", estimate is " << IsEst << endl;
    if (IsEst == PartEst)
      for (size_t i = 0; i < n; ++i)
        *debug << setw(10) << x_[i]
             << setw(10) << pl[i]
             << setw(10) << pu[i]
             << setw(10) << ph[i]
             << endl;
  }
}

void PhaseEquilibrium::SetDependent() const
{
  for (size_t i = 0; i < vc.size(); ++i)
    **st[i_vc[i]] = vc[i]();
  int j = 0;
  for (size_t i = 0; i < sys.size(); ++i)
  {
    if (i_dmf[i] != -1)
    {
      **st[i_dmf[i]] = 1.;
      for (size_t k = 0; k < sys[i].size(); ++k)
      {
        if (j != i_dmf[i])
          **st[i_dmf[i]] -= **st[j];
        ++j;
      }
    }
  }
}

void PhaseEquilibrium::SetPhases()
{
  if (sys.empty())
    throw gError("PhaseEquilibirum: no phases");
  nm_sys.resize(sys.size());
  vx.resize(sys.size());
  size_t nst = 2;
  for (size_t i = 0; i < sys.size(); ++i)
  {
    nst += sys[i].size();
    string name = sys[i].id();
    size_t num = 0;
    size_t first;
    char ch = '1';
    for (size_t j = 0; j < i; ++j)
    {
      if (sys[j].id() == name)
      {
        ++num;
        if (num == 1)
          first = j;
      }
    }
    if (num == 1)
      nm_sys[first] += "_1";
    if (num > 0)
    {
      ch += num;
      name.append("_").append(1, ch);
    }
    nm_sys[i] = name;
  }
  nm_st.resize(nst);
  nm_st2.resize(nst);
  fix.resize(nst);
  st.resize(nst);
  st_in.resize(nst);
  fill(st_in.begin(), st_in.end(), -1);
  st_in.back() = 1.; // p = 1 by default
  fill(fix.begin(), fix.end(), constraint);
  fix.back() = HardConstraint; // for p
  size_t j = 0;
  size_t k = 0;
  i_dmf.resize(sys.size());
  for (size_t i = 0; i < sys.size(); ++i)
  {
    fix[j] = dependent;
    i_dmf[i] = j;
    j += sys[i].size();
    vx[i].resize(sys[i].size());
    for (size_t l = 0; l < sys[i].size(); ++l)
    {
      nm_st[k] = string("x(") + nm_sys[i] + string(",")
        + sys[i].comps()[l].cf() + string(")");
      nm_st2[k] = string("x(") + nm_sys[i] + string(",") 
        + ObjToString(l + 1) + string(")");
      st[k] = &vx[i].ptr(l);
      ++k;
    }
  }
  nm_st[nst - 2] = nm_st2[nst - 2] = "T";
  nm_st[nst - 1] = nm_st2[nst - 1] = "p";
  st[nst - 2] = &Tp.ptrT();
  st[nst - 1] = &Tp.ptrp();
}

void PhaseEquilibrium::SetEq()
{
  vec_formula vf = sys[0].comps();
  for (size_t i = 1; i < sys.size(); ++i)
  {
    const vec_formula &vf2 = sys[i].comps();
    copy(vf2.begin(), vf2.end(), back_insert_iterator<vec_formula>(vf));
  }
  if (debug)
  {
    *debug << "all the formulas ";
    copy(vf.begin(), vf.end(), ostream_iterator<formula>(*debug, " "));
    *debug << endl;
  }
  GetNu ForMat;
  const vec_int &indx = ForMat.reset(vf);
  if (debug)
  {
    *debug << "basis is ";
    copy(indx.begin(), indx.end(), ostream_iterator<int>(*debug, " "));
    *debug << endl;
  }
  size_t N = indx.size();
  size_t M = vf.size() - indx.size();
  IsBasis.resize(vf.size());
  eq.resize(N, M);
  fv.resize(M);
  mu.resize(M);
  mu_basis.resize(N);
  fill(IsBasis.begin(), IsBasis.end(), 0);
  size_t ib = 0;
  size_t ir = 0;
  for (size_t i = 0; i < vf.size(); ++i)
  {
    if (ib < N && i == indx[ib])
    {
      IsBasis[i] = 1;
      ib++;
    }
    else
    {
      IsBasis[i] = 0;
      ForMat.nu(vf[i], &*eq.begin() + ir*N);
      ir++;
    }
  }
  if (debug)
  {
    *debug << "          ";
    vec_int_i it = IsBasis.begin();
    for (size_t i = 0; i < sys.size(); ++i)
      for (size_t j = 0; j < sys[i].size(); ++j)
        if (*it++ == 0)
          *debug << setw(10) << sys[i].comps()[j].cf().c_str();
    *debug << endl;
    size_t k = 0;
    it = IsBasis.begin();
    for (size_t i = 0; i < sys.size(); ++i)
      for (size_t j = 0; j < sys[i].size(); ++j)
        if (*it++)
        {
          *debug << setw(10) << sys[i].comps()[j].cf().c_str();
          for (size_t l = 0; l < eq.NCols(); ++l)
            *debug << setw(10) << eq(k, l);
          *debug << endl;
          ++k;
        }
  }
}

void PhaseEquilibrium::SetVar(const vec_double &lo, const vec_double &up, 
      const vec_convert &vctmp)
{
  Nx = count(fix.begin(), fix.end(), unknown);
  if (Nx > Neq)
    throw gError("PhaseEquilibrium: Nvar > Neq");
  x_.resize(Nx);
  xest.resize(Nx);
  h.resize(Nx);
  u.resize(Nx);
  l.resize(Nx);
  nbd.resize(Nx);
  iun.resize(Nx);
  size_t ivc = count(fix.begin(), fix.end(), functional);
  vc.resize(ivc);
  i_vc.resize(ivc);
  fill(i_dmf.begin(), i_dmf.end(), -1);
  size_t j = 0;
  size_t k = 0;
  size_t ip = 0;
  size_t ix = 0;
  for (size_t i = 0; i < fix.size(); ++i)
  {
    switch (fix[i])
    {
      case unknown:
        xest[j] = st_in[i];
        l[j] = lo[i];
        u[j] = up[i];
        h[j] = u[j] - l[j];
        nbd[j] = 2;
        iun[j] = i;
        *st[i] = &x_[j];
        ++j;
        break;
      case HardConstraint:
      case constraint:
        **st[i] = st_in[i];
        break;
      case dependent:
        if (ip == sys.size())
          throw gError("PhaseEquilbirium: try to set T or P to Dependent Mole Fraction");
        if (i_dmf[ip] != -1)
          throw gError("PhaseEquilibrium: more than one Dependent Mole Fraction");
        i_dmf[ip] = i;
        break;
      case functional:
        vc[k] = vctmp[i];
        i_vc[k] = i;
        break;
    }
    ++ix;
    if (ip < sys.size() && ix == sys[ip].size())
    {
      ix = 0;
      ++ip;
    }
  }
  for (size_t k = 0; k < ivc; ++k)
  {
    vc[k].SetID(nm_st);
    vc[k].SetPtr(st);
  }
  if (debug)
  {
    for (size_t i = 0; i < sys.size(); ++i)
      *debug << sys[i].id() << " ";
    *debug << "          x         l         u         h" << endl;
    for (size_t j = 0; j < xest.size(); ++j)
      *debug << nm_st[iun[j]]
         << setw(10) << xest[j]
         << setw(10) << l[j]
         << setw(10) << u[j]
         << setw(10) << h[j]
         << endl;
  }
}

PhaseEquilibrium& PhaseEquilibrium::operator=(const PhaseEquilibrium &old)
{
  if (this == &old)
    return *this;
// from RefAlgorithm  
  FileName = old.FileName;
  debug = old.debug;
  debug_ = old.debug_;
  
  sys = old.sys;
  nm_sys = old.nm_sys;
  Tp = old.Tp;
  vx = old.vx;
  nm_st = old.nm_st;
  nm_st2 = old.nm_st2;
  fix = old.fix;
  st = old.st;
  size_t k = 0;
  if (st.size())
  {
    for (size_t i = 0; i < sys.size(); ++i)
      for (size_t l = 0; l < sys[i].size(); ++l)
        st[k++] = &vx[i].ptr(l);
    st[st.size() - 2] = &Tp.ptrT();
    st[st.size() - 1] = &Tp.ptrp();
  }
  IsXSet = false;
  solved = false;
  ThrowException = old.ThrowException;
  SaveSolution = old.SaveSolution;
  set = old.set;
  eq = old.eq;
  IsBasis = old.IsBasis;
  fv = old.fv;
  mu = old.mu;
  mu_basis = old.mu_basis;
  prop = old.prop;
  Neq = old.Neq;
  
  x_ = old.x_;
  Nx = old.Nx;
  xest = old.xest;
  h = old.h;
  u = old.u;
  l = old.l;
  nbd = old.nbd;
  iun = old.iun;
  for (size_t i = 0; i < Nx; ++i)
    *st[iun[i]] = &x_[i];
  vc = old.vc;
  for (size_t i = 0; i < vc.size(); ++i)
    vc[i].SetPtr(st);
  i_vc = old.i_vc;
  i_dmf = old.i_dmf;
  prop_in = old.prop_in;
  st_in = old.st_in;
  vd = old.vd;
  return *this;
}

void PhaseEquilibrium::clear()
{
  IsXSet = false;
  solved = false;
  sys.clear();
}

void PhaseEquilibrium::read(const SGML &e)
{
  clear();
  SetDebug(e);
  if (debug)
    set.iprint = 1;
  else
    set.iprint = -1;
  ThrowException = e.FindInt("ThrowException", 1);
  SaveSolution = e.FindInt("SaveSolution", 0);
  SGML el;
  SGML el2;
  parser p(e);
  p.GetSGML(el);
  el.compare("phases");
  {
    parser p2(el);
    phase ph;
    while (!p2.eof())
    {
      p2.GetSGML(el2);
// compatibility      
      if (el2.name.empty())
      {
        ph.find(p2.GetToken());
      }
      else
        ph.read(el2);
      if (ph.size())
        sys.push_back(ph);
    }
  }
  SetPhases();
  SetEq();
  vec_double up(st.size());
  vec_double lo(st.size());
  fill(up.begin(), up.end() - 2, 1.);
  fill(lo.begin(), lo.end() - 2, 0.);
  up[up.size() - 2] = Tmax;
  lo[lo.size() - 2] = Tmin;
  up[up.size() - 1] = 100000;
  lo[lo.size() - 1] = 0.;
  vec_convert vctmp(st.size());
  while (!p.eof())
  {
    p.GetSGML(el);
    StatusOfVariable status;
    double val;
    double lower, upper;
    convert cv;
    if (el.name == "state")
    {
      status = el.FindString("status");
      switch (status)
      {
        case unknown:
          val = el.FindDouble("value", -1.);
          lower = el.FindDouble("lower", HUGE_VAL);
          upper = el.FindDouble("upper", HUGE_VAL);
          break;
        case HardConstraint:
          val = el.FindDouble("value", -1);
          if (val == -1)
            throw gError("PhaseEquilibrium: HardConstraint is not defined");
        case constraint:
          val = el.FindDouble("value", -1);
          break;
        case dependent:
          break;
        case functional:
          {
            parser p3(el);
            p3.GetSGML(el2);
            el2.compare("convert");
            cv.read(el2);
            el.body.erase(0, p3.npos());
          }
          break;
        case value:
          val = el.FindDouble("value", HUGE_VAL);
          break;
      }
    }
// compatibility    
    else if (el.name == "unknown")
    {
      status = unknown;
      val = el.FindDouble("value", -1.);
      lower = el.FindDouble("lower", HUGE_VAL);
      upper = el.FindDouble("upper", HUGE_VAL);
    }
    else if (el.name == "constraint")
    {
      status = HardConstraint;
      val = el.FindDouble("value", -1.);
      if (val == -1)
        throw gError("PhaseEquilibrium: HardConstraint is not defined");
    }
    else if (el.name == "FreeVariable")
    {
      status = constraint;
      val = el.FindDouble("value", -1.);
    }
    else 
      throw gError(string("PhaseEquilbirium: unknown element - ")
          .append(el.name));
    size_t N;
    int ivar = GetOutputIndex(el.body, N);
    if (N != 1)
      throw gError(string("PhaseEquilibrium: wrong variable (N > 1) - ")
          .append(el.body));
    if (status != value && ivar >= 0)
      throw gError(string("PhaseEquilibrium: wrong variable for state - ")
          .append(el.body));
    if (status == value && ivar < 0)
      throw gError(string("PhaseEquilibrium: wrong variable for value - ")
          .append(el.body));
    if (ivar < 0)
      ivar = -(ivar + 1);
    switch (status)
    {
      case unknown:
        fix[ivar] = status;
        st_in[ivar] = val;
        if (lower != HUGE_VAL)
          lo[ivar] = lower;
        if (upper != HUGE_VAL)
          up[ivar] = upper;
        break;
      case HardConstraint:
      case constraint:
        fix[ivar] = status;
        st_in[ivar] = val;
        break;
      case dependent:
        fix[ivar] = status;
        break;
      case functional:
        fix[ivar] = status;
        vctmp[ivar] = cv;
        break;
      case value:
        if (ivar != prop.size())
          throw gError("PhaseEquilibrium: ivar != prop.size()");
        prop.push_back(val);
        break;
    }
  }
  Neq = eq.NCols() + prop.size();
  fv.insert(fv.end(), prop.size(), 0.);
  SetVar(lo, up, vctmp);
  prop_in = prop;
}

void PhaseEquilibrium::WriteBody(ostream& out, size_t shift) const
{
  out << PutTab(shift) << "<phases>" << endl;
  for (size_t i = 0; i < sys.size(); ++i)
    sys[i].write(out, shift + 1);
  out << PutTab(shift) << "</phases>" << endl;
  size_t ivc = 0;
  size_t iun = 0;
  size_t ip = 0;
  size_t ix = 0;
  for (size_t i = 0; i < fix.size(); ++i)
  {
    if (ip == sys.size() || sys[ip].size() != 1)
    {
      out << PutTab(shift) << "<state status=" << fix[i];
      switch (fix[i])
      {
        case unknown:
          if (i == fix.size() - 2)
          {
            if (l[iun] != Tmin)
              out << " lower=" << l[iun];
            if (u[iun] != Tmax)
              out << " upper=" << u[iun];
          }
          else
            out << " lower=" << l[iun] << " upper=" << u[iun];
          if (st_in[i] != -1.)
            out << " value=" << st_in[i];
          out << "> ";
          ++iun;
          break;
        case HardConstraint:
          out << " value=" << st_in[i] << "> ";
          break;
        case constraint:
          if (st_in[i] != -1.)
            out << " value=" << st_in[i];
          out << "> ";
          break;
        case dependent:
          out << "> ";
          break;
        case functional:
          out << ">" << endl;
          vc[ivc++].write(out, shift + 1);
          out << PutTab(shift);
          break;
      }
      out << name(-(i + 1)) << " </state>" << endl;
    }
    ++ix;
    if (ip < sys.size() && ix == sys[ip].size())
    {
      ix = 0;
      ++ip;
    }
  }
  for (size_t i = 0; i < prop.size(); ++i)
  {
    out << PutTab(shift) << "<state status=value";
    if (prop_in[i] != HUGE_VAL)
      out << " value=" << prop_in[i];
    out << "> " << name(i) << " </state>" << endl;
  }
}

int PhaseEquilibrium::GetInputIndex(const string& s) const
{
  string str = RemoveSpaces(s);
  if (str == "SaveSolution")
    return -10000;
  int i = SearchState(str);
  if (i != -1)
  {
    if (fix[i] == unknown || fix[i] == constraint)
      return -(i + 1);
  }
  else
  {
    for (size_t i = 0; i < prop.size(); ++i)
      if (str == vd[i].name)
        return i;
  }
  throw gError(string("PhaseEquilibrium: can not use for input - ")
    .append(str));
}

void PhaseEquilibrium::SetVariable(int i, const double &x) const
{
  if (i == -10000)
    SaveSolution = x;
  else if (i < 0)
  {
    size_t ii = -(i + 1);
    if (fix[ii] != unknown)
      **st[ii] = x;
    else
    {
      int j;
      for (j = 0; j < Nx; ++j)
        if (iun[j] == ii)
        {
          xest[j] = x;
          break;
        }
      if (j == Nx)
        throw gError("PhaseEquilibrium: internal error in SetVariable");
    }
  }
  else
    prop[i] = x;
  solved = false;
  IsXSet = false;
}

const string& PhaseEquilibrium::name(int i) const
{
  if (i == -10000)
  {
    static string str("SaveSolution");
    return str;
  }
  else if (i < 0)
    return nm_st[-(i + 1)];
  else
    return vd[i].name;
}

int PhaseEquilibrium::GetOutputIndex(const string& s, size_t &N) const
{
  string str = RemoveSpaces(s);
  N = 1;
  int i = SearchState(str);
  if (i != -1)
    return -(i + 1);
  for (size_t j = 0; j < vd.size(); ++j)
    if (str == vd[j].name)
    {
      N = vd[j].N;
      return j;
    }
  i = vd.size();
  vd.push_back(descriptor());
  vd[i].name = str;
  if (str == "state(all)")
  {
    vd[i].from = AllState;
    N = vd[i].N = st.size();
    return i;
  }
  size_t ipar = str.find("(");
  string prop;
  string phase;
  string name;
  if (ipar == string::npos)
  {
    prop = str;
    phase = name = "";
  }
  else
  {
    prop = str.substr(0, ipar);
    name = str.substr(ipar + 1, str.rfind(")") - ipar - 1);
    size_t iu = name.find(",");
    if (iu == string::npos)
    {
      phase = name;
      name = "";
    }
    else
    {
      phase = name.substr(0, iu);
      name = name.substr(iu + 1, string::npos);
    }
  }
  if (prop == "fmin")
  {
    vd[i].from = Residual;
    if (name.empty())
    {
      N = vd[i].N = 1;
      vd[i].ix = -2;
    }
    else if (name == "all")
    {
      N = vd[i].N = fv.size();
      vd[i].ix = -1;
    }
    else
    {
      vd[i].ix = atoi(name.c_str());
      if (vd[i].ix < 0)
        throw gError(string("PhaseEquilibrium: wrong output string - ")
            .append(str));
      N = vd[i].N = 1;
    }
    return i;
  }
  size_t iu = prop.find("_");
  string fun;
  if (iu == string::npos)
  {
    fun = prop;
    prop = "";
  }
  else
  {
    fun = prop.substr(0, iu);
    prop.erase(0, iu + 1);
  }
  vd[i].f = fun;
  if (prop.empty())
    vd[i].i = ::full;
  else
    vd[i].i = prop;
  if (phase.empty())
    vd[i].ip = 0;
  else
  {
    vd[i].ip = SearchString(nm_sys, phase);
    if (vd[i].ip == -1)
      throw gError(string("PhaseEquilibrium: no such phase ").append(phase));
  }
  if (name.empty())
  {
    vd[i].from = Property;
    N = vd[i].N = 1;
  }
  else
  {
    vd[i].from = Partial;
    if (name == "all")
    {
      vd[i].ix = -1;
      N = vd[i].N = sys[vd[i].ip].size();
    }
    else
    {
      vd[i].ix = SearchFormula(sys[vd[i].ip].comps(), name);
      if (vd[i].ix == -1)
        throw 
          gError(string("PhaseEquilibrium: not known name in ").append(str));
      N = vd[i].N = 1;
    }
  }
  return i;
}

void PhaseEquilibrium::OutNames(int i, vec_string &vs) const
{
  if (i < 0)
  {
    vs.push_back(nm_st[-(i + 1)]);
    return;
  }
  switch(vd[i].from)
  {
    case Property:
      vs.push_back(vd[i].name);
      return;
    case Partial:
      if (vd[i].ix == -1)
        for (size_t j = 0; j < sys[vd[i].ip].size(); ++j)
        {
          vs.push_back(vd[i].name);
          vs.back().replace(vs.back().find("all"), 3, 
              sys[vd[i].ip].comps()[j].cf());
        }
      else
        vs.push_back(vd[i].name);
      return;
    case AllState:
      vs.insert(vs.end(), nm_st.begin(), nm_st.end());
      return;
    case Residual:
      if (vd[i].ix == -1)
        for (size_t j = 0; j < fv.size(); ++j)
        {
          vs.push_back(vd[i].name);
          vs.back().replace(vs.back().find("all"), 3, 
              ObjToString(j + 1));
        }
      else
        vs.push_back(vd[i].name);
  }
}

void PhaseEquilibrium::OutProperty(int i, double *x_) const
{
  switch(vd[i].from)
  {
    case Property:
      *x_ = sys[vd[i].ip].Z(vd[i].f, vd[i].i, Tp, vx[vd[i].ip]);
      return;
    case Partial:
      if (vd[i].ix != -1)
        *x_ = sys[vd[i].ip].z(vd[i].f, vd[i].i, Tp, vx[vd[i].ip])[vd[i].ix];
      else
      {
        const vec_double &r 
          = sys[vd[i].ip].z(vd[i].f, vd[i].i, Tp, vx[vd[i].ip]);
        copy(r.begin(), r.end(), x_);
      }
      return;
  }
}

void PhaseEquilibrium::out(int i, double *x_) const
{
  solver();
  if (i >= 0)
    switch(vd[i].from)
    {
      case AllState:
        for (size_t ii = 0; ii < st.size(); ++ii)
          *x_++ = **st[ii];
        return;
      case Residual:
        if (vd[i].ix == -2)
          *x_ = fmin;
        else if (vd[i].ix == -1)
          copy(fv.begin(), fv.end(), x_);
        else
          *x_ = fv[vd[i].ix];
        return;
    }
  if (ThrowException && fmin > eps)
    throw CanNotCompute(fmin*penalty);
  if (i < 0)
  {
    *x_ = **st[-(i + 1)];
    return;
  }
  else
    OutProperty(i, x_);
}

