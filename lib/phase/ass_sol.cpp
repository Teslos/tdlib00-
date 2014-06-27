/*
Copyright (C) 1998, 1999, 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                       http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#include <iomanip>
#include <algorithm>
#include <toms.h>
#include <lapack.h>
#include "ass_sol.h"

double AssociatedSolutionBasis::tol = DBL_EPSILON*100.;
double AssSolManyReact::factr = 1000;
double AssSolManyReact::pgtol = sqrt(DBL_EPSILON);
int AssSolManyReact::lbfg_m = 10;

void RegisterAssociatedSolution()
{
  phase::RegisterType("AssociatedSolution", new AssSolOneReact);
  phase::RegisterType("AssociatedSolution", new AssSolManyReact);
  phase::RegisterType("AssociatedSolution", new AssociatedSolution);
#ifdef INCLUDE_VCS
  phase::RegisterType("associated_solution", new associated_solution);
#endif
}

void AssociatedSolution::read(const SGML &e)
{
  AssociatedSolutionBasis *ptr = new AssociatedSolutionBasis;
  ptr->read(e);
  RefPhase *obj;
  if (ptr->R == 1)
    obj = new AssSolOneReact(*ptr, e);
  else
    obj = new AssSolManyReact(*ptr, e);
  delete ptr;
  throw ObjFromProxy(obj);
}

void AssociatedSolutionBasis::read(const SGML &e)
{
  if (e.FindInt("equilibrated") == 0)
    equil = ::ref_ideal;
  else
    equil = ::full;
  parser p(e);
  SGML el;
  ref.read(p.GetSGML(el));
  internal.read(p.GetSGML(el));
  initialize();
}

void AssociatedSolutionBasis::initialize()
{
  C = ref.size();
  N = internal.size();
  if (N <= C)
    throw gError("AssociatedSolutionBasis: N <= C");
  R = N - C;
  clear();
  const vec_formula &ex = ref.comps();
  const vec_formula &in = internal.comps();
  GetNu basis(ex);
  if (basis.NOfComps() != C)
    throw gError("AssociatedSolutionBasis: dependent components in ref");
  basis.FormulaMatrix(in, A);
}

struct n_greater
{
  const vec_double &n;
  n_greater(const vec_double &n_) : n(n_) {}
  bool operator()(int i, int j)
  {
    return n[i] > n[j];
  }
};

size_t AssociatedSolutionBasis::SetOptimalBasis() const
{
  vec_int in(N);
  for (size_t i = 0; i < N; ++i)
    in[i] = i;
  n_greater srt(n); 
  stable_sort(in.begin(), in.end(), srt);
  copy(in.begin(), in.begin() + C, basis.begin());
  copy(in.begin() + C, in.end(), reac.begin());
  size_t ireac = 0;
  bool change = true;
  int info;
  while (change) 
  {
    change = false;
    for (size_t i = 0; i < C; ++i)
      copy(A.begin() + basis[i]*C, A.begin() + basis[i]*C + C, 
          Ac.begin() + i*C);
    dgetf2(C, C, &*Ac.begin(), &*ipiv.begin(), info);
    for (size_t i = 0; i < C; ++i)
      if (fabs(Ac(i, i)) < GetNu::tol)
      {
        if (ireac == R)
          throw gError ("AssociatedSolutionBasis: not enough components in internal solution");
        swap(basis[i], reac[ireac]);
        ireac++;
        change = true;
      }
  }
  for (size_t i = 0; i < R; ++i)
    copy(A.begin() + reac[i]*C, A.begin() + reac[i]*C + C, 
        Z_.begin() + i*C);
  dgetrs(C, C, R, &*Ac.begin(), &*ipiv.begin(), &*Z_.begin());
  for (size_t i = 0; i < Z_.size(); ++i)
    if (fabs(Z_[i]) < GetNu::tol)
      Z_[i] = 0.;
  size_t same = 0;
  for (size_t i = 0; i < C; ++i)
  {
    pair<vec_int_i, vec_int_i> result 
      = equal_range(in.begin(), in.end(), basis[i], srt);
    same += result.second - result.first - 1;
  }
  if (dbg)
  {
    *dbg << "same value as in the basis is " << same << endl;
    const vec_formula &ex = ref.comps();
    const vec_formula &vf = internal.comps();
    *dbg << "Formula Matrix" << endl << "            ";
    for (size_t i = 0; i < N; ++i)
    {
      const string &s = vf[i].cf();
      int k = 12 - s.size();
      for (int j = 0; j < k; ++j)
        *dbg << " ";
      *dbg << s << " ";
    }
    *dbg << endl;
    for (size_t i = 0; i < C; ++i)
    {
      const string &s = ex[i].cf();
      *dbg << s;
      int k = 12 - s.size();
      for (int j = 0; j < k; ++j)
        *dbg << " ";
      for (size_t j = 0; j < N; ++j)
        *dbg << setw(12) << A(i, j) << " ";
      *dbg << endl;
    }
    *dbg << "Optimal Basis" << endl;
    *dbg << "the initial order is ";
    for (size_t i = 0; i < N; ++i)
      *dbg << vf[in[i]].cf() << " ";  
    *dbg << endl;
    *dbg << "Z Matrix" << endl << "            ";
    for (size_t i = 0; i < R; ++i)
    {
      const string &s = vf[reac[i]].cf();
      int k = 12 - s.size();
      for (int j = 0; j < k; ++j)
        *dbg << " ";
      *dbg << s << " ";
    }
    *dbg << endl;
    for (size_t i = 0; i < C; ++i)
    {
      const string &s = vf[basis[i]].cf();
      *dbg << s;
      int k = 12 - s.size();
      for (int j = 0; j < k; ++j)
        *dbg << " ";
      for (size_t j = 0; j < R; ++j)
        *dbg << setw(12) << Z_(i, j) << " ";
      *dbg << endl;
    }
  }
  return same;
}

void AssociatedSolutionBasis::SetStoichiometry() const
{
  for (size_t i = 0; i < R; ++i)
  {
    DelNu[i] = 1.;
    for (size_t j = 0; j < C; ++j)
      DelNu[i] -= Z_(j, i);
  }
  fill(Nu.begin(), Nu.end(), 0.);
  for (size_t i = 0; i < R; ++i)
  {
    for (size_t j = 0; j < C; ++j)
       Nu(basis[j], i) = -Z_(j, i);
    Nu(reac[i], i) = 1.;
  }
  fill(AcInv.begin(), AcInv.end(), 0.);
  for (size_t i = 0; i < C; ++i)
    AcInv(i, i) = 1.;
  dgetrs(C, C, C, &*Ac.begin(), &*ipiv.begin(), &*AcInv.begin());
  if (dbg)
  {
    *dbg << "AcInv" << endl;
    for (size_t i = 0; i < C; ++i)
    {
      for (size_t j = 0; j < C; ++j)
        *dbg << setw(12) << AcInv(i, j) << " ";
      *dbg << endl;
    }
  }
  for (size_t i = 0; i < C; ++i)
  {
    SumAcInv[i] = 0.;
    for (size_t j = 0; j < C; ++j)
      SumAcInv[i] += AcInv(j, i);
  }
}

double AssSolOneReact::f(const double &y) const
{
  eta[0] = y;
  SetMoleFractions();
  const vec_double &z = internal.z(::G, equil, *Tp, vy);
  double Sum = 0.;
  for (size_t i = 0; i < C; ++i)
    Sum -= z[basis[i]]*Z_(i, 0);
  Sum += z[reac[0]];
  if (debug)
  {
    *debug << setw(12) << y << setw(12) << Sum;
    for (size_t i = 0; i < vint.size(); ++i)
      *debug << setw(12) << vint[i];
    *debug << endl;
  }
  return Sum;
}

const vec_double& AssSolOneReact::yeq(const StateTp &Tp_, const StateX &x) const
{
  if (debug)
  {
    *debug << "AssSolOneReact: T is " << Tp_.T() 
      << " p is " << Tp_.p() << " x is ";
    for (size_t i = 0; i < x.size(); ++i)
      *debug << x[i] << " ";
    *debug << endl;
    *debug << "     eta           f          vy" << endl;
  }
  Tp = &Tp_;
  SetNo(x);
  b = HUGE_VAL;
  a = 0.;
  double t;
  for (size_t i = 0; i < C; ++i)
    if (Z_(i, 0) > 0.)
    {
      if (no[i] < 0)
        throw gError("AssSolOneReact: - no[i] < 0");
      t = no[i]/Z_(i, 0);
      if (t < b)
        b = t;
    }
    else if (Z_(i, 0) < 0.)
    {
      t = no[i]/Z_(i, 0);
      if (t > a)
        a = t;
    }
  if (debug)
    *debug << "a is " << a << " b is " << b << endl;
  if (fabs(b - a) < AssociatedSolutionBasis::tol)
  {
    eta[0] = (b + a)/2.;
    SetMoleFractions();
    if (debug)
      *debug << "eta is equal to zero" << endl << endl;
    return vint;
  }
  try
  {
    eta[0] = rroot(makeFunctor((D2D*)0, *this, &AssSolOneReact::f), a, b);
  }
  catch (err_rroot)
  {
    eta[0] = a;
    SetMoleFractions();
    double ga = internal.Z(::G, equil, *Tp, vy); 
    eta[0] = b;
    SetMoleFractions();
    double gb = internal.Z(::G, equil, *Tp, vy); 
    if (gb < ga)
      eta[0] = b;
    else
      eta[0] = a;
    if (debug)
    {
      *debug << "there were no roots: f[" << a 
        << "]=" << ga << " f[" << b << "]=" << gb << endl;
      *debug << "the choice is " << eta[0] << endl;
    }
  }
  catch (...)
  {
    throw gError("AssSolOneReact: not known error - no solution");
  }
  SetMoleFractions();
  if (debug)
    *debug << endl;
  return vint;
}

void AssSolManyReact::LinearProgramming(const StateTp &Tp, 
    const StateX &x) const
{
  int n_ = 2*N + 1;
  int m1 = 0;
  int m2 = C;
  double s;
  matrix ps(N, n_);
  for (size_t i = 0; i < N; ++i)
    for (size_t j = 0; j < N; ++j)
      if (i == j)
        ps(i, j) = 1.;
      else
        ps(i, j) = 0.;
  for (size_t j = 0; j < N; ++j)
    ps(j, N) = 1./double(N);
  for (size_t i = 0; i < N; ++i)
    for (size_t j = 0; j < N; ++j)
      ps(j, N + i + 1) = (ps(j, i) + ps(j, N))/2.;
  Cmatrix A_(C, n_);
  for (size_t i = 0; i < N; ++i)
    for (size_t j = 0; j < C; ++j)
      A_(j, i) = A(j, i);
  for (size_t i = 0; i < N; ++i)
    vy[i] = 0.;
  vec_double c(n_);
  const vec_double &c_ = internal.z(::G, ::ref, Tp, vy); //vy  = 0 is OK
  copy(c_.begin(), c_.end(), c.begin());
  StateX vy_(N);
  for (size_t i = N; i < n_; ++i)
  {
    for (size_t j = 0; j < C; ++j)
    {
      A_(j, i) = 0.;
      for (size_t k = 0; k < N; ++k)
        A_(j, i) += A_(j, k)*ps(k, i);
    }
    for (size_t k = 0; k < N; ++k)
      vy_.ptr(k) = &ps(k, i);
    c[i] = internal.Z(::G, equil, Tp, vy_);
  }
  vec_double b(C);
  for (size_t i = 0; i < C; ++i)
  {
    b[i] = x[i];
    if (b[i] < 0.)
      b[i] = 0.;
  }
// because minit search maximum
  for (size_t i = 0; i < n_; ++i)
    c[i] = -c[i];
  vec_double dsol(m1+m2);
  int ier;  
  vec_double psol(n_);
  lp.zx3lp(A_, b, c, n_, m1, m2, s, psol, dsol, ier);
  fill(n.begin(), n.end(), 0.);
  for (size_t i = 0; i < n_; ++i)
    if (psol[i] > 0)
    {
      for (size_t j = 0; j < N; ++j)
        n[j] += psol[i]*ps(j, i);
    }
  if (debug)
  {
    *debug << "Linear Programming Estimates" << endl;
    *debug << "-Go are ";
    copy (c.begin(), c.end(), ostream_iterator<double>(*debug, " "));
    *debug << endl;
    *debug << "psol are ";
    copy (psol.begin(), psol.end(), ostream_iterator<double>(*debug, " "));
    *debug << endl;
    *debug << "n are ";
    copy (n.begin(), n.end(), ostream_iterator<double>(*debug, " "));
    *debug << endl;
  }
  if (ier > 70)
    throw gError(string("AssSolManyReact: zxlp3 - ")
        .append(ObjToString(ier)));
}

#include <f2c.h>
#include <f2clib.h>

void AssSolManyReact::fun(long n, const StateTp &Tp, double &f) const
{
  if (n != R)
    for (size_t i = 0; i < n; ++i)
      *(vptr[i]) = eta1[i];
  SetMoleFractions();
  const vec_double &z = internal.z(::G, equil, Tp, vy);
  size_t ii = 0;
  for (size_t i = 0; i < R; ++i)
  {
    if (incl[i])
    {
      e[ii] = 0.;
      for (size_t k = 0; k < C; ++k)
        e[ii] -= z[basis[k]]*Z_(k, i);
      e[ii] += z[reac[i]];
      ++ii;
    }
  }
  f = fabs((*ntot))*internal.Z(::G, equil, Tp, vy);
  if (debug)
  {
    *debug << "f is " << f << " G/RT is " 
      << f/global::R/Tp.T() << endl;  
    if (n == R)
    {
      *debug << "eta[" << R << "] is ";
      for (size_t i = 0; i < R; ++i)
        *debug << setw(12) << eta[i];
    }
    else
    {
      *debug << "eta[" << n << "] is ";
      for (size_t i = 0; i < n; ++i)
        *debug << setw(12) << eta1[i];
    }
    *debug << endl;
    *debug << "grad[" << n << "] is ";  
    for (size_t i = 0; i < n; ++i)
      *debug << setw(12) << e[i];
    *debug << endl;
    *debug << "vy[" << N + 1 << "] is ";  
    for (size_t i = 0; i < vint.size(); ++i)
      *debug << setw(12) << vint[i];
    *debug << endl;
  }
}

const vec_double& AssSolManyReact::yeq(const StateTp &Tp, const StateX &x) const
{
  if (debug)
  {
    *debug << "AssSolManyReact: T is " << Tp.T() 
      << " p is " << Tp.p() << " x is ";
    for (size_t i = 0; i < x.size(); ++i)
      *debug << x[i] << " ";
    *debug << endl;
  }
  LinearProgramming(Tp, x);
  size_t change_basis = 0;
  bool finish = false;
  while (!finish && change_basis < 3)
  {
    bool exit;
    size_t iter;
    do 
    {
      iter = SetOptimalBasis();
      SetStoichiometry();
      SetNo(x);
      exit = true;
      for (size_t i = 0; i < C; ++i)
        if (no[i] < 0.)
        {
          if (iter)
          {
            exit = false;
            n[basis[i]] -= DBL_EPSILON;
            break;
          }
          else
            cout << "AssSolManyReact: negative no[" << i << "]=" << no[i] 
              << endl;
        }
    } 
    while (!exit && iter);
    char task[60];
    char csave[60];
    long lsave[4];
    long isave[44];
    double dsave[29];
    double f;
    long m = lbfg_m;
    long n_ = R;
    vec_double l(n_);
    vec_double u(n_);
    vector<long> nbd(n_);
    vec_double wa((2*m + 4)*n_ + 11*m*m + 8*m);
    vector<long> iwa(3*n_);
    long iprint;
    if (debug)
      iprint = 1;
    else
      iprint = -1;
    size_t k = 0;
    for (size_t i = 0; i < R; ++i)
    {
      l[k] = 0.;
      nbd[k] = 1; // only lower bounds
  //    nbd[k] = 2;
      u[k] = HUGE_VAL;
      double t;
      for (size_t j = 0; j < C; ++j)
        if (Z_(j, i) > 0. && (t = no[j]/Z_(j, i)) < u[k])
          u[k] = t;
      if (fabs(u[k] - l[k]) < AssociatedSolutionBasis::tol)
      {
        incl[i] = 0;
        if (n_ == R)
        {
          eta1.resize(k);
          vptr.resize(k);
          for (size_t ii = 0; ii < k; ++ii)
          {
            eta1[ii] = eta[ii];
            vptr[ii] = &eta[ii];
          }
        }
        --n_;
        eta[i] = (u[k] - l[k])/2.;
        if (debug)
          *debug << "eta[" << i << "] is equal to zero " 
            << eta[i] << " l is " << l[k] << " and h is " << u[k] << endl;
      }
      else
      {
        incl[i] = 1;
        double dd = n[reac[i]];
        if (dd <= l[k])
          dd = (u[k] - l[k])/1000.;
        if (n_ == R)
          eta[k] = dd;
        else
        {
          eta1[k] = dd;
          vptr[k] = &eta[i];
        }
        if (debug)
          *debug << "eta[" << i << "] is equal to " << dd
            << " l is " << l[k] << " and h is " << u[k] << endl;
        ++k;
      }
    }
    if (debug)
    {
      *debug << "number of unknowns is " << n_ << endl 
        << "Starting point is" << endl;
      fun(n_, Tp, f);
    }
    if (!n_)
    {
      SetMoleFractions();
      if (debug)
        *debug << endl;
      return vint;
    }
    double *xx;
    if (n_ == R)
      xx = &*eta.begin();
    else
      xx = &*eta1.begin();
    e.resize(n_);
    s_copy(task, "START", 60L, 5L);
    while (1)
    {
      setulb_(&n_, &m, xx, &*l.begin(), &*u.begin(), &*nbd.begin(), &f, 
          &*e.begin(), &factr, &pgtol, &*wa.begin(), &*iwa.begin(), task, 
					&iprint, csave, lsave, isave, dsave);
      if (s_cmp(task, "FG", 2L, 2L) == 0)
      {
        fun(n_, Tp, f);
      }
      else if (s_cmp(task, "NEW_X", 5L, 5L) == 0)
        continue;
      else
        break;
    }
    if (debug)
      *debug << endl;
    finish = true;
    for (size_t i = 0; i < n_; ++i)
      if (fabs(e[i]) > 10.)
      {
        finish = false;
        ++change_basis;
        if (debug)
          *debug << "no convergence - a new try with a new basis" << endl;
        break;
      }
  }
  if (!finish)
    for (size_t i = 0; i < e.size(); ++i)
      if (fabs(e[i]) > 10.)
        cout << "AssSolManyReact: no convergence " << i << "-th gradient is " 
          << e[i] << endl;
  return vint;
}

#ifdef INCLUDE_VCS
void associated_solution::read(const SGML &e_)
{
  clear();
  SGML e = e_;
  ReadComponents(vf, e, ref);
  if (e.FindString("mole_fractions") == "equilibrium")
    cout << "associated_solution: mole_fractions=equilibrium is not supported"
      << endl;
  parser p(e);
  SGML el;
  p.GetSGML(el);
  el.compare("internal_solution");
  {
    parser p2(el);
    p2.GetSGML(e);
    sol.read(e);
  }
  solver.set(vf, sol);
  xsol.resize(sol.size());
  for (size_t i = 0; i < sol.size(); ++i)
    xsol.SetX(i, const_cast<double&>(solver.xi()[i]));
  string FileName = e_.FindString("DebugFile");
  bool debug = e_.FindInt("debug");
  if (phase::debug || debug)
  {
    string s = FileName;
    if (s.empty())
    {
      string id = e_.FindString("id");
      if (id.empty())
      {
        const vec_formula &vf = comps();
        for (size_t i = 0; i < size(); ++i)
        {
          id += vf[i].cf();
          if (i != size() - 1)
            id += "_";
        }
      }
      s = string("debug.phase.") + id;
    }
    SetDebug(s);
  }
  vint.resize(sol.size() + 1);
  vsint.resize(sol.size() + 1);
  vsint[0] = "ntot";
  for (size_t i = 0; i < sol.size(); ++i)
    vsint[i + 1] = string("x(") + sol.comps()[i].cf() + ")";
}

void associated_solution::WriteBody(ostream &out, size_t shift) const
{
  WriteComponents(out, shift);
  ref.write(out, &vf, shift);
  out << PutTab(shift) << "<internal_solution>" << endl;
  sol.write(out, shift + 1);
  out << PutTab(shift) << "</internal_solution>" << endl;
}

double associated_solution::Z(function f, index i, const StateTp &Tp,
                              const StateX &x) const
{
  IntVarsEq(Tp, x);
  double sum = 0.;
  if (i & ::ref)
    sum += ref.Z(f, Tp, x);
  double f1, f2, stepT;
  if ((i & ::mix) == ::mix)
    switch (f)
    {
      case ::G:
        sum += solver.Ntot()*sol.Z(f, ::full, Tp, xsol);
        break;
      case ::H:
        sum += solver.Ntot()*sol.Z(f, ::full, Tp, xsol) - Tp.T()*GnT(Tp);
        break;
      case ::S:
        sum += solver.Ntot()*sol.Z(f, ::full, Tp, xsol) - GnT(Tp);
        break;
      case ::V:
        sum += solver.Ntot()*sol.Z(f, ::full, Tp, xsol) + Gnp(Tp);
        break;
      case ::Cp:
        stepT = 1e-3*(Tp.T() + 1e-3);
        f1 = Z(::H, ::mix, StateTp(Tp.T() - stepT, Tp.p()), x);
        f2 = Z(::H, ::mix, StateTp(Tp.T() + stepT, Tp.p()), x);
        sum += (f2 - f1)/stepT/2.;
        break;
      default:
        throw gError("associated_solution: unknown code for the function");
    }
  return sum;
}
#endif
