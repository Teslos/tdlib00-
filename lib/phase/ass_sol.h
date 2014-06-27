/*
Copyright (C) 1998, 1999, 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                       http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#ifndef __ASS_SOL_H
#define __ASS_SOL_H

#include <iterator>
#include <intvars.h>
#include <simple.h>
#include <float.h>
#include <minit.h>

#ifdef INCLUDE_VCS
#include <../vcs/vcs.h>
#endif

class AssociatedSolution : public RefPhase
{
public:
  AssociatedSolution() {}
  
  virtual AssociatedSolution* clone() const 
    {return new AssociatedSolution(*this);}

  virtual void read(const SGML &e);
  virtual void WriteBody(ostream &out, size_t shift = 0) const
  {
    throw gError("AssociatedSolution: try to write the proxy");
  }
  virtual const vec_formula& comps() const
  {
    throw gError("AssociatedSolution: try to query the proxy");
  }
  virtual double Z(function f, index i, const StateTp &Tp, 
      const StateX &x) const
  {
    throw gError("AssociatedSolution: try to compute the proxy");
  }
};

class AssociatedSolutionBasis
{
  friend class AssociatedSolution;
protected:

  SimpleSolution ref; // C
  SimpleSolution internal; // N
  index equil; //full or ref_ideal
  size_t C;
  size_t N;
  size_t R; // N - C
  
  matrix A; // C*N

  mutable matrix Ac; //C*C
  mutable vec_int basis; // C  -> Ac
  mutable matrix Z_; // C*R -> Ac Z = Ar
  mutable matrix Nu; // N*R (comprises Z and I)
  mutable vec_double DelNu; // R
  mutable vec_int reac; // R   -> Ar (ordinal numbers in N)
  mutable matrix AcInv; // no = AcInv x
  mutable vec_double SumAcInv; // SumAcInv(i) = Sum j AcInv(j, i)
  mutable vec_int ipiv; // C

  mutable vec_int major; // reac = major + minor (ordinal number in R)
  mutable vec_int minor; // 

  mutable matrix mnc; // N*C
  mutable matrix mrc; // R*C
  mutable matrix mnn; // N*N

  mutable vec_double no; // C Ac no = x
  mutable vec_double n; // N moles of internal
  mutable vec_double eta; // R chemical variables
// this is saved to vint  
  mutable StateX vy; // N mole fractions of internal
  mutable double *ntot;
  
  ofstream *dbg;

  void SetNo(const StateX &x) const
  {
    for (size_t i = 0; i < C; ++i)
    {
      no[i] = 0.;
      for (size_t j = 0; j < C; ++j)
        no[i] += AcInv(i, j)*x[j];
      if (no[i] < 0. && no[i] > -DBL_EPSILON*10000.)
        no[i] = 0.;
    }
    if (dbg)
    {
      *dbg << "no are ";
      copy(no.begin(), no.end(), ostream_iterator<double>(*dbg, " "));
      *dbg << endl;
    }
  }
  void SetMoleFractions() const
  {
    *ntot = 0.;
    for (size_t i = 0; i < C; ++i)
    {
      n[basis[i]] = no[i];
      for (size_t j = 0; j < R; ++j)
        n[basis[i]] -= Z_(i, j)*eta[j];
      *ntot += n[basis[i]];
    }
    for (size_t i = 0; i < R; ++i)
    {
      n[reac[i]] = eta[i];
      *ntot += n[reac[i]];
    }
    for (size_t i = 0; i < N; ++i)
      vy[i] = n[i]/(*ntot);
  }
  const matrix& dmudx(const SymMatrix &sm) const
  {
// mnn = dmu_k/dy_l   
    for (size_t l = 0; l < N; ++l)
    {
      double Sum = 0.; // +z[l]-z[l]
      for (size_t m = 0; m < N; ++m)
        Sum -= sm(m, l)*vy[m];
      for (size_t k = 0; k < N; ++k)
        mnn(k, l) = Sum + sm(k, l);
    }
// mnc = Sum_l dmu_k/dy_l dy_l/dx_j   
    for (size_t k = 0; k < N; ++k)
      for (size_t j = 0; j < C; ++j)
      {
        mnc(k, j) = 0;
        size_t ib = 0;
        for (size_t l = 0; l < N; ++l)
        {
          double der = 0.;
          if (ib < C && l == basis[ib])
          {
            der = AcInv(ib, j);
            ++ib;
          }
          mnc(k, j) += mnn(k, l)*(der - vy[l]*SumAcInv[j]);
        }
        mnc(k, j) /= (*ntot);
      }
    return mnc;
  }

  double Zex(function f, const StateTp &Tp, const StateX &x) const
  {
    return ref.Z(f, ::full, Tp, x);
  }
//add to res
  void dZexdx(function f, const StateTp &Tp, const StateX &x, 
      vec_double& res) const
  {
    const vec_double &v = ref.dZdx(f, ::full, Tp, x);
    for (size_t i = 0; i < C; ++i)
      res[i] += v[i];
  }
  
//the next two functions take const vec_double &y from vint 
//do not forget to set vint (for example by calling yeq)
  double Zin(function f, const StateTp &Tp, const StateX &x) const
  {
    return (*ntot)*internal.Z(f, ::full, Tp, vy);
  }
  void dZindx(function f, const StateTp &Tp, const StateX &x, 
      vec_double &res) const
  {
    const vec_double &z = internal.z(f, ::full, Tp, vy);
    for (size_t i = 0; i < C; ++i)
      for (size_t j = 0; j < C; ++j)
        res[i] += z[basis[j]]*AcInv(j, i);
  }
  
public: 
  static double tol;
  
  AssociatedSolutionBasis() : equil(::ref_ideal), dbg(0) {}
  
  void clear()
  {
    A.resize(C, N);
    Ac.resize(C, C);
    basis.resize(C);
    Z_.resize(C, R);
    Nu.resize(N, R);
    DelNu.resize(R);
    reac.resize(R);
    major.clear();
    minor.clear();
    no.resize(C);
    AcInv.resize(C, C);
    SumAcInv.resize(C);
    ipiv.resize(C);
    mnc.resize(N, C);
    mrc.resize(R, C);
    mnn.resize(N, N);
    vy.resize(N);
    n.resize(N);
    eta.resize(R);
  }
  
  void read(const SGML &e);
  void WriteAttributes(ostream &out, size_t shift = 0) const
  {
    out << endl << PutTab(shift) << "equilibrated=";
    if (equil == ::full)
      out << 1;
    else
      out << 0;
  }
  void WriteBody(ostream &out, size_t shift = 0) const
  {
    ref.write(out, shift);
    internal.write(out, shift);
  }
  
  void initialize();
// set Ac and Z 
  size_t SetOptimalBasis() const;
// set DelNu Nu AcInv SumAcInv
  void SetStoichiometry() const;
};

class AssSolOneReact : public OneInternalVariable, 
  public AssociatedSolutionBasis
{
  virtual double Zex(function f, const StateTp &Tp, const StateX &x) const
  {
    return AssociatedSolutionBasis::Zex(f, Tp, x);
  }

//add to res
  virtual void dZexdx(function f, const StateTp &Tp, const StateX &x, 
      vec_double &res) const
  {
    AssociatedSolutionBasis::dZexdx(f, Tp, x, res);
  }

//the next functions take const vec_double &y from vint 
//do not forget to set vint (for example by calling yeq)
  virtual double Zin(function f, const StateTp &Tp, const StateX &x) const
  {
    return AssociatedSolutionBasis::Zin(f, Tp, x);
  }

//add to res
  virtual void dZindx(function f, const StateTp &Tp, const StateX &x, 
      vec_double &res) const
  {
    AssociatedSolutionBasis::dZindx(f, Tp, x, res);
  }
  
  mutable const StateTp *Tp;

  double f(const double &y) const;
  
//result in vint
  virtual const vec_double& yeq(const StateTp &Tp, const StateX &x) const;
  
  virtual double dEqdy(const StateTp &Tp, const StateX &x) const
  {
    const SymMatrix &sm = internal.d2Zdx2(::G, equil, Tp, vy);
    double Sum = 0.;
    for (size_t k = 0; k < N; ++k)
    {
      for (size_t l = 0; l < k; ++l)
        Sum += 2.*sm(k, l)*(Nu(k, 0) - vy[k]*DelNu[0])
          *(Nu(l, 0) - vy[l]*DelNu[0]);
      Sum += sm(k, k)*pow((Nu(k, 0) - vy[k]*DelNu[0]), 2);
    }
    return Sum/(*ntot);
  }
  virtual double dEqdT(const StateTp &Tp, const StateX &x) const
  {
    const vec_double &z = internal.z(::S, equil, Tp, vy);
    double Sum = 0;
    for (size_t i = 0; i < C; ++i)
      Sum -= z[basis[i]]*Z_(i, 0);
    Sum += z[reac[0]];
    return -Sum;
  }
  virtual double dEqdp(const StateTp &Tp, const StateX &x) const
  {
    const vec_double &z = internal.z(::V, equil, Tp, vy);
    double Sum = 0;
    for (size_t i = 0; i < C; ++i)
      Sum -= z[basis[i]]*Z_(i, 0);
    Sum += z[reac[0]];
    return Sum;
  }
//rewrite res
  virtual void dEqdx(const StateTp &Tp, const StateX &x, vec_double &res) const
  {
    const SymMatrix &sm = internal.d2Zdx2(::G, equil, Tp, vy);
    const matrix &mat = dmudx(sm);
    for (size_t j = 0; j < C; ++j)
    {
      res[j] = 0.;
      for (size_t k = 0; k < C; ++k)
        res[j] -= mat(basis[k], j)*Z_(k, 0);
      res[j] += mat(reac[0], j);
    }
  }
  virtual double dZindy(function f, const StateTp &Tp, const StateX &x) const
  {
    double Sum = 0.;
    if (derZin)
    {
      const vec_double &z = internal.z(f, ::full, Tp, vy);
      for (size_t i = 0; i < C; ++i)
      Sum -= z[basis[i]]*Z_(i, 0);
      Sum += z[reac[0]];
    }
    return Sum;
  }

  mutable bool derZin;
  mutable double b;
  mutable double a;
  virtual void CopyInternalVariables(const vec_double &y) const
  {
    if (vint.size() != y.size())
      throw gError("InternalVariables: dimension of y != vint");
    vint = y;
    if ((*ntot)*vy[reac[0]] < a + AssociatedSolutionBasis::tol
        || (*ntot)*vy[reac[0]] > b - AssociatedSolutionBasis::tol)
      derZin = false;
    else
      derZin = true;
  }
public:
  AssSolOneReact() 
  {
    C = N = R = 0;
    clear();
  }
  AssSolOneReact(const AssociatedSolutionBasis &x, const SGML &e) 
    : AssociatedSolutionBasis(x)
  {
    SetDebug(e);
    clear();
  }

  virtual void clear(size_t n_ = 0, size_t m = 0, size_t l = 0)
  {
    if (R > 1)
      throw gError("AssSolOneReact: R > 1");
    dbg = debug;
    OneInternalVariable::clear(C, N + 1);
    ntot = &vint[0];
    vsint[0] = "n";
    for (size_t i = 0; i < N; ++i)
    {
      vy.SetX(i, vint[i + 1]);
      vsint[i + 1] = string("x(") + internal.comps()[i].cf() + ")";
    }
    if (equil == ::full)
      SetIfEqTheSame(true);
    else
      SetIfEqTheSame(false);
    fill(n.begin(), n.end(), 1.);
    SetOptimalBasis();
    SetStoichiometry();
  }
  
  virtual AssSolOneReact* clone() const {return new AssSolOneReact(*this);}
  virtual void read(const SGML &e)
  {
    AssociatedSolutionBasis::read(e);
    SetDebug(e);
    clear();
  }
  virtual void WriteAttributes(ostream &out, size_t shift = 0) const
  {
    AssociatedSolutionBasis::WriteAttributes(out, shift);
    RefInternalVariables::WriteAttributes(out, shift);
  }
  virtual void WriteBody(ostream &out, size_t shift = 0) const
  {
    AssociatedSolutionBasis::WriteBody(out, shift);
  }
  virtual const vec_formula& comps() const
  {
    return ref.comps();
  }
};

class AssSolManyReact : public InternalVariables, 
  public AssociatedSolutionBasis
{
  virtual double Zex(function f, const StateTp &Tp, const StateX &x) const
  {
    return AssociatedSolutionBasis::Zex(f, Tp, x);
  }

//add to res
  virtual void dZexdx(function f, const StateTp &Tp, const StateX &x, 
      vec_double &res) const
  {
    AssociatedSolutionBasis::dZexdx(f, Tp, x, res);
  }

  void LinearProgramming(const StateTp &Tp, const StateX &x) const;
  void fun(long n, const StateTp &Tp, double &f) const;
//result in vint
  virtual const vec_double& yeq(const StateTp &Tp, const StateX &x) const;
  
//the next functions take const vec_double &y from vint 
//do not forget to set vint (for example by calling yeq)
  virtual double Zin(function f, const StateTp &Tp, const StateX &x) const
  {
    return AssociatedSolutionBasis::Zin(f, Tp, x);
  }

//add to res
  virtual void dZindx(function f, const StateTp &Tp, const StateX &x, 
      vec_double &res) const
  {
    AssociatedSolutionBasis::dZindx(f, Tp, x, res);
  }
  
//result in dedy
  virtual matrix& dEqdy(const StateTp &Tp, const StateX &x) const
  {
    const SymMatrix &sm = internal.d2Zdx2(::G, equil, Tp, vy);
    dedy.resize(major.size(), major.size());
    for (size_t i = 0; i < major.size(); ++i)
      for (size_t j = 0; j <= i; ++j)
      {
        dedy(i, j) = 0.;
        for (size_t k = 0; k < N; ++k)
        {
          for (size_t l = 0; l < k; ++l)
            dedy(i, j) += 2.*sm(k, l)*(Nu(k, major[i]) - vy[k]*DelNu[major[i]])
              *(Nu(l, major[j]) - vy[l]*DelNu[major[j]]);
          dedy(i, j) += sm(k, k)*(Nu(k, major[i]) - vy[k]*DelNu[major[i]])
            *(Nu(k, major[j]) - vy[k]*DelNu[major[j]]);
        }
        dedy(j, i) = dedy(i, j) = dedy(i, j)/(*ntot);
      }
    return dedy;
  }
//result in res
  virtual vec_double& dEqdT(const StateTp &Tp, const StateX &x) const
  {
    res.resize(major.size(), 1);
    const vec_double &z = internal.z(::S, equil, Tp, vy);
    for (size_t i = 0; i < major.size(); ++i)
    {
      res[i] = 0.;
      for (size_t k = 0; k < C; ++k)
        res[i] -= z[basis[k]]*Z_(k, major[i]);
      res[i] += z[reac[major[i]]];
      res[i] = -res[i];
    }
    return res;
  }
//result in res
  virtual vec_double& dEqdp(const StateTp &Tp, const StateX &x) const
  {
    res.resize(major.size(), 1);
    const vec_double &z = internal.z(::V, equil, Tp, vy);
    for (size_t i = 0; i < major.size(); ++i)
    {
      res[i] = 0.;
      for (size_t k = 0; k < C; ++k)
        res[i] -= z[basis[k]]*Z_(k, major[i]);
      res[i] += z[reac[major[i]]];
    }
    return res;
  }
//result in res
  virtual matrix& dEqdx(const StateTp &Tp, const StateX &x) const
  {
    res.resize(major.size(), C);
    const SymMatrix &sm = internal.d2Zdx2(::G, equil, Tp, vy);
    const matrix &mat = dmudx(sm);
    for (size_t i = 0; i < major.size(); ++i)
    {
      for (size_t j = 0; j < C; ++j)
      {
        res(i, j) = 0.;
        for (size_t k = 0; k < C; ++k)
          res(i, j) -= mat(basis[k], j)*Z_(k, major[i]);
        res(i, j) += mat(reac[major[i]], j);
      }
    }
    return res;
  }
  
//result in e
  virtual const vec_double& dZindy(function f, const StateTp &Tp, 
      const StateX &x) const 
  {
    const vec_double &z = internal.z(f, ::full, Tp, vy);
    e.resize(major.size());
    for (size_t i = 0; i < major.size(); ++i)
    {
      e[i] = 0.;
      for (size_t k = 0; k < C; ++k)
        e[i] -= z[basis[k]]*Z_(k, major[i]);
      e[i] += z[reac[major[i]]];
    }
    return e;
  }
  virtual void CopyInternalVariables(const vec_double &y) const
  {
    if (vint.size() != y.size())
      throw gError("InternalVariables: dimension of y != vint");
    vint = y;
    major.clear();
    minor.clear();
    for (size_t i = 0; i < R; ++i)
    {
// it is necessary to use values only from y      
      if ((*ntot)*vy[reac[i]] > AssociatedSolutionBasis::tol)
        major.push_back(i);
      else
        minor.push_back(i);
    }
  }
  mutable vec_double eta1;
  mutable vec_double_ptr vptr;
  mutable vec_int incl;
  mutable minit lp;
  
public:
  static int lbfg_m;
  static double factr;
  static double pgtol;

  AssSolManyReact() 
  {
    C = N = R = 0;
    clear();
  } 
  AssSolManyReact(const AssociatedSolutionBasis &x, const SGML &e) 
    : AssociatedSolutionBasis(x)
  {
    SetDebug(e);
    clear();
  }

  virtual void clear(size_t n = 0, size_t m = 0, size_t l = 0)
  {
    lp.debug = dbg = debug;
    InternalVariables::clear(C, N + 1, R);
    ntot = &vint[0];
    vsint[0] = "ntot";
    for (size_t i = 0; i < N; ++i)
    {
      vy.SetX(i, vint[i + 1]);
      vsint[i + 1] = string("x(") + internal.comps()[i].cf() + ")";
    }
    if (equil == ::full)
      SetIfEqTheSame(true);
    else
      SetIfEqTheSame(false);
    eta1.resize(R);
    vptr.resize(R);
    incl.resize(R);
  }
  
  virtual AssSolManyReact* clone() const {return new AssSolManyReact(*this);}
  virtual void read(const SGML &e)
  {
    AssociatedSolutionBasis::read(e);
    SetDebug(e);
    clear();
  }
  virtual void WriteAttributes(ostream &out, size_t shift = 0) const
  {
    AssociatedSolutionBasis::WriteAttributes(out, shift);
    RefInternalVariables::WriteAttributes(out, shift);
  }
  virtual void WriteBody(ostream &out, size_t shift = 0) const
  {
    AssociatedSolutionBasis::WriteBody(out, shift);
  }

  virtual const vec_formula& comps() const
  {
    return ref.comps();
  }
};

#ifdef INCLUDE_VCS
class associated_solution : public RefSimple
{
  SimpleSolution sol;
  FuncTpx ref;
  mutable VCS solver;
  mutable StateX xsol; // set on solver.x()
  mutable vec_double vint;
  vec_string vsint;
  string FileName;
  bool debug;

  double GnT(const StateTp &Tp) const
  {
    const double* b = solver.dT();
    const vec_double &mu = sol.mu(Tp, xsol);
    double Sum = 0.;
    for (size_t i = 0; i < sol.size(); ++i)
      Sum += b[i]*mu[i];
    return Sum;
  }
  double Gnp(const StateTp &Tp) const
  {
    const double* b = solver.dp();
    const vec_double &mu = sol.mu(Tp, xsol);
    double Sum = 0.;
    for (size_t i = 0; i < sol.size(); ++i)
      Sum += b[i]*mu[i];
    return Sum;
  }
  const vec_double& Gnb(index ind, const StateTp &Tp) const
  {
    const matrix& b = solver.db();
    const vec_double &mu = sol.mu(Tp, xsol);
    fill(vz.begin(), vz.end(), 0.);
    if (ind & ::ref)
    {
      dynamic_cast<const Reference*>(ref.pointer())->z(::G, Tp, vz);
    }
    if ((ind & ::mix) == ::mix)
    {
      size_t ne = size();
      for (size_t i = 0; i < ne; ++i)
      {
        size_t m = sol.comps().size();
        for (size_t j = 0; j < m ; ++j)
          vz[i] += mu[j]*b(j, i);
      }
    }
    return vz;
  }
public:
  associated_solution() {}
  virtual ~associated_solution() {}

  virtual void clear()
  {
    sol.clear(); 
    xsol.resize(0); 
    solver.clear(); 
    RefSimple::clear();
  }
  virtual associated_solution* clone() const
    {return new associated_solution(*this);}

  virtual void read(const SGML &e);
  virtual void WriteAttributes(ostream &out, size_t shift = 0) const
  {
    out << endl << PutTab(shift) << "debug=" << debug;
    if (!FileName.empty())
      out << endl << PutTab(shift) << "DebugFile=" << FileName;
  }
  virtual void WriteBody(ostream &out, size_t shift = 0) const;

  void SetDebug(const string &str)
  {
    solver.set_debug_file(str);
  }
  
  virtual size_t NIntVars() const {return vint.size();}
  virtual const vec_string& NamesIntVars() const {return vsint;}
  virtual const vec_double& IntVarsEq(const StateTp &Tp, const StateX &x) const 
  {
    if (!solver.solve(Tp, x))
      cout << "association solution - not solved" << endl;
    vint[0] = solver.Ntot();
    for (size_t i = 0; i < sol.size(); ++i)
      vint[i + 1] = solver.xi()[i];
    return vint;
  }

  virtual double Z(function f, index i, const StateTp &Tp, 
      const StateX &x) const;
  virtual const vec_double& dZdx(function f, index i, const StateTp &Tp, 
      const StateX &x) const
  {
    throw gError("associated_solution: dZdx is not implemented");
  }
  virtual const vec_double& z(function f, index i, const StateTp &Tp, 
      const StateX &x) const
  {
    if (f != ::G)
      throw gError("associated_solution: function not supported");
    IntVarsEq(Tp, x);
    return Gnb(i, Tp);
  }
  virtual const vec_double& z_y(function f, index i, const StateTp &Tp, 
      const StateX &x, const vec_double& y) const
    {return z(f, i, Tp, x);}
};
#endif

#endif
