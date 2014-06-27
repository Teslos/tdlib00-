/*
Copyright (C) 1998, 1999, 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                       http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#include "intvars.h"
#include <lapack.h>
#include <iomanip>
#include <iterator>
#include <float.h>

double RefInternalVariables::step = 10.*pow(DBL_EPSILON, 1./3.);

void RefInternalVariables::DeleteDebug()
{
  if (debug)
  {
    delete debug;
    debug = 0;
  }
}

void RefInternalVariables::SetDebug(const SGML &e)
{
  FileName = e.FindString("DebugFile");
  debug_ = e.FindInt("debug");
  if (phase::debug || debug_)
  {
    DeleteDebug();
    string s = FileName;
    if (s.empty())
    {
      string id = e.FindString("id");
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
    debug = new ofstream(s.c_str());
  }
}

double OneInternalVariable::Z_y(function f, index i, const StateTp &Tp, 
    const StateX &x, const vec_double &y) const
{
  double sum = 0.;
  if (debug)
    *debug << "OneInternalVariable::Z_y function is " << f << " index is " 
      << i << endl;
  if (i & ::ref)
  {
    double z = Zex(f, Tp, x);
    if (debug)
      *debug << "Zex is " << z << endl;
    sum += z;
  }
  if ((i & ::mix) == ::mix)
  {
    CopyInternalVariables(y);
    if (debug)
    {
      *debug << "Internal Variables" << endl;
      copy(y.begin(), y.end(), ostream_iterator<double>(*debug, " "));  
      *debug << endl;
    }
    if (same)
    {
      double z = Zin(f, Tp, x);
      if (debug)
        *debug << "Zin is " << z << endl;
      sum += z;
      switch (f)
      {
        case ::G:
        case ::H:
        case ::S:
        case ::V:
          break;
        case ::Cp:
        {
          double z1 = dZindy(::H, Tp, x);
          double z2 = dydT(Tp, x);
          if (debug)
            *debug << "dZindy(H) is " << z1 << " dydT is " << z2 << endl;
          sum += z1*z2;
          break;
        }
        case ::dVdT:
        {
          double z1 = dZindy(::V, Tp, x);
          double z2 = dydT(Tp, x);
          if (debug)
            *debug << "dZindy(V) is " << z1 << " dydT is " << z2 << endl;
          sum += z1*z2;
          break;
        }
        case ::dVdp:
        {
          double z1 = dZindy(::V, Tp, x);
          double z2 = dydp(Tp, x);
          if (debug)
            *debug << "dZindy(V) is " << z1 << " dydp is " << z2 << endl;
          sum += z1*z2;
          break;
        }
        default:
          throw gError("OneInternalVariable: function is not defined");
      }
    }
    else
      switch (f)
      {
        case ::G:
        case ::H:
        case ::S:
        case ::V:
        {
          double z = Zin(f, Tp, x);
          if (debug)
            *debug << "Zin is " << z << endl;
          sum += z;
          switch (f)
          {
            case ::G:
              break;
            case ::H:
            {
              double z1 = dZindy(::G, Tp, x);
              double z2 = dydT(Tp, x);
              if (debug)
                *debug << "dZindy(G) is " << z1 << " dydT is " << z2 << endl;
              sum += -Tp.T()*z1*z2;
              break;
            }
            case ::S:
            {
              double z1 = dZindy(::G, Tp, x);
              double z2 = dydT(Tp, x);
              if (debug)
                *debug << "dZindy(G) is " << z1 << " dydT is " << z2 << endl;
              sum += -z1*z2;
              break;
            }
            case ::V:
            {
              double z1 = dZindy(::G, Tp, x);
              double z2 = dydp(Tp, x);
              if (debug)
                *debug << "dZindy(G) is " << z1 << " dydp is " << z2 << endl;
              sum += z1*z2;
              break;
            }
          }
          break;
        }
        case ::Cp:
        {
          Z_T num(::H, *this, ::mix, Tp.p(), x);
          double z = FirstDerivative(num, Tp.T(), step);
          if (debug)
            *debug << "FirstDerivative is " << z << endl;
          sum += z;
          break;
        }
        case ::dVdT:
        {
          Z_T num(::V, *this, ::mix, Tp.p(), x);
          double z = FirstDerivative(num, Tp.T(), step);
          if (debug)
            *debug << "FirstDerivative is " << z << endl;
          sum += z;
          break;
        }
        case ::dVdp:
        {
          Z_p num(::V, *this, ::mix, Tp.T(), x);
          double z = FirstDerivative(num, Tp.p(), step);
          if (debug)
            *debug << "FirstDerivative is " << z << endl;
          sum += z;
          break;
        }
        default:
          throw gError("OneInternalVariable: function is not defined");
      }
  }
  if (debug)
    *debug << endl;
  return sum;
}

const vec_double& OneInternalVariable::dZdx_y(function f, index i, 
    const StateTp &Tp, const StateX &x, const vec_double &y) const
{
  if (debug)
    *debug << "OneInternalVariable::dZdx_y function is " << f << " index is " 
      << i << endl;
  if ((i & ::mix) == ::mix)
  {
    CopyInternalVariables(y);
    if (debug)
    {
      *debug << "Internal Variables" << endl;
      copy(y.begin(), y.end(), ostream_iterator<double>(*debug, " "));  
      *debug << endl;
    }
    if (same)
      switch (f)
      {
        case ::G:
          fill(vz.begin(), vz.end(), 0.);
          break;
        case ::H:
        case ::S:
        case ::V:
        {
          dydx(Tp, x, vz);
          if (debug)
          {
            *debug << "dydx is ";
            copy(vz.begin(), vz.end(), ostream_iterator<double>(*debug, " "));  
            *debug << endl;
          }
          double d = dZindy(f, Tp, x);
          if (debug)
            *debug << "dZindy(" << f << ") is " << d << endl;             
          for (size_t i = 0; i < size(); ++i)
            vz[i] *= d;
          break;
        }
        default:
          throw gError("OneInternalVariable: function is not supported");
      }
    else
      switch (f)
      {
        case ::G:
        {
          dydx(Tp, x, vz);
          if (debug)
          {
            *debug << "dydx is ";
            copy(vz.begin(), vz.end(), ostream_iterator<double>(*debug, " "));  
            *debug << endl;
          }
          double d = dZindy(f, Tp, x);
          if (debug)
            *debug << "dZindy(" << f << ") is " << d << endl;             
          for (size_t i = 0; i < size(); ++i)
            vz[i] *= d;
          break;
        }
        case ::H:
        {
          ZA_x num(::H, *this, ::mix, Tp);
          ArrayFirstDerivative(num, x, vz, step);
          if (debug)
          {
            *debug << "ArrayFirstDer is ";
            copy(vz.begin(), vz.end(), ostream_iterator<double>(*debug, " "));  
            *debug << endl;
          }
          break;
        }
        case ::S:
        {
          ZA_x num(::S, *this, ::mix, Tp);
          ArrayFirstDerivative(num, x, vz, step);
          if (debug)
          {
            *debug << "ArrayFirstDer is ";
            copy(vz.begin(), vz.end(), ostream_iterator<double>(*debug, " "));  
            *debug << endl;
          }
          break;
        }
        case ::V:
        {
          ZA_x num(::V, *this, ::mix, Tp);
          ArrayFirstDerivative(num, x, vz, step);
          if (debug)
          {
            *debug << "ArrayFirstDer is ";
            copy(vz.begin(), vz.end(), ostream_iterator<double>(*debug, " "));  
            *debug << endl;
          }
          break;
        }
        default:
          throw gError("OneInternalVariable: function is not supported");
      }
    if (same || f == ::G)
    {
      dZindx(f, Tp, x, vz);
      if (debug)
      {
        *debug << "after dZindx the vector is ";
        copy(vz.begin(), vz.end(), ostream_iterator<double>(*debug, " "));  
        *debug << endl;
      }
    }
  }
  if (i & ::ref)
  {
    dZexdx(f, Tp, x, vz);
    if (debug)
    {
      *debug << "after dZexdx the vector is ";
      copy(vz.begin(), vz.end(), ostream_iterator<double>(*debug, " "));  
      *debug << endl;
    }
  }
  if (debug)
    *debug << endl;
  return vz;
}

//result in res
const vec_double& InternalVariables::dydT(const StateTp &Tp, 
    const StateX &x) const
{
  matrix &m1 = dEqdy(Tp, x);
  vec_double &m2 = dEqdT(Tp, x);
  if (debug)
  {
    *debug << "dedy and res before finding dydT" << endl;
    copy(m1.begin(), m1.end(), ostream_iterator<double>(*debug, " "));  
    *debug << endl;
    copy(m2.begin(), m2.end(), ostream_iterator<double>(*debug, " "));  
    *debug << endl;
  }
  for (size_t i = 0; i < m2.size(); ++i)
    m2[i] = -m2[i];
  dgesv(m1.NRows(), 1, &*m1.begin(), &*ipiv.begin(), &*m2.begin());
  return m2;
}

//result in res
const vec_double& InternalVariables::dydp(const StateTp &Tp, 
    const StateX &x) const
{
  matrix &m1 = dEqdy(Tp, x);
  vec_double &m2 = dEqdp(Tp, x);
  if (debug)
  {
    *debug << "dedy and res before finding dydp" << endl;
    copy(m1.begin(), m1.end(), ostream_iterator<double>(*debug, " "));  
    *debug << endl;
    copy(m2.begin(), m2.end(), ostream_iterator<double>(*debug, " "));  
    *debug << endl;
  }
  for (size_t i = 0; i < m2.size(); ++i)
    m2[i] = -m2[i];
  dgesv(m1.NRows(), 1, &*m1.begin(), &*ipiv.begin(), &*m2.begin());
  return m2;
}

//result in res
const matrix& InternalVariables::dydx(const StateTp &Tp, const StateX &x) const
{
  matrix &m1 = dEqdy(Tp, x);
  matrix &m2 = dEqdx(Tp, x);
  if (debug)
  {
    *debug << "dedy and res before finding dydx" << endl;
    copy(m1.begin(), m1.end(), ostream_iterator<double>(*debug, " "));  
    *debug << endl;
    copy(m2.begin(), m2.end(), ostream_iterator<double>(*debug, " "));  
    *debug << endl;
  }
  for (size_t i = 0; i < m2.size(); ++i)
    m2[i] = -m2[i];
  dgesv(m1.NRows(), m2.NCols(), &*m1.begin(), &*ipiv.begin(), &*m2.begin());
  return m2;
}

double InternalVariables::Z_y(function f, index i, const StateTp &Tp, 
    const StateX &x, const vec_double &y) const
{
  double sum = 0.;
  if (debug)
    *debug << "InternalVariables::Z_y function is " << f << " index is " 
      << i << endl;
  if (i & ::ref)
  {
    double z = Zex(f, Tp, x);
    if (debug)
      *debug << "Zex is " << z << endl;
    sum += z;
  }
  if ((i & ::mix) == ::mix)
  {
    CopyInternalVariables(y);
    if (debug)
    {
      *debug << "Internal Variables" << endl;
      copy(y.begin(), y.end(), ostream_iterator<double>(*debug, " "));  
      *debug << endl;
    }
    if (same)
    {
      double z = Zin(f, Tp, x);
      if (debug)
        *debug << "Zin is " << z << endl;
      sum += z;
      switch (f)
      {
        case ::G:
        case ::H:
        case ::S:
        case ::V:
          break;
        case ::Cp:
        {
          const vec_double &r1 = dZindy(::H, Tp, x);
          const vec_double &r2 = dydT(Tp, x);
          if (debug)
          {
            *debug << "dZindy(H) is ";
            copy(r1.begin(), r1.end(), ostream_iterator<double>(*debug, " "));  
            *debug << endl;
            *debug << "dydT is ";
            copy(r2.begin(), r2.begin() + r1.size(), 
                ostream_iterator<double>(*debug, " ")); 
            *debug << endl;
          }
          for (size_t i = 0; i < r1.size(); ++i)
            sum += r1[i]*r2[i];
          break;
        }
        case ::dVdT:
        {
          const vec_double &r1 = dZindy(::V, Tp, x);
          const vec_double &r2 = dydT(Tp, x);
          if (debug)
          {
            *debug << "dZindy(V) is ";
            copy(r1.begin(), r1.end(), ostream_iterator<double>(*debug, " "));  
            *debug << endl;
            *debug << "dydT is ";
            copy(r2.begin(), r2.begin() + r1.size(), 
                ostream_iterator<double>(*debug, " ")); 
            *debug << endl;
          }
          for (size_t i = 0; i < r1.size(); ++i)
            sum += r1[i]*r2[i];
          break;
        }
        case ::dVdp:
        {
          const vec_double &r1 = dZindy(::V, Tp, x);
          const vec_double &r2 = dydp(Tp, x);
          if (debug)
          {
            *debug << "dZindy(V) is ";
            copy(r1.begin(), r1.end(), ostream_iterator<double>(*debug, " "));  
            *debug << endl;
            *debug << "dydp is ";
            copy(r2.begin(), r2.begin() + r1.size(), 
                ostream_iterator<double>(*debug, " ")); 
            *debug << endl;
          }
          for (size_t i = 0; i < r1.size(); ++i)
            sum += r1[i]*r2[i];
          break;
        }
        default:
          throw gError("InternalVariables: function is not defined");
      }
    }
    else
      switch (f)
      {
        case ::G:
        case ::H:
        case ::S:
        case ::V:
        {
          double z = Zin(f, Tp, x);
          if (debug)
            *debug << "Zin is " << z << endl;
          sum += z;
          switch (f)
          {
            case ::G:
              break;
            case ::H:
            {
              double s = 0;
              const vec_double &r1 = dZindy(::G, Tp, x);
              const vec_double &r2 = dydT(Tp, x);
              if (debug)
              {
                *debug << "dZindy(G) is ";
                copy(r1.begin(), r1.end(), 
                    ostream_iterator<double>(*debug, " ")); 
                *debug << endl;
                *debug << "dydT is ";
                copy(r2.begin(), r2.begin() + r1.size(), 
                    ostream_iterator<double>(*debug, " ")); 
                *debug << endl;
              }
              for (size_t i = 0; i < r1.size(); ++i)
                s += r1[i]*r2[i];
              sum += -Tp.T()*s;
              break;
            }
            case ::S:
            {
              const vec_double &r1 = dZindy(::G, Tp, x);
              const vec_double &r2 = dydT(Tp, x);
              if (debug)
              {
                *debug << "dZindy(G) is ";
                copy(r1.begin(), r1.end(), 
                    ostream_iterator<double>(*debug, " ")); 
                *debug << endl;
                *debug << "dydT is ";
                copy(r2.begin(), r2.begin() + r1.size(), 
                    ostream_iterator<double>(*debug, " ")); 
                *debug << endl;
              }
              for (size_t i = 0; i < r1.size(); ++i)
                sum += -r1[i]*r2[i];
              break;
            }
            case ::V:
            {
              const vec_double &r1 = dZindy(::G, Tp, x);
              const vec_double &r2 = dydp(Tp, x);
              if (debug)
              {
                *debug << "dZindy(G) is ";
                copy(r1.begin(), r1.end(), 
                    ostream_iterator<double>(*debug, " ")); 
                *debug << endl;
                *debug << "dydp is ";
                copy(r2.begin(), r2.begin() + r1.size(), 
                    ostream_iterator<double>(*debug, " ")); 
                *debug << endl;
              }
              for (size_t i = 0; i < r1.size(); ++i)
                sum += r1[i]*r2[i];
              break;
            }
          }
          break;
        }
        case ::Cp:
        {
          Z_T num(::H, *this, ::mix, Tp.p(), x);
          double z = FirstDerivative(num, Tp.T(), step);
          if (debug)
            *debug << "FirstDerivative is " << z << endl;
          sum += z;
          break;
        }
        case ::dVdT:
        {
          Z_T num(::V, *this, ::mix, Tp.p(), x);
          double z = FirstDerivative(num, Tp.T(), step);
          if (debug)
            *debug << "FirstDerivative is " << z << endl;
          sum += z;
          break;
        }
        case ::dVdp:
        {
          Z_p num(::V, *this, ::mix, Tp.T(), x);
          double z = FirstDerivative(num, Tp.p(), step);
          if (debug)
            *debug << "FirstDerivative is " << z << endl;
          sum += z;
          break;
        }
        default:
          throw gError("InternalVariables: function is not defined");
      }
  }
  if (debug)
    *debug << endl;
  return sum;
}

const vec_double& InternalVariables::dZdx_y(function f, index i, 
    const StateTp &Tp, const StateX &x, const vec_double &y) const
{
  if (debug)
    *debug << "InternalVariables::dZdx_y function is " << f << " index is " 
      << i << endl;
  fill(vz.begin(), vz.end(), 0.);
  if ((i & ::mix) == ::mix)
  {
    CopyInternalVariables(y);
    if (debug)
    {
      *debug << "Internal Variables" << endl;
      copy(y.begin(), y.end(), ostream_iterator<double>(*debug, " "));  
      *debug << endl;
    }
    if (same)
      switch (f)
      {
        case ::G:
         break;
        case ::H:
        case ::S:
        case ::V:
        {
          const matrix &m1 = dydx(Tp, x);
          const vec_double &r1 = dZindy(f, Tp, x);
          if (debug)
          {
            *debug << "dydx" << endl;
            for (size_t j = 0; j < r1.size(); ++j)
            {
              for (size_t i = 0; i < size(); ++i)
                *debug << setw(12) << m1(j, i);
              *debug << endl;
            }
            *debug << "dZindy(" << f << ") is ";
            copy(r1.begin(), r1.end(), ostream_iterator<double>(*debug, " "));  
            *debug << endl;
          }
          for (size_t i = 0; i < size(); ++i)
            for (size_t j = 0; j < r1.size(); ++j)
              vz[i] += r1[j]*m1(j, i);
          break;
        }
        default:
          throw gError("InternalVariables: function is not supported");
      }
    else
      switch (f)
      {
        case ::G:
        {
          const matrix &m1 = dydx(Tp, x);
          const vec_double &r1 = dZindy(f, Tp, x);
          if (debug)
          {
            *debug << "dydx" << endl;
            for (size_t j = 0; j < r1.size(); ++j)
            {
              for (size_t i = 0; i < size(); ++i)
                *debug << setw(12) << m1(j, i);
              *debug << endl;
            }
            *debug << "dZindy(" << f << ") is ";
            copy(r1.begin(), r1.end(), ostream_iterator<double>(*debug, " "));  
            *debug << endl;
          }
          for (size_t i = 0; i < size(); ++i)
            for (size_t j = 0; j < r1.size(); ++j)
              vz[i] += r1[j]*m1(j, i);
          break;
        }
        case ::H:
        {
          ZA_x num(::H, *this, ::mix, Tp);
          ArrayFirstDerivative(num, x, vz, step);
          if (debug)
          {
            *debug << "ArrayFirstDer is ";
            copy(vz.begin(), vz.end(), ostream_iterator<double>(*debug, " "));  
            *debug << endl;
          }
          break;
        }
        case ::S:
        {
          ZA_x num(::S, *this, ::mix, Tp);
          ArrayFirstDerivative(num, x, vz, step);
          if (debug)
          {
            *debug << "ArrayFirstDer is ";
            copy(vz.begin(), vz.end(), ostream_iterator<double>(*debug, " "));  
            *debug << endl;
          }
          break;
        }
        case ::V:
        {
          ZA_x num(::V, *this, ::mix, Tp);
          ArrayFirstDerivative(num, x, vz, step);
          if (debug)
          {
            *debug << "ArrayFirstDer is ";
            copy(vz.begin(), vz.end(), ostream_iterator<double>(*debug, " "));  
            *debug << endl;
          }
          break;
        }
        default:
          throw gError("InternalVariables: function is not supported");
      }
    if (same || f == ::G)
    {
      dZindx(f, Tp, x, vz);
      if (debug)
      {
        *debug << "after dZindx the vector is ";
        copy(vz.begin(), vz.end(), ostream_iterator<double>(*debug, " "));  
        *debug << endl;
      }
    }
  }
  if (i & ::ref)
  {
    dZexdx(f, Tp, x, vz);
    if (debug)
    {
      *debug << "after dZexdx vz is ";
      copy(vz.begin(), vz.end(), ostream_iterator<double>(*debug, " "));  
      *debug << endl;
    }
  }
  if (debug)
    *debug << endl;
  return vz;
}

