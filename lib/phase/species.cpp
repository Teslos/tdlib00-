/*
Copyright (C) 1998, 1999, 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                       http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#include <lapack.h>
#include <algorithm>
#include <typeinfo>
#include <float.h>
#include "species.h"
#include "func_tp_imp.h"


Ref_species* species::PTR_NULL;
species::map_species *species::s = 0;
species::Destruct species::clean;

double GetNu::tol = DBL_EPSILON*1000.;

void species::read(const SGML &e)
{
  clear();
  if (!(ptr->id = e.FindString("IDREF")).empty())
  {
    map_species_ci j = s->find(ptr->id);
    if (j == s->end())
      throw gError(string("species: id ").append(ptr->id).append
          (" is not defined"));
    else
      *this = (*j).second;
  }
  else
  {
    istringstream in(e.body);
    in >> ptr->fml;
    SGML el(in);
    if (el.FindString("class") == "func_Tp")
    {
      ptr->f.read(el);
      el.read(in);
    }
    bool FindCoefs = true;
    if (!el.name.empty())
    {
      el.compare("ref_plane");
      species t_r;
      double t_c;
      parser p(el);
      SGML ele;
      while (!p.eof())
      {
        t_c = p.GetDouble();
        if (t_c)
          FindCoefs = false;
        t_r.read(p.GetSGML(ele));
        ptr->c.push_back(t_c);
        ptr->r.push_back(t_r);
      }
    }
    if (!ptr->r.empty() && FindCoefs)
    {
      GetNu A(ptr->r);
      ptr->r = A.basis();
      ptr->c.resize(ptr->r.size());
      A.nu(ptr->fml, &*ptr->c.begin());
      for (int i = ptr->c.size() - 1; i >= 0; --i)
        if (fabs(ptr->c[i]) < GetNu::tol)
        {
          ptr->c.erase(ptr->c.begin() + i);
          ptr->r.erase(ptr->r.begin() + i);
        }
    }
    if (!(ptr->id = e.FindString("id")).empty())
    {
      pair<map_species_i, bool> i;
      i = s->insert(map_species::value_type(ptr->id, *this));
      if (!i.second)
      {
        cout << "species: " << ptr->id << " is already defined - "
          << "set to anonymous" << endl;
        ptr->id.erase();
      }
    }
  }
}

ostream& species::write(ostream& out, size_t shift) const
{
  if (ptr->id.empty() || ptr->p)
  {
    out << PutTab(shift) << "<species";
    if (!ptr->id.empty())
      out << " id=" << ptr->id;
    out << "> " << ptr->fml << endl;
    if (typeid(*ptr->f.pointer()).name() != typeid(null_Tp).name())
      ptr->f.write(out, shift + 1);
    if (!ptr->r.empty())
    {
      out << PutTab(shift + 1) << "<ref_plane>" << endl;
      for (size_t i = 0; i < ptr->r.size(); ++i)
      {
        out << PutTab(shift + 2) << ptr->c[i] << endl;
        ptr->r[i].write(out, shift + 2);
      }
      out << PutTab(shift + 1) << "</ref_plane>" << endl;
    }
    out << PutTab(shift) << "</species>" << endl;
    ptr->p = false;
  }
  else
    out << PutTab(shift) << "<species IDREF=" << ptr->id << "></species>" << endl;
  return out;
}

double species::Z(function fun, index i, const StateTp &Tp) const
{
  double sum = 0.;
  if (i & ::ref)
    for (size_t i = 0; i < ptr->r.size(); ++i)
      sum += ptr->c[i]*ptr->r[i].Z(fun, Tp);
  if (i & ::excess)
    sum += ptr->f.Z(fun, Tp);
  return sum;
}

void GetNu::reset(const vec_species &r)
{
  from_species = true;
  Nel = els.size();
  Nc = min(Nel, r.size());
  size_t left = Nc;
  bs.assign(r.begin(), r.begin() + Nc);
  while (1) 
  {
    SetFormulaMatrix();
    size_t dep = 0;
    for (size_t i = 0; i < Nc; ++i)
      if (fabs(fm[i + i*Nel]) < tol)
      {
        bs.erase(bs.begin() + i - dep);
        ++dep;
        break;
      }
    if (dep)
    {
      size_t end = left + min(dep, r.size() - left);
      bs.insert(bs.end(), r.begin() + left, r.begin() + end);
      left = end;
    }
    else
      break;
  }
}

const vec_int& GetNu::reset(const vec_formula &r)
{
  from_species = false;
  Nel = els.size();
  Nc = min(Nel, r.size());
  size_t left = Nc;
  fbs.assign(r.begin(), r.begin() + Nc);
  indx.resize(Nc);
  for (size_t i = 0; i < Nc; ++i)
    indx[i] = i;
  while (1) 
  {
    SetFormulaMatrix();
    size_t dep = 0;
    for (size_t i = 0; i < Nc; ++i)
      if (fabs(fm[i + i*Nel]) < tol)
      {
        fbs.erase(fbs.begin() + i - dep);
        indx.erase(indx.begin() + i - dep);
        ++dep;
        break;
      }
    if (dep)
    {
      size_t end = left + min(dep, r.size() - left);
      fbs.insert(fbs.end(), r.begin() + left, r.begin() + end);
      for (size_t i = left; i < end; ++i)
        indx.insert(indx.end(), i);
      left = end;
    }
    else
      break;
  }
  return indx;
}

void GetNu::SetFormulaMatrix() const
{
  if (from_species)
    Nc = bs.size();
  else
    Nc = fbs.size();
  fm.resize(Nel*Nc, 0);
  b.resize(Nel);
  ipiv.resize(Nel);
  if (from_species)
    for (size_t j = 0; j < Nc; ++j)
    {
      const vec_double &fv = bs[j].mf().FormulaVector();
      copy(fv.begin(), fv.end(), fm.begin() + j*Nel);
    }
  else
    for (size_t j = 0; j < Nc; ++j)
    {
      const vec_double &fv = fbs[j].mf().FormulaVector();
      copy(fv.begin(), fv.end(), fm.begin() + j*Nel);
    }
    dgetf2(Nel, Nc, &*fm.begin(), &*ipiv.begin(), info);
/*  
  cout << "fmpr" << endl;
  for (size_t i = 0; i < Nel; ++i)
  {
    for (size_t j = 0; j < Nc; ++j)
      cout << setw(8) << fm[i + j*Nel] << " ";
    cout << endl;
  }
*/
}

void GetNu::nu(const formula &frm, double *n) const
{
  if (Nc == 0)
    throw GetNuError("GetNu: basis is empty");
  if (Nel != elem::elements().size())
  {
    Nel = els.size();
    SetFormulaMatrix();
  }
  const vec_double &fv = frm.FormulaVector();
  copy(fv.begin(), fv.end(), b.begin());
  dgetrs(Nel, Nc, 1, &*fm.begin(), &*ipiv.begin(), &*b.begin());
  for (size_t i = 0; i < Nc; ++i)
    if (fabs(b[i]) < tol)
      b[i] = 0.;
  for (size_t i = Nc; i < Nel; ++i)
    if (fabs(b[i]) > tol)
      throw GetNuError(string("GetNu: can not construct ").append(ObjToString(frm)));
  copy(b.begin(), b.begin() + Nc, n);
}

void GetNu::FormulaMatrix(const vec_formula &vf, matrix &A)
{
  if (Nc == 0)
    throw GetNuError("GetNu: basis is empty");
  if (Nel != elem::elements().size())
  {
    Nel = els.size();
    SetFormulaMatrix();
  }
  if (Nel == Nc)
  {
    for (size_t j = 0; j < vf.size(); ++j)
    {
      const vec_double &fv = vf[j].mf().FormulaVector();
      copy(fv.begin(), fv.end(), A.begin() + j*Nel);
    }
    dgetrs(Nel, Nc, vf.size(), &*fm.begin(), &*ipiv.begin(), &*A.begin());
  }
  else
  {
    matrix A_(vf.size(), Nel);
    for (size_t j = 0; j < vf.size(); ++j)
    {
      const vec_double &fv = vf[j].mf().FormulaVector();
      copy(fv.begin(), fv.end(), A_.begin() + j*Nel);
    }
    dgetrs(Nel, Nc, vf.size(), &*fm.begin(), &*ipiv.begin(), &*A_.begin());
    for (size_t j = 0; j < vf.size(); ++j)
    {
      for (size_t i = Nc; i < Nel; ++i)
        if (fabs(A_(j, i)) > tol)
      throw GetNuError(string("GetNu: can not construct ").
          append(ObjToString(vf[j])));
      copy(A_.begin() + j*Nel, A_.begin() + j*Nel + Nc, A.begin() + j*Nc);
    }
  }
  for (size_t i = 0; i < A.size(); ++i)
    if (fabs(A[i]) < tol)
      A[i] = 0.;
}

