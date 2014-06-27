/*
Copyright (C) 1998, 1999, 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                       http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#include "simple.h"
#include <algorithm>
#include <typeinfo>

void RegisterSimple()
{
  phase::RegisterType("pol_solution", new SimpleSolution);
  phase::RegisterType("SimpleSolution", new SimpleSolution);
}

void index_x::init()
{
  for (size_t i = 0; i < m; ++i)
    in[i] = i;
  ilin = 0;
  nm = bin_coef(m, n);
}

void index_x::operator++()
{
  for (size_t i = 0; i < m; ++i)
    if (++(in[m - i - 1]) < n - i)
    {
      for (size_t j = m - i; j < m; ++j)
        in[j] = in[j - 1] + 1;
      ++ilin;
      return;
    }
  ilin = nm;
  return;
}

void SimpleSolution::insert(const FuncTpx &y)
{
  vec_int v = y.GetVars();
  sort(v.begin(), v.end());
  if (v.size() > size())
    throw gError("SimpleSolution::insert - wrong1");
  if (v.back() >= (int)size())
    throw gError("SimpleSolution::insert - wrong2");
  vec_FuncTpx_i i;
  for (i = p.begin(); i != p.end(); ++i)
  {
    if (y.GetIndex() > (*i).GetIndex())
      break;
    else
      continue;
    if (v.size() < (*i).size()) 
      break;
    else if ((v.size() == (*i).size()) && (v < (*i).GetVars()))
      break;
  }
  p.insert(i, y);
}

void SimpleSolution::set(size_t nc, size_t m_max)
{
  if (!(nc > 1 && m_max <= nc)) 
    throw gError("SimpleSolution::set - wrong input");
  clear();
  vf.resize(nc);
  vz.resize(nc);
  vec_int c(nc);
  string id = "_a";
  for (size_t i = 0; i < c.size(); ++i)
  {
    id[1] = '1' + i;
    vf[i] = id.c_str();
    c[i] = i;
  }
  FuncTpx RK;
  RK.create("Reference");
  RK.set(c);
  p.push_back(RK);
  RK.create("IdealMixing");
  RK.set(c);
  p.push_back(RK);
  RK.create("Polynomial");
  for (size_t i = 2; i <= m_max; ++i)
    for (index_x ix(i, nc); ix; ++ix)
    {
      RK.set(ix, func_Tp("Cp_zero"));
      p.push_back(RK);
    }
}

void SimpleSolution::read(const SGML &e_)
{
  clear();
  SGML e = e_;
  FuncTpx t;
  ReadComponents(vf, e, t);
  if (t.size())
    insert(t);
  parser par(e);
  SGML el;
  par.GetSGML(el);
  if (el.name.empty())
  {
// compatibility    
    if (par.GetToken() == "Gid")
    {
      el.name = "IdealMixing";
      t.read(el, &vf);
      insert(t);
      par.GetSGML(el);
    }
  }
  while (!el.name.empty())
  {
    t.read(el, &vf);
    if (t.size() == size())
      t.SetFormalism(::NoChange);
    insert(t);
    par.GetSGML(el);
  }
  mat.resize(size());
}

void SimpleSolution::WriteBody(ostream &out, size_t shift) const
{
  WriteComponents(out, shift);
  for (vec_FuncTpx_ci i = p.begin(); i != p.end(); ++i)
    (*i).write(out, &vf, shift);
}

const IdealMixing& SimpleSolution::FindIdealMixing()
{
  for (vec_FuncTpx_ci i = p.begin(); i != p.end(); ++i)
    if (typeid(*(*i).pointer()).name() == typeid(IdealMixing).name())
      return dynamic_cast<const IdealMixing&>(*(*i).pointer());
  throw gError("SimpleSolution: no IdealMixing term");
}

double SimpleSolution::Z(function f, index idx, const StateTp &Tp,
    const StateX &x) const
{
  double Sum = 0.;
  for (vec_FuncTpx_ci i = p.begin(); i != p.end(); ++i)
    if ((idx & (*i).GetIndex()) == (*i).GetIndex())
      Sum += (*i).Z(f, Tp, x);
 return Sum;
}

const vec_double& SimpleSolution::dZdx(function f, index idx, 
    const StateTp &Tp, const StateX &x) const
{
  fill(vz.begin(), vz.end(), 0.);
  for (vec_FuncTpx_ci i = p.begin(); i != p.end(); ++i)
    if ((idx & (*i).GetIndex()) == (*i).GetIndex())
      (*i).dZdx(f, Tp, x, vz);
  return vz;
}

const SymMatrix& SimpleSolution::d2Zdx2(function f, index idx, 
    const StateTp &Tp, const StateX &x) const
{
  fill(mat.begin(), mat.end(), 0.);
  for (vec_FuncTpx_ci i = p.begin(); i != p.end(); ++i)
    if ((idx & (*i).GetIndex()) == (*i).GetIndex())
      (*i).d2Zdx2(f, Tp, x, mat);
  return mat;
}

