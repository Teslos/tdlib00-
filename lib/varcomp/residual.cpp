/*
Copyright (C) 1999 - 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#include "post_an.h"
#include "residual.h"

map_residual res;

static bool WriteElem = false;
static bool WriteGl = false;

void ReadModelFile(istream &in, const string &str, vec_of &vec)
{
  parser p(in, str);
  SGML e;
  coef cf;
  func_Tp fTp;
  FuncTpx fTpx;
  species sps;
  phase ph;
  algorithm alg;
  residual r;
  output o;
  OutputFile of;
  while (!p.eof())
  {
    p.GetSGML(e);
    const string &cls = e.FindString("class");
    if (e.name == "elements")
    {
      WriteElem = true;
      elem::ReadElements(e);
    }
    else if (e.name == "globals")
    {
      WriteGl = true;
      ReadGlobals(e);
    }
    else if (e.name == "coef")
      cf.read(e);
// compatibility    
    else if (e.name == "PhaseEquilibrium")
      alg.read(e);
    else if (cls == "func_Tp" || e.name == "func_Tp")
      fTp.read(e);
    else if (cls == "FuncTpx" || e.name == "func_Tpx")
      fTpx.read(e);
    else if (e.name == "species")
      sps.read(e);
    else if (cls == "phase" || e.name == "phase")
      ph.read(e);
    else if (cls == "algorithm" || e.name == "algorithm")
      alg.read(e);
    else if (e.name == "residual" || e.name == "PhaseResidual")
    {
      r.read(e);
      if (!r.id().empty())
        res.insert(map_residual::value_type(r.id(), r));
      else
        cout << "anonymous residual - ignored" << endl;
    }
    else if (cls == "output" || e.name == "output")
      o.read(e);
    else if (e.name == "OutputFile")
    {
      of.read(e);
      vec.push_back(of);
    }
    else if (e.name == "")
      break;
    else
    {
      try
      {
        fTp.read(e);
      }
      catch (gError)
      {
        try
        {
          fTpx.read(e);
        }
        catch (gError)
        {
          try
          {
            ph.read(e);
          }
          catch (gError)
          {
            try
            {
              alg.read(e);
            }
            catch (gError)
            {
              try
              {
                o.read(e);
              }
              catch (gError)
              {
                cout << "not known name - " << e.name << " - skipped" << endl;
              }
            }
          }
        }
      }
    }
  }
}

void WriteModelFile(const vec_of &vec)
{
  {
    if (!file_name.empty())
    {
      ofstream out((file_name + ".gl.mod").c_str());
      WriteGlobals(out);
    }
    else
      WriteGlobals(cout);
  }
  if (WriteElem)
  {
    if (!file_name.empty())
    {
      ofstream out((file_name + ".el.mod").c_str());
      elem::WriteElements(out);
    }
    else
      elem::WriteElements(cout);
  }
  {
    if (!file_name.empty())
    {
      ofstream out((file_name + ".sys.mod").c_str());
      out.precision(12);
      phase::WriteAll(out);
    }
    else
      phase::WriteAll(cout);
  }
  {
    if (!file_name.empty())
    {
      ofstream out((file_name + ".alg.mod").c_str());
      algorithm::WriteAll(out);
    }
    else
      algorithm::WriteAll(cout);
  }
  if (!res.empty())
  {
    if (!file_name.empty())
    {
      ofstream out((file_name + ".res.mod").c_str());
      for (map_residual_i i = res.begin(); i != res.end(); ++i)
        (*i).second.write(out);
    }
    else
      for (map_residual_i i = res.begin(); i != res.end(); ++i)
        (*i).second.write(cout);
  }
  if (!vec.empty())
  {
    if (!file_name.empty())
    {
      ofstream out((file_name + ".out.mod").c_str());
      for (vec_of_ci i = vec.begin(); i != vec.end(); ++i)
        (*i).write(out);
    }
    else
      for (vec_of_ci i = vec.begin(); i != vec.end(); ++i)
        (*i).write(cout);
  }
}

void residual::read(const SGML &e)
{
  if (e.name != "PhaseResidual")
    e.compare("residual");
  clear();
  sc = e.FindDouble("ScaleOfX", 1.);
  id_ = e.FindString("id");
  yname = e.FindString("yname");
  if (yname.empty())
    throw gError("residual: no yname");
  xname = e.FindString("xname");
  SGML el;
  parser p(e);
  comp.read(p.GetSGML(el));
}

ostream& residual::write(ostream& out, size_t shift) const
{
  out << PutTab(shift) << "<residual ";
  if (!id_.empty())
    out << "id=" << id_ << endl;
  out << PutTab(shift + 1) << "yname=" << yname << endl;
  if (!xname.empty())
    out << PutTab(shift + 1) << "xname=" << xname << endl;
  out << PutTab(shift + 1) << "ScaleOfX=" << sc << ">" << endl;
  comp.write(out, shift + 1);
  out << PutTab(shift) << "</residual>" << endl;
  return out;
}

double residual::res(double *reg) const
{
  try
  {
    return reg[iy] - f(reg);
  }
  catch (CanNotCompute &t)
  {
    return t.fmin;
  }
}

residual& FindResidual(const string &str)
{
  map_residual_i i = res.find(str);
  if (i == res.end())
    throw gError(string("no such residual: - ").append(str));
  else
    return (*i).second;
}


