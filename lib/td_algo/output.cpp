/*
Copyright (C) 1998, 1999, 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                       http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#include "output.h"
#include <typeinfo> 
#include <fstream> 
#include <iomanip> 

static map_string_string *type;

RefOutput *output::PTR_NULL;
output::map_types *output::types = 0;
output::map_types *output::obj = 0;
output::Destruct output::clean;

void RegisterPostAnalysis();

void (*output::init_ptr[])() = {RegisterPostAnalysis, 0};

const char *comparison::OperationName[] = {"LE", "EQ", "GE", 0};

int OutputFile::width = 12;

istream& RefOutput::read(istream& in)
{
  output::init();
  SGML e(in);
  e.compare((*type)[typeid(*this).name()]);
  read(e);
  return in;
}

ostream& RefOutput::write(ostream& out, size_t shift) const
{
  output::init();
  out << PutTab(shift) << "<" << (*type)[typeid(*this).name()] 
    << " class=output";
  WriteAttributes(out, shift + 1);
  out << ">" << endl;
  WriteBody(out, shift + 1);
  out << PutTab(shift) << "</" << (*type)[typeid(*this).name()] << "> " << endl;
  return out;
}

output::Destruct::~Destruct()
{
  if (output::types)
  {
    delete output::types;
    delete output::obj;
    output::types = 0;
    output::obj = 0;
    delete type;
  }
}

void output::init()
{
  if (!types)
  {
    type = new map_string_string;
    PTR_NULL = new NullOutput;
    types = new map_types;
    obj = new map_types;
    RegisterType("NullOutput", PTR_NULL);
    RegisterType("ComputeOutput", new ComputeOutput);
    for (size_t i = 0; init_ptr[i]; ++i)
      (*init_ptr[i])();
  }
}

void output::RegisterType(const string& name, RefOutput* ptr)
{
  init();
  (*types)[name] = output(ptr);
  (*type)[typeid(*ptr).name()] = name;
}

void output::find(const string &name)
{
  map_types_i i = obj->find(name);
  if (i == obj->end())
  {
    throw gError(string("output: id ").append(ptr->id).append
      (" is not defined"));
  }
  else
    *this = (*i).second;
}

void output::read(const SGML &e)
{
  string id;
  if (!(id = e.FindString("IDREF")).empty())
  {
    find(id);
  }
  else
  {
    create(e.name);
    cow();
    ptr->read(e);
    if (!(ptr->id = e.FindString("id")).empty())
    {
      pair<map_types_i, bool> i;
      i = obj->insert(map_types::value_type(ptr->id, *this));
      if (!i.second)
      {
        cout << "output: " << ptr->id << " is already defined - "
          << "set to anonymous" << endl;
        ptr->id.erase();
      }
    }
  }
}

ostream& output::write(ostream& out, size_t shift) const
{
  if (ptr->id.empty() || ptr->print)
  {
    out << PutTab(shift) << "<" << (*type)[typeid(*ptr).name()] 
      << " class=output";
    if (!ptr->id.empty())
    {
      out << " id=" << ptr->id;
    }
    ptr->WriteAttributes(out, shift + 1);
    out << ">" << endl;
    ptr->WriteBody(out, shift + 1);
    out << PutTab(shift) << "</" << (*type)[typeid(*ptr).name()] << "> " 
      << endl;
    ptr->print = false;
  }
  else
  {
    out << PutTab(shift) << "<output class=output IDREF="
      << ptr->id << "></output>"
      << endl;
  }
  return out;
}

void CycleOutput::ReadStart(const SGML &e)
{
  parser p(e);
  SGML el;
  p.GetSGML(el);
  if (el.name == "compute")
  {
    start.read(el);
    nm_cycle = start.names();
    vds.resize(start.size());
    p.GetSGML(el);
  }
  while (1)
  {
    if (el.name == "var")
    {
      const string &s = el.FindString("name");
      bool exist = false;
      for (size_t i = 0; i < nm_cycle.size(); ++i)
        if (nm_cycle[i] == s)
        {
          cout << "CycleOutput: " << s << " is already defined - ignored" 
            << endl;
          exist = true;
          break;
        }
      if (exist)
      {
        p.GetSGML(el);
        continue;
      }
      nm_cycle.push_back(s);
      vds.push_back(el.FindDouble("value"));
    }
    else
      break;
    p.GetSGML(el);
  }
  vd.resize(vds.size());
}

void CycleOutput::ReadFinish(const SGML &e)
{
  parser p(e);
  SGML el;
  p.GetSGML(el);
  if (el.name == "compute")
  {
    finish.read(el);
    const vec_string& vs = finish.names();
    for (size_t j = 0; j < vs.size(); ++j)
    {
      vif.push_back(-1);
      vof.push_back(comparison("GE"));
      vif.back() = SearchString(nm_cycle, vs[j]);
      if (vif.back() == -1)
        throw gError(string("CycleOutput: no such id in the start - ")
            .append(vs[j]));
    }
    vdf.resize(vif.size());
    p.GetSGML(el);
  }
  iout = jout = -1;
  while (1)
  {
    if (el.name == "var")
    {
      const string &s = el.FindString("name");
      bool exist = false;
      for (size_t i = 0; i < vif.size(); ++i)
        if (nm_cycle[vif[i]] == s)
        {
          vof[i] = el.FindString("operation");
          exist = true;
          break;
        }
      if (exist)
      {
        p.GetSGML(el);
        continue;
      }
      for (size_t i = 0; i < nm_cycle.size(); ++i)
        if (nm_cycle[i] == s)
        {
          vif.push_back(i);
          vdf.push_back(el.FindDouble("value"));
          vof.push_back(el.FindString("operation"));
          exist = true;
          break;
        }
      if (!exist)
        throw gError(string("CycleOutput: no such id in the start - ")
            .append(s));
    }
    else if (el.name == "compare")
    {
      comp = el.FindString("operation");
      string s1 = el.FindString("OutName1");
      string s2 = el.FindString("OutName2");
      for (size_t i = 0; i < nm.size(); ++i)
        if (nm[i] == s1)
          iout = i;
        else if (nm[i] == s2)
          jout = i;
      if (iout == -1)
        throw gError(string("CycleOutput: no such id in the output - ")
            .append(s1));
      if (jout == -1)
        throw gError(string("CycleOutput: no such id in the output - ")
            .append(s2));
    }
    else
      break;
    p.GetSGML(el);
  }
}

void CycleOutput::ReadStep(const SGML &e)
{
  parser p(e);
  SGML el;
  convert conv;
  while (!p.eof())
  {
    p.GetSGML(el);
    if (el.name == "var")
    {
      const string &s = el.FindString("name");
      bool exist = false;
      for (size_t i = 0; i < nm_cycle.size(); ++i)
        if (nm_cycle[i] == s)
        {
          vic.push_back(i);
          istringstream in(el.body);
          in >> conv;
          vc.push_back(conv);
          vc.back().SetID(nm_cycle);
          exist = true;
          break;
        }
      if (!exist)
        throw gError(string("CycleOutput: no such id in the start - ")
            .append(s));
    }
    else
      break;
  }
}

/*
<start>
  <compute> </compute>
  <var name=  value= ></var>
</start>
<finish>
  <compute> </compute>
  <var name= value= operation=></var>
  <compare OutName1= OutName2= operation=></compare>
</finish>
<step>
<var name= ><convert></comvert></var>
</step>
*/

void CycleOutput::WriteCycle(ostream &out, size_t shift) const
{
  out << PutTab(shift) << "<start>" << endl;
  if (start.size())
    start.write(out, shift + 1);
  for (size_t i = start.size(); i < nm_cycle.size(); ++i)
    out << PutTab(shift + 1) << "<var name=" << nm_cycle[i] << " value=" 
      << vds[i] << "></var>" << endl;
  out << PutTab(shift) << "</start>" << endl;
  out << PutTab(shift) << "<finish>" << endl;
  if (finish.size())
    finish.write(out, shift + 1);
  for (size_t i = finish.size(); i < vdf.size(); ++i)
    out << PutTab(shift + 1) << "<var name=" << nm_cycle[vif[i]] << " value=" 
      << vdf[i] << " operation=" << vof[i] << "></var>" << endl;
  if (iout != jout)
    out << PutTab(shift + 1) << "<compare OutName1=" << nm[iout] 
      << " OutName2=" << nm[jout] << " operation=" << comp << "></compare>" 
      << endl;
  out << PutTab(shift) << "</finish>" << endl;
  out << PutTab(shift) << "<step>" << endl;
  for (size_t i = 0; i < vc.size(); ++i)
  {
    out << PutTab(shift + 1) << "<var name=" << nm_cycle[vic[i]] <<">" << endl;
    vc[i].write(out, shift + 2);
    out << PutTab(shift + 1) << "</var>" << endl;
  }
  out << PutTab(shift) << "</step>" << endl;
}

void CycleOutput::est() const
{
  if (start.size())
  {
    try
    {
      const vec_double &r = start.vec();
      copy(r.begin(), r.end(), vds.begin());
    }
    catch (CanNotCompute)
    {
      cout << "CycleOutput: can not cumpute start - ignored" << endl;
      return;
    }
  }
  abort = -1;
  if (finish.size())
  {
    try
    {
      const vec_double &r = finish.vec();
      copy(r.begin(), r.end(), vdf.begin());
    }
    catch (CanNotCompute)
    {
      cout << "CycleOutput: can not cumpute finish - ten outputs" << endl;
      abort = 10;
    }
  }
  SetOnce();
  vd = vds;
  bool going = true;
  while (going)
  {
    if (abort == -1)
    {
      for (size_t i = 0; i < vdf.size(); ++i)
        if (vof[i](vd[vif[i]], vdf[i]))
        {
          vd[vif[i]] = vdf[i];
          going = false;
        }
    }
    else
    {
      if (--abort == 0)
        break;
    }
    try
    {
      vals.insert(vals.end(), NCols(), 0.);
      double *x = &*vals.end() - NCols();
      ComputeCurrentPoint(&*vd.begin(), x);
    }
    catch (CanNotCompute)
    {
//      fill(vals.end() - NCols(), vals.end(), 0./0.);
      fill(vals.end() - NCols(), vals.end(), HUGE_VAL);
    }
    if (iout != jout)
      if (comp(*(vals.end() - NCols() + iout), *(vals.end() - NCols() + jout)))
        break;
    for (size_t i = 0; i < vc.size(); ++i)
      vd[vic[i]] = vc[i](&*vd.begin());
  }
}

void ComputeOutput::read(const SGML &e)
{
  clear();
  if (e.defined("ChangeXY"))
    change = true;
  else
    change = false;
  SGML start, finish, step;
  parser p(e);
  SGML el;
  compute comp;
  while (!p.eof())
  {
    p.GetSGML(el);
    if (el.name == "start")
      start = el;
    else if (el.name == "finish")
      finish = el;
    else if (el.name == "step")
      step = el;
    else if (!el.name.empty())
    {
      comp.read(el);
      nc += comp.size();
      const vec_string &vs = comp.names();
      nm.insert(nm.end(), vs.begin(), vs.end());
      vec.push_back(comp);
    }
    else
      break;
  }
  if (start.name.empty())
  {
    NoCycle = true;
    vals.insert(vals.end(), NCols(), 0.);
  }
  else
  {
    NoCycle = false;
    ReadStart(start);
    ReadFinish(finish);
    ReadStep(step);
    for (vec_compute_ci i = vec.begin(); i != vec.end(); ++i)
    {
      (*i).SetOnceInput(nm_cycle);
      (*i).SetInput(nm_cycle);
    }
  }
}

void ComputeOutput::WriteBody(ostream &out, size_t shift) const
{
  if (!NoCycle)
    WriteCycle(out, shift);
  for (vec_compute_ci i = vec.begin(); i != vec.end(); ++i)
    (*i).write(out, shift);
}

void ComputeOutput::est() const
{
  if (NoCycle)
    ComputeCurrentPoint(0, &*vals.begin());
  else
    CycleOutput::est();
}

void ComputeOutput::ComputeCurrentPoint(double *reg, double *x) const
{
  for (vec_compute_ci i = vec.begin(); i != vec.end(); ++i)
  {
    try
    {
      const vec_double &r = (*i).vec(reg);
      copy(r.begin(), r.end(), x);
      x += (*i).size();
    }
    catch (CanNotCompute)
    {
      for (size_t ii = 0; ii < (*i).size(); ++ii)
//        *x++ = 0./0.;
        *x++ = HUGE_VAL;
    }
  }     
}

void OutputFile::read(const SGML &e)
{
  clear();
  format = e.FindString("format", "file");
  e.compare("OutputFile");
  ext = e.FindString("ext", "axm");
  output o;
  parser p(e);
  SGML el;
  while (!p.eof())
  {
    o.read(p.GetSGML(el));
    vo.push_back(o);
  }
}

ostream& OutputFile::write(ostream& out, size_t shift) const
{
  out << PutTab(shift) << "<OutputFile ext=" << ext 
    << " format=" << format << ">" << endl;
  for (vec_output_ci i = vo.begin(); i != vo.end(); ++i)
    (*i).write(out, shift + 1);
  out << PutTab(shift) << "</OutputFile>" << endl;
  return out;
}

void OutputFile::out(const string &basename) const
{
  for (vec_output_ci i = vo.begin(); i != vo.end(); ++i)
    (*i).est();
  ostream *out_ptr;
  ofstream out;
  if (basename.empty())
    out_ptr = &cout;
  else
  {
    out.open((basename + "." + ext).c_str());
    out_ptr = &out;
  }
  ostream &of = *out_ptr;
  if (format == "axum")
  {
    for (vec_output_ci i = vo.begin(); i != vo.end(); ++i)
      for (size_t j = 0; j < (*i).NCols(); ++j)
        of << setw(width) << (*i).name(j).c_str() << " ";
    of << endl;
    size_t NRows = 0;
    for (vec_output_ci i = vo.begin(); i != vo.end(); ++i)
      if ((*i).NRows() > NRows)
        NRows = (*i).NRows();
    for (size_t i = 0; i < NRows; ++i)
    {
      for (vec_output_ci j = vo.begin(); j != vo.end(); ++j)
      {
        if (i < (*j).NRows())
          for (size_t k = 0; k < (*j).NCols(); ++k)
            of << setw(width) << (*j).value(i, k) << " ";
        else
          for (size_t k = 0; k < (*j).NCols(); ++k)
            of << setw(width) << "miss" << " ";
      }
      of << endl;
    }
  }
  else if (format == "file")
  {
    for (vec_output_ci i = vo.begin(); i != vo.end(); ++i)
    {
      for (size_t j = 0; j < (*i).NCols(); ++j)
        of << setw(width) << (*i).name(j).c_str() << " ";
      of << endl;
      for (size_t j = 0; j < (*i).NRows(); ++j)
      {
        for (size_t k = 0; k < (*i).NCols(); ++k)
          of << setw(width) << (*i).value(j, k) << " ";
        of << endl;
      }
      of << endl;
    }
  }
  else
  {
    of << "plot \\" << endl;
    for (vec_output_ci i = vo.begin(); i != vo.end(); ++i)
    {
      for (size_t j = 1; j < (*i).NCols(); ++j)
      {
        of << "'-' title \"" << (*i).name(j) << "\" ";
        if ((*i).IsLine())
          of << "with line";
        if ((*i).name(j).empty())
          of << " lt -1";
        if (i != vo.end() - 1 || j != (*i).NCols() - 1)
          of << ", \\";
        of << endl;
      }
    }
    for (vec_output_ci j = vo.begin(); j != vo.end(); ++j)
    {
      for (size_t k = 1; k < (*j).NCols(); ++k)
      {
        for (size_t i = 0; i < (*j).NRows(); ++i)
          if (finite((*j).value(i, 0)) && finite((*j).value(i, k)))
          {
            if ((*j).ChangeXY())
              of << (*j).value(i, k) << " " << (*j).value(i, 0) << endl;
            else
              of << (*j).value(i, 0) << " " << (*j).value(i, k) << endl;
          }
        of << "e" << endl;
      }
    }
  }
}

