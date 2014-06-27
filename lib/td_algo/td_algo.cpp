/*
Copyright (C) 1998, 1999, 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                       http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#include <algorithm>
#include <typeinfo>
#include <iterator>
#include "td_algo_imp.h"

static map_string_string *type;

RefAlgorithm *algorithm::PTR_NULL;
algorithm::map_types *algorithm::types = 0;
algorithm::map_types *algorithm::obj = 0;
algorithm::Destruct algorithm::clean;

void RegisterPhaseEquilibrium();

void (*algorithm::init_ptr[])() = {RegisterPhaseEquilibrium, 0};

bool algorithm::debug = false;

void convert::AnalizeId(const string &id)
{
  curr_tok.key = DOUBLE_PTR_PTR;
  curr_tok.x_ptr_ptr = &v[id].ptr;
}

ostream& convert::write(ostream& out, size_t shift) const
{
  out << PutTab(shift) << "<convert";
  if (!name_.empty())
    out << " name=" << name_;
  out << ">" << endl;
  CalculatorWithCoefs::write(out, shift + 1);
  out << PutTab(shift) << "</convert>" << endl;
  return out;
}

void convert::SetID(const vec_string &vs) const
{
  for (varmap_i i = v.begin(); i != v.end(); ++i)
   (*i).second.i = -1;
  varmap_i j;
  for (size_t i = 0; i < vs.size(); ++i)
    if ((j = v.find(vs[i])) != v.end())
      (*j).second.i = i;
  for (varmap_i i = v.begin(); i != v.end(); ++i)
    if ((*i).second.i < 0)
    {
// compatibility     
      if ((*i).first[0] == 'x')
      {
        size_t ii = 0;
        if ((*i).first.length() > 1)
        {
          ii = atoi((*i).first.c_str() + 1);
          if (ii == 0)
            throw gError("convert: ID_UNDEF");
          else
            --ii;
        }
        (*j).second.i = ii;
      }
      else
        throw gError("convert: id is not defined");
    }
}

void RefAlgorithm::DeleteDebug()
{
  if (debug)
  {
    delete debug;
    debug = 0;
  }
}

void RefAlgorithm::SetDebug(const SGML &e)
{
  FileName = e.FindString("DebugFile");
  debug_ = e.FindInt("debug");
  if (algorithm::debug || debug_)
  {
    DeleteDebug();
    string s = FileName;
    if (s.empty())
    {
      string id = e.FindString("id");
      if (id.empty())
        id = "anonymous";
      s = string("debug.algorithm.") + id;
    }
    debug = new ofstream(s.c_str());
  }
}

istream& RefAlgorithm::read(istream& in)
{
  algorithm::init();
  SGML e(in);
  e.compare((*type)[typeid(*this).name()]);
  read(e);
  return in;
}

ostream& RefAlgorithm::write(ostream& out, size_t shift) const
{
  algorithm::init();
  out << PutTab(shift) << "<" << (*type)[typeid(*this).name()] 
    << " class=algorithm";
  WriteAttributes(out, shift + 1);
  out << ">" << endl;
  WriteBody(out, shift + 1);
  out << PutTab(shift) << "</" << (*type)[typeid(*this).name()] << "> " << endl;
  return out;
}

algorithm::Destruct::~Destruct()
{
  if (algorithm::types)
  {
    delete algorithm::types;
    delete algorithm::obj;
    algorithm::types = 0;
    algorithm::obj = 0;
    delete type;
  }
}

void algorithm::init()
{
  if (!types)
  {
    type = new map_string_string;
    PTR_NULL = new PassThrough;
    types = new map_types;
    obj = new map_types;
    RegisterType("PassThrough", PTR_NULL);
    RegisterType("PhaseProperty", new PhaseProperty);
    RegisterType("reaction", new reaction);
    for (size_t i = 0; init_ptr[i]; ++i)
      (*init_ptr[i])();
  }
}

void algorithm::RegisterType(const string& name, RefAlgorithm* ptr)
{
  init();
  (*types)[name] = algorithm(ptr);
  (*type)[typeid(*ptr).name()] = name;
}

void algorithm::find(const string &name)
{
  map_types_i i = obj->find(name);
  if (i == obj->end())
  {
    phase ph;
    try
    {
      ph.find(name);
      create("PhaseProperty");
      cow();
      dynamic_cast<PhaseProperty*>(ptr)->assign(ph);
    }
    catch (gError)
    {
      throw gError(string("algorithm: id ").append(name).append
          (" is not defined"));
    }
  }
  else
    *this = (*i).second;
}

void algorithm::read(const SGML &e)
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
        cout << "algorithm: " << ptr->id << " is already defined - "
          << "set to previous" << endl;
        *this = (*i.first).second;
      }
    }
    else
    {
      cout << "anonymous algorithms are not allowed" << endl;
    }
  }
}

ostream& algorithm::write(ostream& out, size_t shift) const
{
  if (ptr->id.empty() || ptr->print)
  {
    out << PutTab(shift) << "<" << (*type)[typeid(*ptr).name()] 
      << " class=algorithm";
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
    out << PutTab(shift) << "<algorithm class=algorithm IDREF="
      << ptr->id << "></algorithm>"
      << endl;
  }
  return out;
}

void SimpleAlgorithm::read(const SGML &e)
{
  clear();
  parser p(e);
  SGML el;
  p.GetSGML(el);
  el.compare("InputNames");
  parser p2(el);
  while (!p2.eof())
    vs.push_back(p2.GetToken());
  x.resize(vs.size());
}

void SimpleAlgorithm::WriteBody(ostream& out, size_t shift) const
{
  out << PutTab(shift) << "<InputNames> ";
  copy(vs.begin(), vs.end(), ostream_iterator<string>(out, " "));
  out << "</InputNames>" << endl;
}

int SimpleAlgorithm::GetInputIndex(const string& s) const
{
  string str = RemoveSpaces(s);
  int i = SearchString(vs, str);
  if (i == -1)
    throw gError(string("SimpleAlgorithm: no such name - ").append(s));
  return i;
}

const string& SimpleAlgorithm::name(int i) const
{
  if (i >= 0)
    return vs[i];
  else
  {
    static string str("f");
    return str;
  }
}

void InputValue::read(const SGML &e, const algorithm &alg)
{
  iname = e.FindString("name");
  if (iname.empty())
  {
    iy = e.FindInt("NumberOf");
    if (iy == 0)
    {
      if (e.defined("value"))
      {
        x = e.FindDouble("value");
        source = Value;
      }
      else if (e.defined("Value"))
      {
        x = e.FindDouble("Value");
        source = Value;
      }
      else
        source = Convert;
    }
    else
    {
      source = NumberOf;
      --iy;
    }
  }
  else
  {
    iy = -1;
    source = Name;
  }
  if (source == Convert)
  {
    SGML e_ = e;
    parser p(e_);
    SGML el;
    expr.read(p.GetSGML(el));
    e_.body.erase(0, p.npos());
    ind = alg.GetInputIndex(e_.body);
  }
  else
    ind = alg.GetInputIndex(e.body);
}

void InputValue::write(ostream &out, size_t shift, const algorithm &alg, 
    bool once) const
{
  out << PutTab(shift) << "<input";
  if (once)
    out << " once";
  switch (source)
  {
    case Name:
      out << " name=" << iname;
      break;
    case NumberOf:
      out << " NumberOf=" << iy + 1;
      break;
    case Value:
      out << " value=" << x;
  }
  out << "> ";
  if (source == Convert)
  {
    out << endl;
    expr.write(out, shift + 1);
    out << PutTab(shift);
  }
  out << alg.name(ind) << " </input>" << endl;
}

void InputValue::set(const vec_string& vs) const
{
  switch (source)
  {
    case Name:
      iy = SearchString(vs, iname);
      break;
    case Convert:
      expr.SetID(vs);
      break;
  }
}

void InputValue::set(const algorithm &alg, double *reg) const
{
  switch (source)
  {
    case Name:
      if (iy == -1)
        throw gError("InputValue: iy is negative");
    case NumberOf:
      alg.SetVariable(ind, reg[iy]);
      break;
    case Value:
      alg.SetVariable(ind, x);
      break;
    case Convert:
      alg.SetVariable(ind, expr(reg));
      break;
  }
}

void OutputValue::read(const SGML &e, const algorithm &alg)
{
  if (e.defined("noname"))
    noname=true;
  else
    noname = false;
  if (!noname)
    name_ = e.FindString("name");
  else
    name_.erase();
  ind = alg.GetOutputIndex(e.body, N);
}

void OutputValue::write(ostream &out, size_t shift, const algorithm &alg) const
{
  out << PutTab(shift) << "<output";
  if (noname)
    out << " noname";
  else
  {
    if (!name_.empty())
      out << " name=" << name_;
  }
  out << "> " << alg.name(ind) << " </output>" << endl;
}

void OutputValue::names(vec_string &vs, const algorithm &alg) const
{
  if (noname)
  {
    vs.insert(vs.end(), N, "");
  }
  else if (name_.empty())
  {
    alg.OutNames(ind, vs);
  }
  else
  {
    if (N == 1)
      vs.push_back(name_);
    else
    {
      char ch = '1';
      for (size_t i = 0; i < N; ++i)
      {
        vs.push_back(name_ + ch++);
        if (ch == ':')
          ch = 'a';
      }
    }
  }
}

void compute::read(const SGML &e)
{
  clear();
  parser p(e);
  SGML el;
//compatibility 
  if (e.name == "SetPhaseEquilibrium")
  {
    alg.find(e.FindString("PhaseEquilibrium"));
  }
  else
  {
    e.compare("compute");
    alg.read(p.GetSGML(el));
  }
  InputValue i;
  OutputValue o;
  convert c;
  while (!p.eof())
  {
    p.GetSGML(el);
    if (el.name == "input")
    {
      i.read(el, alg);
      if (el.defined("once"))
        once.insert(i);
      else
        in.insert(i);
    }
    else if (el.name == "output")
    {
      o.read(el, alg);
      out.push_back(o);
    }
    else if (el.name == "convert")
    {
      c.read(el);
      con.push_back(c);
    }
    else
      throw gError(string("compute: not known element - ") +
                   el.name);
  }
  v2.resize(con.size());
  size_t n = 0;
  for (vec_output_ci i = out.begin(); i != out.end(); ++i)
    n += (*i).size();
  if (n == 0)
    throw gError("compute: no output");
  v1.resize(n);
  for (vec_output_ci i = out.begin(); i != out.end(); ++i)
    (*i).names(nm_out, alg);
  for (vec_convert_ci i = con.begin(); i != con.end(); ++i)
    nm_con.push_back((*i).name());
  SetOutConvert();
}

ostream& compute::write(ostream& of, size_t shift) const
{
  of << PutTab(shift) << "<compute>" << endl;
  alg.write(of, shift + 1);
  for (set_input_ci i = once.begin(); i != once.end(); ++i)
    (*i).write(of, shift + 1, alg, true);
  for (set_input_ci i = in.begin(); i != in.end(); ++i)
    (*i).write(of, shift + 1, alg);
  for (vec_output_ci i = out.begin(); i != out.end(); ++i)
    (*i).write(of, shift + 1, alg);
  for (vec_convert_ci i = con.begin(); i != con.end(); ++i)
    (*i).write(of, shift + 1);
  of << PutTab(shift) << "</compute>" << endl;
  return of;
}

const double& compute::f(double *reg) const
{
  SetIn(reg);
  if (con.empty())
  {
    out[0].out(alg, &*v1.begin());
    return v1[0];
  }
  else
  {
    SetV1();
    return v2[0] = con[0]();
  }
}

const vec_double& compute::vec(double *reg) const
{
  SetIn(reg);
  SetV1();
  if (con.empty())
    return v1;
  else
  {
    SetV2();
    return v2;
  }
}

int PassThrough::GetInputIndex(const string& s) const
{
  string str = RemoveSpaces(s);
  int i = SearchString(vs, str);
  if (i >= 0)
    return i;
  else
  {
    x.push_back(0.);
    vs.push_back(str);
    return vs.size() - 1;
  }
}

void PhaseProperty::clear()
{
  size_t n = ph.size();
  x.resize(n);
  dmf = 2;
  st.resize(n + 2);
  nm_st.resize(n + 2);
  nm_st2.resize(n + 2);
  nm_st2[0] = nm_st[0] = "T";
  nm_st2[1] = nm_st[1] = "p";
  Tp.SetT(st[0]); // = 1000 from default
  Tp.Setp(st[1]); // = 1 from default
  for (size_t i = 0; i < n; ++i)
  {
    x.SetX(i, st[i + 2] = 0.);
    nm_st[i + 2] = string("x(") + ph.comps()[i].cf() + ")";
    nm_st2[i + 2] = string("x(") + ObjToString(i + 1) + ")";
  }
  vd.clear();
  vd.reserve(5);
  eval = false;
  eval_vint = false;
  vint.resize(ph.NIntVars());
}

void PhaseProperty::read(const SGML &e)
{
  parser p(e);
  SGML el;
  ph.read(p.GetSGML(el));
  clear();
  const string &str = e.FindString("DependentMoleFraction");
  if (str.empty())
    return;
  else if (str == "first")
    dmf = 2;
  else if (str == "last")
    dmf = st.size() - 1;
  else if (str == "no")
    dmf = st.size();
  else
  {
    int ix = SearchFormula(ph.comps(), str);
    if (ix == -1)
    {
      cout << "PhaseProperty: no such component - " << str 
        << " set to default" << endl;
      dmf = 2;
    }
    else
      dmf = ix + 2;
  }
}

void PhaseProperty::WriteAttributes(ostream &out, size_t shift) const
{
  out << endl << PutTab(shift) << "DependentMoleFraction=";
  if (dmf == 2)
    out << "first";
  else if (dmf < st.size() - 2)
    out << dmf - 2;
  else if (dmf == st.size() - 1)
    out << "last";
  else
    out << "no";
}

void PhaseProperty::WriteBody(ostream &out, size_t shift) const
{
  ph.write(out, shift);
}

int PhaseProperty::GetInputIndex(const string& s) const
{
  string str = RemoveSpaces(s);
  int i = SearchState(str);
  if (i != -1 && i != dmf)
    return -(i + 1);
  else
    throw gError(string("PhaseProperty: can not use for input - ")
      .append(str));
}

const string& PhaseProperty::name(int i) const
{
  if (i < 0)
    return nm_st[-(i + 1)];
  else
    return vd[i].name;
}

int PhaseProperty::GetOutputIndex(const string& s, size_t &N) const
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
  vd.back().name = str;
  size_t ipar = str.find("(");
  string prop;
  string name;
  if (ipar == string::npos)
  {
    prop = str;
    name = "";
  }
  else
  {
    prop = str.substr(0, ipar);
    name = str.substr(ipar + 1, str.rfind(")") - ipar - 1);
  }
  if (str == "state(all)")
  {
    vd[i].from = AllState;
    N = vd[i].N = st.size();
  }
  else if (prop == "internal")
  {
    vd[i].from = Internal;
    if (name == "all")
    {
      vd[i].ix = -1;
      N = vd[i].N = ph.NIntVars();
    }
    else
    {
      vd[i].ix = SearchString(ph.NamesIntVars(), name);
      if (vd[i].ix == -1)
        throw gError(string("PhaseProperty: not known name in ").append(str));
      N = vd[i].N = 1;
    }
  }
  else if (prop == "stable")
  {
    vd.back().from = Stability;
    if (name == "all")
    {
      vd[i].ix = -1;
      N = vd[i].N = ph.NIntPhases();
    }
    else
    {
      vd[i].ix = SearchString(ph.NamesIntPhases(), name);
      if (vd[i].ix == -1)
        throw gError(string("PhaseProperty: not known name in ").append(str));
      N = vd[i].N = 1;
    }
  }
  else
  {
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
        N = vd[i].N = ph.size();
      }
      else
      {
        vd[i].ix = SearchFormula(ph.comps(), name);
        if (vd[i].ix == -1)
          throw gError(string("PhaseProperty: not known name in ").append(str));
        N = vd[i].N = 1;
      }
    }
  }
  return i;
}

void PhaseProperty::OutNames(int i, vec_string &vs) const
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
        for (size_t j = 0; j < ph.size(); ++j)
        {
          vs.push_back(vd[i].name);
          vs.back().replace(vs.back().find("all"), 3, ph.comps()[j].cf());
        }
      else
        vs.push_back(vd[i].name);
      return;
    case Internal:
      if (vd[i].ix == -1)
        for (size_t j = 0; j < ph.NIntVars(); ++j)
        {
          vs.push_back(vd[i].name);
          vs.back().replace(vs.back().find("all"), 3, ph.NamesIntVars()[j]);
        }
      else
        vs.push_back(vd[i].name);
      return;
    case Stability:
      if (vd[i].ix == -1)
        for (size_t j = 0; j < ph.NIntPhases(); ++j)
        {
          vs.push_back(vd[i].name);
          vs.back().replace(vs.back().find("all"), 3, ph.NamesIntPhases()[j]);
        }
      else
        vs.push_back(vd[i].name);
      return;
    case AllState:
      vs.insert(vs.end(), nm_st.begin(), nm_st.end());
      return;
  }
}

void PhaseProperty::out(int i, double *x_) const
{
  if (!eval)
  {
    if (dmf < st.size())
    {
      st[dmf] = 1.;
      for (size_t i = 2; i < st.size(); ++i)
        if (i != dmf)
          st[dmf] -= st[i];
    }
    eval = true;
  }
  if (i < 0)
  {
    *x_ = st[-(i + 1)];
    return;
  }
  else
  {
    if (!eval_vint)
    {
      vint = ph.IntVarsEq(Tp, x);
      eval_vint = true;
    }
    switch(vd[i].from)
    {
      case Property:
        *x_ = ph.Z(vd[i].f, vd[i].i, Tp, x, vint);
        return;
      case Partial:
        if (vd[i].ix != -1)
          *x_ = ph.z(vd[i].f, vd[i].i, Tp, x, vint)[vd[i].ix];
        else
        {
          const vec_double &r = ph.z(vd[i].f, vd[i].i, Tp, x, vint);
          copy(r.begin(), r.end(), x_);
        }
        return;
      case Internal:
        if (vd[i].ix != -1)
          *x_ = vint[vd[i].ix];
        else
          copy(vint.begin(), vint.end(), x_);
        return;
      case Stability:
        if (vd[i].ix != -1)
          *x_ = ph.IsStable(vd[i].ix, Tp, x);
        else
          for (size_t j = 0; j < ph.NIntPhases(); ++j)
            *x_++ = ph.IsStable(j, Tp, x);
        return;
      case AllState:
        copy(st.begin(), st.end(), x_);
        return;
    }
  }
}

void reaction::read(const SGML &e)
{
  clear();
  parser p(e);
  SGML el;
  compute alg;
  while (!p.eof())
  {
    p.GetSGML(el);
    alg.read(el);
    reac.push_back(alg);
  }
  OutN = 0;
  for (vec_compute_ci i = reac.begin(); i != reac.end(); ++i)
    OutN += (*i).size();
}

void reaction::WriteBody(ostream& out, size_t shift) const
{
  for (vec_compute_ci i = reac.begin(); i != reac.end(); ++i)
    (*i).write(out, shift + 1);
}

int reaction::GetInputIndex(const string& s) const
{
  string str = RemoveSpaces(s);
  int i = SearchString(vs, str);
  if (i >= 0)
    return i;
  else
  {
    reg.push_back(0.);
    vs.push_back(str);
    set = false;
    return vs.size() - 1;
  }
}

const string& reaction::name(int i) const
{
  if (i >= 0)
    return vs[i];
  else
  {
    static string str("all");
    return str;
  }
}

void reaction::out(int i, double *x) const
{
  if (!set)
  {
    SetInput(vs);
    set = true;
  }
  for (vec_compute_ci i = reac.begin(); i != reac.end(); ++i)
  {
    SetOnce(vs, &*reg.begin());
    const vec_double &r = (*i).vec(&*reg.begin());
    copy(r.begin(), r.end(), x);
    x += r.size();
  }
}

