#include <func_tpx_imp.h>
#include <typeinfo>

size_t bin_coef(size_t n, size_t N)
{
// N!/(N-n)!/n!
  if (n == 0) return 1;
  if ((n > N))
    throw gError("bin_coef - wrong parameters");
  if (N - n < n)
    n = N - n;
  double P = 1;
  for (size_t i = 0; i < n; i++)
    P *= double(N - i)/double(i + 1);
  return size_t(P);
}

const char *FormalismName[] = {"NoChange", "Kohler", "Muggianu", 0};

static map_string_string *type;

RefFuncTpx *FuncTpx::PTR_NULL;
FuncTpx::map_types *FuncTpx::types = 0;
FuncTpx::map_types *FuncTpx::obj = 0;
FuncTpx::Destruct FuncTpx::clean;

void registerTpx()
{
  FuncTpx::RegisterType("IdealMixing", new IdealMixing);
  FuncTpx::RegisterType("Reference", new Reference);
}

void RegisterInteraction();

void (*FuncTpx::init_ptr[])() = {registerTpx, RegisterInteraction, 0};

FuncTpx::Destruct::~Destruct()
{
  if (FuncTpx::types)
  {
    delete FuncTpx::types;
    delete FuncTpx::obj;
    FuncTpx::types = 0;
    FuncTpx::obj = 0;
    delete type;
  }
}

void FuncTpx::init()
{
  if (!types)
  {
    type = new map_string_string;
    PTR_NULL = new NullFuncTpx;
    types = new map_types;
    obj = new map_types;
    RegisterType("NullFuncTpx", PTR_NULL);
    for (size_t i = 0; init_ptr[i]; ++i)
      (*init_ptr[i])();
  }
}

void FuncTpx::RegisterType(const string& name, RefFuncTpx *ptr)
{
  init();
  (*types)[name] = FuncTpx(ptr);
  (*type)[typeid(*ptr).name()] = name;
}

istream& RefFuncTpx::read(istream& in)
{
  FuncTpx::init();
  SGML e(in);
  e.compare((*type)[typeid(*this).name()]);
  read(e);
  return in;
}

ostream& RefFuncTpx::write(ostream& out, size_t shift) const
{
  FuncTpx::init();
  out << PutTab(shift) << "<" << (*type)[typeid(*this).name()] 
    << " class=FuncTpx";
  WriteAttributes(out, 0, shift + 1);
  out << ">" << endl;
  WriteBody(out, 0, shift + 1);
  out << PutTab(shift) << "</" << (*type)[typeid(*this).name()] << "> " << endl;
  return out;
}

void FuncTpx::read(const SGML &e, const vec_formula *f)
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
    ptr->read(e, f);
    if (!(ptr->id = e.FindString("id")).empty())
    {
      pair<map_types_i, bool> i;
      i = obj->insert(map_types::value_type(ptr->id, *this));
      if (!i.second)
      {
        cout << "FuncTpx: " << ptr->id << " is already defined - "
          << "set to anonymous" << endl;
        ptr->id.erase();
      }
    }
  }
}

ostream& FuncTpx::write(ostream& out, const vec_formula *f, 
    size_t shift) const
{
  if (ptr->id.empty() || ptr->print)
  {
    out << PutTab(shift) << "<" << (*type)[typeid(*ptr).name()] 
      << " class=FuncTpx";
    if (!ptr->id.empty())
    {
      out << " id=" << ptr->id;
    }
    ptr->WriteAttributes(out, f, shift + 1);
    out << ">" << endl;
    ptr->WriteBody(out, f, shift + 1);
    out << PutTab(shift) << "</" << (*type)[typeid(*ptr).name()] << "> " << endl;
    ptr->print = false;
  }
  else
  {
    out << PutTab(shift) << "<FuncTpx class=FuncTpx IDREF="
      << ptr->id << "></FuncTpx>"
      << endl;
  }
  return out;
}

void IdealMixing::read(const SGML &e, const vec_formula *f)
{
  clear();
  parser p(e);
  if (f && p.eof())
  {
    vec_int v(f->size());
    set(v);
  }
  else
  {
    int num;
    SGML el;
    while (!p.eof())
    {
      coef cf(1.);
      p.GetToken();
      p.GetSGML(el);
      if (!el.name.empty())
        cf.read(el);
      p.SkipChar('*');
      num = FindFormula(p, f);
      size_t i;
      for (i = 0; i < size(); ++i)
        if (num < vars[i])
          break;
      vars.insert(vars.begin() + i, num);
      vc.insert(vc.begin() + i, cf);
      p.SkipUntil('\n');
    }
  }
}

void IdealMixing::WriteBody(ostream &out, const vec_formula *f, 
    size_t shift) const
{
  for (size_t i = 0; i < size(); ++i)
  {
    out << PutTab(shift) << "+R*T*";
    if (vc[i].x() != 1.)
      out << vc[i] << "*";
    out << "x(";
    WriteFormula(out, i, f);
    out << ")*log(x(";
    WriteFormula(out, i, f);
    out << "))" << endl;
  }
}

double IdealMixing::Z(function f, const StateTp &Tp, const StateX &x) const
{
  double sum = 0.;
  if (f == ::S || f == ::G)
  {
    for (size_t i = 0; i < size(); ++i)
    {
      if (global::neg_log || x[vars[i]] != 0.)
        sum += vc[i].x()*x[vars[i]]*log_z(x[vars[i]]);
    }
    if (f == ::G)
      sum *=global::R*Tp.T();
    else
      sum *=-global::R;
  }
  return sum;
}

void IdealMixing::dZdx(function f, const StateTp &Tp, const StateX &x, 
    vec_double &res) const
{
  if (f == ::S || f == ::G)
  {
    double coef;
    if (f == ::G)
      coef = global::R*Tp.T();
    else
      coef = -global::R;
    for (size_t j = 0; j < size(); ++j)
        res[vars[j]] += vc[j].x()*coef*(log_z(x[vars[j]]) + 1.);
  } 
}

void IdealMixing::d2Zdx2(function f, const StateTp &Tp, const StateX &x,
    SymMatrix &res) const
{
  if (f == ::S || f == ::G)
  {
    double coef;
    if (f == ::G)
      coef = global::R*Tp.T();
    else
      coef = -global::R;
    for (size_t j = 0; j < size(); ++j)
        res(vars[j], vars[j]) += vc[j].x()*coef*dlog_z(x[vars[j]]);
  } 
}

void Reference::read(const SGML &e, const vec_formula *f)
{
  clear();
  parser p(e);
  SGML el;
  species s;
  int num;
  int cur = 0;
  while (!p.eof())
  {
    p.GetSGML(el);
    num = -2;
    if (el.name == "func_x")
    {
      parser p2(el);
      p2.SkipChar('+');
      if (p2.eof())
        num = -1;
      else
        num = FindFormula(p2, f);
      p.GetSGML(el);
    }
    if (el.name == "species")
    {
      s.read(el);
      if (num == -2)
      {
        if (f)
          num = FindFormula(s.cf(), f);
        else
          num = cur++;
      }
      size_t i;
      for (i = 0; i < size(); ++i)
        if (num < vars[i])
          break;
      vars.insert(vars.begin() + i, num);
      sps.insert(sps.begin() + i, s);
    }
    if (el.name.empty())
      break;
  }
}

void Reference::WriteBody(ostream &out, const vec_formula *f, 
    size_t shift) const
{
  for (size_t i = 0; i < size(); ++i)
  {
    out << PutTab(shift) << "<func_x> ";
    if (vars[i] == -1)
      out << "+";
    else
    {
      out << "+x(";
      WriteFormula(out, i, f);
      out << ")*";
    }
    out << " </func_x>" << endl;
    sps[i].write(out, shift);
  }
}

double Reference::Z(function f, const StateTp &Tp, const StateX &x) const
{
  double sum = 0.;
  for (size_t i = 0; i < size(); ++i)
    if (vars[i] == -1)
      sum += sps[i].Z(f, Tp);
    else
      sum += x[vars[i]]*sps[i].Z(f, Tp);
  return sum;
}

void Reference::dZdx(function f, const StateTp &Tp, const StateX &x, 
    vec_double &res) const
{
  for (size_t j = 0; j < size(); ++j)
    if (vars[j] >= 0)
      res[vars[j]] += sps[j].Z(f, Tp);
}

void Reference::z(function f, const StateTp &Tp, vec_double &res) const
{
  double Sum = 0.;
  for (size_t i = 0; i < size(); ++i)
    if (vars[i] == -1)
      Sum += sps[i].Z(f, Tp);
  for (size_t i = 0; i < size(); ++i)
      res[i] += Sum;
  for (size_t i = 0; i < size(); ++i)
    if (vars[i] >= 0)
      res[vars[i]] += sps[i].Z(f, Tp);
}

