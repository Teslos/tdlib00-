#include <../td_algo/phase_eq.h> 
#include "opt.h"
#include <output.h> 
#include <ass_sol.h> 

struct pointer
{
  enum {DOUBLE, BOOL, INT, LONG, STRING} from;
  void *ptr;
  pointer() : ptr(0) {}
  
  pointer(double &x)
  {
    from = DOUBLE;
    ptr = &x;
  }
  pointer(bool &x)
  {
    from = BOOL;
    ptr = &x;
  }
  pointer(int &x)
  {
    from = INT;
    ptr = &x;
  }
  pointer(long &x)
  {
    from = LONG;
    ptr = &x;
  }
  pointer(string &x)
  {
    from = STRING;
    ptr = &x;
  }
  void write(ostream &out)
  {
    switch (from)
    {
      case DOUBLE:
        out << *(double*)ptr;
        break;
      case BOOL:
        out << *(bool*)ptr;
        break;
      case INT:
        out << *(int*)ptr;
        break;
      case LONG:
        out << *(long*)ptr;
        break;
      case STRING:
        out << *(string*)ptr;
        break;
    }
  }
  void assign(const string &str)
  {
    switch (from)
    {
      case DOUBLE:
        *(double*)ptr = atof(str.c_str());
        break;
      case BOOL:
        *(bool*)ptr = atoi(str.c_str());
        break;
      case INT:
        *(int*)ptr = atoi(str.c_str());
        break;
      case LONG:
        *(long*)ptr = atoi(str.c_str());
        break;
      case STRING:
        *(string*)ptr = str;
        break;
    }
  }
};

static map<string, pointer> gl;
static bool init_globals = false;

void InitGlobals()
{
  if (!init_globals)
  {
    gl.insert(map<string, pointer>::value_type("R", global::R));
    gl.insert(map<string, pointer>::value_type("DefaultT", global::T));
    gl.insert(map<string, pointer>::value_type("Defaultp", global::p));
    gl.insert(map<string, pointer>::value_type("FirstDerivativeStep", 
          global::step1));
    gl.insert(map<string, pointer>::value_type("SecondDerivativeStep", 
          global::step2));
    gl.insert(map<string, pointer>::value_type("NegativeLog", global::neg_log));
    gl.insert(map<string, pointer>::value_type("EpsForNegLog", global::eps));
    gl.insert(map<string, pointer>::value_type("DebugPhase", phase::debug));
    gl.insert(map<string, pointer>::value_type("GetNuTolerence", GetNu::tol));
    gl.insert(map<string, pointer>::value_type("OutputWidth", 
          OutputFile::width));
    gl.insert(map<string, pointer>::value_type("DebugAlgorithm", 
          algorithm::debug));
    gl.insert(map<string, pointer>::value_type("PETmin", 
          PhaseEquilibrium::Tmin));
    gl.insert(map<string, pointer>::value_type("PETmax", 
          PhaseEquilibrium::Tmax));
    gl.insert(map<string, pointer>::value_type("PETolerence", 
          PhaseEquilibrium::eps));
    gl.insert(map<string, pointer>::value_type("PEpenalty", 
          PhaseEquilibrium::penalty));
    gl.insert(map<string, pointer>::value_type("AssTolerance", 
          AssociatedSolutionBasis::tol));
    gl.insert(map<string, pointer>::value_type("AssFactr", 
          AssSolManyReact::factr));
    gl.insert(map<string, pointer>::value_type("AssPgtol", 
          AssSolManyReact::pgtol));
    gl.insert(map<string, pointer>::value_type("AssM", 
          AssSolManyReact::lbfg_m));
    gl.insert(map<string, pointer>::value_type("PEm", 
          PhaseEquilibrium::GlSet.m));
    gl.insert(map<string, pointer>::value_type("PEiter", 
          PhaseEquilibrium::GlSet.iter));
    gl.insert(map<string, pointer>::value_type("PEfactr", 
        PhaseEquilibrium::GlSet.factr));
    gl.insert(map<string, pointer>::value_type("PEpgtol", 
        PhaseEquilibrium::GlSet.pgtol));
    gl.insert(map<string, pointer>::value_type("PEsmall", 
        PhaseEquilibrium::GlSet.small));
    gl.insert(map<string, pointer>::value_type("PEstep", 
        PhaseEquilibrium::GlSet.step));
    gl.insert(map<string, pointer>::value_type("IntVarStep", 
          RefInternalVariables::step));
    init_globals = true;
  }
}

void ReadGlobals(const SGML &e)
{
  InitGlobals();
  e.compare("globals");
  map<string, pointer>::iterator j;
  for (map_string_string_ci i = e.attr.begin(); i != e.attr.end(); ++i)
    if ((j = gl.find((*i).first)) != gl.end())
      (*j).second.assign((*i).second);
  SGML el;
  parser p(e);
  while (!p.eof())
  {
    p.GetSGML(el);
    if (el.name == "optimizer")
      optimizer::ReadSGML(el);
  }
}

ostream& WriteGlobals(ostream &out, size_t shift)
{
  InitGlobals();
  out << PutTab(shift) << "<globals" << endl;
  map<string, pointer>::iterator i;
  for (i = gl.begin(); i != gl.end(); ++i)
  {
    out << PutTab(shift + 1) << (*i).first << "=";
    (*i).second.write(out);
    out << endl;
  }
  out << PutTab(shift + 1) << ">" << endl;
  optimizer::WriteSGML(out, shift + 1);
  out << PutTab(shift) << "</globals>" << endl;
  return out;
}

