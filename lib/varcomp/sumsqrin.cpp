/*
Copyright (C) 1994 - 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#include <math.h>
#include "sumsqr.h"
#include "post_an.h"

const char* HELP ="\n\
Evgenii Rudnyi 1994 - 2000 (c) All rights reserved \n\
processing experimental series with systematic errors \n\n\
    program [-options] [-m] models -v values -d datafiles -s setfiles \n\n\
      models - list of .mod files with models, algorithms and residuals \n\
      values - list of .par files with initial estimates \n\
      datafiles - list of .dat files with experimetnal vales \n\
      setfiles - list of .set files with hypotheses \n\
           options \n\
    -i  initialization file [.ini] (assess.ini is always processed) \n\
    -o  output file (default - con) \n\
    -w  new datafile \n\
    -t  write file for tilt-shift plot \n\
    -h  help \n\
    -l  license \n\
    -f  fix all unknowns \n\
    -p  print results only (no iterations results) (default) \n\
    -p1 print the results of big iterations \n\
    -p2 print all the intermediate resutls \n\
    -e  listing without experimental values (default) \n\
    -e1 print used series only \n\
    -e2 print all the series \n\
    -b  start with variance estimates (default) \n\
    -b1 start with parameter estimates \n";

static const char* WRONGOPTION = "\nNo such option \n\n";
static const char* CANTCREATE = "\nCan not create a file - ";
static const char* NOEQ = "\nCan not find appropriate equation - ";

static enum {model, values, data, hypotheses} from;
static bool hide_series = false;
static string outfile, writefile;

void readfile(const string &nm)
{
  if (nm.empty())
    return;
  string name = nm;
  switch (from)
  {
    case model:
    {
      name += ".mod";
      ifstream in(name.c_str(), ios::in);
      if (!in)
      {
        in.open(nm.c_str(), ios::in);
        if (!in)
        {
          cout << "no file - " << name << " - skipped" << endl;
          return;
        }
      }
      ReadModelFile(in, name, vec);
      return;
    }
    case values:
    {
      name += ".par";
      ifstream in(name.c_str(), ios::in);
      if (!in)
      {
        in.open(nm.c_str(), ios::in);
        if (!in)
        {
          cout << "no file - " << name << " - skipped" << endl;
          return;
        }
      }
      coef::InitialValues(in, name);
      return;
    }
    case data:
    {
      name += ".dat";
      ifstream in(name.c_str(), ios::in);
      if (!in)
      {
        in.open(nm.c_str(), ios::in);
        if (!in)
        {
          cout << "no file - " << name << " - skipped" << endl;
          return;
        }
      }
      read_series(in, name);
      return;
    }
    case hypotheses:
    {
      if (!hide_series)
      {
        hide_series = true;
        for (iser = ser.begin(); iser != ser.end(); ++iser)
          (*iser).hide = true;
      }
      name += ".set";
      ifstream in(name.c_str(), ios::in);
      if (!in)
      {
        in.open(nm.c_str(), ios::in);
        if (!in)
        {
          cout << "no file - " << name << " - skipped" << endl;
          return;
        }
      }
      read_setfile(in);
      return;
    }
  }
}

void process_ini(const string &name, int &arg, bool write)
{
  string nm = name + ".ini";
  ifstream in(nm.c_str(), ios::in);
  if (!in && name != "assess")
  {
    in.open(name.c_str(), ios::in);
    if (!in)
    {
      if (write)
        cout << "no file - " << nm << " - skipped" << endl;
      return;
    }
  }
  string str;
  vec_string tok;
  while (in)
  {
    str = "";
    in >> str;
    if (str.empty())
      break;
    else
      tok.push_back(str);
  }
  arg += tok.size();
  vector<char*> argv(tok.size() + 1);
  for (size_t i = 0; i < tok.size(); ++i)
    argv[i + 1] = const_cast<char*>(tok[i].c_str());
  open_streams(tok.size() + 1, &*argv.begin(), false, arg);
}

bool open_streams(int argc, char* argv[], bool ini, int &arg)
{
  from = model;
  string name;
  for (int i = 1; i < argc; i++)
  {
    if (argv[i][0] == '-' || argv[i][0] == '/')
    {
      switch (tolower(argv[i][1]))
      {
        case 'i' :
          name = argv[i] + 2;
          if (name.empty() && ++i < argc)
            name = argv[i];
          if (ini)
            process_ini(name, arg);
          break;
        case 'm' :
          write_model_file = true;
          from = model;
          name = argv[i] + 2;
          readfile(name);
          break;
        case 'v' :
          from = values;
          name = argv[i] + 2;
          readfile(name);
          break;
        case 'd' :
          from = data;
          name = argv[i] + 2;
          readfile(name);
          break;
        case 's' :
          from = hypotheses;
          name = argv[i] + 2;
          readfile(name);
          break;
        case 'o' :
          outfile = argv[i] + 2;
          if (outfile.empty() && ++i < argc)
          {
            outfile = argv[i];
          }
          if (outfile[outfile.size() - 1] == '.')
            outfile.erase(outfile.size() - 1);
          if (outfile.empty())
            cout << "null output file" << endl;
          break;
        case 'w' :
          writefile = argv[i] + 2;
          if (writefile.empty() && ++i < argc)
          {
            writefile = argv[i];
          }
          if (writefile[writefile.size() - 1] == '.')
            writefile.erase(writefile.size() - 1);
          if (writefile.empty())
          {
            cout << "null file for data output" << endl;
            return false;
          }
          break;
        case 't' :
          write_e2 = true;
          break;
        case 'h' : case '?' :
          cout << HELP << endl;
          return false;
        case 'l' :
          cout << LICENSE; 
          return false;
        case 'f' :
          fixmodel = true; 
          break;
        case 'p' :
          switch (argv[i][2])
          {
            case '\0': case '0':
              print = results; 
              break;
            case '1':
              print = big_iter; 
              break;
            case '2':
              print = all_iter; 
              break;
          }
          break;
        case 'e' :
          switch (argv[i][2])
          {
            case '\0': case '0':
              print_series = no_series; 
              break;
            case '1':
              print_series = used_series; 
              break;
            case '2':
              print_series = all_series; 
              break;
          }
          break;
        case 'b' :
          switch (argv[i][2])
          {
            case '\0': case '0':
              start = variances; 
              break;
            case '1':
              start = parameters; 
              break;
          }
          break;
        default:
          cout << WRONGOPTION << HELP; 
          return false;
      }
    }
    else
    {
      name = argv[i];
      readfile(name);
    }
  }
  return true;
}

bool create_output()
{
  coef::EstimateComputed();
  file_name = outfile;
  if (!outfile.empty())
    outfile += ".lst";
  if (!ser.empty())
  {
    if (outfile.empty())
      out_file = &cout;
    else
    {
      off.open(outfile.c_str());
      off.precision(prcsn);
      off.setf(ios::showpoint);
      if (!off)
      {
        cout << CANTCREATE << outfile << endl;
        return false;
      }
      out_file = &off;
    }
    if (!writefile.empty())
    {
      ofstream out((writefile + ".dat").c_str());
      for (iser = ser.begin(); iser != ser.end(); ++iser)
      {
        if ((*iser).hide)
          out << '*';
        out << (*iser).ID << ',' << endl;
        out << (*iser).f.id() << " ";
        for (size_t i = 0; i < (*iser).vs.size(); ++i)
        {
          out << "<var name=" << (*iser).vs[i]
            << " value=" << (*iser).vx[i] << "></var>";
          if (i != (*iser).vs.size() - 1)
            out << endl;
        }
        out << ',' << endl;
        out << (*iser);
        out << endl;
      }
      ofstream out2((writefile + ".axm").c_str());
      errors1(out2, true, true);
      ofstream out3((writefile + ".set").c_str());
      fl_series old = print_series;
      print_series = used_series;
      print_variances(out3, true);
      print_series = old;
    }
  }
  return true;
}

void read_series(istream& df, const string &name)
{
  char ch;
  vec_convert vc;
  series t;
  while (df)
  {
    t.init();
    vc.clear();
    df >> ws;
    if (df.peek() == '*')
    {
      t.hide = true;
      df.ignore();
    }
    if (hide_series)
      t.hide = true;
    df >> ws;
    ch = get_token(df, t.ID);
    if (isspace(ch))
      ch = skip_until(df);
    if (ch != ',')
      continue;
    df >> ws;
    string str;
    ch = get_token(df, str);
    try
    {
      t.f = FindResidual(str);
    }
    catch (gError) 
    {
      cout << NOEQ << t.ID << " - " << str << " -  series ignored" << endl;
      while (ch = df.get(), ch != ';' &&
                            ch != EOF);
      continue;
    }
    parser p(df, name);
    SGML el;
    if (isspace(ch))
      while (df >> ws, ch = df.peek(), ch == '<')
      {
        p.GetSGML(el);
        el.compare("var");
        str = el.FindString("name");
        if (!str.empty())
        {
          t.vs.push_back(str);
          t.vx.push_back(el.FindDouble("value"));
        }
      }
    if (ch != ',')
      ch = skip_until(df);
    if (ch == ',')
      df.ignore();
    else
      continue;
    while (df >> ws, ch = df.peek(), ch == '<')
    {
      vc.push_back(convert());
      vc.back().read(p.GetSGML(el));
    }
    if (ch == ',')
      df.ignore();
    t.read(df, 0, t.f.NameOfX(), t.f.ScaleOfX());
    if (!t)
      continue;
    if (vc.empty())
      ser.push_back(t);
    else
    {
      series t1;
      t1.hide = t.hide;
      t1.ID = t.ID;
      t1.f = t.f;
      t1.vs = t.vs;
      t1.vx = t.vx;
      t1.set(t.Ntot(), vc.size());
      for (size_t i = 0; i < vc.size(); ++i)
      {
        if (vc[i].name().empty())
          throw gError("read_series: no name while converting");
        else
          t1.name(i) = vc[i].name();
      }
      t1.scale() = t.scale();
      t1.NOfX() = SearchString(t1.names(), t1.f.NameOfX());
      for (size_t i = 0; i < vc.size(); ++i)
        vc[i].SetID(t.names());
      for (size_t i = 0; i < t.Ntot(); ++i)
      {
        t1.atr(i) = t.atr(i);
        for (size_t j = 0; j < vc.size(); ++j)
          t1(i, j) = vc[j](t(i));
      }
      t1.set_av();
      ser.push_back(t1);
    }
    ser.back().f.SetInput(ser.back().names());
    ser.back().f.SetOnceInput(ser.back().vs);
  }
}

void read_var_(istream& sf, int& fl, double& var, double* same_var, int& i)
{
  i = 0;
  double tmp = HUGE_VAL;
  sf >> ws;
  if (sf.peek() == '*')
  {
    fl = series::fixed;
    sf.ignore();
  }
  else if (sf.peek() == '#')
  {
    fl = series::same;
    sf.ignore();
    if (isdigit(sf.peek()))
    {
      sf >> i;
      if (i >= nvar) i = 0;
      sf.clear();
    }
  }
  else if (sf.peek() == '%')
  {
    fl = series::own;
    sf.ignore();
  }
  sf >> tmp;
  if (tmp != HUGE_VAL)
  {
    var = tmp*tmp;
    if (fl == series::same)
      same_var[i] = var;
  }
  sf.clear();
}

char read_var(istream& sf, int& fl, double& var, double* same_var, int& i)
{
  read_var_(sf, fl, var, same_var, i);
  char ch = skip_until(sf);
  return ch;
}

char read_var_1(istream& sf, int& fl, double& var, double* same_var,
                                                int& i, double& adjust)
{
  read_var_(sf, fl, var, same_var, i);
  sf >> ws;
  if (sf.peek() == '@')
  {
    sf.ignore();
    double tmp = HUGE_VAL;
    sf >> tmp;
    if (tmp != HUGE_VAL)
      adjust = tmp*tmp;
    sf.clear();
  }
  char ch = skip_until(sf);
  return ch;
}

void read_setfile(istream& sf)
{
  string word;
  bool hide;
  char ch;
  series* tptr;
  while (sf >> ws, sf.peek() != EOF)
  {
    hide = false;
    if (sf.peek() == '*')
    {
      hide = true;
      sf.ignore();
      sf >> ws;
    }
    ch = get_token(sf, word);
    tptr = NULL;
    iser = ser.begin();
    while (iser != ser.end())
    {
      if (word == (*iser).ID)
      {
        tptr = &(*iser);
        tptr->hide = hide;
        break;
      }
      iser++;
    }
    if (!tptr)
    {
      cout << "series " << word << " not found - ignored" << endl;
      if (ch != ';')
        while (ch = sf.get(), ch != ';' &&
                              ch != EOF);
    }
    else
    {
      if (isspace(ch))
        ch = skip_until(sf);
      if (ch != ',')
        continue;
      if (read_var_1(sf, tptr->fl_sr2, tptr->sr2, sr2t, tptr->isr2,
                                                    tptr->adjust) != ',')
        continue;
      if (read_var(sf, tptr->fl_ga, tptr->ga, gat, tptr->igat) != ',')
        continue;
      if (read_var(sf, tptr->fl_gb, tptr->gb, gbt, tptr->igbt) != ',')
        continue;
      while (ch = sf.get(), ch != ';' &&
                            ch != EOF);
    }
  }
}

