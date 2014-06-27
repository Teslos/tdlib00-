/*
Copyright (C) 1998, 1999, 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                       http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#ifndef __CALC_H
#define __CALC_H

#include <general.h>

inline string RemoveSpaces(const string &str)
{
  string out;
  for (string::const_iterator i = str.begin(); i != str.end(); ++i)
    if (!isspace(*i))
      out.append(1, *i);
  return out;
}

class calculator
{
protected:
  enum calculator_key {DOUBLE, DOUBLE_PTR, DOUBLE_PTR_PTR, FUNC1, FUNC2, FUNC3,
                       PLUS, MINUS, MUL, DIV, POW, LP, RP, COMMA, UMIN,
                       END};
  class calculator_entity;
  friend class calculator_entity;
  class calculator_entity
  {
  public:
    calculator_key key;
    union
    {
      double x;
      double *x_ptr;
      double **x_ptr_ptr;
      double (*f1)(double);
      double (*f2)(double, double);
      double (*f3)(double, double, double);
    };

    calculator_entity() {key = END;}
    calculator_entity(double x_) {key = DOUBLE; x = x_;}
    calculator_entity(double *x_ptr_) {key = DOUBLE_PTR; x_ptr = x_ptr_;}
    calculator_entity(double **x_ptr_ptr_)
      {key = DOUBLE_PTR_PTR; x_ptr_ptr = x_ptr_ptr_;}
    calculator_entity(double (*f1_)(double)) {key = FUNC1; f1 = f1_;}
    calculator_entity(double (*f2_)(double, double))
      {key = FUNC2; f2 = f2_;}
    calculator_entity(double (*f3_)(double, double, double))
      {key = FUNC3; f3 = f3_;}
  };

  typedef vector<calculator_entity> vec_entity;
  typedef map<string, calculator_entity, less<string> > map_entity;
  typedef map<char, calculator_key, less<char> > map_key;

  static map_entity *global;
  static map_key *delim_name;
  static bool qty;

  vec_entity tokens;
  string text;
  mutable vec_double stack;

  calculator_entity curr_tok;
  istream* in_ptr;

  void init();

  void get_token() ;
  void expr(bool must);
  void term(bool must);
  void fact(bool must);
  void prim(bool must);

  virtual void AnalizeId(const string &id)
    {throw gError("calculator: ID_UNDEF");}
  virtual string AnalizeSGML(const SGML &e)
    {throw gError("calculator: SGML_UNDEF");}

private:
  class Destruct;
  friend class Destruct;
  class Destruct
  {
  public:
    ~Destruct()
    {
      if (calculator::qty)
      {
        calculator::qty = false;
        delete calculator::global;
        delete calculator::delim_name;
      }
    }
  };
  static Destruct clean;

public:

  calculator() {init();}
  calculator(const calculator& old) {init(); FromString(old.text);}
  calculator(const string &old) {init(); FromString(old);}
  calculator& operator=(const calculator& old)
  {
    if (this != &old) 
      FromString(old.text); 
    return *this;
  }
  ~calculator() {}

  void FromString(const string& old)
  {
    istringstream in(old);
    read(in);
  }
  virtual void clear()
  {
    tokens.clear();
    text.erase();
  }

  istream& read(istream& in)
  {
    clear();
    in_ptr = &in;
    get_token();
    expr(false);
    stack.reserve(10);
    return in;
  }
  ostream& write(ostream& out) const
    {return out << text;}
  const string& expression() const
    {return text;}
  bool empty()
    {return text.empty();}

  double est() const;
};

OVERLOAD_STREAMS(calculator)

#endif

