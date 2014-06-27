/*
Copyright (C) 1998, 1999, 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                       http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#include "calc.h"

calculator::map_entity *calculator::global;
calculator::map_key *calculator::delim_name;
bool calculator::qty = false;
calculator::Destruct calculator::clean;

void calculator::init()
{
  if (!qty)
  {
    qty = true;
    global = new map_entity;
    typedef map_entity::value_type value;
    (*global).insert(value("acos", calculator_entity(acos)));
    (*global).insert(value("abs", calculator_entity(fabs)));
    (*global).insert(value("asin", calculator_entity(asin)));
    (*global).insert(value("atan", calculator_entity(atan)));
    (*global).insert(value("cos", calculator_entity(cos)));
    (*global).insert(value("exp", calculator_entity(exp)));
    (*global).insert(value("log10", calculator_entity(log10)));
    (*global).insert(value("log", calculator_entity(log)));
    (*global).insert(value("sin", calculator_entity(sin)));
    (*global).insert(value("sqrt", calculator_entity(sqrt)));
    (*global).insert(value("tan", calculator_entity(tan)));
    (*global).insert(value("R", calculator_entity(global::R)));

    delim_name = new map_key;
    (*delim_name)['+'] = PLUS;
    (*delim_name)['-'] = MINUS;
    (*delim_name)['*'] = MUL;
    (*delim_name)['/'] = DIV;
    (*delim_name)['^'] = POW;
    (*delim_name)['('] = LP;
    (*delim_name)[')'] = RP;
    (*delim_name)[','] = COMMA;
  }
}

void calculator::get_token()
{
  *in_ptr >> ws;
  char t;
  map_key::iterator i;
  in_ptr->get(t);
  if (!*in_ptr)
  {
    curr_tok.key = END;
    return;
  }
  if ((i = delim_name->find(t)) != delim_name->end())
  {
    text += t;
    curr_tok.key = (*i).second;
    return;
  }
  else if (t == '<')
  {
    in_ptr->putback(t);
    SGML e(*in_ptr);
    if (e.name == "str")
    {
      e.body = RemoveSpaces(e.body);
      text += string("<str> ") + e.body + " </str>";
      AnalizeId(e.body);
    }
    else
      text += AnalizeSGML(e);
    return;
  }
  else if (isdigit(t) || t == '.')
  {
    in_ptr->putback(t);
    *in_ptr >> curr_tok.x;
    curr_tok.key = DOUBLE;
//    char buf[30];
//    text += gcvt(curr_tok.x, 7, buf);
    text += ObjToString(curr_tok.x, 7);
    return;
  }
  else if (t == ';')
  {
    curr_tok.key = END;
    return;
  }
  else if (isalnum(t) || t == '_')
  {
    string id(1, t);
    while (in_ptr->get(t), *in_ptr && (isalnum(t) || t == '_'))
      id += t;
    in_ptr->putback(t);
    map_entity::iterator i;
    text += id;
    if ((i = global->find(id)) != global->end())
      curr_tok = (*i).second;
    else
      AnalizeId(id);
  }
  else
   throw gError("calculator: WRTOK");
}

void calculator::expr(bool must)
{
  calculator_entity tmp;
  term(must);
  for(;;)
    switch ((tmp = curr_tok).key) {
      case PLUS: case MINUS:
        get_token();
        term(true);
        tokens.push_back(tmp);
        break;
			case END: case RP: case COMMA:
        return;
      default:
        throw gError("calculator: NOTERM");
    }
}

void calculator::term(bool must)
{
  calculator_entity tmp;
  fact(must);
  for(;;)
    switch ((tmp = curr_tok).key) {
      case MUL: case DIV:
        get_token();
        fact(true);
        tokens.push_back(tmp);
        break;
			case END: case PLUS: case MINUS: case RP: case COMMA:
  return;
      default:
        throw gError("calculator: NOFACT");
   }
}

void calculator::fact(bool must)
{
  calculator_entity tmp;
  prim(must);
  for(;;)
    switch ((tmp = curr_tok).key) {
      case POW:
        get_token();
        prim(true);
        tokens.push_back(tmp);
        break;
      case END: case PLUS: case MINUS: case MUL: case DIV: case RP: case COMMA:
        return;
      default:
        throw gError("calculator: NOPRIM");
    }
}

void calculator::prim(bool must)
{
  calculator_entity tmp;
  switch (curr_tok.key)
  {
    case DOUBLE:
    case DOUBLE_PTR:
    case DOUBLE_PTR_PTR:
      tokens.push_back(curr_tok);
      get_token();
      return;
    case LP:
      get_token();
      expr(true);
      if (curr_tok.key != RP) throw gError("calculator: NORP");
      get_token();
      return;
    case FUNC1:
      tmp = curr_tok;
      get_token();
      if (curr_tok.key != LP) throw gError("calculator: NOLP");
      get_token();
      expr(true);
      if (curr_tok.key != RP) throw gError("calculator: NORP");
      get_token();
      tokens.push_back(tmp);
      return;
    case FUNC2:
      tmp = curr_tok;
      get_token();
      if (curr_tok.key != LP) throw gError("calculator: NOLP");
      get_token();
      expr(true);
      if (curr_tok.key != COMMA) throw gError("calculator: NOCOMMA");
      get_token();
      expr(true);
      if (curr_tok.key != RP) throw gError("calculator: NORP");
      get_token();
      tokens.push_back(tmp);
      return;
    case FUNC3:
      tmp = curr_tok;
      get_token();
      if (curr_tok.key != LP) throw gError("calculator: NOLP");
      get_token();
      expr(true);
      if (curr_tok.key != COMMA) throw gError("calculator: NOCOMMA");
      get_token();
      expr(true);
      if (curr_tok.key != COMMA) throw gError("calculator: NOCOMMA");
      get_token();
      expr(true);
      if (curr_tok.key != RP) throw gError("calculator: NORP");
      get_token();
      tokens.push_back(tmp);
      return;
    case PLUS:
      get_token();
      prim(true);
      return;
    case MINUS:
      tmp.key = UMIN;
      get_token();
      prim(true);
      tokens.push_back(tmp);
      return;
    case RP:
      throw gError("calculator: NOLP");
    case END:
      if (must) throw gError("calculator: NOPRIM");
      return;
    default:
      throw gError("calculator: NOPRIM");
  }
}

double calculator::est() const
{
  if (tokens.empty()) return 0.;
  double arg2, arg3;
  vec_entity::const_iterator i;
  for (i = tokens.begin(); i != tokens.end(); i++)
  {
    switch ((*i).key)
    {
      case DOUBLE:
        stack.push_back((*i).x);
        break;
      case DOUBLE_PTR:
        stack.push_back(*(*i).x_ptr);
        break;
      case DOUBLE_PTR_PTR:
        stack.push_back(**(*i).x_ptr_ptr);
        break;
      case FUNC1:
        stack.back() = (*(*i).f1)(stack.back());
        break;
      case FUNC2:
        arg2 = stack.back();
        stack.pop_back();
        stack.back() = (*(*i).f2)(stack.back(), arg2);
        break;
      case FUNC3:
        arg3 = stack.back();
        stack.pop_back();
        arg2 = stack.back();
        stack.pop_back();
        stack.back() = (*(*i).f3)(stack.back(), arg2, arg3);
        break;
      case PLUS:
        arg2 = stack.back();
        stack.pop_back();
        stack.back() += arg2;
        break;
      case MINUS:
        arg2 = stack.back();
        stack.pop_back();
        stack.back() -= arg2;
        break;
      case MUL:
        arg2 = stack.back();
        stack.pop_back();
        stack.back() *= arg2;
        break;
      case DIV:
        arg2 = stack.back();
        stack.pop_back();
        stack.back() /= arg2;
        break;
      case POW:
        arg2 = stack.back();
        stack.pop_back();
        stack.back() = pow(stack.back(), arg2);
        break;
      case UMIN:
        stack.back() = -stack.back();
        break;
      default:
        throw gError("calculator: run-time WRTOK");
    }
  }
  if (stack.size() != 1)
    throw gError("calculator: run-time WRONG_TOKENS");
  arg2 = stack.back();
  stack.pop_back();
  return arg2;
}

