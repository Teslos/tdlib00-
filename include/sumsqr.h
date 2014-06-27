/*
Copyright (C) 1994 - 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#ifndef _SUMSQR_H
#define _SUMSQR_H

#include <list>
#include <fstream>
#include <math.h>
#include <toms.h>
#include "data.h"
#include "residual.h"

#include "f2c.h"
#include "f2clib.h"

void init_before_reading();
void init_after_reading();
void post_analysis();
extern char* LICENSE;
const int prcsn = 5;
const int wdth = 12;

const int nvar = 30;
extern double sr2t[nvar];
extern int Nst[nvar];
extern double gat[nvar];
extern int Ngat[nvar];
extern double gbt[nvar];
extern int Ngbt[nvar];

class series : public data
{
public:
  enum var_des {fixed, same, own};

  string ID;
  residual f;
  vec_string vs;
  vec_double vx;
  bool hide;

  double sumi;
  double suma;
  double sumb;

  double sr2, sr2old;
  double adjust;
  int fl_sr2;
  int isr2;
  double ga, gaold;
  int fl_ga;
  int igat;
  double gb, gbold;
  int fl_gb;
  int igbt;

  double SS;

  series()
    {init();}

  void init();
  int operator==(const series& second) const
    {return ID == second.ID;}
  int operator<(const series& second)
    {return ID < second.ID;}
};

extern list<series> ser;
extern list<series>::iterator iser;

enum fl_print {results, big_iter, all_iter};
enum fl_series {no_series, used_series, all_series};
enum fl_start {variances, parameters};
extern fl_print print;
extern fl_series print_series;
extern fl_start start;
extern ostream *out_file;
extern ofstream off;
extern int iss;
extern bool fixmodel;
extern bool write_data_file;
extern bool write_model_file;
extern bool write_e2;
extern int max_iter;

extern string new_data;
extern string file_name;

void print_ss_par_1(double *x, size_t n);
void print_variances(ostream &of, bool print_always_gb = false);
void print_all_points();

void anal_series();
void set_sums_1(series& t);
void set_sums(series& t, double* f, int& i);
void ss(double* x, int m, int n, double* f);

void eval_L(double& L, double& SS);
double gat_solve(double x);
double gbt_solve(double x);
void eval_var_comp();

void process_ini(const string &name, int &arg, bool write = true);
bool open_streams(int agrc, char* argv[], bool ini, int &arg);
bool create_output();
void read_series(istream& df, const string &name);
void read_setfile(istream& sf);

extern vec_string xnm;

#endif
