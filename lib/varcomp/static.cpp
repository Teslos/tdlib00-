/*
Copyright (C) 1994 - 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#include "sumsqr.h"
#include <toms.h>

list<series> ser;
list<series>::iterator iser;

fl_print print = results;
fl_series print_series = no_series;
fl_start start = variances;
int iss = 0;
ostream *out_file;
ofstream off;
bool fixmodel = false;
bool write_data_file = false;
bool write_model_file = false;
bool write_e2 = false;
string new_data;
string file_name;
int max_iter = 5;

double sr2t[nvar];
int Nst[nvar];
double gat[nvar];
int Ngat[nvar];
double gbt[nvar];
int Ngbt[nvar];

vec_string xnm;

