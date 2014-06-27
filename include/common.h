/*
Copyright (C) 1994 - 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#ifndef _COMMON_H
#define _COMMON_H

#include <iostream>
#include <cstdlib>
#include <cstring>

using namespace std;

void warning(const char* txt);

class buffer
{
  char *buf;
  char *cur;
  size_t sz;

public:
  buffer(size_t size = 10240)
    {buf = new char[sz = size]; cur = buf; if (!buf) sz = 0;}
  ~buffer()
    {delete buf;}
  int free()
    {return sz - busy();}
  int busy()
    {return cur - buf;}
  void read(void* ptr)
    {memcpy(ptr, buf, busy()); cur = buf;}
  bool write(void* ptr, size_t n = 1);
  void clear()
    {cur = buf;}
  char* ptr()
    {return buf;}
};

extern buffer Buffer;

#endif
