/*
Copyright (C) 1994 - 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                   http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#include <math.h>
#include "common.h"

buffer Buffer(1024*30);

void warning(const char* txt)
{
  cout << endl << "***WARNING " << txt << endl;
}

bool buffer::write(void* ptr, size_t n)
{
  if (free() < n)
  {
    warning("Buffer is full");
    return false;
  }
  else
  {
    memcpy(cur, ptr, n);
    cur += n;
    return true;
  }
}

