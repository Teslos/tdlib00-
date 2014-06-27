/*
Copyright (C) 1998, 1999, 2000 Evgenii Rudnyi, rudnyi@comp.chem.msu.su
                                       http://www.chem.msu.su/~rudnyi/

This software is a copyrighted work licensed under the terms, described
in the file "FREE_LICENSE". 
*/

#ifndef __LAPACK
#define __LAPACK

#include <stddef.h>
#include <general.h>

class LapackError : public gError
{
public:
  int info;
  LapackError(const string &s, int ier) 
    : gError(s + string(" - info is ").append(ObjToString(ier))), info(ier) {}
};


extern "C" {
  void dgetrf_(int *m, int *n, double *a, int *lda, int *ipiv, int *info);

  void dlaswp_(int *n, double *a, int *lda, int *k1, int *k2, int *ipiv, 
      int *inx);

  void dtrsm_(const char *side, const char *uplo, const char *transa, const char *diag, int *m, int *n, double *alpha, double *a, int* lda, double *b, int* ldb);

  void dpptri_(const char *uplo, int *n, double *ap, int* info);
  void dpptrf_(const char *uplo, int *n, double *ap, int* info);
}


inline void dgetf2(int m, int n, double *a, int *ipiv, int &info)
{
  if (m && n)
    dgetrf_(&m, &n, a, &m, ipiv, &info);
}

inline void dgetrf(int m, int n, double *a, int *ipiv, int &info)
{
  if (m && n)
    dgetrf_(&m, &n, a, &m, ipiv, &info);
}

inline void dgetrs(int m, int n, int nrns, double *a, int *ipiv, double *b)
{
  if (m && n && nrns)
  {
    int one = 1;
    double one_d = 1.;
    dlaswp_(&nrns, b, &m, &one, &n, ipiv, &one);
    dtrsm_("L", "L", "N", "U", &n, &nrns, &one_d, a, &m, b, &m);
    for (size_t i = n; i < m; ++i)
      for (size_t k = 0; k < nrns; ++k)
        for (size_t j = 0; j < n; ++j)
          b[i + m*k] -= a[i + m*j]*b[j + m*k];
    dtrsm_("L", "U", "N", "N", &n, &nrns, &one_d, a, &m, b, &m);
  }
}

inline void dgesv(int n, int nrns, double *a, int *ipiv, double *b)
{
  if (n && nrns)
  {
    int info;
    dgetf2(n, n, a, ipiv, info);
    if (info)
      throw LapackError("Lapack::dgesv", info);
    dgetrs(n, n, nrns, a, ipiv, b);
  }
}

inline void dpptrin(int n, double* a, int &info)
{
  dpptrf_("U", &n, a, &info);
  if (info)
    return;
  dpptri_("U", &n, a, &info);
}
    
#endif
