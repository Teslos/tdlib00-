#include <iomanip>
#include <iostream>
#include <lapack.h>

int main()
{
  int n = 3;
  int info;
  int ipiv[3];
  double a[3][3] = {10.,  1.,  5.,
                     1.,  2., -1.,
                     5., -1.,  5.};
  double b[3][3] = { 1.,  0.,  0.,
                     0.,  1.,  0.,
                     0.,  0.,  1.};
  for (int i = 0; i < n; ++i)
  {
    for (int j = 0; j < n; ++j)
      cout << setw(10) << a[j][i];
    cout << endl;
  }
  cout << endl;
  dgetf2(n, n, (double*)a, ipiv, info);
  for (int i = 0; i < n; ++i)
  {
    for (int j = 0; j < n; ++j)
      cout << setw(10) << a[j][i];
    cout << endl;
  }
  cout << "info is " << info << endl;
  dgetrs(n, n, 3, (double*)a, ipiv, (double*)b);
  for (int i = 0; i < n; ++i)
  {
    for (int j = 0; j < n; ++j)
      cout << setw(10) << b[j][i];
    cout << endl;
  }
  return 0;
}
