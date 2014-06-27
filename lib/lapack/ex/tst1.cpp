#include <iomanip>
#include <iostream>
#include <lapack.h>

int main()
{
  int n = 3;
  int info;
  int ipiv[3];
  double a[3][3] = {10., -3.,  5.,
                    -7.,  2., -1.,
                     0.,  6.,  5.};
  double b[3] = {7., 4., 6.};
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
  dgetrs(n, n, 1, (double*)a, ipiv, b);
  cout << "solution " << setw(10) << b[0]
                      << setw(10) << b[1]
                      << setw(10) << b[2] << endl << endl;
  cout << "Reference results from IBM PC/AT" << endl;
  cout << "number of accurate digits = 5" << endl;
  cout << " solution    0.000000   -1.000000    1.000000" << endl;
  return 0;
}
