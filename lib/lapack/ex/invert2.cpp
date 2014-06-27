#include <iomanip>
#include <iostream>
#include <lapack.h>

int main()
{
  int n = 3;
  int info;
  double a[6] = {10.,  1.,  2., 5., -1., 5.};
  int k = 0;
  for (int i = 0; i < n; ++i)
  {
    for (int j = 0; j <= i; ++j)
      cout << setw(10) << a[k++];
    cout << endl;
  }
  cout << endl;
  dpptrin(n, a, info);
  cout << "info is " << info << endl;
  k = 0;
  for (int i = 0; i < n; ++i)
  {
    for (int j = 0; j <= i; ++j)
      cout << setw(10) << a[k++];
    cout << endl;
  }
  cout << endl;
  return 0;
}
