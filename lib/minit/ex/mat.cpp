#include <iostream>
#include <minit.h>

int main()
{
  Cmatrix A;
  A.resize(3, 6);
  A[1][1] = 2.;
  A[0][4] = 3.;
  A[2][3] = 4.;
  
  for (size_t i = 0; i < 3; ++i)
  {
    for (size_t j = 0; j < 6; ++j)
      cout << A[i][j] << " ";
    cout << endl;
  }
  cout << endl;
  for (size_t i = 0; i < 3; ++i)
  {
    for (size_t j = 0; j < 6; ++j)
      cout << A(i, j) << " ";
    cout << endl;
  }
  cout << endl;
  double b[3][6];
  b[1][1] = 2.;
  b[0][4] = 3.;
  b[2][3] = 4.;
  
  for (size_t i = 0; i < 3; ++i)
  {
    for (size_t j = 0; j < 6; ++j)
      cout << b[i][j] << " ";
    cout << endl;
  }
  cout << endl;
}
