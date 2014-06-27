#include <iomanip>
#include <iostream>
#include <math.h>
#include "toms.h"

void f(long *n, double *x, double *fv, long *iflag)
{
  if (*iflag == 1)
    fv[0] = x[0]*x[1] - pow(x[1], 3) - 1.0;
  else if (*iflag == 2)
    fv[1] = pow(x[0], 2)*x[1] + x[1] - 5.0;
}

int main()
{
  int n = 2;
  double x[2];
  double fv[2];
  x[0] = 1.;
  x[1] = 3.;
//  cout << "guesses (two numbers)? ";
//  cin >> x[0] >> x [1];
  try
  {
    brent1(makeFunctor((VOF*)0, f), n, x, fv);
  }
  catch (err_brent1 &ex)
  {
    cout << ex.message << " " << ex.ier << endl;
  }
  cout << "brent1 results" << endl;
  cout << "solution  = " << setw(10) << x[0] << setw(10) << x[1] << endl;
  cout << "functions = " << setw(10) << fv[0] << setw(10) << fv[1] << endl;
  cout << "reference results from IBM PC/AT" << endl;
  cout << "estimate of solution   2.000000    1.000000" << endl;
  cout << "values of nonlinear functions   0.0000e+00  0.0000e+00" << endl;
  return 0;
}

