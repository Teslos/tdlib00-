#include <iomanip.h>
#include <iostream.h>
#include <math.h>
#include "toms.h"

double f(double *x, long iflag)
{
  double r;
  if (iflag == 1)
    r = x[0]*x[1] - pow(x[1], 3) - 1.0;
  else if (iflag == 2)
    r = pow(x[0], 2)*x[1] + x[1] - 5.0;
  cout << setw(12) << x[0]
       << setw(12) << x[1]
       << setw(6) << iflag
       << setw(12) << r
       << endl;
  return r;
}

int main()
{
  int n = 2;
  double x[2];
  double h[2];
  double fv[2];
  x[0] = 1.;
  x[1] = 0.;
//  cout << "guesses? ";
//  cin >> x[0] >> x [1];
//  cout << "boxes? ";
//  cin >> h[0] >> h [1];
  h[0] = 2.;
  h[1] = 2.;
  try
  {
    cout << chabis(makeFunctor((DL2D*)0, f), n, x, h, fv, 0.1) << endl;
  }
  catch (err_chabis &ex)
  {
    cout << ex.message << " " << ex.ier << endl;
  }
  cout << "chabis results" << endl;
  cout << "solution  = " << setw(10) << x[0] << setw(10) << x[1] << endl;
  cout << "functions = " << setw(10) << fv[0] << setw(10) << fv[1] << endl;
  cout << "reference results from IBM PC/AT" << endl;
  cout << "estimate of solution   2.000000    1.000000" << endl;
  cout << "values of nonlinear functions   0.0000e+00  0.0000e+00" << endl;
  return 0;
}

