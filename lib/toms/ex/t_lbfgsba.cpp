#include <iomanip>
#include <iostream>
#include <math.h>
#include "toms.h"

void SS(long *n, double *x, double *fv, long *iflag)
{
    fv[0] = x[0]*x[1] - pow(x[1], 3) - 1.0;
    fv[1] = pow(x[0], 2)*x[1] + x[1] - 5.0;
}

int main()
{
  int n = 2;
  double x[2];
  double l[2];
  double u[2];
  long nbd[2];
  nbd[0] = 0;
  nbd[1] = 0;
  x[0] = 1.;
  x[1] = 0.;
//  cout << "guesses? ";
//  cin >> x[0] >> x [1];

  double f = lbfgsb_SS(makeFunctor((VOF*)0, SS), n, n, x, l, u, nbd);

  cout << "lbfgsb_SS results" << endl;
  cout << "solution  = " << setw(10) << x[0] << setw(10) << x[1] << endl;
  cout << "SS = " << f << endl;
  cout << "reference results from IBM PC/AT" << endl;
  cout << "estimate of solution   2.000000    1.000000" << endl;
  cout << "values of nonlinear functions   0.0000e+00  0.0000e+00" << endl;
  return 0;
}

