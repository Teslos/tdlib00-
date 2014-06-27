#include <iomanip>
#include <iostream>
#include <math.h>
#include "toms.h"

void SS(long *n, double *x, double *fv, long *iflag)
{
    fv[0] = 2.*(x[0]-3.) + (x[1]-2.)*(x[1]-2.);
    fv[1] = (x[0]-3.)*(x[0]-3.) + 2*(x[1]-2.);
}

int main()
{
  int n = 2;
  double x[2];
  double l[2];
  double u[2];
  long nbd[2];
  nbd[0] = 2;
  nbd[1] = 2;
  l[0] = 5.;
  l[1] = 3.;
  u[0] = 10.;
  u[1] = 10.;
  x[0] = 1.;
  x[1] = 0.;
//  cout << "guesses? ";
//  cin >> x[0] >> x [1];

  double f = lbfgsb_SS(makeFunctor((VOF*)0, SS), n, n, x, l, u, nbd);

  cout << "lbfgsb_SS results" << endl;
  cout << "solution  = " << setw(10) << x[0] << setw(10) << x[1] << endl;
  cout << "SS = " << f << endl;
  return 0;
}

