// this is a modified example from gnussl

#include <iostream>
#include <iomanip>
#include <vector>
#include <time.h>
#include <lapack.h>

int main(int argc,char **argv) {
  if (argc !=2)
  {
    cout << "gauss dim" << endl;
    return 0;
  }
  time_t t1 = time(0);
  int dim=atoi(argv[1]);
  vector<double> a(dim*dim);
  vector<double> b(dim);
  vector<int> ipiv(dim);
  srand(1);              // seed the random # generator with a known value
  double maxr=(double)RAND_MAX;
  for(int r=0; r < dim; r++) {  // set a to a random matrix, i to the identity
    for(int c=0; c < dim; c++) {
      a[r + c*dim] = rand()/maxr;
    }
    b[r] = rand()/maxr;
  }
  vector<double> a1(a);
  vector<double> b1(b);
  int info;
  cout << "matrices allocated and initialised " << time(0) - t1 << endl;
  try {
    time_t t2 = time(0);
    clock_t c2 = clock();
    dgetf2(dim, dim, &*a.begin(), &*ipiv.begin(), info);
    cout << "info is " << info << endl;
    dgetrs(dim, dim, 1, &*a.begin(), &*ipiv.begin(), &*b.begin());
    clock_t c3 = clock();
    cout << "gauss_jordan is over for " << time(0) - t2
      << " " << c3 - c2 << endl;
    double eps = 0.;
    for (int i = 0; i < dim; ++i)
    {
      double sum = 0.;
      for (int j = 0; j < dim; ++j)
        sum += a1[i + j*dim]*b[j];
      eps += fabs(b1[i] - sum);
    }
    cout << "check is " << eps << endl;
  }
  catch (...) {
    cout << "error" << endl;
    exit(0);
  }
  cout << CLK_TCK << endl;
  return 0;
}
