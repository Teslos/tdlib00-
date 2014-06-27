#include <interact.h>
#include <fstream>
#include <iomanip>
#include <algorithm>

const char *HELP =
"\nEvgenii Rudnyi (C) All rights reserved, 2000 \n\n\
test for pol_interaction, usage:\n\n\
interact file_name shift\n";

int main(int argc, char *argv[])
{
  double shift = 0.;
  if (argc == 3)
  {
    shift = atof(argv[2]);
  }
  else if (argc != 2)
  {
    cout << HELP;
    return 0;
  }
  try
  {
    StateX x(2);
    vec_double res(2);
    FuncTpx t1;
    StateTp Tp;
    for (x[1] = 0.; x[1] < 1.01; x[1] += 0.1)
    {
      x[0] = 1. - shift - x[1];
      cout << setw(10) << x[0] << setw(10) << x[1] << endl;
      ifstream in(argv[1]);
      while (in >> ws, in.peek() != EOF)
      {
        in >> t1;
        t1.SetFormalism(::NoChange);
        cout << setw(13) << t1.Z(::G, Tp, x);
        t1.SetFormalism(::Muggianu);
        cout << setw(13) << t1.Z(::G, Tp, x);
        t1.SetFormalism(::Kohler);
        cout << setw(13) << t1.Z(::G, Tp, x);
        t1.SetFormalism(::NoChange);
        res[0] = res[1] = 0.;
        t1.dZdx(::G, Tp, x, res);
        cout << setw(13) << res[1] - res[0];
        t1.SetFormalism(::Muggianu);
        res[0] = res[1] = 0.;
        t1.dZdx(::G, Tp, x, res);
        cout << setw(13) << res[1] - res[0];
        res[0] = res[1] = 0.;
        t1.SetFormalism(::Kohler);
        t1.dZdx(::G, Tp, x, res);
        cout << setw(13) << res[1] - res[0];
        cout << endl;
      }
    }
  }
  catch (gError &t)
  {
    cout << "error: " << t.message << endl;
  }
  return 0;
}
