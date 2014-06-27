#include <phase.h>
#include <fstream>
#include <iomanip>
#include <algorithm>

const char *HELP =
"\nEvgenii Rudnyi (C) All rights reserved, 1998 \n\n\
test for derivative, usage:\n\n\
phs_der model data \n";

int main(int argc, char *argv[])
{
  if (argc != 3)
  {
    cout << HELP;
    return 1;
  }
  try
  {
    ifstream in(argv[1]);
    ifstream data(argv[2]);
    phase t1;
    StateTp Tp;
    Tp.T() = 600.;
    index ind;
    in >> t1;
    string str;
    while (data >> ws, data.peek() != EOF)
    {
      int nc;
      nc = t1.size();
      StateX x(nc);
      if (nc > 0)
      {
        data >> Tp.T();
        x[nc - 1] = 1.;
        for (size_t i = 0; i < nc - 1; ++i)
        {
          data >> x[i];
          x[nc - 1] -= x[i];
        }
        data >> str;
        ind = str;
        cout << setw(10) << Tp.T();
        for (size_t i = 0; i < nc; ++i)
          cout << setw(10) << x[i];
        cout << endl;
      }
      cout << "G" << ind << "=" << setw(10) << t1.Z(::G, ind, Tp, x)
           << " H" << ind << "=" << setw(10) << t1.Z(::H, ind, Tp, x)
           << " S" << ind << "=" << setw(10) << t1.Z(::S, ind, Tp, x)
           << endl;
      cout << "dGdx" << endl;
      vec_double mu = t1.pointer()->RefPhase::dZdx(::G, ind, Tp, x);
      for (size_t i = 0; i < t1.size(); ++i)
        cout << setw(13) << mu[i];
      cout << endl;
      cout << "G" << ind << "=" << setw(10) << t1.Z(::G, ind, Tp, x)
           << " numH" << ind << "=" << setw(10); 
      cout << t1.pointer()->RefPhase::num_H(ind, Tp, x)
           << " numS" << ind << "=" << setw(10);
      cout << t1.pointer()->RefPhase::num_S(ind, Tp, x)
           << endl;
      mu = t1.pointer()->RefPhase::num_dGdx(ind, Tp, x);
      cout << "numerical dGdx" << endl;
      for (size_t i = 0; i < t1.comps().size(); ++i)
        cout << setw(13) << mu[i];
      cout << endl;
      cout << endl;
    }
  }
  catch (gError &t)
  {
    cout << "error: " << t.message << endl;
  }
  return 0;
}
