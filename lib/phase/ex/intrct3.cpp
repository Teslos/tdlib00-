#include <interact.h>
#include <fstream>
#include <iomanip>
#include <algorithm>

const char *HELP =
"\nEvgenii Rudnyi (C) All rights reserved, 2000 \n\n\
test for pol_interaction, usage:\n\n\
interact file_name\n";

int main(int argc, char *argv[])
{
  if (argc != 2)
  {
    cout << HELP;
    return 0;
  }
  try
  {
    StateX x;
    vec_double res;
    SymMatrix mat;
    FuncTpx t1;
    StateTp Tp;
    ifstream in(argv[1]);
    while (in >> ws, in.peek() != EOF)
    {
      in >> t1;
      x.resize(t1.size());
      res.resize(t1.size());
      mat.resize(t1.size());
//      t1.SetFormalism(::NoChange);
//      t1.SetFormalism(::Muggianu);
      t1.SetFormalism(::Kohler);
      for (size_t i = 0; i < t1.size(); ++i)
        x[i] = 1./t1.size();
      if (t1.size())
        x[0] -= 0.1;
      if (t1.size() > 1)
        x[1] -= 0.1;
      cout << "x ";
      for (size_t i = 0; i < t1.size(); ++i)
        cout << x[i] << " ";
      cout << endl;
      cout << "G " << t1.Z(::G, Tp, x) << endl;
      cout << "dGdx ";
      fill(res.begin(), res.end(), 0.);
      t1.dZdx(::G, Tp, x, res);
      copy(res.begin(), res.end(), ostream_iterator<double>(cout, " "));
      cout << endl;
      cout << "d2Gdx2 ";
      fill(mat.begin(), mat.end(), 0.);
      t1.d2Zdx2(::G, Tp, x, mat);
      copy(mat.begin(), mat.end(), ostream_iterator<double>(cout, " "));
      cout << endl;
      for (size_t i = 0; i < t1.size(); ++i)
      {
        double oldx = x[i];
        double step = global::step1*(x[i] + global::step1);
        x[i] += step;
        vec_double res2(t1.size());
        fill(res2.begin(), res2.end(), 0.);
        t1.dZdx(::G, Tp, x, res2);
        cout << "der by x" << i + 1 << endl;
        for (size_t j = 0; j < t1.size(); ++j)
          cout << (res2[j] - res[j])/step << " ";
        cout << endl;
        x[i] = oldx;
      }
    }
  }
  catch (gError &t)
  {
    cout << "error: " << t.message << endl;
  }
  return 0;
}
