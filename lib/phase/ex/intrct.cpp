#include <interact.h>
#include <fstream>
#include <iomanip>
#include <algorithm>

const char *HELP =
"\nEvgenii Rudnyi (C) All rights reserved, 2000 \n\n\
test for pol_interaction, usage:\n\n\
interact file_name \n\
interact file_name create\n";

void create(ofstream &out)
{
  interaction *t1;
  vec_int v(2);
  v[0] = 0;
  v[1] = 1;
  t1 = new RedlichKister;
  t1->set(v, func_Tp("Cp_zero"), 5);
  t1->write(out);
  delete t1;
  t1 = new HochArpshofen;
  t1->set(v, func_Tp("Cp_zero"), 5);
  t1->write(out);
  delete t1;
  t1 = new Polynomial;
  int m, p;
  do
  {
    cout << "order of interaction (m) and power of polynomial (p)?"
         << endl << " (p = 0 to exit) ";
    cin >> m >> p;
    v.resize(m);
    for (int j = 0; j < m; ++j)
      v[j] = j;
    t1->set(v, func_Tp("Cp_zero"), p);
    t1->write(out);
  }
  while (p != 0);
}

int main(int argc, char *argv[])
{
  if (argc == 3)
  {
    ofstream out(argv[1]);
    create(out);
    return 0;
  }
  else if (argc != 2)
  {
    cout << HELP;
    return 0;
  }
  char t;
  try
  {
    ifstream in(argv[1]);
    FuncTpx a;
    StateTp Tp; //Tp.T() = 1000.; Tp.p() = 1
    while (in >> ws, in.peek() != EOF)
    {
      int nc = 0;
      StateX x;
      a.read(in);
      a.write(cout);
      if (a.size() > 0)
      {
        const vec_int &v = a.GetVars();
        for (int i = 0; i < v.size(); ++i)
          if (v[i] > nc)
            nc = v[i];
        x.resize(nc + 1);
        for (int i = 0; i < v.size(); ++i)
          x[v[i]] = double(1)/v.size();
      }
      cout << setw(10) << a.Z(G, Tp, x)
           << setw(10) << a.Z(H, Tp, x)
           << setw(10) << a.Z(S, Tp, x)
           << setw(10) << a.Z(Cp, Tp, x)
           << endl;
      cout << "press Enter to continue" << endl;
      cin.get(t);
    }
  }
  catch (gError &t)
  {
    cout << "error: " << t.message << endl;
  }
  return 0;
}
