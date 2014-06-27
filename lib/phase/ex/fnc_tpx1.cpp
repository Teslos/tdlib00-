#include <func_tpx.h>
#include <fstream.h>
#include <iomanip.h>

const char *HELP =
"\nEvgenii Rudnyi (C) All rights reserved, 2000 \n\n\
test for func_Tp, usage:\n\n\
fnc_Tp file_name \n";


int main(int argc, char *argv[])
{
  if (argc != 2)
  {
    cout << HELP;
    return 0;
  }
  char t;
  try
  {
    ifstream in(argv[1]);
    FuncTpx a;
    StateX x(4);
    x[0] = 0.25;
    x[1] = 0.25;
    x[2] = 0.25;
    x[3] = 0.25;
    StateTp Tp; //Tp.T() = 1000.; Tp.p() = 1
    while (in >> ws, in.peek() != EOF)
    {
      a.read(in);
      a.write(cout);
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
