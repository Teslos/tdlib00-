#include <species.h>
#include <fstream>
#include <iomanip>

const char *HELP =
"\nEvgenii Rudnyi (C) All rights reserved, 1998 \n\n\
test for species, usage:\n\n\
species file_name \n";

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
    species a;
    StateTp Tp; //Tp.T = 1000.; Tp.p = 1
    double H298;
    cout << "   T         Cp        S        H-H298       G " << endl;
    while (in >> ws)
    {
      in >> a;
      cout << a.id() << endl;
      Tp.T() = 298.15;
      cout << setw(10) << Tp.T()
           << setw(10) << a.Z(Cp, Tp)
           << setw(10) << a.Z(S, Tp)
           << setw(10) << (H298 = a.Z(H, Tp))
           << setw(10) << a.Z(G, Tp)
           << endl;
      for (Tp.T() = 300.; Tp.T() < 1501.; Tp.T() += 100.)
        cout << setw(10) << Tp.T()
             << setw(10) << a.Z(Cp, Tp)
             << setw(10) << a.Z(S, Tp)
             << setw(10) << a.Z(H, Tp) - H298
             << setw(10) << a.Z(G, Tp)
             << endl;

    }
  }
  catch (gError &t)
  {
    cout << "error: " << t.message << endl;
  }
  return 0;
}
