#include <species.h>
#include <fstream>
#include <iomanip>

const char *HELP =
"\nEvgenii Rudnyi (C) All rights reserved, 2000 \n\n\
test for species, usage:\n\n\
sps file_name \n";

int main(int argc, char *argv[])
{
  if (argc != 2)
  {
    cout << HELP;
    return 0;
  }
  try
  {
    ifstream in(argv[1]);
    char t;
    species a;
    StateTp Tp; //Tp.T = 1000.; Tp.p = 1
    while (in >> ws, in.peek() != EOF)
    {
      in >> a;
      cout << a << endl;
      cout << setw(10) << a.Z(G, Tp)
           << setw(10) << a.Z(H, Tp)
           << setw(10) << a.Z(S, Tp)
           << setw(10) << a.Z(Cp, Tp)
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
