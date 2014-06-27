#include <func_tp.h>
#include <fstream>
#include <iomanip>

const char *HELP =
"\nEvgenii Rudnyi (C) All rights reserved, 1998 \n\n\
test for func_Tp, usage:\n\n\
fnc_Tp file_name \n\
fnc_Tp file_name create\n";

void create(ofstream &out)
{
  func_Tp a;
  out << a;
  a.create("Cp_zero");
  out << a;
  a.create("Cp_const");
  out << a;
  a.create("Cp_BB2");
  out << a;
  a.create("Cp_BB4");
  out << a;
  a.create("IVT_Tp");
  out << a;
  a.create("SGTE_Tp");
  out << a;
  a.create("ideal_gas");
  out << a;
  a.create("V_const");
  out << a;
  a.create("alpha_const");
  out << a;
  a.create("alpha_kappa_const");
  out << a;
  a.create("alpha_kappa_const2");
  out << a;
  a.create("calc_Tp");
  out << a;
  a.create("complex_Tp");
  out << a;
  a.create("compound_Tp");
  out << a;
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
    func_Tp a;
    StateTp Tp; //Tp.T() = 1000.; Tp.p() = 1
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
