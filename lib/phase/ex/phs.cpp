#include <phase.h>
#include <fstream>
#include <iomanip>
#include <iterator>

const char *HELP =
"\nEvgenii Rudnyi (C) All rights reserved, 1998 \n\n\
test for phase, usage:\n\n\
phase file_name \n";

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
    phase foo;
    while (in >> ws, in.peek() != EOF)
    //main loop continues until the end of the stream in
    {
      in >> foo; //initializing foo by the current phase
      cout << foo;
      size_t N = foo.size();
    //foo.comps() returns vec_species.
    //Its length is equal to the number of components
      if (N == 0)
      {
        cout << "Number of components is equal to zero - skipped" << endl;
        continue;
      }
      cout << "Phase " << foo.id() << ", number of components is " << N << endl;
    //now we need user input for mole fractions
      StateX v(N);
      v[N-1] = 1.;
      if (N > 1)
      {
        cout << "Enter mole fractions for the first " << N - 1 << " components" << endl;
        for (size_t i = 0; i < N; ++i)
        {
          v[i] = 1./N;
        }
/*        for (size_t i = 0; i < N - 1; ++i)
        {
          cout << "x(" <<foo.comps()[i] << ") = ";
          cin >> v[i];
          v[N-1] -= v[i];
        }*/
        cout << "x(" << foo.comps()[N-1] << ") = " << v[N-1] << endl;
    //it would be good to check if mole fractions are feasible
    //but this is left for your homework
      }
    //now we need user input for temperature and pressure
      StateTp Tp;
//      cout << "Enter temperature and pressure ";
//      cin >> Tp.T() >> Tp.p();
    //printing integral thermodynamic properties
      cout << "G = " << foo.G(Tp, v);
      cout << ", Gref = " << foo.Z(::G, ::ref, Tp, v);
      cout << ", Gmix = " << foo.Z(::G, ::mix, Tp, v);
      cout << ", Gideal = " << foo.Z(::G, ::ideal, Tp, v);
      cout << ", Gexcess = " << foo.Z(::G, ::excess, Tp, v);
      cout << endl;
      cout << "H = " << foo.Z(::H, Tp, v);
      cout << ", Href = " << foo.Z(::H, ::ref, Tp, v);
      cout << ", Hmix = " << foo.Z(::H, ::mix, Tp, v);
      cout << ", Hideal = " << foo.Z(::H, ::ideal, Tp, v);
      cout << ", Hexcess = " << foo.Z(::H, ::excess, Tp, v);
      cout << endl;
      cout << "S = " << foo.Z(::S, Tp, v);
      cout << ", Sref = " << foo.Z(::S, ::ref, Tp, v);
      cout << ", Smix = " << foo.Z(::S, ::mix, Tp, v);
      cout << ", Sideal = " << foo.Z(::S, ::ideal, Tp, v);
      cout << ", Sexcess = " << foo.Z(::S, ::excess, Tp, v);
      cout << endl;
    //printing partial properties
      cout << "Chemical potentials" << endl;
      const vec_double &r1 = foo.mu(Tp, v);
      copy(r1.begin(), r1.end(), ostream_iterator<double>(cout, " "));
      cout << endl;
      cout << "Partial enthalpies" << endl;
      const vec_double &r2 = foo.z(::H, Tp, v);
      copy(r2.begin(), r2.end(), ostream_iterator<double>(cout, " "));
      cout << endl;
      cout << "Partial excess entropies" << endl;
      const vec_double &r3 = foo.z(::S, ::excess, Tp, v);
      copy(r3.begin(), r3.end(), ostream_iterator<double>(cout, " "));
      cout << endl;
      cout << "Press Enter to continue" << endl;
      cin.get(t);
    }
  }
  catch (gError &t)
  {
    cout << "error: " << t.message << endl;
  }
  return 0;
}
