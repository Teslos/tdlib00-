#include <simple.h>
#include <fstream.h>
#include <iomanip.h>
#include <algo.h>

const char *HELP =
"\nEvgenii Rudnyi (C) All rights reserved, 2000 \n\n\
test for pol_interaction, usage:\n\n\
smpl file_name \n\
smpl file_name create\n";

void create(ofstream &out)
{
  SimpleSolution t1;
  size_t n;
  do
  {
    cout << "number of components (n)?"
         << endl << " (n = 2 to exit) ";
    cin >> n;
    t1.set(n, n);
    t1.write(out);
  }
  while (n != 2);
}

int main(int argc, char *argv[])
{
  try
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
    ifstream in(argv[1]);
    phase t1;
    StateTp Tp;
    char t;
    while (in >> ws, in.peek() != EOF)
    {
      in >> t1;
      cout << "name " << t1.id() << endl;
      size_t nc = t1.size();
      for (size_t i = 0; i < nc ; ++i)
      cout << t1.comps()[i] << endl;
      cin.get(t);
    }
  }
  catch (gError &t)
  {
    cout << "error: " << t.message << endl;
  }
  return 0;
}
