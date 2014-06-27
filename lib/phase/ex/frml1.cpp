#include "formula.h"
#include <fstream>
#include <iomanip>

const char *HELP =
"\nEvgenii Rudnyi (C) All rights reserved, 2000 \n\n\
test for formula, usage:\n\n\
frml1 file_name \n";

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
    formula a;
    while (in >> ws, in.peek() != EOF)
    {
      in >> a;
      cout << a << endl;
    }
  }
  catch (gError &t)
  {
    cout << "error: " << t.message << endl;
  }
  return 0;
}
