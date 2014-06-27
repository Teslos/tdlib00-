#include <calc.h>

int main(int argc, char *argv[])
{
  if (argc != 2)
    return 0;
  istringstream in(argv[1]);
  try
  {
    calculator c;
    in >> c;
    cout << c.est() << endl;
    cout << c << endl;
  }
  catch (gError &e)
  {
    cout << "Error: " << e.message << endl;
  }
  return 0;
}
