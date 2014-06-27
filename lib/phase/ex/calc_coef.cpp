#include <coef.h>
#include <fstream>

int main(int argc, char *argv[])
{
  if (argc != 2)
    return 0;
  ifstream in(argv[1]);
  try
  {
    CalculatorWithCoefs c;
    in >> c;
    cout << c;
    cout << "result is " << c.est() << endl;
    CalculatorWithCoefs d(c);
    cout << d;
    cout << "result is " << d.est() << endl;
    CalculatorWithCoefs e;
    e = d;
    cout << e;
    cout << "result is " << e.est() << endl;
  }
  catch (gError &e)
  {
    cout << "Error: " << e.message << endl;
  }
  return 0;
}
