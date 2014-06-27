#include <coef.h>
#include <fstream>
#include <iomanip>

int main(int argc, char *argv[])
{
  if (argc != 2)
    return 1;
  try
  {
    ifstream in(argv[1]);
    coef c;
    vec_coef v;
    while (in >> ws, in.peek() != EOF)
    {
      in >> c;
      v.push_back(c);
    }
    for (vec_coef_i i = v.begin(); i != v.end(); ++i)
      cout << setw(10) << (*i).x();
    cout << endl;
    for (vec_coef_i i = v.begin(); i != v.end(); ++i)
      cout << *i << endl;
    cout << "Computed" << endl;
    coef::EstimateComputed();
    cout << "ok" << endl;
    for (vec_coef_i i = v.begin(); i != v.end(); ++i)
      cout << setw(10) << (*i).x();
    cout << endl;
    coef::SetNonPrinted();
    for (vec_coef_i i = v.begin(); i != v.end(); ++i)
      cout << *i << endl;
  }
  catch (gError &t)
  {
    cout << "Error: " << t.message << endl;
  }
  return 0;
}
