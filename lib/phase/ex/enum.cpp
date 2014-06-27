#include <general.h>

int main()
{
  function f;
  cout << f << endl;
  f = ::H;
  cout << f << endl;
  index i;
  cout << i << endl;
  i = ::mix;
  cout << i << endl;
//  while (1)
  {
    string s1("S"), s2("ideal");
    try
    {
//      cin >> s;
      f = s1;
//      cin >> s;
      i = s2;
      cout << f << " " << i << endl;
    }
    catch (gError &e)
    {
      cout << "error - " << e.message << endl;
    }
  }
  return 0;
}
