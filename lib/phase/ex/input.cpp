#include <general.h>
#include <fstream>

int main(int argc, char *argv[])
{
  if (argc != 2)
    return 0;
  ifstream in(argv[1]);
  parser p(in);
  SGML s, s2;
  cout << s << endl;
  cout << "starting file" << endl;
  try
  {
    while (!p.eof())
    {
      p.GetSGML(s);
      if (s.name.empty())
        break;
      cout << s << endl;
      parser p2(s);
      p2.GetSGML(s2);
      cout << s2 << endl;
    }
  }
  catch (gError &t)
  {
    cout << "Error: " << t.message << endl;
  }
  return 0;
}
