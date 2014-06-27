#include <general.h>
#include <fstream.h>

int main(int argc, char *argv[])
{
  if (argc != 2)
    return 0;
  ifstream in(argv[1]);
  parser p(in);
  SGML s;
  cout << s << endl;
  cout << "starting file" << endl;
  try
  {
    while (p)
    {
      p.GetSGML(s);
      cout << s << endl;
    }
  }
  catch (gError &t)
  {
    cout << "Error: " << t.message << endl;
    cout << s << endl;
  }
  return 0;
}
