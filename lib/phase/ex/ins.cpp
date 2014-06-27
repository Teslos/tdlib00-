#include <formula.h>
#include <fstream.h>
#include <algo.h>

int main(int argc, char *argv[])
{
  if (argc != 2)
  {
    cout << "name of file missed" << endl;
    return 0;
  }
  try
  {
    vec_formula vf;
    ifstream in(argv[1]);
    copy(istream_iterator<formula>(in), istream_iterator<formula>(), 
        back_insert_iterator<vec_formula>(vf));
    copy(vf.begin(), vf.end(), ostream_iterator<formula>(cout, " "));
    cout << endl;
  }
  catch (gError &t)
  {
    cout << "error: " << t.message << endl;
  }
}
