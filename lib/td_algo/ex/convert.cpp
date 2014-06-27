#include <td_algo.h>

int main()
{
  try
  {
    convert cv;
    ifstream in("cv.mod");
    cv.read(in);
    cout << cv;

    convert cv1(cv);
  }
  catch (gError &t)
  {
    cout << t.message << endl;
  }
}
