#include <general.h>


int main()
{
	string text = "4.5 -6.8 inf -inf 10.9";
	istringstream in(text);
  double d;
  while (in >> ws)
  {
    in >> CheckInf(d);
    cout << CheckInf(d) << endl;
  }
  return 0;
}
