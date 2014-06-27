#include <general.h>
#include <iostream>
#include <iomanip>

int main()
{
  SGML el;
  el.body = "tt rrr kkk ddd   jkjkjkjkj";
  parser p(el);
  p.GetToken();
	cout << p.token << endl;
  p.GetToken();
	cout << p.token << endl;
  p.GetToken();
	cout << p.token << endl;
  parser p2(el);
  p2.GetToken();
	cout << p2.token << endl;
  p2.GetToken();
	cout << p2.token << endl;
  p2.GetToken();
	cout << p2.token << endl;

}
