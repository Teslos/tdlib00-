
#include <elem.h>
#include <iomanip>
#include <iterator>

void print_elem(elem& t)
{
  cout << "element " << t << ", massa is " << t.mass() << endl;
}

int main()
{
  elem a;
  elem b("Se");
  elem c = b;
  print_elem(a);
  print_elem(b);
  print_elem(c);
  cout << "press RETURN to continue" << endl;
  cin.get();
  elem::WriteElements(cout);
	istringstream in("K Cl");
//  do
  {
//   cout << endl << "enter two elements ";
    in >> a >> b;
    print_elem(a);
    print_elem(b);
    cout << "(first < second) is " << (a < b) << endl;
    cout << "(first = second) is " << (a == b) << endl;
    cout << "(first > second) is " << (a > b) << endl;
  }
//  while (a && b);
  const set_elem &s = elem::elements();
  cout << "all entered elements" << endl;
  copy(s.begin(), s.end(), ostream_iterator<elem>(cout, " "));
  cout << endl;
  return 0;
}
