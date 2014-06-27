#include "formula.h"
#include <iterator>

void print_formula(formula& t)
{
  cout << "formula " << t << ", massa is " << t.mass() << endl;
  cout << "mol_formula "; 
  t.mf().write(cout); 
  cout << endl;
  cout << "number of K atoms is " << t.noa("K") << endl << endl;
}

int main()
{
  formula a;
  print_formula(a);
  formula b("K2SO4");
  a = b;
  print_formula(a);
  istringstream t1("Al2(SO4)3");
  t1 >> a;
  print_formula(a);
  istringstream t2("KSO4-");
  t2 >> b;
  print_formula(b);
  istringstream t3("SO4(e2)  ");
  t3 >> b;
  print_formula(b);
  istringstream t4("C60(e-10)  ");
  t4 >> b;
  print_formula(b);
  istringstream in("YBa2Cu3O6 YBa2Cu3O7 ");
//  do {
    cout << endl;
//    cout << "Enter first formula (spaces are not allowed)" << endl;
    in >> a;
    print_formula(a);
//    cout << "Enter second formula (spaces are not allowed)" << endl;
    in >> b;
    print_formula(b);
    cout << "results of comparisons" << endl;
    cout << "a < b is  " << (a < b) << endl;
//    cout << "a <= b is " << (a <= b) << endl;
    cout << "a == b is " << (a == b) << endl;
    cout << "a != b is " << (a != b) << endl;
		cout << a.cf().size() << endl;
		cout << b.cf().size() << endl;
		string sss = a.cf();
		for (int i = 0; i < sss.size(); ++i)
			cout << i << "x" << sss[i] << "x" << int(sss[i]) << endl; 
//    cout << "a >= b is " << (a >= b) << endl;
//    cout << "a > b is  " << (a > b) << endl;
//  }
//  while (!a.empty() && !b.empty());
  cout << "All entered elements" << endl;
  const set_elem &s = elem::elements();
  copy(s.begin(), s.end(), ostream_iterator<elem>(cout, " "));
  cout << endl;
  return 0;
}
