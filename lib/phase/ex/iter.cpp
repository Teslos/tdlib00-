#include <interact.h>
#include <simple.h>

void check_Polynomial_iterator()
{
  Polynomial t;
  int m = 4, p = 3;
//  for (;p != 0;)
  {
//    cout << "order of interaction (m) and power of polynomial (p)?"
//         << endl << " (p = 0 to exit) ";
//    cin >> m >> p;
    vec_int v(m);
    t.set(v, func_Tp(), p);
    cout << "m_in_p " << t.size() << endl;
    for (Polynomial::iterator ip = t.begin(); !(ip == t.end()); ++ip)
    {
      for (size_t j = 0; j < m; ++j)
        cout << " " << ip[j];
      cout << endl;
    }
  }
}

void check_index_x()
{
  int m = 3, n = 4;
//  for (;n!=2;)
  {
//    cout << "order of interaction (m) and number of components (n) ? "
//         << endl << " (n = 2 to exit) ";
//    cin >> m >> n;
    cout << "m_in_n " << bin_coef(m, n) << endl;
    index_x ix(m, n);
    for (; ix; ++ix)
    {
      for (int j = 0; j < m; ++j)
        cout << " " << ix[j];
      cout << endl;
    }
  }
}

int main()
{
  check_Polynomial_iterator();
  check_index_x();
  return 0;
}

