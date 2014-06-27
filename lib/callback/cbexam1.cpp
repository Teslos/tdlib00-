#include <iostream.h>
#include "callback.hpp"

//do5Times() is a function that takes a functor and invokes it 5 times

void do5Times(const CBFunctor1<int> &doIt)
  {
  for(int i=0;i<5;i++)
    doIt(i);
  }

//Here are some standalone functions

void fred(int i){cout << "fred: " << i<<endl;}
int ethel(long l){cout << "ethel: " << l<<endl;return l;}

//Here is a class with a virtual function, and a derived class

class B{
public:
  virtual void ricky(int i)
     {cout << "B::ricky: " << i<<endl;}
};

class D:public B{
public:
  void ricky(int i)
     {cout << "D::ricky: " << i<<endl;}
};

void main()
  {
  //create a typedef of the functor type to simplify dummy argument
  typedef CBFunctor1<int> *FtorType;

  CBFunctor1<int> ftor; //a functor variable
  //make functor from ptr-to-function
  ftor = makeFunctor((FtorType)0,fred);
  do5Times(ftor);
  //note ethel is not an exact match
  ftor = makeFunctor((FtorType)0,ethel);
  do5Times(ftor);

  //create a D object to be a callback target
  D myD;
  //make functor from object and ptr-to-member-func
  ftor = makeFunctor((FtorType)0,myD,&B::ricky);
  do5Times(ftor);
  }


