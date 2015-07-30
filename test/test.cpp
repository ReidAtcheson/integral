#include <iostream>

#include "../gauss_kronrod.hpp"

using namespace std;
int main(int argc, char** argv){
  int  n     = argc>1 ? atoi(argv[1]) : 2;
  double reltol=1e-3;
  auto f =[n](double x) -> double { return pow(x,n); }; 
  auto intf  = gauss_kronrod(f,reltol);
  cout << "int(f) = " << intf << endl;
  return 0;
}
