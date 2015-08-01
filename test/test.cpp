#include <iostream>

#include "../gauss_kronrod.hpp"

using namespace std;
int main(int argc, char** argv){
  double  n     = argc>1 ? atof(argv[1]) : 2.0;
  double reltol=1e-3;
  auto f =[n](double x) -> double { return sin(n*x)*sin(n*x); }; 
  auto intf  = gauss_kronrod(f,reltol);
  cout << "int(f) = " << intf << endl;
  return 0;
}
