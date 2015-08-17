#include <iostream>
#include <cassert>

#include "../gauss_kronrod.hpp"

bool almost_equal(double,double,double);
double exact_integral(int freq);

using namespace std;
int main(int argc, char** argv){
  valarray<double> reltols = {1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-10,1e-14};
  valarray<int> frequencies= {1,2,3,4,5,10,20,40};
  for(auto r : reltols){
    for(auto f : frequencies){
      auto g =[f](double x) -> double { return sin(f*x)*sin(f*x); };

      cout << endl;
      cout << "TEST: frequency = "<< f << endl;
      cout << "TEST: relative tolerance = " << r << endl;
      auto intf = gauss_kronrod(g,r);



      assert(almost_equal(exact_integral(f),intf,r));

      cout << "SUCCEEDED" << endl;
      cout << endl;

    }
  }
  return 0;
}


double exact_integral(int freq){
  assert(freq>0);
  return 1.0 - sin(freq)*cos(freq)/freq;
}
bool almost_equal(double exact, double approx, double reltol){
  return fabs(exact-approx)<reltol*exact;
}
