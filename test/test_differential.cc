#include <iostream>
#include "numeric.h"
#include "vector.h"
#include "matrix.h"
#include "dense_differential.h"

using std::cout; using std::endl;
using namespace Ariadne;

template<class R, class A, class P>
void henon(R& r, const A& x, const P& p) 
{
  r[0]=p[0]-x[0]*x[0]-p[1]*x[1];
  r[1]=x[0];
}


int main() {

  {
    // Arithmetic
    DenseDifferential<Float> x(2,2,1.0,0);
    DenseDifferential<Float> y(2,2,1.0,1);
    DenseDifferential<Float> z=2*x+y;
    std::cout << x << y << z << std::endl;
    std::cout << z*z << std::endl;
    std::cout << pow(z,2) << std::endl;
    std::cout << pow(z,3) << std::endl;
  }

  {
    // restrict
    double xa[40]={ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,
                    0,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39};
    int pa[2]={0,2};
    array<int> p(pa,pa+2);
    DenseDifferentialVector<Float> x(2,3,3,xa);
    std::cout << x << restrict(x,p) << "\n";
  }

  {
    // restrict
    double xa[40]={ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,
                    0,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39};
    int pa[2]={0,2};
    array<int> p(pa,pa+2);
    DenseDifferentialVector<Float> x(2,3,3,xa);
    std::cout << "restrict:\n" << x << "\n" << restrict(x,p) << "\n\n";
  }

  {
    // expand
    double xa[20]={ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0,11,12,13,14,15,16,17,18,19};
    int pa[2]={1,2};
    array<int> p(pa,pa+2);
    DenseDifferentialVector<Float> x(2,2,3,xa);
    std::cout << "expand:\n" << x << "\n" << expand(x,3,p) << "\n" << restrict(expand(x,3,p),p) << "\n";
  }

  {
    // compose
    Float av[2]={1,2}; Float adv[4]={1,2,3,4};
    Vector<Float> v(2u,av); Matrix<Float> dv(2u,2u,adv);
    std::cout << v << " " << dv << std::endl;

    DenseDifferentialVector<Float> x(v,dv,2);
    DenseDifferentialVector<Float> y(v,dv,2);
    std::cout << x << std::endl;
    std::cout << compose(x,y) << std::endl;
    std::cout << x[0] << std::endl;
    std::cout << std::endl;
  }
    
  {
    // translate
    cout << "translate:"<<endl;
    Float ac[2]={0,1}; Float adv[10]={1,2,3,4,5,6,7,8,9,10};
    Vector<Float> c(2u,ac); 
    Vector<Float> mc=-c;
    DenseDifferentialVector<Float> dv(1u,2u,2u,adv);
    std::cout << "c=" << c << std::endl;
    std::cout << "dv="<< dv << std::endl;
    std::cout << "dv(x+c)" << translate(dv,c) << std::endl;
    std::cout << "dv((x+c)-c)" <<  translate(translate(dv,c),mc) << std::endl;
    std::cout << "dv(x-c)" << translate(dv,mc) << std::endl;
    std::cout << "dv((x+c)-c)" <<  translate(translate(dv,mc),c) << std::endl;
    std::cout << std::endl;
  }

  {
    // evaluate
    cout << "evaluate:"<<endl;
    Float ac[2]={1,2}; Float adv[10]={1,2,3,4,5,6,7,8,9,10};
    Vector<Float> c(2u,ac); 
    DenseDifferentialVector<Float> dv(1u,2u,3u,adv);
    std::cout << "c=" << c << std::endl;
    std::cout << "dv="<< dv << std::endl;
    std::cout << "v(c)" << evaluate(dv,c) << std::endl;
    std::cout << std::endl;
  }

  {
    // inverse
    cout << "inverse:"<<endl;
    Float ac[2]={1,2}; Float adv[12]={0,1,2,3,4,5,0,6,7,8,9,10};
    DenseDifferentialVector<Float> dv(2u,2u,2u,adv);
    std::cout << "dv="<< dv << std::endl;
    std::cout << "inverse(dv)" << inverse(dv) << std::endl;
    std::cout << "inverse(inverse(dv))" << inverse(inverse(dv)) << std::endl;
    std::cout << std::endl;
  }

  {
    // implicit
    cout << "implicit:"<<endl;
    Float adv[20]={0,1,2,3,4,5,6,7,8,9,0,6,7,8,9,10,11,12,13,14};
    //Float adv[20]={0, 3,2,1, 0,0,0,0,0,0, 0, 5,1,1, 0,0,0,0,0,0};
    DenseDifferentialVector<Float> x(2u,3u,2u,adv);
    cout << "x="<< x << endl;
    DenseDifferentialVector<Float> y=implicit(x);
    cout << "implicit(x)" << y << endl;
    DenseDifferentialVector<Float> z(3u,1u,2u);
    z[0]=DenseDifferential<Float>::variable(1u,2u,0.0,0u);
    z[1]=y[0];
    z[2]=y[1];
    cout << "compose(x,i:y)="<<compose(x,z) << endl;
    cout << endl;
  }

}  
