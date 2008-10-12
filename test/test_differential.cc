/***************************************************************************
 *            test_differential.cc
 *
 *  Copyright 2008  Pieter Collins
 * 
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */
 
#include <iostream>
#include "numeric.h"
#include "vector.h"
#include "matrix.h"
#include "dense_differential.h"
#include "sparse_differential.h"
#include "differential_vector.h"

using std::cout; using std::endl; using std::flush;
using namespace Ariadne;

template<class R, class A, class P>
void henon(R& r, const A& x, const P& p) 
{
  r[0]=p[0]-x[0]*x[0]-p[1]*x[1];
  r[1]=x[0];
}

template<class R, class A>
void spiral(R& r, const A& x) 
{
  r[0]=-0.8*x[0]+0.4*x[1]-1.0;
  r[1]=-0.4*x[0]-0.8*x[1];
}

template<class DIFF>
DifferentialVector<DIFF> 
henon(const DifferentialVector<DIFF>& x, const Vector<typename DIFF::ScalarType>& p) 
{
  DifferentialVector<DIFF> r(2,2,x.degree()); henon(r,x,p); return r;
}


template<class DIFF>
int 
test_differential()
{
  typedef typename DIFF::ScalarType ScalarType;
  typedef DIFF DifferentialType;
  typedef DifferentialVector<DIFF> DifferentialVectorType;
  
  {
    //Indexing
    MultiIndex a(4);
    for(uint i=0; i!=100; ++i) { cout << a << "\n"; ++a; }
    
    DifferentialType x(2,3);
    DifferentialType y(2,3);
    a=MultiIndex::zero(2); x[a]=2.0;
    a=MultiIndex::unit(2,0); x[a]=1.0;
    a=MultiIndex::unit(2,1); x[a]=1.0;
    a=MultiIndex::zero(2); y[a]=3.0;
    a=MultiIndex::unit(2,0); y[a]=0.0;
    a=MultiIndex::unit(2,1); y[a]=1.0;
    cout << "x=" << x << endl;
    cout << "y=" << y << endl;
  }

  {
    // Arithmetic
    DifferentialType x=DifferentialType::variable(2,2,1.0,0);
    DifferentialType y=DifferentialType::variable(2,2,1.0,1);
    DifferentialType z=2*x+y;
    std::cout << x << y << z << std::endl;
    std::cout << z*z << std::endl;
    std::cout << pow(z,2) << std::endl;
    std::cout << pow(z,3) << std::endl;
  }

  {
    // reciprocal
    DifferentialType x(2,3);
    x[MultiIndex::zero(2)]=2.0;
    x[MultiIndex::unit(2,0)]=1.0;
    x[MultiIndex::unit(2,1)]=1.0;
    cout << "x=" << x << endl;
    cout << "rec(x)=" << rec(x) << endl;
    cout << "rec(x)*x=" << rec(x)*x << endl;
    cout << "rec(x)*x-1=" << rec(x)*x-1.0 << endl;
    cout << endl;
  }

  {
    // restrict
    double xa[40]={ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,
                    0,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39};
    int pa[2]={0,2};
    array<int> p(pa,pa+2);
    DifferentialVectorType x(2,3,3,xa);
    std::cout << x << restrict(x,p) << "\n";
  }

  {
    // restrict
    double xa[40]={ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,
                    0,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39};
    int pa[2]={0,2};
    array<int> p(pa,pa+2);
    DifferentialVectorType x(2,3,3,xa);
    std::cout << "restrict:\n" << x << "\n" << restrict(x,p) << "\n\n";
  }

  {
    // expand
    double xa[20]={ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0,11,12,13,14,15,16,17,18,19};
    int pa[2]={1,2};
    array<int> p(pa,pa+2);
    DifferentialVectorType x(2,2,3,xa);
    std::cout << "expand:\n" << x << "\n" << expand(x,3,p) << "\n" << restrict(expand(x,3,p),p) << "\n";
  }

  {
    // compose
    Float av[2]={1,2}; Float adv[4]={1,2,3,4};
    Vector<Float> v(2u,av); Matrix<Float> dv(2u,2u,adv);
    std::cout << v << " " << dv << std::endl;

    DifferentialVectorType x=DifferentialVectorType::affine(2,2,2,v,dv);
    DifferentialVectorType y=DifferentialVectorType::affine(2,2,2,v,dv);
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
    DifferentialVectorType dv(1u,2u,2u,adv);
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
    DifferentialVectorType dv(1u,2u,3u,adv);
    std::cout << "c=" << c << std::endl;
    std::cout << "dv="<< dv << std::endl;
    std::cout << "v(c)" << evaluate(dv,c) << std::endl;
    std::cout << std::endl;
  }

  {
    // differentiation
    double a[]={ 1,2,3,4,5,6,7,8,9,10 };
    DifferentialType y(2,3,a);
    cout << "y=" << y << endl;
    cout << "derivative(y,0)=" << derivative(y,0) << endl;
    cout << "antiderivative(derivative(y,0),0)-y=" << antiderivative(derivative(y,0),0)-y << endl;
    cout << "derivative(antiderivative(y,0),0)-y=" << derivative(antiderivative(y,0),0)-y << endl;
    cout << endl;
  }

  {
    // mapping 
    DifferentialVectorType x(2,2,3); 
    x[0][MultiIndex::unit(2,0)]=1; x[1][MultiIndex::unit(2,1)]=1; 
    cout << "x=" << x << endl;
    DifferentialVectorType y(2,2,3);
    Vector<Float> a(2); a[0]=1.5; a[1]=0.375;
    y=henon(x,a);
    cout << "h(x)=" << y << endl;
    x=y; y=henon(x,a); cout << "h(h(x))=" << y << endl;
    y=henon(henon(henon(x,a),a),a);
    cout << endl;
  }

  {
    // inverse
    Vector<Float> a(2); a[0]=1.5; a[1]=0.375;
    Vector<Float> v(2); v[0]=2; v[1]=1;
    DifferentialVectorType w(2,2,3);
    DifferentialVectorType x(2,2,3);
    DifferentialVectorType y(2,2,3);
    DifferentialVectorType z(2,2,3);
    w=DifferentialVectorType::variable(2,2,3,v);
    cout << "w=" << w << endl;
    x=henon(w,a);
    cout << "h(w)=" << w << endl;
    y=inverse(x,x.value());
    cout << "hinv(h(w))=" << y << endl;
    z=inverse(y,y.value());
    cout << "hinv(hinv(h(y)))=" << z << endl;
    cout << "h(hinv(hinv(h(y))))=" << henon(z,y.value()) << endl;

    // inverse
    cout << "inverse:"<<endl;
    Float adv[12]={0,1,2,3,4,5,0,6,7,8,9,10};
    DifferentialVectorType dv(2u,2u,2u,adv);
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
    DifferentialVectorType x(2u,3u,2u,adv);
    cout << "x="<< x << endl;
    DifferentialVectorType y=implicit(x);
    cout << "implicit(x)" << y << endl;
    DifferentialVectorType z(3u,1u,2u);
    z[0]=DifferentialType::variable(1u,2u,0.0,0u);
    z[1]=y[0];
    z[2]=y[1];
    cout << "compose(x,i:y)="<<compose(x,z) << endl;
    cout << endl;
  }

  return 0;
}  

int main() {
  test_differential<DenseDifferential<Float> >();
  test_differential<SparseDifferential<Float> >();
}
