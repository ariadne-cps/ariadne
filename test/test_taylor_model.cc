/***************************************************************************
 *            file
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
 
/*! \file orbit.h
 *  \brief Orbits of dynamic systems
 */
#include <iostream>
#include "numeric.h"
#include "vector.h"
#include "matrix.h"
#include "sparse_differential.h"

using namespace std;

template<class R, class A, class P>
void henon(R& r, const A& x, const P& p) 
{
  r[0]=-(x[0]*x[0])+p[0]-p[1]*x[1];
  r[1]=x[0];
}

template<class R, class A>
void spiral(R& r, const A& x) 
{
  r[0]=-0.8*x[0]+0.4*x[1]-1.0;
  r[1]=-0.4*x[0]-0.8*x[1];
}

template<class X>
TaylorModel<X> henon(const TaylorModel<X>& x, const Vector<X>& p) 
{
  TaylorModel<X> r(2,2,x.degree()x.smoothness()); henon(r,x,p); return r;
}

int main() {
  MultiIndex a(4);
  for(uint i=0; i!=100; ++i) { cout << a << "\n"; ++a; }

  a=MultiIndex(4,2)*2+MultiIndex(4,1); cout << a << "\n";

  TaylorModel<Float> x(2,3);
  SparseDifferential<Float> y(2,3);
  a=MultiIndex(2); x[a]=2.0;
  a=MultiIndex(2,0); x[a]=1.0;
  a=MultiIndex(2,1); x[a]=1.0;
  a=MultiIndex(2); y[a]=3.0;
  a=MultiIndex(2,0); y[a]=0.0;
  a=MultiIndex(2,1); y[a]=1.0;
  cout << x << endl;
  cout << x*x << endl;

  a=MultiIndex(1); 
  SparseDifferential<Float> z(1,3);
  z[a]=3; ++a; z[a]=1;
  cout << z << z*z << z*z*z << endl;

  cout << "rec(2,1.5)=" << flush;
  cout << SparseSeries<Float>::rec(8,1.25) << endl;
  cout << "x=" << x << endl;
  cout << "rec(x)=" << rec(x) << endl;
  cout << "rec(x)*x=" << rec(x)*x << endl;
  cout << "rec(x)*x-1=" << rec(x)*x-1.0 << endl;

  y=rec(x);
  cout << "y=" << y << endl;
  cout << "derivative(y,0)=" << derivative(y,0) << endl;
  cout << "antiderivative(derivative(y,0),0)-y=" << antiderivative(derivative(y,0),0)-y << endl;
  cout << "derivative(antiderivative(y,0),0)-y=" << derivative(antiderivative(y,0),0)-y << endl;

  SparseDifferentialVector<Float> w1(1,2,4);
  w1[0][MultiIndex(2,1)*3]=1.0;
  cout << w1 << std::endl;
  SparseDifferentialVector<Float> w2(2,3,4);
  w2[1][MultiIndex(3,1)*1]=1.0;
  cout << w2 << std::endl;

  cout << compose(w1,w2) << endl;
   
  {
    SparseDifferentialVector<Float> x(2,2,3); 
    x[0][MultiIndex(2,0)]=1; x[1][MultiIndex(2,1)]=1; 
    SparseDifferentialVector<Float> y(2,2,3);
    Vector<Float> a(2); a[0]=1.5; a[1]=0.375;
    y=henon(x,a);
    cout << "h(" << x << ")=" << y << endl;
    x=y; y=henon(x,a); cout << "h(" << x << ")=" << y << endl;
    y=henon(henon(henon(x,a),a),a);
    cout << endl;

    SparseDifferentialVector<Float> p(2,3,3);
    p[0][MultiIndex(3,0)]=1;
    p[1][MultiIndex(3,1)]=1;
    SparseDifferentialVector<Float> i(3,2,3);
    i[0][MultiIndex(2,0)]=1;
    i[1][MultiIndex(2,1)]=1;
    i[2][MultiIndex(2)]=1;
    
  }

  {
    Vector<Float> a(2); a[0]=1.5; a[1]=0.375;
    SparseDifferentialVector<Float> x(2,2,3);
    SparseDifferentialVector<Float> y(2,2,3);
    x[0].value()=1; x[0].gradient(0)=1;
    x[1].value()=2; x[1].gradient(1)=1;
    x=henon(x,a);
    y=inverse(x,x.value());
    cout << "x=" << x << "\ninverse(x)=" << y << endl;
    cout << "inverse(inverse(x))=" << inverse(y,y.value()) << endl;
  }

  {
    SparseDifferentialVector<Float> x(1,2,3);
    x[0].value()=0; x[0][MultiIndex(2,0)*2]=1; x[0][MultiIndex(2,1)*2]=1; 
    cout << "implicit(x)=" << implicit(x,x.value()) << endl;
    cout << endl;
  }


  {
    SparseDifferentialVector<Float> A(Vector<Float>(2),Matrix<Float>(2,2),5);
    
    SparseDifferentialVector<Float> f(2,2,5);
    f[0].value()=-1; f[0].gradient(0)=-0.8; f[0].gradient(1)=-0.4;
    f[1].value()=0; f[1].gradient(0)=-0.4; f[1].gradient(1)=-0.8;
    Vector<Float> x(2);
    x[0]=0.125; x[1]=0.25;
    cout << integrate(f,x,4,2) << std::endl << std::endl;
  }

}
