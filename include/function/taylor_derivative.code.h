/***************************************************************************
 *            taylor_derivative.code.h
 *
 *  Copyright 2007  Alberto Casagrande, Pieter Collins
 *  Email casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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
 
namespace Ariadne {

template<class X> 
void 
Function::compute_composition(TaylorVariable<X>& z, const TaylorVariable<X>& y, const TaylorDerivative<X>& x)
{
  using namespace std;
  //std::cerr << "y=" << y << std::endl;
  //std::cerr << "x=" << x << std::endl;

  assert(z.degree()==x.degree());
  assert(z.degree()==y.degree());
  assert(y.argument_size()==x.result_size());
  size_type d=z.degree();
  size_type ms=x.result_size();
  size_type as=x.argument_size();
  
  TaylorDerivative<X> w=x;
  for(uint i=0; i!=ms; ++i) {
    w[i].value()=0;
  }

  TaylorVariable<X> r(as,d);
  TaylorVariable<X> t(as,d);

  // Use inefficient brute-force approach with lots of storage...
  array< array< TaylorVariable<X> > > val(ms, array< TaylorVariable<X> >(d+1));
  for(uint j=0; j!=ms; ++j) {
    val[j][0]=TaylorVariable<X>::constant(as,d,1.0);
    for(uint k=1; k<=d; ++k) {
      val[j][k]=val[j][k-1]*w[j];
    }
  }
  for(MultiIndex i(ms); i.degree()<=d; ++i) {
    t=y[i];
    for(uint j=0; j!=ms; ++j) {
      t=t*val[j][i[j]];
    }
    r+=t;
  }

  z=r;

  //std::cerr << "z=" << z << std::endl;
}




template<class X> 
void 
Function::compute_composition(TaylorDerivative<X>& z, const TaylorDerivative<X>& y, const TaylorDerivative<X>& x)
{
  using namespace std;
  //std::cerr << "y=" << y << std::endl;
  //std::cerr << "x=" << x << std::endl;

  assert(z.degree()==x.degree());
  assert(z.degree()==y.degree());
  assert(y.argument_size()==x.result_size());
  size_type d=z.degree();
  size_type rs=y.result_size();
  size_type ms=x.result_size();
  size_type as=x.argument_size();
  
  TaylorDerivative<X> w=x;
  for(uint i=0; i!=ms; ++i) {
    w[i].value()=0;
  }

  TaylorDerivative<X> r(rs,as,d);
  TaylorVariable<X> t(as,d);

  // Use inefficient brute-force approach with lots of storage...
  array< array< TaylorVariable<X> > > val(ms, array< TaylorVariable<X> >(d+1));
  for(uint j=0; j!=ms; ++j) {
    val[j][0]=TaylorVariable<X>::constant(as,d,1.0);
    for(uint k=1; k<=d; ++k) {
      val[j][k]=val[j][k-1]*w[j];
    }
  }
  for(MultiIndex j(ms); j.degree()<=d; ++j) {
    t=TaylorVariable<X>::constant(as,d,1.0);
    for(uint k=0; k!=ms; ++k) {
      t=t*val[k][j[k]];
    }
    for(uint i=0; i!=rs; ++i) {
      r[i]+=y[i][j]*t;
    }
  }
  z=r;

  //std::cerr << "z=" << z << std::endl;

  //TaylorDerivative<X> t(y.result_size(),x.argument_size(),0);
  //for(uint n=1; n<=d; ++n) {
    //TaylorDerivative<X> u(y.result_size(),x.argument_size(),n);
    //compute_product(u,t,w);
    //u.value()=y[d-n];
    //t=u;
    //std::cerr << "t[" << n << "]=" << t << std::endl;
  //}
  //z=t;
}


}
