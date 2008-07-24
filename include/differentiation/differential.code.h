/***************************************************************************
 *            differential.code.h
 *
 *  Copyright 2007  Pieter Collins
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
 
#include "linear_algebra/vector.h"
#include "differentiation/multi_index.h"
#include "differentiation/power_series.h"

namespace Ariadne {

namespace {

template<class X>
void 
compute_composition(Differential<X>& z, 
                    const PowerSeries<X>& y, 
                    const Differential<X>& x)
{
  
  size_type as=x.argument_size();
  size_type d=z.degree();

  Differential<X> w=x;
  w.value()=0;
  Differential<X> t(as,d);
  t.value()=y.data()[d];
  for(uint n=1; n<=d; ++n) {
    Differential<X> u(as,d);
    acc(u,t,w);
    t=u+y.data()[d-n];
  }
  z=t;
  return;
}

} // namespace




template<class X> 
Differential<X>::Differential()
  : _argument_size(1), _degree(0), _data(1u,X(0)) 
{
}

template<class X> 
Differential<X>::Differential(size_type a, smoothness_type d)
  : _argument_size(a), _degree(d), _data(compute_polynomial_data_size(a,d),X(0))
{
}



template<class X> 
size_type 
Differential<X>::argument_size() const 
{ 
  return this->_argument_size;
}

template<class X> 
smoothness_type 
Differential<X>::degree() const 
{ 
  return this->_degree;
}

template<class X> 
const X&
Differential<X>::value() const 
{ 
  return this->_data[0];
}

template<class X> 
X&
Differential<X>::value()  
{ 
  return this->_data[0];
}

template<class X> 
const X&
Differential<X>::gradient(size_type j) const 
{ 
  return this->_data[j+1u];
}

template<class X> 
array<X>& 
Differential<X>::data()
{
  return this->_data; 
}

template<class X> 
const array<X>& 
Differential<X>::data() const 
{
  return this->_data; 
}


template<class X> 
X& 
Differential<X>::operator[](const MultiIndex& a) 
{ 
  return this->_data[a.position()]; 
}

template<class X> 
const X& 
Differential<X>::operator[](const MultiIndex& a) const 
{ 
  return this->_data[a.position()]; 
}

template<class X> 
bool
operator<(const Differential<X>& x1, const Differential<X>& x2)
{
  return x1.value() < x2.value();
}

template<class X> 
Differential<X>&
acc(Differential<X>& r, const Differential<X>& x1, const Differential<X>& x2)
{
  ARIADNE_ASSERT(r.argument_size()==x1.argument_size());
  ARIADNE_ASSERT(r.argument_size()==x2.argument_size());
  for(MultiIndex i1(x1.argument_size()); i1.degree() <= std::min(r.degree(),x1.degree()); ++i1) {
    for(MultiIndex i2(x2.argument_size()); i2.degree() <= std::min(x2.degree(),smoothness_type(r.degree()-i1.degree())); ++i2) {
      MultiIndex i0=i1+i2;
      r[i0]+=x1[i1]*x2[i2];
    }
  }
  return r;
}

template<class X> 
Differential<X>&
acc(Differential<X>& r, const X& c, const Differential<X>& x)
{
  ARIADNE_ASSERT(r.argument_size()==x.argument_size());
  size_type n=std::max(r.data().size(),x.data().size());
  for(size_type i=0; i!=n; ++i) {
    r.data()[i]+=c*x.data()[i];
  }
  return r;
}


template<class X> 
Differential<X>&
Differential<X>::assign(const Differential<X>& x)
{
  ARIADNE_ASSERT(this->argument_size()==x.argument_size());
  ARIADNE_ASSERT(this->degree()>=x.degree());
  for(uint i=0; i!=x.data().size(); ++i) {
    this->_data[i]=x.data()[i]; 
  }
  return *this;
}


template<class X> 
Differential<X>&
operator+=(Differential<X>& r, const Differential<X>& x)
{
  ARIADNE_ASSERT(r.argument_size()==x.argument_size());
  ARIADNE_ASSERT(r.degree()==x.degree());
  for(uint i=0; i!=r.data().size(); ++i) {
    r.data()[i]+=x.data()[i];
  }
  //reinterpret_cast<Vector<X>&>(r._data)
  //  += reinterpret_cast<Vector<X>const&>(x._data);
  return r;
}


















template<class X> 
Differential<X> 
compose(const Differential<X>& y, const Differential<X>& x)
{
  assert(y.argument_size()==1);
  PowerSeries<X> t(y.degree(),y.data().begin());
  Differential<X> z(x.argument_size(),std::min(x.degree(),t.degree()));
  compute_composition(z,t,x);
  return z;
}

template<class X> 
Differential<X> 
compose(const PowerSeries<X>& y, const Differential<X>& x)
{
  Differential<X> z(x.argument_size(),std::min(x.degree(),y.degree()));
  compute_composition(z,y,x);
  return z;
}

template<class X>
Differential<X> 
embed(const Differential<X>& x, 
      uint size, uint start)
{  
  assert(start+x.argument_size()<=size);
  Differential<X> r(size,x.degree());
  MultiIndex jr(size);
  for(MultiIndex jx(x.argument_size()); jx.degree()<=x.degree(); ++jx)
  {
    for(uint k=0; k!=x.argument_size(); ++k) {
      jr.set(start+k,jx[k]);
    }
    r[jr]=x[jx];
  }
  return r;
}


template<class X> 
Differential<X> 
derivative(const Differential<X>& x, size_type k)
{
  Differential<X> r(x.argument_size(),x.degree()-1);
  MultiIndex e(x.argument_size());
  e.set(k,1);
  for(MultiIndex j(r.argument_size()); j.degree() <= r.degree(); ++j) {
    r[j]=(j[k]+1)*x[j+e];
  }
  return r;
}


template<class X> 
Differential<X> 
antiderivative(const Differential<X>& x, size_type k)
{
  Differential<X> r(x.argument_size(),x.degree()+1);
  MultiIndex e(x.argument_size());
  e.set(k,1);
  for(MultiIndex j(r.argument_size()); j.degree() <= x.degree(); ++j) {
    r[j+e]=x[j]/(j[k]+1);
  }
  return r;
}





template<class X> 
Differential<X> 
reduce(const Differential<X>& x, const size_type& d)
{
  assert(x.degree()>=d);
  Differential<X> r(x.argument_size(),d);
  for(MultiIndex i(x.argument_size()); i.degree() <= x.degree(); ++i) {
    r[i]=x[i];
  }
}








template<class X> 
std::ostream& 
operator<<(std::ostream& os, const Differential<X>& x) {
  //  return os << "Differential( argument_size=" << x.argument_size() << ", degree=" << x.degree() << ", data=" << x.data() << ")";
  //os << "Differential(";
  os << "V";
  size_type degree=0;
  for(MultiIndex i(x.argument_size()); i.degree()<=x.degree(); ++i) {
    if(i.degree()==0) {
      os << '[';
    } else if(i.degree()==degree) {
      os << ',';
    } else {
      degree=i.degree();
      os << ';';
    }
    os << x[i];
  }
  os << ']';
  //os << ")";
  return os;

//  return os << "Differential( argument_size=" << x.argument_size() << ", degree=" << x.degree() << ", data=" << x.data() << ")";
}


template<class X> 
void
Differential<X>::instantiate() 
{
  const int* n=0;
  X* x=0;
  Vector<X>* v=0;
  PowerSeries<X>* ts=0;
  Differential<X>* tv=0;
  std::ostream* os = 0;

  acc(*tv,*x,*tv);
  acc(*tv,*tv,*tv);

  evaluate(*tv,v->data());
  compose(*ts,*tv);
  derivative(*tv,0u);
  antiderivative(*tv,0u);
 
  operator<(*tv,*tv);

  operator+(*tv);
  operator-(*tv);
  operator+(*tv,*tv);
  operator-(*tv,*tv);
  operator*(*tv,*tv);
  operator/(*tv,*tv);

  operator+=(*tv,*tv);
  operator+=(*tv,*n);
  operator*=(*tv,*n);
  operator+=(*tv,*x);
  operator*=(*tv,*x);
  
  operator<<(*os,*tv);
}


} //namespace Ariadne
