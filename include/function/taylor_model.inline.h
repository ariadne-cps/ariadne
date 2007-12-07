/***************************************************************************
 *            taylor_model.inline.h
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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

#include "sorted_index.h"
#include "multi_index.h"

namespace Ariadne {


template<class R> template<class RR> inline
Function::TaylorModel<R>::TaylorModel(const size_type& rs, const size_type& as, const size_type& d, const size_type& s, const RR* p)
  : _result_size(rs), 
    _argument_size(as), 
    _order(d),
    _smoothness(s),
    _data(p,p+rs*Numeric::bin(d+as,as)) 
{
}

template<class R> inline
Function::TaylorModel<R>::TaylorModel(const TaylorModel<R>& p)
  : _result_size(p.result_size()), 
    _argument_size(p.argument_size()),
    _order(p.order()),
    _smoothness(p.smoothness()),
    _data(p.data())
{
}
        
template<class R> template<class RR> inline
Function::TaylorModel<R>::TaylorModel(const TaylorModel<RR>& p)
  : _result_size(p.result_size()), 
    _argument_size(p.argument_size()),
    _order(p.order()),
    _smoothness(p.smoothness()),
    _data(p.data())
{
}
        
template<class R> inline
Function::TaylorModel<R>&
Function::TaylorModel<R>::operator=(const TaylorModel<R>& p)
{
  if(this!=&p) {
    this->_result_size=p.result_size(); 
    this->_argument_size=p.argument_size();
    this->_order=p.order();
    this->_smoothness=p.smoothness();
    this->_data=p.data();
  }
  return *this;
}
        

template<class R> template<class RR> inline
Function::TaylorModel<R>&
Function::TaylorModel<R>::operator=(const TaylorModel<RR>& p)
{
  this->_result_size=p.result_size(); 
  this->_argument_size=p.argument_size();
  this->_order=p.order();
  this->_smoothness=p.smoothness();
  this->_data=p.data();

  return *this;
}
        

template<class R1,class R2> inline
Function::TaylorModel<typename Numeric::traits<R1,R2>::arithmetic_type>
Function::add(const TaylorModel<R1>& p1, const TaylorModel<R2>& p2)
{
  TaylorModel<typename Numeric::traits<R1,R2>::arithmetic_type> p0;
  add(p0,p1,p2);
  return p0;
}

template<class R1,class R2> inline
Function::TaylorModel<typename Numeric::traits<R1,R2>::arithmetic_type>
Function::mul(const TaylorModel<R1>& p1, const TaylorModel<R2>& p2)
{
  TaylorModel<typename Numeric::traits<R1,R2>::arithmetic_type> p0;
  mul(p0,p1,p2);
  return p0;
}

template<class R1,class R2> inline
Function::TaylorModel<typename Numeric::traits<R1,R2>::arithmetic_type>
Function::operator+(const TaylorModel<R1>& p1, const TaylorModel<R2>& p2)
{
  TaylorModel<typename Numeric::traits<R1,R2>::arithmetic_type> p0;
  add(p0,p1,p2);
  return p0;
}

template<class R1,class R2> inline
Function::TaylorModel<typename Numeric::traits<R1,R2>::arithmetic_type>
Function::operator-(const TaylorModel<R1>& p1, const TaylorModel<R2>& p2)
{
  TaylorModel<typename Numeric::traits<R1,R2>::arithmetic_type> p0;
  sub(p0,p1,p2);
  return p0;
}

template<class R1,class R2> inline
Function::TaylorModel<typename Numeric::traits<R1,R2>::arithmetic_type>
Function::operator*(const TaylorModel<R1>& p1, const TaylorModel<R2>& p2)
{
  TaylorModel<typename Numeric::traits<R1,R2>::arithmetic_type> p0;
  mul(p0,p1,p2);
  return p0;
}

template<class R1,class R2> inline
Function::TaylorModel<typename Numeric::traits<R1,R2>::arithmetic_type>
Function::operator*(const TaylorModel<R1>& p1, const R2& x2)
{
  TaylorModel<typename Numeric::traits<R1,R2>::arithmetic_type> p0(p1);
  Function::scale(p0,x2);
  return p0;
}

template<class R1,class R2> inline
Function::TaylorModel<typename Numeric::traits<R1,R2>::arithmetic_type>
Function::operator/(const TaylorModel<R1>& p1, const R2& x2)
{
  return p1*static_cast<typename Numeric::traits<R1,R2>::arithmetic_type>(1/x2);
}

template<class R1,class R2> inline
Function::TaylorModel<typename Numeric::traits<R1,R2>::arithmetic_type>
Function::operator*(const R1& x1, const TaylorModel<R2>& p2)
{
  return p2*x1;
}

template<class R> inline
Function::TaylorModel<typename Numeric::traits<R>::arithmetic_type>
Function::pow(const TaylorModel<R>& p, const unsigned int& n)
{
  TaylorModel<typename Numeric::traits<R>::arithmetic_type> r;
  pow(r,p,n);
  return r;
}

template<class R1,class R2> inline
Function::TaylorModel<typename Numeric::traits<R1,R2>::arithmetic_type>
Function::compose(const TaylorModel<R1>& p1, const TaylorModel<R2>& p2)
{
  TaylorModel<typename Numeric::traits<R1,R2>::arithmetic_type> p0;
  compose(p0,p1,p2);
  return p0;
}

template<class R> inline
Function::TaylorModel<typename Numeric::traits<R>::arithmetic_type>
Function::derivative(const TaylorModel<R>& p1, const size_type& k)
{
  TaylorModel<typename Numeric::traits<R>::arithmetic_type> p0;
  derivative(p0,p1,k);
  return p0;
}




} // namespace Ariadne
