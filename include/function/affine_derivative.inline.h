/***************************************************************************
 *            affine_derivative.inline.h
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
namespace Function {

template<class X>
class AffineVariableReference
{
 public:
  AffineVariableReference(AffineDerivative<X>& ad, const size_type& i) : _ad(ad), _i(i) { }
  const X& value() const { return _ad._data[_i*(_ad.argument_size()+1u)]; }
  const X& derivative(uint j) const { return _ad._data[_i*(_ad.argument_size()+1u)+j+1u]; }
  void operator=(const AffineVariable<X>& av) { 
    assert(av.argument_size()==_ad.argument_size());
    size_type as=_ad.argument_size();
    _ad._data[_i*(as+1u)]=av.value();
    for(size_type j=0; j!=as; ++j) {
      _ad._data[_i*(as+1u)+j+1u]=av.derivative(j);
    }
  }
private:
  AffineDerivative<X>& _ad; const size_type _i;
};

}
}



namespace Ariadne {
    
template<class X> inline
Function::AffineDerivative<X>::AffineDerivative()
  : _result_size(0), _argument_size(0), _data()
{
}

template<class X> inline
Function::AffineDerivative<X>::AffineDerivative(const uint& rs, const uint& as)
  : _result_size(rs), _argument_size(as), _data(rs*(as+1),0)
{
}

template<class X> template<class XX> inline
Function::AffineDerivative<X>::AffineDerivative(const uint& rs, const uint& as, const XX* ptr)
  : _result_size(rs), _argument_size(as), _data(rs*(as+1),ptr)
{
}

template<class X> template<class XX> inline
Function::AffineDerivative<X>::AffineDerivative(const AffineDerivative<XX>& other)
  : _result_size(other.result_size()), _argument_size(other.argument_size()), _data(other.data())
{
}

template<class X> template<class XX> inline
Function::AffineDerivative<X>&
Function::AffineDerivative<X>::operator=(const AffineDerivative<XX>& other)
{
  this->_result_size=other.result_size();
  this->_argument_size=other.argument_size();
  this->_data=other.data();
  return *this;
}


template<class X> inline
size_type
Function::AffineDerivative<X>::result_size() const
{
  return this->_result_size;
}

template<class X> inline
size_type
Function::AffineDerivative<X>::argument_size() const
{
  return this->_argument_size;
}

template<class X> inline
smoothness_type
Function::AffineDerivative<X>::degree() const
{
  return 1;
}

template<class X> inline
Function::AffineVariable<X>
Function::AffineDerivative<X>::operator[](size_type i) const
{
  return AffineVariable<X>(this->_argument_size,this->_data.begin()+i*(this->_argument_size+1));
}

template<class X> inline
Function::AffineVariableReference<X>
Function::AffineDerivative<X>::operator[](size_type i) 
{
  return AffineVariableReference<X>(*this,i);
}


template<class X1, class X2> inline
bool
Function::operator==(const AffineDerivative<X1>& x1, const AffineDerivative<X2>& x2) 
{
  return x1.argument_size()==x2.argument() && x1.data()==x2.data();
}

template<class X1, class X2> inline
bool
Function::operator!=(const AffineDerivative<X1>& x1, const AffineDerivative<X2>& x2) 
{
  return !(x1==x2);
}

template<class X> inline
std::ostream& 
Function::operator<<(std::ostream& os, const AffineDerivative<X>& x) 
{
  return os << "AffineDerivative("<<x.result_size()<<","<<x.argument_size()<<","<<x.data()<<")";
}




}
