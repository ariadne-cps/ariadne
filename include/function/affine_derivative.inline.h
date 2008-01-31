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
    
template<class X> inline
Function::AffineDerivative<X>::AffineDerivative()
  : _result_size(0), _argument_size(0), _variables()
{
}

template<class X> inline
Function::AffineDerivative<X>::AffineDerivative(const uint& rs, const uint& as)
  : _result_size(rs), _argument_size(as), _variables(rs,AffineVariable<X>(as))
{
}

template<class X> template<class XX> inline
Function::AffineDerivative<X>::AffineDerivative(const uint& rs, const uint& as, const XX* ptr)
  : _result_size(rs), _argument_size(as), _variables(rs,AffineVariable<X>(as))
{
  for(size_type i=0; i!=rs; ++i) {
    _variables[i].data()[0]=*(ptr++);
    for(size_type j=0; j!=as; ++j) {
      _variables[i].data()[j+1u]=*(ptr++);
    }
  }
}

template<class X> template<class XX> inline
Function::AffineDerivative<X>::AffineDerivative(const AffineDerivative<XX>& other)
  : _result_size(other.result_size()), _argument_size(other.argument_size()), _variables(other.variables())
{
}

template<class X> template<class XX> inline
Function::AffineDerivative<X>&
Function::AffineDerivative<X>::operator=(const AffineDerivative<XX>& other)
{
  this->_result_size=other.result_size();
  this->_argument_size=other.argument_size();
  this->_variables=other.variables();
  return *this;
}


template<class X> inline
Function::AffineDerivative<X>
Function::AffineDerivative<X>::constant(size_type as, const LinearAlgebra::Vector<X>& x)
{
  const size_type& rs=x.size();
  AffineDerivative<X> result(rs,as);
  for(size_type i=0; i!=rs; ++i) {
    result[i].data()[0]=x[i];
  }
  return result;
}

template<class X> inline
Function::AffineDerivative<X>
Function::AffineDerivative<X>::variable(const LinearAlgebra::Vector<X>& x)
{
  const size_type& n=x.size();
  AffineDerivative<X> result(n,n);
  for(size_type i=0; i!=n; ++i) {
    result[i].data()[0]=x[i];
    result[i].data()[i+1u]=1;
  }
  return result;
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
void
Function::AffineDerivative<X>::resize(const size_type& rs, const size_type& as) 
{
  this->_result_size=rs;
  this->_argument_size=as;
  this->_data.resize(rs*(as+1u));
}

template<class X> inline
const array< Function::AffineVariable<X> >&
Function::AffineDerivative<X>::variables() const
{
  return this->_variables;
}


template<class X> inline
const Function::AffineVariable<X>&
Function::AffineDerivative<X>::operator[](size_type i) const
{
  return this->_variables[i]; 
}

template<class X> inline
Function::AffineVariable<X>&
Function::AffineDerivative<X>::operator[](size_type i) 
{
  return this->_variables[i]; 
}


template<class X> inline
LinearAlgebra::Vector<X>
Function::AffineDerivative<X>::value() const
{
  LinearAlgebra::Vector<const X> result(this->_result_size);
  for(size_type i=0; i!=this->_result_size; ++i) {
    result[i]=this->_variables[i].value(); 
  }
  return result;
}

template<class X> inline
LinearAlgebra::Matrix<X>
Function::AffineDerivative<X>::jacobian() const
{
  LinearAlgebra::Matrix<X> result(this->_result_size,this->_argument_size);
  for(size_type i=0; i!=this->_result_size; ++i) {
    for(size_type j=0; j!=this->_argument_size; ++j) {
      result[i][j]=this->_variables[i].gradient(j);
    }
  }
  return result;
}

template<class X1, class X2> inline
bool
Function::operator==(const AffineDerivative<X1>& x1, const AffineDerivative<X2>& x2) 
{
  return x1.variables()==x2.variables();
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
  return os << "AffineDerivative("<<x.result_size()<<","<<x.argument_size()<<","<<x.variables()<<")";
}




}
