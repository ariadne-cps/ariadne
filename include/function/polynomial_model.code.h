/***************************************************************************
 *            polynomial_model.code.h
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

#include "polynomial_model.h"

#include <string>
#include <sstream>

#include <list>
#include <set>
#include <vector>
#include <valarray>

#include "../base/stlio.h"
#include "../base/array.h"
#include "../numeric/interval.h"

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

#include "../function/exceptions.h"

namespace Ariadne {

template<class R>
Function::Monomial<R>::Monomial(const std::string& s) 
{
  std::stringstream ss(s);
  ss >> *this;
}

template<class R>
Function::Polynomial<R>::Polynomial(const std::string& s) 
{
  //std::cerr << "Polynomial(\"" << s << "\")" << std::endl;
  std::stringstream ss(s);
  ss >> *this;
}


template<class R>
Function::PolynomialModel<R>::PolynomialModel(const std::string& s) 
{
  //std::cerr << "PolynomialModel(\"" << s << "\")" << std::endl;
  std::stringstream ss(s);
  ss >> *this;
}


template<class R>
size_type
Function::Monomial<R>::degree() const 
{
  size_type result=0;
  for(size_type i=0; i!=this->argument_size(); ++i) {
    result+=this->index(i);
  }
  return result;
}

template<class R>
void
Function::Monomial<R>::_resize(const size_type& n) 
{
  this->_multi_index.resize(n);
}

template<class R>
void
Function::Polynomial<R>::_sort() 
{
  std::sort(this->_terms.begin(), this->_terms.end());
  
  typedef typename std::vector< Monomial<R> >::const_iterator terms_const_iterator;
  std::vector< Monomial<R> > new_terms;
  for(terms_const_iterator i=this->_terms.begin(); i!=this->_terms.end(); ++i) {
    new_terms.push_back(*i);
  }
  this->_terms.swap(new_terms);
}

template<class R>
void
Function::Polynomial<R>::_set_argument_size(const size_type& n) 
{
  this->_argument_size=n;
  for(size_type i=0; i!=this->_terms.size(); ++i) {
    this->_terms[i]._resize(n);
  }
}

template<class R>
size_type
Function::Polynomial<R>::_compute_maximum_term_size() const
{
  size_type result=0;
  for(size_type i=0; i!=this->_terms.size(); ++i) {
    result=std::max(this->_terms[i].argument_size(),result);
  }
  return result;
}


template<class R>
void
Function::PolynomialModel<R>::_set_argument_size(const size_type& n) 
{
  this->_argument_size=n;
  for(size_type i=0; i!=this->_components.size(); ++i) {
    this->_components[i]._set_argument_size(n);
  }
}

template<class R>
size_type
Function::PolynomialModel<R>::_compute_maximum_component_size() const
{
  size_type result=0;
  for(size_type i=0; i!=this->_components.size(); ++i) {
    result=std::max(this->_components[i].argument_size(),result);
  }
  return result;
}



template<class R> 
bool 
Function::operator<(const Monomial<R>& m1, const Monomial<R>& m2)
{
  if(m1.degree() == m2.degree()) { 
    size_type max_arg_dim=std::max(m1.argument_size(),m2.argument_size());
    for(size_type i=0; i!=max_arg_dim; ++i) {
      size_type m1index=(i<m1.argument_size()) ? m1.index(i) : 0u;
      size_type m2index=(i<m2.argument_size()) ? m2.index(i) : 0u;
      if(m1index > m2index) {
        return true;
      }
      else if(m1index < m2index) {
        return false;
      }
    }
    return false;
  }
  else {
    return m1.degree()<m2.degree();
  }
}


template<class R>
typename Function::Monomial<R>::F
Function::Monomial<R>::evaluate(const LinearAlgebra::Vector<F>& s) const 
{
  ARIADNE_CHECK_ARGUMENT_SIZE(*this,s,"Monomial<R>::evaluate(Point<R>)");
  result_type result=_coefficient;
  for(size_type k=0; k!=this->argument_size(); ++k) {
    result *= pow(s[k],_multi_index[k]);
  }
  return result;
}



template<class R>
typename Function::Polynomial<R>::F
Function::Polynomial<R>::evaluate(const LinearAlgebra::Vector<F>& s) const 
{
  //std::cerr << "Polynomial<R>::evaluate(const LinearAlgebra::Vector<R>& s) const " << std::endl;
  ARIADNE_CHECK_ARGUMENT_SIZE(*this,s,"Polynomial<R>::evaluate(Point<R>)");
  F result=0;
  for(size_type j=0; j!=_terms.size(); ++j) {
    result += _terms[j].evaluate(s);
  }
  return result;
}



template<class R>
LinearAlgebra::Vector<typename Function::PolynomialModel<R>::F>
Function::PolynomialModel<R>::evaluate(const LinearAlgebra::Vector<F>& s) const 
{
  ARIADNE_CHECK_ARGUMENT_SIZE(*this,s,"PolynomialModel<R>::apply(Point<R>)");
  LinearAlgebra::Vector<F> result(this->result_size());
  for(size_type i=0; i!=this->result_size(); ++i) {
    result[i] = _components[i].evaluate(s);
  }
  return result;
}


template<class R>
LinearAlgebra::Matrix< typename Function::PolynomialModel<R>::F >
Function::PolynomialModel<R>::jacobian(const LinearAlgebra::Vector<F>& s) const 
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}






template<class R>
std::ostream&
Function::Monomial<R>::write(std::ostream& os) const
{
  const Monomial<R>& m=*this;
  bool leading_term=true;
  if(m.coefficient()==0) {
    return os << "0";
  }
  if(m.degree()==0) {
    return os << m.coefficient();
  }
  if(m.coefficient()==1) {
  }
  else if(m.coefficient()==-1) {
    os << "-";
    leading_term=false;
  }
  else {
    os << m.coefficient();
    leading_term=false;
  }
  
  for(size_type k=0; k!=m.argument_size(); ++k) {
    if(m.index(k) > 0) {
      if(leading_term) {
        leading_term=false;
      }
      else {
        os << " ";
      }
      os << "x_" << k;
      if(m.index(k)>1) {
        os << "^" << m.index(k);
      }
    }
  }
  return os;
}

template<class R>
std::ostream&
Function::Polynomial<R>::write(std::ostream& os) const
{
  const Polynomial<R>& p=*this;
  if(p.number_of_terms()==0) {
    os << "0";
  }
  else {
    for(size_type j=0; j!=p.number_of_terms(); ++j) {
      const Monomial<R>& m = p.term(j);
      if(m.coefficient()>0 && j>0) {
        os << " + " << m;
      }
      else if(m.coefficient()<0 && j>0) {
        os << " - " << -m;
      }
      else {
        os << m;
      }
    }
  }
  return os;
}

template<class R>
std::ostream&
Function::PolynomialModel<R>::write(std::ostream& os) const
{
  const PolynomialModel<R>& pm=*this;
  return os << pm._components;
}




template<class R>
std::istream&
Function::Monomial<R>::read(std::istream& is)
{
  Monomial<R>& m=*this;
  R a;
  std::vector<size_type> powers;
  uint i;
  uint n;
  char c;
  char d;
  
  is >> c;
  if(isdigit(c)) {
    is.putback(c);
    c='+';
  }
  if(c=='+' || c=='-') {
    is >> d;
    if(isdigit(d)) {
      is.putback(d);
      is >> a;
    }
    else {
      is.putback(d);
      a=1;
    }
    if(c=='-') {
      a=-a;
    }
    is >> c;
  }
  else {
    a=R(1);
  }
  
  while(c=='x' || c=='*') {
    if(c=='x') {
      is.putback(c);
    }
    is >> c;
    if(!c=='x') {
      throw InvalidInput(__PRETTY_FUNCTION__);
    }
    is >> c;
    if(!c=='_') {
      throw InvalidInput(__PRETTY_FUNCTION__);
    }
    is >> i;
    is >> c;
    if(c=='^') {
      is >> n;
    }
    else {
      is.putback(c);
      n=1;
    }
    if(powers.size() <= i) {
      powers.resize(i+1);
    }
    powers[i]+=n;
    is >> c;
  }
  is.putback(c);
  
  m._coefficient=a;
  m._multi_index=array<size_type>(powers.begin(),powers.end());
  
  //if(!is) { std::cerr << "Stream left in invalid state\n" << std::endl; }
  return is;
}

template<class R>
std::istream&
Function::Polynomial<R>::read(std::istream& is)
{
  Polynomial<R>& p=*this;
  p=Polynomial<R>(0);
  Monomial<R> m(0);
  char c;
  is >> c;
  while(is && (c=='x' || isdigit(c) || c=='+' || c=='-')) {
    if(c!='+') {
      is.putback(c);
    }
    is >> m;
    p._terms.push_back(m);
    //std::cerr << "'" << c << "'  " << m << std::endl;
    is >> c;
  }
  is.putback(c);
  
  p._set_argument_size(p._compute_maximum_term_size());
  p._sort();
  
  return is;
}

template<class R>
std::istream&
Function::PolynomialModel<R>::read(std::istream& is)
{
  PolynomialModel<R>& pm=*this;
  std::vector< Polynomial<R> > vec;
  is >> vec;
  pm._components=array< Polynomial<R> >(vec.begin(),vec.end());
  pm._set_argument_size(pm._compute_maximum_component_size());
  return is;
}


}
