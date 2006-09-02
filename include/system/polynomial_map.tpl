/***************************************************************************
 *            polynomial_map.tpl
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
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

#include "polynomial_map.h"

#include <string>
#include <sstream>

#include <list>
#include <set>
#include <vector>
#include <valarray>

#include "../utility/stlio.h"
#include "../base/array.h"
#include "../numeric/interval.h"

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"
#include "../linear_algebra/interval_matrix.h"
#include "../linear_algebra/interval_vector.h"

#include "../geometry/point.h"
#include "../geometry/rectangle.h"


namespace Ariadne {
  namespace System {

    template<typename R>
    Monomial<R>::Monomial(const std::string& s) 
    {
      std::stringstream ss(s);
      ss >> *this;
    }
    
    template<typename R>
    Polynomial<R>::Polynomial(const std::string& s) 
    {
      //std::cerr << "Polynomial(\"" << s << "\")" << std::endl;
      std::stringstream ss(s);
      ss >> *this;
    }
    
    
    template<typename R>
    PolynomialMap<R>::PolynomialMap(const std::string& s) 
    {
      //std::cerr << "PolynomialMap(\"" << s << "\")" << std::endl;
      std::stringstream ss(s);
      ss >> *this;
    }
    
    
    template<typename R>
    size_type
    Monomial<R>::degree() const 
    {
      size_type result=0;
      for(dimension_type i=0; i!=this->argument_dimension(); ++i) {
        result+=this->index(i);
      }
      return result;
    }
    
    template<typename R>
    void
    Monomial<R>::_resize(const dimension_type& n) 
    {
      this->_multi_index.resize(n);
    }
    
    template<typename R>
    void
    Polynomial<R>::_sort() 
    {
      std::sort(this->_terms.begin(), this->_terms.end());

      typedef typename std::vector< Monomial<R> >::const_iterator terms_const_iterator;
      std::vector< Monomial<R> > new_terms;
      for(terms_const_iterator i=this->_terms.begin(); i!=this->_terms.end(); ++i) {
        if(!new_terms.empty() && new_terms.back().multi_index() == i->multi_index()) {
          new_terms.back()._coefficient += i->coefficient();
        }
        else if(i->_coefficient!=0) {
          new_terms.push_back(*i);
        }
      }
      this->_terms.swap(new_terms);
    }
    
    template<typename R>
    void
    Polynomial<R>::_set_argument_dimension(const dimension_type& n) 
    {
      this->_argument_dimension=n;
      for(size_type i=0; i!=this->_terms.size(); ++i) {
        this->_terms[i]._resize(n);
      }
    }
    
    template<typename R>
    dimension_type
    Polynomial<R>::_compute_maximum_term_dimension() const
    {
      dimension_type result=0;
      for(size_type i=0; i!=this->_terms.size(); ++i) {
        result=max(this->_terms[i].argument_dimension(),result);
      }
      return result;
    }
    
    
    template<typename R>
    void
    PolynomialMap<R>::_set_argument_dimension(const dimension_type& n) 
    {
      this->_argument_dimension=n;
      for(size_type i=0; i!=this->_components.size(); ++i) {
        this->_components[i]._set_argument_dimension(n);
      }
    }
    
    template<typename R>
    dimension_type
    PolynomialMap<R>::_compute_maximum_component_dimension() const
    {
      dimension_type result=0;
      for(size_type i=0; i!=this->_components.size(); ++i) {
        result=max(this->_components[i].argument_dimension(),result);
      }
      return result;
    }
    
    
    
    template<typename R> 
    bool 
    operator<(const Monomial<R>& m1, const Monomial<R>& m2)
    {
      if(m1.degree() == m2.degree()) { 
        dimension_type max_arg_dim=max(m1.argument_dimension(),m2.argument_dimension());
        for(dimension_type i=0; i!=max_arg_dim; ++i) {
          size_type m1index=(i<m1.argument_dimension()) ? m1.index(i) : 0u;
          size_type m2index=(i<m2.argument_dimension()) ? m2.index(i) : 0u;
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
    
    
    template<typename R>
    R
    Monomial<R>::apply(const Geometry::Point<R>& s) const 
    {
      assert(s.dimension() == this->argument_dimension());
      real_type result=_coefficient;
      for(size_type k=0; k!=this->argument_dimension(); ++k) {
        result *= pow(s[k],_multi_index[k]);
      }
      return result;
    }
    
    template<typename R>
    Interval<R>
    Monomial<R>::apply(const Geometry::Rectangle<R>& r) const 
    {
      assert(r.dimension() == this->argument_dimension());
      Interval<R> result(this->_coefficient);
      for(size_type k=0; k!=this->argument_dimension(); ++k) {
        result *= pow(r[k],_multi_index[k]);
      }
      return result;
    }
    
    
    template<typename R>
    R
    Polynomial<R>::apply(const Geometry::Point<R>& s) const 
    {
      //std::cerr << "Polynomial<R>::apply(const Geometry::Point<R>& s) const " << std::endl;
      assert(s.dimension() == this->argument_dimension());
      real_type result=0;
      for(size_type j=0; j!=_terms.size(); ++j) {
        result += _terms[j].apply(s);
      }
      return result;
    }
    
    template<typename R>
    Interval<R>
    Polynomial<R>::apply(const Geometry::Rectangle<R>& r) const 
    {
      assert(r.dimension() == this->argument_dimension());
      Interval<R> result(R(0));
      for(size_type j=0; j!=_terms.size(); ++j) {
        result += _terms[j].apply(r);
      }
      return result;
    }
    
    
    template<typename R>
    Geometry::Point<R>
    PolynomialMap<R>::apply(const Geometry::Point<R>& s) const 
    {
      assert(s.dimension() == this->argument_dimension());
      state_type result(this->result_dimension());
      for(size_type i=0; i!=this->result_dimension(); ++i) {
        result[i] = _components[i].apply(s);
      }
      return result;
    }
    
    template<typename R>
    Geometry::Rectangle<R>
    PolynomialMap<R>::apply(const Geometry::Rectangle<R>& r) const 
    {
      assert(r.dimension() == this->argument_dimension());
      Geometry::Rectangle<R> result(this->result_dimension());
      for(size_type i=0; i!=this->result_dimension(); ++i) {
        result[i] = _components[i].apply(r);
      }
      return result;
    }
    
    template<typename R>
    LinearAlgebra::Matrix<R>
    PolynomialMap<R>::derivative(const Geometry::Point<R>& s) const 
    {
      //std::cerr << "PolynomialMap<R>::derivative(const Geometry::Point<R>& s) const " << std::endl;
      this->_compute_derivative();
      LinearAlgebra::Matrix<R> result(this->result_dimension(), this->argument_dimension());
      for(size_type i=0; i!=this->result_dimension(); ++i) {
        for(size_type j=0; j!=this->argument_dimension(); ++j) {
          result(i,j)=this->_derivative(i,j).apply(s);
        }
      }
      return result;
    }
    
    
    template<typename R>
    LinearAlgebra::IntervalMatrix<R>
    PolynomialMap<R>::derivative(const Geometry::Rectangle<R>& r) const 
    {
      this->_compute_derivative();
      LinearAlgebra::IntervalMatrix<R> result(this->result_dimension(), this->argument_dimension());
      for(size_type i=0; i!=this->result_dimension(); ++i) {
        for(size_type j=0; j!=this->argument_dimension(); ++j) {
          result(i,j)=this->_derivative(i,j).apply(r);
        }
      }
      return result;
    }
    
    template<typename R>
    const PolynomialMatrix<R>&
    PolynomialMap<R>::derivative() const 
    {
      this->_compute_derivative(); return this->_derivative;
    }
   
  
    template<typename R>
    void
    PolynomialMap<R>::_compute_derivative() const 
    {
      if(this->_derivative.number_of_rows()!=0) {
        return;
      }
      
      this->_derivative=PolynomialMatrix<R>(this->result_dimension(),this->argument_dimension());
      for(dimension_type i=0; i!=this->result_dimension(); ++i) {
        for(dimension_type j=0; j!=this->argument_dimension(); ++j) {
          Polynomial<R>& dp=this->_derivative._matrix(i,j);
          dp._argument_dimension=this->argument_dimension();
          for(size_type k=0; k!=this->_components[i].number_of_terms(); ++k) {
            const Monomial<R>& m=this->_components[i]._terms[k];
            Monomial<R> dm=m;
            dm._coefficient *= m._multi_index[j];
            dm._multi_index[j]-=1;
            dp._terms.push_back(dm);
          }
          dp._sort();
        }
      }
    }
    

    
    template <typename R>
    std::ostream&
    operator<<(std::ostream& os, const Monomial<R>& m) 
    {
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
      
      for(size_type k=0; k!=m.argument_dimension(); ++k) {
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

    template <typename R>
    std::ostream&
    operator<<(std::ostream& os, const Polynomial<R>& p) 
    {
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

    template <typename R>
    std::ostream&
    operator<<(std::ostream& os, const PolynomialMap<R>& p) 
    {
      return os << p._components;
    }
    
    template <typename R>
    std::ostream&
    operator<<(std::ostream& os, const PolynomialMatrix<R>& p) 
    {
      os << "[";
      for(size_type i=0; i!=p.number_of_rows(); ++i) {
        if(i!=0) {
          os << "; ";
        }
        for(size_type j=0; j!=p.number_of_columns(); ++j) {
          if(j!=0) { 
            os << ", ";
          }
          os << p(i,j);
        }
      }
      os << "]";
      return os;
    }
    
    
    
    
    template <typename R>
    std::istream&
    operator>>(std::istream& is, Monomial<R>& m) 
    {
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
          throw std::runtime_error("Invalid format for monomial");
        }
        is >> c;
        if(!c=='_') {
          throw std::runtime_error("Invalid format for monomial");
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

    template <typename R>
    std::istream&
    operator>>(std::istream& is, Polynomial<R>& p) 
    {
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
      
      p._set_argument_dimension(p._compute_maximum_term_dimension());
      p._sort();
      
      return is;
    }
    
    template <typename R>
    std::istream&
    operator>>(std::istream& is, PolynomialMap<R>& pm) 
    {
      std::vector< Polynomial<R> > vec;
      is >> vec;
      pm._components=array< Polynomial<R> >(vec.begin(),vec.end());
      pm._set_argument_dimension(pm._compute_maximum_component_dimension());
      return is;
    }
    
  
  }
}
