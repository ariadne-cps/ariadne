/***************************************************************************
 *            tensor.code.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it Pieter.Collins@cwi.nl
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
 
/*! \file tensor.h
 *  \brief Tensors.
 */

#include "tensor.h"

#include <iostream>
#include <sstream>
#include <string>

#include "../linear_algebra/multi_index.h"
#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

namespace Ariadne {
  namespace LinearAlgebra {

    template<class R>
    void
    Tensor<R>::_instantiate()
    {
      typedef typename Numeric::traits<R>::arithmetic_type F;
      R* s=0; 
      Tensor<R>* T=0; 
      Tensor<F>* qT=0;
      
      *qT = *T + *T;
      *qT = *T * *s;
    }
    
    
    template<class R>
    Tensor<typename Numeric::traits<R>::arithmetic_type>
    operator+(const Tensor<R>& T1, const Tensor<R>& T2)
    {
      typedef typename Numeric::traits<R>::arithmetic_type F;
      if(T1.sizes()!=T2.sizes()) { throw IncompatibleSizes(__PRETTY_FUNCTION__); }
      Tensor<F> result(T1.sizes());
      for(size_type i=0; i!=result.number_of_elements(); ++i) {
        result.begin()[i]=T1.begin()[i]+T2.begin()[i];
      }
      return result;
    }
      
    
    template<class R>
    Tensor<typename Numeric::traits<R>::arithmetic_type>
    operator*(const Tensor<R>& T, const R& s)
    {
      typedef typename Numeric::traits<R>::arithmetic_type F;
      Tensor<F> result(T.sizes());
      for(size_type i=0; i!=result.number_of_elements(); ++i) {
        result.begin()[i]=T.begin()[i]*s;
      }
      return result;
    }
    
    
    
    template<class R>
    void
    SymmetricTensor<R>::_instantiate()
    {
      R* s=0; SymmetricTensor<R>* T=0;
      *T + *T; *T * *s;
    }
    
    
    template<class R>
    SymmetricTensor<typename Numeric::traits<R>::arithmetic_type>
    operator+(const SymmetricTensor<R>& T1, const SymmetricTensor<R>& T2)
    {
      typedef typename Numeric::traits<R>::arithmetic_type F;
      if(!(T1.argument_size()==T2.argument_size() && T1.degree()==T2.degree())) {
        throw IncompatibleSizes(__PRETTY_FUNCTION__); 
      }
      SymmetricTensor<F> result(T1.argument_size(),T1.degree());
      for(size_type i=0; i!=result.number_of_independent_elements(); ++i) {
        result.begin()[i]=T1.begin()[i]+T2.begin()[i];
      }
      return result;
    }
      
    
    template<class R>
    SymmetricTensor<typename Numeric::traits<R>::arithmetic_type>
    operator*(const SymmetricTensor<R>& T, const R& s)
    {
      typedef typename Numeric::traits<R>::arithmetic_type F;
      SymmetricTensor<F> result(T.argument_size(),T.degree());
      for(size_type i=0; i!=result.number_of_elements(); ++i) {
        result.begin()[i]=T.begin()[i]*s;
      }
      return result;
    }
    
    
    
    template<class R>
    void
    DerivativeTensor<R>::_instantiate()
    {
      R* s=0; 
      Vector<R>* v=0; 
      Matrix<R>* A=0; 
      DerivativeTensor<R>* T=0;
      *T + *T;
      *T - *T;
      *T * *s;
      *T * *v;
      *T * *A;
      *T * *T;
    }
    
    
    template<class R>
    DerivativeTensor<R>::DerivativeTensor(const Vector<R>& v) 
      : _res_size(v.size()), _arg_size(1), _degree(0), _elements(v.size())
    {
      for(size_type i=0; i!=this->_res_size; ++i) {
        this->_elements[i]=v(i);
      }
    }
    
    
    template<class R>
    DerivativeTensor<R>::DerivativeTensor(const Matrix<R>& A) 
      : _res_size(A.number_of_rows()), _arg_size(A.number_of_columns()), 
        _degree(1), _elements(A.number_of_rows()*A.number_of_columns())
    {
      for(size_type i=0; i!=this->_res_size; ++i) {
        for(size_type j=0; i!=this->_arg_size; ++j) {
          this->_elements[i*this->_res_size+j]=A(i,j);
        }
      }
    }
    
    
    template<class R>
    DerivativeTensor<typename Numeric::traits<R>::arithmetic_type>
    operator+(const DerivativeTensor<R>& T1, const DerivativeTensor<R>& T2) 
    {
      if(!(T1.result_size()==T2.result_size() && T1.argument_size()!=T2.argument_size() &&
           T1.degree()==T2.degree())) 
      {
        throw IncompatibleSizes(__PRETTY_FUNCTION__);
      }
      DerivativeTensor<typename Numeric::traits<R>::arithmetic_type> result(T1.result_size(),T1.argument_size(),T1.degree());
      for(size_type i=0; i!=result.number_of_independent_elements(); ++i) {
        result.begin()[i]=T1.begin()[i]+T2.begin()[i];
      }
      return result;
    }
      
    
    template<class R>
    DerivativeTensor<typename Numeric::traits<R>::arithmetic_type>
    operator-(const DerivativeTensor<R>& T1, const DerivativeTensor<R>& T2) 
    {
      if(!(T1.result_size()==T2.result_size() && T1.argument_size()==T2.argument_size() &&
           T1.degree()==T2.degree()))
      {
        throw IncompatibleSizes(__PRETTY_FUNCTION__);
      }
      DerivativeTensor<typename Numeric::traits<R>::arithmetic_type> result(T1.result_size(),T1.argument_size(),T1.degree());
      for(size_type i=0; i!=result.number_of_independent_elements(); ++i) {
        result.begin()[i]=T1.begin()[i]-T2.begin()[i];
      }
      return result;
    }
      
    
    template<class R>
    DerivativeTensor<typename Numeric::traits<R>::arithmetic_type>
    operator*(const DerivativeTensor<R>& T, const R& s) 
    {
      DerivativeTensor<typename Numeric::traits<R>::arithmetic_type> result(T.result_size(),T.argument_size(),T.degree());
      for(size_type i=0; i!=result.number_of_independent_elements(); ++i) {
        result.begin()[i]=T.begin()[i]*s;
      }
      return result;
    }
      
      
    template<class R>
    DerivativeTensor<typename Numeric::traits<R>::arithmetic_type>
    operator*(const DerivativeTensor<R>& T, const Vector<R>& v) 
    {
      return operator*(T,DerivativeTensor<R>(v));
    }
      
    
    template<class R>
    DerivativeTensor<typename Numeric::traits<R>::arithmetic_type>
    operator*(const DerivativeTensor<R>& T, const Matrix<R>& A) 
    {
      return operator*(T,DerivativeTensor<R>(A));
    }
      
    
    template<class R>
    DerivativeTensor<typename Numeric::traits<R>::arithmetic_type>
    operator*(const DerivativeTensor<R>& T1, const DerivativeTensor<R>& T2) 
    {
      //std::cerr << __PRETTY_FUNCTION__ << std::endl;
      if(T1.argument_size()!=T2.result_size()) { throw IncompatibleSizes(__PRETTY_FUNCTION__); }
      if(T1.degree()==0) { throw IncompatibleSizes(__PRETTY_FUNCTION__); }
      
      DerivativeTensor<typename Numeric::traits<R>::arithmetic_type> T0(T1.result_size(),T2.argument_size(),T1.degree()+T2.degree()-1);

      MultiIndex m0(T0.argument_size());
      MultiIndex m1(T1.argument_size());
      MultiIndex m2(T2.argument_size());

      m1.set_index(0,T1.degree()-1);
      MultiIndexIterator mi_begin1(m1);
      m1.increment_index(0);
      MultiIndexIterator mi_end1(m1);

      m2.set_index(0,T2.degree());
      MultiIndexIterator mi_begin2(m2);
      m2.increment_index(0);
      MultiIndexIterator mi_end2(m2);

	  unsigned int zn=3;
	  R zx=2;
	  R zy=5;
	  R zr;
	  zr=zn*zx;
	  zr+=zn*(zx*zy);
		
      for(size_type i1=0; i1!=T1.result_size(); ++i1) {
        for(size_type i2=0; i2!=T2.result_size(); ++i2) {
          for(MultiIndexIterator mi_iter2=mi_begin2; mi_iter2!=mi_end2; ++mi_iter2) {
            for(MultiIndexIterator mi_iter1=mi_begin1; mi_iter1!=mi_end1;++mi_iter1) {
              m1=*mi_iter1;
              m2=*mi_iter2;
              m0=m1+m2;
              size_type n=m1.number()*m2.number();
              m1.increment_index(i2);
              //T0(i1,m0)+=n*T1(i1,m1)*T2(i2,m2);
            }
          }
        }
      }
      return T0;
    }

    

    
    template<class R>
    std::ostream&
    Tensor<R>::write(std::ostream& os) const
    {
      return os << "Tensor(...)";
    }
    
    
    template<class R>
    std::ostream&
    SymmetricTensor<R>::write(std::ostream& os) const
    {
      //std::cerr << __PRETTY_FUNCTION__ << std::endl;
      size_type n=this->argument_size();
      size_type d=this->degree();
      
      os << "SymmetricTensor( argument_size=" << n << ", degree=" << d << "\n  elements=[\n";
      for(MultiIndexIterator j(n,d); j!=MultiIndexIterator(n,d+1); ++j) {
        os << "    " << *j << ":" << this->position(*j) << ": " << (*this)(*j) << "; \n"; 
      }
      os << "  ]\n)\n";
      return os;
    }

    
    template<class R>
    std::ostream&
    DerivativeTensor<R>::write(std::ostream& os) const
    {
      //std::cerr << __PRETTY_FUNCTION__ << std::endl;
      size_type m=this->result_size();
      size_type n=this->argument_size();
      size_type d=this->degree();
      
      os << "DerivativeTensor( result_size=" << m << ", argument_size=" << n << ", degree=" << d << "\n  elements=[\n";
      for(size_type i=0; i!=m; ++i) {
        for(MultiIndexIterator j(n,d); j!=MultiIndexIterator(n,d+1); ++j) {
          os << "    " << i << "," << *j << ":" << this->position(i,*j) << ": " << (*this)(i,*j) << "; \n"; 
        }
      }
      os << "  ]\n)\n";
      return os;
    }

  }
}
