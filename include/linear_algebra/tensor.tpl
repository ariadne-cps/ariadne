/***************************************************************************
 *            tensor.tpl
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

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

namespace Ariadne {
  namespace LinearAlgebra {

    template<typename R>
    vector<R>
    product(const tensor<R>& T, const vector<R>& v1, const vector<R>& v2) 
    {
      assert(T.size(1)==v1.size());
      assert(T.size(2)==v2.size());
      
      vector<R> result(T.size(0));
      for(size_type i=0; i!=T.size(0); ++i) {
        for(size_type j=0; j!=T.size(1); ++j) {
          for(size_type k=0; k!=T.size(2); ++k) {
            result(i)+=T(i,j,k)*v1(j)*v2(k);
          }
        }
      }
      return result;
    }
      
    template<typename R>
    matrix<R>
    product(const tensor<R>& T, const vector<R>& v1) 
    {
      assert(T.size(1)==v1.size());
      
      matrix<R> result(T.size(0),T.size(1));
      for(size_type i=0; i!=T.size(0); ++i) {
        for(size_type j=0; j!=T.size(1); ++j) {
          for(size_type k=0; k!=T.size(2); ++k) {
            result(i,k)+=T(i,j,k)*v1(j);
          }
        }
      }
      return result;
    }

    template<typename R>
    tensor<R>
    product(const tensor<R>& T, const matrix<R>& A1) 
    {
      assert(T.size(1)==A1.size1());
      
      tensor<R> result(T.size(0),A1.size2(),T.size(2));
      for(size_type i=0; i!=result.size(0); ++i) {
        for(size_type j=0; j!=result.size(1); ++j) {
          for(size_type k=0; k!=result.size(2); ++k) {
            for(size_type l=0; l!=T.size(1); ++l) {
              result(i,j,k)+=T(i,l,k)*A1(l,j);
            }
          }
        }
      }
      return result;
    }

    template<typename R>
    tensor<R>
    product(const tensor<R>& T, const matrix<R>& A1, const matrix<R>& A2) 
    {
      assert(T.size(1)==A1.size1());
      assert(T.size(2)==A2.size1());
      
      tensor<R> result(T.size(0),A1.size2(),A2.size2());
      for(size_type i=0; i!=T.size(0); ++i) {
        for(size_type j=0; j!=A1.size2(); ++j) {
          for(size_type k=0; k!=A2.size2(); ++k) {
            for(size_type l=0; l!=T.size(1); ++l) {
              for(size_type m=0; k!=T.size(2); ++m) {
                result(i,j,k)+=T(i,l,m)*A1(l,j)*A2(m,k);
              }
            }
          }
        }
      }
      return result;
    }
    
    template<typename R>
    std::ostream&
    operator<< (std::ostream& os, const tensor<R>& T)
    {
      return os;
    }
    
  }
}
