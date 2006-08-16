/***************************************************************************
 *            interval_matrix.tpl
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
 
#include <sstream>
#include <string>

#include "declarations.h"


#include "interval_matrix.h"

#include "../numeric/arithmetic.h"
#include "../numeric/interval.h"
#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"
#include "../linear_algebra/interval_vector.h"
#include "../linear_algebra/lu_matrix.h"


namespace Ariadne {
  namespace LinearAlgebra {

    template<typename R>
    IntervalMatrix<R> 
    IntervalMatrix<R>::zero(const size_type r, const size_type c) {
      return IntervalMatrix<R>(r,c); 
    }
    
    template<typename R>
    IntervalMatrix<R> 
    IntervalMatrix<R>::identity(const size_type n) {
      IntervalMatrix<R> result(n,n); 
      for(size_type i=0; i!=n; ++i) { result(i,i)=Interval<R>(1); } 
      return result;
    }
    
  
    
      

    template <typename R>
    IntervalVector<R> 
    prod(const Matrix<R>& A, const IntervalVector<R>& v) {
      IntervalVector<R> result(A.size1());
      for (size_type i=0; i!=result.size(); ++i) {
        result(i)=R(0);
        for (size_type j=0; j!=v.size(); ++j) {
          result(i)+=A(i,j)*v(j);
        }
      }
      return result;
    }
      
    template <typename R>
    IntervalVector<R> 
    prod(const IntervalMatrix<R>& A, const Vector<R>& v) {
      IntervalVector<R> result(A.size1());
      for (size_type i=0; i!=result.size(); ++i) {
        result(i)=R(0);
        for (size_type j=0; j!=v.size(); ++j) {
          result(i)+=A(i,j)*v(j);
        }
      }
      return result;
    }
      
    template <typename R>
    IntervalVector<R> 
    prod(const IntervalMatrix<R>& A, const IntervalVector<R>& v) {
      IntervalVector<R> result(A.size1());
      for (size_type i=0; i!=result.size(); ++i) {
        result(i)=R(0);
        for (size_type j=0; j!=v.size(); ++j) {
          result(i)+=A(i,j)*v(j);
        }
      }
      return result;
    }
      
    template <typename R>
    IntervalMatrix<R> 
    prod(const IntervalMatrix<R>& A, const Matrix<R>& B) {
      IntervalMatrix<R> result(A.size1(),B.size2());
      for (size_type i=0; i!=A.size1(); ++i) {
        for (size_type j=0; j!=B.size2(); ++j) {
          result(i,j)=R(0);
          for (size_type k=0; k!=A.size2(); ++k) {
            result(i,j)+=A(i,k)*B(k,j);
          }
        }
      }
      return result;
    }
      
    template <typename R>
    IntervalMatrix<R> 
    prod(const Matrix<R>& A, const IntervalMatrix<R>& B) {
      IntervalMatrix<R> result(A.size1(),B.size2());
      for (size_type i=0; i!=A.size1(); ++i) {
        for (size_type j=0; j!=B.size2(); ++j) {
          result(i,j)=R(0);
          for (size_type k=0; k!=A.size2(); ++k) {
            result(i,j)+=A(i,k)*B(k,j);
          }
        }
      }
      return result;
    }
    
    template <typename R>
    IntervalMatrix<R> 
    prod(const IntervalMatrix<R>& A, const IntervalMatrix<R>& B) {
      IntervalMatrix<R> result(A.size1(),B.size2());
      for (size_type i=0; i!=A.size1(); ++i) {
        for (size_type j=0; j!=B.size2(); ++j) {
          result(i,j)=R(0);
          for (size_type k=0; k!=A.size2(); ++k) {
            result(i,j)+=A(i,k)*B(k,j);
          }
        }
      }
      return result;
    }

    template<typename R>
    IntervalMatrix<R>
    fprod(const Matrix<typename numerical_traits<R>::field_extension_type>& A, 
         const IntervalMatrix<R>& B) 
    {
      /* FIXME: Write this code! */
      return IntervalMatrix<R>(A.size1(),B.size2());
    }

    template<typename R>
    Matrix<R>
    IntervalMatrix<R>::centre() const
    {
      const IntervalMatrix<R>& A=*this;      
      Matrix<R> result(A.size1(),A.size2());
      for(size_type i=0; i!=A.size1(); ++i) {
        for(size_type j=0; j!=A.size2(); ++j) {
          result(i,j)=A(i,j).centre();
        }
      }
      return result;
    }

    
    template<typename R>
    R
    IntervalMatrix<R>::radius_norm() const
    {
      const IntervalMatrix<R>& A=*this;
      R diameter=0;
      for(size_type i=0; i!=A.size1(); ++i) {
        R row_sum=0;
        for(size_type j=0; j!=A.size2(); ++j) {
          row_sum+=A(i,j).length();
        }
        diameter=max(diameter,row_sum);
      }
      return diameter/2;
    }

    template<typename R>
    IntervalVector<R>
    IntervalMatrix<R>::radius_row_sum() const
    { 
      IntervalVector<R> result(this->size1());
      const IntervalMatrix<R>& self=*this;
      for(dimension_type i=0; i!=self.size1(); ++i) {
        R radius=0;
        for(dimension_type j=0; j!=self.size2(); ++j) {
          radius+=self(i,j).length();
        }
        radius /= 2;
        result[i]=Interval<R>(-radius,radius);
      }
      return result;
    }


    template<typename R>
    R
    IntervalMatrix<R>::upper_norm() const
    {
      return this->norm().upper();
    }
        
    template<typename R>
    R
    IntervalMatrix<R>::upper_log_norm() const
    {
      return this->log_norm().upper();
    }
        
/*
    template<typename R>
    IntervalMatrix<R>
    IntervalMatrix<R>::inverse() const
    {
      LUMatrix< Interval<R> > lu(*this);
      IntervalMatrix<R> luinv=lu.inverse();
      return luinv;
     
      throw std::domain_error("IntervalMatrix<R>::inverse() const not implemented.");
    }
*/
    
    
    template<typename R>
    Matrix<R>
    over_approximation(const IntervalMatrix<R>& A)
    {
      
      typedef typename numerical_traits<R>::field_extension_type F;

      assert(A.size1()==A.size2());
      dimension_type n=A.size1();
      
      IntervalMatrix<F> AF(n,n);
      for(size_type i=0; i!=n; ++i) {
        for(size_type j=0; j!=n; ++j) {
          AF(i,j)=Interval<F>(A(i,j).lower(),A(i,j).upper());
        }
      }

      Matrix<R> Amid(n,n);
      for(size_type i=0; i!=n; ++i) {
        for(size_type j=0; j!=n; ++j) {
          Amid(i,j)=(A(i,j).upper()+A(i,j).lower())/2;
        }
      }
      
      R Aaccuracy=upper_norm(IntervalMatrix<R> (A-Amid));
      R approx_accuracy=Aaccuracy/1024;
      
      Matrix<F> Amidinv=inverse(Amid);
      IntervalMatrix<F> E=Amidinv*AF;
      
      F excess=norm(E).upper();
      R approx_excess=Ariadne::approximate<R>(excess,approx_accuracy)+approx_accuracy;
      
      return approx_excess*Amid;
    }
    
    template<typename R>
    IntervalMatrix<R>
    approximate(const Matrix<typename numerical_traits<R>::field_extension_type>& A, const R& e)
    {
      R abserr=e/(2*A.number_of_columns());
      Interval<R> err(-abserr,abserr);
      IntervalMatrix<R> result(A.size1(),A.size2());
      for(size_type i=0; i!=result.size1(); ++i) {
        for(size_type j=0; j!=result.size2(); ++j) {
          const R& Aij=A(i,j);
          result(i,j)=err+Aij;
        }
      }
      return result;
    }
      
     
/*
    template<typename R>
    IntervalMatrix<R>
    identity_IntervalMatrix(const size_type& n) 
    {
      IntervalMatrix<R> result(n,n);
      for(size_type i=0; i!=n; ++i) {
        result(i,i)=R(1);
      }
      return result;
    }
*/
    
 
    
  }
}
