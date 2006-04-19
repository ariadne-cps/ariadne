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
 
#include "interval_matrix.h"

#include "../numeric/arithmetic.h"
#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"
#include "../linear_algebra/interval_vector.h"


namespace Ariadne {
  namespace LinearAlgebra {

    template <typename R>
    std::ostream&
    operator<<(std::ostream& os, const IntervalMatrix<R>& A)
    {
      if(A.size1()==0 || A.size2()==0) {
        return os << "[ ]";
      }
      
      for(uint i=0; i!=A.size1(); ++i) {
        for(uint j=0; j!=A.size2(); ++j) {
          os << (j==0 ? (i==0 ? "[ " : "; ") : ",");
          os << A(i,j);
        }
      }
      os << " ]";
      return os;
    }
     
    template <>
    std::ostream&
    operator<<(std::ostream& os, const IntervalMatrix<Rational>& A)
    {
      if(A.size1()==0 || A.size2()==0) {
        return os << "[ ]";
      }
      
      for(uint i=0; i!=A.size1(); ++i) {
        for(uint j=0; j!=A.size2(); ++j) {
          os << (j==0 ? (i==0 ? "[ " : "; ") : ",");
          double l=Ariadne::convert_to<double>(A(i,j).lower());
          double u=Ariadne::convert_to<double>(A(i,j).upper());
          os << Ariadne::Interval<double>(l,u);
        }
      }
      os << " ]";
      return os;
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
    IntervalMatrix<R>::IntervalMatrix(const Matrix<R>& A, const R& r)
      : Base(A.size1(),A.size2()) 
    { 
      for(size_type i=0; i!=A.size1(); ++i) {
        for(size_type j=0; j!=A.size2(); ++j) {
          Base::operator()(i,j)=Interval<R>(A(i,j)-r,A(i,j)+r);
        }
      }
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
    IntervalMatrix<R>::radius() const
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
    Interval<R>
    IntervalMatrix<R>::norm() const 
    {
      const IntervalMatrix<R>& A=*this;
      R lower_bound=0;
      R upper_bound=0;
      for(size_type i=0; i!=A.size1(); ++i) {
        R lower_row_sum=0;
        R upper_row_sum=0;
        for(size_type j=0; j!=A.size2(); ++j) {
          if(!(A(i,j).lower()<=0 && A(i,j).upper()>=0)) {
            lower_row_sum+=min( abs(A(i,j).lower()), abs(A(i,j).upper()) );
          }
          upper_row_sum+=max( abs(A(i,j).lower()), abs(A(i,j).upper()) );
        }
        lower_bound=max(lower_bound,lower_row_sum);
        upper_bound=max(upper_bound,upper_row_sum);
      }
      return Interval<R>(lower_bound,upper_bound);
    }
        
    template<typename R>
    R
    IntervalMatrix<R>::upper_norm() const
    {
      const IntervalMatrix<R>& A=*this;
      R upper_bound=0;
      for(size_type i=0; i!=A.size1(); ++i) {
         R upper_row_sum=0;
        for(size_type j=0; j!=A.size2(); ++j) {
          upper_row_sum+=max( abs(A(i,j).lower()), abs(A(i,j).upper()) );
        }
        upper_bound=max(upper_bound,upper_row_sum);
      }
      return upper_bound;
    }
        
    template<typename R>
    R
    IntervalMatrix<R>::upper_log_norm() const
    {
      const IntervalMatrix<R>& A=*this;
      R upper_bound=0;
      for(size_type i=0; i!=A.size1(); ++i) {
         R upper_row_sum=A(i,i).upper();
        for(size_type j=0; j!=A.size2(); ++j) {
          if(i!=j) {
           upper_row_sum+=max( abs(A(i,j).lower()), abs(A(i,j).upper()) );
          }
         }
        upper_bound=max(upper_bound,upper_row_sum);
      }
      return upper_bound;
    }
        

    

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
      R err=e/2;
      IntervalMatrix<R> result(A.size1(),A.size2());
      for(size_type i=0; i!=result.size1(); ++i) {
        for(size_type j=0; j!=result.size1(); ++j) {
          R lower=Ariadne::approximate<R>(A(i,j),err)-err;
          R upper=Ariadne::approximate<R>(A(i,j),err)+err;
          result(i,j)=Interval<R>(lower,upper);
        }
      }
      return result;
    }
      
    template<typename R>
    IntervalMatrix<R>
    identity_IntervalMatrix(const size_type& n) 
    {
      IntervalMatrix<R> result(n,n);
      for(size_type i=0; i!=n; ++i) {
        result(i,i)=1;
      }
      return result;
    }
    
    template<typename R>
    IntervalMatrix<R>
    exp(const IntervalMatrix<R>& A) 
    {
      assert(A.size1()==A.size2());
      R err=A.radius()/65536;
      if(err==0) {
        err=A.upper_norm()/65536;
        err/=65536;
        err/=65536;
      }
            
      IntervalMatrix<R> result=identity_IntervalMatrix<R>(A.size1())+A;
      IntervalMatrix<R> term=A;
      unsigned int n=1;
      while(term.upper_norm()>err) {
        n=n+1;
        term=(term*A)/Interval<R>(R(n));
        result+=term;
      }
      term=Interval<R>(-1,1)*term;
      result+=term;
      
      return result;
    }
    
    
  }
}
