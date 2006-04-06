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
 
#include "../linear_algebra/interval_matrix.h"


namespace Ariadne {
  namespace LinearAlgebra {

    template <>
    std::ostream&
    operator<<(std::ostream& os, const interval_matrix<Float64>& A)
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
    operator<<(std::ostream& os, const interval_matrix<Dyadic>& A)
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
    operator<<(std::ostream& os, const interval_matrix<Rational>& A)
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
    interval_vector<R> 
    prod(const matrix<R>& A, const interval_vector<R>& v) {
      interval_vector<R> result(A.size1());
      for (size_type i=0; i!=result.size(); ++i) {
        result(i)=Interval<R>(0);
        for (size_type j=0; j!=v.size(); ++j) {
          result(i)+=A(i,j)*v(j);
        }
      }
      return result;
    }
      
    template <typename R>
    interval_vector<R> 
    prod(const interval_matrix<R>& A, const vector<R>& v) {
      interval_vector<R> result(A.size1());
      for (size_type i=0; i!=result.size(); ++i) {
        result(i)=Interval<R>(0);
        for (size_type j=0; j!=v.size(); ++j) {
          result(i)+=A(i,j)*v(j);
        }
      }
      return result;
    }
      
    template <typename R>
    interval_vector<R> 
    prod(const interval_matrix<R>& A, const interval_vector<R>& v) {
      interval_vector<R> result(A.size1());
      for (size_type i=0; i!=result.size(); ++i) {
        result(i)=Interval<R>(0);
        for (size_type j=0; j!=v.size(); ++j) {
          result(i)+=A(i,j)*v(j);
        }
      }
      return result;
    }
      
    template <typename R>
    interval_matrix<R> 
    prod(const interval_matrix<R>& A, const matrix<R>& B) {
      interval_matrix<R> result(A.size1(),B.size2());
      for (size_type i=0; i!=A.size1(); ++i) {
        for (size_type j=0; j!=B.size2(); ++j) {
          result(i,j)=Interval<R>(0);
          for (size_type k=0; k!=A.size2(); ++k) {
            result(i,j)+=A(i,k)*B(k,j);
          }
        }
      }
      return result;
    }
      
    template <typename R>
    interval_matrix<R> 
    prod(const matrix<R>& A, const interval_matrix<R>& B) {
      interval_matrix<R> result(A.size1(),B.size2());
      for (size_type i=0; i!=A.size1(); ++i) {
        for (size_type j=0; j!=B.size2(); ++j) {
          result(i,j)=Interval<R>(0);
          for (size_type k=0; k!=A.size2(); ++k) {
            result(i,j)+=A(i,k)*B(k,j);
          }
        }
      }
      return result;
    }
    
    template <typename R>
    interval_matrix<R> 
    prod(const interval_matrix<R>& A, const interval_matrix<R>& B) {
      interval_matrix<R> result(A.size1(),B.size2());
      for (size_type i=0; i!=A.size1(); ++i) {
        for (size_type j=0; j!=B.size2(); ++j) {
          result(i,j)=Interval<R>(0);
          for (size_type k=0; k!=A.size2(); ++k) {
            result(i,j)+=A(i,k)*B(k,j);
          }
        }
      }
      return result;
    }

    template<typename R>
    interval_matrix<R>
    fprod(const matrix<typename numerical_traits<R>::field_extension_type>& A, 
         const interval_matrix<R>& B) 
    {
      /* FIXME: Write this code! */
      return interval_matrix<R>(A.size1(),B.size2());
    }
          
    template<typename R>
    matrix<R>
    centre(const interval_matrix<R>& A)
    {
      dimension_type m=A.size1();
      dimension_type n=A.size2();
      
      matrix<R> result(m,n);
      for(size_type i=0; i!=m; ++i) {
        for(size_type j=0; j!=n; ++j) {
          result(i,j)=(A(i,j).lower()+A(i,j).upper())/2;
        }
      }
      return result;
    }

    template<typename R>
    R
    radius(const interval_matrix<R>& A)
    {
      R diameter=0;
      for(size_type i=0; i!=A.size1(); ++i) {
        R row_sum=0;
        for(size_type j=0; j!=A.size2(); ++j) {
          row_sum+=A(i,j).upper()-A(i,j).lower();
        }
        diameter=max(diameter,row_sum);
      }
      return diameter/2;
    }



    template<typename R>
    Interval<R>
    norm(const interval_matrix<R>& A) 
    {
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
    upper_norm(const interval_matrix<R>& A) 
    {
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
    upper_log_norm(const interval_matrix<R>& A) 
    {
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
    matrix<R>
    over_approximation(const interval_matrix<R>& A)
    {
      
      typedef typename numerical_traits<R>::field_extension_type F;

      assert(A.size1()==A.size2());
      dimension_type n=A.size1();
      
      interval_matrix<F> AF(n,n);
      for(size_type i=0; i!=n; ++i) {
        for(size_type j=0; j!=n; ++j) {
          AF(i,j)=Interval<F>(A(i,j).lower(),A(i,j).upper());
        }
      }

      matrix<R> Amid(n,n);
      for(size_type i=0; i!=n; ++i) {
        for(size_type j=0; j!=n; ++j) {
          Amid(i,j)=(A(i,j).upper()+A(i,j).lower())/2;
        }
      }
      
      R Aaccuracy=upper_norm(interval_matrix<R> (A-Amid));
      R approx_accuracy=Aaccuracy/1024;
      
      matrix<F> Amidinv=inverse(Amid);
      interval_matrix<F> E=Amidinv*AF;
      
      F excess=norm(E).upper();
      R approx_excess=Ariadne::approximate<R>(excess,approx_accuracy)+approx_accuracy;
      
      return approx_excess*Amid;
    }
    
    interval_matrix<Dyadic>
    approximate(const matrix<Rational>& A, const Dyadic& e)
    {
      Dyadic err=e/2;
      interval_matrix<Dyadic> result(A.size1(),A.size2());
      for(size_type i=0; i!=result.size1(); ++i) {
        for(size_type j=0; j!=result.size1(); ++j) {
          Dyadic lower=Ariadne::approximate<Dyadic>(A(i,j),err)-err;
          Dyadic upper=Ariadne::approximate<Dyadic>(A(i,j),err)+err;
          result(i,j)=Interval<Dyadic>(lower,upper);
        }
      }
      return result;
    }
      
      
  }
}
