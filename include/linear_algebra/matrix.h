/***************************************************************************
 *            matrix.h
 *
 *  Mon May  3 12:31:15 2004
 *  Copyright  2004-6  Alberto Casagrande, Pieter Collins
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
 
/*! \file matrix.h
 *  \brief Matrices and Matrix operations.
 */

#ifndef ARIADNE_MATRIX_H
#define ARIADNE_MATRIX_H

#include <iosfwd>

#include "../base/types.h"
#include "../base/array.h"
#include "../numeric/numerical_traits.h"
#include "../numeric/integer.h"
#include "../numeric/interval.h"

#include "../linear_algebra/vector.h"
#include "../linear_algebra/exceptions.h"

namespace Ariadne {
  namespace LinearAlgebra {

    /*!\brief Base class for all matrix expressions. */
    template<class E>
    class MatrixExpression 
    {
     public:
      /*!\brief Convert \a *this to a reference to E. */
      E& operator() () { return static_cast<E&>(*this); }
      /*!\brief Convert \a *this to a constant reference to E. */
      const E& operator() () const { return static_cast<const E&>(*this); }
    };
    


    /* Proxy for a matrix row in operator[] */
    template<class Mx>
    struct MatrixRow
    {
      MatrixRow(Mx& A, const size_type& i) : _mx(A), _i(i) { }
      typename Mx::value_type& operator[](const size_type& j) { return _mx(_i,j); }
      Mx& _mx; const size_type _i;
    };
    
    /*!\ingroup LinearAlgebra
     * \brief A matrix over \a R. 
     *
     * \internal
     * Agreed on number_of_rows() and number_of_columns() for size functions as 
     * this is clearest and best English.
     */
    template<class R>
    class Matrix : public MatrixExpression< Matrix<R> > 
    {
     private:
      size_type _nr;
      size_type _nc;
      array<R> _array;
     private:
      typedef typename Numeric::traits<R>::arithmetic_type F;
     public:
      /*! \brief The type of real number stored by the matrix. */
      typedef R value_type;
       
      /*! \brief Construct a 0 by 0 matrix. */
      explicit Matrix() : _nr(0), _nc(0), _array() { }
      /*! \brief Construct an \a r by \a c matrix, all of whose entries are zero. */
      explicit Matrix(const size_type& r, const size_type& c) : _nr(r), _nc(c), _array(r*c,static_cast<R>(0)) { }
      /*! \brief Construct an \a r by \a c matrix from the array beginning at \a ptr. */
      explicit Matrix(const size_type& nr, const size_type& nc,const R* ptr)
        : _nr(nr), _nc(nc), _array(ptr,ptr+nr*nc) { }
      /*! \brief Construct an \a r by \a c matrix from the array beginning at \a ptr, 
       *  incrementing the input row elements by \a ri and input columns by \a ci. */
      explicit Matrix(const size_type& nr, const size_type& nc, 
             const R* ptr, const size_type& ri, const size_type& ci=1)
        : _nr(nr), _nc(nc), _array(nr*nc) 
      { 
        for(size_type i=0; i!=nr; ++i) { 
          for(size_type j=0; j!=nc; ++j) { 
            (*this)(i,j)=ptr[i*ri+j*ci];
          } 
        } 
      }

      /*! \brief Convert from a matrix expression. */
      template<class E> Matrix(const MatrixExpression<E>& A)
        : _nr(A().number_of_rows()), _nc(A().number_of_columns()), _array(_nc*_nr)
      { 
        const E& mx=A();
        for(size_type i=0; i!=_nr; ++i) { 
          for(size_type j=0; j!=_nc; ++j) { 
            (*this)(i,j)=mx(i,j);
          } 
        } 
      }
        
      /*! \brief Construct from a string literal of the form "[a11,a12,...,a1n; a21,a22,...,a2n;...;am1,am2,...amn]". */
      explicit Matrix(const std::string& s);

      /*! \brief Copy constructor. */
      Matrix(const Matrix<R>& A) : _nr(A._nr), _nc(A._nc), _array(A._array) { }
      /*! \brief Copy assignment operator. */
      Matrix<R>& operator=(const Matrix<R>& A) {
        if(this!=&A) { this->_nr=A._nr; this->_nc=A._nc; this->_array=A._array; } return *this; }
      /*! \brief Assignment operator. */
      template<class E> Matrix<R>& operator=(const MatrixExpression<E>& A) {
        resize(A().number_of_rows(),A().number_of_columns());
        for(size_type i=0; i!=this->number_of_rows(); ++i) {
          for(size_type j=0; j!=this->number_of_columns(); ++j) {
            (*this)(i,j)=A()(i,j); } } 
        return *this;
      }
      
      /*! \brief The equality operator. */
      bool operator==(const Matrix<R>& A) const { 
        return this->_nr==A._nr && this->_array==A._array; }
      /*! \brief The inequality operator. */
      bool operator!=(const Matrix<R>& A) const {
        return !(*this==A); }

      /*! A reference to the array holding the element data. */
      array<R>& data() { return this->_array; }
      /*! A constant reference to the array holding the element data. */
      const array<R>& data() const { return this->_array; }
      
      /*! \brief Resize the matrix. */
      void resize(const size_type& nr, const size_type nc) {
        if(nr*nc!=this->_array.size()) { _array.resize(nr*nc); }
        this->_nr=nr; this->_nc=nc;
      }
      
      /*! \brief The number of rows of the matrix. */
      size_type number_of_rows() const { return this->_nr; }
      /*! \brief The number of columns of the matrix. */
      size_type number_of_columns() const { return this->_nc; }
      /*! \brief A two-element array giving the number of rows and number of columns of the matrix. */
      array<size_type,2u> size() const { return array<size_type,2u>(this->_nr,this->_nc); }

      /*! \brief A constant reference to \a i,\a j th element. */
      const R& operator() (const size_type& i, const size_type& j) const { 
        return this->_array[i*this->_nc+j]; }
      /*! \brief A reference to \a i,\a j th element. */
      R& operator() (const size_type& i, const size_type& j) { 
        return this->_array[i*this->_nc+j]; }
      
      /*! \brief Subscripting operator; A[i][j] returns a reference to the (i,j)th element. */
      MatrixRow< Matrix<R> > operator[](const size_type& i) {
        return MatrixRow< Matrix<R> >(*this,i); }
        
      /*! \brief Constant subscripting operator; A[i][j] returns the (i,j)th element. */
      //MatrixRow< const Matrix<R> > operator[](const size_type& i) const {
      //  return MatrixRow< const Matrix<R> >(*this,i); }
        
      /*! \brief An \a r by \a c matrix, all of whose entries are zero. */
      static Matrix<R> zero(const size_type r, const size_type c) {
        return Matrix<R>(r,c); }
        
      /*! \brief An \a r by \a c matrix, all of whose entries are one. */
      static Matrix<R> one(const size_type r, const size_type c) {
        R one(1); return Matrix<R>(r,c,&one,0,0); }
        
      /*! \brief The \a n by \a n identity matrix. */
      static Matrix<R> identity(const size_type n) {
        Matrix<R> result(n,n); 
        for(size_type i=0; i!=n; ++i) { result(i,i)=R(1); } 
        return result;
      }
      
      /*! \brief The \a i th row. */
      VectorSlice<const R> row(const size_type& i) const {
        return VectorSlice<const R>(this->_nc,this->begin()+i*_nc,1u); }
        
      /*! \brief The \a j th column. */
      VectorSlice<const R> column(const size_type& j) const {
        return VectorSlice<const R>(this->_nr,this->begin()+j,this->_nc); }
        
      /*! \brief Concatenate the columns of two matrices. */
      static Matrix<R> concatenate_columns(const Matrix<R>& A1, const Matrix<R>& A2);
      
      /*! \brief The operator norm with respect to the supremum norm on vectors. 
       * Equal to the supremum over all rows of the sum of absolute values.
       */
      typename Numeric::traits<R>::arithmetic_type norm() const;

      /*! \brief The logarithmic norm. */
      typename Numeric::traits<R>::arithmetic_type log_norm() const;
      
      /*! \brief True if the matrix is singular. */
      bool singular() const;
      /*! \brief The determinant of the matrix. */
      typename Numeric::traits<R>::arithmetic_type determinant() const;
      
      /*! \brief The transposed matrix. */
      Matrix<R> transpose() const;

      /*! \brief The inverse of the matrix. */
      Matrix<typename Numeric::traits<R>::arithmetic_type> inverse() const;
      /*! \brief The solution of the linear equation \f$ Ax=b\f$. */
      Vector<typename Numeric::traits<R>::arithmetic_type> solve(const Vector<R>& b) const;

      /*! \brief A pointer to the first element of the array of values. */
      //R* begin() { return &(*this)(0,0); }
      R* begin() { return this->data().begin(); }
      /*! \brief A constant pointer to the first element of the array of values. */
      //const R* begin() const { return const_cast< Matrix<R>* >(this)->begin(); }
      const R* begin() const { return this->data().begin(); }

      /*! \brief The increment to the array element needed to move one row down. */
      size_type row_increment() const { return this->number_of_columns(); }
      /*! \brief The increment to the array element needed to move one column right. */
      size_type column_increment() const { return 1u; }
      
      /*! \brief Write to an output stream . */
      std::ostream& write(std::ostream& os) const;
      /*! \brief Read from an input stream . */
      std::istream& read(std::istream& is);

#ifdef DOXYGEN
      typedef Numeric::traits<R>::arithmetic_type AT;
      /*! \brief The additive inverse of the matrix \a A. */
      friend Matrix<AT> operator-<>(const Matrix<R>& A);
      /*! \brief The sum of \a A1 and \a A2. */
      friend Matrix<AT> operator+<>(const Matrix<R>& A1, const Matrix<R>& A2);
      /*! \brief The difference of \a A1 and \a A2. */
      friend Matrix<AT> operator-<>(const Matrix<R>& A1, const Matrix<R>& A2);
      /*! \brief The scalar product of \a A by \a s. */
      friend Matrix<AT> operator*<>(const R& s, const Matrix<R>& A);
      /*! \brief The scalar product of \a A by \a s. */
      friend Matrix<AT> operator*<>(const Matrix<R>& A, const R& s);
      /*! \brief The scalar product of \a A by the reciprocal of \a s. */
      friend Matrix<RAT> operator/<>(const Matrix<R>& A, const R& s);
      /*! \brief The product of matrix \a A1 with matrix \a A2. */
      friend Matrix<AT> operator*<>(const Matrix<R>& A1, const Matrix<R>& A2);
      /*! \brief The product of matrix \a A with vector \a v. */
      friend Vector<AT> operator*<>(const Matrix<R>& A, const Vector<R>& v);
      /*! \brief The product of vector \a v with matrix \a A. */
      friend Vector<AT> operator*<>(const Vector<R>& v, const Matrix<R>& A);

      /*! \brief The inverse of \a A. */
      friend Matrix<AT> inverse<>(const Matrix<R>& A);
      /*! \brief Solve the linear system \f$Ax=b\f$. */
      friend Matrix<AT> solve<>(const Matrix<R>& A, const Vector<R>& b);

      /*! \brief Catenate the columns of \a A1 with those of \a A2. */
      friend Matrix<AT> concatenate_columns<>(const Matrix<R>& A1, const Matrix<R>& A2);
      /*! \brief Checks if the matrices have the same dimensions. */
      friend bool have_same_dimensions<>(const Matrix<R>& A1, const Matrix<R>& A2);
#endif 

    };
  


    /*!\ingroup LinearAlgebra
     * \brief A slice through a matrix with equally spaced row and column increments. 
     */
    template<class R>
    class MatrixSlice : public MatrixExpression< MatrixSlice<R> >
    {
      typedef typename Numeric::traits<R>::arithmetic_type F;
     public:
      typedef R value_type;

      MatrixSlice(const size_type& nr, const size_type& nc, R* ptr, const size_type& rinc, const size_type& cinc=1u)
        : _number_of_rows(nr), _number_of_columns(nc), _begin(ptr), _row_increment(rinc), _column_increment(cinc) { }
      MatrixSlice(const size_type& nr, const size_type& nc, R* ptr)
        : _number_of_rows(nr), _number_of_columns(nc), _begin(ptr), _row_increment(nc), _column_increment(1u) { }
      MatrixSlice(const Matrix<R>& m)
        : _number_of_rows(m.number_of_rows()), _number_of_columns(m.number_of_constraints()), _begin(m.begin()),
          _row_increment(m.row_increment()), _column_increment(m.column_increment()) { }

      size_type number_of_rows() const { return this->_number_of_rows; }
      size_type number_of_columns() const { return this->_number_of_columns; }
      array<size_type,2u> size() const { return array<size_type,2>(this->_nr,this->_nc); }
      const R* begin() const { return this->_begin; }
      R* begin() { return this->_begin; }
      size_type row_increment() const { return this->_row_increment; }
      size_type column_increment() const { return this->_column_increment; }
     
      const R& operator() (const size_type& i, const size_type& j) const { 
        return this->_begin[i*this->_row_increment+j*this->_column_increment]; }
      R& operator() (const size_type& i, const size_type& j) { 
        return this->_begin[i*this->_row_increment+j*this->_column_increment]; }
     
      VectorSlice<const R> column(const size_type& j) const {
        return VectorSlice<const R>(this->number_of_rows(),this->begin()+j*this->column_increment(),this->row_increment()); }
        
      F norm() const { return Matrix<R>(*this).norm(); }
      Matrix<F> inverse() const { return Matrix<R>(*this).inverse(); }
      F determinant() const { return Matrix<R>(*this).determinant(); }

      MatrixSlice<R>& operator=(const R& x) {
        for(size_type i=0; i!=this->number_of_rows(); ++i) { 
          for(size_type j=0; j!=this->number_of_columns(); ++j) { 
            (*this)(i,j)=x; } }
        return *this;
      }

      template<class E> MatrixSlice<R>& operator=(const MatrixExpression< E >& m) {
        ARIADNE_CHECK_MATRIX_EQUAL_SIZES(*this,m(),"MatrixSlice& MatrixSlice::operator=(MatrixExpression)");
        const E& e=m(); 
        for(size_type i=0; i!=this->number_of_rows(); ++i) { 
          for(size_type j=0; j!=this->number_of_columns(); ++j) { 
            (*this)(i,j)=e(i,j); } }
        return *this;
      }

      std::ostream& write(std::ostream& os) const { return Matrix<R>(*this).write(os); }
     private:
      size_type _number_of_rows;
      size_type _number_of_columns;
      R* _begin;
      size_type _row_increment;
      size_type _column_increment;
    };
    


    template<class R>
    inline 
    bool
    contains_value(const Matrix< Numeric::Interval<R> >& iA, const Matrix<R>& A) 
    {
      ARIADNE_CHECK_MATRIX_EQUAL_SIZES(iA,A,"bool contains_value(Matrix<Interval>,Matrix<Real>)");
      for(size_type i=0; i!=A.number_of_rows(); ++i) {
        for(size_type j=0; j!=A.number_of_columns(); ++j) {
          if(!Numeric::contains_value(iA(i,j),A(i,j))) {
            return false;
          }
        }
      }
      return true;
    }


    template<class R>
    inline 
    Matrix<R>
    approximate_value(const Matrix< Numeric::Interval<R> >& im) 
    {
      Matrix<R> result(im.number_of_rows(),im.number_of_columns());
      for(size_type i=0; i!=im.number_of_rows(); ++i) {
        for(size_type j=0; j!=im.number_of_columns(); ++j) {
          result(i,j)=approximate_value(im(i,j));
        }
      }
      return result;
    }

    template<class R1,class R2>
    inline 
    Matrix<R1>
    approximate(const Matrix<R2>& im) 
    {
      Matrix<R1> result(im.number_of_rows(),im.number_of_columns());
      for(size_type i=0; i!=im.number_of_rows(); ++i) {
        for(size_type j=0; j!=im.number_of_columns(); ++j) {
          result(i,j)=Numeric::conv_approx<R1>(im(i,j));
        }
      }
      return result;
    }

    template<class R>
    inline
    Vector< Numeric::Interval<R> >
    radius_row_sum(const Matrix< Numeric::Interval<R> >& im) 
    { 
      Vector< Numeric::Interval<R> > result(im.number_of_rows());
      for(dimension_type i=0; i!=im.number_of_rows(); ++i) {
        R radius=0;
        for(dimension_type j=0; j!=im.number_of_columns(); ++j) {
          radius=add_up(radius,im(i,j).length());
        }
        radius = div_up(radius,static_cast<R>(2));
        result[i]=Numeric::Interval<R>(-radius,radius);
      }
      return result;
    }

    
    template<class R>
    inline
    bool
    have_same_dimensions(const Matrix<R>& A1, const Matrix<R>& A2) 
    {
      return A1.number_of_rows()==A2.number_of_rows() 
          && A1.number_of_columns()==A2.number_of_columns();
    }
    
    
    template<class R> inline
    Matrix<R>
    operator-(const Matrix<R>& A) 
    {
      Matrix<R> result(A.number_of_rows(),A.number_of_columns());
      for(size_type i=0; i!=A.number_of_rows(); ++i) {
        for(size_type j=0; j!=A.number_of_columns(); ++j) {
          result(i,j)=-A(i,j);
        }
      }
      return result;
    }
 
    template<class R> inline
    Matrix<R>
    operator-(const MatrixSlice<R>& A) 
    {
      Matrix<R> result(A.number_of_rows(),A.number_of_columns());
      for(size_type i=0; i!=A.number_of_rows(); ++i) {
        for(size_type j=0; j!=A.number_of_columns(); ++j) {
          result(i,j)=-A(i,j);
        }
      }
      return result;
    }
 
    template<class R1, class R2> inline
    Matrix<R1>&
    operator+=(Matrix<R1>& A1, const Matrix<R2>& A2) 
    {
      for(size_type i=0; i!=A1.number_of_rows(); ++i) {
        for(size_type j=0; j!=A1.number_of_columns(); ++j) {
          A1(i,j)+=A2(i,j);
        }
      }
      return A1;
    }
    
    template<class R1, class R2> inline
    Matrix<typename Numeric::traits<R1,R2>::arithmetic_type>
    operator+(const Matrix<R1>& A1, const Matrix<R2>& A2) 
    {
      typedef typename Numeric::traits<R1,R2>::arithmetic_type R3;
      Matrix<R3> result(A1.number_of_rows(),A1.number_of_columns());
      for(size_type i=0; i!=A1.number_of_rows(); ++i) {
        for(size_type j=0; j!=A1.number_of_columns(); ++j) {
          result(i,j)=A1(i,j)+A2(i,j);
        }
      }
      return result;
    }
    
    template<class R1, class R2> inline
    Matrix<typename Numeric::traits<R1,R2>::arithmetic_type>
    operator-(const Matrix<R1>& A1, const Matrix<R2>& A2) 
    {
      typedef typename Numeric::traits<R1,R2>::arithmetic_type R3;
      Matrix<R3> result(A1.number_of_rows(),A1.number_of_columns());
      for(size_type i=0; i!=A1.number_of_rows(); ++i) {
        for(size_type j=0; j!=A1.number_of_columns(); ++j) {
          result(i,j)=A1(i,j)-A2(i,j);
        }
      }
      return result;
    }
    
    template<class R1, class R2> inline
    Matrix<typename Numeric::traits<R1,R2>::arithmetic_type>
    operator*(const R1& s, const Matrix<R2>& A) 
    {
      typedef typename Numeric::traits<R1,R2>::arithmetic_type R3;
      Matrix<R3> result(A.number_of_rows(),A.number_of_columns());
      for(size_type i=0; i!=A.number_of_rows(); ++i) {
        for(size_type j=0; j!=A.number_of_columns(); ++j) {
          result(i,j)=s*A(i,j);
        }
      }
      return result;
    }
    
    template<class R1, class R2> inline
    Matrix<typename Numeric::traits<R1,R2>::arithmetic_type>
    operator*(const Matrix<R1>& A, const R2& s) 
    {
      return s*A;
    }
    
    template<class R1, class R2> inline
    Matrix<typename Numeric::traits<R1,R2>::arithmetic_type>
    operator/(const Matrix<R1>& A, const R2& s) 
    {
      return static_cast<typename Numeric::traits<R2>::arithmetic_type>(static_cast<R2>(1)/s)*A;
    }
    
    template<class R1, class R2> inline
    Vector<typename Numeric::traits<R1,R2>::arithmetic_type>
    operator*(const Matrix<R1>& A, const Vector<R2>& v) 
    {
      typedef typename Numeric::traits<R1,R2>::arithmetic_type R3;
      Vector<R3> result(A.number_of_rows());
      for(size_type i=0; i!=A.number_of_rows(); ++i) {
        for(size_type j=0; j!=A.number_of_columns(); ++j) {
          result(i)+=A(i,j)*v(j);
        }
      }
      return result;
    }
    
    template<class R1, class R2> inline
    Vector<typename Numeric::traits<R1,R2>::arithmetic_type>
    operator*(const Vector<R1>& v, const Matrix<R2>& A) 
    {
      typedef typename Numeric::traits<R1,R2>::arithmetic_type R3;
      Vector<R3> result(A.number_of_columns());
      for(size_type i=0; i!=A.number_of_rows(); ++i) {
        for(size_type j=0; j!=A.number_of_columns(); ++j) {
          result(j)+=v(i)*A(i,j);
        }
      }
      return result;
    }
    
    template<class R1, class R2> inline
    Matrix<typename Numeric::traits<R1,R2>::arithmetic_type>
    operator*(const Matrix<R1>& A1, const Matrix<R2>& A2) 
    {
      typedef typename Numeric::traits<R1,R2>::arithmetic_type R3;
      Matrix<R3> result(A1.number_of_rows(),A2.number_of_columns());
      for(size_type i=0; i!=A1.number_of_rows(); ++i) {
        for(size_type j=0; j!=A2.number_of_columns(); ++j) {
          for(size_type k=0; k!=A1.number_of_columns(); ++k) {
            result(i,j)+=A1(i,k)*A2(k,j);
          }
        }
      }
      return result;
    }
    
    template<class R1, class R2> inline
    Matrix<typename Numeric::traits<R1,R2>::arithmetic_type>
    operator*(const Matrix<R1>& A1, const MatrixSlice<R2>& A2) 
    {
      return A1*Matrix<R2>(A2);
    }
    
    template<class R1, class R2> inline
    Matrix<typename Numeric::traits<R1,R2>::arithmetic_type>
    operator*(const MatrixSlice<R1>& A1, const Matrix<R2>& A2) 
    {
      return Matrix<R1>(A1)*A2;
    }
    
    template<class R1, class R2> inline
    Matrix<typename Numeric::traits<R1,R2>::arithmetic_type>
    operator*(const MatrixSlice<R1>& A1, const MatrixSlice<R2>& A2) 
    {
      return Matrix<R1>(A1)*Matrix<R2>(A2);
    }
    
    
 
    template<class R1, class R2> inline
    Matrix<typename Numeric::traits<R1,R2>::arithmetic_type>
    operator*(const R1& s, const MatrixSlice<R2>& A) 
    {
      typedef typename Numeric::traits<R1,R2>::arithmetic_type R3;
      Matrix<R3> result(A.number_of_rows(),A.number_of_columns());
      for(size_type i=0; i!=A.number_of_rows(); ++i) {
        for(size_type j=0; j!=A.number_of_columns(); ++j) {
          result(i,j)=s*A(i,j);
        }
      }
      return result;
    }
    
    template<class R1, class R2> inline
    Matrix<typename Numeric::traits<R1,R2>::arithmetic_type>
    operator*(const MatrixSlice<R1>& A, const R2& s) 
    {
      return s*A;
    }
    
    template<class R1, class R2> inline
    Vector<typename Numeric::traits<R1,R2>::arithmetic_type>
    operator*(const MatrixSlice<R1>& A, const Vector<R2>& v) 
    {
      typedef typename Numeric::traits<R1,R2>::arithmetic_type R3;
      Vector<R3> result(A.number_of_rows());
      for(size_type i=0; i!=A.number_of_rows(); ++i) {
        for(size_type j=0; j!=A.number_of_columns(); ++j) {
          result(i)+=A(i,j)*v(j);
        }
      }
      return result;
    }
    

    template<class R1, class R2> inline
    Matrix<typename Numeric::traits<R1,R2>::arithmetic_type>
    outer_product(const Vector<R1>& v1, const Vector<R2>& v2) 
    {
      typedef typename Numeric::traits<R1,R2>::arithmetic_type R3;
      Matrix<R3> result(v1.size(),v2.size());
      for(size_type i=0; i!=v1.size(); ++i) {
        for(size_type j=0; j!=v2.size(); ++j) {
          result(i,j)=v1(i)*v2(j);
        }
      }
      return result;
    }
    




    template<class R>
    inline
    Vector<typename Numeric::traits<R>::arithmetic_type>
    row_norms(const Matrix<R>& A) 
    {
      Vector<typename Numeric::traits<R>::arithmetic_type> result(A.number_of_rows());
      for(size_type i=0; i!=A.number_of_rows(); ++i) {
        for(size_type j=0; j!=A.number_of_columns(); ++j) {
          result(i)+=abs(A(i,j));
        }
      }
      return result;
    }
    
    template<class R>
    inline
    typename Numeric::traits<R>::arithmetic_type
    norm(const Matrix<R>& A) 
    {
      return A.norm();
    }
    
    template<class R>
    inline
    typename Numeric::traits<R>::arithmetic_type
    norm(const MatrixSlice<R>& A) 
    {
      return A.norm();
    }
    
    template<class R>
    inline
    typename Numeric::traits<R>::arithmetic_type
    log_norm(const Matrix<R>& A)
    {
      return A.log_norm();
    }

    template<class R>
    inline
    Matrix<typename Numeric::traits<R>::arithmetic_type> 
    inverse(const Matrix<R>& A)
    {
      return A.inverse();
    }
    
    template<class R>
    inline
    Vector<typename Numeric::traits<R>::arithmetic_type> 
    solve(const Matrix<R>& A, const Vector<R>& v)
    {
      return A.solve(v);
    }
    
    template<class R>
    inline
    Matrix<R>
    concatenate_columns(const Matrix<R>& A1, const Matrix<R>& A2) {
      return Matrix<R>::concatenate_columns(A1,A2);
    }
    
    template<class R>
    inline
    Matrix<R>
    concatenate_columns(const MatrixSlice<R>& A1, const MatrixSlice<R>& A2) {
      return concatenate_columns(Matrix<R>(A1),Matrix<R>(A2));
    }
    
    
    /*! \brief A matrix \f$A\f$ such that for all zonotopes \f$Z\f$, \f$AZ\subset \overline{\underline{A}}\f$. */
    template<class R>
    Matrix<R>
    over_approximation(const Matrix< Numeric::Interval<R> >& A)
    {
      if(A.number_of_rows()!=A.number_of_columns()) {
        ARIADNE_THROW(NotImplemented,"Matrix<Real> over_approximation(Matrix<Interval>)","A="<<A<<" (only implemented for square matrices)"); 
      }
      dimension_type n=A.number_of_rows();
      
      Matrix<R> Amid(n,n);
      for(size_type i=0; i!=n; ++i) {
        for(size_type j=0; j!=n; ++j) {
          Amid(i,j)=(A(i,j).upper()+A(i,j).lower())/2;
        }
      }
      Matrix< Numeric::Interval<R> > I=LinearAlgebra::inverse(Matrix< Numeric::Interval<R> >(Amid))*A;
      
      R excess=LinearAlgebra::norm(I).upper();
      
      // FIXME: Outer bound on multiplication
      return excess*Amid;
    }


    template<class R>
    std::ostream&
    operator<<(std::ostream& os, const MatrixSlice<R>& A) {
      return A.write(os);
    }

    template<class R>
    std::ostream&
    operator<<(std::ostream& os, const Matrix<R>& A) {
      return A.write(os);
    }

    template<class R>
    std::istream&
    operator>>(std::istream& is, Matrix<R>& A) {
      return A.read(is); 
    }
   

  }
}


#endif /* ARIADNE_MATRIX_H */
