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

#ifndef _ARIADNE_MATRIX_H
#define _ARIADNE_MATRIX_H

#include <iosfwd>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include "../declarations.h"
#include "../numeric/integer.h"
#include "../numeric/interval.h"

namespace Ariadne {
  namespace LinearAlgebra {

    using boost::numeric::ublas::identity_matrix;
    using boost::numeric::ublas::herm;
    using boost::numeric::ublas::matrix_row;
    using boost::numeric::ublas::matrix_column;
  
    /*! \ingroup LinearAlgebra
     *  \brief A matrix over \a R. 
     */
    template<class R>
    class Matrix : public boost::numeric::ublas::matrix<R> 
    {
      typedef boost::numeric::ublas::matrix<R> _Base;
      typedef boost::numeric::ublas::matrix<R> _boost_matrix;
      typedef typename Numeric::traits<R>::arithmetic_type F;
     public:
      /*! \brief Construct a 0 by 0 matrix. */
      Matrix() : _Base() { }
      /*! \brief Construct an \a r by \a c matrix, all of whose entries are zero. */
      Matrix(const size_type& r, const size_type& c) : _Base(r,c) { }
      /*! \brief Construct an \a r by \a c matrix from the array beginning at \a ptr, 
       *  incrementing the input row elements by \a ri and input columns by \a ci. */
      Matrix(const size_type& nr, const size_type& nc, 
             const R* ptr, const size_type& ri, const size_type& ci=1) : _Base(nr,nc) 
      { 
        for(size_type i=0; i!=nr; ++i) { 
          for(size_type j=0; j!=nc; ++j) { 
            _Base::operator()(i,j)=ptr[i*ri+j*ci];
          } 
        } 
      }

      template<class E> Matrix(const boost::numeric::ublas::matrix_expression<E>& A) : _Base(A) { }
            
      /*! \brief Construct from a string literal of the form "[a11,a12,...,a1n; a21,a22,...,a2n;...;am1,am2,...amn]". */
      explicit Matrix(const std::string& s);

#ifdef DOXYGEN
      /*! \brief Copy constructor. */
      Matrix(const Matrix<R>& A);
      /*! \brief Copy assignment operator. */
      Matrix<R>& operator=(const Matrix<R>& A);
#endif
      
      /*! \brief The equality operator. */
      bool operator==(const Matrix<R>& A) const;
      /*! \brief The inequality operator. */
      bool operator!=(const Matrix<R>& A) const;

    
     
      /*! \brief The number of elements in the \a d th dimension. */
      size_type size(const size_type& d) const {
        if(d==0) { return number_of_rows(); } 
        else if(d==1) { return number_of_columns(); }
        else { assert(false); } 
      }
      /*! \brief The number of rows of the matrix. */
      size_type number_of_rows() const { return _Base::size1(); }
      /*! \brief The number of columns of the matrix. */
      size_type number_of_columns() const { return _Base::size2(); }
      
      /*! \brief A constant reference to \a i,\a j th element. */
      const R& operator() (const size_type& i, const size_type& j) const { 
        return this->_Base::operator()(i,j); }
      /*! \brief A reference to \a i,\a j th element. */
      R& operator() (const size_type& i, const size_type& j) { 
        return this->_Base::operator()(i,j); }
      
      // /*! \brief A slice through the \a i th row. */
      //VectorSlice<R> operator[] (const size_type& i) {
      //  return VectorSlice<R>(this->number_of_columns(),this->begin()+i*this->column_increment(),this->row_increment()); }
        
      /*! \brief An \a r by \a c matrix, all of whose entries are zero. */
      static Matrix<R> zero(const size_type r, const size_type c) {
        return Matrix<R>(r,c); }
      /*! \brief The \a n by \a n identity matrix. */
      static Matrix<R> identity(const size_type n) {
        Matrix<R> result(n,n); 
        for(size_type i=0; i!=n; ++i) { result(i,i)=R(1); } 
        return result;
      }
      
      /*! \brief The \a i th row. */
      matrix_row< Matrix<R> > row(const size_type& i) {
        return boost::numeric::ublas::row(*this,i); }
        
      matrix_column< Matrix<R> > column(const size_type& j) {
        return boost::numeric::ublas::column(*this,j); }
        
      /*! \brief Concatenate the columns of two matrices. */
      static Matrix<R> concatenate_columns(const Matrix<R>& A1, const Matrix<R>& A2);
      
      /*! \brief The operator norm with respect to the supremum norm on vectors. 
       * Equal to the supremum over all rows of the sum of absolute values.
       */
      F norm() const;

      /*! \brief The logarithmic norm. */
      F log_norm() const;
      
      /*! \brief True if the matrix is singular. */
      bool singular() const;
      /*! \brief The determinant of the matrix. */
      F determinant() const;
      
      /*! \brief The transposed matrix. */
      Matrix<R> transpose() const;

      /*! \brief The inverse of the matrix. */
      Matrix<F> inverse() const;
      /*! \brief The solution of the linear equation \f$ Ax=b\f$. */
      Vector<F> solve(const Vector<R>& b) const;

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
      /*! \brief The additive inverse of the matrix \a A. */
      friend Matrix<R> operator-<>(const Matrix<R>& A);
      /*! \brief The sum of \a A1 and \a A2. */
      friend Matrix<R> operator+<>(const Matrix<R>& A1, const Matrix<R>& A2);
      /*! \brief The difference of \a A1 and \a A2. */
      friend Matrix<R> operator-<>(const Matrix<R>& A1, const Matrix<R>& A2);
      /*! \brief The scalar product of \a A by \a s. */
      friend Matrix<R> operator*<>(const R& s, const Matrix<R>& A);
      /*! \brief The scalar product of \a A by \a s. */
      friend Matrix<R> operator*<>(const Matrix<R>& A, const R& s);
      /*! \brief The scalar product of \a A by the reciprocal of \a s. */
      friend Matrix<R> operator/<>(const Matrix<R>& A, const R& s);
      /*! \brief The product of matrix \a A1 with matrix \a A2. */
      friend Matrix<R> operator*<>(const Matrix<R>& A1, const Matrix<R>& A2);
      /*! \brief The product of matrix \a A with vector \a v. */
      friend Vector<R> operator*<>(const Matrix<R>& A, const Vector<R>& v);

      /*! \brief The inverse of \a A. */
      friend Matrix<R> inverse<>(const Matrix<R>& A);
      /*! \brief Solve the linear system \f$Ax=b\f$. */
      friend Matrix<R> solve<>(const Matrix<R>& A, const Vector<R>& b);

      /*! \brief Catenate the columns of \a A1 with those of \a A2. */
      friend Matrix<R> concatenate_columns<>(const Matrix<R>& A1, const Matrix<R>& A2);
      /*! \brief Checks if the matrices have the same dimensions. */
      friend bool have_same_dimensions<>(const Matrix<R>& A1, const Matrix<R>& A2);
#endif 

    };
  
    /*!\ingroup LinearAlgebra
     * \brief A slice through a matrix with equally spaced row and column increments. 
     */
    template<class R>
    class MatrixSlice : public boost::numeric::ublas::matrix_expression< MatrixSlice<R> >
    {
     public:
      MatrixSlice(const size_type& nr, const size_type& nc, R* ptr, const size_type& rinc, const size_type& cinc=1u)
        : _number_of_rows(nr), _number_of_columns(nc), _begin(ptr), _row_increment(rinc), _column_increment(cinc) { }
      MatrixSlice(const size_type& nr, const size_type& nc, R* ptr)
        : _number_of_rows(nr), _number_of_columns(nc), _begin(ptr), _row_increment(nc), _column_increment(1u) { }
      MatrixSlice(const Matrix<R>& m)
        : _number_of_rows(m.number_of_rows()), _number_of_columns(m.number_of_constraints()), _begin(m.begin()),
          _row_increment(m.row_increment()), _column_increment(m.column_increment()) { }

      size_type size1() const { return this->_number_of_rows; }
      size_type size2() const { return this->_number_of_columns; }
      
      size_type number_of_rows() const { return this->_number_of_rows; }
      size_type number_of_columns() const { return this->_number_of_columns; }
      const R* begin() const { return this->_begin; }
      R* begin() { return this->_begin; }
      size_type row_increment() const { return this->_row_increment; }
      size_type column_increment() const { return this->_column_increment; }
     
      const R& operator() (const size_type& i, const size_type& j) const { 
        return this->_begin[i*this->_row_increment+j*this->_column_increment]; }
      R& operator() (const size_type& i, const size_type& j) { 
        return this->_begin[i*this->_row_increment+j*this->_column_increment]; }
     
      template<class E> MatrixSlice<R>& operator=(const boost::numeric::ublas::matrix_expression< E > m) {
        const E& e=m(); assert(this->number_of_rows()==e.size1() && this->number_of_columns()==e.size2());
        for(size_type i=0; i!=this->number_of_rows(); ++i) { 
          for(size_type j=0; j!=this->number_of_columns(); ++j) { 
            (*this)(i,j)=e(i,j); } }
        return *this;
      }
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
    contains_value(const Matrix< Interval<R> >& iA, const Matrix<R>& A) 
    {
      assert(A.number_of_rows()==iA.number_of_rows() 
        && A.number_of_columns()==iA.number_of_columns());
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
    approximate_value(const Matrix< Interval<R> >& im) 
    {
      Matrix<R> result(im.number_of_rows(),im.number_of_columns());
      for(size_type i=0; i!=im.number_of_rows(); ++i) {
        for(size_type j=0; j!=im.number_of_columns(); ++j) {
          result(i,j)=approximate_value(im(i,j));
        }
      }
      return result;
    }

    template<class R>
    inline
    Vector< Interval<R> >
    radius_row_sum(const Matrix< Interval<R> >& im) 
    { 
      Vector< Interval<R> > result(im.number_of_rows());
      for(dimension_type i=0; i!=im.number_of_rows(); ++i) {
        R radius=0;
        for(dimension_type j=0; j!=im.number_of_columns(); ++j) {
          radius=add_up(radius,im(i,j).length());
        }
        radius = div_up(radius,static_cast<R>(2));
        result[i]=Interval<R>(-radius,radius);
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
    
    template<class R>
    inline
    Vector<R>
    operator*(const Matrix<R>& A, const Vector<R>& B) {
      return boost::numeric::ublas::prod(A,B);
    }
           
    template<class R>
    inline
    Vector< Interval<R> >
    operator*(const Matrix<R>& A, const Vector< Interval<R> >& B) {
      return boost::numeric::ublas::prod(Matrix< Interval<R> >(A),B);
    }
           
    template<class R>
    inline
    Vector< Interval<R> >
    operator*(const Matrix< Interval<R> >& A, const Vector<R>& B) {
      return boost::numeric::ublas::prod(A,Vector< Interval<R> >(B));
    }
           
    
    template<class R>
    inline
    Matrix<R>
    operator*(const Matrix<R>& A, const Matrix<R>& B) {
      return boost::numeric::ublas::prod(A,B);
    }
    
    template<class R>
    inline
    Matrix< Interval<R> >
    operator*(const Matrix<R>& A, const Matrix< Interval<R> >& B) {
      return boost::numeric::ublas::prod(Matrix< Interval<R> >(A),B);
    }
    
    template<class R>
    inline
    Matrix< Interval<R> >
    operator*(const Matrix< Interval<R> >& A, const Matrix<R>& B) {
      return boost::numeric::ublas::prod(A,Matrix< Interval<R> >(B));
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
    
    
    /*! \brief A matrix \f$A\f$ such that for all zonotopes \f$Z\f$, \f$AZ\subset \overline{\underline{A}}\f$. */
    template<class R>
    Matrix<R>
    over_approximation(const Matrix< Interval<R> >& A)
    {
      assert(A.number_of_rows()==A.number_of_columns());
      dimension_type n=A.number_of_rows();
      
      Matrix<R> Amid(n,n);
      for(size_type i=0; i!=n; ++i) {
        for(size_type j=0; j!=n; ++j) {
          Amid(i,j)=(A(i,j).upper()+A(i,j).lower())/2;
        }
      }
      Matrix< Interval<R> > I=LinearAlgebra::inverse(Matrix< Interval<R> >(Amid))*A;
      
      R excess=LinearAlgebra::norm(I).upper();
      
      // FIXME: Outer bound on multiplication
      return excess*Amid;
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

    template<>
    class Matrix<int> : public boost::numeric::ublas::matrix<int> 
    {
     public:
      Matrix()
        : boost::numeric::ublas::matrix<int>() { }
      Matrix(const size_type& nr, const size_type& nc)
        : boost::numeric::ublas::matrix<int>(nr,nc) { }
    };
             
    
/*        
        
    
    
    template<class R>
    bool independent_rows(Matrix<R> A);

    
    template<class R>
    bool 
    equivalent_columns(const Matrix<R> &A, 
                       const size_type &A_col, 
                       const Matrix<R> &B, 
                       const size_type &B_col);
    
  
    template<class R>
    size_type 
    find_first_not_null_in_col(const Matrix<R> &A, 
                               const size_type &col);
    
    template<class R>
    Matrix<R> 
    remove_null_columns_but_one(const Matrix<R> &A);
    
    template<class R>
    void 
    remove_null_columns(const Matrix<R>& A, 
                        array<size_type>& row, 
                        array<size_type>& col);

    template<class R>
    Matrix<R> 
    compute_space(const Matrix<R>& SA, 
                  array<size_type>& row,
                  const array<size_type>& col);
*/

  }
}


#endif /* _ARIADNE_MATRIX_H */
