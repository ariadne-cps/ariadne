/***************************************************************************
 *            matrix.inline.h
 *
 *  Copyright  2004-7  Alberto Casagrande, Pieter Collins
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
 

namespace Ariadne {

template<class R> inline
LinearAlgebra::Matrix<R>::Matrix()
  : _nr(0), _nc(0), _array() 
{ 
}


template<class R> inline
LinearAlgebra::Matrix<R>::Matrix(const size_type& r, const size_type& c)
  : _nr(r), _nc(c), _array(r*c,static_cast<R>(0)) 
{ 
}


template<class R> template<class RR> inline
LinearAlgebra::Matrix<R>::Matrix(const size_type& nr, const size_type& nc, const RR* ptr)
  
  : _nr(nr), _nc(nc), _array(ptr,ptr+nr*nc) 
{ 
}

template<class R> template<int NN, class RR> inline
LinearAlgebra::Matrix<R>::Matrix(const size_type& nr, const size_type& nc, const RR ary[][NN])
  
  : _nr(nr), _nc(nc), _array(nr*nc) 
{ 
  for(size_type i=0; i!=nr; ++i) { 
    for(size_type j=0; j!=nc; ++j) { 
      (*this)(i,j)=static_cast<RR>(ary[i][j]);
    } 
  } 
}


template<class R> template<class RR> inline
LinearAlgebra::Matrix<R>::Matrix(const size_type& nr, const size_type& nc, 
                                 const RR* ptr, const size_type& ri, const size_type& ci)
  : _nr(nr), _nc(nc), _array(nr*nc) 
{ 
  for(size_type i=0; i!=nr; ++i) { 
    for(size_type j=0; j!=nc; ++j) { 
      (*this)(i,j)=ptr[i*ri+j*ci];
    } 
  } 
}


template<class R> template<class E> inline
LinearAlgebra::Matrix<R>::Matrix(const MatrixExpression<E>& A)
  : _nr(A().number_of_rows()), _nc(A().number_of_columns()), _array(_nc*_nr)
{ 
  const E& mx=A();
  for(size_type i=0; i!=_nr; ++i) { 
    for(size_type j=0; j!=_nc; ++j) { 
      (*this)(i,j)=mx(i,j);
    } 
  } 
}


template<class R> inline
LinearAlgebra::Matrix<R>::Matrix(const Matrix<R>& A)
  : _nr(A._nr), _nc(A._nc), _array(A._array) 
{ 
}


template<class R> inline
LinearAlgebra::Matrix<R>& 
LinearAlgebra::Matrix<R>::operator=(const Matrix<R>& A) 
{
  if(this!=&A) { 
    this->_nr=A._nr; 
    this->_nc=A._nc; 
    this->_array=A._array;
  }
  return *this; 
}


template<class R> template<class E> inline
LinearAlgebra::Matrix<R>& 
LinearAlgebra::Matrix<R>::operator=(const MatrixExpression<E>& A) 
{
  resize(A().number_of_rows(),A().number_of_columns());
  for(size_type i=0; i!=this->number_of_rows(); ++i) {
    for(size_type j=0; j!=this->number_of_columns(); ++j) {
      (*this)(i,j)=A()(i,j); 
    }
  } 
  return *this;
}


template<class R> inline
bool 
LinearAlgebra::Matrix<R>::operator==(const Matrix<R>& A) const 
{ 
  return this->_nr==A._nr && this->_array==A._array;
}


template<class R> inline
bool 
LinearAlgebra::Matrix<R>::operator!=(const Matrix<R>& A) const 
{
  return !(*this==A);
}


template<class R> inline
array<R>& 
LinearAlgebra::Matrix<R>::data() 
{ 
  return this->_array; 
}


template<class R> inline
const array<R>& 
LinearAlgebra::Matrix<R>::data() const 
{ 
  return this->_array;
}


template<class R> inline
void 
LinearAlgebra::Matrix<R>::resize(const size_type& nr, const size_type nc) 
{
  if(nr*nc!=this->_array.size()) {
    _array.resize(nr*nc);
  }
  this->_nr=nr;
  this->_nc=nc;
}


template<class R> inline
size_type 
LinearAlgebra::Matrix<R>::number_of_rows() const 
{ 
  return this->_nr; 
}


template<class R> inline
size_type 
LinearAlgebra::Matrix<R>::number_of_columns() const 
{ 
  return this->_nc; 
}


template<class R> inline
array<size_type,2u> 
LinearAlgebra::Matrix<R>::size() const 
{ 
  return array<size_type,2u>(this->_nr,this->_nc); 
}


template<class R> inline
const R& 
LinearAlgebra::Matrix<R>::operator() (const size_type& i, const size_type& j) const 
{ 
  return this->_array[i*this->_nc+j]; 
}


template<class R> inline
R& 
LinearAlgebra::Matrix<R>::operator() (const size_type& i, const size_type& j) 
{ 
  return this->_array[i*this->_nc+j]; 
}


template<class R> inline
LinearAlgebra::MatrixRow< LinearAlgebra::Matrix<R> > 
LinearAlgebra::Matrix<R>::operator[](const size_type& i) 
{
  return MatrixRow< Matrix<R> >(*this,i); 
}        


template<class R> inline
LinearAlgebra::Matrix<R>
LinearAlgebra::Matrix<R>:: zero(const size_type r, const size_type c) 
{
  return Matrix<R>(r,c); 
}


template<class R> inline
LinearAlgebra::Matrix<R> 
LinearAlgebra::Matrix<R>::one(const size_type r, const size_type c) 
{
  R one(1); return Matrix<R>(r,c,&one,0,0); 
}


template<class R> inline
LinearAlgebra::Matrix<R> 
LinearAlgebra::Matrix<R>::identity(const size_type n) 
{
  Matrix<R> result(n,n); 
  for(size_type i=0; i!=n; ++i) { 
    result(i,i)=R(1); 
  } 
  return result;
}


template<class R> inline
LinearAlgebra::VectorSlice<const R>
LinearAlgebra::Matrix<R>::row(const size_type& i) const 
{
  return VectorSlice<const R>(this->_nc,this->begin()+i*_nc,1u); 
}


template<class R> inline
LinearAlgebra::VectorSlice<const R> 
LinearAlgebra::Matrix<R>::column(const size_type& j) const 
{
  return VectorSlice<const R>(this->_nr,this->begin()+j,this->_nc); 
}


template<class R> inline
R* 
LinearAlgebra::Matrix<R>::begin() 
{
  return this->data().begin(); 
}


template<class R> inline
const R* 
LinearAlgebra::Matrix<R>::begin() const 
{
  return this->data().begin(); 
}


template<class R> inline
size_type 
LinearAlgebra::Matrix<R>::row_increment() const 
{
  return this->number_of_columns(); 
}


template<class R> inline
size_type 
LinearAlgebra::Matrix<R>::column_increment() const 
{ 
  return 1u; 
}






template<class R> inline
LinearAlgebra::MatrixSlice<R>::MatrixSlice(const size_type& nr, const size_type& nc, R* ptr, const size_type& rinc, const size_type& cinc)
  
  : _number_of_rows(nr), _number_of_columns(nc), _begin(ptr), _row_increment(rinc), _column_increment(cinc) 
{ 
}


template<class R> inline
LinearAlgebra::MatrixSlice<R>::MatrixSlice(const size_type& nr, const size_type& nc, R* ptr)
  
  : _number_of_rows(nr), _number_of_columns(nc), _begin(ptr), _row_increment(nc), _column_increment(1u) 
{ 
}


template<class R> inline
LinearAlgebra::MatrixSlice<R>::MatrixSlice(const Matrix<R>& m)
  
  : _number_of_rows(m.number_of_rows()), _number_of_columns(m.number_of_constraints()), _begin(m.begin()),
    _row_increment(m.row_increment()), _column_increment(m.column_increment()) 
{ 
}


template<class R> inline
size_type
LinearAlgebra::MatrixSlice<R>::number_of_rows() const 
{
  return this->_number_of_rows;
}

template<class R> inline
size_type
LinearAlgebra::MatrixSlice<R>::number_of_columns() const { 
  return this->_number_of_columns;
}


template<class R> inline
array<size_type,2u>
LinearAlgebra::MatrixSlice<R>::size() const 
{
  return array<size_type,2>(this->_nr,this->_nc); }
template<class R> inline
const R*
LinearAlgebra::MatrixSlice<R>::begin() const 
{
  return this->_begin;
}


template<class R> inline
R*
LinearAlgebra::MatrixSlice<R>::begin() 
{
  return this->_begin;
}


template<class R> inline
size_type
LinearAlgebra::MatrixSlice<R>::row_increment() const 
{
  return this->_row_increment;
}


template<class R> inline
size_type
LinearAlgebra::MatrixSlice<R>::column_increment() const 
{
  return this->_column_increment;
}

template<class R> inline
const R&
LinearAlgebra::MatrixSlice<R>::operator() (const size_type& i, const size_type& j) const 
{
  return this->_begin[i*this->_row_increment+j*this->_column_increment];
}


template<class R> inline
R&
LinearAlgebra::MatrixSlice<R>::operator() (const size_type& i, const size_type& j) 
{
  return this->_begin[i*this->_row_increment+j*this->_column_increment]; 
}


template<class R> inline
LinearAlgebra::VectorSlice<const R>
LinearAlgebra::MatrixSlice<R>::column(const size_type& j) const 
{
  return VectorSlice<const R>(this->number_of_rows(),this->begin()+j*this->column_increment(),this->row_increment());
}


template<class R> inline
LinearAlgebra::MatrixSlice<R>&
LinearAlgebra::MatrixSlice<R>::operator=(const R& x) 
{
  for(size_type i=0; i!=this->number_of_rows(); ++i) { 
    for(size_type j=0; j!=this->number_of_columns(); ++j) { 
      (*this)(i,j)=x;
    } 
  }
  return *this;
}


template<class R> template<class E> inline
LinearAlgebra::MatrixSlice<R>&
LinearAlgebra::MatrixSlice<R>::operator=(const MatrixExpression< E >& m) 
{
  ARIADNE_CHECK_MATRIX_EQUAL_SIZES(*this,m(),"MatrixSlice& MatrixSlice::operator=(MatrixExpression)");
  const E& e=m(); 
  for(size_type i=0; i!=this->number_of_rows(); ++i) { 
    for(size_type j=0; j!=this->number_of_columns(); ++j) { 
      (*this)(i,j)=e(i,j);
    }
  }
  return *this;
}


template<class R> inline
std::ostream&
LinearAlgebra::MatrixSlice<R>::write(std::ostream& os) const 
{
  return Matrix<R>(*this).write(os);
}




template<class R1,class R2> inline 
LinearAlgebra::Matrix<R1>
LinearAlgebra::approximation(const Matrix<R2>& im) 
{
  Matrix<R1> result(im.number_of_rows(),im.number_of_columns());
  for(size_type i=0; i!=im.number_of_rows(); ++i) {
    for(size_type j=0; j!=im.number_of_columns(); ++j) {
      result(i,j)=Numeric::conv_approx<R1>(im(i,j));
    }
  }
  return result;
}


template<class R> inline 
LinearAlgebra::Matrix<R>
LinearAlgebra::midpoint(const Matrix< Numeric::Interval<R> >& im) 
{
  Matrix<R> result(im.number_of_rows(),im.number_of_columns());
  for(size_type i=0; i!=im.number_of_rows(); ++i) {
    for(size_type j=0; j!=im.number_of_columns(); ++j) {
      result(i,j)=midpoint(im(i,j));
    }
  }
  return result;
}


template<class R> inline 
bool
LinearAlgebra::encloses(const Matrix< Numeric::Interval<R> >& iA, const Matrix<R>& A) 
{
  ARIADNE_CHECK_MATRIX_EQUAL_SIZES(iA,A,"bool encloses(Matrix<Interval>,Matrix<Real>)");
  for(size_type i=0; i!=A.number_of_rows(); ++i) {
    for(size_type j=0; j!=A.number_of_columns(); ++j) {
      if(!Numeric::encloses(iA(i,j),A(i,j))) {
        return false;
      }
    }
  }
  return true;
}


template<class R> inline 
bool
LinearAlgebra::refines(const Matrix< Numeric::Interval<R> >& iA1, const Matrix< Numeric::Interval<R> >& iA2) 
{
  ARIADNE_CHECK_MATRIX_EQUAL_SIZES(iA1,iA2,"bool refines(Matrix<Interval>,Matrix<Interval>)");
  for(size_type i=0; i!=iA1.number_of_rows(); ++i) {
    for(size_type j=0; j!=iA1.number_of_columns(); ++j) {
      if(!Numeric::refines(iA1(i,j),iA2(i,j))) {
        return false;
      }
    }
  }
  return true;
}



template<class R> inline
LinearAlgebra::Vector< Numeric::Interval<R> >
LinearAlgebra::radius_row_sum(const Matrix< Numeric::Interval<R> >& im) 
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


template<class R> inline
bool
LinearAlgebra::have_same_dimensions(const Matrix<R>& A1, const Matrix<R>& A2) 
{
  return A1.number_of_rows()==A2.number_of_rows() 
    && A1.number_of_columns()==A2.number_of_columns();
}


template<class R> inline 
LinearAlgebra::Matrix<R>
LinearAlgebra::operator-(const Matrix<R>& A) 
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
LinearAlgebra::Matrix<R>
LinearAlgebra::operator-(const MatrixSlice<R>& A) 
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
LinearAlgebra::Matrix<R1>&
LinearAlgebra::operator+=(Matrix<R1>& A1, const Matrix<R2>& A2) 
{
  for(size_type i=0; i!=A1.number_of_rows(); ++i) {
    for(size_type j=0; j!=A1.number_of_columns(); ++j) {
      A1(i,j)+=A2(i,j);
    }
  }
  return A1;
}


template<class R1, class R2> inline 
LinearAlgebra::Matrix<R1>&
LinearAlgebra::operator-=(Matrix<R1>& A1, const Matrix<R2>& A2) 
{
  for(size_type i=0; i!=A1.number_of_rows(); ++i) {
    for(size_type j=0; j!=A1.number_of_columns(); ++j) {
      A1(i,j)-=A2(i,j);
    }
  }
  return A1;
}


template<class R1, class R2>  inline
LinearAlgebra::Matrix<typename Numeric::traits<R1,R2>::arithmetic_type>
LinearAlgebra::operator+(const Matrix<R1>& A1, const Matrix<R2>& A2) 
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


template<class R1, class R2>  inline
LinearAlgebra::Matrix<typename Numeric::traits<R1,R2>::arithmetic_type>
LinearAlgebra::operator-(const Matrix<R1>& A1, const Matrix<R2>& A2) 
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


template<class R1, class R2>  inline
LinearAlgebra::Matrix<typename Numeric::traits<R1,R2>::arithmetic_type>
LinearAlgebra::operator*(const R1& s, const Matrix<R2>& A) 
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

template<class R1, class R2>  inline
LinearAlgebra::Matrix<typename Numeric::traits<R1,R2>::arithmetic_type>
LinearAlgebra::operator*(const Matrix<R1>& A, const R2& s) 
{
  return s*A;
}

template<class R1, class R2>  inline
LinearAlgebra::Matrix<typename Numeric::traits<R1,R2>::arithmetic_type>
LinearAlgebra::operator/(const Matrix<R1>& A, const R2& s) 
{
  return static_cast<typename Numeric::traits<R2>::arithmetic_type>(static_cast<R2>(1)/s)*A;
}

template<class R1, class R2>  inline
LinearAlgebra::Vector<typename Numeric::traits<R1,R2>::arithmetic_type>
LinearAlgebra::operator*(const Matrix<R1>& A, const Vector<R2>& v) 
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

template<class R1, class R2>  inline
LinearAlgebra::Vector<typename Numeric::traits<R1,R2>::arithmetic_type>
LinearAlgebra::operator*(const Vector<R1>& v, const Matrix<R2>& A) 
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

template<class R1, class R2>  inline
LinearAlgebra::Matrix<typename Numeric::traits<R1,R2>::arithmetic_type>
LinearAlgebra::operator*(const Matrix<R1>& A1, const Matrix<R2>& A2) 
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

template<class R1, class R2>  inline
LinearAlgebra::Matrix<typename Numeric::traits<R1,R2>::arithmetic_type>
LinearAlgebra::operator*(const Matrix<R1>& A1, const MatrixSlice<R2>& A2) 
{
  return A1*Matrix<R2>(A2);
}

template<class R1, class R2>  inline
LinearAlgebra::Matrix<typename Numeric::traits<R1,R2>::arithmetic_type>
LinearAlgebra::operator*(const MatrixSlice<R1>& A1, const Matrix<R2>& A2) 
{
  return Matrix<R1>(A1)*A2;
}

template<class R1, class R2>  inline
LinearAlgebra::Matrix<typename Numeric::traits<R1,R2>::arithmetic_type>
LinearAlgebra::operator*(const MatrixSlice<R1>& A1, const MatrixSlice<R2>& A2) 
{
  return Matrix<R1>(A1)*Matrix<R2>(A2);
}



template<class R1, class R2>  inline
LinearAlgebra::Matrix<typename Numeric::traits<R1,R2>::arithmetic_type>
LinearAlgebra::operator*(const R1& s, const MatrixSlice<R2>& A) 
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

template<class R1, class R2> 
LinearAlgebra::Matrix<typename Numeric::traits<R1,R2>::arithmetic_type>
LinearAlgebra::operator*(const MatrixSlice<R1>& A, const R2& s) 
  {
  return s*A;
}

template<class R1, class R2>  inline
LinearAlgebra::Vector<typename Numeric::traits<R1,R2>::arithmetic_type>
LinearAlgebra::operator*(const MatrixSlice<R1>& A, const Vector<R2>& v) 
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


template<class R1, class R2>  inline
LinearAlgebra::Matrix<typename Numeric::traits<R1,R2>::arithmetic_type>
LinearAlgebra::outer_product(const Vector<R1>& v1, const Vector<R2>& v2) 
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





template<class R> inline
typename Numeric::traits<R>::arithmetic_type
LinearAlgebra::norm(const Matrix<R>& A) 
{
  return sup_norm(A);
}


template<class R> inline
LinearAlgebra::Vector<typename Numeric::traits<R>::arithmetic_type>
LinearAlgebra::row_norms(const Matrix<R>& A) 
{
  Vector<typename Numeric::traits<R>::arithmetic_type> result(A.number_of_rows());
  for(size_type i=0; i!=A.number_of_rows(); ++i) {
    for(size_type j=0; j!=A.number_of_columns(); ++j) {
      result(i)+=abs(A(i,j));
    }
  }
  return result;
}

template<class R> inline
typename Numeric::traits<R>::arithmetic_type
LinearAlgebra::norm(const MatrixSlice<R>& A) 
{
  return sup_norm(A);
}


template<class R> inline
LinearAlgebra::Matrix<R>
LinearAlgebra::concatenate_columns(const MatrixSlice<R>& A1, const MatrixSlice<R>& A2) {
  return concatenate_columns(Matrix<R>(A1),Matrix<R>(A2));
}



template<class R> inline
std::ostream&
LinearAlgebra::operator<<(std::ostream& os, const MatrixSlice<R>& A) {
  return A.write(os);
}

template<class R> inline
std::ostream&
LinearAlgebra::operator<<(std::ostream& os, const Matrix<R>& A) {
  return A.write(os);
}

template<class R> inline
std::istream&
LinearAlgebra::operator>>(std::istream& is, Matrix<R>& A) {
  return A.read(is); 
}


} // namespace Ariadne
