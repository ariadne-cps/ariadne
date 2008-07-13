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
Matrix<R>::Matrix()
  : _nr(0), _nc(0), _array() 
{ 
}


template<class R> inline
Matrix<R>::Matrix(const size_type& r, const size_type& c)
  : _nr(r), _nc(c), _array(r*c,static_cast<R>(0)) 
{ 
}


template<class R> template<class RR> inline
Matrix<R>::Matrix(const size_type& nr, const size_type& nc, const RR* ptr)
  
  : _nr(nr), _nc(nc), _array(ptr,ptr+nr*nc) 
{ 
}

template<class R> template<int NN, class RR> inline
Matrix<R>::Matrix(const size_type& nr, const size_type& nc, const RR ary[][NN])
  
  : _nr(nr), _nc(nc), _array(nr*nc) 
{ 
  for(size_type i=0; i!=nr; ++i) { 
    for(size_type j=0; j!=nc; ++j) { 
      (*this)(i,j)=static_cast<RR>(ary[i][j]);
    } 
  } 
}


template<class R> template<class RR> inline
Matrix<R>::Matrix(const size_type& nr, const size_type& nc, 
                                 const RR* ptr, const size_type& ri, const size_type& ci)
  : _nr(nr), _nc(nc), _array(nr*nc) 
{ 
  for(size_type i=0; i!=nr; ++i) { 
    for(size_type j=0; j!=nc; ++j) { 
      (*this)(i,j)=ptr[i*ri+j*ci];
    } 
  } 
}


template<class R> template<class RR> inline
Matrix<R>::Matrix(const Matrix<RR>& A)
  : _nr(A.number_of_rows()), _nc(A.number_of_columns()), _array(_nc*_nr)
{ 
  for(size_type i=0; i!=_nr; ++i) { 
    for(size_type j=0; j!=_nc; ++j) { 
      (*this)(i,j)=R(A(i,j));
    } 
  } 
}

template<class R> template<class E> inline
Matrix<R>::Matrix(const MatrixExpression<E>& A)
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
Matrix<R>::Matrix(const Matrix<R>& A)
  : _nr(A._nr), _nc(A._nc), _array(A._array) 
{ 
}


template<class R> inline
Matrix<R>& 
Matrix<R>::operator=(const Matrix<R>& A) 
{
  if(this!=&A) { 
    this->_nr=A._nr; 
    this->_nc=A._nc; 
    this->_array=A._array;
  }
  return *this; 
}


template<class R> template<class E> inline
Matrix<R>& 
Matrix<R>::operator=(const MatrixExpression<E>& A) 
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
Matrix<R>::operator==(const Matrix<R>& A) const 
{ 
  return this->_nr==A._nr && this->_array==A._array;
}


template<class R> inline
bool 
Matrix<R>::operator!=(const Matrix<R>& A) const 
{
  return !(*this==A);
}


template<class R> inline
array<R>& 
Matrix<R>::data() 
{ 
  return this->_array; 
}


template<class R> inline
const array<R>& 
Matrix<R>::data() const 
{ 
  return this->_array;
}


template<class R> inline
void 
Matrix<R>::resize(const size_type& nr, const size_type nc) 
{
  if(nr*nc!=this->_array.size()) {
    _array.resize(nr*nc);
  }
  this->_nr=nr;
  this->_nc=nc;
}


template<class R> inline
size_type 
Matrix<R>::number_of_rows() const 
{ 
  return this->_nr; 
}


template<class R> inline
size_type 
Matrix<R>::number_of_columns() const 
{ 
  return this->_nc; 
}


template<class R> inline
size_type 
Matrix<R>::row_size() const 
{ 
  return this->_nr; 
}


template<class R> inline
size_type 
Matrix<R>::column_size() const 
{ 
  return this->_nc; 
}

template<class R> inline
array<size_type,2u> 
Matrix<R>::size() const 
{ 
  return array<size_type,2u>(this->_nr,this->_nc); 
}



template<class R> inline
MatrixSlice<R>
Matrix<R>::operator() (const Slice& i, const Slice& j)  
{ 
  R* begin=this->_array.begin()+this->row_increment()*i.start()+this->column_increment()*j.start();
  return MatrixSlice<R>(i.size(),j.size(),begin,
                        this->row_increment()*i.stride(), this->column_increment()*j.stride());
}

template<class R> inline
MatrixSlice<const R> 
Matrix<R>::operator() (const Slice& i, const Slice& j) const
{ 
  const R* begin=this->_array.begin()+this->row_increment()*i.start()+this->column_increment()*j.start();
  return MatrixSlice<const R>(i.size(),j.size(),begin,
                              this->row_increment()*i.stride(), this->column_increment()*j.stride());
}


template<class R> inline
const R& 
Matrix<R>::operator() (const size_type& i, const size_type& j) const 
{ 
  return this->_array[i*this->_nc+j]; 
}

template<class R> inline
R& 
Matrix<R>::operator() (const size_type& i, const size_type& j) 
{ 
  return this->_array[i*this->_nc+j]; 
}


template<class R> inline
MatrixRow< Matrix<R> > 
Matrix<R>::operator[](const size_type& i) 
{
  return MatrixRow< Matrix<R> >(*this,i); 
}        

template<class R> inline
MatrixRow< const Matrix<R> > 
Matrix<R>::operator[](const size_type& i) const
{
  return MatrixRow< const Matrix<R> >(*this,i); 
}        


template<class R> inline
Matrix<R>
Matrix<R>:: zero(const size_type r, const size_type c) 
{
  return Matrix<R>(r,c); 
}


template<class R> inline
Matrix<R> 
Matrix<R>::one(const size_type r, const size_type c) 
{
  R one(1); return Matrix<R>(r,c,&one,0,0); 
}


template<class R> inline
Matrix<R> 
Matrix<R>::identity(const size_type n) 
{
  Matrix<R> result(n,n); 
  for(size_type i=0; i!=n; ++i) { 
    result(i,i)=R(1); 
  } 
  return result;
}


template<class R> inline
VectorSlice<const R>
Matrix<R>::row(const size_type& i) const 
{
  return VectorSlice<const R>(this->number_of_columns(),this->begin()+i*this->row_increment(),this->column_increment()); 
}


template<class R> inline
VectorSlice<const R> 
Matrix<R>::column(const size_type& j) const 
{
  return VectorSlice<const R>(this->number_of_rows(),this->begin()+j*this->column_increment(),this->row_increment());
}

template<class R> inline
VectorSlice<R>
Matrix<R>::row(const size_type& i)  
{
  return VectorSlice<R>(this->number_of_columns(),this->begin()+i*this->row_increment(),this->column_increment()); 
}


template<class R> inline
VectorSlice<R> 
Matrix<R>::column(const size_type& j)  
{
  return VectorSlice<R>(this->number_of_rows(),this->begin()+j*this->column_increment(),this->row_increment());
}


template<class R> inline
R* 
Matrix<R>::begin() 
{
  return this->data().begin(); 
}


template<class R> inline
const R* 
Matrix<R>::begin() const 
{
  return this->data().begin(); 
}


template<class R> inline
size_type 
Matrix<R>::row_increment() const 
{
  return this->number_of_columns(); 
}


template<class R> inline
size_type 
Matrix<R>::column_increment() const 
{ 
  return 1u; 
}









template<class R> inline 
Matrix<R>
midpoint(const Matrix< Interval<R> >& im) 
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
Matrix<R>
approximation(const Matrix<R>& A) 
{
  return A;
}

template<class R> inline 
Matrix<R>
approximation(const Matrix< Interval<R> >& iA) 
{
  return midpoint(iA);
}

template<class R> inline 
bool
encloses(const Matrix< Interval<R> >& iA, const Matrix<R>& A) 
{
  ARIADNE_CHECK_MATRIX_EQUAL_SIZES(iA,A,"bool encloses(Matrix<Interval>,Matrix<Real>)");
  for(size_type i=0; i!=A.number_of_rows(); ++i) {
    for(size_type j=0; j!=A.number_of_columns(); ++j) {
      if(!encloses(iA(i,j),A(i,j))) {
        return false;
      }
    }
  }
  return true;
}


template<class R> inline 
bool
refines(const Matrix< Interval<R> >& iA1, const Matrix< Interval<R> >& iA2) 
{
  ARIADNE_CHECK_MATRIX_EQUAL_SIZES(iA1,iA2,"bool refines(Matrix<Interval>,Matrix<Interval>)");
  for(size_type i=0; i!=iA1.number_of_rows(); ++i) {
    for(size_type j=0; j!=iA1.number_of_columns(); ++j) {
      if(!refines(iA1(i,j),iA2(i,j))) {
        return false;
      }
    }
  }
  return true;
}



template<class R> inline
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


template<class R> inline
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
Matrix<R1>&
operator-=(Matrix<R1>& A1, const Matrix<R2>& A2) 
{
  for(size_type i=0; i!=A1.number_of_rows(); ++i) {
    for(size_type j=0; j!=A1.number_of_columns(); ++j) {
      A1(i,j)-=A2(i,j);
    }
  }
  return A1;
}


template<class R1, class R2>  inline
Matrix<typename traits<R1,R2>::arithmetic_type>
operator+(const Matrix<R1>& A1, const Matrix<R2>& A2) 
{
  typedef typename traits<R1,R2>::arithmetic_type R3;
  Matrix<R3> result(A1.number_of_rows(),A1.number_of_columns());
  for(size_type i=0; i!=A1.number_of_rows(); ++i) {
    for(size_type j=0; j!=A1.number_of_columns(); ++j) {
      result(i,j)=A1(i,j)+A2(i,j);
    }
  }
  return result;
}


template<class R1, class R2>  inline
Matrix<typename traits<R1,R2>::arithmetic_type>
operator-(const Matrix<R1>& A1, const Matrix<R2>& A2) 
{
  typedef typename traits<R1,R2>::arithmetic_type R3;
  Matrix<R3> result(A1.number_of_rows(),A1.number_of_columns());
  for(size_type i=0; i!=A1.number_of_rows(); ++i) {
    for(size_type j=0; j!=A1.number_of_columns(); ++j) {
      result(i,j)=A1(i,j)-A2(i,j);
    }
  }
  return result;
}


template<class R1, class R2>  inline
Matrix<typename traits<R1,R2>::arithmetic_type>
operator*(const R1& s, const Matrix<R2>& A) 
{
  typedef typename traits<R1,R2>::arithmetic_type R3;
  Matrix<R3> result(A.number_of_rows(),A.number_of_columns());
  for(size_type i=0; i!=A.number_of_rows(); ++i) {
    for(size_type j=0; j!=A.number_of_columns(); ++j) {
      result(i,j)=s*A(i,j);
    }
  }
  return result;
}

template<class R1, class R2>  inline
Matrix<typename traits<R1,R2>::arithmetic_type>
operator*(const Matrix<R1>& A, const R2& s) 
{
  return s*A;
}

template<class R1, class R2>  inline
Matrix<typename traits<R1,R2>::arithmetic_type>
operator/(const Matrix<R1>& A, const R2& s) 
{
  return static_cast<typename traits<R2>::arithmetic_type>(static_cast<R2>(1)/s)*A;
}

template<class R1, class R2>  inline
Vector<typename traits<R1,R2>::arithmetic_type>
operator*(const Matrix<R1>& A, const Vector<R2>& v) 
{
  typedef typename traits<R1,R2>::arithmetic_type R3;
  Vector<R3> result(A.number_of_rows());
  for(size_type i=0; i!=A.number_of_rows(); ++i) {
    for(size_type j=0; j!=A.number_of_columns(); ++j) {
      result(i)+=A(i,j)*v(j);
    }
  }
  return result;
}

template<class R1, class R2>  inline
Vector<typename traits<R1,R2>::arithmetic_type>
operator*(const Vector<R1>& v, const Matrix<R2>& A) 
{
  typedef typename traits<R1,R2>::arithmetic_type R3;
  Vector<R3> result(A.number_of_columns());
  for(size_type i=0; i!=A.number_of_rows(); ++i) {
    for(size_type j=0; j!=A.number_of_columns(); ++j) {
      result(j)+=v(i)*A(i,j);
    }
  }
  return result;
}

template<class R1, class R2>  inline
Matrix<typename traits<R1,R2>::arithmetic_type>
operator*(const Matrix<R1>& A1, const Matrix<R2>& A2) 
{
  typedef typename traits<R1,R2>::arithmetic_type R3;
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
Matrix<typename traits<R1,R2>::arithmetic_type>
operator*(const Matrix<R1>& A1, const MatrixSlice<R2>& A2) 
{
  return A1*Matrix<R2>(A2);
}

template<class R1, class R2>  inline
Matrix<typename traits<R1,R2>::arithmetic_type>
operator*(const MatrixSlice<R1>& A1, const Matrix<R2>& A2) 
{
  return Matrix<R1>(A1)*A2;
}

template<class R1, class R2>  inline
Matrix<typename traits<R1,R2>::arithmetic_type>
operator*(const MatrixSlice<R1>& A1, const MatrixSlice<R2>& A2) 
{
  return Matrix<R1>(A1)*Matrix<R2>(A2);
}



template<class R1, class R2>  inline
Matrix<typename traits<R1,R2>::arithmetic_type>
operator*(const R1& s, const MatrixSlice<R2>& A) 
{
  typedef typename traits<R1,R2>::arithmetic_type R3;
  Matrix<R3> result(A.number_of_rows(),A.number_of_columns());
  for(size_type i=0; i!=A.number_of_rows(); ++i) {
    for(size_type j=0; j!=A.number_of_columns(); ++j) {
      result(i,j)=s*A(i,j);
    }
  }
  return result;
}

template<class R1, class R2> 
Matrix<typename traits<R1,R2>::arithmetic_type>
operator*(const MatrixSlice<R1>& A, const R2& s) 
  {
  return s*A;
}

template<class R1, class R2>  inline
Vector<typename traits<R1,R2>::arithmetic_type>
operator*(const MatrixSlice<R1>& A, const Vector<R2>& v) 
{
  typedef typename traits<R1,R2>::arithmetic_type R3;
  Vector<R3> result(A.number_of_rows());
  for(size_type i=0; i!=A.number_of_rows(); ++i) {
    for(size_type j=0; j!=A.number_of_columns(); ++j) {
      result(i)+=A(i,j)*v(j);
    }
  }
  return result;
}


template<class R1, class R2>  inline
Matrix<typename traits<R1,R2>::arithmetic_type>
outer_product(const Vector<R1>& v1, const Vector<R2>& v2) 
{
  typedef typename traits<R1,R2>::arithmetic_type R3;
  Matrix<R3> result(v1.size(),v2.size());
  for(size_type i=0; i!=v1.size(); ++i) {
    for(size_type j=0; j!=v2.size(); ++j) {
      result(i,j)=v1(i)*v2(j);
    }
  }
  return result;
}





template<class R> inline
typename traits<R>::arithmetic_type
norm(const Matrix<R>& A) 
{
  return sup_norm(A);
}


template<class R> inline
Vector<typename traits<R>::arithmetic_type>
row_norms(const Matrix<R>& A) 
{
  Vector<typename traits<R>::arithmetic_type> result(A.number_of_rows());
  for(size_type i=0; i!=A.number_of_rows(); ++i) {
    for(size_type j=0; j!=A.number_of_columns(); ++j) {
      result(i)+=abs(A(i,j));
    }
  }
  return result;
}

template<class R> inline
typename traits<R>::arithmetic_type
norm(const MatrixSlice<R>& A) 
{
  return sup_norm(A);
}


template<class R> inline
Matrix<R>
concatenate_columns(const MatrixSlice<R>& A1, const MatrixSlice<R>& A2) {
  return concatenate_columns(Matrix<R>(A1),Matrix<R>(A2));
}



template<class R> inline
std::ostream&
operator<<(std::ostream& os, const Matrix<R>& A) {
  return A.write(os);
}

template<class R> inline
std::istream&
operator>>(std::istream& is, Matrix<R>& A) {
  return A.read(is); 
}



} // namespace Ariadne
