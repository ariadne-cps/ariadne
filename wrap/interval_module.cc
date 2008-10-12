/***************************************************************************
 *            file
 *
 *  Copyright 2008  Pieter Collins
 * 
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
 
/*! \file orbit.h
 *  \brief Orbits of dynamic systems
 */
#include "numeric.h"
#include "vector.h"
#include "matrix.h"

#include <boost/python.hpp>

using namespace boost::python;

using namespace Ariadne;

void read(Interval& ivl, const boost::python::object& obj) {
  boost::python::extract<boost::python::list> lx(obj);
  boost::python::extract<boost::python::tuple> tx(obj);
  boost::python::extract<double> dx(obj);
  boost::python::extract<Interval> ix(obj);

  if(lx.check()) {
    boost::python::list lst=lx();
    assert(boost::python::len(lst)==2);
    double l=boost::python::extract<double>(lst[0]);
    double u=boost::python::extract<double>(lst[1]);
    ivl=Interval(l,u);
  } else if(tx.check()) {
    boost::python::tuple tpl=tx();
    assert(boost::python::len(tpl)==2);
    double l=boost::python::extract<double>(tpl[0]);
    double u=boost::python::extract<double>(tpl[1]);
    ivl=Interval(l,u);
  } else if(dx.check()) {
    ivl=static_cast<Interval>(dx());
  } else if(ix.check()) {
    ivl=ix();
  }
}


/*
template<class X> class Vector
{
  uint _size; X* _ptr;
 public:
  ~Vector() { delete[] _ptr; }
  Vector(const boost::python::object& obt);
  Vector(uint n) : _size(n), _ptr(new X[n]) { 
    for(uint i=0; i!=_size; ++i) { _ptr[i]=0; } }
  Vector(uint n, const X* p) : _size(n), _ptr(new X[n]) { 
    for(uint i=0; i!=_size; ++i) { _ptr[i]=p[i]; } }
  Vector(const Vector<X>& v) : _size(v._size), _ptr(new X[_size]) { 
    for(uint i=0; i!=_size; ++i) { _ptr[i]=v.begin()[i]; } }
  Vector<X>& operator=(const Vector<X>& v) { 
    if(this!=&v) { delete[] _ptr; _size=v._size; _ptr=new X[_size]; 
      for(uint i=0; i!=_size; ++i) { _ptr[i]=v.begin()[i]; } } return *this; }
  void resize(uint n) { 
    delete[] _ptr; _size=n; _ptr=new X[n]; for(uint i=0; i!=n; ++i) { _ptr[i]=0; } }
  uint size() const { return _size; }
  const X* begin() const { return _ptr; }
  const X& get(uint i) const { return _ptr[i]; }
  template<class T> void set(uint i, const T& x) { _ptr[i] = x; }
  const X& operator[](uint i) const { return _ptr[i]; }
  X& operator[](uint i) { return _ptr[i]; }
};

template<class X> Vector<X> operator-(const Vector<X>& v) {
  Vector<X> r(v.size()); for(uint i=0; i!=r.size(); ++i) { r[i]=-v[i]; } return r; }
template<class X> Vector<X> operator+(const Vector<X>& v1, const Vector<X>& v2) {
  Vector<X> r(v1.size()); for(uint i=0; i!=r.size(); ++i) { r[i]=v1[i]+v2[i]; } return r; }
template<class X> Vector<X> operator-(const Vector<X>& v1, const Vector<X>& v2) {
  Vector<X> r(v1.size()); for(uint i=0; i!=r.size(); ++i) { r[i]=v1[i]-v2[i]; } return r; }
template<class X> Vector<X> operator*(const X& s, const Vector<X>& v) {
  Vector<X> r(v.size()); for(uint i=0; i!=r.size(); ++i) { r[i]=s*v[i]; } return r; }
template<class X> Vector<X> operator*(const Vector<X>& v, const X& s) {
  Vector<X> r(v.size()); for(uint i=0; i!=r.size(); ++i) { r[i]=v[i]*s; } return r; }
template<class X> Vector<X> operator/(const Vector<X>& v, const X& s) {
  Vector<X> r(v.size()); for(uint i=0; i!=r.size(); ++i) { r[i]=v[i]/s; } return r; }

template<class X> Vector<X> operator*(const double& s, const Vector<X>& v) {
  Vector<X> r(v.size()); for(uint i=0; i!=r.size(); ++i) { r[i]=s*v[i]; } return r; }
template<class X> Vector<X> operator*(const Vector<X>& v, const double& s) {
  Vector<X> r(v.size()); for(uint i=0; i!=r.size(); ++i) { r[i]=v[i]*s; } return r; }
template<class X> Vector<X> operator/(const Vector<X>& v, const double& s) {
  Vector<X> r(v.size()); for(uint i=0; i!=r.size(); ++i) { r[i]=v[i]/s; } return r; }

template<class X> std::ostream& operator<<(std::ostream& os, const Vector<X>& v) {
  for(uint i=0; i!=v.size(); ++i) { os << (i==0?"[":",") << v[i]; } return os << "]";
}
 
template<class X> Vector<X>::Vector(const boost::python::object& obj) {
  boost::python::list elements=boost::python::extract<boost::python::list>(obj);
  int n=boost::python::len(elements);
  this->_size=n; 
  this->_ptr=new X[n];
  for(int i=0; i!=n; ++i) {
    read(this->_ptr[i], elements[i]);
  }
}


template<class X> class Matrix
{
  uint _rows; uint _cols; X* _ptr;
 public:
  ~Matrix() { delete[] _ptr; }
  Matrix(const boost::python::object& obj);
  Matrix(uint r, uint c) : _rows(r), _cols(c), _ptr(new X[r*c]) { 
    for(uint i=0; i!=r*c; ++i) { _ptr[i]=0; } }
  Matrix(uint r, uint c, const X* p) : _rows(r), _cols(c), _ptr(new X[r*c]) { 
    for(uint i=0; i!=r*c; ++i) { _ptr[i]=p[i]; } }
  Matrix(const Matrix<X>& A) : _rows(A._rows), _cols(A._cols), _ptr(new X[_rows*_cols]) { 
    for(uint i=0; i!=_rows*_cols; ++i) { _ptr[i]=A.begin()[i]; } }
  Matrix<X>& operator=(const Matrix<X>& A) { 
    if(this!=&A) { delete[] _ptr; _rows=A._rows; _cols=A._cols; _ptr=new X[_rows*_cols]; 
      for(uint i=0; i!=_rows*_cols; ++i) { _ptr[i]=A.begin()[i]; } } return *this; }
  void resize(uint r, uint c) { delete[] _ptr; _rows=r; _cols=c; _ptr=new X[r*c]; 
    for(int i=0; i!=r*c; ++i) { _ptr[i]=0; } }
  uint rows() const { return _rows; } 
  uint columns() const { return _cols; } 
  const X* begin() const { return _ptr; }
  X* begin() { return _ptr; }
  const X& get(uint i, uint j) const { return _ptr[i*_cols+j]; }
  template<class T> void set(uint i, uint j, const T& x) { _ptr[i*_cols+j] = x; }
  const X* operator[](uint i) const { return _ptr+i*_cols; }
  X* operator[](uint i) { return _ptr+i*_cols; }
};

template<class X> std::ostream& operator<<(std::ostream& os, const Matrix<X>& a) {
  for(uint i=0; i!=a.rows(); ++i) { os << (i==0?"[":";"); 
    for(uint j=0; j!=a.columns(); ++j) { os << (j==0?"[":",") << a.get(i,j); } } return os << "]";
}
 
template<class X> Matrix<X> operator-(const Matrix<X>& a) {
  Matrix<X> r(a.rows(),a.columns()); for(uint i=0; i!=a.rows()*a.columns(); ++i) { r.begin()[i]=-a.begin()[i]; } return r; }
template<class X> Matrix<X> operator+(const Matrix<X>& a1, const Matrix<X>& a2) {
  Matrix<X> r(a1.rows(),a1.columns()); for(uint i=0; i!=r.rows()*r.columns(); ++i) { r.begin()[i]=a1.begin()[i]+a2.begin()[i]; } return r; }
template<class X> Matrix<X> operator-(const Matrix<X>& a1, const Matrix<X>& a2) {
  Matrix<X> r(a1.rows(),a1.columns()); for(uint i=0; i!=r.rows()*r.columns(); ++i) { r.begin()[i]=a1.begin()[i]-a2.begin()[i]; } return r; }
template<class X> Matrix<X> operator*(const X& s, const Matrix<X>& a) {
  Matrix<X> r(a.rows(),a.columns()); for(uint i=0; i!=a.rows()*a.columns(); ++i) { r.begin()[i]=s*a.begin()[i]; } return r; }
template<class X> Matrix<X> operator*(const Matrix<X>& a, const X& s) {
  Matrix<X> r(a.rows(),a.columns()); for(uint i=0; i!=a.rows()*a.columns(); ++i) { r.begin()[i]=a.begin()[i]*s; } return r; }
template<class X> Matrix<X> operator/(const Matrix<X>& a, const X& s) {
  Matrix<X> r(a.rows(),a.columns()); for(uint i=0; i!=a.rows()*a.columns(); ++i) { r.begin()[i]=a.begin()[i]/s; } return r; }

template<class X> Matrix<X> operator*(const double& s, const Matrix<X>& a) {
  Matrix<X> r(a.rows(),a.columns()); for(uint i=0; i!=a.rows()*a.columns(); ++i) { r.begin()[i]=s*a.begin()[i]; } return r; }
template<class X> Matrix<X> operator*(const Matrix<X>& a, const double& s) {
  Matrix<X> r(a.rows(),a.columns()); for(uint i=0; i!=a.rows()*a.columns(); ++i) { r.begin()[i]=a.begin()[i]*s; } return r; }
template<class X> Matrix<X> operator/(const Matrix<X>& a, const double& s) {
  Matrix<X> r(a.rows(),a.columns()); for(uint i=0; i!=a.rows()*a.columns(); ++i) { r.begin()[i]=a.begin()[i]/s; } return r; }

template<class X> Vector<X> operator*(const Matrix<X>& a, const Vector<X>& v) {
  Vector<X> r(a.rows()); 
  for(uint i=0; i!=a.rows(); ++i) { 
    for(uint j=0; j!=a.columns(); ++j) {
      r[i]+=a[i][j]*v[j]; 
    }
  }
  return r;
}

template<class X> Matrix<X> operator*(const Matrix<X>& a1, const Matrix<X>& a2) {
  Matrix<X> r(a1.rows(),a2.columns()); 
  for(uint i=0; i!=a1.rows(); ++i) { 
    for(uint j=0; j!=a2.columns(); ++j) {
      for(uint k=0; k!=a1.columns(); ++k) {
        r[i][j]+=a1[i][k]*a2[k][j];
      }
    }
  }
  return r;
}

template<class X> Matrix<X>::Matrix(const boost::python::object& obj) {
  boost::python::list elements=boost::python::extract<boost::python::list>(obj);
  int m=boost::python::len(elements);
  boost::python::list row=boost::python::extract<boost::python::list>(elements[0]);
  int n=boost::python::len(row);
  this->_rows=m;
  this->_cols=n;
  this->_ptr=new X[_rows*_cols];
  for(int i=0; i!=m; ++i) {
    row=boost::python::extract<boost::python::list>(elements[i]);
    if(boost::python::len(row)!=n) {
      throw std::runtime_error("Matrix with rows of different sizes");
    }
    for(int j=0; j!=n; ++j) {
      read(this->_ptr[i*_cols+j],row[j]);
    }
  }
}

*/

typedef Interval Ivl;
typedef Vector<Interval> Vec;
typedef Matrix<Interval> Mx;



void export_interval() 
{
    using boost::python::class_;
    using boost::python::init;
    using boost::python::self;
    using boost::python::return_value_policy;
    using boost::python::copy_const_reference;
    using boost::python::def;

    //typedef Interval (*IFUN)(const Interval&);
    //typedef Interval (*IZFUN)(const Interval&, int n);
    //typedef Interval (*INFUN)(const Interval&, unsigned int n);
    
    typedef Interval (*IFUN)(Interval);
    typedef Interval (*IZFUN)(Interval, int n);
    typedef Interval (*INFUN)(Interval, unsigned int n);
    
    def("down",&down);
    def("up",&up);

    class_< Interval > interval_class("Interval");
    interval_class.def(init<double,double>());
    interval_class.def(init<double>());
    interval_class.def(-self);
    interval_class.def(self + self);
    interval_class.def(self - self);
    interval_class.def(self * self);
    interval_class.def(self / self);
    interval_class.def(self + double());
    interval_class.def(self - double());
    interval_class.def(self * double());
    interval_class.def(self / double());
    interval_class.def(double() + self);
    interval_class.def(double() - self);
    interval_class.def(double() * self);
    interval_class.def(double() / self);
    interval_class.def("lower", &Interval::lower, return_value_policy<copy_const_reference>());
    interval_class.def("upper", &Interval::upper, return_value_policy<copy_const_reference>());
    interval_class.def(boost::python::self_ns::str(self));

    def("med", (IFUN) &med);
    def("rad", (IFUN) &rad);
    def("diam", (IFUN) &diam);

    def("trunc", (INFUN) &trunc, "truncate to n binary digits");
    def("abs", (IFUN) &abs, "interval absolute value function");
    def("pow",  (IZFUN) &pow, "interval power function");
    def("sqr", (IFUN) &sqr, "interval square function");
    def("sqrt", (IFUN) &sqrt);
    def("exp", (IFUN) &exp);
    def("log", (IFUN) &log);
    /*
    def("sin", isin);
    def("cos", icos);
    def("tan", itan);
    def("asin", iasin);
    def("acos", iacos);
    def("atan", iatan);
    def("sinh", isinh);
    def("cosh", icosh);
    def("tanh", itanh);
    def("asinh", iasinh);
    def("acosh", iacosh);
    def("atanh", iatanh);
    */
}

    

template<class X>
void export_vector()
{
    typedef Vector<X> V;

    class_< Vector<X> > vector_class("Vector",init<int>());
    //vector_class.def(init<const boost::python::object&>());
    vector_class.def("__len__", &Vector<X>::size);
    vector_class.def("__setitem__", (void(Vector<X>::*)(uint,const double&)) &Vector<X>::set);
    vector_class.def("__setitem__", (void(Vector<X>::*)(uint,const X&)) &Vector<X>::set);
    vector_class.def("__getitem__", &Vector<X>::get, return_value_policy<copy_const_reference>());
    vector_class.def(-self);
    vector_class.def(self + self);
    vector_class.def(self - self);
    vector_class.def(X() * self);
    vector_class.def(self * X());
    vector_class.def(self / X());
    vector_class.def(double() * self);
    vector_class.def(self * double());
    vector_class.def(self / double());
    vector_class.def(self += self);
    vector_class.def(self -= self);
    vector_class.def(self *= double());
    vector_class.def(self *= X());
    vector_class.def(self /= double());
    vector_class.def(self /= X());
    vector_class.def(boost::python::self_ns::str(self));
}


template<class X>
void export_matrix()
{
    class_< Matrix<X> > matrix_class("Matrix",init<int,int>());
    //matrix_class.def(init<const boost::python::object&>());
    matrix_class.def("rows", &Matrix<X>::row_size);
    matrix_class.def("columns", &Matrix<X>::column_size);
    matrix_class.def("row_size", &Matrix<X>::row_size);
    matrix_class.def("column_size", &Matrix<X>::column_size);
    matrix_class.def("__setitem__", (void(Matrix<X>::*)(uint,uint,const double&)) &Matrix<X>::set);
    matrix_class.def("__setitem__", (void(Matrix<X>::*)(uint,uint,const X&)) &Matrix<X>::set);
    matrix_class.def("__getitem__", &Matrix<X>::get, return_value_policy<copy_const_reference>());
    matrix_class.def(-self);        // __neg__
    matrix_class.def(self + self);  // __add__
    matrix_class.def(self - self);  // __sub__
    matrix_class.def(X() * self);  // __mul__
    matrix_class.def(self * X());  // __mul__
    matrix_class.def(self / X());  // __div__
    matrix_class.def(double() * self);  // __mul__
    matrix_class.def(self * double());  // __mul__
    matrix_class.def(self / double());  // __div__
    matrix_class.def("__mul__",(Vector<X>(*)(const Matrix<X>&,const Vector<X>&))&operator*);
    matrix_class.def("__mul__",(Matrix<X>(*)(const Matrix<X>&,const Matrix<X>&))&operator*);
    //matrix_class.def("inverse", &inverse<X>);
    //matrix_class.def("determinant", &Matrix::determinant);
    //matrix_class.def("transpose", &Matrix::transpose);
    //matrix_class.def("solve", &Matrix::solve);
    matrix_class.def(boost::python::self_ns::str(self));    // __str__
}



BOOST_PYTHON_MODULE(interval)
{
  export_interval();
  export_vector<Interval>();
  export_matrix<Interval>();
}

      

