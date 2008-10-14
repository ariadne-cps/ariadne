/***************************************************************************
 *            function.h
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
 
/*! \file function.h
 *  \brief Concrete function types.
 */
#ifndef ARIADNE_FUNCTION_H
#define ARIADNE_FUNCTION_H

#include <iosfwd>
#include <iostream>
#include "function_interface.h"
#include "vector.h"

namespace Ariadne {

static const int SMOOTH=255;

// A wrapper for transformations  
// This class is for internal use only; we can easily specify parameter types.
template<class T> 
class FunctionBase
  : public FunctionInterface 
{
  T t;
 protected:
  FunctionBase() : t() { }
  template<class S> FunctionBase(const S& s) : t(s) { }
  template<class S1, class S2> FunctionBase(const S1& s1, const S2& s2) : t(s1,s2) { }
  template<class S1, class S2, class S3> FunctionBase(const S1& s1, const S2& s2, const S3& s3) : t(s1,s2,s3) { }
 public:
  virtual FunctionBase<T>* clone() const { return new FunctionBase<T>(*this); }
  virtual uint result_size() const { return t.result_size(); }
  virtual uint argument_size() const { return t.argument_size(); }
  virtual ushort smoothness() const { return t.smoothness(); }
  virtual Vector<Float> evaluate(const Vector<Float>& x) const {
    Vector<Float> r; t.compute(r,x); return r; } 
  virtual Vector<Interval> evaluate(const Vector<Interval>& x) const {
    Vector<Interval> r; t.compute(r,x); return r; }                          
  virtual Matrix<Float> jacobian(const Vector<Float>& x) const {
    return this->_expansion(x,1u).jacobian(); }                         
  virtual Matrix<Interval> jacobian(const Vector<Interval>& x) const {
    return this->_expansion(x,1u).jacobian(); }                         
  virtual DifferentialVector< SparseDifferential<Float> > expansion(const Vector<Float>& x, const ushort& s) const {
    return this->_expansion(x,s); } 
  virtual DifferentialVector< SparseDifferential<Interval> > expansion(const Vector<Interval>& x, const ushort& s) const {
    return this->_expansion(x,s); }
  virtual std::ostream& write(std::ostream& os) const  {
    return os << "Function"; }
 private:
  template<class X> DifferentialVector< SparseDifferential<X> > _expansion(const Vector<X>& x, const ushort& s) const {
    const uint& rs=t.result_size();
    const uint& as=t.argument_size();
    DifferentialVector< SparseDifferential<X> > dx(as,as,s);
    DifferentialVector< SparseDifferential<X> > dr(rs,as,s);
    for(uint i=0; i!=as; ++i) { dx[i]=x[i]; }
    for(uint i=0; i!=as; ++i) { dx[i][i+1]=1; }
    t.compute(dr,dx);
    return dr;
  }                                             
};


template<uint RS, uint AS, uint PS=0u, uint SM=255u>
struct FunctionData 
{
  const uint result_size() const { return RS; }
  const uint argument_size() const { return AS; }
  const uint parameter_size() const { return PS; }
  const uint smoothness() const { return SM; }
};
  
// FIXME: Make interval computations use interval version of parameters!
template<class T> struct FunctionWrapper
  : public T
{
  FunctionWrapper() : T(), _p(Vector<Float>(T::parameter_size())) { }
  FunctionWrapper(const Vector<Float>& p) : T(), _p(p) { ARIADNE_ASSERT(p.size()==T::parameter_size()); }
  FunctionWrapper(const Vector<Interval>& p) : T(), _p(midpoint(p)) { ARIADNE_ASSERT(p.size()==T::parameter_size()); }
  template<class R, class A> void compute(R& r, const A& x) const { this->T::compute(r,x,_p); }
 private:
  Vector<Float> _p;
};

template<class T> class Function
  : public FunctionBase< FunctionWrapper<T> >
{
 public:
  Function() : FunctionBase< FunctionWrapper<T> >() { }
  Function(const Vector<Float>& p) : FunctionBase< FunctionWrapper<T> >(p) { }
  Function(const Vector<Interval>& p) : FunctionBase< FunctionWrapper<T> >(p) { }
};



struct ConstantTransformation 
{
  ConstantTransformation(const Vector<Float>& c, uint as)
    : _as(as), _c(c) { }
  const uint result_size() const { return _c.size(); }
  const uint argument_size() const { return _as; }
  const int smoothness() const { return SMOOTH; }
  template<class R, class A>
  void compute(R& r, const A& x) const {
    for(uint i=0; i!=result_size(); ++i) { r[i]=_c[i]; } }
 private:
  uint _as;
  Vector<Float> _c;
};

struct IdentityTransformation 
{
  IdentityTransformation(uint n)
    : _n(n) { }
  const uint result_size() const { return _n; }
  const uint argument_size() const { return _n; }
  const int smoothness() const { return SMOOTH; }
  template<class R, class A>
  void compute(R& r, const A& x) const {
    for(uint i=0; i!=result_size(); ++i) { r[i]=x[i]; } }
 private:
  uint _n;
};

struct AffineTransformation 
{
  AffineTransformation(const Matrix<Float>& A, const Vector<Float>& b)
    : _A(A), _b(b) { ARIADNE_ASSERT(A.row_size()==b.size()); }
  const uint result_size() const { return _A.row_size(); }
  const uint argument_size() const { return _A.column_size(); }
  const int smoothness() const { return SMOOTH; }
  
  template<class R, class A>
  void compute(R& r, const A& x) const {
    for(uint i=0; i!=result_size(); ++i) {
      r[i]=_b[i]; 
      for(uint j=0; j!=argument_size(); ++j) {
        r[i]+=_A[i][j]*x[j];
      }
    }
  }
 private:
  Matrix<Float> _A;
  Vector<Float> _b;
};



//! A constant function \f$ x\mapsto c\f$ from \f$R^m\f$ to \f$R^n\f$.
class ConstantFunction
  : public FunctionBase<ConstantTransformation>
{
  //! A constant function with value \f$c\in\R^n\f$ and domain \f$\R^m\f$.
  ConstantFunction(const Vector<Float>& c, uint m) 
    : FunctionBase<ConstantTransformation>(c,m) { }

};


//! The identity function \f$ x\mapsto x\f$ in \f$\R^n\f$.
class IdentityFunction
  : public FunctionBase<IdentityTransformation>
{
  //! Construct the identity function in dimension \a n.
  IdentityFunction(uint n) 
    : FunctionBase<IdentityTransformation>(n) { }
};

//! An affine function \f$x\mapsto Ax+b\f$. 
class AffineFunction
  : public FunctionBase<AffineTransformation>
{
  //! Construct an affine function from the matrix \a A and vector \a b.
  AffineFunction(const Matrix<Float>& A, const Vector<Float>& b) 
    : FunctionBase<AffineTransformation>(A,b) { }
};


} // namespace Ariadne

#endif
