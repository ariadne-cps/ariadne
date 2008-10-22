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
#include "matrix.h"
#include "sparse_differential.h"
#include "differential_vector.h"

namespace Ariadne {

static const int SMOOTH=255;

// A wrapper for transformations  
// This class is for internal use only; we can easily specify parameter types.
template<class T> 
class FunctionBase
  : public FunctionInterface,
    public T
{
 protected:
  FunctionBase() : T() { }
  template<class S> FunctionBase(const S& s) : T(s) { }
  template<class S1, class S2> FunctionBase(const S1& s1, const S2& s2) : T(s1,s2) { }
  template<class S1, class S2, class S3> FunctionBase(const S1& s1, const S2& s2, const S3& s3) : T(s1,s2,s3) { }
 public:
  virtual FunctionBase<T>* clone() const { return new FunctionBase<T>(*this); }
  virtual uint result_size() const { return this->T::result_size(); }
  virtual uint argument_size() const { return this->T::argument_size(); }
  virtual ushort smoothness() const { return this->T::smoothness(); }
  virtual Vector<Float> evaluate(const Vector<Float>& x) const {
    Vector<Float> r(this->T::result_size()); this->T::compute(r,x); return r; } 
  virtual Vector<Interval> evaluate(const Vector<Interval>& x) const {
    Vector<Interval> r(this->T::result_size()); this->T::compute(r,x); return r; }                          
  virtual Matrix<Float> jacobian(const Vector<Float>& x) const {
    return this->_expansion(x,1u).get_jacobian(); }                         
  virtual Matrix<Interval> jacobian(const Vector<Interval>& x) const {
    return this->_expansion(x,1u).get_jacobian(); }                         
  virtual DifferentialVector< SparseDifferential<Float> > expansion(const Vector<Float>& x, const ushort& s) const {
    return this->_expansion(x,s); } 
  virtual DifferentialVector< SparseDifferential<Interval> > expansion(const Vector<Interval>& x, const ushort& s) const {
    return this->_expansion(x,s); }
  // TODO: Find a better way for writing functions which can handle transformations which may not have a
  // write() method or operator<<.
  virtual std::ostream& write(std::ostream& os) const  {
    return os << "Function( result_size="<<this->result_size()<<", argument_size="<<this->argument_size()<<" )"; }
 private:
  template<class X> DifferentialVector< SparseDifferential<X> > _expansion(const Vector<X>& x, const ushort& s) const {
    const uint rs=this->T::result_size();
    const uint as=this->T::argument_size();
    DifferentialVector< SparseDifferential<X> > dx(as,as,s);
    DifferentialVector< SparseDifferential<X> > dr(rs,as,s);
    for(uint i=0; i!=as; ++i) { dx[i]=x[i]; }
    for(uint i=0; i!=as; ++i) { dx[i][i]=1; }
    this->T::compute(dr,dx);
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
  const Vector<Float>& parameters() const { return this->_p; }
  template<class R, class A> void compute(R& r, const A& x) const { this->T::compute(r,x,_p); }
 private:
  Vector<Float> _p;
};

//! \brief A wrapper for converting templated C++ functions to %Ariadne functions.
//!
//! Given a C++ class T with a (static) template method 
//!   <code>template<class R, class A, class P> compute(R& r, const A& a, const P& p);</code>
//! the type <code>Function<T></code> is an Ariadne function defined by \f$r=f(a)\f$. 
//! The constructor for Function<T> takes a Vector<Float> argument which is used for \a p.
//! 
//! The class T must also define meta-data <c>result_size(), argument_size(), parameter_size()
//! and smoothness()</c>. These are most easily defined by inheriting from the 
//! <tt>FunctionData<RS,AS,PS,SM=SMOOTH></tt> class.
//! 
//! The constant \a SMOOTH is used for an arbitrarily-differentiable function.
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
  const Vector<Float>& c() const { return _c; }
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

struct ProjectionTransformation 
{
  ProjectionTransformation(const array<uint>& p, uint as) 
    : _as(as), _p(p) 
  { for(uint i=0; i!=_p.size(); ++i) { ARIADNE_ASSERT p[i]<as; } } 
  ProjectionTransformation(const Range& rng, uint as) 
    : _as(as), _p(rng.size()) 
  { ARIADNE_ASSERT range.start()+range.size()<=as;
    for(uint i=0; i!=_p.size(); ++i) { _p[i]=rng.start()+i; } }
  const array<uint>& p() const { return _p; }
  const uint result_size() const { return _p.size(); }
  const uint argument_size() const { return _as; }
  const int smoothness() const { return SMOOTH; }
  template<class R, class A>
  void compute(R& r, const A& x) const {
    for(uint i=0; i!=result_size(); ++i) { r[i]=x[p[i]]; } }
 private:
  uint _as;
  array<uint> _p;
};

struct ScalingTransformation 
{
  ScalingTransformation(const Vector<Float>& origin, 
                        const Vector<Float>& lengths)
    : _o(origin), _l(lengths) { ARIADNE_ASSERT(origin.size()==lengths.size()); }
  explicit ScalingTransformation(const Vector<Interval>& range)
    : _o(midpoint(range)), _l(range.size()) { for(uint i=0; i!=_l.size(); ++i) { _l[i]=range[i].radius(); } }
  const Vector<Float>& origin() const { return _o; }
  const Vector<Float>& lengths() const { return _l; }
  const uint result_size() const { return _l.size(); }
  const uint argument_size() const { return _l.size(); }
  const int smoothness() const { return SMOOTH; }
  template<class R, class A>
  void compute(R& r, const A& x) const {
    for(uint i=0; i!=result_size(); ++i) {
      r[i]=_o[i]+_l[i]*x[i]; 
    }
  }
 private:
  Vector<Float> _o;
  Vector<Float> _l;
};


struct AffineTransformation 
{
  AffineTransformation(const Matrix<Float>& A, const Vector<Float>& b)
    : _A(A), _b(b) { ARIADNE_ASSERT(A.row_size()==b.size()); }
  const Matrix<Float>& A() const { return _A; }
  const Vector<Float>& b() const { return _b; }
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
 public:
 //! A constant function with value \f$c\in\R^n\f$ and domain \f$\R^m\f$.
  ConstantFunction(const Vector<Float>& c, uint m) 
    : FunctionBase<ConstantTransformation>(c,m) { }
  std::ostream& write(std::ostream& os) const {
    return os << "ConstantFunction( argument_size=" << this->argument_size() 
              << ", c=" << this->c() << " )"; }

};


//! A projection function \f$ x'_i= x_{p(i)}\f$.
class ProjectionFunction
  : public FunctionBase<IdentityTransformation>
{
 public:
  //! Construct the identity function in dimension \a n.
  ProjectionFunction(const array<uint>& p, uint as) 
    : FunctionBase<ProjectionTransformation>(_p,as) { }
  std::ostream& write(std::ostream& os) const {
    return os << "ProjectionFunction( p=" << this->p() << " )"; }
};


//! The identity function \f$ x\mapsto x\f$ in \f$\R^n\f$.
class IdentityFunction
  : public FunctionBase<IdentityTransformation>
{
 public:
  //! Construct the identity function in dimension \a n.
  IdentityFunction(uint n) 
    : FunctionBase<IdentityTransformation>(n) { }
  std::ostream& write(std::ostream& os) const {
    return os << "IdentityFunction( size=" << this->result_size() << " )"; }
};


//! An scaling function \f$x_i' = o_i+l_ix_i\f$. 
class ScalingFunction
  : public FunctionBase<ScalingTransformation>
{
 public:
  //! Construct an affine function from the matrix \a A and vector \a b.
  ScalingFunction(const Vector<Float>& o, const Vector<Float>& l) 
    : FunctionBase<ScalingTransformation>(o,l) { }
  explicit ScalingFunction(const Vector<Interval>& bx) 
    : FunctionBase<ScalingTransformation>(bx) { }
  std::ostream& write(std::ostream& os) const {
    return os << "ScalingFunction( o=" << this->origin() << ", l=" << this->lengths() << " )"; }
};


//! An affine function \f$x\mapsto Ax+b\f$. 
class AffineFunction
  : public FunctionBase<AffineTransformation>
{
 public:
  //! Construct an affine function from the matrix \a A and vector \a b.
  AffineFunction(const Matrix<Float>& A, const Vector<Float>& b) 
    : FunctionBase<AffineTransformation>(A,b) { }
  std::ostream& write(std::ostream& os) const {
    return os << "AffineFunction( A=" << this->A() << ", b=" << this->b() << " )"; }
};


} // namespace Ariadne

#endif
