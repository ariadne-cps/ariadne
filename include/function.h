#ifndef ARIADNE_FUNCTION_H
#define ARIADNE_FUNCTION_H

#include <iosfwd>
#include <iostream>
#include "function_interface.h"
#include "vector.h"

namespace Ariadne {


template<class T> 
class Function : public FunctionInterface 
{
  T t; Vector<Float> p;
 public:
  Function(const Vector<Float>& parameters) : t(), p(parameters) { }
  virtual Function<T>* clone() const { return new Function<T>(*this); }
  virtual uint result_size() const { return t.result_size; }
  virtual uint argument_size() const { return t.argument_size; }
  virtual ushort smoothness() const { return t.smoothness; }
  virtual Vector<Float> evaluate(const Vector<Float>& x) const {
    Vector<Float> r; t.compute(r,x,p); return r; } 
  virtual Vector<Interval> evaluate(const Vector<Interval>& x) const {
    Vector<Interval> r; t.compute(r,x,p); return r; }                          
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
    const uint& rs=t.result_size;
    const uint& as=t.argument_size;
    DifferentialVector< SparseDifferential<X> > dx(as,as,s);
    DifferentialVector< SparseDifferential<X> > dr(rs,as,s);
    for(uint i=0; i!=as; ++i) { dx[i]=x[i]; }
    for(uint i=0; i!=as; ++i) { dx[i][i+1]=1; }
    t.compute(dr,dx,p);
    return dr;
  }                                             \
};


} // namespace Ariadne

#endif
