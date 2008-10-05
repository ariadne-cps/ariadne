#ifndef ARIADNE_SET_H
#define ARIADNE_SET_H

#include <iosfwd>

#include "tribool.h"

namespace Ariadne {

class Interval;
template<class X> class Vector;

class SetInterface {
 public:
  virtual ~SetInterface() { };
  virtual SetInterface* clone() const = 0;
  virtual uint dimension() const = 0;
  virtual tribool disjoint(const Vector<Interval>& bx) const = 0;
  virtual std::ostream& write(std::ostream& os) const = 0;
};

inline std::ostream& operator<<(std::ostream& os, const SetInterface& s) {
  return s.write(os); 
}


} // namespace Ariadne



#define ARIADNE_BUILD_FUNCTION(Nm,f,rs,as,np,sm)        \
  class Nm                                              \
    : public FunctionInterface                          \
  {                                                     \
  public:                                               \
    explicit Nm(const Vector<Float>& p) : _p(p) { }     \
    virtual Nm* clone() const { return new Nm(*this); }                 \
    virtual uint result_size() const { return rs; }                     \
    virtual uint argument_size() const { return as; }                   \
    virtual ushort smoothness() const { return sm; }                    \
    virtual Vector<Interval> evaluate(const Vector<Interval>& x, const ushort& s) const { Vector<Interval> r; f(r,x,_p); return r; } \
    virtual Vector<Interval> evaluate(const Vector<Interval>& x) const { Vector<Interval> r; f(r,x,_p); return r; } \
    virtual Matrix<Interval> jacobian(const Vector<Interval>& x) const { return this->_expansion(x,1u).jacobian(); } \
    virtual SparseDifferentialVector<Float> expansion(const Vector<Float>& x, const ushort& s) const { return this->_expansion(x,s); } \
    virtual SparseDifferentialVector<Interval> expansion(const Vector<Interval>& x, const ushort& s) const { return this->_expansion(x,s); } \
    virtual std::ostream& write(std::ostream& os) const  { return os << "Nm"; } \
  private:                                                              \
    template<class X> SparseDifferentialVector<X> _expansion(const Vector<X>& x, const ushort& s) const { \
      const Vector<X>& p=this->_p;                                      \
      SparseDifferentialVector<X> dx(as,as,s);                                \
      SparseDifferentialVector<X> dr(rs,as,s);                                \
      for(uint i=0; i!=as; ++i) { dx[i]=x[i]; }                         \
      for(uint i=0; i!=as; ++i) { dx[i][i+1]=1; }                       \
      f(dr,dx,p);                                                       \
      return dr;                                                        \
    }                                                                   \
  private:                                                              \
    Vector<Float> _p;                                                   \
  };                                                                    



#endif
