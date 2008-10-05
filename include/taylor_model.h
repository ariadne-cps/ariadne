#ifndef ARIADNE_TAYLOR_MODEL_H
#define ARIADNE_TAYLOR_MODEL_H

#include <iosfwd>
#include "numeric.h"
#include "vector.h"
#include "sparse_differential.h"

namespace Ariadne {

template<class X> class Vector;
template<class X> class Matrix;

class FunctionInterface;

class MultiIndex;
template<class X> class SparseDifferential;
template<class X> class SparseDifferentialVector;

class TaylorModel;

TaylorModel recentre(const TaylorModel&, const Vector<Interval>& bx, const Vector<Float>& pt);
TaylorModel truncate(const TaylorModel&, const Vector<Interval>&, uint, uint);

TaylorModel compose(const TaylorModel&, const TaylorModel&);
TaylorModel inverse(const TaylorModel&, const Vector<Float>&);
TaylorModel implicit(const TaylorModel&, const Vector<Float>&);
TaylorModel derivative(const TaylorModel&, uint);




/* \brief A taylor_model with multivalued output, using a den.
 */
class TaylorModel {
  typedef Float R;
  typedef Interval I;
  typedef Vector<R> Point;
  typedef Vector<I> Box;
 public:
  /*! \brief Default constructor constructs a Taylor model of order zero with no arguments and no result variables. */
  TaylorModel();
  /*! \brief The zero Taylor model in \a as variables with size \a rs image, order \a o and smoothness \a s, defined on the whole space with centre at the origin. */
  TaylorModel(uint rs, uint as, ushort o, ushort s);
  
  /*! \brief Construct from a domain, centre, and two derivative expansions, one for the centre and one over the entire domain. */
  TaylorModel(const Box& domain, const Point& centre, 
              const SparseDifferentialVector<I>& centre_derivatives, 
              const SparseDifferentialVector<I>& domain_derivatives);
  
  /*! \brief Construct from a domain, centre, an order and a function. */
  TaylorModel(const Box& domain, const Point& centre,
              ushort order, ushort smoothness,
              const FunctionInterface& function);
  
  
  /*! \brief Equality operator. */
  bool operator==(const TaylorModel& p) const;
  /*! \brief Inequality operator. */
  bool operator!=(const TaylorModel& p) const;
  
  // Data access
  /*! \brief The data used to define the centre of the Taylor model. */
  const SparseDifferentialVector<I>& centre_derivatives() const;
  /*! \brief The bounds on the derivative values over the domain of the Taylor model. */
  const SparseDifferentialVector<I>& domain_derivatives() const;
  
  // Data access
  /*! \brief The order of the Taylor model. */
  ushort order() const;
  /*! \brief The smoothness of the function. */
  ushort smoothness() const;
  /*! \brief The size of the argument. */
  uint argument_size() const;
  /*! \brief The size of the result. */
  uint result_size() const;
  
  /// Resizing
  void resize(uint rs, uint as, ushort d, ushort s);

  /*! \brief The domain of validity of the Taylor model. */
  Box domain() const;
  /*! \brief The centre of the derivative expansion. */
  Vector<R> centre() const;
  /*! \brief The range of values the Taylor model can take. */
  Box range() const;
  
  /*! \brief Evaluate the Taylor model at the point \a x. */
  Vector<I> evaluate(const Vector<I>& x) const;
  Vector<I> evaluate(const Vector<R>& x) const;
  
  /*! \brief Compute the derivate of the map at a point. */
  Matrix<I> jacobian(const Vector<I>& s) const;
  Matrix<I> jacobian(const Vector<R>& s) const;
  
  /*! \brief Truncate to a model of lower order and/or smoothness, possibly on a different domain. */
  TaylorModel truncate(const Box& domain, const Vector<R>& centre, 
                       ushort order, ushort smoothness) const;
  
  /*! \brief The zero Taylor model with result size \a rs and argument size \a as. */
  static TaylorModel zero(uint rs, uint as);
  /*! \brief The unit Taylor model with result size 1 and argument size \a as. */
  static TaylorModel one(uint as);
  /*! \brief The constant Taylor model with result size 1 and argument size \a as. */
  static TaylorModel constant(uint as, const R& c);
  
  /*! \brief Write to an output stream. */
  std::ostream& write(std::ostream& os) const;
  
#ifdef DOXYGEN
  /*! \brief Addition. */
  friend template<class R> TaylorModel operator+(const TaylorModel<R1>&, const TaylorModel<R2>&);
  /*! \brief Subtraction. */
  friend template<class R> TaylorModel operator-(const TaylorModel<R1>&, const TaylorModel<R2>&);
  /*! \brief Multiplication. At least one argument must be scalar-valued. */
  friend template<class R> TaylorModel operator*(const TaylorModel<R1>&, const TaylorModel<R2>&);
  
  /*! \brief Multiplication by a scalar. */
  friend template<class R> TaylorModel operator*(const R1&, const TaylorModel<R2>&);
  /*! \brief Multiplication by a scalar. */
  friend template<class R> TaylorModel operator*(const TaylorModel<R1>&, const R2&);
  /*! \brief Division by a scalar. */
  friend template<class R> TaylorModel operator/(const TaylorModel<R1>&, const R2&);
  
  /*! \brief Composition \f$p\circ q(x)=p(q(x))\f$. */
  friend template<class R> TaylorModel compose(const TaylorModel<R1>&, const TaylorModel<R2>&);
  /*! \brief Power of a scalar Taylor model. */
  friend template<class R> TaylorModel pow(const TaylorModel& p, const unsigned int& n);
  /*! \brief Derivative with respect to variable \a k. */
  friend template<class R> TaylorModel derivative(const TaylorModel&, uint k);
  /*! \brief Truncate within \a r to a Taylor model of order at most \a d, putting the error into terms of order \a s. */
  friend template<class R> TaylorModel truncate(const TaylorModel& p, const Rectangle& bb, uint d, uint s);
#endif
 private:
  static void instantiate();
  array< array<I> > _powers(const Vector<I>&) const;
  void _compute_jacobian() const;
  void _set_argument_size(uint n);
  uint _compute_maximum_component_size() const;
 private:
  friend TaylorModel recentre(const TaylorModel&, const Box& bx, const Vector<R>&);
  friend TaylorModel inverse(const TaylorModel&, const Vector<R>&);
  //template<class R> friend void add(TaylorModel&,const TaylorModel&,const TaylorModel&);
  //template<class R> friend void sub(TaylorModel&,const TaylorModel&,const TaylorModel&);
  //template<class R> friend void mul(TaylorModel&,const TaylorModel&,const TaylorModel&);
  //template<class R> friend void div(TaylorModel&,const TaylorModel&,const TaylorModel&);
  //template<class R> friend void compose(TaylorModel&,const TaylorModel&,const TaylorModel&);
  //template<class R> friend void scale(TaylorModel&,const R&);
 private:
  /* Domain of definition. */
  Vector<I> _domain;
  /* The centre of the derivative expansion. */
  Vector<R> _centre;
  uint _smoothness; 
  SparseDifferentialVector<I> _centre_derivatives;
  SparseDifferentialVector<I> _domain_derivatives;
};


TaylorModel operator*(const TaylorModel& tm, const Interval& ivl);
TaylorModel operator*(const TaylorModel& tm1, const TaylorModel& tm2);

TaylorModel add(const TaylorModel&, const TaylorModel&);
TaylorModel sub(const TaylorModel&, const TaylorModel&);

TaylorModel combine(const TaylorModel&, const TaylorModel&);
TaylorModel join(const TaylorModel&, const TaylorModel&);

std::ostream& operator<<(std::ostream&, const TaylorModel&);

} // namespace Ariadne

#endif /* ARIADNE_TAYLOR_MODEL_H */
