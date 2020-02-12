/***************************************************************************
 *            function/function.hpp
 *
 *  Copyright  2008-20  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

/*! \file function/function.hpp
 *  \brief Built-in and user functions and expressions
 */

#ifndef ARIADNE_FUNCTION_HPP
#define ARIADNE_FUNCTION_HPP

#include "../config.hpp"

#include <cstdarg>
#include <iosfwd>
#include <iostream>

#include "../utility/declarations.hpp"
#include "../utility/macros.hpp"
#include "../utility/pointer.hpp"
#include "../utility/container.hpp"
#include "../utility/metaprogramming.hpp"

#include "../function/function_interface.hpp"

#include "../numeric/numeric.hpp"
#include "../algebra/vector.hpp"
#include "../algebra/covector.hpp"
#include "../algebra/differential.hpp"
#include "../function/domain.hpp"

namespace Ariadne {

template<class P, class D> struct VectorFunctionElementReference;

template<class T> class Variable;
typedef Variable<Real> RealVariable;
template<class T> class Expression;
typedef Expression<Real> RealExpression;
template<class P, class D, class C> class FunctionExpression;

//! \brief Named (static) constructors for constant and coordinate functions.
template<class P>
class FunctionConstructors {
    static_assert(Or<IsSame<P,ApproximateTag>,IsSame<P,ValidatedTag>,IsSame<P,EffectiveTag>>::value,"P must be an information level/paradigm.");
    typedef Number<P> Y;
  public:
    typedef Y NumericType;

    //@{
    //! \name Named (static) constructors
    static ScalarFunction<P,BoxDomainType> zero(SizeType as); //!< \brief %Constant scalar function zero in \a as arguments.
    static ScalarFunction<P,BoxDomainType> constant(SizeType as, NumericType c); //!< \brief %Constant scalar function \a c in \a as arguments.
    static ScalarFunction<P,BoxDomainType> coordinate(SizeType as, SizeType j); //!< \brief %Coordinate function taking \a j -th index in \a as arguments.
    static List<ScalarFunction<P,BoxDomainType>> coordinates(SizeType ns); //!< \brief All coordinates over \a ns variables.
    static VectorFunction<P,BoxDomainType> zeros(SizeType rs, SizeType as); //!< \brief The zero function over \a as variables returning a result of size \a rs.
    static VectorFunction<P,BoxDomainType> constant(SizeType as, Vector<NumericType> c); //!< \brief %Constant vector function \a c in \a as arguments.
    static VectorFunction<P,BoxDomainType> identity(SizeType ns); //!< \brief The identity over \a ns variables.

    static ScalarFunction<P,BoxDomainType> zero(BoxDomainType dom); //!< \brief %Constant scalar function zero over domain \a dom.
    static ScalarFunction<P,BoxDomainType> constant(BoxDomainType dom, NumericType c); //!< \brief %Constant function \a c over domain \a dom.
    static ScalarFunction<P,BoxDomainType> coordinate(BoxDomainType dom, SizeType j); //!< \brief %Coordinate function taking \a j -th index over domain \a dom.
    static List<ScalarFunction<P,BoxDomainType>> coordinates(BoxDomainType dom); //!< \brief %List of all coordinate functions over domain \a dom.
    static VectorFunction<P,BoxDomainType> zeros(SizeType rs, BoxDomainType dom); //!< \brief Zero function over domain \a dom returning a result of size \a rs.
    static VectorFunction<P,BoxDomainType> constant(BoxDomainType dom, Vector<NumericType> c); //!< \brief Constant function over domain \a dom  taking values \a c.
    static VectorFunction<P,BoxDomainType> identity(BoxDomainType dom); //!< \brief Identity function over domain \a dom.

    static ScalarFunction<P,IntervalDomainType> zero(); //!< \brief %Constant scalar univariate function zero.
    static ScalarFunction<P,IntervalDomainType> constant(NumericType c); //!< \brief %Constant scalar univariate function with value \a c.
    static ScalarFunction<P,IntervalDomainType> coordinate(); //!< \brief %Scalar univariate coordinate function.
    static VectorFunction<P,IntervalDomainType> zeros(SizeType rs); //!< \brief %Constant vector univariate function zero returning a result of size \a rs.
    static VectorFunction<P,IntervalDomainType> constant(Vector<NumericType> c); //!< \brief Constant univariate function taking values \a c.
    static ScalarFunction<P,IntervalDomainType> identity(); //!< \brief %Scalar univariate identity function.

    static ScalarFunction<P,IntervalDomainType> zero(IntervalDomainType dom); //!< \brief %Constant scalar function zero over domain \a dom.
    static ScalarFunction<P,IntervalDomainType> constant(IntervalDomainType dom, NumericType c); //!< \brief %Constant scalar function \a c over domain \a dom.
    static ScalarFunction<P,IntervalDomainType> coordinate(IntervalDomainType dom, SizeOne j); //!< \brief %Coordinate function over domain \a dom.
    static ScalarFunction<P,IntervalDomainType> coordinate(IntervalDomainType dom); //!< \brief %Coordinate function over domain \a dom.
    static VectorFunction<P,IntervalDomainType> zeros(SizeType rs, IntervalDomainType dom); //!< \brief %Zero vector function over domain \a dom returning a result of size \a rs.
    static VectorFunction<P,IntervalDomainType> constant(IntervalDomainType dom, Vector<NumericType> c); //!< \brief Constant function over domain \a dom  taking values \a c.
    static ScalarFunction<P,IntervalDomainType> identity(IntervalDomainType dom); //!< \brief %Identity function over domain \a dom.

    //@}

};

template<class P, class X> using EvaluateType = decltype(declval<ScalarMultivariateFunctionInterface<P>>()._evaluate(declval<Vector<X>>()));

//! \brief Function operations which depend on the scalar/vector nature of the arguments and result.
template<class P, class D, class C> class FunctionFacade {
};

//! \brief Function operations for scalar univariate functions.
template<class P> class FunctionFacade<P,IntervalDomainType,IntervalDomainType> {
  public:
    //! \relates Function \brief The slope (derivative) of the function at the point \a x.
    template<class X> Scalar<EvaluateType<P,X>> slope(X const& x) const;
    //! \relates Function \brief The derivative of the function at the point \a x.
    FunctionExpression<P,IntervalDomainType,IntervalDomainType> operator() (const RealVariable& x) const;
};

//! \brief Function operations for vector univariate functions.
template<class P> class FunctionFacade<P,IntervalDomainType,BoxDomainType> {
  public:
    //! \relates Function \brief The tangent vector of derivatives of the function at the point \a x.
    template<class X> Vector<EvaluateType<P,X>> tangent(X const& x) const;
    FunctionExpression<P,IntervalDomainType,BoxDomainType> operator() (const RealVariable& x) const;
};

//! \brief Function operations for scalar multivariate functions.
template<class P> class FunctionFacade<P,BoxDomainType,IntervalDomainType> {
    typedef Number<P> Y;
  public:
    //! \relates Function \brief The gradient covector of derivatives of the function at the point \a x.
    template<class X> Covector<EvaluateType<P,X>> gradient(Vector<X> const& x) const;
    FunctionExpression<P,BoxDomainType,IntervalDomainType> operator() (const Vector<RealVariable>& x) const;
    //FunctionExpression<P,BoxDomainType,IntervalDomainType> operator() (const Vector<RealExpression>& x) const;
};

//! \brief Function operations for vector multivariate functions.
template<class P> class FunctionFacade<P,BoxDomainType,BoxDomainType> {
    typedef Number<P> Y;
  public:
    //! \relates Function \brief The Jacobian matrix of derivatives of the function at the point \a x.
    template<class X> Matrix<EvaluateType<P,X>> jacobian(Vector<X> const& x) const;
    FunctionExpression<P,BoxDomainType,BoxDomainType> operator() (const Vector<RealVariable>& x) const;
};

template<class P, class D, class C> class DeclareFunctionOperations;
template<class P, class D> class DeclareFunctionOperations<P,D,IntervalDomainType>
    : DeclareElementaryAlgebraOperations<Function<P,D,IntervalDomainType>,Number<P>> { };
template<class P, class D> class DeclareFunctionOperations<P,D,BoxDomainType>
    : DeclareVectorAlgebraOperators<Function<P,D,BoxDomainType>,Function<P,D,IntervalDomainType>,Vector<Number<P>>,Number<P>> { };

//! \brief A mixin which provides Scalar/Vector function methods and friends.
template<class P, class D, class C> class DispatchFunctionOperations;
template<class P, class D> class DispatchFunctionOperations<P,D,IntervalDomainType>
    : DispatchElementaryAlgebraOperations<Function<P,D,IntervalDomainType>,Number<P>>
{
};
template<class P, class D> class DispatchFunctionOperations<P,D,BoxDomainType>
    : DeclareVectorAlgebraOperators<Function<P,D,BoxDomainType>,Function<P,D,IntervalDomainType>,Vector<Number<P>>,Number<P>> { };

//! \ingroup FunctionModule
//! \brief A generic function which can be evaluated over the number type \a X,  \f$f:\X^n\rightarrow\X^m\f$.
template<class P, class D, class C>
class Function
    : public FunctionConstructors<P>
    , public FunctionFacade<P,D,C>
//    , public DeclareFunctionOperations<P,D,C>
    , public DispatchFunctionOperations<P,D,C>
{
    static_assert(IsStronger<P,ApproximateTag>::value,"P must be an information level/paradigm.");
    typedef Number<P> Y;
  protected:
    SharedPointer< const FunctionInterface<P,D,C> > _ptr;
  public:
    typedef P InformationTag; //!< The type of information (Effective, Validated or Approximate) provided by the function implementation.
    typedef P Paradigm;
    typedef D DomainType; //!< The type of the domain.
    typedef C CodomainType; //!< The type of the codomain.
    typedef Number<P> NumericType; //!< The numeric type required to construct a constant scalar function.
    typedef typename C::DimensionType ResultSizeType; //!< The type used to descibe the size of an element of the codomain.
    typedef typename D::DimensionType ArgumentSizeType; //!< The type used to descibe the size of an element of the domain.
    typedef typename D::DimensionType ArgumentIndexType; //!< The type used to descibe an index into an element of the domain.

    //! \brief The type of an argument to the function whose scalar type is \a Y.
    template<class Y> using Argument = typename ElementTraits<D>::template Type<Y>;
    //! \brief The type for the result of calling the function whose scalar type is \a Y.
    template<class Y> using Result = typename ElementTraits<C>::template Type<Y>;


    //@{
    //! \name User constructors.

    explicit Function(EuclideanDomain dom); //!< \deprecated
    explicit Function(ResultSizeType rs, EuclideanDomain dom); //!< \deprecated

    //! \brief Construct the zero function with the given domain \a dom and default codomain (Zero-dimensional for a vector function, one-dimensional for a scalar function).
    explicit Function(DomainType dom);
    //! \brief Construct the zero function with the given domain \a dom and number of result variables \a rs.
    explicit Function(ResultSizeType rs, DomainType dom);

    //! \brief Construct a function with the given domain \a dom and number of result variables \a rs.
    explicit Function(DomainType dom, Result<Formula<Y>>const& e);
    //! \brief Construct a function each of whose components are \a sf with the given domain \a dom.
    explicit Function(ResultSizeType rs, ScalarFunction<P,D> sf);

    //! \brief Create a vector function from an initializer list of scalar functions.
    Function(InitializerList<ScalarFunction<P,D>> const& lsf);
    //! \brief Create a vector function from a list of scalar functions.
    Function(List<ScalarFunction<P,D>> const& lsf);
    //! \brief Create a vector function from a vector of scalar functions.
    Function(Vector<ScalarFunction<P,D>> const& lsf);
    //@}

    //@{
    //! \name Prototype constructors.

    //! \brief Construct a zero scalar function with the same domain.
    ScalarFunction<P,D> create_zero() const { return ScalarFunction<P,D>::zero(this->domain()); }
    //! \brief Construct a scalar constant function with value \a c and the same domain.
    ScalarFunction<P,D> create_constant(NumericType c) const { return ScalarFunction<P,D>::constant(this->domain(),c); }
    //! \brief Construct the scalar coordinate function for the index \a j.
    ScalarFunction<P,D> create_coordinate(ArgumentIndexType j) const { return ScalarFunction<P,D>::coordinate(this->domain(),j); }
    //! \brief Construct the vector constant function with values \a c.
    VectorFunction<P,D> create_constant(Vector<NumericType> c) const { return VectorFunction<P,D>::constant(this->domain(),c); }
    //@}

    //@{
    //! \name Handle-interface methods.
    Function(); //!< \brief Create an invalid (null) function.
    explicit Function(FunctionInterface<P,D,C>* p) : _ptr(p) { } //!< \brief Capture a newly-allocated function pointer.
    explicit Function(SharedPointer<FunctionInterface<P,D,C>> p) : _ptr(p) { } //!< \brief Construct from a managed pointer.
    Function(const FunctionInterface<P,D,C>& t) : _ptr(t._clone()) { }  //!< \brief Clone from a reference.

    //! \brief Assign from a reference.
    Function<P,D,C>& operator=(const FunctionInterface<P,D,C>& f) {
        _ptr=std::shared_ptr< FunctionInterface<P,D,C> >(f._clone()); return *this; }

    //! \brief Return a managed (shared) pointer to the underlying interface class.
    SharedPointer< const FunctionInterface<P,D,C> > managed_pointer() const  { return _ptr; }
    //! \brief Return a raw (builtin) pointer to the underlying interface class.
    const FunctionInterface<P,D,C>* raw_pointer() const  { return _ptr.operator->(); }
    //! \brief Return a reference to the underlying interface class.
    const FunctionInterface<P,D,C>& reference() const  { return _ptr.operator*(); }
    //! \brief Convert to a reference to the underlying interface class.
    operator const FunctionInterface<P,D,C>& () const { return _ptr.operator*(); }

    //@}

    //@{
    //! \name Conversions and assignment.

    //! \brief Convert from a function class specifying more information.
    template<class PP, EnableIf<IsStronger<PP,P>> =dummy>
    Function(const Function<PP,D,C>& f)
        : _ptr(std::dynamic_pointer_cast< const FunctionInterface<P,D,C> >(f.managed_pointer())) { }
    //! \brief Assign from a function class specifying more information.
    template<class PP, EnableIf<IsStronger<PP,P>> =dummy>
        Function<P,D,C>& operator=(Result<NumericType> const& c); // { return *this=this->create_constant(c); }

    //! \brief Set equal to the constant value \a c.
    Function<P,D,C>& operator=(const Result<NumericType>& c) {
        return (*this)=this->create_constant(c); }

    //@}

    //@{
    //! \name Query domain and codomain.

    //! \brief The domain of the function.
    DomainType domain() const {
        return this->reference().domain(); }
    //! \brief The codomain of the function.
    CodomainType codomain() const {
        return this->reference().codomain(); }
    //! \brief The number of scalar variables in the argument.
    ArgumentSizeType argument_size() const {
        return this->reference().argument_size(); }
    //! \brief The number of scalar variables in the result.
    ResultSizeType result_size() const {
        return this->reference().result_size(); }
    //@}

    //@{
    //! \name Call/evaluate

    //! \brief Call the function on an argument of concrete scalar type \a X.
    template<class X> auto operator() (const Argument<X>& x) const -> decltype(this->reference()._evaluate(x)) {
        return this->reference()._evaluate(x); }
#ifdef DOXYGEN
    //! \brief Call the function on an argument of concrete scalar type \a X.
    friend template<class X> auto evaluate (const Function<P,D,C>& f,const Argument<X>& x) const -> decltype(f.reference()._evaluate(x)) {
        return f(x); }
#endif
//! \brief Call the function on an argument of concrete scalar type \a X. \deprecated
    template<class X> auto evaluate(const Argument<X>& x) const -> decltype(this->reference()._evaluate(x)) {
        return this->reference()._evaluate(x); }
    //@}

    friend VectorFunction<P,D> operator*(ScalarFunction<P,D> const&, Vector<Y> const&);

    //@{
    //! \name Differential function operations.

    //! \brief The derivative of the function with respect to the \a k -th variable.
    Function<P,D,C> derivative(ElementIndexType<D> k) const {
        return Function<P,D,C>(this->reference()._derivative(k)); }
    //! \brief The derivative of the function \a f with respect to the \a k -th variable.
    friend Function<P,D,C> derivative(Function<P,D,C> const& f, ElementIndexType<D> k) {
        return f.derivative(k); }

    //! \brief The differential (partial derivatives) of the function at the point \a x, computed to degree \a d.
    template<class X> decltype(auto) differential(const Argument<X>& x, DegreeType d) const {
        return this->_ptr->_evaluate(Differential<EvaluateType<P,X>>::identity(d,x)); }
    //@}

    //@{
    //! \name Vector (indexing) operations.

    //! \brief Set the \a i -th component of a vector function to \a f.
    Void set(SizeType i, ScalarFunction<P,D> f);
    //! \brief Get the \a i -th component of a vector function.
    Function<P,D,IntervalDomainType> get(SizeType i) const;
    //! \brief Get the \a i -th component of a vector function.
    const Function<P,D,IntervalDomainType> operator[](SizeType i) const;
    //! \brief Get the components of a specified by \a rng function.
    const Function<P,D,BoxDomainType> operator[](Range rng) const;
    //! \brief A reference to the \a i -th component of a vector function.
    VectorFunctionElementReference<P,D> operator[](SizeType i);

    //@}

    //@{
    //! \name Input/output operations.

    //! \brief Write to an output stream.
    friend OutputStream& operator<<(OutputStream& os, Function<P,D,C> const& f) { f._ptr->_write(os); return os; }

    //@}
};

template<class P, class D> struct AlgebraOperations<ScalarFunction<P,D>,Number<P>>
{
    using F=ScalarFunction<P,D>;
    using C=Number<P>;
    static ScalarFunction<P,D> apply(BinaryElementaryOperator op, ScalarFunction<P,D> const& f1, ScalarFunction<P,D> const& f2);
    static ScalarFunction<P,D> apply(UnaryElementaryOperator op, ScalarFunction<P,D> const& f);
    static ScalarFunction<P,D> apply(BinaryElementaryOperator op, ScalarFunction<P,D> const& f1, Number<P> const& c2);
    static ScalarFunction<P,D> apply(BinaryElementaryOperator op, Number<P> const& c1, ScalarFunction<P,D> const& f2);
    static ScalarFunction<P,D> apply(GradedElementaryOperator op, ScalarFunction<P,D> const& f, Int n);
};

template<class P, class D, class C> inline OutputStream&
operator<<(OutputStream& os, const Function<P,D,C>& f) {
    return f._write(os); }

template<class P, class C, class X> inline decltype(auto)
evaluate(const Function<P,IntervalDomainType,C>& f, const Scalar<X>& x) {
    return f(x); }

template<class P, class C, class X> inline decltype(auto)
evaluate(const Function<P,BoxDomainType,C>& f, const Vector<X>& x) {
    return f(x); }


template<class P, class C, class X> inline decltype(auto)
differential(const Function<P,IntervalDomainType,C>& f, const X& x, DegreeType d) {
    return f.differential(x,d); }

template<class P, class C, class X> inline decltype(auto)
differential(const Function<P,BoxDomainType,C>& f, const Vector<X>& x, DegreeType d) {
    return f.differential(x,d); }

template<class P, class D, class C, EnableIf<IsSame<ElementSizeType<D>,SizeOne>> =dummy> inline
Function<P,D,C> derivative(Function<P,D,C> const& f) {
    return f.derivative(SizeOne()); }

template<class P, class X> Scalar<EvaluateType<P,X>>
slope(const ScalarUnivariateFunction<P>& f, const Scalar<X>& x) {
    return differential(f,x,1u).gradient()[0]; }

template<class P, class X> Vector<EvaluateType<P,X>>
tangent(const VectorUnivariateFunction<P>& f, const Scalar<X>& x) {
    return column(differential(f,x,1u).jacobian(),0u); }

template<class P, class X> Covector<EvaluateType<P,X>>
gradient(const ScalarMultivariateFunction<P>& f, const Vector<X>& x) {
    return differential(f,x,1u).gradient(); }

template<class P, class X> Matrix<EvaluateType<P,X>>
jacobian(const VectorMultivariateFunction<P>& f, const Vector<X>& x) {
    return differential(f,x,1u).jacobian(); }

template<class P> template<class X> EvaluateType<P,X>
FunctionFacade<P,IntervalDomainType,IntervalDomainType>::slope(X const& x) const {
    return Ariadne::slope(static_cast<Function<P,IntervalDomainType,IntervalDomainType>const&>(*this),x);
}

template<class P> template<class X> Vector<EvaluateType<P,X>>
FunctionFacade<P,IntervalDomainType,BoxDomainType>::tangent(X const& x) const {
    return Ariadne::tangent(static_cast<Function<P,IntervalDomainType,BoxDomainType>const&>(*this),x);
}

template<class P> template<class X> Covector<EvaluateType<P,X>>
FunctionFacade<P,BoxDomainType,IntervalDomainType>::gradient(Vector<X> const& x) const {
    return Ariadne::gradient(static_cast<Function<P,BoxDomainType,IntervalDomainType>const&>(*this),x);
}

template<class P> template<class X> Matrix<EvaluateType<P,X>>
FunctionFacade<P,BoxDomainType,BoxDomainType>::jacobian(Vector<X> const& x) const {
    return Ariadne::jacobian(static_cast<Function<P,BoxDomainType,BoxDomainType>const&>(*this),x);
}



EffectiveScalarUnivariateFunction compose(const EffectiveScalarUnivariateFunction& f, const EffectiveScalarUnivariateFunction& g);
EffectiveScalarUnivariateFunction compose(const EffectiveScalarMultivariateFunction& f, const EffectiveVectorUnivariateFunction& g);
EffectiveVectorUnivariateFunction compose(const EffectiveVectorUnivariateFunction& f, const EffectiveScalarUnivariateFunction& g);
EffectiveVectorUnivariateFunction compose(const EffectiveVectorMultivariateFunction& f, const EffectiveVectorUnivariateFunction& g);
EffectiveScalarMultivariateFunction compose(const EffectiveScalarUnivariateFunction& f, const EffectiveScalarMultivariateFunction& g);
EffectiveScalarMultivariateFunction compose(const EffectiveScalarMultivariateFunction& f, const EffectiveVectorMultivariateFunction& g);
EffectiveVectorMultivariateFunction compose(const EffectiveVectorUnivariateFunction& f, const EffectiveScalarMultivariateFunction& g);
EffectiveVectorMultivariateFunction compose(const EffectiveVectorMultivariateFunction& f, const EffectiveVectorMultivariateFunction& g);

ValidatedScalarUnivariateFunction compose(const ValidatedScalarUnivariateFunction& f, const ValidatedScalarUnivariateFunction& g);
ValidatedScalarUnivariateFunction compose(const ValidatedScalarMultivariateFunction& f, const ValidatedVectorUnivariateFunction& g);
ValidatedVectorUnivariateFunction compose(const ValidatedVectorUnivariateFunction& f, const ValidatedScalarUnivariateFunction& g);
ValidatedVectorUnivariateFunction compose(const ValidatedVectorMultivariateFunction& f, const ValidatedVectorUnivariateFunction& g);
ValidatedScalarMultivariateFunction compose(const ValidatedScalarUnivariateFunction& f, const ValidatedScalarMultivariateFunction& g);
ValidatedScalarMultivariateFunction compose(const ValidatedScalarMultivariateFunction& f, const ValidatedVectorMultivariateFunction& g);
ValidatedVectorMultivariateFunction compose(const ValidatedVectorUnivariateFunction& f, const ValidatedScalarMultivariateFunction& g);
ValidatedVectorMultivariateFunction compose(const ValidatedVectorMultivariateFunction& f, const ValidatedVectorMultivariateFunction& g);

ApproximateScalarUnivariateFunction compose(const ApproximateScalarUnivariateFunction& f, const ApproximateScalarUnivariateFunction& g);
ApproximateScalarUnivariateFunction compose(const ApproximateScalarMultivariateFunction& f, const ApproximateVectorUnivariateFunction& g);
ApproximateVectorUnivariateFunction compose(const ApproximateVectorUnivariateFunction& f, const ApproximateScalarUnivariateFunction& g);
ApproximateVectorUnivariateFunction compose(const ApproximateVectorMultivariateFunction& f, const ApproximateVectorUnivariateFunction& g);
ApproximateScalarMultivariateFunction compose(const ApproximateScalarUnivariateFunction& f, const ApproximateScalarMultivariateFunction& g);
ApproximateScalarMultivariateFunction compose(const ApproximateScalarMultivariateFunction& f, const ApproximateVectorMultivariateFunction& g);
ApproximateVectorMultivariateFunction compose(const ApproximateVectorUnivariateFunction& f, const ApproximateScalarMultivariateFunction& g);
ApproximateVectorMultivariateFunction compose(const ApproximateVectorMultivariateFunction& f, const ApproximateVectorMultivariateFunction& g);


ValidatedScalarMultivariateFunction& operator*=(ValidatedScalarMultivariateFunction& sf, const ExactNumber& c);
EffectiveVectorMultivariateFunction operator*(const EffectiveNumber& c, const EffectiveVectorMultivariateFunction& vf);

EffectiveScalarMultivariateFunction embed(SizeType as1, const EffectiveScalarMultivariateFunction& f2, SizeType as3);
EffectiveVectorMultivariateFunction embed(SizeType as1, const EffectiveVectorMultivariateFunction& f2, SizeType as3);

EffectiveVectorMultivariateFunction join(const EffectiveScalarMultivariateFunction& f1, const EffectiveScalarMultivariateFunction& f2);
EffectiveVectorMultivariateFunction join(const EffectiveScalarMultivariateFunction& f1, const EffectiveVectorMultivariateFunction& f2);
EffectiveVectorMultivariateFunction join(const EffectiveVectorMultivariateFunction& f1, const EffectiveScalarMultivariateFunction& f2);
EffectiveVectorMultivariateFunction join(const EffectiveVectorMultivariateFunction& f1, const EffectiveVectorMultivariateFunction& f2);

EffectiveScalarMultivariateFunction lie_derivative(const EffectiveScalarMultivariateFunction& g, const EffectiveVectorMultivariateFunction& f);
EffectiveVectorMultivariateFunction lie_derivative(const EffectiveVectorMultivariateFunction& g, const EffectiveVectorMultivariateFunction& f);

Formula<Real> make_formula(const EffectiveScalarMultivariateFunction& f);
Vector<Formula<Real>> make_formula(const EffectiveVectorMultivariateFunction& f);
//RealExpression evaluate(EffectiveScalarMultivariateFunction const& f, Vector<RealVariable> const& vars);



ValidatedVectorMultivariateFunction join(const ValidatedScalarMultivariateFunction& f1, const ValidatedScalarMultivariateFunction& f2);
ValidatedVectorMultivariateFunction join(const ValidatedScalarMultivariateFunction& f1, const ValidatedVectorMultivariateFunction& f2);
ValidatedVectorMultivariateFunction join(const ValidatedVectorMultivariateFunction& f1, const ValidatedScalarMultivariateFunction& f2);
ValidatedVectorMultivariateFunction join(const ValidatedVectorMultivariateFunction& f1, const ValidatedVectorMultivariateFunction& f2);


ApproximateVectorMultivariateFunction join(const ApproximateScalarMultivariateFunction& f1, const ApproximateScalarMultivariateFunction& f2);
ApproximateVectorMultivariateFunction join(const ApproximateScalarMultivariateFunction& f1, const ApproximateVectorMultivariateFunction& f2);
ApproximateVectorMultivariateFunction join(const ApproximateVectorMultivariateFunction& f1, const ApproximateScalarMultivariateFunction& f2);
ApproximateVectorMultivariateFunction join(const ApproximateVectorMultivariateFunction& f1, const ApproximateVectorMultivariateFunction& f2);


//! \brief A reference into an element of a vector function.
template<class P, class D>
struct VectorFunctionElementReference
    : DispatchFunctionOperations<P,D,IntervalDomainType>
{
    template<class Y> using Argument = typename ElementTraits<D>::template Type<Y>;
    typedef IntervalDomainType SC; typedef BoxDomainType VC;
    typedef VectorFunctionElementReference<P,D> SelfType;
    typedef Function<P,D,SC> ElementType;
    typedef typename ElementType::NumericType NumericType;
    Function<P,D,VC>& _vf; SizeType _i;
    VectorFunctionElementReference<P,D>(Function<P,D,VC>& vf, SizeType i) : _vf(vf), _i(i) { }
    template<class WP> operator Function<WP,D,SC> () const;
    Void operator=(const Function<P,D,SC>& sf);
    VectorFunctionElementReference<P,D>& operator=(const VectorFunctionElementReference<P,D>& sfr);
    D domain() const;
    template<class XX> XX evaluate(const Argument<XX> & x) const;
    template<class XX> XX operator()(const Argument<XX>& x) const;
};

template<class P, class D> template<class WP> inline
VectorFunctionElementReference<P,D>::operator Function<WP,D,SC> () const {
    return this->_vf.get(this->_i);
}

template<class P, class D> inline VectorFunctionElementReference<P,D> make_element_reference(ScalarFunction<P,D>& sf, SizeType i) {
    throw std::runtime_error(""); }
template<class P, class D> inline VectorFunctionElementReference<P,D> make_element_reference(VectorFunction<P,D>& vf, SizeType i) {
    return VectorFunctionElementReference<P,D>(vf,i); }

template<class P, class D, class C> inline VectorFunctionElementReference<P,D> Function<P,D,C>::operator[](SizeType i) {
    return make_element_reference(*this,i); }

template<class P, class D, class C> inline ScalarFunction<P,D> const Function<P,D,C>::operator[](SizeType i) const {
    return this->get(i); }

template<class P, class D> inline OutputStream& operator<<(OutputStream& os, const VectorFunctionElementReference<P,D>& vfe) {
    return  os << static_cast< Function<P,D,IntervalDomainType> >(vfe); }

template<class P, class D> inline Void VectorFunctionElementReference<P,D>::operator=(const Function<P,D,SC>& sf) {
    _vf.set(_i,sf); }
template<class P, class D> inline VectorFunctionElementReference<P,D>& VectorFunctionElementReference<P,D>::operator=(const VectorFunctionElementReference<P,D>& sfr) {
    _vf.set(_i,static_cast<Function<P,D,SC>>(sfr)); return *this; }
template<class P, class D> inline D VectorFunctionElementReference<P,D>::domain() const {
    return _vf.domain(); }
template<class P, class D> template<class XX> inline XX VectorFunctionElementReference<P,D>::evaluate(const Argument<XX> & x) const {
    return static_cast<Function<P,D,SC>>(*this).evaluate(x); }
template<class P, class D> template<class XX> inline XX VectorFunctionElementReference<P,D>::operator()(const Argument<XX> & x) const {
    return static_cast<Function<P,D,SC>>(*this).evaluate(x); }



UpperIntervalType evaluate_range(ScalarMultivariateFunction<ValidatedTag>const& f, const Vector<UpperIntervalType>& x);
Vector<UpperIntervalType> evaluate_range(VectorMultivariateFunction<ValidatedTag>const& f, const Vector<UpperIntervalType>& x);
Vector<Differential<UpperIntervalType>> derivative_range(VectorMultivariateFunction<ValidatedTag>const& f, const Vector<Differential<UpperIntervalType>>& x);
Covector<UpperIntervalType> gradient_range(ValidatedScalarMultivariateFunction const& f, const Vector<UpperIntervalType>& x);
Matrix<UpperIntervalType> jacobian_range(ValidatedVectorMultivariateFunction const& f, const Vector<UpperIntervalType>& x);



/*
inline Matrix<UpperIntervalType> jacobian(VectorMultivariateFunction<ValidatedTag>const& f, const Vector<UpperIntervalType>& x) {
    return static_cast<Matrix<UpperIntervalType>>(f.jacobian(reinterpret_cast<Vector<ValidatedNumericType>const&>(x))); }
inline Matrix<UpperIntervalType> jacobian(VectorMultivariateFunction<ValidatedTag>const& f, const Vector<ExactInterval>& x) {
    return static_cast<Matrix<UpperIntervalType>>(f.jacobian(reinterpret_cast<Vector<ValidatedNumericType>const&>(x))); }
inline Matrix<UpperIntervalType> jacobian_range(VectorMultivariateFunction<ValidatedTag>const& f, const Vector<UpperIntervalType>& x) {
    return static_cast<Matrix<UpperIntervalType>>(f.jacobian(reinterpret_cast<Vector<ValidatedNumericType>const&>(x))); }

// FIXME: Needed to override templated gradient and jacobian
inline Covector<UpperIntervalType> gradient(ScalarMultivariateFunction<EffectiveTag>const& f, const Vector<UpperIntervalType>& x) {
    return static_cast<Covector<UpperIntervalType>>(gradient(f,reinterpret_cast<Vector<ValidatedNumericType>const&>(x))); }
inline Covector<UpperIntervalType> gradient(ScalarMultivariateFunction<EffectiveTag>const& f, const Vector<ExactInterval>& x) {
    return gradient(f,static_cast<Vector<UpperIntervalType>>(x)); }
inline Matrix<UpperIntervalType> jacobian(VectorMultivariateFunction<EffectiveTag>const& f, const Vector<UpperIntervalType>& x) {
    return static_cast<Matrix<UpperIntervalType>>(f.jacobian(reinterpret_cast<Vector<ValidatedNumericType>const&>(x))); }
inline Matrix<UpperIntervalType> jacobian(VectorMultivariateFunction<EffectiveTag>const& f, const Vector<ExactInterval>& x) {
    return jacobian(f,static_cast<Vector<UpperIntervalType>>(x)); }
*/

//! \brief A class for constructing generic functions with a specific implementation.
//! \sa Function
template<class P> class FunctionFactory {
    SharedPointer< const FunctionFactoryInterface<P> > _ptr;
  public:
    FunctionFactory(const FunctionFactoryInterface<P>& ref) : _ptr(ref.clone()) { }
    FunctionFactory(const FunctionFactoryInterface<P>* ptr) : _ptr(ptr) { }
    FunctionFactory(SharedPointer< const FunctionFactoryInterface<P> > ptr) : _ptr(ptr) { }
    //! \brief Create a scalar function equal to \a f over domain \a dom.
    inline ScalarMultivariateFunction<P> create(const BoxDomainType& dom, const ScalarMultivariateFunctionInterface<P>& f) const;
    //! \brief Create a vector function equal to \a f over domain \a dom.
    inline VectorMultivariateFunction<P> create(const BoxDomainType& dom, const VectorMultivariateFunctionInterface<P>& f) const;
    //! \brief Create a scalar zero function over domain \a dom.
    inline ScalarMultivariateFunction<P> create_zero(const BoxDomainType& dom) const;
    //! \brief Create the identity function over domain \a dom.
    inline VectorMultivariateFunction<P> create_identity(const BoxDomainType& dom) const;
    //! \brief Write to an output stream.
    friend OutputStream& operator<<(OutputStream& os, const FunctionFactory<P>& factory);
};

typedef FunctionFactory<ValidatedTag> ValidatedFunctionFactory;

template<>
class FunctionFactory<ValidatedTag>
{
    SharedPointer< const FunctionFactoryInterface<ValidatedTag> > _ptr;
  public:
    FunctionFactory(const FunctionFactoryInterface<ValidatedTag>& ref) : _ptr(ref.clone()) { }
    FunctionFactory(const FunctionFactoryInterface<ValidatedTag>* ptr) : _ptr(ptr) { }
    FunctionFactory(SharedPointer< const FunctionFactoryInterface<ValidatedTag> > ptr) : _ptr(ptr) { }
    inline ValidatedScalarMultivariateFunction create(const BoxDomainType& d, const ValidatedScalarMultivariateFunctionInterface& f) const;
    inline ValidatedVectorMultivariateFunction create(const BoxDomainType& d, const ValidatedVectorMultivariateFunctionInterface& f) const;
    inline ValidatedScalarMultivariateFunction create_zero(const BoxDomainType& d) const;
    inline ValidatedVectorMultivariateFunction create_identity(const BoxDomainType& d) const;
    friend OutputStream& operator<<(OutputStream& os, const FunctionFactory<ValidatedTag>& factory);
};

inline ValidatedScalarMultivariateFunction FunctionFactoryInterface<ValidatedTag>::create(const BoxDomainType& domain, const ValidatedScalarMultivariateFunctionInterface& function) const {
    return ValidatedScalarMultivariateFunction(SharedPointer<ValidatedScalarMultivariateFunctionInterface>(this->_create(domain,function))); }
inline ValidatedVectorMultivariateFunction FunctionFactoryInterface<ValidatedTag>::create(const BoxDomainType& domain, const ValidatedVectorMultivariateFunctionInterface& function) const {
    return ValidatedVectorMultivariateFunction(SharedPointer<ValidatedVectorMultivariateFunctionInterface>(this->_create(domain,function))); }
inline ValidatedScalarMultivariateFunction FunctionFactoryInterface<ValidatedTag>::create_zero(const BoxDomainType& domain) const {
    return this->create(domain,EffectiveScalarMultivariateFunction::zero(domain.dimension())); }
inline ValidatedVectorMultivariateFunction FunctionFactoryInterface<ValidatedTag>::create_identity(const BoxDomainType& domain) const {
    return this->create(domain,EffectiveVectorMultivariateFunction::identity(domain.dimension())); }

inline ValidatedScalarMultivariateFunction FunctionFactory<ValidatedTag>::create(const BoxDomainType& domain, const ValidatedScalarMultivariateFunctionInterface& function) const {
    return this->_ptr->create(domain,function); }
inline ValidatedVectorMultivariateFunction FunctionFactory<ValidatedTag>::create(const BoxDomainType& domain, const ValidatedVectorMultivariateFunctionInterface& function) const {
    return this->_ptr->create(domain,function); }
inline ValidatedScalarMultivariateFunction FunctionFactory<ValidatedTag>::create_zero(const BoxDomainType& domain) const {
    return this->_ptr->create_zero(domain); }
inline ValidatedVectorMultivariateFunction FunctionFactory<ValidatedTag>::create_identity(const BoxDomainType& domain) const {
    return this->_ptr->create_identity(domain); }

inline OutputStream& operator<<(OutputStream& os, const ValidatedFunctionFactory& factory) {
    factory._ptr->_write(os); return os; }

} // namespace Ariadne

#endif
