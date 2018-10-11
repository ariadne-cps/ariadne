/***************************************************************************
 *            scaled_function_patch.hpp
 *
 *  Copyright 2008-17  Pieter Collins
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

/*! \file scaled_function_patch.hpp
 *  \brief Over-approximations of functions on box domains.
 */

#ifndef ARIADNE_SCALED_FUNCTION_PATCH_HPP
#define ARIADNE_SCALED_FUNCTION_PATCH_HPP

#include <iosfwd>
#include "../utility/container.hpp"
#include "../utility/exceptions.hpp"
#include "../utility/declarations.hpp"
#include "../numeric/numeric.hpp"
#include "../algebra/vector.hpp"
#include "../function/taylor_model.hpp"

#include "../algebra/operations.hpp"
#include "../function/function_interface.hpp"
#include "../function/function_mixin.hpp"
#include "../function/function_model.hpp"
#include "../function/function_model_mixin.hpp"

// FIXME: Added to prevent compilation error in Clang++-5.0. Should not be necessary.
#include "../function/formula.hpp"

namespace Ariadne {

template<class X> class Polynomial;

template<class T> using NumericType = typename T::NumericType;
template<class T> using FunctionType = typename T::FunctionType;
template<class T> using GenericType = typename T::GenericType;

template<class M> using ScalarFunctionType = typename M::ScalarFunctionType;
template<class M> using VectorFunctionType = typename M::VectorFunctionType;
template<class M> using ScalarFunctionModelType = ScalarFunctionModel<Paradigm<M>,BoxDomainType,PrecisionType<M>,ErrorPrecisionType<M>>;
template<class M> using VectorFunctionModelType = VectorFunctionModel<Paradigm<M>,BoxDomainType,PrecisionType<M>,ErrorPrecisionType<M>>;

template<class M> class ScaledFunctionPatch;
template<class M> using ScalarScaledFunctionPatch = ScaledFunctionPatch<M>;
template<class M> class VectorScaledFunctionPatch;
template<class M> class VectorScaledFunctionPatchElementReference;

inline FloatDPApproximation convert_error_to_bounds(const PositiveFloatDPApproximation& e) { return FloatDPApproximation(0.0); }
inline FloatDPBounds convert_error_to_bounds(const PositiveFloatDPUpperBound& e) { return FloatDPBounds(-e.raw(),+e.raw()); }
inline FloatDPBounds convert_error_to_bounds(const FloatDPError& e) { return FloatDPBounds(-e.raw(),+e.raw()); }

inline FloatMPApproximation convert_error_to_bounds(const PositiveFloatMPApproximation& e) { return FloatMPApproximation(0.0,e.precision()); }
inline FloatMPBounds convert_error_to_bounds(const PositiveFloatMPUpperBound& e) { return FloatMPBounds(-e.raw(),+e.raw()); }
inline FloatMPBounds convert_error_to_bounds(const FloatMPError& e) { return FloatMPBounds(-e.raw(),+e.raw()); }

/*
template<class X> X operator+(const X& x1, const GenericType<X>& y2) { return x1+x1.create(y2); }
template<class X> X operator-(const X& x1, const GenericType<X>& y2) { return x1-x1.create(y2); }
template<class X> X operator*(const X& x1, const GenericType<X>& y2) { return x1*x1.create(y2); }
template<class X> X operator/(const X& x1, const GenericType<X>& y2) { return x1/x1.create(y2); }
template<class X> X operator+(const GenericType<X>& y1, const X& x2) { return x2.create(y1)+x2; }
template<class X> X operator-(const GenericType<X>& y1, const X& x2) { return x2.create(y1)-x2; }
template<class X> X operator*(const GenericType<X>& y1, const X& x2) { return x2.create(y1)*x2; }
template<class X> X operator/(const GenericType<X>& y1, const X& x2) { return x2.create(y1)/x2; }
*/

// FIXME: Try to abstract away this template
template<class R, class F, class X> class CanEvaluate {
    template<class RR, class FF, class XX, class=decltype(declval<RR>()=evaluate(declval<FF>(),declval<XX>()))> static True test(int);
    template<class RR, class FF, class XX> static False test(...);
  public:
    static const bool value = decltype(test<R,F,X>(1))::value;
};
template<class R, class F, class X> class CanCall {
    template<class RR, class FF, class XX, class=decltype(declval<RR>()=declval<FF>()(declval<XX>()))> static True test(int);
    template<class RR, class FF, class XX> static False test(...);
  public:
    static const bool value = decltype(test<R,F,X>(1))::value;
};

template<class M> class ScaledFunctionPatchFactory;
template<class M> class ScaledFunctionPatchCreator;

template<class M> struct AlgebraOperations<ScaledFunctionPatch<M>> {
    typedef typename M::NumericType NumericType;
  public:
    template<class OP> static ScaledFunctionPatch<M> apply(OP op, ScaledFunctionPatch<M> const& f);
    template<class OP> static ScaledFunctionPatch<M> apply(OP op, ScaledFunctionPatch<M> const& f1, ScaledFunctionPatch<M> const& f2);
    template<class OP> static ScaledFunctionPatch<M> apply(OP op, ScaledFunctionPatch<M> const& f1, NumericType const& c2);
    template<class OP> static ScaledFunctionPatch<M> apply(OP op, NumericType const& c1, ScaledFunctionPatch<M> const& f2);
    static ScaledFunctionPatch<M> apply(Pow op, ScaledFunctionPatch<M> const& f, Int n);
};

/*! \ingroup FunctionModelSubModule
 *  \brief A ValidatedScalarMultivariateTaylorFunctionModelDP is a type of FunctionModel in which a the restriction of a scalar function \f$f:\R^n\rightarrow\R\f$ on a domain \f$D\f$ is approximated by polynomial \f$p\f$ with uniform error \f$e\f$.
 *
 * Formally, a ValidatedScalarMultivariateTaylorFunctionModelDP is a triple \f$(D,p,e)\f$ representing a set of continuous functions \f$\mathrm{T}(D,p,e)\f$ by
 * \f[ \mathrm{T}(D,p,e) = \{ f:\R^n\rightarrow \R \mid \sup_{x\in D}|f(x)-p(x)| \leq e \} . \f]
 * Note that there is no need for the functions \f$f\f$ to be themselves polynomial, and that no information is given
 * about the values of \f$f\f$ outside of \f$D\f$. Information about the derivatives of \f$f\f$ is also unavailable.
 * However, integrals of \f$f\f$ can be computed.
 *
 * Internally, the polynomial \f$p\f$ is represented as the composition \f$p=m\circ s^{-1}\f$,
 * where \f$m:[-1,+1]^n\rightarrow\R\f$ and \f$s:[-1,+1]^n\rightarrow D\f$ is a scaling function,
 * \f$s_i(y_i)=(a_i+b_i)/2+(b_i-a_i)y_i/2\f$ where \f$D_i=[a_i,b_i]\f$ is the \f$i^\textrm{th}\f$ subinterval of \f$D\f$.
 *
 * When solving algebraic equations by iterative Newton-like methods, it is necessary to compute the derivatives of \f$f\f$.
 * For these applications, it suffices to compute the derivative of \f$p\f$, since only a uniform approximation to the solution is required.
 *
 * Finding exact bounds for the range of \f$p\f$ over \f$D\f$ is an NP-complete problem,
 * for but there are a number of techniques available.
 *
 * \sa Expansion, TaylorModel, ValidatedVectorMultivariateTaylorFunctionModelDP, TaylorConstrainedImageSet.
 */
template<class M> class ScaledFunctionPatch
    : public ScalarMultivariateFunctionModelMixin<ScaledFunctionPatch<M>, typename M::Paradigm, BoxDomainType, typename M::PrecisionType, typename M::ErrorPrecisionType>
    , public DispatchTranscendentalAlgebraOperations<ScaledFunctionPatch<M>, NumericType<M>>
    , public ProvideConcreteGenericArithmeticOperators<ScaledFunctionPatch<M>, ScalarMultivariateFunction<typename M::Paradigm>>
    , public DispatchConcreteGenericAlgebraNumberOperations<ScaledFunctionPatch<M>,NumericType<M>,Number<typename M::Paradigm>>
{
    typedef BoxDomainType D;
    typedef IntervalDomainType C;
    typedef typename M::Paradigm P;
    typedef typename M::RawFloatType F;
    typedef typename M::PrecisionType PR;
    typedef typename M::ErrorPrecisionType PRE;
  public:
    typedef D DomainType;
    typedef M ModelType;
    typedef typename ModelType::CodomainType CodomainType;
    typedef typename ModelType::RangeType RangeType;
    typedef typename ModelType::ExpansionType ExpansionType;
    typedef typename ModelType::CoefficientType CoefficientType;
    typedef typename ModelType::ErrorType ErrorType;
    typedef typename ModelType::NumericType NumericType;
    typedef typename ModelType::Paradigm Paradigm;
    typedef typename ModelType::PrecisionType PrecisionType;
    typedef typename ModelType::NormType NormType;
    typedef ScaledFunctionPatch<M> ScaledFunctionPatchType;
    typedef ScalarMultivariateFunction<Paradigm> FunctionType;
    typedef ScalarMultivariateFunction<Paradigm> GenericType;
    typedef Number<Paradigm> GenericNumericType;
    typedef typename M::PropertiesType PropertiesType;

    template<class Y> using Argument = typename ElementTraits<D>::template Type<Y>;
    template<class Y> using Result = ElementTraits<C>::template Type<Y>;
  private:
    static const CoefficientType _zero;
    DomainType _domain;
    ModelType _model;
  public:

    //@{
    //! \name Constructors and destructors.
    //! \brief Default constructor.
    explicit ScaledFunctionPatch();
    //! \brief Construct a ScaledFunctionPatch<M> over the domain \a d.
    //explicit ScaledFunctionPatch(const DomainType& d);
    explicit ScaledFunctionPatch(const DomainType& d, PropertiesType prp);
    //! \brief Construct a ScaledFunctionPatch<M> over the domain \a d, based on the scaled model \a m.
    explicit ScaledFunctionPatch(const DomainType& d, const ModelType& m);

    explicit ScaledFunctionPatch(const BoxDomainType& d, const Expansion<MultiIndex,FloatValue<PR>>& p, const FloatError<PR>& e, const Sweeper<RawFloat<PR>>& prp);
    explicit ScaledFunctionPatch(const BoxDomainType& d, const Expansion<MultiIndex,RawFloat<PR>>& p, const RawFloat<PR>& e, const Sweeper<RawFloat<PR>>& prp);

    explicit ScaledFunctionPatch(const ScalarFunctionModelType<M>& f);
    ScaledFunctionPatch& operator=(const ScalarFunctionModelType<M>& f);

    //! \brief Construct a ScaledFunctionPatch over the domain \a d from the function \a f.
    explicit ScaledFunctionPatch(const DomainType& d, const ScalarFunctionType<M>& f, PropertiesType prp);
    //@}

    //@{
    //! \name Assignment to constant values.
    //! \brief Set equal to a constant, keeping the same number of arguments.
    ScaledFunctionPatch<M>& operator=(const NumericType& c) { this->_model=c; return *this; }
    //! \brief Set equal to a constant, keeping the same number of arguments.
    ScaledFunctionPatch<M>& operator=(const GenericNumericType c) { this->_model=c; return *this; }
    //@}

    //@{
    //! \name Named constructors.
    //! \brief Construct a zero function over domain \a d.
    static ScaledFunctionPatch<M> zero(const DomainType& d, PropertiesType prp);
    //! \brief Construct a constant quantity in \a as independent variables.
    static ScaledFunctionPatch<M> constant(const DomainType& d, const NumericType& c, PropertiesType prp);
    //! \brief Construct the coordinate \f$x_{j}\f$ over the domain \a d.
    static ScaledFunctionPatch<M> coordinate(const DomainType& d, SizeType j, PropertiesType prp);
    //! \brief Construct a constant quantity in \a as independent variables with value zero and uniform error \a 1
    static ScaledFunctionPatch<M> unit_ball(const DomainType& d, PropertiesType prp);
    //! \brief Construct the coordinate \f$x_{j}\f$ over the domain \a d.
    static ScaledFunctionPatch<M> projection(const BoxDomainType& d, SizeType j, PropertiesType prp);
    //! \brief Return the vector function of coordinates in the range \a js over domain \a d.
    static VectorScaledFunctionPatch<M> projection(const BoxDomainType& d, Range js, PropertiesType prp);
    //! \brief Construct the identity function over the domain \a d.
    static VectorScaledFunctionPatch<M> identity(const DomainType& d, PropertiesType prp);

    //! \brief Construct the quantity \f$c+\sum g_jx_j\f$ over the domain \a d. // DEPRECATED
    static ScaledFunctionPatch<M> affine(const DomainType& d, const CoefficientType& c, const Vector<CoefficientType>& g, PropertiesType prp);
    //! \brief Construct the quantity \f$c+\sum g_jx_j \pm e\f$ over domain \a d. // DEPRECATED
    static ScaledFunctionPatch<M> affine(const DomainType& d, const CoefficientType& x, const Vector<CoefficientType>& g, const ErrorType& e, PropertiesType prp) ;

    //! \brief Return the vector of constants with interval values \a c over domain \a d.
    static Vector<ScaledFunctionPatch<M>> constants(const DomainType& d, const Vector<NumericType>& c, PropertiesType prp);
    //! \brief Return the vector of variables with values \a x over domain \a d.
    static Vector<ScaledFunctionPatch<M>> coordinates(const DomainType& d, PropertiesType prp);
    //! \brief Return the vector of variables in the range \a imin to \a imax with values \a x over domain \a d.
    static Vector<ScaledFunctionPatch<M>> coordinates(const DomainType& d, SizeType imin, SizeType imax, PropertiesType prp);
    //@}

    //@{
    //! \name Prototype constructors.
    friend ScaledFunctionPatchCreator<M> factory(ScaledFunctionPatch<M>const& f) {
        return ScaledFunctionPatchCreator<M>(f.domain(),f.properties()); }
    //! \brief Construct a zero function over the same domain with the same computational properties.
    ScaledFunctionPatch<M> create_zero() const;
    //@}

    //@{
    //! \name Data access
    //! \brief The domain of the quantity.
    const DomainType domain() const { return this->_domain; }
    //! \brief The scaled expansion over a unit box with error bound.
    const ModelType& model() const { return this->_model; }
    //! \brief A reference to the scaled expansion over a unit box with error bound.
    ModelType& model() { return this->_model; }

    //! \brief An over-approximation to the range of the quantity; not necessarily tight.
    const CodomainType codomain() const { return this->_model.codomain(); }
    //! \brief The scaled expansion over a unit box.
    const ExpansionType& expansion() const { return this->_model.expansion(); }
    //! \brief The error of the expansion over the domain.
    const ErrorType& error() const { return this->_model.error(); }
    //! \brief The accuracy parameter used to control approximation of the function model.
    PropertiesType properties() const { return this->_model.properties(); }
    //! \brief The precision of the numbers used.
    PrecisionType precision() const { return this->_model.precision(); }
    //! \brief A reference to the expansion.
    ExpansionType& expansion() { return this->_model.expansion(); }
    //! \brief A reference to the error of the expansion over the domain.
    ErrorType& error() { return this->_model.error(); }

    //! \brief The constant term in the expansion.
    const CoefficientType& value() const { return this->_model.value(); }
    //! \brief The gradient at the centre of the domain.
    const CoefficientType gradient_value(SizeType i) const {
        // FIXME: Cannot be guaranteed to be exact
        Bounds<F> radius(this->_domain[i].radius(),this->model().precision()); Value<F> gradient=this->_model.gradient_value(i);
        return cast_exact(gradient/radius); }

    //! \brief A polynomial representation.
    Polynomial<Bounds<F>> polynomial() const;
    //! \brief A multivalued function equal to the model on the domain.
    ScalarFunctionType<M> function() const;
    //! \brief Cast to a generic function.
    ScalarFunctionType<M> generic() const;

    //! \brief Set the error of the expansion.
    Void set_error(const ErrorType& ne) { this->_model.set_error(ne); }
    //! \brief Set the constant term in the expansion.
    Void set_value(const CoefficientType& c) { this->_model.set_value(c); }

    //! \brief The coefficient of the term in $x^a$.
    const CoefficientType& operator[](const MultiIndex& a) const { return this->_model[a]; }
    //! \brief A read/write reference to the coefficient of the term in $x^a$.
    CoefficientType& operator[](const MultiIndex& a) { return this->_model[a]; }

    //! \brief The number of variables in the argument of the quantity.
    SizeType argument_size() const { return this->_model.argument_size(); }
    //! \brief The maximum degree of terms in the expansion expansion.
    DegreeType degree() const { return this->_model.degree(); }
    //! \brief The number of nonzero terms in the expansion expansion.
    SizeType number_of_nonzeros() const { return this->_model.number_of_nonzeros(); }
    //@}

    //@{
    //! \name Comparison operators.
    Bool operator==(const ScaledFunctionPatch<M>& tv) const;
    //! \brief Inequality operator.
    Bool operator!=(const ScaledFunctionPatch<M>& tv) const { return !(*this==tv); }
    //@}

    //@{
    //! \name Function operations.
    //! \brief An over-approximation to the range of the function.
    RangeType range() const { return this->_model.range(); }
    //! \brief Evaluate the function at the point \a x.
    FloatBounds<PR> operator()(const Vector<FloatBounds<PR>>& x) const;
    FloatBounds<PR> operator()(const Vector<FloatValue<PR>>& x) const;
    FloatApproximation<PR> operator()(const Vector<FloatApproximation<PR>>& x) const;
    ValidatedNumber operator()(const Vector<ValidatedNumber>& x) const;

    //! \brief Compute an approximation to gradient derivative of the function at the point \a x.
    Covector<NumericType> gradient(const Vector<NumericType>& x) const;
    //@}

    //@{
    //! \name Simplification operations.
   //! \brief Remove all terms whose coefficient has magnitude
    //! lower than the cutoff threshold of the quantity.
    ScaledFunctionPatch<M>& simplify() { this->_model.simplify(); return *this; }
    //! \brief Remove all terms whose degree is higher than \a deg or
    //! whose coefficient has magnitude less than \a eps.
    ScaledFunctionPatch<M>& simplify(const PropertiesType& prp) { this->_model.simplify(prp); return *this; }
    //@}

    //@{
    //! \name Accuracy parameters.
    //! \copydoc TaylorModel::set_properties()
    Void set_properties(const PropertiesType& prp) { this->_model.set_properties(prp); }
    //@}

    //@{
    //! \name Non-arithmetic operations.
    //! \brief Restrict to a subdomain.
    Void restrict(const DomainType& d);
    //@}

    //@{
    //! \name Stream input/output operators.
    //! \brief Write to an output stream.
    OutputStream& write(OutputStream& os) const;
    //! \brief Write a full representation to an output stream.
    OutputStream& repr(OutputStream& os) const;
    //! \brief Write to an output stream.
    friend OutputStream& operator<<(OutputStream& os, const ScaledFunctionPatch<M>& x) {
        return x.write(os); }
    //@}

  public:
    Void clobber() { this->_model.clobber(); }
  private:
    friend class TaylorFunctionFactory;
    friend class FunctionMixin<ScaledFunctionPatch<M>, P, D,C>;
    friend class FunctionModelMixin<ScaledFunctionPatch<M>, P, D, C, PR>;
  public:
    template<class X, EnableIf<CanCall<X,M,Vector<X>>> =dummy> Void _compute(X& r, const Vector<X>& a) const;
    template<class X, DisableIf<CanCall<X,M,Vector<X>>> =dummy> Void _compute(X& r, const Vector<X>& a) const;
  private:
    ScaledFunctionPatch<M>* _derivative(SizeType j) const;
    ScaledFunctionPatch<M>* _clone() const;
    ScaledFunctionPatch<M>* _create() const;
    virtual ScaledFunctionPatchFactory<M>* _factory() const;
  public:
    using ScalarMultivariateFunctionModelMixin<ScaledFunctionPatch<M>, typename M::Paradigm, BoxDomainType, typename M::PrecisionType, typename M::ErrorPrecisionType>::_apply;
    template<class OP> static ScaledFunctionPatch<M> _apply(OP op, ScaledFunctionPatch<M> const& f);
    template<class OP> static ScaledFunctionPatch<M> _apply(OP op, ScaledFunctionPatch<M> const& f1, ScaledFunctionPatch<M> const& f2);

    friend ScaledFunctionPatch<M> max(const ScaledFunctionPatch<M>& f1, const ScaledFunctionPatch<M>& f2) { return _apply(Max(),f1,f2); }
    friend ScaledFunctionPatch<M> min(const ScaledFunctionPatch<M>& f1, const ScaledFunctionPatch<M>& f2) { return _apply(Min(),f1,f2); }
    friend ScaledFunctionPatch<M> abs(const ScaledFunctionPatch<M>& f) { return _apply(Abs(),f); }

    friend VectorScaledFunctionPatch<M> join(const ScaledFunctionPatch<M>& f1, const ScaledFunctionPatch<M>& f2);
    friend VectorScaledFunctionPatch<M> combine(const ScaledFunctionPatch<M>& f1, const ScaledFunctionPatch<M>& f2);

    //! \brief Restrict to a subdomain.
    friend ScaledFunctionPatch<M> restriction(const ScaledFunctionPatch<M>& f, const BoxDomainType& dom) {
        if(not(subset(dom,f.domain()))) { ARIADNE_THROW(DomainException,"restiction(ScaledFunctionPatch<M>,BoxDomainType)","f="<<f<<", dom="<<dom); }
        return unchecked_compose(f,ScaledFunctionPatch<M>::identity(dom,f.properties())); }


    //! \brief Restrict a component to a subdomain
    // To scale from a model on [a,b] to a model on [c,d], use scale factor s=(d-c)/(b-a)
    // and translation t=((c+d)-(a+b))/(b-a)
    // Because we are scaling the model on [-1,+1], this is not the same as
    // the mapping taking [a,b] to [c,d]
    friend ScaledFunctionPatch<M> partial_restriction(const ScaledFunctionPatch<M>& f, SizeType k, const IntervalDomainType& ivl) {
        BoxDomainType dom=f.domain(); dom[k]=ivl; return restriction(f,dom); }

    //! \brief Pre-compose with a projection from the Cartesian product of the given domains.
    friend ScaledFunctionPatch<M> embed(const BoxDomainType& dom1, const ScaledFunctionPatch<M>& f2,const BoxDomainType& dom3) {
        return ScaledFunctionPatch<M>(product(dom1,f2.domain(),dom3),embed(dom1.size(),f2.model(),dom3.size())); }

    //! \brief Split the domain into two.
    friend Pair<ScaledFunctionPatch<M>,ScaledFunctionPatch<M>> split(const ScaledFunctionPatch<M>& f, SizeType j) {
        Pair<ModelType,ModelType> models=split(f.model(),j); Pair<BoxDomainType,BoxDomainType> subdomains=split(f.domain(),j);
        return make_pair(ScaledFunctionPatch<M>(subdomains.first,models.first), ScaledFunctionPatch<M>(subdomains.second,models.second)); }

/*
    friend ScaledFunctionPatch<M>& operator+=(ScaledFunctionPatch<M>& f, const NumericType& c) { f.model()+=c; return f; }
    friend ScaledFunctionPatch<M>& operator*=(ScaledFunctionPatch<M>& f, const NumericType& c) { f.model()*=c; return f; }
    friend ScaledFunctionPatch<M> operator-(const ScaledFunctionPatch<M>& f) { return apply(Neg(),f); }
    friend ScaledFunctionPatch<M> operator+(const ScaledFunctionPatch<M>& f1, const ScaledFunctionPatch<M>& f2) { return apply(Add(),f1,f2); }
    friend ScaledFunctionPatch<M> operator-(const ScaledFunctionPatch<M>& f1, const ScaledFunctionPatch<M>& f2) { return apply(Sub(),f1,f2); }
    friend ScaledFunctionPatch<M> operator*(const ScaledFunctionPatch<M>& f1, const ScaledFunctionPatch<M>& f2) { return apply(Mul(),f1,f2); }
    friend ScaledFunctionPatch<M> operator/(const ScaledFunctionPatch<M>& f1, const ScaledFunctionPatch<M>& f2) { return apply(Div(),f1,f2); }
*/
    friend ScaledFunctionPatch<M> derivative(const ScaledFunctionPatch<M>& f, SizeType k) {
        return ScaledFunctionPatch<M>(f.domain(),derivative(f.model(),k)/rad(f.domain()[k])); }
    friend ScaledFunctionPatch<M> antiderivative(const ScaledFunctionPatch<M>& f, SizeType k) {
        return ScaledFunctionPatch<M>(f.domain(),antiderivative(f.model(),k)*rad(f.domain()[k])); }
    friend ScaledFunctionPatch<M> antiderivative(const ScaledFunctionPatch<M>& f, SizeType k, const NumericType& c) {
        ARIADNE_ASSERT(k<f.argument_size()); ARIADNE_ASSERT(decide(contains(f.domain()[k],c)));
        ScaledFunctionPatch<M> g = antiderivative(f,k);
        VectorScaledFunctionPatch<M> h ( VectorScaledFunctionPatch<M>::identity(f.domain(),f.properties()) );
        h[k] = ScaledFunctionPatch<M>::constant(f.domain(),c,f.properties());
        return g-compose(g,h); }
    friend ScaledFunctionPatch<M> antiderivative(const ScaledFunctionPatch<M>& f, SizeType k, const GenericNumericType& c) {
        return antiderivative(f,k,NumericType(c,f.precision())); }
        
    friend ScaledFunctionPatch<M> partial_evaluate(const ScaledFunctionPatch<M>& f, SizeType k, const NumericType& c) {
        ARIADNE_ASSERT(decide(contains(f.domain()[k],c)));
        return ScaledFunctionPatch<M>(remove(f.domain(),k),partial_evaluate(f.model(),k,unscale(c,f.domain()[k]))); }
    friend ScaledFunctionPatch<M> partial_evaluate(const ScaledFunctionPatch<M>& f, SizeType k, const GenericNumericType& c) {
        ARIADNE_ASSERT(decide(contains(f.domain()[k],c)));
        return partial_evaluate(f,k,NumericType(c,f.precision())); }
    friend NumericType evaluate(const ScaledFunctionPatch<M>& f, const Vector<NumericType>& x) {
        if(!definitely(contains(f.domain(),x))) { ARIADNE_THROW(DomainException,"evaluate(f,x) with f="<<f<<", x="<<x,"x is not an element of f.domain()="<<f.domain()); }
        return unchecked_evaluate(f,x); }
    friend NumericType evaluate(const ScaledFunctionPatch<M>& f, const Vector<GenericNumericType>& x) {
        return evaluate(f,Vector<NumericType>(x,f.precision())); }
    friend NumericType unchecked_evaluate(const ScaledFunctionPatch<M>& f, const Vector<NumericType>& x) {
        return evaluate(f.model(),unscale(x,f.domain())); }
    friend NumericType unchecked_evaluate(const ScaledFunctionPatch<M>& f, const Vector<GenericNumericType>& x) {
        return unchecked_evaluate(f,Vector<NumericType>(x,f.precision())); }

    friend NormType norm(const ScaledFunctionPatch<M>& f) {
        return norm(f.model()); }
    friend NormType distance(const ScaledFunctionPatch<M>& f1, const ScaledFunctionPatch<M>& f2) {
        return norm(f1-f2); }
    friend NormType distance(const ScaledFunctionPatch<M>& f1, const ScalarMultivariateFunction<P>& f2) {
        return distance(f1,f1.create(f2)); }

    friend Polynomial<NumericType> polynomial(const ScaledFunctionPatch<M>& tfn) { return tfn.polynomial(); }

};

template<class M> template<class X, EnableIf<CanCall<X,M,Vector<X>>>> Void ScaledFunctionPatch<M>::_compute(X& r, const Vector<X>& a) const {
    r = this->_model(unscale(a,this->_domain));
}
template<class M> template<class X, DisableIf<CanCall<X,M,Vector<X>>>> Void ScaledFunctionPatch<M>::_compute(X& r, const Vector<X>& a) const {
    assert(false);
}

template<class FP1, class FP2> Void check_function_patch_domain(String const& op_str, const FP1& f1, const FP2& f2) {
    ARIADNE_ASSERT_MSG(!is_empty(intersection(f1.domain(),f2.domain())),
                    op_str<<"((Vector)ScaledFunctionPatch<M> f1, (Vector)ScaledFunctionPatch<M> f2) with f1="<<f1<<" f2="<<f2<<
                    ": domains are disjoint");
}
template<class M> template<class OP> ScaledFunctionPatch<M> ScaledFunctionPatch<M>::_apply(OP op, ScaledFunctionPatch<M> const& f) {
    return ScaledFunctionPatch<M>(f.domain(),op(f.model()));
}
template<class M> template<class OP> ScaledFunctionPatch<M> ScaledFunctionPatch<M>::_apply(OP op, ScaledFunctionPatch<M> const& f1, ScaledFunctionPatch<M> const& f2) {
    check_function_patch_domain(to_str(op),f1,f2);
    if(f1.domain()==f2.domain()) {
        return ScaledFunctionPatch<M>(f1.domain(),op(f1.model(),f2.model()));
    } else {
        BoxDomainType domain=intersection(f1.domain(),f2.domain());
        return ScaledFunctionPatch<M>(domain,op(restriction(f1,domain).model(),restriction(f2,domain).model()));
    }
}


template<class M> template<class OP> ScaledFunctionPatch<M> AlgebraOperations<ScaledFunctionPatch<M>>::apply(OP op, ScaledFunctionPatch<M> const& f1, ScaledFunctionPatch<M> const& f2) {
    assert(f1.domain()==f2.domain()); return ScaledFunctionPatch<M>(f1.domain(),op(f1.model(),f2.model()));
}
template<class M> template<class OP> ScaledFunctionPatch<M> AlgebraOperations<ScaledFunctionPatch<M>>::apply(OP op, ScaledFunctionPatch<M> const& f) {
    return ScaledFunctionPatch<M>(f.domain(),op(f.model()));
}
template<class M> template<class OP> ScaledFunctionPatch<M> AlgebraOperations<ScaledFunctionPatch<M>>::apply(OP op, ScaledFunctionPatch<M> const& f1, typename M::NumericType const& c2) {
    return ScaledFunctionPatch<M>(f1.domain(),op(f1.model(),c2));
}
template<class M> template<class OP> ScaledFunctionPatch<M> AlgebraOperations<ScaledFunctionPatch<M>>::apply(OP op, typename M::NumericType const& c1, ScaledFunctionPatch<M> const& f2) {
    return ScaledFunctionPatch<M>(f2.domain(),op(c1,f2.model()));
}
template<class M> ScaledFunctionPatch<M> AlgebraOperations<ScaledFunctionPatch<M>>::apply(Pow op, ScaledFunctionPatch<M> const& f, Int n) {
    return ScaledFunctionPatch<M>(f.domain(),op(f.model(),n));
}

//! \brief Test if the function models have the same representation.
template<class M> Bool same(const ScaledFunctionPatch<M>& tv1, const ScaledFunctionPatch<M>& tv2) {
    return tv1.domain()==tv2.domain() && same(tv1.model(),tv2.model());
}
//! \brief Test if the quantity is a better approximation than \a t throughout the domain.
template<class M> Bool refines(const ScaledFunctionPatch<M>& tv1, const ScaledFunctionPatch<M>& tv2) {
    if(tv1.domain()==tv2.domain()) { return refines(tv1.model(),tv2.model()); }
    if(subset(tv2.domain(),tv1.domain())) { return refines(restriction(tv1,tv2.domain()).model(),tv2.model()); }
    else { return false; }
}
//! \brief Test if the function models are inconsistent with representing the same exact function.
template<class M> Bool inconsistent(const ScaledFunctionPatch<M>& tv1, const ScaledFunctionPatch<M>& tv2) {
    if(tv1.domain()==tv2.domain()) {
        return inconsistent(tv1.model(),tv2.model());
    } else {
        BoxDomainType domain=intersection(tv1.domain(),tv2.domain());
        return inconsistent(restriction(tv1,domain).model(),restriction(tv2,domain).model());
    }
}
//! \brief Compute an over-approximation to the common refinement of \a x1 and \a x2.
template<class M> ScaledFunctionPatch<M> refinement(const ScaledFunctionPatch<M>& tv1, const ScaledFunctionPatch<M>& tv2) {
    ARIADNE_ASSERT(tv1.domain()==tv2.domain());
    return ScaledFunctionPatch<M>(tv1.domain(),refinement(tv1.model(),tv2.model()));
}
//! \brief Remove the error term.
template<class M> ScaledFunctionPatch<M> midpoint(const ScaledFunctionPatch<M>& f) {
    M tm=f.model();
    tm.set_error(0u);
    return ScaledFunctionPatch<M>(f.domain(),tm);
}








/*! \ingroup FunctionModelSubModule
 *  \brief A Taylor function model with multivalued codomain built from the TaylorModel class.
 *
 *  See also TaylorModel, ScaledFunctionPatch<M>, ValidatedVectorMultivariateTaylorFunctionModelDP.
 */
template<class M> class VectorScaledFunctionPatch
    : public VectorMultivariateFunctionModelMixin<VectorScaledFunctionPatch<M>,typename M::Paradigm,BoxDomainType,typename M::PrecisionType,typename M::ErrorPrecisionType>
{
    friend class VectorScaledFunctionPatchElementReference<M>;
    typedef BoxDomainType D;
    typedef BoxDomainType C;
    typedef typename M::Paradigm P;
    typedef typename M::RawFloatType F;
    typedef typename M::PrecisionType PR;
    typedef typename M::ErrorPrecisionType PRE;
  public:
    typedef D DomainType;
    typedef M ModelType;
    typedef P Paradigm;
    typedef PR PrecisionType;
    typedef PRE ErrorPrecisionType;
    typedef Box<typename ModelType::CodomainType> CodomainType;
    typedef Box<typename ModelType::RangeType> RangeType;
    typedef typename ModelType::ExpansionType ExpansionType;
    typedef typename ModelType::CoefficientType CoefficientType;
    typedef typename ModelType::ErrorType ErrorType;
    typedef typename ModelType::NumericType NumericType;
    typedef Number<Paradigm> GenericNumericType;
    typedef typename ModelType::NormType NormType;
    typedef typename M::PropertiesType PropertiesType;

    typedef FloatApproximation<PR> ApproximateNumericType;
    typedef FloatBounds<PR> ValidatedNumericType;
    typedef FloatValue<PR> ExactNumericType;

    template<class Y> using Argument = typename ElementTraits<D>::template Type<Y>;
    template<class Y> using Result = ElementTraits<C>::template Type<Y>;

    //! \brief Default constructor constructs a Taylor model of order zero with no arguments and no result variables.
    VectorScaledFunctionPatch<M>();

    //! \brief Construct the zero vector function over an unspecified domain.
    explicit VectorScaledFunctionPatch<M>(SizeType result_size);

    //! \brief Construct from a result size and a domain.
    VectorScaledFunctionPatch<M>(SizeType result_size, const BoxDomainType& domain, PropertiesType properties);

    //! \brief Construct a vector function all of whose components are the same.
    VectorScaledFunctionPatch<M>(SizeType result_size, const ScaledFunctionPatch<M>& scalar_function);

    //! \brief Construct from a domain and the expansion.
   VectorScaledFunctionPatch<M>(const BoxDomainType& domain,
                                const Vector<Expansion<MultiIndex,FloatValue<PR>>>& expansion,
                                PropertiesType properties);

    //! \brief Construct from a domain, and expansion and errors.
    VectorScaledFunctionPatch<M>(const BoxDomainType& domain,
                                 const Vector<Expansion<MultiIndex,FloatValue<PR>>>& expansion,
                                 const Vector<FloatError<PR>>& error,
                                 PropertiesType properties);

    //! \brief Construct from a domain, and expansion and errors.
    VectorScaledFunctionPatch<M>(const BoxDomainType& domain,
                                const Vector<Expansion<MultiIndex,RawFloat<PR>>>& expansion,
                                const Vector<RawFloat<PR>>& error,
                                PropertiesType properties);

    //! \brief Construct from a domain, and expansion and errors.
    VectorScaledFunctionPatch<M>(const BoxDomainType& domain,
                                 const Vector<Expansion<MultiIndex,RawFloat<PR>>>& expansion,
                                 PropertiesType properties);

    //! \brief Construct from a domain and the models.
    explicit VectorScaledFunctionPatch<M>(const BoxDomainType& domain, const Vector< ModelType >& variables);

    //! \brief Construct from a \a domain, a \a function, and \a properties determining the accuracy.
    VectorScaledFunctionPatch<M>(const BoxDomainType& domain,
                                 const VectorFunctionType<M>& function,
                                 const PropertiesType& properties);

    //! \brief Construct from a vector of scalar Taylor functions.
    explicit VectorScaledFunctionPatch<M>(const Vector<ScaledFunctionPatch<M>>& components);

    //! \brief Construct from a list of scalar Taylor functions.
    explicit VectorScaledFunctionPatch<M>(const List<ScaledFunctionPatch<M>>& components);

    //! \brief Construct from an initializer list of scalar Taylor functions.
    VectorScaledFunctionPatch<M>(InitializerList<ScaledFunctionPatch<M>> components);

    //! \brief Construct from a vector expression.
    template<class E> explicit VectorScaledFunctionPatch<M>(const VectorExpression<E>& ve);

    explicit VectorScaledFunctionPatch<M> (const VectorFunctionModelType<M>& f);
    VectorScaledFunctionPatch<M>& operator=(const VectorFunctionModelType<M>& f);

    //! \brief Equality operator.
    Bool operator==(const VectorScaledFunctionPatch<M>& p) const;
    //! \brief Inequality operator.
    Bool operator!=(const VectorScaledFunctionPatch<M>& p) const;

    // Data access
    //! \brief The properties used to control approximation of the function model.
    PropertiesType properties() const;
    //! \brief The precision of the numbers used.
    PrecisionType precision() const;
    //! \brief Set the properties used to control approximation of the function model.
    Void set_properties(PropertiesType prp);
    //! \brief The data used to define the domain of the Taylor model.
    const BoxDomainType domain() const;
    //! \brief A rough bound for the range of the function.
    const BoxDomainType codomain() const;
    //! \brief The centre of the Taylor model.
    const Vector<CoefficientType> centre() const;
    //! \brief The range of the Taylor model.
    const RangeType range() const;
    //! \brief The data used to define the Taylor models.
    const Vector<ModelType>& models() const;
    Vector<ModelType>& models();
    //! \brief The data used to define the centre of the Taylor models.
    const Vector<Expansion<MultiIndex,FloatValue<PR>>> expansions() const;

    //! \brief The \a i<sup>th</sup> Taylor model used to define the function.
    const ModelType& model(SizeType i) const;
    //! \brief The \a i<sup>th</sup> Taylor model used to define the function.
    ModelType& model(SizeType i);

    //! \brief The size of the argument.
    SizeType argument_size() const;
    //! \brief The size of the result.
    SizeType result_size() const;

    // For compatibility with Vector.
    //! \brief The number of scalar function components.
    SizeType size() const;
    //! \brief A null function compatible with the elements.
    ScaledFunctionPatch<M> zero_element() const;
    //! \brief Get the \a ith Taylor variable
    ScaledFunctionPatch<M> get(SizeType i) const;
    //! \brief Set the \a ith Taylor variable
    Void set(SizeType i, const ScaledFunctionPatch<M>& te);
    //! \brief The \a ith Taylor variable
    ScaledFunctionPatch<M> const operator[](SizeType i) const;
    //! \brief The \a ith Taylor variable
    VectorScaledFunctionPatchElementReference<M> operator[](SizeType i);


    //! \brief Evaluate the Taylor model at the point \a x.
    Vector<ValidatedNumericType> operator()(const Vector<ValidatedNumericType>& x) const;
    Vector<ApproximateNumericType> operator()(const Vector<ApproximateNumericType>& x) const;
    Vector<ValidatedNumericType> operator()(const Vector<ExactNumericType>& x) const;
    Vector<ValidatedNumber> operator()(const Vector<ValidatedNumber>& x) const;
    //! \brief Compute an approximation to Jacobian derivative of the Taylor model sat the point \a x.
    Matrix<NumericType> jacobian(const Vector<NumericType>& x) const;

    //! \brief Simplify the representation.
    VectorScaledFunctionPatch<M>& simplify();
    //! \brief Simplify the representation as specified by \a properties.
    VectorScaledFunctionPatch<M>& simplify(const PropertiesType& properties);
    //! \brief Set the error to zero.
    Void clobber();

    friend ScaledFunctionPatchCreator<M> factory(VectorScaledFunctionPatch<M>const& f) {
        return ScaledFunctionPatchCreator<M>(f.domain(),f.properties()); }

    //! \brief The constant Taylor model with range \a r and argument domain \a d.
    static VectorScaledFunctionPatch<M> constant(const BoxDomainType& d, const Vector<NumericType>& r, PropertiesType prp);
    //! \brief The identity Taylor model on domain \a d.
    static VectorScaledFunctionPatch<M> identity(const BoxDomainType& d, PropertiesType prp);
    //! \brief Return the vector of variables in the range with values \a x over domain \a d.
    static VectorScaledFunctionPatch<M> projection(const BoxDomainType& d, Range js, PropertiesType prp);

    //! \brief Convert to an interval polynomial.
    Vector<Polynomial<FloatBounds<PR>>> polynomials() const;
    //! \brief The vector of roundoff/truncation errors of each component.
    Vector<ErrorType> const errors() const;
    //! \brief The maximum roundoff/truncation error of the components.
    ErrorType const error() const;
    //! \brief A multivalued function equal to the model on the domain.
    VectorFunctionType<M> function() const;
    //! \brief Cast to a generic function.
    VectorFunctionType<M> generic() const;

    //! \brief Truncate terms higher than \a bd.
    VectorScaledFunctionPatch<M>& truncate(const MultiIndexBound& bd);
    //! \brief Restrict to a subdomain.
    Void restrict(const BoxDomainType& d);
    //! \brief Adjoin a scalar function.
    Void adjoin(const ScaledFunctionPatch<M>& sf);

    //! \brief Write to an output stream.
    OutputStream& write(OutputStream& os) const;

    //! \brief Write a full representation to an output stream.
    OutputStream& repr(OutputStream& os) const;

  private:
    Void _compute_jacobian() const;
    Void _set_argument_size(SizeType n);
    SizeType _compute_maximum_component_size() const;
    virtual ScalarScaledFunctionPatch<M>* _get(SizeType i) const { return new ScaledFunctionPatch<M>(this->_domain,this->_models[i]); }
    virtual VectorScaledFunctionPatch<M>* _clone() const;
    virtual VectorScaledFunctionPatch<M>* _create() const;
    virtual ScaledFunctionPatchFactory<M>* _factory() const;
  private:
    friend class VectorFunctionMixin<VectorScaledFunctionPatch<M>,P,BoxDomainType>;
    friend class TaylorFunctionFactory;
  public:
    template<class X, EnableIf<CanCall<X,M,Vector<X>>> =dummy> Void _compute(Vector<X>& r, const Vector<X>& a) const;
    template<class X, DisableIf<CanCall<X,M,Vector<X>>> =dummy> Void _compute(Vector<X>& r, const Vector<X>& a) const;
  private:
    /* Domain of definition. */
    BoxDomainType _domain;
    Vector< ModelType > _models;

  public:
    //! \brief Compute the function \f$(f \oplus g)(x)=(f(x),g(x))\f$.
    friend VectorScaledFunctionPatch<M> join(const VectorScaledFunctionPatch<M>& f1, const ScaledFunctionPatch<M>& f2) {
        ARIADNE_ASSERT_MSG(f1.domain()==f2.domain(),"f1="<<f1<<", f2="<<f2);
        return VectorScaledFunctionPatch<M>(f1.domain(),join(f1.models(),f2.model()));
    }
    friend VectorScaledFunctionPatch<M> join(const VectorScaledFunctionPatch<M>& f, const VectorScaledFunctionPatch<M>& g) {
        ARIADNE_ASSERT(f.domain()==g.domain());
        return VectorScaledFunctionPatch<M>(f.domain(),join(f.models(),g.models()));
    }
    friend VectorScaledFunctionPatch<M> join(const ScaledFunctionPatch<M>& f1, const ScaledFunctionPatch<M>& f2) {
        ARIADNE_ASSERT(f1.domain()==f2.domain());
        return VectorScaledFunctionPatch<M>(f1.domain(),{f1.model(),f2.model()});
    }
    friend VectorScaledFunctionPatch<M> join(const ScaledFunctionPatch<M>& f1, const VectorScaledFunctionPatch<M>& f2) {
        ARIADNE_ASSERT(f1.domain()==f2.domain());
        return VectorScaledFunctionPatch<M>(f1.domain(),join(f1.model(),f2.models()));
    }

    //! \brief Compute the function \f$(f\otimes g)(x,y)=(f(x),g(y))\f$.
    friend VectorScaledFunctionPatch<M> combine(const ScaledFunctionPatch<M>& f1, const ScaledFunctionPatch<M>& f2) {
        return VectorScaledFunctionPatch<M>(product(f1.domain(),f2.domain()),combine(Vector<M>(1u,f1.model()),{f2.model()}));
    }
    friend VectorScaledFunctionPatch<M> combine(const ScaledFunctionPatch<M>& f1, const VectorScaledFunctionPatch<M>& f2) {
        return VectorScaledFunctionPatch<M>(product(f1.domain(),f2.domain()),combine({f1.model()},f2.models()));
    }
    friend VectorScaledFunctionPatch<M> combine(const VectorScaledFunctionPatch<M>& f1, const ScaledFunctionPatch<M>& f2) {
        return VectorScaledFunctionPatch<M>(product(f1.domain(),f2.domain()),combine(f1.models(),{f2.model()}));
    }
    friend VectorScaledFunctionPatch<M> combine(const VectorScaledFunctionPatch<M>& f1, const VectorScaledFunctionPatch<M>& f2) {
        return VectorScaledFunctionPatch<M>(product(f1.domain(),f2.domain()),combine(f1.models(),f2.models()));
    }

    friend VectorScaledFunctionPatch<M> embed(const VectorScaledFunctionPatch<M>& f, const IntervalDomainType& d) {
        return embed(BoxDomainType(),f,BoxDomainType(1u,d));
    }
    friend VectorScaledFunctionPatch<M> embed(const VectorScaledFunctionPatch<M>& f, const BoxDomainType& d) {
        return embed(BoxDomainType(),f,d);
    }
    friend VectorScaledFunctionPatch<M> embed(const BoxDomainType& d, const VectorScaledFunctionPatch<M>& f) {
        return embed(d,f,BoxDomainType());
    }
    friend VectorScaledFunctionPatch<M> embed(const BoxDomainType& d1, const VectorScaledFunctionPatch<M>& f, const BoxDomainType& d2) {
        return VectorScaledFunctionPatch<M>(product(d1,f.domain(),d2),embed(d1.size(),f.models(),d2.size()));
    }

    friend VectorScaledFunctionPatch<M> restriction(const VectorScaledFunctionPatch<M>& f, const BoxDomainType& d) {
        ARIADNE_ASSERT_MSG(subset(d,f.domain()),"Cannot restriction "<<f<<" to non-sub-domain "<<d);
        if(d==f.domain()) { return f; }
        VectorScaledFunctionPatch<M> r(f.result_size(),d,f.properties());
        for(SizeType i=0; i!=r.result_size(); ++i) {
            r.set(i,restriction(f[i],d));
        }
        return r;
    }
    friend VectorScaledFunctionPatch<M> partial_restriction(const VectorScaledFunctionPatch<M>& tf, SizeType k, const IntervalDomainType& d) {
        VectorScaledFunctionPatch<M> r(tf.result_size(),tf.domain(),tf.properties());
        for(SizeType i=0; i!=tf.result_size(); ++i) {
            r[i]=partial_restriction(tf[i],k,d);
        }
        return r;
    }
    friend VectorScaledFunctionPatch<M> restriction(const VectorScaledFunctionPatch<M>& tf, SizeType k, const IntervalDomainType& d) {
        return partial_restriction(tf,k,d);
    }
    friend Pair<VectorScaledFunctionPatch<M>,VectorScaledFunctionPatch<M>> split(const VectorScaledFunctionPatch<M>& tf, SizeType j) {
        Pair<Vector<M>,Vector<M>> models=split(tf.models(),j);
        Pair<BoxDomainType,BoxDomainType> subdomains=split(tf.domain(),j);
        return make_pair(VectorScaledFunctionPatch<M>(subdomains.first,models.first),
                        VectorScaledFunctionPatch<M>(subdomains.second,models.second));

    }

    friend VectorScaledFunctionPatch<M>& operator+=(VectorScaledFunctionPatch<M>& f, const VectorScaledFunctionPatch<M>& g) {
        ARIADNE_ASSERT(f.result_size()==g.result_size());
        ARIADNE_ASSERT(subset(f.domain(),g.domain()));
        ARIADNE_ASSERT(f.domain()==g.domain());
        f.models()+=g.models();
        return f;
    }
    friend VectorScaledFunctionPatch<M>& operator-=(VectorScaledFunctionPatch<M>& f, const VectorScaledFunctionPatch<M>& g) {
        ARIADNE_ASSERT(f.result_size()==g.result_size());
        ARIADNE_ASSERT(subset(f.domain(),g.domain()));
        ARIADNE_ASSERT(f.domain()==g.domain());
        f.models()+=g.models();
        return f;
    }
    friend VectorScaledFunctionPatch<M>& operator+=(VectorScaledFunctionPatch<M>& f, const Vector<NumericType>& c) {
        ARIADNE_ASSERT(f.result_size()==c.size());
        f.models()+=c;
        return f;
    }
    friend VectorScaledFunctionPatch<M>& operator-=(VectorScaledFunctionPatch<M>& f, const Vector<NumericType>& c) {
        ARIADNE_ASSERT(f.result_size()==c.size());
        f.models()-=c;
        return f;
    }
    friend VectorScaledFunctionPatch<M>& operator*=(VectorScaledFunctionPatch<M>& f, const NumericType& c) {
        f.models()*=c;
        return f;
    }
    friend VectorScaledFunctionPatch<M>& operator/=(VectorScaledFunctionPatch<M>& f, const NumericType& c) {
        f.models()/=c;
        return f;
    }

    template<class OP> static VectorScaledFunctionPatch<M> _apply(OP op, const VectorScaledFunctionPatch<M>& f1, const VectorScaledFunctionPatch<M>& f2);
    template<class OP> static VectorScaledFunctionPatch<M> _apply(OP op, const ScaledFunctionPatch<M>& f1, const VectorScaledFunctionPatch<M>& f2);
    template<class OP> static VectorScaledFunctionPatch<M> _apply(OP op, const VectorScaledFunctionPatch<M>& f1, const ScaledFunctionPatch<M>& f2);
    friend VectorScaledFunctionPatch<M> operator+(const VectorScaledFunctionPatch<M>& f1, const VectorScaledFunctionPatch<M>& f2) {
        return _apply(Plus(),f1,f2); }
    friend VectorScaledFunctionPatch<M> operator-(const VectorScaledFunctionPatch<M>& f1, const VectorScaledFunctionPatch<M>& f2) {
        return _apply(Minus(),f1,f2); }
    friend VectorScaledFunctionPatch<M> operator*(const ScaledFunctionPatch<M>& f1, const VectorScaledFunctionPatch<M>& f2) {
        return _apply(Times(),f1,f2); }
    friend VectorScaledFunctionPatch<M> operator*(const VectorScaledFunctionPatch<M>& f1, const ScaledFunctionPatch<M>& f2) {
        return _apply(Times(),f1,f2); }
    friend VectorScaledFunctionPatch<M> operator/(const VectorScaledFunctionPatch<M>& f1, const ScaledFunctionPatch<M>& f2) {
        return _apply(Divides(),f1,f2); }


    friend VectorScaledFunctionPatch<M> operator+(const VectorScaledFunctionPatch<M>& f) {
        return VectorScaledFunctionPatch<M>(f.domain(),Vector<M>(+f.models())); }
    friend VectorScaledFunctionPatch<M> operator-(const VectorScaledFunctionPatch<M>& f) {
        return VectorScaledFunctionPatch<M>(f.domain(),Vector<M>(-f.models())); }
    friend VectorScaledFunctionPatch<M> operator*(const NumericType& c, const VectorScaledFunctionPatch<M>& f) {
        return VectorScaledFunctionPatch<M>(f.domain(),Vector<M>(f.models()*c)); }
    friend VectorScaledFunctionPatch<M> operator*(const VectorScaledFunctionPatch<M>& f, const NumericType& c) {
        return VectorScaledFunctionPatch<M>(f.domain(),Vector<M>(f.models()*c)); }
    friend VectorScaledFunctionPatch<M> operator/(const VectorScaledFunctionPatch<M>& f, const NumericType& c) {
        return VectorScaledFunctionPatch<M>(f.domain(),Vector<M>(f.models()/c)); }
    friend VectorScaledFunctionPatch<M> operator+(const VectorScaledFunctionPatch<M>& f, const Vector<NumericType>& c) {
        return VectorScaledFunctionPatch<M>(f.domain(),Vector<M>(f.models()+c)); }
    friend VectorScaledFunctionPatch<M> operator-(const VectorScaledFunctionPatch<M>& f, const Vector<NumericType>& c) {
        return VectorScaledFunctionPatch<M>(f.domain(),Vector<M>(f.models()-c)); }

    friend VectorScaledFunctionPatch<M> operator*(const Matrix<NumericType>& A, const VectorScaledFunctionPatch<M>& f) {
        ARIADNE_PRECONDITION(A.column_size()==f.size());
        return VectorScaledFunctionPatch<M>(f.domain(),A*f.models());
    }
    friend VectorScaledFunctionPatch<M> operator+(const VectorMultivariateFunction<P>& f1, const VectorScaledFunctionPatch<M>& tf2) {
        return VectorScaledFunctionPatch<M>(tf2.domain(),f1,tf2.properties())+tf2; }
    friend VectorScaledFunctionPatch<M> operator-(const VectorMultivariateFunction<P>& f1, const VectorScaledFunctionPatch<M>& tf2) {
        return VectorScaledFunctionPatch<M>(tf2.domain(),f1,tf2.properties())-tf2; }
    friend VectorScaledFunctionPatch<M> operator*(const ScalarMultivariateFunction<P>& f1, const VectorScaledFunctionPatch<M>& tf2) {
        return ScaledFunctionPatch<M>(tf2.domain(),f1,tf2.properties())*tf2; }
    friend VectorScaledFunctionPatch<M> operator*(const VectorMultivariateFunction<P>& f1, const ScaledFunctionPatch<M>& tf2) {
        return VectorScaledFunctionPatch<M>(tf2.domain(),f1,tf2.properties())*tf2; }
    friend VectorScaledFunctionPatch<M> operator/(const VectorMultivariateFunction<P>& f1, const ScaledFunctionPatch<M>& tf2) {
        return VectorScaledFunctionPatch<M>(tf2.domain(),f1,tf2.properties())/tf2; }
    friend VectorScaledFunctionPatch<M> operator+(const VectorScaledFunctionPatch<M>& tf1, const VectorMultivariateFunction<P>& f2) {
        return tf1+VectorScaledFunctionPatch<M>(tf1.domain(),f2,tf1.properties()); }
    friend VectorScaledFunctionPatch<M> operator-(const VectorScaledFunctionPatch<M>& tf1, const VectorMultivariateFunction<P>& f2) {
        return tf1-VectorScaledFunctionPatch<M>(tf1.domain(),f2,tf1.properties()); }
    friend VectorScaledFunctionPatch<M> operator*(const ScaledFunctionPatch<M>& tf1, const VectorMultivariateFunction<P>& f2) {
        return tf1*VectorScaledFunctionPatch<M>(tf1.domain(),f2,tf1.properties()); }
    friend VectorScaledFunctionPatch<M> operator*(const VectorScaledFunctionPatch<M>& tf1, const ScalarMultivariateFunction<P>& f2) {
        return tf1*ScaledFunctionPatch<M>(tf1.domain(),f2,tf1.properties()); }
    friend VectorScaledFunctionPatch<M> operator/(const VectorScaledFunctionPatch<M>& tf1, const ScalarMultivariateFunction<P>& f2) {
        return tf1/ScaledFunctionPatch<M>(tf1.domain(),f2,tf1.properties()); }

    friend VectorScaledFunctionPatch<M> partial_evaluate(const VectorScaledFunctionPatch<M>& tf, SizeType k, const NumericType& c) {
        ARIADNE_ASSERT(decide(contains(tf.domain()[k],c)));
        return VectorScaledFunctionPatch<M>(remove(tf.domain(),k),partial_evaluate(tf.models(),k,unscale(c,tf.domain()[k])));
    }
    friend VectorScaledFunctionPatch<M> partial_evaluate(const VectorScaledFunctionPatch<M>& tf, SizeType k, const GenericNumericType& c) {
        return partial_evaluate(tf,k,NumericType(c,tf.precision())); }
        
    friend Vector<NumericType> evaluate(const VectorScaledFunctionPatch<M>& f, const Vector<NumericType>& x) {
        if(!definitely(contains(f.domain(),x))) {
            ARIADNE_THROW(DomainException,"evaluate(f,x) with f="<<f<<", x="<<x,"x is not a subset of f.domain()="<<f.domain());
        }
        return unchecked_evaluate(f,x);
    }
    friend Vector<GenericNumericType> evaluate(const VectorScaledFunctionPatch<M>& f, const Vector<GenericNumericType>& x) {
        return evaluate(f,Vector<NumericType>(x,f.precision()));
    }
    friend Vector<NumericType> unchecked_evaluate(const VectorScaledFunctionPatch<M>& f, const Vector<NumericType>& x) {
        return evaluate(f.models(),unscale(x,f.domain()));
    }
    friend Vector<GenericNumericType> unchecked_evaluate(const VectorScaledFunctionPatch<M>& f, const Vector<GenericNumericType>& x) {
        return unchecked_evaluate(f,Vector<NumericType>(x,f.precision()));
    }

    friend ScaledFunctionPatch<M> compose(const ScalarMultivariateFunction<P>& g, const VectorScaledFunctionPatch<M>& f) {
        return ScaledFunctionPatch<M>(f.domain(),g.evaluate(f.models()));
    }
    friend ScaledFunctionPatch<M> compose(const ScaledFunctionPatch<M>& g, const VectorScaledFunctionPatch<M>& f) {
        if(!subset(f.codomain(),g.domain())) {
            ARIADNE_THROW(DomainException,"compose(g,f) with g="<<g<<", f="<<f,"f.codomain()="<<f.codomain()<<" is not a subset of g.domain()="<<g.domain());
        }
        return unchecked_compose(g,f);
    }
    friend ScaledFunctionPatch<M> unchecked_compose(const ScaledFunctionPatch<M>& g, const VectorScaledFunctionPatch<M>& f) {
        return ScaledFunctionPatch<M>(f.domain(),compose(g.model(),unscale(f.models(),g.domain())));
    }


    friend VectorScaledFunctionPatch<M> compose(const VectorMultivariateFunction<P>& g, const VectorScaledFunctionPatch<M>& f) {
        return VectorScaledFunctionPatch<M>(f.domain(),g.evaluate(f.models()));
    }
    friend VectorScaledFunctionPatch<M> compose(const VectorMultivariateFunctionModel<P,PR,PRE>& g, const VectorScaledFunctionPatch<M>& f) {
        return VectorScaledFunctionPatch<M>(f.domain(),g.evaluate(f.models()));
    }
    friend VectorScaledFunctionPatch<M> compose(const VectorScaledFunctionPatch<M>& g, const VectorScaledFunctionPatch<M>& f) {
        if(!subset(f.codomain(),g.domain())) {
            ARIADNE_THROW(DomainException,"compose(g,f) with g="<<g<<", f="<<f,"f.codomain()="<<f.codomain()<<" is not a subset of g.domain()="<<g.domain());
        }
        return unchecked_compose(g,f);
    }
    friend VectorScaledFunctionPatch<M> unchecked_compose(const VectorScaledFunctionPatch<M>& g, const VectorScaledFunctionPatch<M>& f) {
        return VectorScaledFunctionPatch<M>(f.domain(),compose(g.models(),unscale(f.models(),g.domain())));
    }



    friend VectorScaledFunctionPatch<M> derivative(const VectorScaledFunctionPatch<M>& f, SizeType k) {
        ARIADNE_ASSERT_MSG(k<f.argument_size(),"f="<<f<<"\n f.argument_size()="<<f.argument_size()<<" k="<<k);
        VectorScaledFunctionPatch<M> g=f;
        for(SizeType i=0; i!=g.size(); ++i) {
            g[i]=derivative(f[i],k);
        }
        return g;
    }
    friend VectorScaledFunctionPatch<M> antiderivative(const VectorScaledFunctionPatch<M>& f, SizeType k) {
        ARIADNE_ASSERT_MSG(k<f.argument_size(),"f="<<f<<"\n f.argument_size()="<<f.argument_size()<<" k="<<k);
        VectorScaledFunctionPatch<M> g=f;
        for(SizeType i=0; i!=g.size(); ++i) {
            g[i]=antiderivative(f[i],k);
        }
        return g;
    }
    friend VectorScaledFunctionPatch<M> antiderivative(const VectorScaledFunctionPatch<M>& f, SizeType k, NumericType c) {
        ARIADNE_ASSERT_MSG(k<f.argument_size(),"f="<<f<<"\n f.argument_size()="<<f.argument_size()<<" k="<<k);
        VectorScaledFunctionPatch<M> g=f;
        for(SizeType i=0; i!=g.size(); ++i) {
            g[i]=antiderivative(f[i],k,c);
        }
        return g;
    }
    friend VectorScaledFunctionPatch<M> antiderivative(const VectorScaledFunctionPatch<M>& f, SizeType k, GenericNumericType c) {
        return antiderivative(f,k,NumericType(c,f.precision())); 
    }
    friend NormType norm(const VectorScaledFunctionPatch<M>& f) {
        NormType res=norm(f.zero_element());
        for(SizeType i=1; i!=f.result_size(); ++i) {
            res=max(res,norm(f[i]));
        }
        return res;
    }
    NormType distance(const VectorScaledFunctionPatch<M>& f1, const VectorScaledFunctionPatch<M>& f2) {
        return norm(f1-f2);
    }
    NormType distance(const VectorScaledFunctionPatch<M>& f1, const VectorMultivariateFunction<P>& f2) {
        return distance(f1,VectorScaledFunctionPatch<M>(f1.domain(),f2,f1.properties()));
    }

    friend OutputStream& operator<<(OutputStream& os, const VectorScaledFunctionPatch<M>& p) {
        return p.write(os);
    }


    friend Vector< Polynomial<NumericType> > polynomials(const VectorScaledFunctionPatch<M>& tfn) {
        return tfn.polynomials();
    }

};

template<class M> template<class X, EnableIf<CanCall<X,M,Vector<X>>>> Void VectorScaledFunctionPatch<M>::_compute(Vector<X>& r, const Vector<X>& a) const {
    ARIADNE_DEBUG_ASSERT_MSG(r.size()==this->result_size(),"\nr="<<r<<"\nf="<<(*this)<<"\n");
    Vector<X> sa=Ariadne::unscale(a,this->_domain);
    for(SizeType i=0; i!=r.size(); ++i) {
        r[i]=this->_models[i](sa);
    }
}
template<class M> template<class X, DisableIf<CanCall<X,M,Vector<X>>>> Void VectorScaledFunctionPatch<M>::_compute(Vector<X>& r, const Vector<X>& a) const {
    assert(false);
}

template<class M> template<class OP>
VectorScaledFunctionPatch<M> VectorScaledFunctionPatch<M>::_apply(OP op, const VectorScaledFunctionPatch<M>& f1, const VectorScaledFunctionPatch<M>& f2) {
    check_function_patch_domain(to_str(op),f1,f2);
    if(f1.domain()==f2.domain()) {
        return VectorScaledFunctionPatch<M>(f1.domain(),Vector<ModelType>(op(f1.models(),f2.models())));
    } else {
        BoxDomainType new_domain=intersection(f1.domain(),f2.domain());
        return op(restriction(f1,new_domain),restriction(f2,new_domain));
    }
}
template<class M> template<class OP>
VectorScaledFunctionPatch<M> VectorScaledFunctionPatch<M>::_apply(OP op, const ScaledFunctionPatch<M>& f1, const VectorScaledFunctionPatch<M>& f2) {
    check_function_patch_domain(to_str(op),f1,f2);
    if(f1.domain()==f2.domain()) {
        return VectorScaledFunctionPatch<M>(f1.domain(),Vector<ModelType>(op(f1.model(),f2.models())));
    } else {
        BoxDomainType new_domain=intersection(f1.domain(),f2.domain());
        return op(restriction(f1,new_domain),restriction(f2,new_domain));
    }
}
template<class M> template<class OP>
VectorScaledFunctionPatch<M> VectorScaledFunctionPatch<M>::_apply(OP op, const VectorScaledFunctionPatch<M>& f1, const ScaledFunctionPatch<M>& f2) {
    check_function_patch_domain(to_str(op),f1,f2);
    if(f1.domain()==f2.domain()) {
        return VectorScaledFunctionPatch<M>(f1.domain(),Vector<ModelType>(op(f1.models(),f2.model())));
    } else {
        BoxDomainType new_domain=intersection(f1.domain(),f2.domain());
        return op(restriction(f1,new_domain),restriction(f2,new_domain));
    }
}

template<class M> template<class E> VectorScaledFunctionPatch<M>::VectorScaledFunctionPatch(const VectorExpression<E>& ve)
    : _domain(), _models(ve().size(),ve().zero_element().model())
{
    if(ve().size()!=0) { this->_domain=ve().zero_element().domain(); }
    for(SizeType i=0; i!=ve().size(); ++i) { this->set(i,ve()[i]); }
}


//! \brief Test if the function models have the same representation.
template<class M> Bool same(const VectorScaledFunctionPatch<M>& f1, const VectorScaledFunctionPatch<M>& f2) {
    if(!same(f1.domain(),f2.domain())) { return false; }
    for(SizeType i=0; i!=f1.result_size(); ++i) {
        if(!same(f1[i],f2[i])) { return false; }
    }
    return true;
}
//! \brief Tests if a function \a f refines another function \a g.
//! To be a refinement, the domain of \a f must contain the domain of \a g.
template<class M> Bool refines(const VectorScaledFunctionPatch<M>& f1, const VectorScaledFunctionPatch<M>& f2) {
    ARIADNE_ASSERT(f1.result_size()==f2.result_size());
    for(SizeType i=0; i!=f1.result_size(); ++i) {
        if(!refines(f1[i],f2[i])) { return false; }
    }
    return true;
}
template<class M> Bool inconsistent(const VectorScaledFunctionPatch<M>& f1, const VectorScaledFunctionPatch<M>& f2) {
    ARIADNE_ASSERT(f1.result_size()==f2.result_size());
    for(SizeType i=0; i!=f1.result_size(); ++i) {
        if(inconsistent(f1[i],f2[i])) { return true; }
    }
    return false;
}
template<class M> VectorScaledFunctionPatch<M> refinement(const VectorScaledFunctionPatch<M>& f1, const VectorScaledFunctionPatch<M>& f2) {
    ARIADNE_ASSERT(f1.result_size()==f2.result_size());
    VectorScaledFunctionPatch<M> r(f1.result_size());
    for(SizeType i=0; i!=r.result_size(); ++i) {
        r[i]=refinement(f1[i],f2[i]);
    }
    return r;
}

// Sanitised output
template<class M> OutputStream& operator<<(OutputStream& os, const Representation<ScaledFunctionPatch<M>>& repr) {
    return repr.pointer->repr(os); }
template<class M> OutputStream& operator<<(OutputStream& os, const Representation<VectorScaledFunctionPatch<M>>& repr) {
    return repr.pointer->repr(os); }

template<class F> struct ModelRepresentation { const F* pointer; double threshold; };
template<class F> ModelRepresentation<F> model_representation(const F& f, double prpt) {
    ModelRepresentation<F> r={&f,prpt}; return r; }
template<class M> OutputStream& operator<<(OutputStream&,const ModelRepresentation<ScaledFunctionPatch<M>>&);
template<class M> OutputStream& operator<<(OutputStream&,const ModelRepresentation<VectorScaledFunctionPatch<M>>&);

template<class F> struct PolynomialRepresentation { const F* pointer; double threshold; List<String> names; };
template<class F> PolynomialRepresentation<F> polynomial_representation(const F& f, double prpt) {
    PolynomialRepresentation<F> r={&f,prpt}; return r; }
template<class F> PolynomialRepresentation<F> polynomial_representation(const F& f, double prpt, const List<String>& names) {
    PolynomialRepresentation<F> r={&f,prpt,names}; return r; }
template<class M> OutputStream& operator<<(OutputStream&,const PolynomialRepresentation<ScaledFunctionPatch<M>>&);
template<class M> OutputStream& operator<<(OutputStream&,const PolynomialRepresentation<VectorScaledFunctionPatch<M>>&);


template<class M> OutputStream& operator<<(OutputStream& os, const PolynomialRepresentation<VectorScaledFunctionPatch<M>>& repr) {
    const VectorScaledFunctionPatch<M>& function = *repr.pointer;
    os << "[";
    for(SizeType i=0; i!=function.result_size(); ++i) {
        if(i!=0) { os << ","; }
        os << polynomial_representation(function[i],repr.threshold,repr.names);
    }
    return os << "]";
}

template<class M> OutputStream& operator<<(OutputStream& os, const PolynomialRepresentation<List<ScaledFunctionPatch<M>>>& repr) {
    const List<ScaledFunctionPatch<M>>& functions = *repr.pointer;
    os << "[";
    for(SizeType i=0; i!=functions.size(); ++i) {
        if(i!=0) { os << ","; }
        os << polynomial_representation(functions[i],repr.threshold);
    }
    return os << "]";
}

#ifdef ARIADNE_UNDEF

template<class M> Bool check(const Vector<ScaledFunctionPatch<M>>& tv)
{
    for(SizeType i=0; i!=tv.size(); ++i) {
        if(tv.zero_element().domain()!=tv[i].domain()) { return false; }
    }
    return true;
}

template<class M> Vector<Expansion<MultiIndex,typename M::CoefficientType>> expansion(const Vector<ScaledFunctionPatch<M>>& x)
{
    Vector< Expansion<MultiIndex,typename M::CoefficientType> > r(x.size());
    for(SizeType i=0; i!=x.size(); ++i) {
        r[i]=x[i].expansion();
    }
    return r;
}

template<class M> Vector<typename M::ErrorType> error(const Vector<ScaledFunctionPatch<M>>& x)
{
    Vector<typename M::ErrorType> r(x.size());
    for(SizeType i=0; i!=x.size(); ++i) {
        r[i]=x[i].error();
    }
    return r;
}

template<class M> Vector<typename M::CoefficientType> value(const Vector<ScaledFunctionPatch<M>>& x)
{
    Vector<typename M::CoefficientType> r(x.size());
    for(SizeType i=0; i!=x.size(); ++i) {
        r[i]=x[i].value();
    }
    return r;
}

template<class M> Vector<typename M::RangeType> ranges(const Vector<ScaledFunctionPatch<M>>& x)
{
    Vector<typename M::RangeType> r(x.size());
    for(SizeType i=0; i!=x.size(); ++i) {
        r[i]=x[i].range();
    }
    return r;
}

#endif // ARIADNE_UNDEF

template<class M> class VectorScaledFunctionPatchElementReference
    : public DispatchTranscendentalAlgebraOperations<ScaledFunctionPatch<M>, NumericType<M>>
    , public ProvideConcreteGenericArithmeticOperators<ScaledFunctionPatch<M>, ScalarMultivariateFunction<typename M::Paradigm>>
    , public DispatchConcreteGenericAlgebraNumberOperations<ScaledFunctionPatch<M>,NumericType<M>,Number<typename M::Paradigm>>
{
    typedef M ModelType;
    typedef typename M::NumericType NumericType;
    typedef typename M::NormType NormType;
    typedef VectorScaledFunctionPatchElementReference<M> SelfType;
    typedef ScaledFunctionPatch<M> FunctionType;
 public:
    VectorScaledFunctionPatchElementReference(VectorScaledFunctionPatch<M>& c, SizeType i) : _c(&c), _i(i) { }
    operator ScalarScaledFunctionPatch<M> () const { return this->_c->get(this->_i); }
    Void operator=(const VectorScaledFunctionPatchElementReference<M>& x) { this->_c->set(this->_i,x._c->get(x._i)); }
    Void operator=(const ScalarScaledFunctionPatch<M>& x) { this->_c->set(this->_i,x); }
    ScalarScaledFunctionPatch<M> element() const { return this->_c->get(this->_i); }
    BoxDomainType const& domain() const { return this->_c->domain(); }
    const ModelType& model() const { return this->_c->_models[this->_i]; }
    ErrorType error() const { return this->_c->_models[this->_i].error(); }
    Void set_error(const ErrorType& e) { this->_c->_models[this->_i].set_error(e); }
    Void add_error(const ErrorType& e) { this->_c->_models[this->_i].set_error(this->_c->_models[this->_i].error()+e); }
    Void simplify() { this->_c->_models[this->_i].simplify(); }
    template<class X> X operator()(const Vector<X>& x) const { return this->_c->get(this->_i).operator()(x); }
    friend OutputStream& operator<<(OutputStream& os, const VectorScaledFunctionPatchElementReference<M>& f) { return os<<f.element(); }
  public:
    friend FunctionType operator+(SelfType e1, SelfType e2) { return e1.element()+e2.element(); }
    friend FunctionType operator+(NumericType c, SelfType e) { return c+e.element(); }
    friend FunctionType operator+(SelfType e, NumericType c) { return e.element()+c; }
    friend FunctionType operator+(FunctionType f, SelfType e) { return f+e.element(); }
    friend FunctionType operator+(SelfType e, FunctionType f) { return e.element()+f; }
    friend FunctionType operator*(SelfType e1, SelfType e2) { return e1.element()*e2.element(); }
    friend FunctionType operator*(NumericType c, SelfType e) { return c*e.element(); }
    friend FunctionType operator*(SelfType e, NumericType c) { return e.element()*c; }
    friend FunctionType operator*(FunctionType f, SelfType e) { return f*e.element(); }
    friend FunctionType operator*(SelfType e, FunctionType f) { return e.element()*f; }
  private:
    VectorScaledFunctionPatch<M>* _c; SizeType _i;
};


template<class M> class ScaledFunctionPatchFactory
    : public FunctionModelFactoryMixin<ScaledFunctionPatchFactory<M>, ValidatedTag, typename M::PrecisionType, typename M::ErrorPrecisionType>
{
    typedef BoxDomainType D;
    typedef IntervalDomainType SD;

    typedef typename M::Paradigm P;
    typedef typename M::PrecisionType PR;
    typedef typename M::ErrorPrecisionType PRE;
  public:
    typedef P Paradigm;
    typedef PR PrecisionType;
    typedef PRE ErrorPrecisionType;
    typedef typename M::PropertiesType PropertiesType;
    typedef BoxDomainType DomainType;

    explicit ScaledFunctionPatchFactory<M>(PropertiesType properties) : _properties(properties) { }
    PropertiesType properties() const { return this->_properties; }

    CanonicalNumericType<P,PR,PRE> create(const Number<P>& number) const;
    ScalarScaledFunctionPatch<M> create(const BoxDomainType& domain, const ScalarFunctionInterface<P,D>& function) const;
    VectorScaledFunctionPatch<M> create(const BoxDomainType& domain, const VectorFunctionInterface<P,D>& function) const;

    ScaledFunctionPatch<M> create_zero(const DomainType& domain) const;
    ScaledFunctionPatch<M> create_constant(const DomainType& domain, Number<P> const& value) const;
    ScalarFunctionModel<P,D,PR,PRE> create_constant(const BoxDomainType& domain, const CanonicalNumericType<P,PR,PRE>& value) const {
        return create_constant(domain,Number<P>(value)); };
    ScaledFunctionPatch<M> create_coordinate(const DomainType& domain, SizeType index) const;
    VectorScaledFunctionPatch<M> create_zeros(SizeType result_size, const DomainType& domain) const;
    VectorScaledFunctionPatch<M> create_constants(const DomainType& domain, Vector<Number<P>> const& values) const;
    VectorScaledFunctionPatch<M> create_identity(const DomainType& domain) const;
    ScalarScaledFunctionPatch<M> create_identity(const IntervalDomainType& domain) const { return this->create_coordinate(BoxDomainType(1u,domain),0u); };
    CanonicalNumericType<P,PR,PRE> create_number(const Number<P>& number) const { return this->create(number); }
    friend OutputStream& operator<<(OutputStream& os, ScaledFunctionPatchFactory<M> const& factory) {
        return os << "ScaledFunctionPatchFactory( properties=" << factory._properties << " )"; }
  private:
    PropertiesType _properties;
};

template<class M> class ScaledFunctionPatchCreator
    : public FunctionModelCreator<ScaledFunctionPatchFactory<M>,BoxDomainType>
{
  public:
    typedef BoxDomainType DomainType;
    typedef typename M::PropertiesType PropertiesType;
    explicit ScaledFunctionPatchCreator<M>(DomainType domain, PropertiesType properties)
        : FunctionModelCreator<ScaledFunctionPatchFactory<M>,DomainType>(domain,ScaledFunctionPatchFactory<M>(properties)) { }
    PropertiesType properties() const { return this->_factory.properties(); }
};


template<class M> ScaledFunctionPatch<M>* ScaledFunctionPatch<M>::_clone() const {
    return new ScaledFunctionPatch<M>(*this);
}

template<class M> VectorScaledFunctionPatch<M>* VectorScaledFunctionPatch<M>::_clone() const {
    return new VectorScaledFunctionPatch<M>(*this);
}

} // namespace Ariadne

#endif // ARIADNE_FUNCTION_PATCH_HPP
