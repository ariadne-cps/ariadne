/***************************************************************************
 *            function_patch.h
 *
 *  Copyright 2008-15  Pieter Collins
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

/*! \file function_patch.h
 *  \brief Over-approximations of functions on box domains.
 */

#ifndef ARIADNE_FUNCTION_PATCH_H
#define ARIADNE_FUNCTION_PATCH_H

#include <iosfwd>
#include "utility/container.h"
#include "utility/exceptions.h"
#include "utility/declarations.h"
#include "numeric/numeric.h"
#include "algebra/vector.h"
#include "function/taylor_model.h"

#include "algebra/operations.h"
#include "function/function_interface.h"
#include "function/function_mixin.h"
#include "function/function_model.h"
#include "function/function_model_mixin.h"

namespace Ariadne {

template<class X> class Polynomial;

template<class T> using NumericType = typename T::NumericType;
template<class T> using FunctionType = typename T::FunctionType;
template<class T> using GenericType = typename T::GenericType;

template<class M> using ScalarFunctionType = typename M::ScalarFunctionType;
template<class M> using VectorFunctionType = typename M::VectorFunctionType;
template<class M> using ScalarFunctionModelType = ScalarFunctionModel<Paradigm<M>,PrecisionType<M>,ErrorPrecisionType<M>>;
template<class M> using VectorFunctionModelType = VectorFunctionModel<Paradigm<M>,PrecisionType<M>,ErrorPrecisionType<M>>;

template<class M> class FunctionPatch;
template<class M> using ScalarFunctionPatch = FunctionPatch<M>;
template<class M> class VectorFunctionPatch;
template<class M> class VectorFunctionPatchElementReference;

inline Float64Approximation convert_error_to_bounds(const PositiveFloat64Approximation& e) { return Float64Approximation(0.0); }
inline Float64Bounds convert_error_to_bounds(const PositiveFloat64UpperBound& e) { return Float64Bounds(-e.raw(),+e.raw()); }
inline Float64Bounds convert_error_to_bounds(const Float64Error& e) { return Float64Bounds(-e.raw(),+e.raw()); }

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

template<class M> class FunctionPatchFactory;
template<class M> class FunctionPatchCreator;


/*! \ingroup FunctionModelSubModule
 *  \brief A ScalarTaylorFunction is a type of FunctionModel in which a the restriction of a scalar function \f$f:\R^n\rightarrow\R\f$ on a domain \f$D\f$ is approximated by polynomial \f$p\f$ with uniform error \f$e\f$.
 *
 * Formally, a ScalarTaylorFunction is a triple \f$(D,p,e)\f$ representing a set of continuous functions \f$\mathrm{T}(D,p,e)\f$ by
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
 * \sa Expansion, TaylorModel, VectorTaylorFunction, TaylorConstrainedImageSet.
 */
template<class M> class FunctionPatch
    : public ScalarFunctionModelMixin<FunctionPatch<M>, typename M::Paradigm, typename M::PrecisionType, typename M::ErrorPrecisionType>
    , public DispatchSymbolicAlgebraOperations<FunctionPatch<M>, NumericType<M>>
    , public ProvideConcreteGenericArithmeticOperators<FunctionPatch<M>, ScalarFunction<typename M::Paradigm>>
    , public DispatchConcreteGenericAlgebraNumberOperations<FunctionPatch<M>,NumericType<M>,Number<typename M::Paradigm>>
{
    typedef typename M::Paradigm P;
    typedef typename M::RawFloatType F;
    typedef typename M::PrecisionType PR;
    typedef typename M::ErrorPrecisionType PRE;
  public:
    typedef BoxDomain DomainType;
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
    typedef FunctionPatch<M> FunctionPatchType;
    typedef ScalarFunction<Paradigm> FunctionType;
    typedef ScalarFunction<Paradigm> GenericType;
    typedef Number<Paradigm> GenericNumericType;
    typedef typename M::PropertiesType PropertiesType;
  private:
    static const CoefficientType _zero;
    DomainType _domain;
    ModelType _model;
  public:

    //@{
    //! \name Constructors and destructors.
    //! \brief Default constructor.
    explicit FunctionPatch();
    //! \brief Construct a FunctionPatch<M> over the domain \a d.
    //explicit FunctionPatch(const DomainType& d);
    explicit FunctionPatch(const DomainType& d, PropertiesType prp);
    //! \brief Construct a FunctionPatch<M> over the domain \a d, based on the scaled model \a m.
    explicit FunctionPatch(const DomainType& d, const ModelType& m);

    explicit FunctionPatch(const ExactBoxType& d, const Expansion<FloatValue<PR>>& p, const FloatError<PR>& e, const Sweeper<RawFloat<PR>>& prp);
    explicit FunctionPatch(const ExactBoxType& d, const Expansion<RawFloat<PR>>& p, const RawFloat<PR>& e, const Sweeper<RawFloat<PR>>& prp);

    explicit FunctionPatch(const ScalarFunctionModelType<M>& f);
    FunctionPatch& operator=(const ScalarFunctionModelType<M>& f);

    //! \brief Construct a FunctionPatch over the domain \a d from the function \a f.
    explicit FunctionPatch(const DomainType& d, const ScalarFunctionType<M>& f, PropertiesType prp);
    //@}

    //@{
    //! \name Assignment to constant values.
    //! \brief Set equal to an interval constant, keeping the same number of arguments.
    FunctionPatch<M>& operator=(const NumericType& c) { this->_model=c; return *this; }
    //@}

    //@{
    //! \name Named constructors.
    //! \brief Construct a zero function over domain \a d.
    static FunctionPatch<M> zero(const DomainType& d, PropertiesType prp);
    //! \brief Construct a constant quantity in \a as independent variables.
    static FunctionPatch<M> constant(const DomainType& d, const NumericType& c, PropertiesType prp);
    //! \brief Construct the quantity \f$x_j\f$ over the domain \a d.
    static FunctionPatch<M> coordinate(const DomainType& d, SizeType j, PropertiesType prp);
    //! \brief Construct a constant quantity in \a as independent variables with value zero and uniform error \a 1
    static FunctionPatch<M> unit_ball(const DomainType& d, PropertiesType prp);
    //! \brief Construct the quantity \f$x_j\f$ over the domain \a d.
    static VectorFunctionPatch<M> identity(const DomainType& d, PropertiesType prp);

    //! \brief Construct the quantity \f$c+\sum g_jx_j\f$ over the domain \a d. // DEPRECATED
    static FunctionPatch<M> affine(const DomainType& d, const CoefficientType& c, const Vector<CoefficientType>& g, PropertiesType prp);
    //! \brief Construct the quantity \f$c+\sum g_jx_j \pm e\f$ over domain \a d. // DEPRECATED
    static FunctionPatch<M> affine(const DomainType& d, const CoefficientType& x, const Vector<CoefficientType>& g, const ErrorType& e, PropertiesType prp) ;

    //! \brief Return the vector of constants with interval values \a c over domain \a d.
    static Vector<FunctionPatch<M>> constants(const DomainType& d, const Vector<NumericType>& c, PropertiesType prp);
    //! \brief Return the vector of variables with values \a x over domain \a d.
    static Vector<FunctionPatch<M>> coordinates(const DomainType& d, PropertiesType prp);
    //! \brief Return the vector of variables in the range \a imin to \a imax with values \a x over domain \a d.
    static Vector<FunctionPatch<M>> coordinates(const DomainType& d, SizeType imin, SizeType imax, PropertiesType prp);
    //@}

    //@{
    //! \name Prototype constructors.
    friend FunctionPatchCreator<M> factory(FunctionPatch<M>const& f) {
        return FunctionPatchCreator<M>(f.domain(),f.properties()); }
    //! \brief Construct a zero function over the same domain with the same computational properties.
    FunctionPatch<M> create_zero() const;
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
        FloatBounds<PR> radius(this->_domain[i].radius(),this->model().precision()); return cast_exact(this->_model.gradient_value(i)/radius); }

    //! \brief A polynomial representation.
    Polynomial<FloatBounds<PR>> polynomial() const;
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
    Bool operator==(const FunctionPatch<M>& tv) const;
    //! \brief Inequality operator.
    Bool operator!=(const FunctionPatch<M>& tv) const { return !(*this==tv); }
    //@}

    //@{
    //! \name Function operations.
    //! \brief An over-approximation to the range of the function.
    RangeType range() const { return this->_model.range(); }
    //! \brief Evaluate the function at the point \a x.
    FloatBounds<PR> operator()(const Vector<FloatBounds<PR>>& x) const;
    FloatBounds<PR> operator()(const Vector<FloatValue<PR>>& x) const;
    FloatApproximation<PR> operator()(const Vector<FloatApproximation<PR>>& x) const;

    //! \brief Compute an approximation to gradient derivative of the function at the point \a x.
    Covector<NumericType> gradient(const Vector<NumericType>& x) const;
    //@}

    //@{
    //! \name Simplification operations.
   //! \brief Remove all terms whose coefficient has magnitude
    //! lower than the cutoff threshold of the quantity.
    FunctionPatch<M>& simplify() { this->_model.simplify(); return *this; }
    //! \brief Remove all terms whose degree is higher than \a deg or
    //! whose coefficient has magnitude less than \a eps.
    FunctionPatch<M>& simplify(const PropertiesType& prp) { this->_model.simplify(prp); return *this; }
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
    friend OutputStream& operator<<(OutputStream& os, const FunctionPatch<M>& x) {
        return x.write(os); }
    //@}

  public:
    Void clobber() { this->_model.clobber(); }
  private:
    friend class TaylorFunctionFactory;
    friend class ScalarFunctionMixin<FunctionPatch<M>, P>;
    friend class ScalarFunctionModelMixin<FunctionPatch<M>, P>;
  public:
    template<class X, EnableIf<CanCall<X,M,Vector<X>>> =dummy> Void _compute(X& r, const Vector<X>& a) const;
    template<class X, DisableIf<CanCall<X,M,Vector<X>>> =dummy> Void _compute(X& r, const Vector<X>& a) const;
  public:
    template<class OP> static FunctionPatch<M> _create(OP op, FunctionPatch<M> const& fp);
    template<class OP> static FunctionPatch<M> _create(OP op, FunctionPatch<M> const& fp1, FunctionPatch<M> const& fp2);
    template<class OP> static FunctionPatch<M> _create(OP op, FunctionPatch<M> const& fp1, NumericType const& c2);
    template<class OP> static FunctionPatch<M> _create(OP op, NumericType const& c1, FunctionPatch<M> const& fp2);
    static FunctionPatch<M> _create(Pow op, FunctionPatch<M> const& fp, Int n);
  private:
    FunctionPatch<M>* _derivative(SizeType j) const;
    FunctionPatch<M>* _clone() const;
    FunctionPatch<M>* _create() const;
    virtual FunctionPatchFactory<M>* _factory() const;
  public:
    VectorFunctionPatch<M> join(const FunctionPatch<M>& f1, const FunctionPatch<M>& f2);
    VectorFunctionPatch<M> combine(const FunctionPatch<M>& f1, const FunctionPatch<M>& f2);

    //! \brief Restrict to a subdomain.
    friend FunctionPatch<M> restriction(const FunctionPatch<M>& fp, const ExactBoxType& dom) {
        if(not(subset(dom,fp.domain()))) { ARIADNE_THROW(DomainException,"restiction(FunctionPatch<M>,ExactBoxType)","fp="<<fp<<", dom="<<dom); }
        return unchecked_compose(fp,FunctionPatch<M>::identity(dom,fp.properties()));
    }
    //! \brief Restrict a component to a subdomain
    // To scale from a model on [a,b] to a model on [c,d], use scale factor s=(d-c)/(b-a)
    // and translation t=((c+d)-(a+b))/(b-a)
    // Because we are scaling the model on [-1,+1], this is not the same as
    // the mapping taking [a,b] to [c,d]
    friend FunctionPatch<M> partial_restriction(const FunctionPatch<M>& fp, SizeType k, const ExactIntervalType& ivl) {
        ExactBoxType dom=fp.domain(); dom[k]=ivl; return restriction(fp,dom);
    }
    //! \brief Extend over a larger domain. Only possible if the larger domain is only larger where the smaller domain is a singleton.
    //! The extension is performed keeping \a x constant over the new coordinates. // DEPRECATED
    friend FunctionPatch<M> extension(const FunctionPatch<M>& fp, const ExactBoxType& dom) {
        return unchecked_compose(fp,FunctionPatch<M>::identity(dom,fp.properties()));
    }
    //! \brief Pre-compose with a projection from the Cartesian product of the given domains.
    friend FunctionPatch<M> embed(const ExactBoxType& dom1, const FunctionPatch<M>& tv2,const ExactBoxType& dom3) {
        return FunctionPatch<M>(product(dom1,tv2.domain(),dom3),embed(dom1.size(),tv2.model(),dom3.size())); }
    //! \brief Split the domain into two.
    friend Pair<FunctionPatch<M>,FunctionPatch<M>> split(const FunctionPatch<M>& tv, SizeType j) {
        Pair<ModelType,ModelType> models={split(tv.model(),j,SplitPart::LOWER),split(tv.model(),j,SplitPart::LOWER)};
        Pair<ExactBoxType,ExactBoxType> subdomains=split(tv.domain(),j);
        return make_pair(FunctionPatch<M>(subdomains.first,models.first),
                         FunctionPatch<M>(subdomains.second,models.second));
    }

    template<class OP> static FunctionPatch<M> _apply(OP op, FunctionPatch<M> const& f);
    template<class OP> static FunctionPatch<M> _apply(OP op, FunctionPatch<M> const& f1, FunctionPatch<M> const& f2);

/*
    friend FunctionPatch<M>& operator+=(FunctionPatch<M>& f, const NumericType& c) { f.model()+=c; return f; }
    friend FunctionPatch<M>& operator*=(FunctionPatch<M>& f, const NumericType& c) { f.model()*=c; return f; }
    friend FunctionPatch<M> operator-(const FunctionPatch<M>& f) { return apply(Neg(),f); }
    friend FunctionPatch<M> operator+(const FunctionPatch<M>& f1, const FunctionPatch<M>& f2) { return apply(Add(),f1,f2); }
    friend FunctionPatch<M> operator-(const FunctionPatch<M>& f1, const FunctionPatch<M>& f2) { return apply(Sub(),f1,f2); }
    friend FunctionPatch<M> operator*(const FunctionPatch<M>& f1, const FunctionPatch<M>& f2) { return apply(Mul(),f1,f2); }
    friend FunctionPatch<M> operator/(const FunctionPatch<M>& f1, const FunctionPatch<M>& f2) { return apply(Div(),f1,f2); }
*/
    friend FunctionPatch<M> max(const FunctionPatch<M>& f1, const FunctionPatch<M>& f2) { return _apply(Max(),f1,f2); }
    friend FunctionPatch<M> min(const FunctionPatch<M>& f1, const FunctionPatch<M>& f2) { return _apply(Min(),f1,f2); }
    friend FunctionPatch<M> abs(const FunctionPatch<M>& f) { return _apply(Abs(),f); }

    friend FunctionPatch<M> antiderivative(const FunctionPatch<M>& f, SizeType k) {
        return FunctionPatch<M>(f.domain(),antiderivative(f.model(),k)*rad(f.domain()[k])); }
    friend FunctionPatch<M> antiderivative(const FunctionPatch<M>& f, SizeType k, const NumericType& c) {
        ARIADNE_ASSERT(k<f.argument_size());
        ARIADNE_ASSERT(decide(contains(f.domain()[k],c)));
        FunctionPatch<M> g = antiderivative(f,k);
        VectorFunctionPatch<M> h ( VectorFunctionPatch<M>::identity(f.domain(),f.properties()) );
        h[k] = FunctionPatch<M>::constant(f.domain(),c,f.properties());
        return g-compose(g,h);
    }
    friend FunctionPatch<M> derivative(const FunctionPatch<M>& f, SizeType k) {
        return FunctionPatch<M>(f.domain(),derivative(f.model(),k)*rec(rad(f.domain()[k]))); }


    friend FunctionPatch<M> partial_evaluate(const FunctionPatch<M>& te, SizeType k, const NumericType& c) {
        // Scale c to domain
        const SizeType as=te.argument_size();
        ARIADNE_ASSERT(k<as);
        const ExactBoxType& domain=te.domain();
        const ExactIntervalType& dk=domain[k];
        NumericType sc=(c-med(dk))/rad(dk);

        ExactBoxType new_domain(as-1);
        for(SizeType i=0; i!=k; ++i) { new_domain[i]=domain[i]; }
        for(SizeType i=k; i!=as-1; ++i) { new_domain[i]=domain[i+1]; }

        M new_model=partial_evaluate(te.model(),k,sc);

        return FunctionPatch<M>(new_domain,new_model);
    }
    friend NumericType evaluate(const FunctionPatch<M>& f, const Vector<NumericType>& x) {
        if(!definitely(contains(f.domain(),x))) {
            ARIADNE_THROW(DomainException,"evaluate(f,x) with f="<<f<<", x="<<x,"x is not an element of f.domain()="<<f.domain());
        }
        return unchecked_evaluate(f,x);
    }
    friend NumericType unchecked_evaluate(const FunctionPatch<M>& f, const Vector<NumericType>& x) {
        return evaluate(f.model(),unscale(x,f.domain()));
    }


    friend NormType norm(const FunctionPatch<M>& f) {
        return norm(f.model());
    }
    friend NormType distance(const FunctionPatch<M>& f1, const FunctionPatch<M>& f2) {
        return norm(f1-f2);
    }
    friend NormType distance(const FunctionPatch<M>& f1, const ScalarFunction<P>& f2) {
        return distance(f1,FunctionPatch<M>(f1.domain(),f2,f1.properties()));
    }


    friend Polynomial<NumericType> polynomial(const FunctionPatch<M>& tfn) {
        return tfn.polynomial();
    }

};

template<class M> template<class X, EnableIf<CanCall<X,M,Vector<X>>>> Void FunctionPatch<M>::_compute(X& r, const Vector<X>& a) const {
    r = this->_model(unscale(a,this->_domain));
}
template<class M> template<class X, DisableIf<CanCall<X,M,Vector<X>>>> Void FunctionPatch<M>::_compute(X& r, const Vector<X>& a) const {
    assert(false);
}

template<class FP1, class FP2> Void check_function_patch_domain(String const& op_str, const FP1& fp1, const FP2& fp2) {
    ARIADNE_ASSERT_MSG(!is_empty(intersection(fp1.domain(),fp2.domain())),
                    op_str<<"((Vector)FunctionPatch<M> fp1, (Vector)FunctionPatch<M> fp2) with fp1="<<fp1<<" fp2="<<fp2<<
                    ": domains are disjoint");
}
template<class M> template<class OP> FunctionPatch<M> FunctionPatch<M>::_apply(OP op, FunctionPatch<M> const& f) {
    return FunctionPatch<M>(f.domain(),op(f.model()));
}
template<class M> template<class OP> FunctionPatch<M> FunctionPatch<M>::_apply(OP op, FunctionPatch<M> const& f1, FunctionPatch<M> const& f2) {
    check_function_patch_domain(to_str(op),f1,f2);
    if(f1.domain()==f2.domain()) {
        return FunctionPatch<M>(f1.domain(),op(f1.model(),f2.model()));
    } else {
        ExactBoxType domain=intersection(f1.domain(),f2.domain());
        return FunctionPatch<M>(domain,op(restriction(f1,domain).model(),restriction(f2,domain).model()));
    }
}


template<class M> template<class OP> FunctionPatch<M> FunctionPatch<M>::_create(OP op, FunctionPatch<M> const& fp1, FunctionPatch<M> const& fp2) {
    assert(fp1.domain()==fp2.domain()); return FunctionPatch<M>(fp1.domain(),op(fp1.model(),fp2.model()));
}
template<class M> template<class OP> FunctionPatch<M> FunctionPatch<M>::_create(OP op, FunctionPatch<M> const& fp) {
    return FunctionPatch<M>(fp.domain(),op(fp.model()));
}
template<class M> template<class OP> FunctionPatch<M> FunctionPatch<M>::_create(OP op, FunctionPatch<M> const& fp1, typename M::NumericType const& c2) {
    return FunctionPatch<M>(fp1.domain(),op(fp1.model(),c2));
}
template<class M> template<class OP> FunctionPatch<M> FunctionPatch<M>::_create(OP op, typename M::NumericType const& c1, FunctionPatch<M> const& fp2) {
    return FunctionPatch<M>(fp2.domain(),op(c1,fp2.model()));
}
template<class M> FunctionPatch<M> FunctionPatch<M>::_create(Pow op, FunctionPatch<M> const& fp, Int n) {
    return FunctionPatch<M>(fp.domain(),op(fp.model(),n));
}

//! \brief Test if the function models have the same representation.
template<class M> Bool same(const FunctionPatch<M>& tv1, const FunctionPatch<M>& tv2) {
    return tv1.domain()==tv2.domain() && same(tv1.model(),tv2.model());
}
//! \brief Test if the quantity is a better approximation than \a t throughout the domain.
template<class M> Bool refines(const FunctionPatch<M>& tv1, const FunctionPatch<M>& tv2) {
    if(tv1.domain()==tv2.domain()) { return refines(tv1.model(),tv2.model()); }
    if(subset(tv2.domain(),tv1.domain())) { return refines(restriction(tv1,tv2.domain()).model(),tv2.model()); }
    else { return false; }
}
//! \brief Test if the function models are inconsistent with representing the same exact function.
template<class M> Bool inconsistent(const FunctionPatch<M>& tv1, const FunctionPatch<M>& tv2) {
    if(tv1.domain()==tv2.domain()) {
        return inconsistent(tv1.model(),tv2.model());
    } else {
        ExactBoxType domain=intersection(tv1.domain(),tv2.domain());
        return inconsistent(restriction(tv1,domain).model(),restriction(tv2,domain).model());
    }
}
//! \brief Compute an over-approximation to the common refinement of \a x1 and \a x2.
template<class M> FunctionPatch<M> refinement(const FunctionPatch<M>& tv1, const FunctionPatch<M>& tv2) {
    ARIADNE_ASSERT(tv1.domain()==tv2.domain());
    return FunctionPatch<M>(tv1.domain(),refinement(tv1.model(),tv2.model()));
}
//! \brief Remove the error term.
template<class M> FunctionPatch<M> midpoint(const FunctionPatch<M>& f) {
    M tm=f.model();
    tm.set_error(0u);
    return FunctionPatch<M>(f.domain(),tm);
}








/*! \ingroup FunctionModelSubModule
 *  \brief A Taylor function model with multivalued codomain built from the TaylorModel class.
 *
 *  See also TaylorModel, FunctionPatch<M>, VectorTaylorFunction.
 */
template<class M> class VectorFunctionPatch
    : public VectorFunctionModelMixin<VectorFunctionPatch<M>,typename M::Paradigm,typename M::PrecisionType,typename M::ErrorPrecisionType>
{
    friend class VectorFunctionPatchElementReference<M>;
    typedef typename M::Paradigm P;
    typedef typename M::RawFloatType F;
    typedef typename M::PrecisionType PR;
    typedef typename M::ErrorPrecisionType PRE;
  public:
    typedef ExactBoxType DomainType;
    typedef M ModelType;
    typedef Box<typename ModelType::CodomainType> CodomainType;
    typedef Box<typename ModelType::RangeType> RangeType;
    typedef typename ModelType::ExpansionType ExpansionType;
    typedef typename ModelType::CoefficientType CoefficientType;
    typedef typename ModelType::ErrorType ErrorType;
    typedef typename ModelType::NumericType NumericType;
    typedef typename ModelType::NormType NormType;
    typedef typename M::PropertiesType PropertiesType;
    typedef P Paradigm;
    typedef PR PrecisionType;
    typedef PRE ErrorPrecisionType;

    typedef FloatApproximation<PR> ApproximateNumericType;
    typedef FloatBounds<PR> ValidatedNumericType;
    typedef FloatValue<PR> ExactNumericType;


    //! \brief Default constructor constructs a Taylor model of order zero with no arguments and no result variables.
    VectorFunctionPatch<M>();

    //! \brief Construct the zero vector function over an unspecified domain.
    explicit VectorFunctionPatch<M>(SizeType result_size);

    //! \brief Construct from a result size and a domain.
    VectorFunctionPatch<M>(SizeType result_size, const ExactBoxType& domain, PropertiesType properties);

    //! \brief Construct a vector function all of whose components are the same.
    VectorFunctionPatch<M>(SizeType result_size, const FunctionPatch<M>& scalar_function);

    //! \brief Construct from a domain and the expansion.
    VectorFunctionPatch<M>(const ExactBoxType& domain,
                         const Vector<Expansion<FloatValue<PR>>>& expansion,
                         PropertiesType properties);

    //! \brief Construct from a domain, and expansion and errors.
    VectorFunctionPatch<M>(const ExactBoxType& domain,
                         const Vector<Expansion<FloatValue<PR>>>& expansion,
                         const Vector<FloatError<PR>>& error,
                         PropertiesType properties);

    //! \brief Construct from a domain, and expansion and errors.
    VectorFunctionPatch<M>(const ExactBoxType& domain,
                         const Vector<Expansion<RawFloat<PR>>>& expansion,
                         const Vector<RawFloat<PR>>& error,
                         PropertiesType properties);

    //! \brief Construct from a domain, and expansion and errors.
    VectorFunctionPatch<M>(const ExactBoxType& domain,
                         const Vector<Expansion<RawFloat<PR>>>& expansion,
                         PropertiesType properties);

    //! \brief Construct from a domain and the models.
    explicit VectorFunctionPatch<M>(const ExactBoxType& domain, const Vector< ModelType >& variables);

    //! \brief Construct from a \a domain, a \a function, and \a properties determining the accuracy.
    VectorFunctionPatch<M>(const ExactBoxType& domain,
                         const VectorFunctionType<M>& function,
                         const PropertiesType& properties);

    //! \brief Construct from a vector of scalar Taylor functions.
    explicit VectorFunctionPatch<M>(const Vector<FunctionPatch<M>>& components);

    //! \brief Construct from a list of scalar Taylor functions.
    explicit VectorFunctionPatch<M>(const List<FunctionPatch<M>>& components);

    //! \brief Construct from an initializer list of scalar Taylor functions.
    VectorFunctionPatch<M>(InitializerList<FunctionPatch<M>> components);

    //! \brief Construct from a vector expression.
    template<class E> explicit VectorFunctionPatch<M>(const VectorExpression<E>& ve);

    explicit VectorFunctionPatch<M> (const VectorFunctionModelType<M>& f);
    VectorFunctionPatch<M>& operator=(const VectorFunctionModelType<M>& f);

    //! \brief Equality operator.
    Bool operator==(const VectorFunctionPatch<M>& p) const;
    //! \brief Inequality operator.
    Bool operator!=(const VectorFunctionPatch<M>& p) const;

    // Data access
    //! \brief The properties used to control approximation of the function model.
    PropertiesType properties() const;
    //! \brief Set the properties used to control approximation of the function model.
    Void set_properties(PropertiesType prp);
    //! \brief The data used to define the domain of the Taylor model.
    const ExactBoxType domain() const;
    //! \brief A rough bound for the range of the function.
    const ExactBoxType codomain() const;
    //! \brief The centre of the Taylor model.
    const Vector<CoefficientType> centre() const;
    //! \brief The range of the Taylor model.
    const RangeType range() const;
    //! \brief The data used to define the Taylor models.
    const Vector<ModelType>& models() const;
    Vector<ModelType>& models();
    //! \brief The data used to define the centre of the Taylor models.
    const Vector<Expansion<FloatValue<PR>>> expansions() const;

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
    FunctionPatch<M> zero_element() const;
    //! \brief Get the \a ith Taylor variable
    FunctionPatch<M> get(SizeType i) const;
    //! \brief Set the \a ith Taylor variable
    Void set(SizeType i, const FunctionPatch<M>& te);
    //! \brief The \a ith Taylor variable
    FunctionPatch<M> const operator[](SizeType i) const;
    //! \brief The \a ith Taylor variable
    VectorFunctionPatchElementReference<M> operator[](SizeType i);


    //! \brief Evaluate the Taylor model at the point \a x.
    Vector<ValidatedNumericType> operator()(const Vector<ValidatedNumericType>& x) const;
    Vector<ApproximateNumericType> operator()(const Vector<ApproximateNumericType>& x) const;
    Vector<ValidatedNumericType> operator()(const Vector<ExactNumericType>& x) const;
    //! \brief Compute an approximation to Jacobian derivative of the Taylor model sat the point \a x.
    Matrix<NumericType> jacobian(const Vector<NumericType>& x) const;

    //! \brief Simplify the representation.
    VectorFunctionPatch<M>& simplify();
    //! \brief Simplify the representation as specified by \a properties.
    VectorFunctionPatch<M>& simplify(const PropertiesType& properties);
    //! \brief Set the error to zero.
    Void clobber();

    friend FunctionPatchCreator<M> factory(VectorFunctionPatch<M>const& f) {
        return FunctionPatchCreator<M>(f.domain(),f.properties()); }

    //! \brief The constant Taylor model with range \a r and argument domain \a d.
    static VectorFunctionPatch<M> constant(const ExactBoxType& d, const Vector<NumericType>& r, PropertiesType prp);
    //! \brief The identity Taylor model on domain \a d.
    static VectorFunctionPatch<M> identity(const ExactBoxType& d, PropertiesType prp);
    //! \brief Return the vector of variables in the range with values \a x over domain \a d.
    static VectorFunctionPatch<M> projection(const ExactBoxType& d, SizeType imin, SizeType imax, PropertiesType prp);

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
    VectorFunctionPatch<M>& truncate(const MultiIndexBound& bd);
    //! \brief Restrict to a subdomain.
    Void restrict(const ExactBoxType& d);
    //! \brief Adjoin a scalar function.
    Void adjoin(const FunctionPatch<M>& sf);

    //! \brief Write to an output stream.
    OutputStream& write(OutputStream& os) const;

    //! \brief Write a full representation to an output stream.
    OutputStream& repr(OutputStream& os) const;

  private:
    Void _compute_jacobian() const;
    Void _set_argument_size(SizeType n);
    SizeType _compute_maximum_component_size() const;
    Void _resize(SizeType rs, SizeType as, ushort d, ushort s);
    virtual ScalarFunctionPatch<M>* _get(SizeType i) const { return new FunctionPatch<M>(this->_domain,this->_models[i]); }
    virtual VectorFunctionPatch<M>* _clone() const;
    virtual VectorFunctionPatch<M>* _create() const;
    virtual FunctionPatchFactory<M>* _factory() const;
  private:
    friend class VectorFunctionMixin<VectorFunctionPatch<M>,P>;
    friend class TaylorFunctionFactory;
  public:
    template<class X, EnableIf<CanCall<X,M,Vector<X>>> =dummy> Void _compute(Vector<X>& r, const Vector<X>& a) const;
    template<class X, DisableIf<CanCall<X,M,Vector<X>>> =dummy> Void _compute(Vector<X>& r, const Vector<X>& a) const;
  private:
    /* Domain of definition. */
    ExactBoxType _domain;
    Vector< ModelType > _models;

  public:
    //! \brief Compute the function \f$(f \oplus g)(x)=(f(x),g(x))\f$.
    friend VectorFunctionPatch<M> join(const VectorFunctionPatch<M>& f1, const FunctionPatch<M>& f2) {
        ARIADNE_ASSERT_MSG(f1.domain()==f2.domain(),"f1="<<f1<<", f2="<<f2);
        return VectorFunctionPatch<M>(f1.domain(),join(f1.models(),f2.model()));
    }
    friend VectorFunctionPatch<M> join(const VectorFunctionPatch<M>& f, const VectorFunctionPatch<M>& g) {
        ARIADNE_ASSERT(f.domain()==g.domain());
        return VectorFunctionPatch<M>(f.domain(),join(f.models(),g.models()));
    }
    friend VectorFunctionPatch<M> join(const FunctionPatch<M>& f1, const FunctionPatch<M>& f2) {
        ARIADNE_ASSERT(f1.domain()==f2.domain());
        return VectorFunctionPatch<M>(f1.domain(),{f1.model(),f2.model()});
    }
    friend VectorFunctionPatch<M> join(const FunctionPatch<M>& f1, const VectorFunctionPatch<M>& f2) {
        ARIADNE_ASSERT(f1.domain()==f2.domain());
        return VectorFunctionPatch<M>(f1.domain(),join(f1.model(),f2.models()));
    }

    //! \brief Compute the function \f$(f\otimes g)(x,y)=(f(x),g(y))\f$.
    friend VectorFunctionPatch<M> combine(const FunctionPatch<M>& f1, const FunctionPatch<M>& f2) {
        return VectorFunctionPatch<M>(product(f1.domain(),f2.domain()),combine(Vector<M>(1u,f1.model()),{f2.model()}));
    }
    friend VectorFunctionPatch<M> combine(const FunctionPatch<M>& f1, const VectorFunctionPatch<M>& f2) {
        return VectorFunctionPatch<M>(product(f1.domain(),f2.domain()),combine({f1.model()},f2.models()));
    }
    friend VectorFunctionPatch<M> combine(const VectorFunctionPatch<M>& f1, const FunctionPatch<M>& f2) {
        return VectorFunctionPatch<M>(product(f1.domain(),f2.domain()),combine(f1.models(),{f2.model()}));
    }
    friend VectorFunctionPatch<M> combine(const VectorFunctionPatch<M>& f1, const VectorFunctionPatch<M>& f2) {
        return VectorFunctionPatch<M>(product(f1.domain(),f2.domain()),combine(f1.models(),f2.models()));
    }

    friend VectorFunctionPatch<M> embed(const VectorFunctionPatch<M>& f, const ExactIntervalType& d) {
        return embed(ExactBoxType(),f,ExactBoxType(1u,d));
    }
    friend VectorFunctionPatch<M> embed(const VectorFunctionPatch<M>& f, const ExactBoxType& d) {
        return embed(ExactBoxType(),f,d);
    }
    friend VectorFunctionPatch<M> embed(const ExactBoxType& d, const VectorFunctionPatch<M>& f) {
        return embed(d,f,ExactBoxType());
    }
    friend VectorFunctionPatch<M> embed(const ExactBoxType& d1, const VectorFunctionPatch<M>& f, const ExactBoxType& d2) {
        return VectorFunctionPatch<M>(product(d1,f.domain(),d2),embed(d1.size(),f.models(),d2.size()));
    }

    friend VectorFunctionPatch<M> restriction(const VectorFunctionPatch<M>& f, const ExactBoxType& d) {
        ARIADNE_ASSERT_MSG(subset(d,f.domain()),"Cannot restriction "<<f<<" to non-sub-domain "<<d);
        if(d==f.domain()) { return f; }
        VectorFunctionPatch<M> r(f.result_size(),d,f.properties());
        for(SizeType i=0; i!=r.result_size(); ++i) {
            r.set(i,restriction(f[i],d));
        }
        return r;
    }
    friend VectorFunctionPatch<M> partial_restriction(const VectorFunctionPatch<M>& tf, SizeType k, const ExactIntervalType& d) {
        VectorFunctionPatch<M> r(tf.result_size(),tf.domain(),tf.properties());
        for(SizeType i=0; i!=tf.result_size(); ++i) {
            r[i]=partial_restriction(tf[i],k,d);
        }
        return r;
    }
    friend VectorFunctionPatch<M> restriction(const VectorFunctionPatch<M>& tf, SizeType k, const ExactIntervalType& d) {
        return partial_restriction(tf,k,d);
    }
    friend Pair<VectorFunctionPatch<M>,VectorFunctionPatch<M>> split(const VectorFunctionPatch<M>& tf, SizeType j) {
        typedef M ModelType;
        Pair<Vector<ModelType>,Vector<ModelType>> models=split(tf.models(),j);
        Pair<ExactBoxType,ExactBoxType> subdomains=split(tf.domain(),j);
        return make_pair(VectorFunctionPatch<M>(subdomains.first,models.first),
                        VectorFunctionPatch<M>(subdomains.second,models.second));

    }

    friend VectorFunctionPatch<M>& operator+=(VectorFunctionPatch<M>& f, const VectorFunctionPatch<M>& g) {
        ARIADNE_ASSERT(f.result_size()==g.result_size());
        ARIADNE_ASSERT(subset(f.domain(),g.domain()));
        ARIADNE_ASSERT(f.domain()==g.domain());
        f.models()+=g.models();
        return f;
    }
    friend VectorFunctionPatch<M>& operator-=(VectorFunctionPatch<M>& f, const VectorFunctionPatch<M>& g) {
        ARIADNE_ASSERT(f.result_size()==g.result_size());
        ARIADNE_ASSERT(subset(f.domain(),g.domain()));
        ARIADNE_ASSERT(f.domain()==g.domain());
        f.models()+=g.models();
        return f;
    }
    friend VectorFunctionPatch<M>& operator+=(VectorFunctionPatch<M>& f, const Vector<NumericType>& c) {
        ARIADNE_ASSERT(f.result_size()==c.size());
        f.models()+=c;
        return f;
    }
    friend VectorFunctionPatch<M>& operator-=(VectorFunctionPatch<M>& f, const Vector<NumericType>& c) {
        ARIADNE_ASSERT(f.result_size()==c.size());
        f.models()-=c;
        return f;
    }
    friend VectorFunctionPatch<M>& operator*=(VectorFunctionPatch<M>& f, const NumericType& c) {
        f.models()*=c;
        return f;
    }
    friend VectorFunctionPatch<M>& operator/=(VectorFunctionPatch<M>& f, const NumericType& c) {
        f.models()/=c;
        return f;
    }

    template<class OP> static VectorFunctionPatch<M> _apply(OP op, const VectorFunctionPatch<M>& f1, const VectorFunctionPatch<M>& f2);
    template<class OP> static VectorFunctionPatch<M> _apply(OP op, const FunctionPatch<M>& f1, const VectorFunctionPatch<M>& f2);
    friend VectorFunctionPatch<M> operator+(const VectorFunctionPatch<M>& f1, const VectorFunctionPatch<M>& f2) {
        return _apply(Plus(),f1,f2); }
    friend VectorFunctionPatch<M> operator-(const VectorFunctionPatch<M>& f1, const VectorFunctionPatch<M>& f2) {
        return _apply(Minus(),f1,f2); }
    friend VectorFunctionPatch<M> operator*(const FunctionPatch<M>& f1, const VectorFunctionPatch<M>& f2) {
        return _apply(Times(),f1,f2); }
    friend VectorFunctionPatch<M> operator*(const VectorFunctionPatch<M>& f1, const FunctionPatch<M>& f2) {
        return _apply(Times(),f1,f2); }
    friend VectorFunctionPatch<M> operator/(const VectorFunctionPatch<M>& f1, const FunctionPatch<M>& f2) {
        return _apply(Divides(),f1,f2); }


    friend VectorFunctionPatch<M> operator-(const VectorFunctionPatch<M>& f) {
        return VectorFunctionPatch<M>(f.domain(),Vector<M>(-f.models())); }
    friend VectorFunctionPatch<M> operator*(const NumericType& c, const VectorFunctionPatch<M>& f) {
        return VectorFunctionPatch<M>(f.domain(),Vector<M>(f.models()*c)); }
    friend VectorFunctionPatch<M> operator*(const VectorFunctionPatch<M>& f, const NumericType& c) {
        return VectorFunctionPatch<M>(f.domain(),Vector<M>(f.models()*c)); }
    friend VectorFunctionPatch<M> operator/(const VectorFunctionPatch<M>& f, const NumericType& c) {
        return VectorFunctionPatch<M>(f.domain(),Vector<M>(f.models()/c)); }
    friend VectorFunctionPatch<M> operator+(const VectorFunctionPatch<M>& f, const Vector<NumericType>& c) {
        return VectorFunctionPatch<M>(f.domain(),Vector<M>(f.models()+c)); }
    friend VectorFunctionPatch<M> operator-(const VectorFunctionPatch<M>& f, const Vector<NumericType>& c) {
        return VectorFunctionPatch<M>(f.domain(),Vector<M>(f.models()-c)); }

    friend VectorFunctionPatch<M> operator*(const Matrix<NumericType>& A, const VectorFunctionPatch<M>& f) {
        ARIADNE_PRECONDITION(A.column_size()==f.size());
        return VectorFunctionPatch<M>(f.domain(),A*f.models());
    }
    friend VectorFunctionPatch<M> operator+(const VectorFunction<P>& f1, const VectorFunctionPatch<M>& tf2) {
        return VectorFunctionPatch<M>(tf2.domain(),f1,tf2.properties())+tf2; }
    friend VectorFunctionPatch<M> operator-(const VectorFunction<P>& f1, const VectorFunctionPatch<M>& tf2) {
        return VectorFunctionPatch<M>(tf2.domain(),f1,tf2.properties())-tf2; }
    friend VectorFunctionPatch<M> operator*(const ScalarFunction<P>& f1, const VectorFunctionPatch<M>& tf2) {
        return FunctionPatch<M>(tf2.domain(),f1,tf2.properties())*tf2; }
    friend VectorFunctionPatch<M> operator*(const VectorFunction<P>& f1, const FunctionPatch<M>& tf2) {
        return VectorFunctionPatch<M>(tf2.domain(),f1,tf2.properties())*tf2; }
    friend VectorFunctionPatch<M> operator/(const VectorFunction<P>& f1, const FunctionPatch<M>& tf2) {
        return VectorFunctionPatch<M>(tf2.domain(),f1,tf2.properties())/tf2; }
    friend VectorFunctionPatch<M> operator+(const VectorFunctionPatch<M>& tf1, const VectorFunction<P>& f2) {
        return tf1+VectorFunctionPatch<M>(tf1.domain(),f2,tf1.properties()); }
    friend VectorFunctionPatch<M> operator-(const VectorFunctionPatch<M>& tf1, const VectorFunction<P>& f2) {
        return tf1-VectorFunctionPatch<M>(tf1.domain(),f2,tf1.properties()); }
    friend VectorFunctionPatch<M> operator*(const FunctionPatch<M>& tf1, const VectorFunction<P>& f2) {
        return tf1*VectorFunctionPatch<M>(tf1.domain(),f2,tf1.properties()); }
    friend VectorFunctionPatch<M> operator*(const VectorFunctionPatch<M>& tf1, const ScalarFunction<P>& f2) {
        return tf1*FunctionPatch<M>(tf1.domain(),f2,tf1.properties()); }
    friend VectorFunctionPatch<M> operator/(const VectorFunctionPatch<M>& tf1, const ScalarFunction<P>& f2) {
        return tf1/FunctionPatch<M>(tf1.domain(),f2,tf1.properties()); }



    friend VectorFunctionPatch<M> partial_evaluate(const VectorFunctionPatch<M>& tf, SizeType k, const NumericType& c) {
        // Scale c to domain
        const SizeType as=tf.argument_size();
        ARIADNE_ASSERT(k<as);
        const Vector<ExactIntervalType>& domain=tf.domain();
        const ExactIntervalType& dk=domain[k];
        NumericType sc=(c-med(dk))/rad(dk);

        Vector<ExactIntervalType> new_domain(as-1);
        for(SizeType i=0; i!=k; ++i) { new_domain[i]=domain[i]; }
        for(SizeType i=k; i!=as-1; ++i) { new_domain[i]=domain[i+1]; }

        Vector<M> new_models=partial_evaluate(tf.models(),k,sc);

        return VectorFunctionPatch<M>(new_domain,new_models);
    }
    friend Vector<NumericType> evaluate(const VectorFunctionPatch<M>& f, const Vector<NumericType>& x) {
        if(!definitely(contains(f.domain(),x))) {
            ARIADNE_THROW(DomainException,"evaluate(f,x) with f="<<f<<", x="<<x,"x is not a subset of f.domain()="<<f.domain());
        }
        return unchecked_evaluate(f,x);
    }
    friend Vector<NumericType> unchecked_evaluate(const VectorFunctionPatch<M>& f, const Vector<NumericType>& x) {
        return evaluate(f.models(),unscale(x,f.domain()));
    }

    friend FunctionPatch<M> compose(const ScalarFunction<P>& g, const VectorFunctionPatch<M>& f) {
        return FunctionPatch<M>(f.domain(),g.evaluate(f.models()));
    }
    friend FunctionPatch<M> compose(const FunctionPatch<M>& g, const VectorFunctionPatch<M>& f) {
        if(!subset(f.codomain(),g.domain())) {
            ARIADNE_THROW(DomainException,"compose(g,f) with g="<<g<<", f="<<f,"f.codomain()="<<f.codomain()<<" is not a subset of g.domain()="<<g.domain());
        }
        return unchecked_compose(g,f);
    }
    friend FunctionPatch<M> unchecked_compose(const FunctionPatch<M>& g, const VectorFunctionPatch<M>& f) {
        return FunctionPatch<M>(f.domain(),compose(g.model(),unscale(f.models(),g.domain())));
    }


    friend VectorFunctionPatch<M> compose(const VectorFunction<P>& g, const VectorFunctionPatch<M>& f) {
        return VectorFunctionPatch<M>(f.domain(),g.evaluate(f.models()));
    }
    friend VectorFunctionPatch<M> compose(const VectorFunctionPatch<M>& g, const VectorFunctionPatch<M>& f) {
        if(!subset(f.codomain(),g.domain())) {
            ARIADNE_THROW(DomainException,"compose(g,f) with g="<<g<<", f="<<f,"f.codomain()="<<f.codomain()<<" is not a subset of g.domain()="<<g.domain());
        }
        return unchecked_compose(g,f);
    }
    friend VectorFunctionPatch<M> unchecked_compose(const VectorFunctionPatch<M>& g, const VectorFunctionPatch<M>& f) {
        return VectorFunctionPatch<M>(f.domain(),compose(g.models(),unscale(f.models(),g.domain())));
    }



    friend VectorFunctionPatch<M> derivative(const VectorFunctionPatch<M>& f, SizeType k) {
        ARIADNE_ASSERT_MSG(k<f.argument_size(),"f="<<f<<"\n f.argument_size()="<<f.argument_size()<<" k="<<k);
        VectorFunctionPatch<M> g=f;
        for(SizeType i=0; i!=g.size(); ++i) {
            g[i]=derivative(f[i],k);
        }
        return g;
    }
    friend VectorFunctionPatch<M> antiderivative(const VectorFunctionPatch<M>& f, SizeType k) {
        ARIADNE_ASSERT_MSG(k<f.argument_size(),"f="<<f<<"\n f.argument_size()="<<f.argument_size()<<" k="<<k);
        VectorFunctionPatch<M> g=f;
        for(SizeType i=0; i!=g.size(); ++i) {
            g[i]=antiderivative(f[i],k);
        }
        return g;
    }
    friend VectorFunctionPatch<M> antiderivative(const VectorFunctionPatch<M>& f, SizeType k, NumericType c) {
        ARIADNE_ASSERT_MSG(k<f.argument_size(),"f="<<f<<"\n f.argument_size()="<<f.argument_size()<<" k="<<k);
        VectorFunctionPatch<M> g=f;
        for(SizeType i=0; i!=g.size(); ++i) {
            g[i]=antiderivative(f[i],k,c);
        }
        return g;
    }

    friend NormType norm(const VectorFunctionPatch<M>& f) {
        NormType res=norm(f.zero_element());;
        for(SizeType i=1; i!=f.result_size(); ++i) {
            res=max(res,norm(f[i]));
        }
        return res;
    }
    NormType distance(const VectorFunctionPatch<M>& f1, const VectorFunctionPatch<M>& f2) {
        return norm(f1-f2);
    }
    NormType distance(const VectorFunctionPatch<M>& f1, const VectorFunction<P>& f2) {
        return distance(f1,VectorFunctionPatch<M>(f1.domain(),f2,f1.properties()));
    }

    friend OutputStream& operator<<(OutputStream& os, const VectorFunctionPatch<M>& p) {
        return p.write(os);
    }


    friend Vector< Polynomial<NumericType> > polynomials(const VectorFunctionPatch<M>& tfn) {
        return tfn.polynomials();
    }

};

template<class M> template<class X, EnableIf<CanCall<X,M,Vector<X>>>> Void VectorFunctionPatch<M>::_compute(Vector<X>& r, const Vector<X>& a) const {
    ARIADNE_DEBUG_ASSERT_MSG(r.size()==this->result_size(),"\nr="<<r<<"\nf="<<(*this)<<"\n");
    Vector<X> sa=Ariadne::unscale(a,this->_domain);
    for(SizeType i=0; i!=r.size(); ++i) {
        r[i]=this->_models[i](sa);
    }
}
template<class M> template<class X, DisableIf<CanCall<X,M,Vector<X>>>> Void VectorFunctionPatch<M>::_compute(Vector<X>& r, const Vector<X>& a) const {
    assert(false);
}

template<class M> template<class OP>
VectorFunctionPatch<M> VectorFunctionPatch<M>::_apply(OP op, const VectorFunctionPatch<M>& f1, const VectorFunctionPatch<M>& f2) {
    check_function_patch_domain(to_str(op),f1,f2);
    if(f1.domain()==f2.domain()) {
        return VectorFunctionPatch<M>(f1.domain(),Vector<ModelType>(op(f1.models(),f2.models())));
    } else {
        ExactBoxType new_domain=intersection(f1.domain(),f2.domain());
        return op(restriction(f1,new_domain),restriction(f2,new_domain));
    }
}
template<class M> template<class OP>
VectorFunctionPatch<M> VectorFunctionPatch<M>::_apply(OP op, const FunctionPatch<M>& f1, const VectorFunctionPatch<M>& f2) {
    check_function_patch_domain(to_str(op),f1,f2);
    if(f1.domain()==f2.domain()) {
        return VectorFunctionPatch<M>(f1.domain(),Vector<ModelType>(op(f1.models(),f2.models())));
    } else {
        ExactBoxType new_domain=intersection(f1.domain(),f2.domain());
        return op(restriction(f1,new_domain),restriction(f2,new_domain));
    }
}

template<class M> template<class E> VectorFunctionPatch<M>::VectorFunctionPatch(const VectorExpression<E>& ve)
    : _domain(), _models(ve().size(),ve().zero_element().model())
{
    if(ve().size()!=0) { this->_domain=ve().zero_element().domain(); }
    for(SizeType i=0; i!=ve().size(); ++i) { this->set(i,ve()[i]); }
}


//! \brief Test if the function models have the same representation.
template<class M> Bool same(const VectorFunctionPatch<M>& f1, const VectorFunctionPatch<M>& f2) {
    if(!same(f1.domain(),f2.domain())) { return false; }
    for(SizeType i=0; i!=f1.result_size(); ++i) {
        if(!same(f1[i],f2[i])) { return false; }
    }
    return true;
}
//! \brief Tests if a function \a f refines another function \a g.
//! To be a refinement, the domain of \a f must contain the domain of \a g.
template<class M> Bool refines(const VectorFunctionPatch<M>& f1, const VectorFunctionPatch<M>& f2) {
    ARIADNE_ASSERT(f1.result_size()==f2.result_size());
    for(SizeType i=0; i!=f1.result_size(); ++i) {
        if(!refines(f1[i],f2[i])) { return false; }
    }
    return true;
}
template<class M> Bool inconsistent(const VectorFunctionPatch<M>& f1, const VectorFunctionPatch<M>& f2) {
    ARIADNE_ASSERT(f1.result_size()==f2.result_size());
    for(SizeType i=0; i!=f1.result_size(); ++i) {
        if(inconsistent(f1[i],f2[i])) { return true; }
    }
    return false;
}
template<class M> VectorFunctionPatch<M> refinement(const VectorFunctionPatch<M>& f1, const VectorFunctionPatch<M>& f2) {
    ARIADNE_ASSERT(f1.result_size()==f2.result_size());
    VectorFunctionPatch<M> r(f1.result_size());
    for(SizeType i=0; i!=r.result_size(); ++i) {
        r[i]=refinement(f1[i],f2[i]);
    }
    return r;
}

// Sanitised output
template<class M> OutputStream& operator<<(OutputStream& os, const Representation<FunctionPatch<M>>& repr) {
    return repr.pointer->repr(os); }
template<class M> OutputStream& operator<<(OutputStream& os, const Representation<VectorFunctionPatch<M>>& repr) {
    return repr.pointer->repr(os); }

template<class F> struct ModelRepresentation { const F* pointer; double threshold; };
template<class F> ModelRepresentation<F> model_representation(const F& f, double prpt) {
    ModelRepresentation<F> r={&f,prpt}; return r; }
template<class M> OutputStream& operator<<(OutputStream&,const ModelRepresentation<FunctionPatch<M>>&);
template<class M> OutputStream& operator<<(OutputStream&,const ModelRepresentation<VectorFunctionPatch<M>>&);

template<class F> struct PolynomialRepresentation { const F* pointer; double threshold; List<String> names; };
template<class F> PolynomialRepresentation<F> polynomial_representation(const F& f, double prpt) {
    PolynomialRepresentation<F> r={&f,prpt}; return r; }
template<class F> PolynomialRepresentation<F> polynomial_representation(const F& f, double prpt, const List<String>& names) {
    PolynomialRepresentation<F> r={&f,prpt,names}; return r; }
template<class M> OutputStream& operator<<(OutputStream&,const PolynomialRepresentation<FunctionPatch<M>>&);
template<class M> OutputStream& operator<<(OutputStream&,const PolynomialRepresentation<VectorFunctionPatch<M>>&);


template<class M> OutputStream& operator<<(OutputStream& os, const PolynomialRepresentation<VectorFunctionPatch<M>>& repr) {
    const VectorFunctionPatch<M>& function = *repr.pointer;
    os << "[";
    for(SizeType i=0; i!=function.result_size(); ++i) {
        if(i!=0) { os << ","; }
        os << polynomial_representation(function[i],repr.threshold,repr.names);
    }
    return os << "]";
}

template<class M> OutputStream& operator<<(OutputStream& os, const PolynomialRepresentation<List<FunctionPatch<M>>>& repr) {
    const List<FunctionPatch<M>>& functions = *repr.pointer;
    os << "[";
    for(SizeType i=0; i!=functions.size(); ++i) {
        if(i!=0) { os << ","; }
        os << polynomial_representation(functions[i],repr.threshold);
    }
    return os << "]";
}

#ifdef ARIADNE_UNDEF

template<class M> Bool check(const Vector<FunctionPatch<M>>& tv)
{
    for(SizeType i=0; i!=tv.size(); ++i) {
        if(tv.zero_element().domain()!=tv[i].domain()) { return false; }
    }
    return true;
}

template<class M> Vector<Expansion<typename M::CoefficientType>> expansion(const Vector<FunctionPatch<M>>& x)
{
    Vector< Expansion<typename M::CoefficientType> > r(x.size());
    for(SizeType i=0; i!=x.size(); ++i) {
        r[i]=x[i].expansion();
    }
    return r;
}

template<class M> Vector<typename M::ErrorType> error(const Vector<FunctionPatch<M>>& x)
{
    Vector<typename M::ErrorType> r(x.size());
    for(SizeType i=0; i!=x.size(); ++i) {
        r[i]=x[i].error();
    }
    return r;
}

template<class M> Vector<typename M::CoefficientType> value(const Vector<FunctionPatch<M>>& x)
{
    Vector<typename M::CoefficientType> r(x.size());
    for(SizeType i=0; i!=x.size(); ++i) {
        r[i]=x[i].value();
    }
    return r;
}

template<class M> Vector<typename M::RangeType> ranges(const Vector<FunctionPatch<M>>& x)
{
    Vector<typename M::RangeType> r(x.size());
    for(SizeType i=0; i!=x.size(); ++i) {
        r[i]=x[i].range();
    }
    return r;
}

#endif // ARIADNE_UNDEF

template<class M> class VectorFunctionPatchElementReference
    : public DispatchSymbolicAlgebraOperations<FunctionPatch<M>, NumericType<M>>
    , public ProvideConcreteGenericArithmeticOperators<FunctionPatch<M>, ScalarFunction<typename M::Paradigm>>
    , public DispatchConcreteGenericAlgebraNumberOperations<FunctionPatch<M>,NumericType<M>,Number<typename M::Paradigm>>
{
    typedef M ModelType;
    typedef typename M::NumericType NumericType;
    typedef typename M::NormType NormType;
    typedef VectorFunctionPatchElementReference<M> SelfType;
    typedef FunctionPatch<M> FunctionType;
 public:
    VectorFunctionPatchElementReference(VectorFunctionPatch<M>& c, SizeType i) : _c(&c), _i(i) { }
    operator ScalarFunctionPatch<M> () const { return this->_c->get(this->_i); }
    Void operator=(const VectorFunctionPatchElementReference<M>& x) { this->_c->set(this->_i,x._c->get(x._i)); }
    Void operator=(const ScalarFunctionPatch<M>& x) { this->_c->set(this->_i,x); }
    ScalarFunctionPatch<M> element() const { return this->_c->get(this->_i); }
    ExactBoxType const& domain() const { return this->_c->domain(); }
    const ModelType& model() const { return this->_c->_models[this->_i]; }
    ErrorType error() const { return this->_c->_models[this->_i].error(); }
    Void set_error(const ErrorType& e) { this->_c->_models[this->_i].set_error(e); }
    Void simplify() { this->_c->_models[this->_i].simplify(); }
    template<class X> X operator()(const Vector<X>& x) const { return this->_c->get(this->_i).operator()(x); }
    friend OutputStream& operator<<(OutputStream& os, const VectorFunctionPatchElementReference<M>& f) { return os<<f.element(); }
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
    VectorFunctionPatch<M>* _c; SizeType _i;
};


template<class M> class FunctionPatchFactory
    : public FunctionModelFactoryMixin<FunctionPatchFactory<M>, ValidatedTag, typename M::PrecisionType, typename M::ErrorPrecisionType>
{
    typedef typename M::Paradigm P;
    typedef typename M::PrecisionType PR;
    typedef typename M::ErrorPrecisionType PRE;
  public:
    typedef P Paradigm;
    typedef PR PrecisionType;
    typedef PRE ErrorPrecisionType;
    typedef typename M::PropertiesType PropertiesType;
    typedef BoxDomain DomainType;

    explicit FunctionPatchFactory<M>(PropertiesType properties) : _properties(properties) { }
    PropertiesType properties() const { return this->_properties; }

    CanonicalNumericType<P,PR,PRE> create(const Number<P>& number) const;
    ScalarFunctionPatch<M> create(const ExactBoxType& domain, const ScalarFunctionInterface<P>& function) const;
    VectorFunctionPatch<M> create(const ExactBoxType& domain, const VectorFunctionInterface<P>& function) const;

    FunctionPatch<M> create_zero(const DomainType& domain) const;
    FunctionPatch<M> create_constant(const DomainType& domain, Number<P> const& value) const;
    ScalarFunctionModel<P,PR,PRE> create_constant(const ExactBoxType& domain, const CanonicalNumericType<P,PR,PRE>& value) const {
        return create_constant(domain,Number<P>(value)); };
    FunctionPatch<M> create_coordinate(const DomainType& domain, SizeType index) const;
    VectorFunctionPatch<M> create_zeros(SizeType result_size, const DomainType& domain) const;
    VectorFunctionPatch<M> create_constants(const DomainType& domain, Vector<Number<P>> const& values) const;
    VectorFunctionPatch<M> create_identity(const DomainType& domain) const;
    ScalarFunctionPatch<M> create_identity(const ExactIntervalType& domain) const { return this->create_coordinate(ExactBoxType(1u,domain),0u); };
    CanonicalNumericType<P,PR,PRE> create_number(const Number<P>& number) const { return this->create(number); }
    friend OutputStream& operator<<(OutputStream& os, FunctionPatchFactory<M> const& factory) {
        return os << "FunctionPatchFactory( properties=" << factory._properties << " )"; }
  private:
    PropertiesType _properties;
};

template<class M> class FunctionPatchCreator
    : public FunctionModelCreator<FunctionPatchFactory<M>>
{
  public:
    typedef BoxDomain DomainType;
    typedef typename M::PropertiesType PropertiesType;
    explicit FunctionPatchCreator<M>(DomainType domain, PropertiesType properties)
        : FunctionModelCreator<FunctionPatchFactory<M>>(domain,FunctionPatchFactory<M>(properties)) { }
    PropertiesType properties() const { return this->_factory.properties(); }
};



} // namespace Ariadne

#endif // ARIADNE_FUNCTION_PATCH_H
