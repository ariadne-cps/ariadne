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

#include "algebra/algebra_operations.h"
#include "function/function_interface.h"
#include "function/function_mixin.h"
#include "function/function_model.h"

namespace Ariadne {

template<class X> class Polynomial;

template<class T> using NumericType = typename T::NumericType;
template<class T> using FunctionType = typename T::FunctionType;
template<class T> using GenericType = typename T::GenericType;
template<class T> using ScalarFunctionType = typename T::ScalarFunctionType;
template<class T> using VectorFunctionType = typename T::VectorFunctionType;

template<class M> class FunctionPatch;
template<class M> using ScalarFunctionPatch = FunctionPatch<M>;
template<class M> class VectorFunctionPatch;
template<class M> class VectorFunctionPatchElementReference;

// Declare basic algebra operators
template<class M> struct IsAlgebra<FunctionPatch<M>> : True { };
template<class M> struct IsSymbolicAlgebra<FunctionPatch<M>> : True { };

// Declare basic algebra operators
template<class A> A operator+(const A& ca1, const GenericType<A>& ga2) { return ca1+ca1.create(ga2); }
template<class A> A operator-(const A& ca1, const GenericType<A>& ga2) { return ca1-ca1.create(ga2); }
template<class A> A operator*(const A& ca1, const GenericType<A>& ga2) { return ca1*ca1.create(ga2); }
template<class A> A operator/(const A& ca1, const GenericType<A>& ga2) { return ca1/ca1.create(ga2); }
template<class A> A operator+(const GenericType<A>& ga1, const A& ca2) { return ca2.create(ga1)+ca2; }
template<class A> A operator-(const GenericType<A>& ga1, const A& ca2) { return ca2.create(ga1)+ca2; }
template<class A> A operator*(const GenericType<A>& ga1, const A& ca2) { return ca2.create(ga1)+ca2; }
template<class A> A operator/(const GenericType<A>& ga1, const A& ca2) { return ca2.create(ga1)+ca2; }

inline ApproximateFloat64 convert_error_to_bounds(const PositiveApproximateFloat64& e) { return ApproximateFloat64(0.0); }
inline ValidatedFloat64 convert_error_to_bounds(const PositiveUpperFloat64& e) { return ValidatedFloat64(-e.raw(),+e.raw()); }



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
    : public ScalarFunctionModelMixin<FunctionPatch<M>, typename M::Paradigm>
{
    typedef typename M::Paradigm P;
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
    typedef ScalarFunction<Paradigm> FunctionType;
    typedef ScalarFunction<Paradigm> GenericType;
  private:
    static const CoefficientType _zero;
    DomainType _domain;
    ModelType _model;
  public:

    //@{
    /*! \name Constructors and destructors. */
    //! \brief Default constructor.
    explicit FunctionPatch();
    //! \brief Construct a FunctionPatch<M> over the domain \a d.
    //explicit FunctionPatch(const DomainType& d);
    explicit FunctionPatch(const DomainType& d, Sweeper swp);
    //! \brief Construct a FunctionPatch<M> over the domain \a d, based on the scaled model \a m.
    explicit FunctionPatch(const DomainType& d, const ModelType& m);

    explicit FunctionPatch(const ExactBox& d, const Expansion<ExactFloat64>& p, const ErrorFloat64& e, const Sweeper& swp);
    explicit FunctionPatch(const ExactBox& d, const Expansion<RawFloat64>& p, const RawFloat64& e, const Sweeper& swp);

    explicit FunctionPatch(const ScalarFunctionModel<ValidatedTag>& f);
    FunctionPatch& operator=(const ScalarFunctionModel<ValidatedTag>& f);

    //! \brief Construct a FunctionPatch over the domain \a d from the function \a f.
    explicit FunctionPatch(const DomainType& d, const ScalarFunctionType<M>& f, Sweeper swp);
    //@}

    //@{
    /*! \name Assignment to constant values. */
    //! \brief Set equal to an interval constant, keeping the same number of arguments.
    FunctionPatch<M>& operator=(const NumericType& c) { this->_model=c; return *this; }
    //@}

    //@{
    /*! \name Named constructors. */
    //! \brief Construct a zero function over domain \a d.
    static FunctionPatch<M> zero(const DomainType& d, Sweeper swp);
    //! \brief Construct a constant quantity in \a as independent variables.
    static FunctionPatch<M> constant(const DomainType& d, const NumericType& c, Sweeper swp);
    //! \brief Construct the quantity \f$x_j\f$ over the domain \a d.
    static FunctionPatch<M> coordinate(const DomainType& d, SizeType j, Sweeper swp);
    //! \brief Construct the quantity \f$x_j\f$ over the domain \a d.
    static VectorFunctionPatch<M> identity(const DomainType& d, Sweeper swp);

    //! \brief Construct the quantity \f$c+\sum g_jx_j\f$ over the domain \a d. // DEPRECATED
    static FunctionPatch<M> affine(const DomainType& d, const CoefficientType& c, const Vector<CoefficientType>& g, Sweeper swp);
    //! \brief Construct the quantity \f$c+\sum g_jx_j \pm e\f$ over domain \a d. // DEPRECATED
    static FunctionPatch<M> affine(const DomainType& d, const CoefficientType& x, const Vector<CoefficientType>& g, const ErrorType& e, Sweeper swp) ;

    //! \brief Return the vector of constants with interval values \a c over domain \a d.
    static Vector<FunctionPatch<M>> constants(const DomainType& d, const Vector<NumericType>& c, Sweeper swp);
    //! \brief Return the vector of variables with values \a x over domain \a d.
    static Vector<FunctionPatch<M>> coordinates(const DomainType& d, Sweeper swp);
    //! \brief Return the vector of variables in the range \a imin to \a imax with values \a x over domain \a d.
    static Vector<FunctionPatch<M>> coordinates(const DomainType& d, SizeType imin, SizeType imax, Sweeper swp);
    //@}

    //@{
    /*! \name Prototype constructors. */
    FunctionPatch<M> create_zero() const;
    FunctionPatch<M> create_constant(NumericType const& c) const;
    FunctionPatch<M> create_coordinate(SizeType j) const;
    FunctionPatch<M> create(GenericType const& f) const;
    //@}

    //@{
    /*! \name Data access */
    //! \brief The domain of the quantity.
    const DomainType domain() const { return this->_domain; }
    //! \brief The scaled expansion over a unit box with error bound.
    const ModelType& model() const { return this->_model; }
    //! \brief A reference to the scaled expansion over a unit box with error bound.
    ModelType& model() { return this->_model; }

    //! \brief An over-approximation to the range of the quantity; not necessarily tight.
    const CodomainType codomain() const { this->_model.codomain(); }
    //! \brief The scaled expansion over a unit box.
    const ExpansionType& expansion() const { return this->_model.expansion(); }
    //! \brief The error of the expansion over the domain.
    const ErrorType& error() const { return this->_model.error(); }
    /*! \brief The accuracy parameter used to control approximation of the Taylor function. */
    Sweeper sweeper() const { return this->_model.sweeper(); }
    //! \brief A reference to the expansion.
    ExpansionType& expansion() { return this->_model.expansion(); }
    //! \brief A reference to the error of the expansion over the domain.
    ErrorType& error() { return this->_model.error(); }

    //! \brief The constant term in the expansion.
    const CoefficientType& value() const { return this->_model.value(); }
    //! \brief The gradient at the centre of the domain.
    const CoefficientType gradient_value(SizeType i) const { return cast_exact(this->_model.gradient_value(i)/this->_domain[i].radius()); }

    //! \brief A polynomial representation.
    Polynomial<ValidatedFloat64> polynomial() const;
    //! \brief A multivalued function equal to the model on the domain.
    ScalarFunctionType<M> function() const;

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
    /*! \name Comparison operators. */
    //! \brief Equality operator. Tests equality of representation, including error term.
    Bool operator==(const FunctionPatch<M>& tv) const;
    //! \brief Inequality operator.
    Bool operator!=(const FunctionPatch<M>& tv) const { return !(*this==tv); }
    //@}

    //@{
    /*! \name Function operations. */
    //! \brief An over-approximation to the range of the function.
    UpperInterval range() const { return this->_model.range(); }
    //! \brief Evaluate the function at the point \a x.
    ValidatedNumber operator()(const Vector<ValidatedNumber>& x) const;
    ValidatedNumber operator()(const Vector<ExactNumber>& x) const;
    ApproximateNumber operator()(const Vector<ApproximateNumber>& x) const;

    /*! \brief Compute an approximation to gradient derivative of the function at the point \a x. */
    Covector<NumericType> gradient(const Vector<NumericType>& x) const;
    //@}

    //@{
    /*! \name Simplification operations. */
   //! \brief Remove all terms whose coefficient has magnitude
    //! lower than the cutoff threshold of the quantity.
    FunctionPatch<M>& sweep() { this->_model.sweep(); return *this; }
    //! \brief Remove all terms whose degree is higher than \a deg or
    //! whose coefficient has magnitude less than \a eps.
    FunctionPatch<M>& sweep(const SweeperInterface& swp) { this->_model.sweep(swp); return *this; }
    //@}

    //@{
    /*! \name Accuracy parameters. */
    //! \copydoc TaylorModel<ValidatedFloat64>::set_sweeper()
    Void set_sweeper(const Sweeper& swp) { this->_model.set_sweeper(swp); }
    //@}

    //@{
    /*! \name Non-arithmetic operations. */
    //! \brief Restrict to a subdomain.
    Void restrict(const DomainType& d);


    //@}

    /*! \name Arithmetic operations. */
    //@}

    //@{
    /*! \name Stream input/output operators. */
    //! \brief Write to an output stream.
    OutputStream& write(OutputStream& os) const;
    /*! \brief Write a full representation to an output stream. */
    OutputStream& repr(OutputStream& os) const;
    //! \brief Write to an output stream.
    friend OutputStream& operator<<(OutputStream& os, const FunctionPatch<M>& x) {
        return x.write(os); }
    //@}

  public:
    Void clobber() { this->_model.clobber(); }
  private:
    friend class TaylorFunctionFactory;
    friend class ScalarFunctionMixin<FunctionPatch<M>, ValidatedTag>;
    friend class ScalarFunctionModelMixin<FunctionPatch<M>, ValidatedTag>;
  public:
    template<class T> Void _compute(T& r, const Vector<T>& a) const {
        typedef typename T::NumericType R;
        r=Ariadne::horner_evaluate(this->_model.expansion(),Ariadne::unscale(a,this->_domain))
            + convert_error_to_bounds(this->_model.error());
    }
  private:
    FunctionPatch<M>* _derivative(SizeType j) const;
    FunctionPatch<M>* _clone() const;
    FunctionPatch<M>* _create() const;
    VectorFunctionModelInterface<Paradigm>* _create_identity() const;
    VectorFunctionModelInterface<Paradigm>* _create_vector(SizeType i) const;
};

//! \brief Restrict to a subdomain.
template<class M> FunctionPatch<M> restriction(const FunctionPatch<M>& x, const ExactBox& d);
//! \brief Extend over a larger domain. Only possible if the larger domain is only larger where the smaller domain is a singleton.
//! The extension is performed keeping \a x constant over the new coordinates. // DEPRECATED
template<class M> FunctionPatch<M> extension(const FunctionPatch<M>& x, const ExactBox& d);

//! \brief Test if the quantity is a better approximation than \a t throughout the domain.
template<class M> Bool refines(const FunctionPatch<M>& x1, const FunctionPatch<M>& x2);
//! \brief Test if the function models are inconsistent with representing the same exact function.
template<class M> Bool inconsistent(const FunctionPatch<M>& x1, const FunctionPatch<M>& x2);
//! \brief Compute an over-approximation to the common refinement of \a x1 and \a x2.
template<class M> FunctionPatch<M> refinement(const FunctionPatch<M>& x1, const FunctionPatch<M>& x2);


// Normed algebra
template<class M> NormType norm(const FunctionPatch<M>& f);

// Differential algebra
template<class M> FunctionPatch<M> antiderivative(const FunctionPatch<M>& x, SizeType k, ExactFloat64 c);
template<class M> FunctionPatch<M> antiderivative(const FunctionPatch<M>& x, SizeType k);
template<class M> FunctionPatch<M> derivative(const FunctionPatch<M>& x, SizeType k);

// Function algebra
template<class M> FunctionPatch<M> embed(const ExactBox& dom1, const FunctionPatch<M>& tv2,const ExactBox& dom3);
// Set the value of the \a kth variable to c
template<class M> FunctionPatch<M> partial_evaluate(const FunctionPatch<M>& f, SizeType k, const NumericType<M>& c);
// Evaluate a scalar Taylor function on a vector.
template<class M> NumericType<M> unchecked_evaluate(const FunctionPatch<M>&, const Vector<NumericType<M>>&);
template<class M> NumericType<M> evaluate(const FunctionPatch<M>&, const Vector<NumericType<M>>&);

// Compose with an function.
template<class M> FunctionPatch<M> compose(const FunctionType<M>& x, const VectorFunctionPatch<M>& y);
template<class M> FunctionPatch<M> unchecked_compose(const FunctionPatch<M>&, const VectorFunctionPatch<M>&);

// Split the variable over two domains, subdividing along the independent variable j.
template<class M> Pair<FunctionPatch<M>,FunctionPatch<M>> split(const FunctionPatch<M>& x, SizeType j);


template<class M> FunctionPatch<M>& operator+=(FunctionPatch<M>& f, const NumericType<M>& c);
template<class M> FunctionPatch<M>& operator*=(FunctionPatch<M>& f, const NumericType<M>& c);

template<class M> FunctionPatch<M> operator-(const FunctionPatch<M>& f);
template<class M> FunctionPatch<M> operator+(const FunctionPatch<M>& f1, const FunctionPatch<M>& f2);
template<class M> FunctionPatch<M> operator*(const FunctionPatch<M>& f1, const FunctionPatch<M>& f2);

template<class M, class OP> FunctionPatch<M> apply(OP op, FunctionPatch<M> const& f);

// Remove the error term
template<class M> FunctionPatch<M> midpoint(const FunctionPatch<M>& x);






/*! \ingroup FunctionModelSubModule
 *  \brief A Taylor function model with multivalued codomain built from the TaylorModel class.
 *
 *  See also TaylorModel, FunctionPatch<M>, VectorTaylorFunction.
 */
template<class M> class VectorFunctionPatch
    : public VectorFunctionModelMixin<VectorFunctionPatch<M>,typename M::Paradigm>
{
    friend class VectorFunctionPatchElementReference<M>;
    typedef typename M::Paradigm P;
  public:
    typedef ExactBox DomainType;
    typedef M ModelType;
    typedef Box<typename ModelType::CodomainType> CodomainType;
    typedef Box<typename ModelType::RangeType> RangeType;
    typedef typename ModelType::ExpansionType ExpansionType;
    typedef typename ModelType::CoefficientType CoefficientType;
    typedef typename ModelType::ErrorType ErrorType;
    typedef typename ModelType::NumericType NumericType;
    typedef P Paradigm;

    /*! \brief Default constructor constructs a Taylor model of order zero with no arguments and no result variables. */
    VectorFunctionPatch<M>();

    /*! \brief Construct the zero vector function over an unspecified domain. */
    explicit VectorFunctionPatch<M>(SizeType result_size);

    /*! \brief Construct from a result size and a domain. */
    VectorFunctionPatch<M>(SizeType result_size, const ExactBox& domain, Sweeper swp);

    /*! \brief Construct a vector function all of whose components are the same. */
    VectorFunctionPatch<M>(SizeType result_size, const FunctionPatch<M>& scalar_function);

    /*! \brief Construct from a domain and the expansion. */
    VectorFunctionPatch<M>(const ExactBox& domain,
                         const Vector<Expansion<ExactFloat64>>& expansion,
                         Sweeper swp);

    /*! \brief Construct from a domain, and expansion and errors. */
    VectorFunctionPatch<M>(const ExactBox& domain,
                         const Vector<Expansion<ExactFloat64>>& expansion,
                         const Vector<ErrorFloat64>& error,
                         Sweeper swp);

    /*! \brief Construct from a domain, and expansion and errors. */
    VectorFunctionPatch<M>(const ExactBox& domain,
                         const Vector<Expansion<RawFloat64>>& expansion,
                         const Vector<RawFloat64>& error,
                         Sweeper swp);

    /*! \brief Construct from a domain, and expansion and errors. */
    VectorFunctionPatch<M>(const ExactBox& domain,
                         const Vector<Expansion<RawFloat64>>& expansion,
                         Sweeper swp);

    /*! \brief Construct from a domain and the models. */
    explicit VectorFunctionPatch<M>(const ExactBox& domain, const Vector< ModelType >& variables);

    /*! \brief Construct from a domain, a function, and a sweeper determining the accuracy. */
    VectorFunctionPatch<M>(const ExactBox& domain,
                         const VectorFunctionType<M>& function,
                         const Sweeper& sweeper);

    /*! \brief Construct from a vector of scalar Taylor functions. */
    explicit VectorFunctionPatch<M>(const Vector<FunctionPatch<M>>& components);

    /*! \brief Construct from a list of scalar Taylor functions. */
    explicit VectorFunctionPatch<M>(const List<FunctionPatch<M>>& components);

    /*! \brief Construct from an initializer list of scalar Taylor functions. */
    VectorFunctionPatch<M>(InitializerList<FunctionPatch<M>> components);

    /*! \brief Construct from a vector expression. */
    template<class E> explicit VectorFunctionPatch<M>(const VectorExpression<E>& ve);

    explicit VectorFunctionPatch<M> (const VectorFunctionModel<ValidatedTag>& f);
    VectorFunctionPatch<M>& operator=(const VectorFunctionModel<ValidatedTag>& f);

    /*! \brief Equality operator. */
    Bool operator==(const VectorFunctionPatch<M>& p) const;
    /*! \brief Inequality operator. */
    Bool operator!=(const VectorFunctionPatch<M>& p) const;

    // Data access
    /*! \brief The sweeper used to control approximation of the Taylor function. */
    Sweeper sweeper() const;
    /*! \brief Set the sweeper used to control approximation of the Taylor function. */
    Void set_sweeper(Sweeper swp);
    /*! \brief The data used to define the domain of the Taylor model. */
    const ExactBox domain() const;
    /*! \brief A rough bound for the range of the function. */
    const ExactBox codomain() const;
    /*! \brief The centre of the Taylor model. */
    const Vector<CoefficientType> centre() const;
    /*! \brief The range of the Taylor model. */
    const UpperBox range() const;
    /*! \brief The data used to define the Taylor models. */
    const Vector<ModelType>& models() const;
    Vector<ModelType>& models();
    /*! \brief The data used to define the centre of the Taylor models. */
    const Vector<Expansion<ExactFloat64>> expansions() const;

    /*! \brief The \a i<sup>th</sup> Taylor model used to define the function. */
    const ModelType& model(SizeType i) const;
    /*! \brief The \a i<sup>th</sup> Taylor model used to define the function. */
    ModelType& model(SizeType i);

    /*! \brief The size of the argument. */
    SizeType argument_size() const;
    /*! \brief The size of the result. */
    SizeType result_size() const;

    // For compatibility wit Vector.
    SizeType size() const { return this->result_size(); }
    /*! \brief Get the \a ith Taylor variable */
    FunctionPatch<M> get(SizeType i) const;
    /*! \brief Set the \a ith Taylor variable */
    Void set(SizeType i, const FunctionPatch<M>& te);
    /*! \brief The \a ith Taylor variable */
    FunctionPatch<M> const operator[](SizeType i) const;
    /*! \brief The \a ith Taylor variable */
    VectorFunctionPatchElementReference<M> operator[](SizeType i);


    /*! \brief Evaluate the Taylor model at the point \a x. */
    Vector<ValidatedNumber> operator()(const Vector<ValidatedNumber>& x) const;
    Vector<ApproximateNumber> operator()(const Vector<ApproximateNumber>& x) const;
    Vector<ValidatedNumber> operator()(const Vector<ExactNumber>& x) const;
    /*! \brief Compute an approximation to Jacobian derivative of the Taylor model sat the point \a x. */
    Matrix<NumericType> jacobian(const Vector<NumericType>& x) const;

    //! \brief Remove all terms whose coefficient has magnitude
    //! lower than the cutoff threshold of the quantity.
    VectorFunctionPatch<M>& sweep();
    //! \brief Remove all terms as specified by \a sweeper.
    VectorFunctionPatch<M>& sweep(const SweeperInterface& sweeper);
    /*! \brief Set the error to zero. */
    Void clobber();

    /*! \brief The constant Taylor model with range \a r and argument domain \a d. */
    static VectorFunctionPatch<M> constant(const ExactBox& d, const Vector<NumericType>& r, Sweeper swp);
    /*! \brief The identity Taylor model on domain \a d. */
    static VectorFunctionPatch<M> identity(const ExactBox& d, Sweeper swp);
    //! \brief Return the vector of variables in the range with values \a x over domain \a d.
    static VectorFunctionPatch<M> projection(const ExactBox& d, SizeType imin, SizeType imax, Sweeper swp);

    /*! \brief Convert to an interval polynomial. */
    Vector<Polynomial<ValidatedFloat64>> polynomials() const;
    /*! \brief The vector of roundoff/truncation errors of each component. */
    Vector<ErrorType> const errors() const;
    /*! \brief The maximum roundoff/truncation error of the components. */
    ErrorType const error() const;
    //! \brief A multivalued function equal to the model on the domain.
    VectorFunctionType<M> function() const;

    /*! \brief Truncate terms higher than \a bd. */
    VectorFunctionPatch<M>& truncate(const MultiIndexBound& bd);
    /*! \brief Restrict to a subdomain. */
    Void restrict(const ExactBox& d);
    //! \brief Adjoin a scalar function.
    Void adjoin(const FunctionPatch<M>& sf);

    /*! \brief Write to an output stream. */
    OutputStream& write(OutputStream& os) const;

    /*! \brief Write a full representation to an output stream. */
    OutputStream& repr(OutputStream& os) const;

  private:
    Void _compute_jacobian() const;
    Void _set_argument_size(SizeType n);
    SizeType _compute_maximum_component_size() const;
    Void _resize(SizeType rs, SizeType as, ushort d, ushort s);
    virtual ScalarFunctionPatch<M>* _get(SizeType i) const { return new FunctionPatch<M>(this->_domain,this->_models[i]); }
    virtual VectorFunctionPatch<M>* _clone() const;
    virtual VectorFunctionPatch<M>* _create() const;
    virtual VectorFunctionPatch<M>* _create_identity() const;
    virtual ScalarFunctionPatch<M>* _create_zero() const;
  private:
    friend class VectorFunctionMixin<VectorFunctionPatch<M>,ValidatedTag>;
    friend class TaylorFunctionFactory;
  public:
    template<class X> Void _compute(Vector<X>& r, const Vector<X>& a) const;
  private:
    /* Domain of definition. */
    ExactBox _domain;
    Vector< ModelType > _models;
};


/*! \brief Inplace addition. */
template<class M> VectorFunctionPatch<M>& operator+=(VectorFunctionPatch<M>& f, const VectorFunctionPatch<M>& g);
/*! \brief Inplace subtraction. */
template<class M> VectorFunctionPatch<M>& operator-=(VectorFunctionPatch<M>& f, const VectorFunctionPatch<M>& g);
/*! \brief Inplace addition. */
template<class M> VectorFunctionPatch<M>& operator+=(VectorFunctionPatch<M>& f, const Vector<NumericType<M>>& c);
/*! \brief Inplace subtraction. */
template<class M> VectorFunctionPatch<M>& operator-=(VectorFunctionPatch<M>& f, const Vector<NumericType<M>>& c);
/*! \brief Inplace scalar multiplication. */
template<class M> VectorFunctionPatch<M>& operator*=(VectorFunctionPatch<M>& f, const NumericType<M>& c);
/*! \brief Inplace scalar division. */
template<class M> VectorFunctionPatch<M>& operator/=(VectorFunctionPatch<M>& f, const NumericType<M>& c);

/*! \brief Negation. */
template<class M> VectorFunctionPatch<M> operator-(const VectorFunctionPatch<M>& f);
/*! \brief Addition. */
template<class M> VectorFunctionPatch<M> operator+(const VectorFunctionPatch<M>& f1, const VectorFunctionPatch<M>& f2);
/*! \brief Subtraction. */
template<class M> VectorFunctionPatch<M> operator-(const VectorFunctionPatch<M>& f1, const VectorFunctionPatch<M>& f2);
/*! \brief Multiplication. */
template<class M> VectorFunctionPatch<M> operator*(const FunctionPatch<M>& f1, const VectorFunctionPatch<M>& f2);
/*! \brief Multiplication. */
template<class M> VectorFunctionPatch<M> operator*(const VectorFunctionPatch<M>& f1, const FunctionPatch<M>& f2);
/*! \brief Division. */
template<class M> VectorFunctionPatch<M> operator/(const VectorFunctionPatch<M>& f1, const FunctionPatch<M>& f2);

/*! \brief Addition of a constant. */
template<class M> VectorFunctionPatch<M> operator+(const VectorFunctionPatch<M>& f, const Vector<NumericType<M>>& c);
/*! \brief Subtraction of a constant. */
template<class M> VectorFunctionPatch<M> operator-(const VectorFunctionPatch<M>& f, const Vector<NumericType<M>>& c);
/*! \brief Multiplication by a scalar. */
template<class M> VectorFunctionPatch<M> operator*(const NumericType<M>& c, const VectorFunctionPatch<M>& f);
/*! \brief Multiplication by a scalar. */
template<class M> VectorFunctionPatch<M> operator*(const VectorFunctionPatch<M>& f, const NumericType<M>& c);
/*! \brief Division by a scalar. */
template<class M> VectorFunctionPatch<M> operator/(const VectorFunctionPatch<M>& f, const NumericType<M>& c);
/*! \brief Multiplication by a matrix. */
template<class M> VectorFunctionPatch<M> operator*(const Matrix<ExactNumber>& A, const VectorFunctionPatch<M>& f);
/*! \brief Multiplication by a matrix. */
template<class M> VectorFunctionPatch<M> operator*(const Matrix<NumericType<M>>& A, const VectorFunctionPatch<M>& f);

template<class M> VectorFunctionPatch<M> operator+(const VectorFunctionType<M>& f1, const VectorFunctionPatch<M>& tf2);
template<class M> VectorFunctionPatch<M> operator-(const VectorFunctionType<M>& f1, const VectorFunctionPatch<M>& tf2);
template<class M> VectorFunctionPatch<M> operator*(const ScalarFunctionType<M>& f1, const VectorFunctionPatch<M>& tf2);
template<class M> VectorFunctionPatch<M> operator*(const VectorFunctionType<M>& f1, const FunctionPatch<M>& tf2);
template<class M> VectorFunctionPatch<M> operator/(const VectorFunctionType<M>& f1, const FunctionPatch<M>& tf2);
template<class M> VectorFunctionPatch<M> operator+(const VectorFunctionPatch<M>& tf1, const VectorFunctionType<M>& f2);
template<class M> VectorFunctionPatch<M> operator-(const VectorFunctionPatch<M>& tf1, const VectorFunctionType<M>& f2);
template<class M> VectorFunctionPatch<M> operator*(const FunctionPatch<M>& tf1, const VectorFunctionType<M>& f2);
template<class M> VectorFunctionPatch<M> operator*(const VectorFunctionPatch<M>& tf1, const ScalarFunctionType<M>& f2);
template<class M> VectorFunctionPatch<M> operator/(const VectorFunctionPatch<M>& tf1, const ScalarFunctionType<M>& f2);

//! \brief Composition \f$f\circ g(x)=f(g(x))\f$.
template<class M> FunctionPatch<M> compose(const FunctionType<M>& f, const VectorFunctionPatch<M>& g);
//! \brief Composition \f$f\circ g(x)=f(g(x))\f$.
template<class M> VectorFunctionPatch<M> compose(const VectorFunctionType<M>& f, const VectorFunctionPatch<M>& g);
//! \brief Composition \f$f\circ g(x)=f(g(x))\f$.
template<class M> FunctionPatch<M> compose(const FunctionPatch<M>& f, const VectorFunctionPatch<M>& g);
//! \brief Composition \f$f\circ g(x)=f(g(x))\f$.
template<class M> VectorFunctionPatch<M> compose(const VectorFunctionPatch<M>& f, const VectorFunctionPatch<M>& g);

//! \brief Weak derivative of \a f with respect to variable \a k.
template<class M> VectorFunctionPatch<M> derivative(const VectorFunctionPatch<M>& f, SizeType k);
//! \brief Antiderivative of \a f with respect to variable \a k.
template<class M> VectorFunctionPatch<M> antiderivative(const VectorFunctionPatch<M>& f, SizeType k);
//! \brief Antiderivative of \a f with respect to variable \a k, taking value \c 0 when \a x[k]=c.
template<class M> VectorFunctionPatch<M> antiderivative(const VectorFunctionPatch<M>& f, SizeType k, ExactFloat64 c);

template<class M> NormType norm(const VectorFunctionPatch<M>& f);
template<class M> NormType distance(const VectorFunctionPatch<M>& f1, const VectorFunctionPatch<M>& f2);
template<class M> NormType distance(const VectorFunctionPatch<M>& f1, const VectorFunctionType<M>& f2);

//! \brief Restrict the function \a f to a subdomain \a d.
template<class M> VectorFunctionPatch<M> restriction(const VectorFunctionPatch<M>& f, const ExactBox& d);
//! \brief Restrict the function \a f to a larger domain \a d.
template<class M> VectorFunctionPatch<M> extension(const VectorFunctionPatch<M>& f, const ExactBox& d);

template<class M> VectorFunctionPatch<M> embed(const ExactBox& d1, const VectorFunctionPatch<M>& tv2,const ExactBox& d3);

//! \brief Tests if a function \a f refines another function \a g.
//! To be a refinement, the domain of \a f must contain the domain of \a g.
template<class M> Bool refines(const VectorFunctionPatch<M>& f, const VectorFunctionPatch<M>& g);
template<class M> Bool inconsistent(const VectorFunctionPatch<M>&, const VectorFunctionPatch<M>&);
template<class M> VectorFunctionPatch<M> refinement(const VectorFunctionPatch<M>&, const VectorFunctionPatch<M>&);

//! \brief Compute the function \f$(f \oplus g)(x)=(f(x),g(x))\f$.
template<class M> VectorFunctionPatch<M> join(const VectorFunctionPatch<M>& f, const VectorFunctionPatch<M>& g);
template<class M> VectorFunctionPatch<M> join(const VectorFunctionPatch<M>& f, const FunctionPatch<M>& g);
template<class M> VectorFunctionPatch<M> join(const FunctionPatch<M>& f, const FunctionPatch<M>& g);
template<class M> VectorFunctionPatch<M> join(const FunctionPatch<M>& f, const VectorFunctionPatch<M>& g);
//! \brief Compute the function \f$(f\otimes g)(x,y)=(f(x),g(y))\f$.
template<class M> VectorFunctionPatch<M> combine(const VectorFunctionPatch<M>& f, const VectorFunctionPatch<M>& g);
template<class M> VectorFunctionPatch<M> combine(const VectorFunctionPatch<M>& f, const FunctionPatch<M>& g);
template<class M> VectorFunctionPatch<M> combine(const FunctionPatch<M>& f, const VectorFunctionPatch<M>& g);
template<class M> VectorFunctionPatch<M> combine(const FunctionPatch<M>& f, const FunctionPatch<M>& g);

template<class M> VectorFunctionPatch<M> partial_evaluate(const VectorFunctionPatch<M>& f, SizeType k, const NumericType<M>& c);

template<class M> Vector<NumericType<M>> unchecked_evaluate(const VectorFunctionPatch<M>&, const Vector<NumericType<M>>&);
template<class M> FunctionPatch<M> unchecked_compose(const FunctionPatch<M>&, const VectorFunctionPatch<M>&);
template<class M> VectorFunctionPatch<M> unchecked_compose(const VectorFunctionPatch<M>&, const VectorFunctionPatch<M>&);

// Split the domain into halves along the \a j<sup>th</sup> coordinate.
template<class M> Pair<VectorFunctionPatch<M>,VectorFunctionPatch<M>> split(const VectorFunctionPatch<M>& x, SizeType j);

template<class M> OutputStream& operator<<(OutputStream&, const VectorFunctionPatch<M>&);

// Conversion operatations
template<class M> Polynomial<NumericType<M>> polynomial(const FunctionPatch<M>& tfn);
template<class M> Vector< Polynomial<NumericType<M>> > polynomials(const VectorFunctionPatch<M>& tfn);


// Sanitised output
template<class M> OutputStream& operator<<(OutputStream&, const Representation<FunctionPatch<M>>&);
template<class M> OutputStream& operator<<(OutputStream&, const Representation<VectorFunctionPatch<M>>&);
template<class F> struct ModelRepresentation { const F* pointer; double threshold; };
template<class F> ModelRepresentation<F> model_representation(const F& f, double swpt) {
    ModelRepresentation<F> r={&f,swpt}; return r; }
template<class M> OutputStream& operator<<(OutputStream&,const ModelRepresentation<FunctionPatch<M>>&);
template<class M> OutputStream& operator<<(OutputStream&,const ModelRepresentation<VectorFunctionPatch<M>>&);
template<class F> struct PolynomialRepresentation { const F* pointer; double threshold; List<String> names; };
template<class F> PolynomialRepresentation<F> polynomial_representation(const F& f, double swpt) {
    PolynomialRepresentation<F> r={&f,swpt}; return r; }
template<class F> PolynomialRepresentation<F> polynomial_representation(const F& f, double swpt, const List<String>& names) {
    PolynomialRepresentation<F> r={&f,swpt,names}; return r; }
template<class M> OutputStream& operator<<(OutputStream&,const PolynomialRepresentation<FunctionPatch<M>>&);
template<class M> OutputStream& operator<<(OutputStream&,const PolynomialRepresentation<VectorFunctionPatch<M>>&);


template<class M> template<class E> VectorFunctionPatch<M>::VectorFunctionPatch(const VectorExpression<E>& ve)
    : _domain(), _models(ve().size(),ve().zero_element().model())
{
    if(ve().size()!=0) { this->_domain=ve().zero_element().domain(); }
    for(SizeType i=0; i!=ve().size(); ++i) { this->set(i,ve()[i]); }
}

template<class M> class VectorFunctionPatchElementReference
{
    typedef M ModelType;
    typedef typename M::NumericType NumericType;
    typedef VectorFunctionPatchElementReference<M> SelfType;
    typedef FunctionPatch<M> FunctionType;
 public:
    VectorFunctionPatchElementReference(VectorFunctionPatch<M>& c, SizeType i) : _c(&c), _i(i) { }
    operator ScalarFunctionPatch<M> () const { return this->_c->get(this->_i); }
    Void operator=(const VectorFunctionPatchElementReference<M>& x) { this->_c->set(this->_i,x._c->get(x._i)); }
    Void operator=(const ScalarFunctionPatch<M>& x) { this->_c->set(this->_i,x); }
    ScalarFunctionPatch<M> element() const { return this->_c->get(this->_i); }
    ExactBox const& domain() const { return this->_c->domain(); }
    const ModelType& model() const { return this->_c->_models[this->_i]; }
    ErrorType error() const { return this->_c->_models[this->_i].error(); }
    Void set_error(const ErrorType& e) { this->_c->_models[this->_i].set_error(e); }
    Void sweep() { this->_c->_models[this->_i].sweep(); }
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





// Non-member template code
template<class M> FunctionPatch<M>::FunctionPatch(const ScalarFunctionModel<ValidatedTag>& f) {
     ARIADNE_ASSERT_MSG(dynamic_cast<const FunctionPatch<M>*>(f._ptr.operator->())," f="<<f);
     *this=dynamic_cast<const FunctionPatch<M>&>(*f._ptr);
}

// To scale from a model on [a,b] to a model on [c,d], use scale factor s=(d-c)/(b-a)
// and translation t=((c+d)-(a+b))/(b-a)
// Because we are scaling the model on [-1,+1], this is not the same as
// the mapping taking [a,b] to [c,d]
template<class M> FunctionPatch<M> partial_restriction(const FunctionPatch<M>& fp, SizeType k, const ExactInterval& ivl) {
    ExactBox dom=fp.domain(); dom[k]=ivl; return restriction(fp,dom);
}

template<class M> FunctionPatch<M> restriction(const FunctionPatch<M>& fp, const ExactBox& dom) {
    if(not(subset(dom,fp.domain()))) { ARIADNE_THROW(DomainException,"restiction(FunctionPatch<M>,ExactBox)","fp="<<fp<<", dom="<<dom); }
    return unchecked_compose(fp,FunctionPatch<M>::identity(dom,fp.sweeper()));
}

template<class M> FunctionPatch<M> extension(const FunctionPatch<M>& fp, const ExactBox& dom) {
    return unchecked_compose(fp,FunctionPatch<M>::identity(dom,fp.sweeper()));
}



template<class M, class OP> FunctionPatch<M> apply(OP op, FunctionPatch<M> const& f) {
    return FunctionPatch<M>(f.domain(),op(f.model()));
}

template<class M, class OP> FunctionPatch<M> apply(OP op, FunctionPatch<M> const& f1, FunctionPatch<M> const& f2) {
    if(f1.domain()==f2.domain()) {
        return FunctionPatch<M>(f1.domain(),op(f1.model(),f2.model()));
    } else {
        ExactBox domain=intersection(f1.domain(),f2.domain());
        return FunctionPatch<M>(domain,op(restriction(f1,domain).model(),restriction(f2,domain).model()));
    }
}


template<class M> FunctionPatch<M>& operator+=(FunctionPatch<M>& f, const NumericType<M>& c) {
    f.model()+=c; return f; }

template<class M> FunctionPatch<M>& operator*=(FunctionPatch<M>& f, const NumericType<M>& c) {
    f.model()*=c; return f;
}

template<class M> FunctionPatch<M> operator-(const FunctionPatch<M>& f) {
    return apply(Neg(),f); }

template<class M> FunctionPatch<M> operator+(const FunctionPatch<M>& f1, const FunctionPatch<M>& f2) {
    return apply(Add(),f1,f2); }

template<class M> FunctionPatch<M> operator-(const FunctionPatch<M>& f1, const FunctionPatch<M>& f2) {
    return apply(Sub(),f1,f2); }

template<class M> FunctionPatch<M> operator*(const FunctionPatch<M>& f1, const FunctionPatch<M>& f2) {
    return apply(Mul(),f1,f2); }

template<class M> FunctionPatch<M> max(const FunctionPatch<M>& f1, const FunctionPatch<M>& f2) {
    return apply(Max(),f1,f2); }

template<class M> FunctionPatch<M> min(const FunctionPatch<M>& f1, const FunctionPatch<M>& f2) {
    return apply(Min(),f1,f2); }

template<class M> FunctionPatch<M> abs(const FunctionPatch<M>& f) {
    return apply(Abs(),f); }


template<class M> FunctionPatch<M> antiderivative(const FunctionPatch<M>& x, SizeType k) {
    ValidatedNumber sf=rad_val(x.domain()[k]);
    return FunctionPatch<M>(x.domain(),antiderivative(x.model(),k)*sf); }

template<class M> FunctionPatch<M> derivative(const FunctionPatch<M>& x, SizeType k) {
    ValidatedNumber sf=rec(rad_val(x.domain()[k]));
    return FunctionPatch<M>(x.domain(),derivative(x.model(),k)*sf); }

template<class M> FunctionPatch<M> embed(const ExactBox& dom1, const FunctionPatch<M>& tv2,const ExactBox& dom3) {
    return FunctionPatch<M>(product(dom1,tv2.domain(),dom3),embed(dom1.size(),tv2.model(),dom3.size())); }

template<class M> NumericType<M> evaluate(const FunctionPatch<M>& f, const Vector<NumericType<M>>& x) {
    if(!definitely(contains(f.domain(),x))) {
        ARIADNE_THROW(DomainException,"evaluate(f,x) with f="<<f<<", x="<<x,"x is not an element of f.domain()="<<f.domain());
    }
    return unchecked_evaluate(f,x);
}

template<class M> NumericType<M> unchecked_evaluate(const FunctionPatch<M>& f, const Vector<NumericType<M>>& x) {
    return evaluate(f.model(),unscale(x,f.domain()));
}


template<class M> FunctionPatch<M> compose(const FunctionType<M>& g, const VectorFunctionPatch<M>& f) {
    return FunctionPatch<M>(f.domain(),g.evaluate(f.models()));
}

template<class M> FunctionPatch<M> compose(const FunctionPatch<M>& g, const VectorFunctionPatch<M>& f) {
    if(!subset(f.codomain(),g.domain())) {
        ARIADNE_THROW(DomainException,"compose(g,f) with g="<<g<<", f="<<f,"f.codomain()="<<f.codomain()<<" is not a subset of g.domain()="<<g.domain());
    }
    return unchecked_compose(g,f);
}

template<class M> FunctionPatch<M> unchecked_compose(const FunctionPatch<M>& g, const VectorFunctionPatch<M>& f) {
    return FunctionPatch<M>(f.domain(),compose(g.model(),unscale(f.models(),g.domain())));
}



template<class M> FunctionPatch<M> partial_evaluate(const FunctionPatch<M>& te, SizeType k, const NumericType<M>& c) {
    // Scale c to domain
    const SizeType as=te.argument_size();
    ARIADNE_ASSERT(k<as);
    const ExactBox& domain=te.domain();
    const ExactInterval& dk=domain[k];
    ValidatedNumber sc=(c-med_val(dk))/rad_val(dk);

    ExactBox new_domain(as-1);
    for(SizeType i=0; i!=k; ++i) { new_domain[i]=domain[i]; }
    for(SizeType i=k; i!=as-1; ++i) { new_domain[i]=domain[i+1]; }

    M new_model=partial_evaluate(te.model(),k,sc);

    return FunctionPatch<M>(new_domain,new_model);
}



template<class M> FunctionPatch<M> antiderivative(const FunctionPatch<M>& f, SizeType k, ExactFloat64 c) {
    ARIADNE_ASSERT(k<f.argument_size());
    ARIADNE_ASSERT(contains(f.domain()[k],c));

    FunctionPatch<M> g = antiderivative(f,k);
    VectorFunctionPatch<M> h = VectorFunctionPatch<M>::identity(g.domain(),g.sweeper());
    h[k] = FunctionPatch<M>::constant(g.domain(),c,g.sweeper());

    return g-compose(g,h);
}





template<class M> Pair<FunctionPatch<M>,FunctionPatch<M>> split(const FunctionPatch<M>& tv, SizeType j) {
    typedef M ModelType;
    Pair<ModelType,ModelType> models={split(tv.model(),j,SplitPart::LOWER),split(tv.model(),j,SplitPart::LOWER)};
    Pair<ExactBox,ExactBox> subdomains=split(tv.domain(),j);
    return make_pair(FunctionPatch<M>(subdomains.first,models.first),
                     FunctionPatch<M>(subdomains.second,models.second));

}

template<class M> Bool refines(const FunctionPatch<M>& tv1, const FunctionPatch<M>& tv2) {
    if(tv1.domain()==tv2.domain()) { return refines(tv1.model(),tv2.model()); }
    if(subset(tv2.domain(),tv1.domain())) { return refines(restriction(tv1,tv2.domain()).model(),tv2.model()); }
    else { return false; }
}

template<class M> Bool inconsistent(const FunctionPatch<M>& tv1, const FunctionPatch<M>& tv2) {
    if(tv1.domain()==tv2.domain()) {
        return inconsistent(tv1.model(),tv2.model());
    } else {
        ExactBox domain=intersection(tv1.domain(),tv2.domain());
        return inconsistent(restriction(tv1,domain).model(),restriction(tv2,domain).model());
    }
}

template<class M> FunctionPatch<M> refinement(const FunctionPatch<M>& tv1, const FunctionPatch<M>& tv2) {
    ARIADNE_ASSERT(tv1.domain()==tv2.domain());
    return FunctionPatch<M>(tv1.domain(),refinement(tv1.model(),tv2.model()));
}

template<class M> ErrorFloat64 norm(const FunctionPatch<M>& f) {
    return norm(f.model());
}

template<class M> ErrorFloat64 distance(const FunctionPatch<M>& f1, const FunctionPatch<M>& f2) {
    return norm(f1-f2);
}

template<class M> ErrorFloat64 distance(const FunctionPatch<M>& f1, const ValidatedScalarFunction& f2) {
    return distance(f1,FunctionPatch<M>(f1.domain(),f2,f1.sweeper()));
}


template<class M> FunctionPatch<M> midpoint(const FunctionPatch<M>& f) {
    M tm=f.model();
    tm.set_error(0u);
    return FunctionPatch<M>(f.domain(),tm);
}

template<class M> OutputStream& operator<<(OutputStream& os, const FunctionPatch<M>& tf) {
    return tf.write(os);
}

template<class M> OutputStream& operator<<(OutputStream& os, const Representation<FunctionPatch<M>>& tf) {
    return tf.pointer->repr(os);
}






#ifdef ARIADNE_UNDEF

template<class M> Bool check(const Vector<FunctionPatch<M>>& tv)
{
    for(SizeType i=0; i!=tv.size(); ++i) {
        if(tv.zero_element().domain()!=tv[i].domain()) { return false; }
    }
    return true;
}

template<class M> Vector<Expansion<ExactFloat64>> expansion(const Vector<FunctionPatch<M>>& x)
{
    Vector< Expansion<ExactFloat64> > r(x.size());
    for(SizeType i=0; i!=x.size(); ++i) {
        r[i]=x[i].expansion();
    }
    return r;
}

template<class M> Vector<ErrorFloat64> error(const Vector<FunctionPatch<M>>& x)
{
    Vector<ErrorFloat64> r(x.size());
    for(SizeType i=0; i!=x.size(); ++i) {
        r[i]=x[i].error();
    }
    return r;
}

template<class M> Vector<ExactFloat64> value(const Vector<FunctionPatch<M>>& x)
{
    Vector<ExactFloat64> r(x.size());
    for(SizeType i=0; i!=x.size(); ++i) {
        r[i]=x[i].value();
    }
    return r;
}

template<class M> Vector<UpperInterval> ranges(const Vector<FunctionPatch<M>>& x)
{
    Vector<UpperInterval> r(x.size());
    for(SizeType i=0; i!=x.size(); ++i) {
        r[i]=x[i].range();
    }
    return r;
}

#endif // ARIADNE_UNDEF





template<class M> VectorFunctionPatch<M> join(const VectorFunctionPatch<M>& f1, const FunctionPatch<M>& f2)
{
    ARIADNE_ASSERT_MSG(f1.domain()==f2.domain(),"f1="<<f1<<", f2="<<f2);
    return VectorFunctionPatch<M>(f1.domain(),join(f1.models(),f2.model()));
}

template<class M> VectorFunctionPatch<M> join(const VectorFunctionPatch<M>& f, const VectorFunctionPatch<M>& g)
{
    ARIADNE_ASSERT(f.domain()==g.domain());
    return VectorFunctionPatch<M>(f.domain(),join(f.models(),g.models()));
}

template<class M> VectorFunctionPatch<M> join(const FunctionPatch<M>& f1, const FunctionPatch<M>& f2)
{
    ARIADNE_ASSERT(f1.domain()==f2.domain());
    return VectorFunctionPatch<M>(f1.domain(),{f1.model(),f2.model()});
}

template<class M> VectorFunctionPatch<M> join(const FunctionPatch<M>& f1, const VectorFunctionPatch<M>& f2)
{
    ARIADNE_ASSERT(f1.domain()==f2.domain());
    return VectorFunctionPatch<M>(f1.domain(),join(f1.model(),f2.models()));
}

template<class M> VectorFunctionPatch<M> combine(const FunctionPatch<M>& f1, const FunctionPatch<M>& f2)
{
    return VectorFunctionPatch<M>(product(f1.domain(),f2.domain()),combine(Vector<M>(1u,f1.model()),{f2.model()}));
}

template<class M> VectorFunctionPatch<M> combine(const FunctionPatch<M>& f1, const VectorFunctionPatch<M>& f2)
{
    return VectorFunctionPatch<M>(product(f1.domain(),f2.domain()),combine({f1.model()},f2.models()));
}

template<class M> VectorFunctionPatch<M> combine(const VectorFunctionPatch<M>& f1, const FunctionPatch<M>& f2)
{
    return VectorFunctionPatch<M>(product(f1.domain(),f2.domain()),combine(f1.models(),{f2.model()}));
}

template<class M> VectorFunctionPatch<M> combine(const VectorFunctionPatch<M>& f1, const VectorFunctionPatch<M>& f2)
{
    return VectorFunctionPatch<M>(product(f1.domain(),f2.domain()),combine(f1.models(),f2.models()));
}


template<class M> VectorFunctionPatch<M> embed(const VectorFunctionPatch<M>& f, const ExactInterval& d)
{
    return embed(ExactBox(),f,ExactBox(1u,d));
}

template<class M> VectorFunctionPatch<M> embed(const VectorFunctionPatch<M>& f, const ExactBox& d)
{
    return embed(ExactBox(),f,d);
}

template<class M> VectorFunctionPatch<M> embed(const ExactBox& d, const VectorFunctionPatch<M>& f)
{
    return embed(d,f,ExactBox());
}

template<class M> VectorFunctionPatch<M> embed(const ExactBox& d1, const VectorFunctionPatch<M>& f, const ExactBox& d2)
{
    return VectorFunctionPatch<M>(product(d1,f.domain(),d2),embed(d1.size(),f.models(),d2.size()));
}

template<class M> VectorFunctionPatch<M> restriction(const VectorFunctionPatch<M>& f, const ExactBox& d)
{
    ARIADNE_ASSERT_MSG(subset(d,f.domain()),"Cannot restriction "<<f<<" to non-sub-domain "<<d);
    if(d==f.domain()) { return f; }
    VectorFunctionPatch<M> r(f.result_size(),d,f.sweeper());
    for(SizeType i=0; i!=r.result_size(); ++i) {
        r.set(i,restriction(f[i],d));
    }
    return r;
}

template<class M> Pair<VectorFunctionPatch<M>,VectorFunctionPatch<M>> split(const VectorFunctionPatch<M>& tf, SizeType j)
{
    typedef M ModelType;
    Pair< Vector<ModelType>,Vector<ModelType> > models=split(tf.models(),j);
    Pair<ExactBox,ExactBox> subdomains=split(tf.domain(),j);
    return make_pair(VectorFunctionPatch<M>(subdomains.first,models.first),
                     VectorFunctionPatch<M>(subdomains.second,models.second));

}

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

template<class M> VectorFunctionPatch<M>& operator+=(VectorFunctionPatch<M>& f, const VectorFunctionPatch<M>& g)
{
    ARIADNE_ASSERT(f.result_size()==g.result_size());
    ARIADNE_ASSERT(subset(f.domain(),g.domain()));
    ARIADNE_ASSERT(f.domain()==g.domain());
    f.models()+=g.models();
    return f;
}

template<class M> VectorFunctionPatch<M>& operator-=(VectorFunctionPatch<M>& f, const VectorFunctionPatch<M>& g)
{
    ARIADNE_ASSERT(f.result_size()==g.result_size());
    ARIADNE_ASSERT(subset(f.domain(),g.domain()));
    ARIADNE_ASSERT(f.domain()==g.domain());
    f.models()+=g.models();
    return f;
}

template<class M> VectorFunctionPatch<M>& operator+=(VectorFunctionPatch<M>& f, const Vector<NumericType<M>>& c)
{
    ARIADNE_ASSERT(f.result_size()==c.size());
    f.models()+=c;
    return f;
}

template<class M> VectorFunctionPatch<M>& operator-=(VectorFunctionPatch<M>& f, const Vector<NumericType<M>>& c)
{
    ARIADNE_ASSERT(f.result_size()==c.size());
    f.models()-=c;
    return f;
}

template<class M> VectorFunctionPatch<M>& operator*=(VectorFunctionPatch<M>& f, const NumericType<M>& c)
{
    f.models()*=c;
    return f;
}

template<class M> VectorFunctionPatch<M>& operator/=(VectorFunctionPatch<M>& f, const NumericType<M>& c)
{
    f.models()/=c;
    return f;
}


template<class M> VectorFunctionPatch<M> operator+(const VectorFunctionPatch<M>& f1, const VectorFunctionPatch<M>& f2)
{
    typedef M ModelType;
    ARIADNE_ASSERT_MSG(!empty(intersection(f1.domain(),f2.domain())),
                       "operator+(VectorFunctionPatch<M> f1, VectorFunctionPatch<M> f2) with f1="<<f1<<" f2="<<f2<<
                       ": domains are disjoint");
    if(f1.domain()==f2.domain()) {
        return VectorFunctionPatch<M>(f1.domain(),Vector<ModelType>(f1.models()+f2.models()));
    } else {
        ExactBox new_domain=intersection(f1.domain(),f2.domain());
        return operator+(restriction(f1,new_domain),restriction(f2,new_domain));
    }
}


template<class M> VectorFunctionPatch<M> operator-(const VectorFunctionPatch<M>& f1, const VectorFunctionPatch<M>& f2)
{
    typedef M ModelType;
    ARIADNE_ASSERT(!is_empty(intersection(f1.domain(),f2.domain())));
    if(f1.domain()==f2.domain()) {
        return VectorFunctionPatch<M>(f1.domain(),Vector<ModelType>(f1.models()-f2.models()));
    } else {
        ExactBox new_domain=intersection(f1.domain(),f2.domain());
        return operator-(restriction(f1,new_domain),restriction(f2,new_domain));
    }
}

template<class M> VectorFunctionPatch<M> operator*(const FunctionPatch<M>& f1, const VectorFunctionPatch<M>& f2)
{
    typedef M ModelType;
    ARIADNE_ASSERT(!is_empty(intersection(f1.domain(),f2.domain())));
    if(f1.domain()==f2.domain()) {
        return VectorFunctionPatch<M>(f1.domain(),Vector<ModelType>(f1.model()*f2.models()));
    } else {
        ExactBox new_domain=intersection(f1.domain(),f2.domain());
        return operator*(restriction(f1,new_domain),restriction(f2,new_domain));
    }
}

template<class M> VectorFunctionPatch<M> operator*(const VectorFunctionPatch<M>& f1, const FunctionPatch<M>& f2)
{
    typedef M ModelType;
    ARIADNE_ASSERT(!empty(intersection(f1.domain(),f2.domain())));
    if(f1.domain()==f2.domain()) {
        return VectorFunctionPatch<M>(f1.domain(),Vector<ModelType>(f1.models()*f2.model()));
    } else {
        ExactBox new_domain=intersection(f1.domain(),f2.domain());
        return operator*(restriction(f1,new_domain),restriction(f2,new_domain));
    }
}

template<class M> VectorFunctionPatch<M> operator/(const VectorFunctionPatch<M>& f1, const FunctionPatch<M>& f2)
{
    return f1 * rec(f2);
}



template<class M> VectorFunctionPatch<M> operator-(const VectorFunctionPatch<M>& f)
{
    return VectorFunctionPatch<M>(f.domain(),Vector<M>(-f.models()));
}

template<class M> VectorFunctionPatch<M> operator*(const NumericType<M>& c, const VectorFunctionPatch<M>& f)
{
    return VectorFunctionPatch<M>(f.domain(),Vector<M>(f.models()*c));
}

template<class M> VectorFunctionPatch<M> operator*(const VectorFunctionPatch<M>& f, const NumericType<M>& c)
{
    return VectorFunctionPatch<M>(f.domain(),Vector<M>(f.models()*c));
}

template<class M> VectorFunctionPatch<M> operator/(const VectorFunctionPatch<M>& f, const NumericType<M>& c)
{
    return VectorFunctionPatch<M>(f.domain(),Vector<M>(f.models()/c));
}

template<class M> VectorFunctionPatch<M> operator+(const VectorFunctionPatch<M>& f, const Vector<NumericType<M>>& c)
{
    return VectorFunctionPatch<M>(f.domain(),Vector<M>(f.models()+c));
}

template<class M> VectorFunctionPatch<M> operator-(const VectorFunctionPatch<M>& f, const Vector<NumericType<M>>& c)
{
    return VectorFunctionPatch<M>(f.domain(),Vector<M>(f.models()-c));
}

template<class M> VectorFunctionPatch<M> operator*(const Matrix<Float64>& A, const VectorFunctionPatch<M>& f)
{
    ARIADNE_PRECONDITION(A.column_size()==f.size());
    Vector<M> models(A.row_size(),M(f.argument_size(),f.sweeper()));
    for(SizeType i=0; i!=A.row_size(); ++i) {
        for(SizeType j=0; j!=A.column_size(); ++j) {
            models[i] += A.get(i,j) * f.model(j);
        }
    }
    return VectorFunctionPatch<M>(f.domain(),models);
}

template<class M> VectorFunctionPatch<M> operator*(const Matrix<NumericType<M>>& A, const VectorFunctionPatch<M>& f)
{
    ARIADNE_PRECONDITION(A.column_size()==f.size());
    Vector<M> models(A.row_size(),M(f.argument_size(),f.sweeper()));
    for(SizeType i=0; i!=A.row_size(); ++i) {
        for(SizeType j=0; j!=A.column_size(); ++j) {
            models[i] += A.get(i,j) * f.model(j);
        }
    }
    return VectorFunctionPatch<M>(f.domain(),models);
}

template<class M> VectorFunctionPatch<M> operator+(const ValidatedVectorFunction& f1, const VectorFunctionPatch<M>& tf2) {
    return VectorFunctionPatch<M>(tf2.domain(),f1,tf2.sweeper())+tf2; }
template<class M> VectorFunctionPatch<M> operator-(const ValidatedVectorFunction& f1, const VectorFunctionPatch<M>& tf2) {
    return VectorFunctionPatch<M>(tf2.domain(),f1,tf2.sweeper())-tf2; }
template<class M> VectorFunctionPatch<M> operator*(const ValidatedScalarFunction& f1, const VectorFunctionPatch<M>& tf2) {
    return FunctionPatch<M>(tf2.domain(),f1,tf2.sweeper())*tf2; }
template<class M> VectorFunctionPatch<M> operator*(const ValidatedVectorFunction& f1, const FunctionPatch<M>& tf2) {
    return VectorFunctionPatch<M>(tf2.domain(),f1,tf2.sweeper())*tf2; }
template<class M> VectorFunctionPatch<M> operator/(const ValidatedVectorFunction& f1, const FunctionPatch<M>& tf2) {
    return VectorFunctionPatch<M>(tf2.domain(),f1,tf2.sweeper())/tf2; }
template<class M> VectorFunctionPatch<M> operator+(const VectorFunctionPatch<M>& tf1, const ValidatedVectorFunction& f2) {
    return tf1+VectorFunctionPatch<M>(tf1.domain(),f2,tf1.sweeper()); }
template<class M> VectorFunctionPatch<M> operator-(const VectorFunctionPatch<M>& tf1, const ValidatedVectorFunction& f2) {
    return tf1-VectorFunctionPatch<M>(tf1.domain(),f2,tf1.sweeper()); }
template<class M> VectorFunctionPatch<M> operator*(const FunctionPatch<M>& tf1, const ValidatedVectorFunction& f2) {
    return tf1*VectorFunctionPatch<M>(tf1.domain(),f2,tf1.sweeper()); }
template<class M> VectorFunctionPatch<M> operator*(const VectorFunctionPatch<M>& tf1, const ValidatedScalarFunction& f2) {
    return tf1*FunctionPatch<M>(tf1.domain(),f2,tf1.sweeper()); }
template<class M> VectorFunctionPatch<M> operator/(const VectorFunctionPatch<M>& tf1, const ValidatedScalarFunction& f2) {
    return tf1/FunctionPatch<M>(tf1.domain(),f2,tf1.sweeper()); }






template<class M> VectorFunctionPatch<M> partial_evaluate(const VectorFunctionPatch<M>& tf, SizeType k, const NumericType<M>& c)
{
    // Scale c to domain
    const SizeType as=tf.argument_size();
    ARIADNE_ASSERT(k<as);
    const Vector<ExactInterval>& domain=tf.domain();
    const ExactInterval& dk=domain[k];
    NumericType<M> sc=(c-med_val(dk))/rad_val(dk);

    Vector<ExactInterval> new_domain(as-1);
    for(SizeType i=0; i!=k; ++i) { new_domain[i]=domain[i]; }
    for(SizeType i=k; i!=as-1; ++i) { new_domain[i]=domain[i+1]; }

    Vector<M> new_models=partial_evaluate(tf.models(),k,sc);

    return VectorFunctionPatch<M>(new_domain,new_models);
}


template<class M> VectorFunctionPatch<M> partial_restriction(const VectorFunctionPatch<M>& tf, SizeType k, const ExactInterval& d)
{
    VectorFunctionPatch<M> r(tf.result_size(),tf.domain(),tf.sweeper());
    for(SizeType i=0; i!=tf.result_size(); ++i) {
        r[i]=partial_restriction(tf[i],k,d);
    }
    return r;
}

template<class M> VectorFunctionPatch<M> restriction(const VectorFunctionPatch<M>& tf, SizeType k, const ExactInterval& d)
{
    return partial_restriction(tf,k,d);
}


template<class M> Vector<ValidatedNumber> evaluate(const VectorFunctionPatch<M>& f, const Vector<ValidatedNumber>& x) {
    if(!definitely(contains(f.domain(),x))) {
        ARIADNE_THROW(DomainException,"evaluate(f,x) with f="<<f<<", x="<<x,"x is not a subset of f.domain()="<<f.domain());
    }
    return unchecked_evaluate(f,x);
}

template<class M> Vector<NumericType<M>> unchecked_evaluate(const VectorFunctionPatch<M>& f, const Vector<NumericType<M>>& x) {
    return evaluate(f.models(),unscale(x,f.domain()));
}

template<class M> VectorFunctionPatch<M> compose(const VectorFunctionType<M>& g, const VectorFunctionPatch<M>& f) {
    return VectorFunctionPatch<M>(f.domain(),g.evaluate(f.models()));
}

template<class M> VectorFunctionPatch<M> compose(const VectorFunctionPatch<M>& g, const VectorFunctionPatch<M>& f)
{
    if(!subset(f.codomain(),g.domain())) {
        ARIADNE_THROW(DomainException,"compose(g,f) with g="<<g<<", f="<<f,"f.codomain()="<<f.codomain()<<" is not a subset of g.domain()="<<g.domain());
    }
    return unchecked_compose(g,f);
}


template<class M> VectorFunctionPatch<M> unchecked_compose(const VectorFunctionPatch<M>& g, const VectorFunctionPatch<M>& f)
{
    return VectorFunctionPatch<M>(f.domain(),compose(g.models(),unscale(f.models(),g.domain())));
}



template<class M> VectorFunctionPatch<M> derivative(const VectorFunctionPatch<M>& f, SizeType k)
{
    ARIADNE_ASSERT_MSG(k<f.argument_size(),"f="<<f<<"\n f.argument_size()="<<f.argument_size()<<" k="<<k);
    ValidatedNumber fdomkrad=rad_val(f.domain()[k]);
    VectorFunctionPatch<M> g=f;
    for(SizeType i=0; i!=g.size(); ++i) {
        g[i]=derivative(f[i],k);
    }
    return g;
}

template<class M> VectorFunctionPatch<M> antiderivative(const VectorFunctionPatch<M>& f, SizeType k)
{
    ARIADNE_ASSERT_MSG(k<f.argument_size(),"f="<<f<<"\n f.argument_size()="<<f.argument_size()<<" k="<<k);
    ValidatedNumber fdomkrad=rad_val(f.domain()[k]);
    VectorFunctionPatch<M> g=f;
    for(SizeType i=0; i!=g.size(); ++i) {
        g.models()[i].antidifferentiate(k);
        g.models()[i]*=fdomkrad;
    }
    return g;
}

template<class M> VectorFunctionPatch<M> antiderivative(const VectorFunctionPatch<M>& f, SizeType k, ExactNumber c)
{
    ARIADNE_ASSERT_MSG(k<f.argument_size(),"f="<<f<<"\n f.argument_size()="<<f.argument_size()<<" k="<<k);
    ValidatedNumber fdomkrad=rad_val(f.domain()[k]);
    VectorFunctionPatch<M> g=f;
    for(SizeType i=0; i!=g.size(); ++i) {
        g[i]=antiderivative(f[i],k,c);
    }
    return g;
}






template<class M> ErrorFloat64 norm(const VectorFunctionPatch<M>& f) {
    ErrorFloat64 res=0u;
    for(SizeType i=0; i!=f.result_size(); ++i) {
        res=max(res,norm(f[i]));
    }
    return res;
}

template<class M> ErrorFloat64 distance(const VectorFunctionPatch<M>& f1, const VectorFunctionPatch<M>& f2) {
    return norm(f1-f2);
}

template<class M> ErrorFloat64 distance(const VectorFunctionPatch<M>& f1, const ValidatedVectorFunction& f2) {
    return distance(f1,VectorFunctionPatch<M>(f1.domain(),f2,f1.sweeper()));
}


template<class M> OutputStream& operator<<(OutputStream& os, const Representation<VectorFunctionPatch<M>>& repr)
{
    return repr.pointer->repr(os);
}

template<class M> OutputStream& operator<<(OutputStream& os, const PolynomialRepresentation<VectorFunctionPatch<M>>& repr)
{
    const VectorFunctionPatch<M>& function = *repr.pointer;
    os << "[";
    for(SizeType i=0; i!=function.result_size(); ++i) {
        if(i!=0) { os << ","; }
        os << polynomial_representation(function[i],repr.threshold,repr.names);
    }
    return os << "]";
}

template<class M> OutputStream& operator<<(OutputStream& os, const PolynomialRepresentation< List<FunctionPatch<M>> >& repr)
{
    const List<FunctionPatch<M>>& functions = *repr.pointer;
    os << "[";
    for(SizeType i=0; i!=functions.size(); ++i) {
        if(i!=0) { os << ","; }
        os << polynomial_representation(functions[i],repr.threshold);
    }
    return os << "]";
}



template<class M> OutputStream& operator<<(OutputStream& os, const VectorFunctionPatch<M>& p)
{
    return p.write(os);
}

template<class M> Polynomial<NumericType<M>> polynomial(const FunctionPatch<M>& tfn) {
    return tfn.polynomial();
}

template<class M> Vector< Polynomial<NumericType<M>> > polynomials(const VectorFunctionPatch<M>& tfn) {
    return tfn.polynomials();
}


} // namespace Ariadne

#endif // ARIADNE_FUNCTION_PATCH_H
