/***************************************************************************
 *            differential.hpp
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

/*! \file differential.hpp
 *  \brief Differential algebra variables with a sparse representation.
 */

#ifndef ARIADNE_DIFFERENTIAL_HPP
#define ARIADNE_DIFFERENTIAL_HPP

#include <map>

#include "../utility/macros.hpp"
#include "../utility/array.hpp"
#include "../numeric/number.decl.hpp"
#include "../numeric/float.decl.hpp"
#include "../algebra/vector.hpp"
#include "../algebra/matrix.hpp"
#include "../algebra/multi_index.hpp"
#include "../algebra/series.hpp"
#include "../algebra/expansion.hpp"
#include "../algebra/operations.hpp"

#include "../algebra/univariate_differential.hpp"

namespace Ariadne {

template<class X> class Vector;
template<class X> class Matrix;
template<class X> class Series;

template<class I, class X> class Expansion;
template<class X> class Differential;
template<class X> class UnivariateDifferential;
template<class X> class NonAssignableDifferential;

typedef Differential<FloatDP> FloatDifferential;
typedef Differential<FloatDPApproximation> FloatApproximationDifferential;
typedef Differential<FloatDPBounds> FloatBoundsDifferential;
typedef Differential<UpperIntervalType> UpperIntervalDifferentialType;

template<class X> class DifferentialFactory {
    typedef PrecisionType<X> PR;
    PR _pr;
  public:
    DifferentialFactory<X>(PR const& pr) : _pr(pr) { }

    template<class Y, EnableIf<IsConstructible<X,Y,PR>> =dummy> X create(Y const& y) { return X(y,_pr); }
    template<class Y> UnivariateDifferential<X> create(UnivariateDifferential<Y> const&);
    template<class Y> Differential<X> create(Differential<Y> const&);
    Differential<X>const& create(NonAssignableDifferential<X> const& dx) { return dx; }
};
template<class X> inline DifferentialFactory<X> differential_factory(X const& x) { return DifferentialFactory<X>(x.precision()); }
template<class X> inline DifferentialFactory<X> factory(UnivariateDifferential<X> const& dx) {
    return DifferentialFactory<X>(dx.value().precision()); }
template<class X> inline DifferentialFactory<X> factory(Differential<X> const& dx) {
    return DifferentialFactory<X>(dx.value().precision()); }




//! \ingroup DifferentiationModule
//! \brief A class representing the partial derivatives of a scalar quantity
//! depending on multiple arguments.
//!
//! Based on a power series Expansion, centred on the point at which the partial derivatives are
//! being evaluated.
//!
//! \invariant The expansion is sorted using graded_sort(). In particular, the
//! total degree of the terms is increasing, and the linear terms appear in coordinate order.
template<class X>
class Differential
    : public DispatchTranscendentalAlgebraOperations<Differential<X>,X>
    , public ProvideConcreteGenericArithmeticOperators<Differential<X>>
{
    static_assert(!IsSame<X,FloatDPValue>::value,"");
    typedef Differential<X> SelfType;

    static const DegreeType MAX_DEGREE=65535;

    SortedExpansion<MultiIndex,X,GradedIndexLess> _expansion;
    DegreeType _degree;
  public:
    //! \brief The type used for a power series.
    typedef UnivariateDifferential<X> SeriesType;
    //! \brief The comparison used to sort the coefficients.
    typedef GradedLess IndexComparisonType;
    //! \brief The type of used to represent numbers in the coefficient.
    typedef SortedExpansion<MultiIndex,X,GradedIndexLess> ExpansionType;
    //! \brief The kind of information provided by the concrete values.
    typedef typename X::Paradigm Paradigm;
    //! \brief The type of used to represent numbers in the coefficient.
    typedef typename X::NumericType NumericType;
    //! \brief The type of used to represent the coefficients.
    typedef X ValueType;
    //! \brief The type of an Iterator through (index,coefficient) pairs..
    typedef typename ExpansionType::Iterator Iterator;
    //! \brief The type of a constant Iterator through (index,coefficient) pairs..
    typedef typename ExpansionType::ConstIterator ConstIterator;

    //! \brief Default constructor constructs a differential with degree zero in no variables.
    Differential() = delete;
    //! \brief Constructs a differential with degree \a deg in \a as variables, and zero element based on \a prs.
    template<class... PRS, EnableIf<IsConstructible<X,PRS...>> =dummy> explicit Differential(SizeType as, DegreeType deg, PRS... prs)
        : Differential(as,deg,X(prs...)) { }
    explicit Differential(SizeType as, DegreeType deg, X const& z);
    //! \brief Construct a differential from a mapping giving a coefficient for a finite number of multi-indices.
    explicit Differential(const Map<MultiIndex,X>& map, DegreeType deg);
    //! \brief Construct a differential of degree \a deg from the power-series expansion \a e.
    //! Terms in \a e of degree higher than \a deg are truncated
    explicit Differential(const Expansion<MultiIndex,X>& e, DegreeType deg);
    //! \brief Construct a differential of degree \a deg from an initializer list of (index,coefficient) pairs.
    explicit Differential(SizeType as, DegreeType deg, InitializerList< Pair<InitializerList<DegreeType>,X> > lst);
    //! \brief Construct a differential of degree \a deg from an initializer list of (index,coefficient) pairs.
    template<class... PRS, EnableIf<IsConstructible<X,Dbl,PRS...>> =dummy>
        explicit Differential(SizeType as, DegreeType deg, InitializerList< Pair<InitializerList<DegreeType>,Dbl> > lst, PRS... pr);
    //! \brief Construct a dense differential of degree \a deg from an initializer list of coefficient arranged in graded lexicographic order of the index.
    template<class... PRS, EnableIf<IsConstructible<X,Dbl,PRS...>> =dummy>
        explicit Differential(SizeType as, DegreeType deg, InitializerList<Dbl> lst, PRS... prs);

    //! \brief Conversion constructor from a different numerical type.
    template<class Y, class... PRS, EnableIf<IsConstructible<X,Y,PRS...>> =dummy> Differential(const Differential<Y>& dy, PRS... prs);

    //! \brief Set the differential equal to a constant, without changing the degree or number of arguments.
    Differential<X>& operator=(const X& c);
    template<class W, EnableIf<IsAssignable<X,W>> =dummy>
        Differential<X>& operator=(const W& c) { X xc=nul(this->value()); xc=c; return (*this)=xc; }

    //! \brief A constant differential of degree \a deg in \a as arguments with value \a c.
    static Differential<X> constant(SizeType as, DegreeType deg, const X& c);
    //! \brief A differential of degree \a deg in \a as arguments representing the quantity \f$v+x_j\f$.
    static Differential<X> variable(SizeType as, DegreeType deg, const X& v, SizeType j);
    //! \brief A differential of degree \a deg in \a as arguments representing a quantity with value \a v and gradient \a g.
    static Differential<X> affine(SizeType as, DegreeType d, const X& v, const Covector<X>& g);
    static Differential<X> affine(DegreeType d, const X& v, const Covector<X>& g);

    //! \brief A vector of \a rs constant differentials of degree \a deg in \a as arguments with values \a c[i].
    static Vector< Differential<X> > constants(SizeType rs, SizeType as, DegreeType deg, const Vector<X>& c);
    static Vector< Differential<X> > constants(SizeType as, DegreeType deg, const Vector<X>& c);
    //! \brief A vector of \a rs differentials of degree \a deg in \a as arguments with values \f$v_i+dx_i\f$.
    //! \pre \a rs == \a as == c.size().
    static Vector< Differential<X> > variables(SizeType rs, SizeType as, DegreeType deg, const Vector<X>& v);
    static Vector< Differential<X> > variables(DegreeType deg, const Vector<X>& v);
    //! \brief A vector of \a rs affine differentials of degree \a deg in \a as arguments with values \a v[i] and gradients \a G[i][j].
    static Vector<Differential<X>> affine(SizeType rs, SizeType as, DegreeType deg, const Vector<X>& v, const Matrix<X>& G);
    static Vector<Differential<X>> affine(DegreeType deg, const Vector<X>& v, const Matrix<X>& G);

//    static UnivariateDifferential<X> identity(DegreeType deg, const Vector<X>& v);
    static Differential<X> identity(DegreeType deg, const X& v); // FIXME: Should use UnivariateDifferentia; required for derivative of UnivariateFunction.
    static Vector<Differential<X>> identity(DegreeType deg, const Vector<X>& v);

    //! \brief Equality operator. Tests equality of representation, so comparing two differentials which are mathematically equal may return false if the structural zeros are different.
    EqualityType<X> operator==(const Differential<X>& other) const;
    //! \brief Inequality operator.
    InequalityType<X> operator!=(const Differential<X>& other) const;


    //! \brief The number of independent variables.
    SizeType argument_size() const;
    //! \brief The maximum degree of the stored terms.
    DegreeType degree() const;
    //! \brief A zero with the precision argument of the numbers.
    X zero_coefficient() const;
    //! \brief The internal representation of the polynomial expansion.
    const Expansion<MultiIndex,X>& expansion() const;
    //! \brief A reference to the internal representation of the polynomial expansion.
    Expansion<MultiIndex,X>& expansion();
    //! \brief The value of the differential i.e. the coefficient of \f$1\f$.
    const X& value() const;
    //! \brief The coefficient of \f$x_j\f$.
    const X& gradient(SizeType j) const;
    //! \brief The vector of coefficients of \f$x_j\f$.
    Covector<X> gradient() const;
    //! \brief The Hessian matrix.
    //! \note Note the the components of the Hessian matrix are \em half those of the values indexed by the differential.
    //! This is because the differential stores the coefficients of the Taylor expansion, rather than the derivatives themselves.
    Matrix<X> hessian() const;
    Matrix<X> half_hessian() const;

    //! \brief A reference to the coefficient of \f$x_j\f$.
    X& operator[](const SizeType& j);
    //! \brief The coefficient of \f$x_j\f$.
    const X& operator[](const SizeType& j) const;
    //! \brief The coefficient of \f$x^a\f$ i.e. \f$\prod x_j^{a_j}\f$.
    X& operator[](const MultiIndex& a);
    //! \brief A reference to the coefficient of \f$x^a\f$.
    const X& operator[](const MultiIndex& a) const;

    //! \brief Set the degree to be equal to \a d.
    Void set_degree(DegreeType d);
    //! \brief Set the coefficient of the constant term to \a c.
    Void set_value(const X& c);
    //! \brief Set the coefficient of the term in \f$x_j\f$ to \a d.
    Void set_gradient(SizeType j, const X& d);

    //! \brief An Iterator to the first structural nonzero.
    Iterator begin();
    //! \brief An Iterator to past-the-end structural nonzero.
    Iterator end();
    //! \brief A constant Iterator to the first structural nonzero.
    ConstIterator begin() const;
    //! \brief A constant Iterator to past-the-end structural nonzero.
    ConstIterator end() const;

    //! \brief The zero element of the differential algebra.
    Differential<X> create() const;
    Differential<X> create_zero() const;
    Differential<X> create_constant(const NumericType& c) const;
    Differential<X> create_constant(Int c) const;
    //! \brief Set all coefficients to zero.
    Void clear();
    //! \brief Remove all terms with coefficient \f$0\f$.
    Void cleanup();
    //! \brief Check that differential is sorted and all terms have degree less than maximum degree.
    Void check() const;

  public:
    //! \brief The derivatives of the composition of functions.
    template<class XX> friend Differential<XX> compose(UnivariateDifferential<XX> const& dx, Differential<XX> const& dy);
    //! \brief The derivatives of the composition of functions.
    template<class XX> friend Differential<XX> compose(Differential<XX> const& dx, Vector<Differential<XX>> const& dy);
    //! \brief The derivatives of the composition of functions.
    template<class XX> friend Vector<Differential<X>> compose(Vector<Differential<X>> const& dx, Vector<Differential<X>> const& dy);
    //! \brief Compute the differential of the derivative.
    template<class XX> friend Differential<XX> derivative(Differential<XX> const& dx, SizeType k);
    //! \brief Compute an antiderivative with respect to the variable \a i.
    template<class XX> friend Differential<XX> antiderivative(Differential<XX> const& dx, SizeType k);

    friend X value(const Differential<X>& x) { return x.value(); }
    friend Covector<X> gradient(const Differential<X>& x) { return x.gradient(); }
    friend Matrix<X> hessian(const Differential<X>& x) { return x.hessian(); }
    friend Matrix<X> half_hessian(const Differential<X>& x) { return x.half_hessian(); }

    friend OutputStream& operator<<(OutputStream& os, Differential<X> const& dx) { return dx._write(os); }
  public:
    friend Differential<X> min(const Differential<X>& x1, const Differential<X>& x2) {
        return AlgebraOperations<Differential<X>>::_min(x1,x2); }
    friend Differential<X> max(const Differential<X>& x1, const Differential<X>& x2) {
        return AlgebraOperations<Differential<X>>::_max(x1,x2); }
    friend Differential<X> abs(const Differential<X>& x) {
        return AlgebraOperations<Differential<X>>::_abs(x); }

    friend Differential<X> compose(UnivariateDifferential<X> const& x, Differential<X> const& y) {
        return Differential<X>::_compose(x,y); }
    friend Differential<X> compose(Differential<X> const& x, Vector<Differential<X>> const& y) {
        return Differential<X>::_compose(x,y); }
    friend Differential<X> derivative(Differential<X> const& x, SizeType k) { return
        Differential<X>::_derivative(x,k); }
    friend Differential<X> antiderivative(Differential<X> const& x, SizeType k) { return
        Differential<X>::_antiderivative(x,k); }

  public:
    static Differential<X> _derivative(Differential<X> const& dx, SizeType k);
    static Differential<X> _antiderivative(Differential<X> const& dx, SizeType k);
    static Differential<X> _compose(Series<X> const& f, Differential<X> const& dx) {
        return _compose(UnivariateDifferential<X>(dx.degree(),f),dx); }
    static Differential<X> _compose(UnivariateDifferential<X> const&, Differential<X> const&);
    static Differential<X> _compose(Differential<X> const&, Vector<Differential<X>> const&);
    static Vector<Differential<X>> _compose(Vector<Differential<X>> const&, Vector<Differential<X>> const&);

    static Matrix<X> _jacobian(Vector<Differential<X>> const&);
    OutputStream& _write(OutputStream& os) const;
};

template<class X, class Y, EnableIf<IsFloat<X>> =dummy, EnableIf<IsGenericNumericType<Y>> =dummy>
decltype(auto) operator+(Differential<X> const& x, Y const& y) { return x+factory(x).create(y); }

template<class X> struct AlgebraOperations<Differential<X>> : GradedAlgebraOperations<Differential<X>> {
    template<class OP> static Differential<X> apply(OP op, Differential<X> dx) {
        return compose(UnivariateDifferential<X>(op,dx.degree(),dx.value()),dx); }
    static Differential<X> apply(Pos op, Differential<X> dx);
    static Differential<X> apply(Neg op, Differential<X> dx);
    static Differential<X> apply(Add op, Differential<X> const& dx1, Differential<X> const& dx2);
    static Differential<X> apply(Sub op, Differential<X> const& dx1, Differential<X> const& dx2);
    static Differential<X> apply(Mul op, Differential<X> const& dx1, Differential<X> const& dx2);
    static Differential<X> apply(Div op, Differential<X> const& dx1, Differential<X> const& dx2);
    static Differential<X> apply(Pow op, Differential<X> const& dx, Int n);
    static Differential<X> apply(Add op, Differential<X> dx, X const& c);
    static Differential<X> apply(Mul op, Differential<X> dx, X const& c);
    static Differential<X> apply(Min op, Differential<X> const& dx1, Differential<X> const& dx2);
    static Differential<X> apply(Max op, Differential<X> const& dx1, Differential<X> const& dx2);
    static Differential<X> apply(Abs op, Differential<X> const& dx);
};


template<class X> template<class Y, class... PRS, EnableIf<IsConstructible<X,Y,PRS...>>>
Differential<X>::Differential(const Differential<Y>& x, PRS... prs)
    : _expansion(x.expansion(),prs...), _degree(x.degree())
{
}

template<class X> template<class... PRS, EnableIf<IsConstructible<X,Dbl,PRS...>>>
Differential<X>::Differential(SizeType as, DegreeType deg, InitializerList< Pair<InitializerList<DegreeType>,Dbl> > lst, PRS... prs)
    : Differential<X>(Expansion<MultiIndex,X>(lst,prs...),deg)
{ }

template<class X> template<class... PRS, EnableIf<IsConstructible<X,Dbl,PRS...>>>
Differential<X>::Differential(SizeType as, DegreeType deg, InitializerList<Dbl> lst, PRS... prs) : _expansion(as,X(prs...)), _degree(deg) {
    auto iter = lst.begin();
    for(MultiIndex j(as); j.degree()<=deg; ++j) {
        Dbl const& y=*iter; if(y!=0.0) { _expansion.append(j,X(y,prs...)); } ++iter; }
    this->cleanup();
}



template<class X> Vector<Differential<X>> compose(Vector<Differential<X>> const&, Vector<Differential<X>> const&);
template<class X> Vector<Differential<X>> derivative(Vector<Differential<X>> const&, SizeType k);
template<class X> Vector<Differential<X>> antiderivative(Vector<Differential<X>> const&, SizeType k);
template<class X> Vector<X> value(const Vector<Differential<X>>& x);
template<class X> Matrix<X> jacobian(Vector<Differential<X>> const&);

template<class X> Vector<Differential<X>> lie_derivative(const Vector<Differential<X> >& df, const Vector<Differential<X> >& dg);




template<class X>
class NonAssignableDifferential
    : public Differential<X>
{
  public:
    NonAssignableDifferential<X>& operator=(const Differential<X>& other) {
        ARIADNE_PRECONDITION_MSG( this->argument_size()==other.argument_size(),
                                  "\n  self="<<*this<<"\n  other="<<other<<"\n"<<std::flush  );
        ARIADNE_PRECONDITION( this->degree()==other.degree() );
        this->Differential<X>::operator=(other); return *this; }
    NonAssignableDifferential<X>& operator=(const X& c) {
        this->Differential<X>::operator=(c); return *this; }
    NonAssignableDifferential<X>& operator+=(const Differential<X>& x) {
        static_cast<Differential<X>&>(*this) += x; return *this; }
    // FIXME: Operators below should not be necessary
    friend Differential<X> operator+(const NonAssignableDifferential<X>& dx1, const NonAssignableDifferential<X>& dx2) {
        return static_cast<Differential<X>const&>(dx1) + static_cast<Differential<X>const&>(dx2); }
    friend Differential<X> operator*(const NonAssignableDifferential<X>& dx1, const NonAssignableDifferential<X>& dx2) {
        return static_cast<Differential<X>const&>(dx1) * static_cast<Differential<X>const&>(dx2); }
};

template<class X> class DifferentialCharacteristics {
    SizeType _argument_size; DegreeType _degree; X _zero_coefficient;
  public:
    DifferentialCharacteristics() : _argument_size(0u), _degree(0u), _zero_coefficient(0) { }
    DifferentialCharacteristics(SizeType as, DegreeType d) : _argument_size(as), _degree(d), _zero_coefficient(0) { }
    template<class PR> DifferentialCharacteristics(SizeType as, DegreeType d, PR pr) : _argument_size(as), _degree(d), _zero_coefficient(pr) { }
    DifferentialCharacteristics(const Differential<X>& d) : _argument_size(d.argument_size()), _degree(d.degree()), _zero_coefficient(d.zero_coefficient()) { }
    SizeType argument_size() const { return this->_argument_size; }
    DegreeType degree() const { return this->_degree; }
    X zero_coefficient() const { return this->_zero_coefficient; }
};

//! \brief A class representing the derivatives of a vector quantity depending on multiple arguments.
template<class X>
class Vector< Differential<X> >
    : public VectorExpression< Vector< Differential<X> > >
    , ProvideVectorOperations
{
    //BOOST_CONCEPT_ASSERT((DifferentialVectorConcept<DifferentialVector<X> >));
  public:
    DifferentialCharacteristics<X> _chars;
    Array< Differential<X> > _ary;
  public:
    // The type of the class
    typedef Vector< Differential<X> > SelfType;
    // The type used for accessing elements
    typedef SizeType IndexType;
    // The value stored in the vector.
    typedef Differential<X> ValueType;
    // The type used for scalars.
    typedef Differential<X> ScalarType;
    // The type used for scalars.
    typedef X NumericType;

    Vector() = delete;
    Vector(SizeType rs) = delete;
    template<class... PRS, EnableIf<IsConstructible<X,PRS...>> =dummy> Vector(SizeType rs, SizeType as, DegreeType d, PRS... prs)
        : Vector<Differential<X>>(rs,as,d,X(prs...)) { }
    Vector(SizeType rs, SizeType as, DegreeType d, X const& z);
    Vector(SizeType rs, const Differential<X>& sd);
    Vector(InitializerList<Differential<X>> const& lst);
    Vector(Array<Differential<X>> ary);
    template<class E> Vector(const VectorExpression<E>& ve);
    template<class E> Vector< Differential<X> >& operator=(const VectorExpression<E>& ve);
    template<class G, EnableIf<IsInvocableReturning<Differential<X>,G,SizeType>> =dummy> Vector(SizeType n, G const& g);
    template<class Y, class... PRS, EnableIf<IsConstructible<X,Y,PRS...>> =dummy> Vector(const Vector<Differential<Y>>& dv, PRS... prs);

    template<class... PRS, EnableIf<IsConstructible<X,Dbl,PRS...>> =dummy>
        Vector(SizeType rs, SizeType as, DegreeType d, InitializerList<InitializerList<Pair<InitializerList<DegreeType>,Dbl>>> lst, PRS... prs);
    template<class... PRS, EnableIf<IsConstructible<X,Dbl,PRS...>> =dummy>
        explicit Vector<Differential<X>>(SizeType rs, SizeType as, DegreeType deg, InitializerList<InitializerList<Dbl>> lst, PRS... prs);

    const Differential<X>& operator[](SizeType i) const { return this->_ary[i]; }
    NonAssignableDifferential<X>& operator[](SizeType i) { return static_cast<NonAssignableDifferential<X>&>(_ary[i]); }

    const Differential<X> zero_element() const { return Differential<X>::constant(this->argument_size(),this->degree(),this->zero_coefficient()); }

    NonAssignableDifferential<X>& at(SizeType i) { return static_cast<NonAssignableDifferential<X>&>(_ary[i]); }
    const Differential<X>& at(SizeType i) const { return this->_ary[i]; }
    const Differential<X>& get(SizeType i) const { return this->_ary[i]; }
    Void set(SizeType i, const Differential<X>& x) {
        ARIADNE_PRECONDITION(i<this->size()); ARIADNE_PRECONDITION(this->argument_size()==x.argument_size()); this->_ary[i]=x; }

    SizeType size() const { return this->_ary.size(); }
    SizeType result_size() const { return this->size(); }
    SizeType argument_size() const { return this->_chars.argument_size(); }
    DegreeType degree() const { return this->_chars.degree(); }
    X zero_coefficient() const { return this->_chars.zero_coefficient(); }
    const Array<Differential<X>>& array() const { return _ary; }

    Vector<X> value() const;
    Matrix<X> jacobian() const;

    Void set_value(const Vector<X>& c);

    //! \brief Set the degree to be equal to \a d.
    Void set_degree(DegreeType d);

    friend Vector<X> value(const Vector<Differential<X>>& x) { return x.value(); }
    friend Matrix<X> jacobian(const Vector<Differential<X>>& x) { return x.jacobian(); }
    //! \brief Compute the differential of the composiition \f$f\circ g\f$.
    friend Vector<Differential<X>> compose(const Vector<Differential<X>>& x, Vector<Differential<X>> const& y) { return _compose(x,y); }
    //! \brief Compute the differential of the composiition \f$f\circ g\f$.
    friend Vector<Differential<X>> derivative(const Vector<Differential<X>>& x, SizeType k) { return _derivative(x,k); }
    //! \brief Compute the differential of the antiderivative \f$ \int^{x_k} f(x_0,\ldots,x_{k-1},\xi,x_{k+1},\ldots) d\xi\f$.
    friend Vector<Differential<X>> antiderivative(const Vector<Differential<X>>& x, SizeType k) { return _antiderivative(x,k); }
    //! \brief Compute the differential of the Lie derivative \f$L_fg(x) = f(x) \cdot \nabla{g}(x)\f$.
    friend Vector<Differential<X>> lie_derivative(const Vector<Differential<X> >& df, const Vector<Differential<X> >& dg) { return _lie_derivative(df,dg); }

    //! \brief Compute the differential of the solutions \f$y=h(x)\f$ to \f$f(x,y)=0\f$. Requires \f$f(x0,y0)=0\f$.
    friend Vector<Differential<X>> solve(const Vector<Differential<X> >& df, const Vector<X>& y0) { return _solve(df,y0); }
    //! \brief Compute the differential of the flow \f$ \phi_{,t}(x0,t) = f(\phi(x0,t)) \f$ at \f$t=0\f$, with respect to the independent variables of
    friend Vector<Differential<X>> flow(const Vector<Differential<X> >& df, const Vector<X>& x0) { return _flow(df,x0); }

  public:
    static Vector<Differential<X>> _compose(Vector<Differential<X>> const& x, Vector<Differential<X>> const& y);
    static Vector<Differential<X>> _derivative(Vector<Differential<X>> const& x, SizeType k);
    static Vector<Differential<X>> _antiderivative(Vector<Differential<X>> const& x, SizeType k);
    static Vector<Differential<X>> _lie_derivative(const Vector<Differential<X> >& df, const Vector<Differential<X> >& dg);
    static Vector<Differential<X>> _solve(const Vector<Differential<X> >& df, const Vector<X>& y0);
    static Vector<Differential<X>> _flow(const Vector<Differential<X> >& df, const Vector<X>& x0);
};

template<class X> template<class Y, class... PRS, EnableIf<IsConstructible<X,Y,PRS...>>>
Vector<Differential<X>>::Vector(const Vector<Differential<Y>>& dv, PRS... prs)
    : _chars(dv.argument_size(),dv.degree(),X(prs...)), _ary(dv._ary,prs...) {
}

template<class X> template<class... PRS, EnableIf<IsConstructible<X,Dbl,PRS...>>>
Vector<Differential<X>>::Vector(SizeType rs, SizeType as, DegreeType d, InitializerList<InitializerList<Dbl>> lst, PRS... prs)
    : _chars(as,d), _ary(rs,Differential<X>(as,d)) {
    auto iter=lst.begin();
    for(SizeType i=0; i!=rs; ++i) { _ary[i]=Differential<X>(as,d,*iter,prs...); ++iter; }
}

template<class X> template<class... PRS, EnableIf<IsConstructible<X,Dbl,PRS...>>>
Vector<Differential<X>>::Vector(SizeType rs, SizeType as, DegreeType d, InitializerList<InitializerList<Pair<InitializerList<DegreeType>,Dbl>>> lst, PRS... prs)
    : _chars(as,d), _ary(rs,Differential<X>(as,d)) {
    auto iter=lst.begin();
    for(SizeType i=0; i!=rs; ++i) { _ary[i]=Differential<X>(as,d,*iter,prs...); ++iter; }
}

template<class X> template<class G, EnableIf<IsInvocableReturning<Differential<X>,G,SizeType>>>
Vector<Differential<X>>::Vector(SizeType n, G const& g)
    : _chars(g(0).argument_size(),g(0).degree()), _ary(n,g) { }

template<class X> template<class E>
Vector<Differential<X>>::Vector(const VectorExpression<E>& ve) : _ary(ve().size(),Uninitialised()) {
    for(SizeType i=0; i!=_ary.size(); ++i) { new (&_ary[i]) Differential<X>(ve()[i]); }
    ARIADNE_ASSERT(_ary.size()!=0); _chars=DifferentialCharacteristics<X>(_ary[0]);
}

template<class X> template<class E>
Vector<Differential<X>>& Vector<Differential<X>>::operator=(const VectorExpression<E>& ve) {
    _ary.resize(ve().size()); for(SizeType i=0; i!=_ary.size(); ++i) { static_cast<Differential<X>&>(_ary[i])=ve()[i]; }
    ARIADNE_ASSERT(_ary.size()!=0); _chars=DifferentialCharacteristics<X>(_ary[0]); return *this;
}





} //namespace Ariadne

#endif // ARIADNE_DIFFERENTIAL_HPP
