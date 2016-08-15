/***************************************************************************
 *            differential.h
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

/*! \file differential.h
 *  \brief Differential algebra variables with a sparse representation.
 */

#ifndef ARIADNE_DIFFERENTIAL_H
#define ARIADNE_DIFFERENTIAL_H

#include <map>

#include "utility/macros.h"
#include "utility/array.h"
#include "numeric/number.decl.h"
#include "numeric/float.decl.h"
#include "algebra/vector.h"
#include "algebra/matrix.h"
#include "algebra/multi_index.h"
#include "algebra/series.h"
#include "algebra/expansion.h"
#include "algebra/operations.h"
#include <boost/concept_check.hpp>

namespace Ariadne {

template<class X> class Vector;
template<class X> class Matrix;
template<class X> class Series;

template<class X> class Expansion;
template<class X> class Differential;

typedef Differential<Float64> FloatDifferential;
typedef Differential<Float64Approximation> FloatApproximationDifferential;
typedef Differential<Float64Bounds> FloatBoundsDifferential;
typedef Differential<UpperIntervalType> UpperIntervalDifferentialType;


template<class X> using EqualityType = decltype(declval<X>()==declval<X>());
template<class X> using InequalityType = decltype(declval<X>()!=declval<X>());




//! \ingroup DifferentiationModule
//! \brief Arbitrary-order derivatives with respect to a single argument.
template<class X> class UnivariateDifferential
    : public DispatchTranscendentalAlgebraOperations<UnivariateDifferential<X>,X>
{
    Array<X> _ary;
  public:
    typedef X NumericType;
    typedef typename X::Paradigm Paradigm;
    typedef UnivariateDifferential<X> SelfType;
    typedef typename Array<X>::Iterator Iterator;
    typedef typename Array<X>::ConstIterator ConstIterator;

    UnivariateDifferential();
    UnivariateDifferential(DegreeType d);
    UnivariateDifferential(DegreeType d, InitializerList<X> e);
    template<class XX> UnivariateDifferential(DegreeType d, XX const* p);
    UnivariateDifferential(DegreeType d, Series<X> const& s); // explicit

    static SelfType constant(DegreeType d, const NumericType& c);
    static SelfType variable(DegreeType d, const NumericType& c);

    SelfType create_zero() const;
    SelfType create_constant(const NumericType& c) const;
    SelfType create_variable(const NumericType& c) const;

    SizeType argument_size() const;
    DegreeType degree() const;
    const Array<X>& array() const;
    Expansion<X>& array();
    const X& operator[](SizeType k) const;
    X& operator[](SizeType k);

    SelfType& operator=(const NumericType& c);

    SelfType& operator+=(const SelfType& x);
    SelfType& operator-=(const SelfType& x);
    SelfType& operator*=(const SelfType& x);
    SelfType& operator+=(const NumericType& c);
    SelfType& operator*=(const NumericType& c);

    SelfType apply(Rec) const;

    X value() const;
    X gradient() const;
    X hessian() const;

    Void clear();
    OutputStream& write(OutputStream& os) const;

    friend OutputStream& operator<<(OutputStream& os, UnivariateDifferential<X> const& dx) {
        return dx.write(os); }
    static UnivariateDifferential<X> _compose(Series<X> const& f, Differential<X> const& dx);
};

template<class X> template<class XX> UnivariateDifferential<X>::UnivariateDifferential(DegreeType d, XX const* ptr)
    : _ary(d+1,X(0))
{
    std::copy(ptr,ptr+d+1,_ary.begin());
}



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
    , public DeclareMixedArithmeticOperators<Differential<X>,Int>
{
    typedef Differential<X> SelfType;

    static const DegreeType MAX_DEGREE=65535;
    static const X _zero;
    static const X _one;

    DegreeType _degree;
    SortedExpansion<X,GradedIndexLess> _expansion;
  public:
    //! \brief The type used for a power series.
    typedef UnivariateDifferential<X> SeriesType;
    //! \brief The comparison used to sort the coefficients.
    typedef GradedLess IndexComparisonType;
    //! \brief The type of used to represent numbers in the coefficient.
    typedef SortedExpansion<X,GradedIndexLess> ExpansionType;
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
    //! \brief Constructs a differential with degree \a deg in \a as variables.
    explicit Differential(SizeType as, DegreeType deg);
    //! \brief Construct a differential from a mapping giving a coefficient for a finite number of multi-indices.
    explicit Differential(const Map<MultiIndex,X>& map, DegreeType deg);
    //! \brief Construct a differential of degree \a deg from the power-series expansion \a e.
    //! Terms in \a e of degree higher than \a deg are truncated
    explicit Differential(const Expansion<X>& e, DegreeType deg);
    //! \brief Construct a differential of degree \a deg from an initializer list list of (index,coefficient) pairs.
    explicit Differential(SizeType as, DegreeType deg, InitializerList< Pair<InitializerList<DegreeType>,X> > lst);
    //! \brief Construct a differential of degree \a deg from an initializer list list of (index,coefficient) pairs.
    template<class PR, EnableIf<IsConstructible<X,Dbl,PR>> =dummy>
        explicit Differential(SizeType as, DegreeType deg, InitializerList< Pair<InitializerList<DegreeType>,Dbl> > lst, PR pr);

    //! \brief Construct a dense differential of degree \a deg in \a as variables from a list of coefficients beginning at \a ptr.
    template<class Y, class... PRS, EnableIf<IsConstructible<X,Y,PRS...>> =dummy> Differential(SizeType as, DegreeType deg, const Y* ptr, PRS... prs);
    //! \brief Conversion constructor from a different numerical type.
    template<class Y, class... PRS, EnableIf<IsConstructible<X,Y,PRS...>> =dummy> Differential(const Differential<Y>& dy, PRS... prs);


    //! \brief Set the differential equal to a constant, without changing the degree or number of arguments.
    Differential<X>& operator=(const X& c);

    //! \brief A constant differential of degree \a deg in \a as arguments with value \a c.
    static Differential<X> constant(SizeType as, DegreeType deg, const X& c);
    //! \brief A differential of degree \a deg in \a as arguments representing the quantity \f$v+x_j\f$.
    static Differential<X> variable(SizeType as, DegreeType deg, const X& v, SizeType j);
    //! \brief A vector of \a rs constant differentials of degree \a deg in \a as arguments with values \a c[j].
    static Vector< Differential<X> > constants(SizeType rs, SizeType as, DegreeType deg, const Vector<X>& c);
    //! \brief A vector of \a rs differentials of degree \a deg in \a as arguments with values \f$v_i+x_i\f$.
    //! \pre \a rs == \a as == c.size().
    static Vector< Differential<X> > variables(SizeType rs, SizeType as, DegreeType deg, const Vector<X>& x);

    //! \brief \brief A vector of differentials of degree \a deg in \a as arguments with values \f$c_i+x_i\f$.
    static Vector< Differential<X> > variables(DegreeType deg, const Vector<X>& x);

    //TODO: Complete UnivariateDifferential
    //static UnivariateDifferential<X> identity(DegreeType deg, const X& x);
    static Differential<X> identity(DegreeType deg, const X& x);
    static Vector< Differential<X> > identity(DegreeType deg, const Vector<X>& x);

    //! \brief Equality operator. Tests equality of representation, so comparing two differentials which are mathematically equal may return false if the structural zeros are different.
    EqualityType<X> operator==(const Differential<X>& other) const;
    //! \brief Inequality operator.
    InequalityType<X> operator!=(const Differential<X>& other) const;


    //! \brief The number of independent variables.
    SizeType argument_size() const;
    //! \brief The maximum degree of the stored terms.
    DegreeType degree() const;
    //! \brief The internal representation of the polynomial expansion.
    const Expansion<X>& expansion() const;
    //! \brief A reference to the internal representation of the polynomial expansion.
    Expansion<X>& expansion();
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

template<class X> struct AlgebraOperations<Differential<X>> : GradedAlgebraOperations<Differential<X>> {
    static Differential<X> _pos(Differential<X> dx);
    static Differential<X> _neg(Differential<X> dx);
    static Differential<X> _add(Differential<X> dx, X const& c);
    static Differential<X> _mul(Differential<X> dx, X const& c);
    static Differential<X> _add(Differential<X> const& dx1, Differential<X> const& dx2);
    static Differential<X> _sub(Differential<X> const& dx1, Differential<X> const& dx2);
    static Differential<X> _mul(Differential<X> const& dx1, Differential<X> const& dx2);
    static Differential<X> _div(Differential<X> const& dx1, Differential<X> const& dx2);
    static Differential<X> _rec(Differential<X> const& dx);
    static Differential<X> _min(Differential<X> const& dx1, Differential<X> const& dx2);
    static Differential<X> _max(Differential<X> const& dx1, Differential<X> const& dx2);
    static Differential<X> _abs(Differential<X> const& dx);
};

template<class X> template<class Y, class... PRS, EnableIf<IsConstructible<X,Y,PRS...>>>
Differential<X>::Differential(SizeType as, DegreeType deg, const Y* ptr, PRS... prs) : _expansion(as), _degree(deg) {
    for(MultiIndex j(as); j.degree()<=deg; ++j) {
        Y const& y=*ptr;
        if(!decide(y==0)) {
            _expansion.append(j,X(y,prs...));
        }
        ++ptr;
    }
    this->cleanup();
}

template<class X> template<class Y, class... PRS, EnableIf<IsConstructible<X,Y,PRS...>>>
Differential<X>::Differential(const Differential<Y>& x, PRS... prs)
    : _expansion(x.expansion(),prs...), _degree(x.degree())
{
}




template<class X> Vector<Differential<X>> compose(Vector<Differential<X>> const&, Vector<Differential<X>> const&);
template<class X> Vector<Differential<X>> derivative(Vector<Differential<X>> const&, SizeType k);
template<class X> Vector<Differential<X>> antiderivative(Vector<Differential<X>> const&, SizeType k);
template<class X> Vector<X> value(const Vector<Differential<X>>& x);
template<class X> Matrix<X> jacobian(Vector<Differential<X>> const&);

template<class X> Vector<Differential<X>> lie_derivative(const Vector<Differential<X> >& df, const Vector<Differential<X> >& dg);





template<class X>
struct NonAssignableDifferential
    : public Differential<X>
{
    NonAssignableDifferential<X>& operator=(const Differential<X>& other) {
        ARIADNE_PRECONDITION_MSG( this->argument_size()==other.argument_size(),
                                  "\n  self="<<*this<<"\n  other="<<other<<"\n"<<std::flush  );
        ARIADNE_PRECONDITION( this->degree()==other.degree() );
        this->Differential<X>::operator=(other); return *this; }
    NonAssignableDifferential<X>& operator=(const X& c) {
        this->Differential<X>::operator=(c); return *this; }
    NonAssignableDifferential<X>& operator+=(const Differential<X>& x) {
        static_cast<Differential<X>&>(*this) += x; return *this; }
};

class DifferentialCharacteristics {
    uint _argument_size; uint _degree;
  public:
    DifferentialCharacteristics() : _argument_size(0u), _degree(0u) { };
    DifferentialCharacteristics(uint as, uint d) : _argument_size(as), _degree(d) { };
    template<class X> DifferentialCharacteristics(const Differential<X>& d) : _argument_size(d.argument_size()), _degree(d.degree()) { };
    uint argument_size() const { return this->_argument_size; }
    uint degree() const { return this->_degree; }
};

//! \brief A class representing the derivatives of a vector quantity depending on multiple arguments.
template<class X>
class Vector< Differential<X> >
    : public VectorExpression< Vector< Differential<X> > >
    , ProvideVectorOperations
{
    //BOOST_CONCEPT_ASSERT((DifferentialVectorConcept<DifferentialVector<X> >));
  public:
    DifferentialCharacteristics _chars;
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
    Vector(SizeType rs, SizeType as, DegreeType d);
    Vector(SizeType rs, const Differential<X>& sd);
    Vector(SizeType rs, const Differential<X>* p);
    template<class Y, class... PRS, EnableIf<IsConstructible<X,Y,PRS...>> =dummy> Vector(const Vector<Differential<Y>>& dv, PRS... prs);
    template<class Y, class... PRS, EnableIf<IsConstructible<X,Y,PRS...>> =dummy> Vector(SizeType rs, SizeType as, DegreeType d, const Y* ptr, PRS... prs);
    Vector(SizeType rs, SizeType as, DegreeType d,const Vector<X>& v, const Matrix<X>& A);
    template<class E> Vector(const VectorExpression<E>& ve);
    template<class E> Vector< Differential<X> >& operator=(const VectorExpression<E>& ve);


    const Differential<X>& operator[](SizeType i) const { return this->_ary[i]; }
    NonAssignableDifferential<X>& operator[](SizeType i) { return static_cast<NonAssignableDifferential<X>&>(_ary[i]); }

    const Differential<X> zero_element() const { return Differential<X>(this->argument_size(),this->degree()); }

    NonAssignableDifferential<X>& at(SizeType i) { return static_cast<NonAssignableDifferential<X>&>(_ary[i]); }
    const Differential<X>& get(SizeType i) const { return this->_ary[i]; }
    Void set(SizeType i, const Differential<X>& x) {
        ARIADNE_PRECONDITION(i<this->size()); ARIADNE_PRECONDITION(this->argument_size()==x.argument_size()); this->_ary[i]=x; }

    SizeType size() const { return this->_ary.size(); }
    SizeType result_size() const { return this->size(); }
    SizeType argument_size() const { return this->_chars.argument_size(); }
    DegreeType degree() const { return this->_chars.degree(); }

    Vector<X> value() const;
    Matrix<X> jacobian() const;

    Void set_value(const Vector<X>& c);

    static Vector<Differential<X>> constant(SizeType rs, SizeType as, DegreeType d, const Vector<X>& c);
    static Vector<Differential<X>> variable(SizeType rs, SizeType as, DegreeType d, const Vector<X>& x);
    static Vector<Differential<X>> affine(SizeType rs, SizeType as, DegreeType d, const Vector<X>& b, const Matrix<X>& A);

    friend Vector<X> value(const Vector<Differential<X>>& x) { return x.value(); }
    friend Matrix<X> jacobian(const Vector<Differential<X>>& x) { return x.jacobian(); }
    friend Vector<Differential<X>> compose(const Vector<Differential<X>>& x, Vector<Differential<X>> const& y) { return _compose(x,y); }
    friend Vector<Differential<X>> derivative(const Vector<Differential<X>>& x, SizeType k) { return _derivative(x,k); }
    friend Vector<Differential<X>> antiderivative(const Vector<Differential<X>>& x, SizeType k) { return _antiderivative(x,k); }
    friend Vector<Differential<X>> lie_derivative(const Vector<Differential<X> >& df, const Vector<Differential<X> >& dg) { return _lie_derivative(df,dg); }

  public:
    static Vector<Differential<X>> _compose(Vector<Differential<X>> const& x, Vector<Differential<X>> const& y);
    static Vector<Differential<X>> _derivative(Vector<Differential<X>> const& x, SizeType k);
    static Vector<Differential<X>> _antiderivative(Vector<Differential<X>> const& x, SizeType k);
    static Vector<Differential<X>> _lie_derivative(const Vector<Differential<X> >& df, const Vector<Differential<X> >& dg);
};

template<class X> template<class Y, class... PRS, EnableIf<IsConstructible<X,Y,PRS...>>>
Vector<Differential<X>>::Vector(const Vector<Differential<Y>>& dv, PRS... prs)
    : _chars(dv._chars), _ary(dv._ary,prs...) {
}

template<class X> template<class Y, class... PRS, EnableIf<IsConstructible<X,Y,PRS...>>>
Vector<Differential<X>>::Vector(SizeType rs, SizeType as, DegreeType d, const Y* ptr, PRS... prs) {
    for(SizeType i=0; i!=rs; ++i) { for(MultiIndex j(as); j.degree()<=d; ++j) { if(*ptr!=0) { (*this)[i][j]=X(*ptr,prs...); } ++ptr; } }
}

template<class X> template<class E>
Vector<Differential<X>>::Vector(const VectorExpression<E>& ve) : _ary(ve().size(),Uninitialised()) {
    for(SizeType i=0; i!=_ary.size(); ++i) { new (&_ary[i]) Differential<X>(ve()[i]); }
    ARIADNE_ASSERT(_ary.size()!=0); _chars=DifferentialCharacteristics(_ary[0]);
}

template<class X> template<class E>
Vector<Differential<X>>& Vector<Differential<X>>::operator=(const VectorExpression<E>& ve) {
    _ary.resize(ve().size()); for(SizeType i=0; i!=_ary.size(); ++i) { static_cast<Differential<X>&>(_ary[i])=ve()[i]; }
    ARIADNE_ASSERT(_ary.size()!=0); _chars=DifferentialCharacteristics(_ary[0]); return *this;
}

template<class X> template<class PR, EnableIf<IsConstructible<X,Dbl,PR>>>
Differential<X>::Differential(SizeType as, DegreeType deg, InitializerList< Pair<InitializerList<DegreeType>,Dbl> > lst, PR pr)
    : Differential<X>(Expansion<X>(lst,pr),deg)
{ }




} //namespace Ariadne

#endif // ARIADNE_DIFFERENTIAL_H
