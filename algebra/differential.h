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
#include "algebra/algebra_operations.h"
#include <boost/concept_check.hpp>

namespace Ariadne {

template<class X> class Vector;
template<class X> class Matrix;
template<class X> class Series;

template<class X> class Expansion;
template<class X> class Differential;

typedef Differential<Float> FloatDifferential;
typedef Differential<ApproximateFloat> ApproximateFloatDifferential;
typedef Differential<ValidatedFloat> ValidatedFloatDifferential;
typedef Differential<UpperInterval> UpperIntervalDifferential;


template<class X> using EqualityType = decltype(declval<X>()==declval<X>());
template<class X> using InequalityType = decltype(declval<X>()!=declval<X>());




//! \ingroup DifferentiationModule
//! \brief Arbitrary-order derivatives with respect to a single argument.
template<class X> class UnivariateDifferential
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
{
    typedef Differential<X> SelfType;

    static const DegreeType MAX_DEGREE=65535;
    static const X _zero;
    static const X _one;

    DegreeType _degree;
    SortedExpansion<X,GradedKeyLess> _expansion;
  public:
    //! \brief The type used for a power series.
    typedef UnivariateDifferential<X> SeriesType;
    //! \brief The comparison used to sort the coefficients.
    typedef GradedKeyLess ComparisonType;
    //! \brief The type of used to represent numbers in the coefficient.
    typedef SortedExpansion<X,GradedKeyLess> ExpansionType;
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
    explicit Differential();
    //! \brief Constructs a differential with degree \a deg in \a as variables.
    explicit Differential(SizeType as, DegreeType deg);
    //! \brief Construct a differential from a mapping giving a coefficient for a finite number of multi-indices.
    explicit Differential(const Map<MultiIndex,X>& map, DegreeType deg);
    //! \brief Construct a differential of degree \a deg from the power-series expansion \a e.
    //! Terms in \a e of degree higher than \a deg are truncated
    explicit Differential(const Expansion<X>& e, DegreeType deg);
    //! \brief Construct a differential of degree \a deg from an initializer list list of (index,coefficient) pairs.
    explicit Differential(SizeType as, DegreeType deg, InitializerList< PairType<InitializerList<Int>,X> > lst);

    //! \brief Construct a dense differential of degree \a deg in \a as variables from a list of coefficients beginning at \a ptr.
    template<class XX> Differential(SizeType as, DegreeType deg, const XX* ptr);
    //! \brief Conversion constructor from a different numerical type.
    template<class XX> Differential(const Differential<XX>& x);


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
    //! \brief Inplace addition of a constant.
    friend Differential<X>& operator+=(Differential<X>& x, const NumericType& c) { return Differential<X>::_iadd(x,c); }
    //! \brief Inplace multiplication by a constant.
    friend Differential<X>& operator*=(Differential<X>& x, const NumericType& c) { return Differential<X>::_imul(x,c); }

    //! \brief Unary plus.
    friend Differential<X> operator+(Differential<X> x) { return Differential<X>::_pos(std::move(x)); }
    //! \brief Negation.
    friend Differential<X> operator-(Differential<X> x) { return Differential<X>::_neg(std::move(x)); }

    //! \brief Addition.
    friend Differential<X> operator+(const Differential<X>& x1, const Differential<X>& x2) {
        return Differential<X>::_add(x1,x2); }
    //! \brief Subtraction.
    friend Differential<X> operator-(const Differential<X>& x1, const Differential<X>& x2) {
        return Differential<X>::_sub(x1,x2); }
    //! \brief Multiplication.
    friend Differential<X> operator*(const Differential<X>& x1, const Differential<X>& x2) {
        return Differential<X>::_mul(x1,x2); }
    //! \brief Division.
    friend Differential<X> operator/(const Differential<X>& x1, const Differential<X>& x2) {
        return Differential<X>::_div(x1,x2); }

    friend Differential<X> min(const Differential<X>& x1, const Differential<X>& x2) {
        return Differential<X>::_min(x1,x2); }
    friend Differential<X> max(const Differential<X>& x1, const Differential<X>& x2) {
        return Differential<X>::_max(x1,x2); }
    friend Differential<X> abs(const Differential<X>& x) {
        return Differential<X>::_abs(x); }

  private:
    static Differential<X> _pos(Differential<X> dx);
    static Differential<X> _neg(Differential<X> dx);
    static Differential<X>& _iadd(Differential<X>& dx, X const& c);
    static Differential<X>& _imul(Differential<X>& dx, X const& c);
    static Differential<X> _add(Differential<X> const& dx1, Differential<X> const& dx2);
    static Differential<X> _sub(Differential<X> const& dx1, Differential<X> const& dx2);
    static Differential<X> _mul(Differential<X> const& dx1, Differential<X> const& dx2);
    static Differential<X> _div(Differential<X> const& dx1, Differential<X> const& dx2);
    static Differential<X> _min(Differential<X> const& dx1, Differential<X> const& dx2);
    static Differential<X> _max(Differential<X> const& dx1, Differential<X> const& dx2);
    static Differential<X> _abs(Differential<X> const& dx);
    static Differential<X> _derivative(Differential<X> const& dx, SizeType k);
    static Differential<X> _antiderivative(Differential<X> const& dx, SizeType k);
    static Differential<X> _compose(UnivariateDifferential<X> const&, Differential<X> const&);
    static Differential<X> _compose(Differential<X> const&, Vector<Differential<X>> const&);
    static Vector<Differential<X>> _compose(Vector<Differential<X>> const&, Vector<Differential<X>> const&);

    static Matrix<X> _jacobian(Vector<Differential<X>> const&);
    OutputStream& _write(OutputStream& os) const;
};

template<class X> struct IsAlgebra<Differential<X>> : True { };
template<class X> struct IsGradedAlgebra<Differential<X>> : True { };

template<class X> template<class XX> Differential<X>::Differential(SizeType as, DegreeType deg, const XX* ptr) : _expansion(as), _degree(deg) {
    for(MultiIndex j(as); j.degree()<=deg; ++j) {
        XX const& x=*ptr;
        if(!decide(x==0)) {
            _expansion.append(j,*ptr);
        }
        ++ptr;
    }
    this->cleanup();
}

template<class X> template<class XX> Differential<X>::Differential(const Differential<XX>& x)
    : _expansion(x.expansion()), _degree(x.degree())
{
}

template<class X> Differential<X> compose(UnivariateDifferential<X> const& x, Differential<X> const& y) { return Differential<X>::_compose(x,y); }
template<class X> Differential<X> compose(Differential<X> const& x, Vector<Differential<X>> const& y) { return Differential<X>::_compose(x,y); }
template<class X> Differential<X> derivative(Differential<X> const& x, SizeType k) { return Differential<X>::_derivative(x,k); }
template<class X> Differential<X> antiderivative(Differential<X> const& x, SizeType k) { return Differential<X>::_antiderivative(x,k); }

template<class X> inline Differential<X> compose(Series<X> const& x, Differential<X> const& y) {
    return compose(UnivariateDifferential<X>(y.degree(),x),y); }




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
        ARIADNE_PRECONDITION( this->argument_size()==other.argument_size()  );
        ARIADNE_PRECONDITION( this->degree()==other.degree() );
        this->Differential<X>::operator=(other); return *this; }
    NonAssignableDifferential<X>& operator=(const X& c) {
        this->Differential<X>::operator=(c); return *this; }
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

    Vector() : _chars(), _ary(0) { }
    Vector(SizeType rs) : _chars(), _ary(rs) { }
    Vector(SizeType rs, SizeType as, DegreeType d) : _chars(as,d), _ary(rs) {
        for(SizeType i=0; i!=rs; ++i) { this->_ary[i]=Differential<X>(as,d); } }
    Vector(SizeType rs, const Differential<X>& sd) : _chars(sd), _ary(rs) {
        for(SizeType i=0; i!=rs; ++i) { this->_ary[i]=sd; } }
    Vector(SizeType rs, const Differential<X>* p) : _chars(), _ary(rs) {
        ARIADNE_ASSERT(rs>0); _chars=DifferentialCharacteristics(p[0]); for(SizeType i=0; i!=rs; ++i) { this->_ary[i]=p[i]; } }
    template<class XX> Vector(const Vector< Differential<XX> > dv) : _chars(dv._chars), _ary(dv._ary) { }
    template<class XX> Vector(SizeType rs, SizeType as, DegreeType d, const XX* ptr);
    Vector(SizeType rs, SizeType as, DegreeType d,const Vector<X>& v, const Matrix<X>& A);
    template<class E> Vector(const VectorExpression<E>& ve) : _ary(ve().size()) {
        for(SizeType i=0; i!=_ary.size(); ++i) { static_cast<Differential<X>&>(_ary[i])=ve()[i]; }
        ARIADNE_ASSERT(_ary.size()!=0); _chars=DifferentialCharacteristics(_ary[0]); }
    template<class E> Vector< Differential<X> >& operator=(const VectorExpression<E>& ve) {
        _ary.resize(ve().size()); for(SizeType i=0; i!=_ary.size(); ++i) { static_cast<Differential<X>&>(_ary[i])=ve()[i]; }
        ARIADNE_ASSERT(_ary.size()!=0); _chars=DifferentialCharacteristics(_ary[0]); return *this; }


    const Differential<X>& operator[](SizeType i) const { return this->_ary[i]; }
    NonAssignableDifferential<X>& operator[](SizeType i) { return static_cast<NonAssignableDifferential<X>&>(_ary[i]); }

    const Differential<X> zero_element() const { return Differential<X>(this->argument_size(),this->degree()); }

    NonAssignableDifferential<X>& at(SizeType i) { return static_cast<NonAssignableDifferential<X>&>(_ary[i]); }
    const Differential<X>& get(SizeType i) const { return this->_ary[i]; }
    Void set(SizeType i, const Differential<X>& x) {
        ARIADNE_PRECONDITION(i<this->size());
        ARIADNE_PRECONDITION(this->argument_size()==x.argument_size());
        this->_ary[i]=x;
    }

    SizeType size() const { return this->_ary.size(); }
    SizeType result_size() const { return this->size(); }
    SizeType argument_size() const { return this->_chars.argument_size(); }
    DegreeType degree() const { return this->_chars.degree(); }

    Vector<X> value() const {
        Vector<X> r(this->result_size());
        for(SizeType i=0; i!=r.size(); ++i) { r[i]=(*this)[i].value(); } return r; }
    Matrix<X> jacobian() const { Matrix<X> r(this->result_size(),this->argument_size());
        for(SizeType i=0; i!=r.row_size(); ++i) { for(SizeType j=0; j!=r.column_size(); ++j) { r[i][j]=(*this)[i].gradient(j); } } return r; }

    Void set_value(const Vector<X>& c) {
        ARIADNE_ASSERT(this->result_size()==c.size());
        for(SizeType i=0; i!=c.size(); ++i) { (*this)[i].set_value(c[i]); } }

    Vector<Differential<X>> constant(SizeType rs, SizeType as, DegreeType d, const Vector<X>& c) {
        ARIADNE_ASSERT(c.size()==rs);
        Vector< Differential<X> > result(rs,as,d);
        for(SizeType i=0; i!=rs; ++i) { result[i]=c[i]; }
        return result;
    }

    Vector<Differential<X>> variable(SizeType rs, SizeType as, DegreeType d, const Vector<X>& x) {
        ARIADNE_ASSERT(x.size()==rs);
        Vector< Differential<X> > result(rs,as,d);
        for(SizeType i=0; i!=rs; ++i) { result[i]=x[i]; result[i][i]=X(1.0); }
        return result;
    }

    Vector<Differential<X>> affine(SizeType rs, SizeType as, DegreeType d, const Vector<X>& b, const Matrix<X>& A) {
        ARIADNE_ASSERT(b.size()==rs);
        ARIADNE_ASSERT(A.row_size()==rs);
        ARIADNE_ASSERT(A.column_size()==as);
        Vector< Differential<X> > result(rs,as,d);
        for(SizeType i=0; i!=rs; ++i) {
            result[i]=b[i];
            for(SizeType j=0; j!=as; ++j) {
                result[i][j]=A[i][j];
            }
        }
        return result;
    }

};


template<class X> template<class XX>
Vector<Differential<X>>::Vector(SizeType rs, SizeType as, DegreeType d, const XX* ptr)
    : _chars(as,d), _ary(rs,Differential<X>(as,d))
{
    for(SizeType i=0; i!=rs; ++i) {
        for(MultiIndex j(as); j.degree()<=d; ++j) {
            if(*ptr!=0) { (*this)[i][j]=*ptr; } ++ptr; } }
}

template<class X>
Vector<Differential<X>>::Vector(SizeType rs, SizeType as, DegreeType d,
                                  const Vector<X>& v, const Matrix<X>& A)
    :  _chars(as,d), _ary(rs,Differential<X>(as,d))
{
    ARIADNE_ASSERT(rs==v.size());
    ARIADNE_ASSERT(rs==A.row_size());
    ARIADNE_ASSERT(as==A.column_size());
    for(SizeType i=0; i!=this->result_size(); ++i) {
        (*this)[i]=v[i];
        for(SizeType j=0; j!=this->argument_size(); ++j) {
            const X& x=A[i][j];
            if(x!=0) { (*this)[i][j]=x; }
        }
    }
}

template<class X>
Vector<X>
value(const Vector<Differential<X>>& x)
{
    Vector<X> r(x.size());
    for(SizeType i=0; i!=x.size(); ++i) {
        r[i]=x[i].value();
    }
    return r;
}

template<class X>
Matrix<X>
jacobian(const Vector<Differential<X>>& x)
{
    if(x.size()==0) { return Matrix<X>(); }
    for(SizeType i=1; i!=x.size(); ++i) {
        ARIADNE_ASSERT(x[i].argument_size()==x[0].argument_size());
    }

    Matrix<X> r(x.size(),x[0].argument_size());
    for(SizeType i=0; i!=x.size(); ++i) {
        for(SizeType j=0; j!=x[0].argument_size(); ++j) {
            r[i][j]=x[i].gradient(j);

        }
    }
    return r;
}


template<class X> Vector<Differential<X>> compose(Vector<Differential<X>> const& x, Vector<Differential<X>> const& y) {
    Vector<Differential<X>> r(x.size(),y.zero_element());
    for(SizeType i=0; i!=r.size(); ++i) { r[i]=compose(x[i],y); }
    return r;
}

template<class X>
Vector<Differential<X>>
derivative(const Vector<Differential<X>>& x, SizeType k) {
    Vector<Differential<X>> r(x.size(), Differential<X>(x.argument_size(),min(x.degree(),1u)-1u));
    for(SizeType i=0; i!=r.size(); ++i) { r[i]=derivative(x[i],k); }
    return r;
}

template<class X>
Vector<Differential<X>>
antiderivative(const Vector<Differential<X>>& x, SizeType k) {
    Vector<Differential<X>> r(x.size(), Differential<X>(x.argument_size(),x.degree()+1));
    for(SizeType i=0; i!=r.size(); ++i) { r[i]=antiderivative(x[i],k); }
    return r;
}




template<class X>
Vector<Differential<X>>
lie_derivative(const Vector<Differential<X> >& df, const Vector<Differential<X> >& dg)
{
    Vector<Differential<X>> r(df.result_size(),df.argument_size(),df.degree()-1);
    Differential<X> t(df.argument_size(), df.degree()-1);
    MultiIndex a; X c;
    for(SizeType i=0; i!=df.result_size(); ++i) {
        Expansion<X> const& dfi_expansion = df[i].expansion();
        Expansion<X>& t_expansion = t.expansion();
        for(SizeType j=0; j!=df.argument_size(); ++j) {
            for(typename Expansion<X>::ConstIterator iter=dfi_expansion.begin(); iter!=dfi_expansion.end(); ++iter) {
                if(iter->key()[j]!=0) {
                    a=iter->key();
                    c=iter->data()*SizeType(a[j]);
                    a[j]-=1;
                    t_expansion.append(a,c);
                }
            }
            r[i]+=t*dg[j];
            t.clear();
        }
    }
    return r;
}




} //namespace Ariadne

#endif // ARIADNE_DIFFERENTIAL_H
