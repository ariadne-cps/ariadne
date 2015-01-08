/***************************************************************************
 *            differential.h
 *
 *  Copyright 2008  Pieter Collins
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
#include "numeric/float.decl.h"
#include "algebra/vector.h"
#include "algebra/matrix.h"
#include "algebra/multi_index.h"
#include "algebra/series.h"
#include "algebra/expansion.h"
#include <boost/concept_check.hpp>

namespace Ariadne {

class ExactInterval;

template<class X> class Vector;
template<class X> class Matrix;
template<class X> class Series;

template<class X> class Expansion;
template<class X> class Differential;
template<class X> class Vector<Differential<X>>;

typedef Differential<Float> FloatDifferential;
typedef Differential<ExactInterval> ExactIntervalDifferential;
typedef Differential<UpperInterval> UpperIntervalDifferential;
typedef Vector< Differential<Float> > FloatDifferentialVector;
typedef Vector< Differential<ExactInterval> > ExactIntervalDifferentialVector;
typedef Vector< Differential<UpperInterval> > UpperIntervalDifferentialVector;

template<class X> Differential<X>& operator+=(Differential<X>& x, const Differential<X>& y);
template<class X> Differential<X>& operator-=(Differential<X>& x, const Differential<X>& y);

template<class X, class R> Differential<X>& operator+=(Differential<X>& x, const R& c);
template<class X, class R> Differential<X>& operator-=(Differential<X>& x, const R& c);
template<class X, class R> Differential<X>& operator*=(Differential<X>& x, const R& c);
template<class X, class R> Differential<X>& operator/=(Differential<X>& x, const R& c);

template<class X> Differential<X> operator+(const Differential<X>& x);
template<class X> Differential<X> operator-(const Differential<X>& x);
template<class X> Differential<X> operator+(const Differential<X>& x, const Differential<X>& y);
template<class X> Differential<X> operator-(const Differential<X>& x, const Differential<X>& y);
template<class X> Differential<X> operator*(const Differential<X>& x, const Differential<X>& y);
template<class X> Differential<X> operator/(const Differential<X>& x, const Differential<X>& y);

template<class X> Differential<X> neg(const Differential<X>& x);
template<class X> Differential<X> rec(const Differential<X>& x);
template<class X> Differential<X> pow(const Differential<X>& x, int n);
template<class X> Differential<X> sqr(const Differential<X>& x);
template<class X> Differential<X> sqrt(const Differential<X>& x);
template<class X> Differential<X> exp(const Differential<X>& x);
template<class X> Differential<X> log(const Differential<X>& x);
template<class X> Differential<X> sin(const Differential<X>& x);
template<class X> Differential<X> cos(const Differential<X>& x);
template<class X> Differential<X> tan(const Differential<X>& x);

template<class X, class Y> Y evaluate(const Differential<X>& y, const Vector<Y>& z);
template<class X> Differential<X> compose(const Series<X>& x, const Differential<X>& y);
template<class X> Differential<X> derivative(const Differential<X>& x, uint i);
template<class X> Differential<X> antiderivative(const Differential<X>& x, uint i);

template<class X> Differential<X> compose(const Differential<X>&, const Vector< Differential<X> >&);
template<class X> Vector< Differential<X> > compose(const Vector< Differential<X> >&, const Vector< Differential<X> >&);



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
    static const uint MAX_DEGREE=65535;
    static const X _zero;
    static const X _one;
  public:
    //! \brief The type of used to represent numbers in the coefficient.
    typedef typename X::NumericType NumericType;
    //! \brief The type of used to represent the coefficients.
    typedef X ValueType;
    //! \brief The type of an Iterator through (index,coefficient) pairs..
    typedef typename Expansion<X>::Iterator Iterator;
    //! \brief The type of a constant Iterator through (index,coefficient) pairs..
    typedef typename Expansion<X>::ConstIterator ConstIterator;

    //! \brief Default constructor constructs a differential with degree zero in no variables.
    explicit Differential() : _expansion(0), _degree(0) { }
    //! \brief Constructs a differential with degree \a deg in \a as variables.
    explicit Differential(uint as, uint deg) : _expansion(as),_degree(deg) { }
    //! \brief Construct a differential from a mapping giving a coefficient for a finite number of multi-indices.
    explicit Differential(const std::map<MultiIndex,X>& map, uint deg) : _expansion(map), _degree(deg) { this->cleanup(); }
    //! \brief Construct a differential of degree \a deg from the power-series expansion \a e.
    //! Terms in \a e of degree higher than \a deg are truncated
    explicit Differential(const Expansion<X>& e, uint deg) : _expansion(e.argument_size()),_degree(deg) {
        for(typename Expansion<X>::ConstIterator iter=e.begin(); iter!=e.end(); ++iter) {
            if(iter->key().degree()<=deg) { this->_expansion.append(iter->key(),iter->data()); } } this->cleanup(); }
    //! \brief Construct a differential of degree \a deg from an initializer list list of (index,coefficient) pairs.
    explicit Differential(unsigned int as, unsigned int deg, InitializerList< std::pair<InitializerList<int>,X> > lst);

    //! \brief Construct a dense differential of degree \a deg in \a as variables from a list of coefficients beginning at \a ptr.
    template<class XX> Differential(uint as, uint deg, const XX* ptr) : _expansion(as), _degree(deg) {
        for(MultiIndex j(as); j.degree()<=deg; ++j) { XX const& x=*ptr; if(!decide(x==0)) { _expansion.append(j,*ptr); } ++ptr; } this->cleanup(); }
    //! \brief Conversion constructor from a different numerical type.
    template<class XX> Differential(const Differential<XX>& x)
        : _expansion(x.expansion()), _degree(x.degree()) { }

    //! \brief Set the differential equal to a constant, without changing the degree or number of arguments.
    Differential<X>& operator=(const X& c) {
        this->_expansion.clear(); this->_expansion.append(MultiIndex(this->argument_size()),c); return *this; }

    //! \brief A constant differential of degree \a deg in \a as arguments with value \a c.
    static Differential<X> constant(uint as, uint deg, const X& c) {
        Differential<X> r(as,deg); r._expansion.append(MultiIndex(as),c); return r; }
    //! \brief A differential of degree \a deg in \a as arguments representing the quantity \f$v+x_j\f$.
    static Differential<X> variable(uint as, uint deg, const X& v, uint j) {
        Differential<X> r(as,deg);
        MultiIndex a(as);
        r._expansion.append(a,v);
        a[j]=1; r._expansion.append(a,_one); return r; }

    //! \brief A vector of \a rs constant differentials of degree \a deg in \a as arguments with values \a c[j].
    static Vector< Differential<X> > constants(uint rs, uint as, uint deg, const Vector<X>& c) {
        ARIADNE_ASSERT(c.size()==rs);
        Vector< Differential<X> > result(rs,Differential(as,deg));
        for(uint i=0; i!=rs; ++i) { result[i]=c[i]; }
        return result;
    }
    //! \brief A vector of \a rs differentials of degree \a deg in \a as arguments with values \f$v_i+x_i\f$.
    //! \pre \a rs == \a as == c.size().
    static Vector< Differential<X> > variables(uint rs, uint as, uint deg, const Vector<X>& x) {
        ARIADNE_ASSERT(x.size()==rs);  ARIADNE_ASSERT(as==x.size());
        return variables(deg,x);
    }

    //! \brief \brief A vector of differentials of degree \a deg in \a as arguments with values \f$c_i+x_i\f$.
    static Vector< Differential<X> > variables(uint deg, const Vector<X>& x) {
        Vector< Differential<X> > result(x.size(),Differential<X>(x.size(),deg));
        MultiIndex a(x.size());
        for(uint i=0; i!=x.size(); ++i) {
            result[i]._expansion.append(a,x[i]);
            a[i]=1; result[i]._expansion.append(a,_one); a[i]=0;
        }
        return result;
    }

    //! \brief Equality operator. Tests equality of representation, so comparing two differentials which are mathematically equal may return false if the structural zeros are different.
    bool operator==(const Differential<X>& other) const;

    //! \brief Inequality operator.
    bool operator!=(const Differential<X>& other) const {
        return !(*this==other); }

    //! \brief The number of independent variables.
    uint argument_size() const { return this->_expansion.argument_size(); }
    //! \brief The maximum degree of the stored terms.
    uint degree() const { return this->_degree; }
    //! \brief The internal representation of the polynomial expansion.
    const Expansion<X>& expansion() const { return this->_expansion; }
    //! \brief A reference to the internal representation of the polynomial expansion.
    Expansion<X>& expansion() { return this->_expansion; }
    //! \brief The value of the differential i.e. the coefficient of \f$1\f$.
    const X& value() const { return this->operator[](MultiIndex(this->argument_size())); }
    //! \brief The coefficient of \f$x_j\f$.
    const X& gradient(uint j) const { return this->operator[](MultiIndex::unit(this->argument_size(),j)); }
    //! \brief The vector of coefficients of \f$x_j\f$.
    Vector<X> gradient() const { Vector<X> g(this->argument_size());
        for(uint j=0; j!=g.size(); ++j) { g[j]=this->gradient(j); } return g; }
    //! \brief The Hessian matrix.
    //! \note Note the the components of the Hessian matrix are \em half those of the values indexed by the differential.
    //! This is because the differential stores the coefficients of the Taylor expansion, rather than the derivatives themselves.
    Matrix<X> hessian() const;


    //! \brief A reference to the coefficient of \f$x_j\f$.
    X& operator[](const uint& j) {
        return this->operator[](MultiIndex::unit(this->argument_size(),j)); }
    //! \brief The coefficient of \f$x_j\f$.
    const X& operator[](const uint& j) const {
        return this->operator[](MultiIndex::unit(this->argument_size(),j)); }
    //! \brief The coefficient of \f$x^a\f$ i.e. \f$\prod x_j^{a_j}\f$.
    X& operator[](const MultiIndex& a) {
        ARIADNE_ASSERT_MSG(a.number_of_variables()==this->argument_size()," d="<<*this<<", a="<<a);
        return this->_expansion._at(a,GradedKeyLess()); }
    //! \brief A reference to the coefficient of \f$x^a\f$.
    const X& operator[](const MultiIndex& a) const {
        ARIADNE_ASSERT_MSG(a.number_of_variables()==this->argument_size()," d="<<*this<<", a="<<a);
        ConstIterator iter=this->_expansion.find(a);
        if(iter==this->_expansion.end()) { return _zero; } else { return iter->data(); } }

    //! \brief Set the degree to be equal to \a d.
    void set_degree(uint d) { this->_degree = d; }
    //! \brief Set the coefficient of the constant term to \a c.
    void set_value(const X& c) { this->operator[](MultiIndex(this->argument_size()))=c; }
    //! \brief Set the coefficient of the term in \f$x_j\f$ to \a d.
    void set_gradient(uint j, const X& d) { this->operator[](MultiIndex::unit(this->argument_size(),j))=d; }

    //! \brief An Iterator to the first structural nonzero.
    Iterator begin() { return this->_expansion.begin(); }
    //! \brief An Iterator to past-the-end structural nonzero.
    Iterator end() { return this->_expansion.end(); }
    //! \brief A constant Iterator to the first structural nonzero.
    ConstIterator begin() const { return this->_expansion.begin(); }
    //! \brief A constant Iterator to past-the-end structural nonzero.
    ConstIterator end() const { return this->_expansion.end(); }

    //! \brief The zero element of the differential algebra.
    Differential<X> create() const { return Differential<X>(this->argument_size(),this->degree()); }
    Differential<X> create_zero() const { return Differential<X>(this->argument_size(),this->degree()); }
    //! \brief Set all coefficients to zero.
    void clear() { this->_expansion.clear(); }
    //! \brief Remove all terms with coefficient \f$0\f$.
    void cleanup() { this->_expansion.graded_sort(); this->_expansion.combine_terms(); this->_expansion.remove_zeros(); }
    //! \brief Check that differential is sorted and all terms have degree less than maximum degree.
    void check() const;


    //! \brief Inplace addition.
    template<class XX> friend Differential<XX>& operator+=(Differential<XX>& x, const Differential<XX>& y);
    //! \brief Inplace subtraction.
    template<class XX> friend Differential<XX>& operator-=(Differential<XX>& x, const Differential<XX>& y);

    //! \brief Inplace addition of a constant.
    template<class XX, class RR> friend Differential<XX>& operator+=(Differential<XX>& x, const RR& c);
    //! \brief Inplace subtraction of a constant.
    template<class XX, class RR> friend Differential<XX>& operator-=(Differential<XX>& x, const RR& c);
    //! \brief Inplace multiplication by constant.
    template<class XX, class RR> friend Differential<XX>& operator*=(Differential<XX>& x, const RR& c);
    //! \brief Inplace division by a constant.
    template<class XX, class RR> friend Differential<XX>& operator/=(Differential<XX>& x, const RR& c);

    //! \brief Unary plus.
    friend Differential<X> operator+<>(const Differential<X>& x);
    //! \brief Unary minus.
    friend Differential<X> operator-<>(const Differential<X>& x);
    //! \brief Addition.
    friend Differential<X> operator+<>(const Differential<X>& x, const Differential<X>& y);
    //! \brief Subtraction.
    friend Differential<X> operator-<>(const Differential<X>& x, const Differential<X>& y);
    //! \brief Multiplication.
    friend Differential<X> operator*<>(const Differential<X>& x, const Differential<X>& y);
    //! \brief Division.
    friend Differential<X> operator/<>(const Differential<X>& x, const Differential<X>& y);

    //! \brief Negation.
    friend Differential<X> neg<>(const Differential<X>& x);
    //! \brief Reciprocal.
    friend Differential<X> rec<>(const Differential<X>& x);
    //! \brief Integer power.
    friend Differential<X> pow<>(const Differential<X>& x, int n);
    //! \brief Square.
    friend Differential<X> sqr<>(const Differential<X>& x);
    //! \brief Square root.
    friend Differential<X> sqrt<>(const Differential<X>& x);
    //! \brief Exponential function.
    friend Differential<X> exp<>(const Differential<X>& x);
    //! \brief Natural logarithm.
    friend Differential<X> log<>(const Differential<X>& x);
    //! \brief Sine function.
    friend Differential<X> sin<>(const Differential<X>& x);
    //! \brief Cosine function.
    friend Differential<X> cos<>(const Differential<X>& x);
    //! \brief Tangent function.
    friend Differential<X> tan<>(const Differential<X>& x);

#ifdef DOXYGEN
    //! \brief Compose by a power series in one variable.
    friend Differential<X> compose<>(const Series<X>& x, const Differential<X>& y);
    //! \brief Compose differentials at a point.
    friend Differential<X> compose<>(const Differential<X>& x, const Vector< Differential<X> >& y);
    //! \brief Compute the differential of the derivative.
    friend Differential<X> derivative<>(const Differential<X>& x, uint i);
    //! \brief Compute an antiderivative with respect to the variable \a i.
    friend Differential<X> antiderivative<>(const Differential<X>& x, uint i);
#endif
  private:
    Expansion<X> _expansion;
    uint _degree;
  private:
    //BOOST_CONCEPT_ASSERT((DifferentialConcept< Differential<X> >));
};

template<class X1, class X2> struct Arithmetic< X1,Differential<X2> > {
    typedef Differential<typename Arithmetic<X1,X2>::ResultType> ResultType; };
template<class X1, class X2> struct Arithmetic< Differential<X1>,X2 > {
    typedef Differential<typename Arithmetic<X1,X2>::ResultType> ResultType; };
template<class X1, class X2> struct Arithmetic< Differential<X1>,Differential<X2> > {
    typedef Differential<typename Arithmetic<X1,X2>::ResultType> ResultType; };

template<class X>
const X Differential<X>::_zero=X(0);

template<class X>
const X Differential<X>::_one=X(1);

template<class X>
Differential<X>::Differential(unsigned int as, unsigned int deg,
                              InitializerList< std::pair<InitializerList<int>,X> > lst)
    : _expansion(as,lst), _degree(deg)
{
    this->cleanup();
}

template<class X>
void
Differential<X>::check() const
{
    for(typename Differential<X>::ConstIterator iter=this->begin(); iter!=this->end(); ++iter) {
        ARIADNE_ASSERT_MSG(iter->key().degree()<=this->degree(), *this);
        typename Differential<X>::ConstIterator next = iter; ++next;
        ARIADNE_ASSERT_MSG(graded_less(iter->key(),next->key()),"Error in ordering Differential "<<this->expansion());
    }
}

/*
template<class X>
Differential<X>& Differential<X>::operator+=(const Differential<X>& x)
{
    for(ConstIterator iter=x._expansion.begin(); iter!=x._expansion.end(); ++iter) {
        this->_expansion[iter->key()]+=static_cast<const X&>(iter->data());
    }
    return *this;
}

template<class X>
Differential<X>& Differential<X>::operator-=(const Differential<X>& x)
{
    for(ConstIterator iter=x._expansion.begin(); iter!=x._expansion.end(); ++iter) {
        this->_expansion[iter->key()]-=static_cast<const X&>(iter->data());
    }
    return *this;
}

template<class X> template<class R>
Differential<X>& Differential<X>::operator+=(const R& c)
{
    this->_expansion[MultiIndex(this->argument_size())]+=c; return *this;
}

template<class X> template<class R>
Differential<X>& Differential<X>::operator-=(const R& c)
{
    this->_expansion[MultiIndex(this->argument_size())]-=c; return *this;
}

template<class X> template<class R>
Differential<X>& Differential<X>::operator*=(const R& c)
{
    if(c==0) {
        X zero=this->_expansion.begin()->data(); zero*=0;
        this->_expansion.clear();
        this->_expansion[MultiIndex(this->argument_size())]=zero;
    } else {
        for(Iterator iter=this->_expansion.begin(); iter!=this->_expansion.end(); ++iter) {
            static_cast<X&>(iter->data())*=c;
        }
    }
    return *this;
}


template<class X> template<class R>
Differential<X>& Differential<X>::operator/=(const R& c)
{
    for(Iterator iter=this->_expansion.begin(); iter!=this->_expansion.end(); ++iter) {
        static_cast<X&>(iter->data())/=c;
    }
    return *this;
}
*/


template<class X> Matrix<X> Differential<X>::hessian() const {
    ARIADNE_PRECONDITION(this->degree()>=2);
    Matrix<X> H(this->argument_size(),this->argument_size());
    ConstIterator iter=this->begin();
    while(iter!=this->end() && iter->key().degree()<=1) { ++iter; }
    uint i=0; uint j=1;
    while(iter!=this->end() && iter->key().degree()<=2) {
        const MultiIndex& a=iter->key(); const X& c=iter->data();
        while(a[i]==0) { ++i; j=i+1; }
        if(a[i]==2) { H[i][i]=c*2; }
        else { while(a[j]==0) { ++j; } H[i][j]=c*2; H[j][i]=c*2; }
        ++iter;
    }
    return H;
}

template<class X> bool Differential<X>::operator==(const Differential<X>& other) const {
    Differential<X> const& self=*this;
    if(self.argument_size()!=other.argument_size()) { return false; }
    ConstIterator self_iter=self.begin(); ConstIterator other_iter=other.begin();
    while(self_iter!=self.end() && other_iter!=other.end()) {
        if(self_iter->data()==0) {
            ++self_iter;
        } else if(other_iter->data()==0) {
            ++other_iter;
        } else {
            if(self_iter->key()!=other_iter->key() || self_iter->data()!=other_iter->data()) { return false; }
            ++self_iter; ++other_iter;
        }
    }
    while(self_iter!=self.end()) {
        if(self_iter->data()!=0) { return false; } ++self_iter;
    }
    while(other_iter!=other.end()) {
        if(other_iter->data()!=0) { return false; } ++other_iter;
    }
    return true;
}



template<class X> Differential<X> operator+(const Differential<X>& x, const Differential<X>& y);
template<class X> Differential<X> operator-(const Differential<X>& x, const Differential<X>& y);



template<class X>
Differential<X>& operator+=(Differential<X>& x, const Differential<X>& y)
{
    x=x+y;
    return x;
}

template<class X>
Differential<X>& operator-=(Differential<X>& x, const Differential<X>& y)
{
    x=x-y;
    return x;
}

template<class X, class R>
Differential<X>& operator+=(Differential<X>& x, const R& c)
{
    MultiIndex a(x.argument_size());
    if(x.expansion().empty()) {
        x.expansion().append(a,c);
    } else if(x.begin()->key()!=a) {
        x.expansion().prepend(a,c);
    } else {
        x.begin()->data()+=c;
    }
    return x;
}

template<class X, class R>
Differential<X>& operator-=(Differential<X>& x, const R& c)
{
    MultiIndex a(x.argument_size());
    if(x.expansion().empty()) {
        x.expansion().append(a,-c);
    } else if(x.begin()->key()!=a) {
        x.expansion().prepend(a,-c);
    } else {
        x.begin()->data()-=c;
    }
    return x;
}

template<class X, class R>
Differential<X>& operator*=(Differential<X>& x, const R& c)
{
    typedef typename Differential<X>::Iterator Iterator;
    if(c==static_cast<X>(0)) {
        x.clear();
    } else {
        for(Iterator iter=x.begin(); iter!=x.end(); ++iter) {
            static_cast<X&>(iter->data())*=c;
        }
    }
    return x;
}


template<class X, class R>
Differential<X>& operator/=(Differential<X>& x, const R& c)
{
    typedef typename Differential<X>::Iterator Iterator;
    for(Iterator iter=x.begin(); iter!=x.end(); ++iter) {
        static_cast<X&>(iter->data())/=static_cast<X>(c);
    }
    return x;
}


template<class X>
Differential<X> operator+(const Differential<X>& x)
{
    Differential<X> r(x.argument_size(),x.degree());
    r.expansion().reserve(x.expansion().number_of_nonzeros());
    for(typename Differential<X>::ConstIterator iter=x.begin(); iter!=x.end(); ++iter) {
        r.expansion().append(iter->key(), +iter->data());
    }
    return r;
}

template<class X>
Differential<X> operator-(const Differential<X>& x)
{
    Differential<X> r(x.argument_size(),x.degree());
    r.expansion().reserve(x.expansion().number_of_nonzeros());
    for(typename Differential<X>::ConstIterator iter=x.begin(); iter!=x.end(); ++iter) {
        r.expansion().append(iter->key(), -iter->data());
    }
    return r;
}


template<class X, class R>
EnableIfNumeric<R,Differential<X> >
operator+(const Differential<X>& x, const R& c)
{
    Differential<X> r(x); r+=X(c); return r;
}


template<class X, class R>
EnableIfNumeric<R,Differential<X> >
operator+(const R& c, const Differential<X>& x)
{
    Differential<X> r(x); r+=X(c); return r;
}

template<class X, class R>
EnableIfNumeric<R,Differential<X> >
operator-(const Differential<X>& x, const R& c)
{
    Differential<X> r(x); r-=X(c); return r;
}

template<class X, class R>
EnableIfNumeric<R,Differential<X> >
operator-(const R& c, const Differential<X>& x)
{
    Differential<X> r(-x); r+=X(c); return r;
}

template<class X, class R>
EnableIfNumeric<R,Differential<X> >
operator*(const Differential<X>& x, const R& c)
{
    Differential<X> r(x); r*=X(c); return r;
}

template<class X, class R>
EnableIfNumeric<R,Differential<X> >
operator*(const R& c, const Differential<X>& x)
{
    Differential<X> r(x); r*=X(c); return r;
}

template<class X, class R>
EnableIfNumeric<R,Differential<X> >
operator/(const Differential<X>& x, const R& c)
{
    Differential<X> r(x); r/=X(c); return r;
}

template<class X, class R>
EnableIfNumeric<R,Differential<X> >
operator/(const R& c, const Differential<X>& x)
{
    Differential<X> r(rec(x)); r*=X(c); return r;
}








template<class X>
Differential<X> operator+(const Differential<X>& x, const Differential<X>& y)
{
    ARIADNE_ASSERT_MSG(x.argument_size()==y.argument_size(),"x="<<x<<" y="<<y);
    Differential<X> r(x.argument_size(),std::min(x.degree(),y.degree()));
    typename Differential<X>::ConstIterator xiter=x.begin();
    typename Differential<X>::ConstIterator yiter=y.begin();
    // No need to check if maximum degree has been reached below,
    // since if one Iterator is above the maximum degree, the other is at the end.
    while(xiter!=x.end() && yiter!=y.end()) {
        if(xiter->key()==yiter->key()) {
            r.expansion().append(xiter->key(),xiter->data()+yiter->data());
            ++xiter; ++yiter;
        } else if(graded_less(xiter->key(),yiter->key())) {
            r.expansion().append(xiter->key(),xiter->data());
            ++xiter;
        } else {
            r.expansion().append(yiter->key(),yiter->data());
            ++yiter;
        }
    }
    while(xiter!=x.end() && xiter->key().degree()<=r.degree()) {
        r.expansion().append(xiter->key(),xiter->data());
        ++xiter;
    }
    while(yiter!=y.end() && yiter->key().degree()<=r.degree()) {
        r.expansion().append(yiter->key(),yiter->data());
        ++yiter;
    }
    //std::cerr<<"x="<<x<<" y="<<y<<" x+y="<<r<<"\n";
    return r;
}

template<class X>
Differential<X> operator-(const Differential<X>& x, const Differential<X>& y)
{
    ARIADNE_ASSERT_MSG(x.argument_size()==y.argument_size(),"x="<<x<<" y="<<y);
    Differential<X> r(x.argument_size(),std::min(x.degree(),y.degree()));
    typename Differential<X>::ConstIterator xiter=x.begin();
    typename Differential<X>::ConstIterator yiter=y.begin();
    // No need to check if maximum degree has been reached below,
    // since if one Iterator is above the maximum degree, the other is at the end.
    while(xiter!=x.end() && yiter!=y.end()) {
        if(xiter->key()==yiter->key()) {
            r.expansion().append(xiter->key(),xiter->data()-yiter->data());
            ++xiter; ++yiter;
        } else if(graded_less(xiter->key(),yiter->key())) {
            r.expansion().append(xiter->key(),xiter->data());
            ++xiter;
        } else {
            r.expansion().append(yiter->key(),-yiter->data());
            ++yiter;
        }
    }
    while(xiter!=x.end() && xiter->key().degree()<=r.degree()) {
        r.expansion().append(xiter->key(),xiter->data());
        ++xiter;
    }
    while(yiter!=y.end() && yiter->key().degree()<=r.degree()) {
        r.expansion().append(yiter->key(),-yiter->data());
        ++yiter;
    }
    //std::cerr<<"x="<<x<<" y="<<y<<" x-y="<<r<<"\n";
    return r;
}

template<class X>
Differential<X> operator*(const Differential<X>& x, const Differential<X>& y)
{
    typedef typename Differential<X>::ConstIterator ConstIterator;
    ARIADNE_ASSERT_MSG(x.argument_size()==y.argument_size(),"x="<<x<<" y="<<y);
    Differential<X> r(x.argument_size(),std::min(x._degree,y._degree));
    MultiIndex a(x.argument_size());
    X c(0);
    for(ConstIterator xiter=x._expansion.begin(); xiter!=x._expansion.end(); ++xiter) {
        if(xiter->key().degree()>r.degree()) { break; }
        for(ConstIterator yiter=y._expansion.begin(); yiter!=y._expansion.end(); ++yiter) {
            if(xiter->key().degree()+yiter->key().degree()>r.degree()) { break; }
            a=xiter->key()+yiter->key();
            c=static_cast<const X&>(xiter->data())*static_cast<const X&>(yiter->data());
            r._expansion.append(a,c);
        }
    }
    r.cleanup();
    //std::cerr<<"x="<<x<<" y="<<y<<" x*y="<<r<<"\n";
    return r;
}

template<class X>
Differential<X> operator/(const Differential<X>& x, const Differential<X>& y)
{
    ARIADNE_ASSERT_MSG(x.argument_size()==y.argument_size(),"x="<<x<<" y="<<y);
    return x*rec(y);
}

template<class X>
Differential<X>
min(const Differential<X>& x1, const Differential<X>& x2)
{
    ARIADNE_ASSERT_MSG(x1.argument_size()==x2.argument_size(),"x1="<<x1<<" x2="<<x2);
    if(x1.value()==x2.value()) {
        ARIADNE_THROW(std::runtime_error,"min(Differential<X> x1, Differential<X> x2)","x1[0]==x2[0]");
    }
    return x1.value()<x2.value() ? x1 : x2;
}


template<class X>
Differential<X>
max(const Differential<X>& x1,const Differential<X>& x2)
{
    ARIADNE_ASSERT_MSG(x1.argument_size()==x2.argument_size(),"x1="<<x1<<" x2="<<x2);
    if(x1.value()==x2.value()) {
        ARIADNE_THROW(std::runtime_error,"max(Differential<X> x1, Differential<X> x2)","x1[0]==x2[0]");
    }
    return x1.value()>x2.value() ? x1 : x2;
}

template<class X>
Differential<X>
abs(const Differential<X>& x)
{
    // FIXME: Maybe need different code for validated and approximate paradigms
    if(decide(x.value()==0)) {
        ARIADNE_THROW(std::runtime_error,"abs(Differential<X> x)","x[0]==0");
    }
    return decide(x.value()>0) ? pos(x) : neg(x);
}


template<class X>
Differential<X>
pos(const Differential<X>& x)
{
    return x;
}

template<class X>
Differential<X>
neg(const Differential<X>& x)
{
    return -x;
}

template<class X>
Differential<X> rec(const Differential<X>& x)
{
    return compose(Series<X>::rec(x.degree(),x.value()),x);
}

template<class X>
Differential<X> sqr(const Differential<X>& x)
{
    return pow(x,2);
}

template<class X>
Differential<X> pow(const Differential<X>& x, int n)
{
    return compose(Series<X>::pow(x.degree(),x.value(),n),x);
}

template<class X>
Differential<X> sqrt(const Differential<X>& x)
{
    return compose(Series<X>::sqrt(x.degree(),x.value()),x);
}

template<class X>
Differential<X> exp(const Differential<X>& x)
{
    return compose(Series<X>::exp(x.degree(),x.value()),x);
}

template<class X>
Differential<X> log(const Differential<X>& x)
{
    return compose(Series<X>::log(x.degree(),x.value()),x);
}

template<class X>
Differential<X> sin(const Differential<X>& x)
{
    return compose(Series<X>::sin(x.degree(),x.value()),x);
}

template<class X>
Differential<X> cos(const Differential<X>& x)
{
    return compose(Series<X>::cos(x.degree(),x.value()),x);
}

template<class X>
Differential<X> tan(const Differential<X>& x)
{
    return compose(Series<X>::tan(x.degree(),x.value()),x);
}



template<class X, class Y>
Y
evaluate(const Differential<X>& y, const Vector<Y>& x)
{
    //std::cerr<<ARIADNE_PRETTY_FUNCTION<<std::endl;
    ARIADNE_ASSERT_MSG(y.argument_size()==x.size(), "y="<<y<<" x="<<x);
    //std::cerr << "y=" << y << std::endl;
    //std::cerr << "x=" << x << std::endl;
    uint d=y.degree();
    uint ms=x.size();
    ARIADNE_ASSERT(d>=1);

    Y zero = x.zero_element();
    Y one = zero; one+=1;

    // Use inefficient brute-force approach with lots of storage...
    Array< Array< Y > > val(ms, Array< Y >(d+1));
    for(uint j=0; j!=ms; ++j) {
        val[j][0]=one;
        val[j][1]=x[j];
        for(uint k=2; k<=d; ++k) {
            val[j][k]=val[j][k-1]*x[j];
        }
    }

    Y r(zero);
    for(typename Differential<X>::ConstIterator iter=y.begin();
        iter!=y.end(); ++iter)
        {
            const MultiIndex& j=iter->key();
            const X& c=iter->data();
            Y t=one;
            for(uint k=0; k!=ms; ++k) {
                t=t*val[k][j[k]];
            }
            t*=c;
            r+=t;
            //std::cerr<<" j="<<j<<" c="<<c<<" r="<<r<<std::endl;
        }
    return r;
}


template<class X>
Differential<X> compose(const Series<X>& x, const Differential<X>& y)
{
    uint as=y.argument_size();
    uint d=std::min(x.degree(),y.degree());

    Differential<X> w=y;
    if(w.begin()->key().degree()==0) { w.begin()->data()=0; }
    Differential<X> r(as,d);
    r[MultiIndex(as)]=x[d];
    for(uint n=1; n<=d; ++n) {
        r=r*w;
        r+=x[d-n];
    }
    return r;
}


template<class X>
Vector<X>
gradient(const Differential<X>& x)
{
    Vector<X> r(x.argument_size());
    for(uint j=0; j!=x.argument_size(); ++j) {
        r[j]=x.gradient(j);
    }
    return r;
}


template<class X>
Differential<X> derivative(const Differential<X>& x, uint i)
{
    if(x.degree()==0) { return Differential<X>(x.argument_size(),0u); }
    Differential<X> r(x.argument_size(), x.degree()-1);
    MultiIndex a(x.argument_size());
    for(typename Differential<X>::ConstIterator iter=x.begin(); iter!=x.end(); ++iter) {
        a=iter->key();
        unsigned int n=a[i];
        if(n!=0) {
            const X& xc=x[a];
            --a[i];
            X& rc=r[a];
            rc=xc*n;
        }
    }
    return r;
}

template<class X>
Differential<X> antiderivative(const Differential<X>& x, uint i)
{
    Differential<X> r(x.argument_size(), x.degree()+1);
    MultiIndex a(x.argument_size());
    MultiIndex ra=MultiIndex(x.argument_size());
    MultiIndex ai=MultiIndex::unit(x.argument_size(),i);
    for(typename Differential<X>::ConstIterator iter=x.begin(); iter!=x.end(); ++iter) {
        a=iter->key();
        const X& xc=x[a];
        ++a[i];
        unsigned int n=a[i];
        X& rc=r[a];
        rc=xc/n;
    }
    return r;
}


template<class X>
OutputStream& operator<<(OutputStream& os, const Differential<X>& x)
{
    Expansion<X> e=x.expansion();
    //e.graded_sort();
    os << "SD("<<x.argument_size()<<","<<x.degree()<<"){";
    for(typename Expansion<X>::ConstIterator iter=e.begin(); iter!=e.end(); ++iter) {
        if(iter!=e.begin()) { os << ","; } os << " ";
        for(uint i=0; i!=e.argument_size(); ++i) {
            if(i!=0) { os << ","; }
            os << uint(iter->key()[i]);
        }
        os << ":" << X(iter->data());
    }
    return os << " }";
}


//! Embed the arguments in a space of dimension \a size, starting at position \a start.
template<class X>
Differential<X>
embed(uint before_size, const Differential<X>& x,
      uint after_size)
{
    return Differential<X>(x.degree(),x.expansion().embed(before_size,after_size));
}





template<class X>
struct NonAssignableDifferential
    : public Differential<X>
{
    NonAssignableDifferential<X>& operator=(const Differential<X>& other) {
        //ARIADNE_PRECONDITION(this->degree()==other.degree());
        //ARIADNE_PRECONDITION(this->argument_size()==other.argument_size());
        ARIADNE_ASSERT( this->argument_size()==other.argument_size()  );
        ARIADNE_ASSERT( this->degree()==other.degree() );
        this->Differential<X>::operator=(other); return *this;
    }
    NonAssignableDifferential<X>& operator=(const X& c) {
        this->Differential<X>::operator=(c); return *this;
    }
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

/*! \brief A class representing the derivatives of a vector quantity depending on multiple arguments. */
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
    typedef uint IndexType;
    // The value stored in the vector.
    typedef Differential<X> ValueType;
    // The type used for scalars.
    typedef Differential<X> ScalarType;
    // The type used for scalars.
    typedef X NumericType;

    Vector() : _chars(), _ary(0) { }
    Vector(uint rs) : _chars(), _ary(rs) { }
    Vector(uint rs, uint as, uint d) : _chars(as,d), _ary(rs) {
        for(uint i=0; i!=rs; ++i) { this->_ary[i]=Differential<X>(as,d); } }
    Vector(uint rs, const Differential<X>& sd) : _chars(sd), _ary(rs) {
        for(uint i=0; i!=rs; ++i) { this->_ary[i]=sd; } }
    Vector(uint rs, const Differential<X>* p) : _chars(), _ary(rs) {
        ARIADNE_ASSERT(rs>0); _chars=DifferentialCharacteristics(p[0]); for(uint i=0; i!=rs; ++i) { this->_ary[i]=p[i]; } }
    template<class XX> Vector(const Vector< Differential<XX> > dv) : _chars(dv._chars), _ary(dv._ary) { }
    template<class XX> Vector(uint rs, uint as, uint d, const XX* ptr);
    Vector(uint rs, uint as, uint d,const Vector<X>& v, const Matrix<X>& A);
    template<class E> Vector(const VectorExpression<E>& ve) : _ary(ve().size()) {
        for(uint i=0; i!=_ary.size(); ++i) { static_cast<Differential<X>&>(_ary[i])=ve()[i]; }
        ARIADNE_ASSERT(_ary.size()!=0); _chars=DifferentialCharacteristics(_ary[0]); }
    template<class E> Vector< Differential<X> >& operator=(const VectorExpression<E>& ve) {
        _ary.resize(ve().size()); for(uint i=0; i!=_ary.size(); ++i) { static_cast<Differential<X>&>(_ary[i])=ve()[i]; }
        ARIADNE_ASSERT(_ary.size()!=0); _chars=DifferentialCharacteristics(_ary[0]); return *this; }


    const Differential<X>& operator[](SizeType i) const { return this->_ary[i]; }
    NonAssignableDifferential<X>& operator[](SizeType i) { return static_cast<NonAssignableDifferential<X>&>(_ary[i]); }

    const Differential<X> zero_element() const { return Differential<X>(this->argument_size(),this->degree()); }

    NonAssignableDifferential<X>& at(SizeType i) { return static_cast<NonAssignableDifferential<X>&>(_ary[i]); }
    const Differential<X>& get(SizeType i) const { return this->_ary[i]; }
    void set(SizeType i, const Differential<X>& x) {
        ARIADNE_PRECONDITION(i<this->size());
        ARIADNE_PRECONDITION(this->argument_size()==x.argument_size());
        this->_ary[i]=x;
    }

    uint size() const { return this->_ary.size(); }
    uint result_size() const { return this->size(); }
    uint argument_size() const { return this->_chars.argument_size(); }
    uint degree() const { return this->_chars.degree(); }

    Vector<X> value() const {
        Vector<X> r(this->result_size());
        for(uint i=0; i!=r.size(); ++i) { r[i]=(*this)[i].value(); } return r; }
    Matrix<X> jacobian() const { Matrix<X> r(this->result_size(),this->argument_size());
        for(uint i=0; i!=r.row_size(); ++i) { for(uint j=0; j!=r.column_size(); ++j) { r[i][j]=(*this)[i].gradient(j); } } return r; }

    void set_value(const Vector<X>& c) {
        ARIADNE_ASSERT(this->result_size()==c.size());
        for(uint i=0; i!=c.size(); ++i) { (*this)[i].set_value(c[i]); } }

    static Vector< Differential<X> > constant(uint rs, uint as, uint d, const Vector<X>& c) {
        ARIADNE_ASSERT(c.size()==rs);
        Vector< Differential<X> > result(rs,as,d);
        for(uint i=0; i!=rs; ++i) { result[i]=c[i]; }
        return result;
    }

    static Vector< Differential<X> > variable(uint rs, uint as, uint d, const Vector<X>& x) {
        ARIADNE_ASSERT(x.size()==rs);
        Vector< Differential<X> > result(rs,as,d);
        for(uint i=0; i!=rs; ++i) { result[i]=x[i]; result[i][i]=X(1.0); }
        return result;
    }

    static Vector< Differential<X> > affine(uint rs, uint as, uint d, const Vector<X>& b, const Matrix<X>& A) {
        ARIADNE_ASSERT(b.size()==rs);
        ARIADNE_ASSERT(A.row_size()==rs);
        ARIADNE_ASSERT(A.column_size()==as);
        Vector< Differential<X> > result(rs,as,d);
        for(uint i=0; i!=rs; ++i) {
            result[i]=b[i];
            for(uint j=0; j!=as; ++j) {
                result[i][j]=A[i][j];
            }
        }
        return result;
    }

};


template<class X> template<class XX>
Vector< Differential<X> >::Vector(uint rs, uint as, uint d, const XX* ptr)
    : _chars(as,d), _ary(rs,Differential<X>(as,d))
{
    for(uint i=0; i!=rs; ++i) {
        for(MultiIndex j(as); j.degree()<=d; ++j) {
            if(*ptr!=0) { (*this)[i][j]=*ptr; } ++ptr; } }
}

template<class X>
Vector< Differential<X> >::Vector(uint rs, uint as, uint d,
                                  const Vector<X>& v, const Matrix<X>& A)
    :  _chars(as,d), _ary(rs,Differential<X>(as,d))
{
    ARIADNE_ASSERT(rs==v.size());
    ARIADNE_ASSERT(rs==A.row_size());
    ARIADNE_ASSERT(as==A.column_size());
    for(uint i=0; i!=this->result_size(); ++i) {
        (*this)[i]=v[i];
        for(uint j=0; j!=this->argument_size(); ++j) {
            const X& x=A[i][j];
            if(x!=0) { (*this)[i][j]=x; }
        }
    }
}


/*! \brief Compose the series of \a x and the series of \a y, assuming that \a x is centred at the value of \a y. The value of \a y is therefore unused by this method. */
template<class X, class Y>
Vector<Y>
evaluate(const Vector< Differential<X> >& x,
         const Vector<Y>& y)
{
    Vector<Y> r(x.result_size());
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=evaluate(x[i],y);
    }
    return r;
}

/*! \brief Compose the series of \a x and the series of \a y, assuming that \a x is centred at the value of \a y. The value of \a y is therefore unused by this method. */
template<class X>
Differential<X>
compose(const Differential<X>& x,
        const Vector< Differential<X> >& y)
{
    Vector< Differential<X> >& ync=const_cast< Vector< Differential<X> >&>(y);
    Vector<X> yv(y.size());
    for(uint i=0; i!=ync.result_size(); ++i) { yv[i]=ync[i].value(); ync[i].set_value(0); }
    Differential<X> r=evaluate(x,ync);
    for(uint i=0; i!=ync.result_size(); ++i) { ync[i].set_value(yv[i]); }
    return r;
}


/*! \brief Compose the series of \a x and the series of \a y, assuming that \a x is centred at the value of \a y. The value of \a y is therefore unused by this method. */
template<class X>
Vector< Differential<X> >
compose(const Vector< Differential<X> >& x,
        const Vector< Differential<X> >& y)
{
    ARIADNE_ASSERT(x.degree()==y.degree());
    //std::cerr<<"compose(DV x, DV y)\n x="<<x<<"\n y="<<y<<std::endl;
    Vector< Differential<X> >& ync=const_cast< Vector< Differential<X> >&>(y);
    Vector<X> yv(y.size());
    for(uint i=0; i!=ync.result_size(); ++i) { yv[i]=ync[i].value(); ync[i].set_value(0); }
    Vector< Differential<X> > r(x.size(),y.argument_size(),y.degree());
    for(uint i=0; i!=x.result_size(); ++i) { r[i]=evaluate(x[i],y); }
    for(uint i=0; i!=ync.result_size(); ++i) { ync[i].set_value(yv[i]); }
    return r;
}

template<class X>
Vector< Differential<X> >
lie_derivative(const Vector<Differential<X> >& df, const Vector<Differential<X> >& dg)
{
    Vector< Differential<X> > r(df.result_size(),df.argument_size(),df.degree()-1);
    Differential<X> t(df.argument_size(), df.degree()-1);
    for(uint i=0; i!=df.result_size(); ++i) {
        Expansion<X> const& dfi_expansion = df[i].expansion();
        Expansion<X>& t_expansion = t.expansion();
        for(uint j=0; j!=df.argument_size(); ++j) {
            for(typename Expansion<X>::ConstIterator iter=dfi_expansion.begin(); iter!=dfi_expansion.end(); ++iter) {
                if(iter->key()[j]!=0) {
                    t_expansion.append(iter->key(),iter->data());
                    t_expansion.back().data()*=t_expansion.back().key()[j];
                    t_expansion.back().key()[j]-=1;
                }
            }
            r[i]+=t*dg[j];
            t.clear();
        }
    }
    return r;
}

template<class X>
Vector< Differential<X> >
antiderivative(const Vector< Differential<X> >& x, uint k) {
    Vector< Differential<X> > r(x.size(), Differential<X>(x.argument_size(),x.degree()+1));
    for(uint i=0; i!=r.size(); ++i) { r[i]=antiderivative(x[i],k); }
    return r;
}

template<class X>
Vector<X>
value(const Vector< Differential<X> >& x)
{
    Vector<X> r(x.size());
    for(uint i=0; i!=x.size(); ++i) {
        r[i]=x[i].value();
    }
    return r;
}

template<class X>
Matrix<X>
jacobian(const Vector< Differential<X> >& x)
{
    if(x.size()==0) { return Matrix<X>(); }
    for(uint i=1; i!=x.size(); ++i) {
        ARIADNE_ASSERT(x[i].argument_size()==x[0].argument_size());
    }

    Matrix<X> r(x.size(),x[0].argument_size());
    for(uint i=0; i!=x.size(); ++i) {
        for(uint j=0; j!=x[0].argument_size(); ++j) {
            r[i][j]=x[i].gradient(j);

        }
    }
    return r;
}




} //namespace Ariadne

#endif // ARIADNE_DIFFERENTIAL_H
