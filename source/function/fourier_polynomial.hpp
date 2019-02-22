/***************************************************************************
 *            fourier_polynomial.hpp
 *
 *  Copyright 2008-18  Pieter Collins
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

/*! \file fourier_polynomial.hpp
 *  \brief Fourier functions on a bounded domain with a sparse representation.
 */

#ifndef ARIADNE_FOURIER_POLYNOMIAL_HPP
#define ARIADNE_FOURIER_POLYNOMIAL_HPP

#include <map>

#include "utility/macros.hpp"
#include "utility/typedefs.hpp"
#include "utility/array.hpp"
#include "utility/pointer.hpp"
#include "algebra/vector.hpp"
#include "algebra/covector.hpp"
#include "algebra/multi_index.hpp"
#include "algebra/expansion.hpp"
#include "algebra/operations.hpp"

namespace Ariadne {

struct IndexLess;
template<class X> class Complex;

template<class X> class UnivariateFourierPolynomial;
template<class X> class MultivariateFourierPolynomial;

template<class X> class OffsetArray {
    typedef std::ptrdiff_t IndexType;
    IndexType _l;
    IndexType _u;
    X* _p;
  public:
    ~OffsetArray() { _p+=_l; delete[] _p; _l=0; _p=nullptr; }
    OffsetArray(IndexType l, IndexType u, X const& z)
        : _l(l), _u(u), _p((new X[u-l])-l) { this->fill(z); }

    OffsetArray(OffsetArray<X> const& a) : _l(a._l), _u(a._u), _p((new X[_u-_l])-_l) {
        for (IndexType i=_l; i!=_u; ++i) { _p[i]=a._p[i]; } }
    OffsetArray(OffsetArray<X>&& a) : _l(a._l), _u(a._u), _p(a._p) { a._p=nullptr; }
    OffsetArray& operator=(OffsetArray<X> const& a) { if(this!=&a) { if (_u-_l != a._u-a._l) { delete[] _p-_l; _p=(new X[_u-_l])-_l; }
        this->_l=a._l; this->_u=a._u; for (IndexType i=_l; i!=_u; ++i) { _p[i]=a._p[i]; } } return *this; }
    OffsetArray& operator=(OffsetArray<X>&& a) { if(this!=&a) { this->_l=a._l; this->_u=a._u; this->_p=a._p; a._l=0; a._p=nullptr; } return *this; }

    SizeType size() const { return _u-_l; }
    IndexType lower_index() const { return _l; }
    IndexType upper_index() const { return _l; }
    X& operator[](IndexType i) { return _p[i]; }
    X const& operator[](IndexType i) const { return _p[i]; }
};






struct BiUniIndex {
    short int _a;
  public:
    BiUniIndex(short int a) : _a(a) { }
    operator short int() const { return _a; }
    DegreeType degree() const { return abs(_a); }
    SizeOne size() const { return SizeOne(); }
    BiUniIndex& operator--() { --_a; return *this; }
    BiUniIndex& operator++() { ++_a; return *this; }
};

SizeOne argument_size_of(UniformList<BiUniIndex> const&) { return SizeOne(); }

template<> struct IndexTraits<BiUniIndex> {
    typedef SizeOne SizeOfType;
    typedef IndexZero IndexIntoType;
    typedef short int InitializerType;
    typedef String NameType;
    template<class X> using Argument = Scalar<X>;
};

template<class X> class UnivariateFourierPolynomial;

//! \brief A \f$2\pi\f$ complex-valued periodic function defined by its Fourier coefficients
template<class X> class UnivariateFourierPolynomial<Complex<X>>
    : public DispatchAlgebraOperations<UnivariateFourierPolynomial<Complex<X>>,Complex<X>>
{
    typedef Complex<X> CX;
    typedef typename X::Paradigm P;
    typedef typename X::PrecisionType PR;
    SortedExpansion<BiUniIndex,Complex<X>,IndexLess> _terms;
  public:
    using FourierPolynomialType = UnivariateFourierPolynomial<Complex<X>>;

    typedef P Paradigm;
    typedef PR PrecisionType;

    typedef Complex<X> NumericType;

    UnivariateFourierPolynomial(PR pr) : _terms(SizeOne(),pr) { }
    UnivariateFourierPolynomial(BiUniIndex dmin, BiUniIndex dmax, std::function<CX(BiUniIndex)> const& g) : _terms(SizeOne(),nul(g(0))) {
        for(BiUniIndex i=dmin; i<=dmax; ++i) { _terms.append(i,g(i)); } }

    UnivariateFourierPolynomial(std::deque<CX> const& cs);

    static FourierPolynomialType constant(Number<P> y, PR pr);
    static FourierPolynomialType constant(CX const& c);
    static FourierPolynomialType basis(BiUniIndex k, PR pr);

    UnivariateFourierPolynomial<CX> create_constant(Number<P> y);

    DegreeType degree() const { return std::max(_terms.front().index().degree(),_terms.back().index().degree()); }
    PR precision() const { return _terms.zero_coefficient().precision(); }
    CX const& zero_coefficient() const { return _terms.zero_coefficient(); }
    CX const& operator[] (SizeType i) const { return _terms[i]; }
    Expansion<BiUniIndex,CX> const& terms() const { return this->_terms; }

    MagType<X> sup_norm() const;
    MagType<X> two_norm() const;
    friend MagType<X> sup_norm(FourierPolynomialType const& fp) { return fp.sup_norm(); }
    friend MagType<X> two_norm(FourierPolynomialType const& fp) { return fp.two_norm(); }

    friend decltype(auto) operator==(FourierPolynomialType const& fp1, FourierPolynomialType const& fp2) {
        return fp1._terms==fp2._terms; }

  private: // FIXME: Put these concrete-generic operations in proper place
    template<class Y, EnableIf<IsAssignable<X,Y>> =dummy>
        friend FourierPolynomialType operator+(FourierPolynomialType p, const Y& c) {
            CX xc=p.zero_coefficient(); xc=c; return p+xc; }
    template<class Y, EnableIf<IsAssignable<X,Y>> =dummy>
        friend FourierPolynomialType operator-(FourierPolynomialType p, const Y& c) {
            CX xc=p.zero_coefficient(); xc=c; return p-xc; }
    template<class Y, EnableIf<IsAssignable<X,Y>> =dummy>
        friend FourierPolynomialType operator*(const Y& c, FourierPolynomialType p) {
            CX xc=p.zero_coefficient(); xc=c; return xc*p; }
    template<class Y, EnableIf<IsAssignable<X,Y>> =dummy>
        friend FourierPolynomialType operator*(FourierPolynomialType p, const Y& c) {
            CX xc=p.zero_coefficient(); xc=c; return p*xc; }
    template<class Y, EnableIf<IsAssignable<X,Y>> =dummy>
        friend FourierPolynomialType operator/(FourierPolynomialType p, const Y& c) {
            CX xc=p.zero_coefficient(); xc=c; return p/xc; }

  public:
    friend OutputStream& operator<<(OutputStream& os, FourierPolynomialType const& fp) { return fp._write(os); }

    CX operator() (X const& x) const;
    template<class Y> friend decltype(auto) evaluate(FourierPolynomialType const& p, Y const& y) { return _evaluate(p,y); }
  private: public:
    static FourierPolynomialType apply(Neg, FourierPolynomialType fp);
    static FourierPolynomialType apply(Add, FourierPolynomialType const& fp1, FourierPolynomialType const& fp2);
    static FourierPolynomialType apply(Sub, FourierPolynomialType const& fp1, FourierPolynomialType const& fp2);
    static FourierPolynomialType apply(Mul, FourierPolynomialType const& fp1, FourierPolynomialType const& fp2);
    static FourierPolynomialType apply(Add, FourierPolynomialType fp1, NumericType const& s2);
    static FourierPolynomialType apply(Mul, FourierPolynomialType fp1, NumericType const& s2);
    static FourierPolynomialType apply(Sqr, FourierPolynomialType fp);

    OutputStream& _write(OutputStream& os) const;

  public:
    template<class Y> static auto _evaluate(FourierPolynomialType const& f, Y const& v) -> ArithmeticType<Complex<X>,Y>;
};

template<class X> struct AlgebraOperations<UnivariateFourierPolynomial<X>,X> : public UnivariateFourierPolynomial<X> { };


template<class X> template<class Y> auto
UnivariateFourierPolynomial<Complex<X>>::_evaluate(UnivariateFourierPolynomial<Complex<X>> const& f, Y const& x) -> ArithmeticType<Complex<X>,Y> {
    typedef ArithmeticType<Complex<X>,Y> R;
    R r=nul(x);
    Complex<X> i(0,1);
    for (auto iter=f._terms.begin(); iter!=f._terms.end(); ++iter) {
        Int k=iter->index(); r+=iter->coefficient()*exp(i*k*x);
    }
    return r;
}


/*

template<class X> class MultivariateFourierPolynomial
    : public DispatchAlgebraOperations<MultivariateFourierPolynomial<X>,X>
{
    typedef typename X::Paradigm P;
    typedef typename X::PrecisionType PR;
    Expansion<MultiIndex,X> _terms;
  public:
    using FourierPolynomialType = MultivariateFourierPolynomial<X>;

    typedef P Paradigm;
    typedef PR PrecisionType;
    typedef X NumericType;
    typedef typename Expansion<MultiIndex,X>::ConstIterator ConstIterator;

    MultivariateFourierPolynomial(SizeType as, PR pr) : _terms(as,X(pr)) { }

    static FourierPolynomialType constant(SizeType as, Number<P> y, PR pr);
    static FourierPolynomialType constant(SizeType as, X const& c);
    static FourierPolynomialType coordinate(SizeType as, SizeType i, PR pr);
    static FourierPolynomialType basis(SizeType as, SizeType i, SizeType k, PR pr);

    FourierPolynomialType create_constant(Number<P> y);

    SizeType argument_size() const { return this->_terms.argument_size(); }
    PR precision() const { return this->_terms.zero_coefficient().precision(); }
    X const& zero_coefficient() const { return _terms.zero_coefficient(); }

    MagType<X> sup_norm() const;
    friend MagType<X> sup_norm(FourierPolynomialType const& fp) { return fp.sup_norm(); }

    friend FourierPolynomialType operator*(Number<P> const& s1, FourierPolynomialType fp2) { return X(s1,fp2.precision()) * fp2; }
    friend OutputStream& operator<<(OutputStream& os, FourierPolynomialType const& fp) { return fp._write(os); }

    X operator() (Vector<X> const& x) const;
    friend X evaluate(FourierPolynomialType const& f, Vector<X> const& x) { return f(x); }
  private: public:
    static FourierPolynomialType apply(Neg, FourierPolynomialType fp);
    static FourierPolynomialType apply(Sqr, FourierPolynomialType fp);
    static FourierPolynomialType apply(Add, FourierPolynomialType const& fp1, FourierPolynomialType const& fp2);
    static FourierPolynomialType apply(Sub, FourierPolynomialType const& fp1, FourierPolynomialType const& fp2);
    static FourierPolynomialType apply(Mul, FourierPolynomialType const& fp1, FourierPolynomialType const& fp2);
    static FourierPolynomialType apply(Add, FourierPolynomialType fp1, Scalar<X> const& s2);
    static FourierPolynomialType apply(Mul, FourierPolynomialType fp1, Scalar<X> const& s2);
    OutputStream& _write(OutputStream& os) const;

    static X _eval(ConstIterator& from, ConstIterator end, Vector<X> const& s, SizeType k);
    static FourierPolynomialType _mul_from(ConstIterator& from1, ConstIterator end1, ConstIterator& from2, ConstIterator end2, SizeType k);
};

template<class X> struct AlgebraOperations<MultivariateFourierPolynomial<X>,X> : public MultivariateFourierPolynomial<X> { };

*/

} // namespace Ariadne

#endif // ARIADNE_FOURIER_POLYNOMIAL_HPP
