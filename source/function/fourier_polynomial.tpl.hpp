/***************************************************************************
 *            fourier_polynomial.tpl.hpp
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
 *  GNU Library General Public License for more detai1ls.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "numeric/numeric.hpp"

#include <iomanip>
#include <limits>

#include "numeric/rounding.hpp"
#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"
#include "algebra/covector.hpp"
#include "algebra/matrix.hpp"
#include "algebra/expansion.hpp"
#include "algebra/series.hpp"
#include "algebra/differential.hpp"
#include "function/fourier_polynomial.hpp"
#include "function/taylor_series.hpp"
#include "function/function.hpp"
#include "utility/exceptions.hpp"

#include "algebra/expansion.inl.hpp"
#include "algebra/expansion.tpl.hpp"
#include "algebra/evaluate.tpl.hpp"
#include "algebra/algebra_operations.tpl.hpp"

#include "algebra/multi_index-noaliasing.hpp"
#include "function/function_mixin.hpp"
#include "algebra/vector.hpp"

namespace Ariadne {

/*
FloatDP operator+(FloatDP x1, FloatDP x2);
FloatDP operator-(FloatDP x1, FloatDP x2);
FloatDP operator*(FloatDP x1, FloatDP x2);
FloatDP operator/(FloatDP x1, FloatDP x2);
FloatDP& operator+=(FloatDP& x1, FloatDP x2);
FloatDP& operator-=(FloatDP& x1, FloatDP x2);
FloatDP& operator*=(FloatDP& x1, FloatDP x2);
FloatDP& operator/=(FloatDP& x1, FloatDP x2);
FloatMP operator+(FloatMP const& x1, FloatMP const& x2);
FloatMP operator-(FloatMP const& x1, FloatMP const& x2);
FloatMP operator*(FloatMP const& x1, FloatMP const& x2);
FloatMP operator/(FloatMP const& x1, FloatMP const& x2);
FloatMP& operator+=(FloatMP& x1, FloatMP const& x2);
FloatMP& operator-=(FloatMP& x1, FloatMP const& x2);
FloatMP& operator*=(FloatMP& x1, FloatMP const& x2);
FloatMP& operator/=(FloatMP& x1, FloatMP const& x2);
*/

template<class X> auto UnivariateFourierPolynomial<Complex<X>>::constant(Complex<X> const& x) -> UnivariateFourierPolynomial<Complex<X>> {
    UnivariateFourierPolynomial<Complex<X>> fpr(x.precision());
    fpr._terms.append(0u,x);
    return fpr;
}

template<class X> auto UnivariateFourierPolynomial<Complex<X>>::constant(Number<P> c, PR pr) -> UnivariateFourierPolynomial<Complex<X>> {
    UnivariateFourierPolynomial<Complex<X>> fpr(pr);
    fpr._terms.append(0u,X(c,pr));
    return fpr;
}

template<class X> auto UnivariateFourierPolynomial<Complex<X>>::basis(BiUniIndex k, PR pr) -> UnivariateFourierPolynomial<Complex<X>> {
    UnivariateFourierPolynomial<Complex<X>> fpr(pr);
    fpr._terms.append(k,X(1,pr));
    return fpr;
}


template<class X> auto UnivariateFourierPolynomial<Complex<X>>::create_constant(Number<P> c) -> UnivariateFourierPolynomial<Complex<X>> {
    return UnivariateFourierPolynomial<Complex<X>>::constant(c,this->precision());
}


template<class X> auto
UnivariateFourierPolynomial<Complex<X>>::two_norm() const -> MagType<X> {
    MagType<X> r = mag(this->zero_coefficient());
    for(auto term : this->_terms) { r+=sqr(mag(term.coefficient())); }
    return sqrt(r);
}

template<class X> auto
UnivariateFourierPolynomial<Complex<X>>::sup_norm() const -> MagType<X> {
    MagType<X> r = mag(this->zero_coefficient());
    for(auto term : this->_terms) { r+=mag(term.coefficient()); }
    return r;
}


template<class X> auto
UnivariateFourierPolynomial<Complex<X>>::apply(Neg, UnivariateFourierPolynomial<Complex<X>> fp) -> UnivariateFourierPolynomial<Complex<X>> {
    for(auto term : fp._terms) { term.coefficient() = -term.coefficient(); }
    return fp;
}

template<class X> auto
UnivariateFourierPolynomial<Complex<X>>::apply(Add, UnivariateFourierPolynomial<Complex<X>> fp, NumericType const& s) -> UnivariateFourierPolynomial<Complex<X>> {
    auto iter = fp._terms.begin();
    while (iter!=fp._terms.end()) {
        if (iter->index()==0) { iter->coefficient()+=s; return fp; }
        if (iter->index()>0) { fp._terms.insert(BiUniIndex(0),s); return fp; }
        ++iter;
    }
   fp._terms.append(0u,s);
   return fp;
}

template<class X> auto
UnivariateFourierPolynomial<Complex<X>>::apply(Mul, UnivariateFourierPolynomial<Complex<X>> fp, NumericType const& s) -> UnivariateFourierPolynomial<Complex<X>> {
    for(auto term : fp._terms) { term.coefficient() *= s; }
    return fp;
}

template<class X> auto
UnivariateFourierPolynomial<Complex<X>>::apply(Sqr, UnivariateFourierPolynomial<Complex<X>> fp) -> UnivariateFourierPolynomial<Complex<X>> {
    return fp*fp;
}

/*
template<class X> template<class Y, EnableIf<IsConvertible<SumType<X,Y>,X>>> auto
UnivariateFourierPolynomial<Complex<X>>::_add(UnivariateFourierPolynomial<Complex<X>> fp, Y const& s) -> UnivariateFourierPolynomial<Complex<X>> {
    if (fp._terms.empty()) { fp._terms.append(0u,fp._terms.zero_coefficient()+s); }
    else if (fp._terms.front().index()==0) { fp._terms.front().coefficient()+=s; }
    else { assert(false); }
    return fp;
}

template<class X> template<class Y, EnableIf<IsConvertible<DifferenceType<X,Y>,X>>> auto
UnivariateFourierPolynomial<Complex<X>>::_sub(UnivariateFourierPolynomial<Complex<X>> fp, Y const& s) -> UnivariateFourierPolynomial<Complex<X>> {
    if (fp._terms.empty()) { fp._terms.append(0u,fp._terms.zero_coefficient()-s); }
    else if (fp._terms.front().index()==0) { fp._terms.front().coefficient()-=s; }
    else { assert(false); }
    return fp;
}

template<class X> template<class Y, EnableIf<IsConvertible<ProductType<X,Y>,X>>> auto
UnivariateFourierPolynomial<Complex<X>>::_mul(UnivariateFourierPolynomial<Complex<X>> fp, Y const& s) -> UnivariateFourierPolynomial<Complex<X>> {
    for(auto term : fp._terms) { term.coefficient() *= s; }
    return fp;
}

template<class X> template<class Y, EnableIf<IsConvertible<QuotientType<X,Y>,X>>> auto
UnivariateFourierPolynomial<Complex<X>>::_div(UnivariateFourierPolynomial<Complex<X>> fp, Y const& s) -> UnivariateFourierPolynomial<Complex<X>> {
    for(auto term : fp._terms) { term.coefficient() = term.coefficient()/s; }
    return fp;
}
*/

template<class X> auto
UnivariateFourierPolynomial<Complex<X>>::apply(Add, UnivariateFourierPolynomial<Complex<X>> const& fp1, UnivariateFourierPolynomial<Complex<X>> const& fp2) -> UnivariateFourierPolynomial<Complex<X>> {
    UnivariateFourierPolynomial<Complex<X>> fpr(min(fp1.precision(),fp2.precision()));
    auto iter1=fp1._terms.begin();
    auto iter2=fp2._terms.begin();
    while(iter1!=fp1._terms.end() && iter2!=fp2._terms.end()) {
        if(iter1->index()==iter2->index()) {
            fpr._terms.append(iter1->index(),iter1->coefficient()+iter2->coefficient());
            ++iter1; ++iter2;
        } else if(iter1->index()<iter2->index()) {
            fpr._terms.append(iter1->index(),iter1->coefficient());
            ++iter1;
        } else { // iter1->index()>iter2->index()
            fpr._terms.append(iter2->index(),iter2->coefficient());
            ++iter2;
        }
    }
    while(iter1!=fp1._terms.end()) {
        fpr._terms.append(iter1->index(),iter1->coefficient());
        ++iter1;
    }
    while(iter2!=fp2._terms.end()) {
        fpr._terms.append(iter2->index(),iter2->coefficient());
        ++iter2;
    }
    return fpr;
}

template<class X> auto
UnivariateFourierPolynomial<Complex<X>>::apply(Sub, UnivariateFourierPolynomial<Complex<X>> const& fp1, UnivariateFourierPolynomial<Complex<X>> const& fp2) -> UnivariateFourierPolynomial<Complex<X>> {
    return add(fp1,neg(fp2));
}

inline short int cast_signed(unsigned short int m) { return static_cast<short int>(m); }
inline int cast_signed(unsigned int m) { return static_cast<int>(m); }

inline unsigned int cast_unsigned(int n) { return static_cast<unsigned int>(n); }
inline unsigned short int cast_unsigned(short int n) { return static_cast<unsigned short int>(n); }


template<class I, class X> Map<I,X> add(Map<I,X> const& a1, Map<I,X> const& a2) {
    Map<I,X> r;
    auto iter1 = a1.begin(); auto iter2=a2.begin();
    while (iter1!=a1.end() && iter2!=a2.end()) {
        if (iter1.index()==iter2.index()) {
            r[iter1.index()] = iter1.coefficient()+iter2.coefficient();
            ++iter1; ++iter2;
        } else if (iter1.index()<iter2.index()) {
            r[iter1.index()] = iter1.coefficient();
            ++iter1;
        } else {
            r[iter2.index()] = iter2.coefficient();
            ++iter2;
        }
    }
    while (iter1!=a1.end()) {
        r[iter1.index()] = iter1.coefficient();
        ++iter1;
    }
    while (iter2!=a2.end()) {
        r[iter2.index()] = iter2.coefficient();
        ++iter2;
    }
    return r;
}

template<class I, class X> Map<I,X> mul(Map<I,X> const& a1, Map<I,X> const& a2) {
    Map<I,X> r;
    for (auto t1 : a1) {
        for (auto t2 : a2) {
            r[t1.index()+t2.index()] += t1.coefficient()*t2.coefficient();
        }
    }
    return r;
}

template<class X> auto
UnivariateFourierPolynomial<Complex<X>>::apply(Mul, UnivariateFourierPolynomial<Complex<X>> const& fp1, UnivariateFourierPolynomial<Complex<X>> const& fp2) -> UnivariateFourierPolynomial<Complex<X>> {
    typedef short int I; typedef unsigned short int J;
    PR pr=min(fp1.precision(),fp2.precision());
    X z(pr);
    I l1=fp1._terms.front().index(); I u1=fp1._terms.back().index();
    I l2=fp2._terms.front().index(); I u2=fp2._terms.back().index();
    I l=l1+l2;
    J n1=static_cast<J>(u1-l1); J n2=static_cast<J>(u2-l2); J n=n1+n2+1u;
    assert(l1<=u1); assert(l2<=u2);
    Array<CX> buf(n,z);
    for (auto term1 : fp1._terms) {
        for (auto term2 : fp2._terms) {
            CX coefficient_product=term1.coefficient()*term2.coefficient();
            buf[static_cast<SizeType>((term1.index()+term2.index())-l)]+=coefficient_product;
        }
    }
    UnivariateFourierPolynomial<Complex<X>> fpr(pr);
    for (SizeType i=0; i!=buf.size(); ++i) {
        fpr._terms.append((J)i+l,buf[i]);
    }
    return fpr;
}

template<class X> auto
UnivariateFourierPolynomial<Complex<X>>::operator() (X const& x) const -> Complex<X> {
    return _evaluate(*this,x);
}

template<class X> auto
UnivariateFourierPolynomial<Complex<X>>::_write (OutputStream& os) const -> OutputStream& {
    for (auto term : this->_terms) {
        String c = to_str(term.coefficient());
        if (c[0]!='+' && c[0]!='-') { os << "+"; }
        os << "(" << c << ")" << "*exp(i*" << term.index() << "*x)";
    }
    return os;
}




/*

template<class X> auto MultivariateFourierPolynomial<X>::constant(SizeType as, X const& c) -> MultivariateFourierPolynomial<X> {
    MultivariateFourierPolynomial<X> fpr(as,c.precision());
    fpr._terms.append(MultiIndex::zero(as),c);
    return fpr;
}

template<class X> auto MultivariateFourierPolynomial<X>::constant(SizeType as, Number<P> c, PR pr) -> MultivariateFourierPolynomial<X> {
    MultivariateFourierPolynomial<X> fpr(as,pr);
    fpr._terms.append(MultiIndex::zero(as),X(c,pr));
    return fpr;
}

template<class X> auto MultivariateFourierPolynomial<X>::coordinate(SizeType as, SizeType i, PR pr) -> MultivariateFourierPolynomial<X> {
    MultivariateFourierPolynomial<X> fpr(as,pr);
    fpr._terms.append(MultiIndex::unit(as,i),X(1,pr));
    return fpr;
}

template<class X> auto MultivariateFourierPolynomial<X>::basis(SizeType as, SizeType i, SizeType k, PR pr) -> MultivariateFourierPolynomial<X> {
    MultivariateFourierPolynomial<X> fpr(as,pr);
    fpr._terms.append(MultiIndex::unit(as,i)*k,X(1,pr));
    return fpr;
}


template<class X> auto MultivariateFourierPolynomial<X>::create_constant(Number<P> c) -> MultivariateFourierPolynomial<X> {
    return MultivariateFourierPolynomial<X>::constant(this->argument_size(),c,this->precision());
}


template<class X> auto
MultivariateFourierPolynomial<X>::sup_norm() const -> MagType<X> {
    MagType<X> r = mag(this->zero_coefficient());
    for(auto term : this->_terms) { r+=mag(term.coefficient()); }
    return r;
}

template<class X> auto
MultivariateFourierPolynomial<X>::apply(Neg, MultivariateFourierPolynomial<X> fp) -> MultivariateFourierPolynomial<X> {
    for(auto term : fp._terms) { term.coefficient() = -term.coefficient(); }
    return fp;
}

template<class X> auto
MultivariateFourierPolynomial<X>::apply(Add, MultivariateFourierPolynomial<X> fp, Scalar<X> const& s) -> MultivariateFourierPolynomial<X> {
    if (fp._terms.empty()) { fp._terms.append(MultiIndex::zero(fp.argument_size()),s); }
    else if (fp._terms.front().index().degree()==0) { fp._terms.front().coefficient()+=s; }
    else { assert(false); }
    return fp;
}

template<class X> auto
MultivariateFourierPolynomial<X>::apply(Mul, MultivariateFourierPolynomial<X> fp, Scalar<X> const& s) -> MultivariateFourierPolynomial<X> {
    for(auto term : fp._terms) { term.coefficient() *= s; }
    return fp;
}

template<class X> auto
MultivariateFourierPolynomial<X>::apply(Add, MultivariateFourierPolynomial<X> const& fp1, MultivariateFourierPolynomial<X> const& fp2) -> MultivariateFourierPolynomial<X> {
    ARIADNE_ASSERT(fp1.argument_size()==fp2.argument_size());
    MultivariateFourierPolynomial<X> fpr(fp1.argument_size(),min(fp1.precision(),fp2.precision()));
    auto iter1=fp1._terms.begin();
    auto iter2=fp2._terms.begin();
    while(iter1!=fp1._terms.end() && iter2!=fp2._terms.end()) {
        if(iter1->index()==iter2->index()) {
            fpr._terms.append(iter1->index(),iter1->coefficient()+iter2->coefficient());
            ++iter1; ++iter2;
        } else if(reverse_lexicographic_less(iter1->index(),iter2->index())) {
            fpr._terms.append(iter1->index(),iter1->coefficient());
            ++iter1;
        } else { // iter1->index()>iter2->index()
            fpr._terms.append(iter2->index(),iter2->coefficient());
            ++iter2;
        }
    }
    while(iter1!=fp1._terms.end()) {
        fpr._terms.append(iter1->index(),iter1->coefficient());
        ++iter1;
    }
    while(iter2!=fp2._terms.end()) {
        fpr._terms.append(iter2->index(),iter2->coefficient());
        ++iter2;
    }
    return fpr;
}

template<class X> auto
MultivariateFourierPolynomial<X>::apply(Sub, MultivariateFourierPolynomial<X> const& fp1, MultivariateFourierPolynomial<X> const& fp2) -> MultivariateFourierPolynomial<X> {
    return add(fp1,neg(fp2));
}

template<class X> auto
MultivariateFourierPolynomial<X>::apply(Mul, MultivariateFourierPolynomial<X> const& fp1, MultivariateFourierPolynomial<X> const& fp2) -> MultivariateFourierPolynomial<X> {
    ConstIterator iter1=fp1._terms.begin(); ConstIterator iter2=fp2._terms.begin();
    ConstIterator end1=fp1._terms.end(); ConstIterator end2=fp2._terms.end();
    return _mul_from(iter1,end1,iter2,end2,0u);
}


template<class X> auto
MultivariateFourierPolynomial<X>::apply(Sqr, MultivariateFourierPolynomial<X> fp) -> MultivariateFourierPolynomial<X> {
    return fp*fp;
}

template<class X> struct PartialFourierPolynomial {
    typedef typename MultivariateFourierPolynomial<X>::ConstIterator ConstIterator;
    ConstIterator _from, _to; SizeType _k;
    PartialFourierPolynomial(ConstIterator from, SizeType k, ConstIterator end)
        : _from(from), _to(from), _k(k) { while(_to!=end && equal_up_to(_from->index(),_to->index(),k)) { ++_to; } }
};

template<class X> class FourierValues {
    mutable List<X> _lst;
  public:
    FourierValues(X const& x) : _lst() { _lst.reserve(16); _lst.append(X(1,x.precision())); _lst.append(x); }
    X const& operator[] (SizeType i) const {
        for(SizeType j=_lst.size(); j<=i; ++j) { _lst.append(2*_lst[1]*_lst[j-1]-_lst[j-2]); } return _lst[i]; }
};

inline Boolean equal_up_to(MultiIndex const& a1, MultiIndex const& a2, SizeType k) {
    for (SizeType i=0; i!=k; ++i) { if (a1[i]!=a2[i]) { return false; } } return true;
}

template<class X> auto
MultivariateFourierPolynomial<X>::operator() (Vector<Scalar<X>> const& x) const -> Scalar<X> {
    ConstIterator from=this->_terms.begin();
    return _eval(from,this->_terms.end(),x,0u);
}

template<class X> auto
MultivariateFourierPolynomial<X>::_eval(ConstIterator& from, ConstIterator end, Vector<Scalar<X>> const& x, SizeType k) -> Scalar<X> {
    FourierValues<X> v(x[k]);
    MultiIndex a_front=from->index();
    X r(x.zero_element());
    if(k+1u==x.size()) {
        while (from!=end && equal_up_to(a_front,from->index(),k)) {
            r+=from->coefficient()*v[from->index()[k]];
            ++from;
        }
    } else {
        while (from!=end && equal_up_to(a_front,from->index(),k)) {
            auto a=from->index();
            auto c=_eval(from,end,x,k+1);
            r+=c*v[a[k]];
        }
    }
    return r;
}



template<class X> auto
MultivariateFourierPolynomial<X>::_mul_from(ConstIterator& from1, ConstIterator end1, ConstIterator& from2, ConstIterator end2, SizeType k) -> MultivariateFourierPolynomial<X> {
    MultiIndexConstReference a1=from1->index();
    MultiIndexConstReference a2=from2->index();
    SizeType n=a1.size();
    PR pr=max(from1->coefficient().precision(),from2->coefficient().precision());
    X z(pr);
    if (k==n) {
        ConstIterator& curr1=from1;
        ConstIterator& curr2=from2;
        MultiIndex a(n);
        X c=curr1->coefficient()*curr2->coefficient();
        MultivariateFourierPolynomial<X> r(n,pr);
        if (not same(c,z)) { r._terms.append(a,c); }
        ++curr1; ++curr2;
        return r;
    } else {
        MultivariateFourierPolynomial<X> zm(n,pr);
        DegreeType i1max=0u; for (auto curr1=from1; curr1!=end1 && equal_up_to(curr1->index(),a1,k); ++curr1) { i1max=curr1->index()[k]; }
        DegreeType i2max=0u; for (auto curr2=from2; curr2!=end2 && equal_up_to(curr2->index(),a2,k); ++curr2) { i2max=curr2->index()[k]; }
        DegreeType imax=i1max+i2max;
        List<MultivariateFourierPolynomial<X>> buf(imax+1u,zm);
        ConstIterator curr1=from1;
        ConstIterator curr2=from2;
        while (from1!=end1 && equal_up_to(from1->index(),a1,k)) {
            curr2=from2;
            while (curr2!=end2 && equal_up_to(curr2->index(),a2,k)) {
                curr1=from1;
                DegreeType a1k=curr1->index()[k]; DegreeType a2k=curr2->index()[k];
                MultivariateFourierPolynomial<X> product = _mul_from(curr1,end1,curr2,end2,k+1u);
                buf[abs(a1k-a2k)]+=product;
                buf[static_cast<DegreeType>(a1k+a2k)]+=product;
            }
            from1=curr1;
        }
        from2=curr2;
        MultivariateFourierPolynomial<X> r(n,pr);
        MultiIndex a(n);
        for(DegreeType i=0; i<=imax; ++i) {
            for(auto term : buf[i]._terms) {
                a=term.index();
                a[k]=i;
                r._terms.append(a,term.coefficient()/2);
            }
        }
        return r;
    }
}



template<class X> auto
MultivariateFourierPolynomial<X>::_write (OutputStream& os) const -> OutputStream& {
    for (auto term_iter = this->_terms.begin(); term_iter!=this->_terms.end(); ++term_iter) {
        auto const& term = *term_iter;
        String cs = to_str(term.coefficient());
        if (cs[0]!='+' && cs[0]!='-' && term_iter!=this->_terms.begin()) { os << "+"; }
        os << cs;
        for (SizeType i=0; i!=this->argument_size(); ++i) { if (term.index()[i]!=0) { os << "*T" << term.index()[i] << "(x" << i << ")"; } }
    }
    return os;
}

*/

} //namespace Ariadne


