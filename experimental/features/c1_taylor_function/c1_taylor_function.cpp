/***************************************************************************
 *            c1_taylor_function.cpp
 *
 *  Copyright 2008--17  Pieter Collins
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

#include "../function/functional.hpp"
#include "../config.hpp"

#include <iostream>
#include <iomanip>

#include "../numeric/float.hpp"
#include "../numeric/rational.hpp"

#include "../utility/macros.hpp"
#include "../utility/exceptions.hpp"
#include "../numeric/numeric.hpp"
#include "../algebra/vector.hpp"
#include "../algebra/multi_index.hpp"
#include "../algebra/expansion.hpp"

#include "c1_taylor_function.hpp"

#define VOLATILE ;

namespace Ariadne {

static const char* plusminus = u8"\u00B1";

template<class X, class Y>
static Y polynomial_evaluate(const std::vector<X>& f, const Y& x)
{
    if(f.size()>=2) {
        Nat i=f.size()-2;
        Y r=x*f[i+1]+f[i];
        while(i!=0) {
            i=i-1;
            r=r*x;
            r+=f[i];
        }
        return r;
    } else if(f.size()==1) {
        return x*X(0)+f[0];
    } else {
        return x*X(0);
    }
}


C1TaylorSeries::C1TaylorSeries()
    : _coefficients(1,FloatDP(0))
    , _zero_error(0)
    , _uniform_error(0)
    , _derivative_error(0)
{
}

C1TaylorSeries::C1TaylorSeries(Nat d)
    : _coefficients(d+1,FloatDP(0))
    , _zero_error(0)
    , _uniform_error(0)
    , _derivative_error(0)
{
}

C1TaylorSeries C1TaylorSeries::constant(FloatDP c) {
    C1TaylorSeries result(1u);
    result._coefficients[0]=FloatDP(c);
    return result;
}

C1TaylorSeries C1TaylorSeries::coordinate() {
    C1TaylorSeries result(1u);
    result._coefficients[1]=1;
    return result;
}

ExactIntervalType C1TaylorSeries::domain() const {
    return ExactIntervalType(-1,+1);
}

Nat C1TaylorSeries::degree() const {
    return this->_coefficients.size()-1;
}

#define ARIADNE_BOUNDS_INTERVAL_SUM
#if defined ARIADNE_MIDPOINT_INTERVAL_SUM
C1TaylorSeries& operator+=(C1TaylorSeries& f, ExactIntervalType ic) {
    FloatDP::set_rounding_upward();
    FloatDP& fv=f._coefficients[0];
    FloatDP c=ic.midpoint();
    FloatDP::set_rounding_upward();
    VOLATILE FloatDP fvu=fv+c;
    VOLATILE FloatDP mfvl=(-fv)-c;
    FloatDP e=(fvu+mfvl)/2;
    e+=max(ic.upper()-c,c-ic.lower());
    f._zero_error+=(fvu+mfvl)/2;
    f._uniform_error+=(fvu+mfvl)/2;
    FloatDP::set_rounding_to_nearest();
    fv+=c;
    ARIADNE_ASSERT_MSG(f._zero_error>=0,"f="<<f<<" c="<<c);
    return f;
}

#elif defined ARIADNE_BOUNDS_INTERVAL_SUM
C1TaylorSeries& operator+=(C1TaylorSeries& f, ValidatedNumericType ic) {
    FloatDP::set_rounding_upward();
    FloatDP& fv=f._coefficients[0];
    FloatDP::set_rounding_upward();
    VOLATILE FloatDP fvu=fv+ic.upper().raw();
    VOLATILE FloatDP mfvl=(-fv)-ic.lower().raw();
    FloatDP e=(fvu+mfvl)/2;
    f._zero_error+=e;
    f._uniform_error+=e;
    FloatDP::set_rounding_to_nearest();
    fv=(fvu-mfvl)/2;
    return f;
}
#endif

C1TaylorSeries& operator*=(C1TaylorSeries& f, ValidatedNumericType ic) {
    FloatDP& fze=f._zero_error;
    FloatDP& fue=f._uniform_error;
    FloatDP& fde=f._derivative_error;
    const FloatDP c=ic.value().raw();
    const FloatDP ac=max(abs(ic.lower().raw()),abs(ic.upper().raw()));

    FloatDP::set_rounding_upward();
    const FloatDP rc=(ic.upper().raw()-ic.lower().raw())/2;

    fze*=ac;
    fue*=ac;
    fde*=ac;
    std::cerr<<"WARNING: operator*=(C1TaylorSeries&,ValidatedNumericType): Mistake in errors\n";
    {
        FloatDP& fv=f._coefficients[0];
        VOLATILE FloatDP fvu=fv*c;
        VOLATILE FloatDP mfvl=(-fv)*c;
        const FloatDP e=(fvu+mfvl)/2;
        fze+=e;
        fue+=e;
    }

    for(Nat i=1; i!=f._coefficients.size(); ++i) {
        FloatDP::set_rounding_upward();
        FloatDP& fv=f._coefficients[i];
        VOLATILE FloatDP fvu=fv*c;
        VOLATILE FloatDP mfvl=(-fv)*c;
        const FloatDP e=(fvu+mfvl)/2;
        fue+=e;
        fde+=i*e;
    };

    FloatDP::set_rounding_to_nearest();
    for(Nat i=0; i!=f._coefficients.size(); ++i) {
        f._coefficients[i]*=c;
    }

    ARIADNE_ASSERT_MSG(f._zero_error>=0,"f="<<f<<" c="<<c);
    return f;

}

C1TaylorSeries operator+(C1TaylorSeries f1, C1TaylorSeries f2) {
    C1TaylorSeries r(std::max(f1.degree(),f2.degree()));

    const std::vector<FloatDP>& f1a=f1._coefficients;
    const std::vector<FloatDP>& f2a=f2._coefficients;
    std::vector<FloatDP>& f0a=r._coefficients;

    FloatDP::set_rounding_upward();
    VOLATILE FloatDP vu=f1a[0]+f2a[0];
    VOLATILE FloatDP mvl=(-f1a[0])-f2a[0];
    r._zero_error=(vu+mvl)/2;
    r._uniform_error=(vu+mvl)/2;
    for(SizeType i=1u; i!=std::min(f1a.size(),f2a.size()); ++i) {
        vu=f1a[i]+f2a[i];
        mvl=(-f1a[i])-f2a[i];
        r._uniform_error+=(vu+mvl)/2;
        r._derivative_error+=i*((vu+mvl)/2);
    }
    r._zero_error+=(f1._zero_error+f2._zero_error);
    r._uniform_error+=f1._uniform_error+f2._uniform_error;
    r._derivative_error+=f1._derivative_error+f2._derivative_error;

    FloatDP::set_rounding_to_nearest();
    for(Nat i=0; i!=std::min(f1a.size(),f2a.size()); ++i) {
        f0a[i]=f1a[i]+f2a[i];
    }
    for(Nat i=std::min(f1a.size(),f2a.size()); i!=f1a.size(); ++i) {
        f0a[i]=f1a[i];
    }
    for(Nat i=std::min(f1a.size(),f2a.size()); i!=f2a.size(); ++i) {
        f0a[i]=f2a[i];
    }
    return r;
}

inline FloatDP abssum(std::vector<FloatDP> const& a) {
    FloatDP s=0;
    for(Nat i=0; i!=a.size(); ++i) {
        s+=abs(a[i]);
    }
    return s;
}

inline FloatDP indabssum(std::vector<FloatDP> const& a) {
    ARIADNE_DEBUG_ASSERT(a.size()>=1);
    FloatDP s=0;
    for(Nat i=1; i!=a.size(); ++i) {
        s+=i*abs(a[i]);
    }
    return s;
}

C1TaylorSeries operator*(C1TaylorSeries f1, C1TaylorSeries f2) {
    C1TaylorSeries fr(f1.degree()+f2.degree());
    // std::cerr<<"d0="<<fr.degree()<<", d1="<<f1.degree()<<", d2="<<f2.degree()<<"\n";

    const std::vector<FloatDP>& f1a=f1._coefficients;
    const std::vector<FloatDP>& f2a=f2._coefficients;
    std::vector<FloatDP>& fra=fr._coefficients;

    const FloatDP& f1ze0=f1._zero_error;
    const FloatDP& f2ze0=f2._zero_error;
    FloatDP& frze0=fr._zero_error;
    const FloatDP& f1e0=f1._uniform_error;
    const FloatDP& f2e0=f2._uniform_error;
    FloatDP& fre0=fr._uniform_error;
    const FloatDP& f1e1=f1._derivative_error;
    const FloatDP& f2e1=f2._derivative_error;
    FloatDP& fre1=fr._derivative_error;

    FloatDP::set_rounding_upward();

    FloatDP f1sa0=abssum(f1a);
    FloatDP f2sa0=abssum(f2a);
    FloatDP f1sa1=indabssum(f1a);
    FloatDP f2sa1=indabssum(f2a);

    VOLATILE FloatDP vu;
    VOLATILE FloatDP mvl;

    vu=f1a[0]*f2a[0];
    mvl=(-f1a[0])*f2a[0];
    frze0=(vu+mvl)/2;
    fre0=(vu+mvl)/2;
    // std::cerr<<"ir="<<0<<", i1="<<0<<", i2="<<0<<"\n";
    for(Nat ir=1; ir!=f1a.size()+f2a.size()-1; ++ir) {
        vu=0.0;
        mvl=0.0;
        for(Nat i1=std::max(0,Int(ir)-Int(f2a.size()-1)); i1!=std::min(Int(ir+1),Int(f1a.size())); ++i1) {
            Nat i2=ir-i1;
            // std::cerr<<"ir="<<ir<<", i1="<<i1<<", i2="<<i2<<"\n";
            ARIADNE_DEBUG_ASSERT(i2<f2a.size());
            vu+=f1a[i1]*f2a[i2];
            mvl+=(-f1a[i1])*f2a[i2];
        }
        fre0+=((vu+mvl)/2);
        fre1+=ir*((vu+mvl)/2);
    }

    frze0+=f1ze0*f2a[0]+f1a[0]*f2ze0+f1ze0*f2ze0;

    fre0+=f1e0*f2sa0+f1sa0*f2e0+f1e0*f2e0;

    fre1 += ( (f1e1*f2sa0 + f1sa1*f2e0 + f1e1*f2e0)
             + (f1e0*f2sa1 + f1sa0*f2e1 + f1e0*f2e1) );

    FloatDP::set_rounding_to_nearest();
    for(Nat i1=0; i1!=f1a.size(); ++i1) {
        for(Nat i2=0; i2!=f2a.size(); ++i2) {
            Nat i0=i1+i2;
            fra[i0]+=f1a[i1]*f2a[i2];
        }
    }

    return fr;
}


C1TaylorSeries compose(C1TaylorSeries f, C1TaylorSeries g) {
    Nat i=f.degree();
    C1TaylorSeries r=C1TaylorSeries::constant(FloatDP(f._coefficients[i]));
    while (i!=0) {
        i=i-i;
        r=r*g;
        r+=ValidatedNumericType(f._coefficients[i]);
    }

    FloatDP::set_rounding_upward();
    r._zero_error+=f._zero_error;
    r._uniform_error+=f._uniform_error;
    r._derivative_error+=f._derivative_error;
    FloatDP::set_rounding_to_nearest();

    return r;
}

UpperIntervalType evaluate(C1TaylorSeries f, UpperIntervalType x) {
    Nat i=f.degree();
    UpperIntervalType r=UpperIntervalType(FloatDP(f._coefficients[i]));
    while (i!=0) {
        i=i-i;
        r*=x;
        r+=UpperIntervalType(FloatDP(f._coefficients[i]));
    }
    if(f._zero_error+FloatDP(x.centre())*f._derivative_error < f._uniform_error) {
        r+=UpperIntervalType(-f._zero_error,+f._zero_error);
        r+=UpperIntervalType(x)*UpperIntervalType(-f._derivative_error,+f._derivative_error);
    } else {
        r+=UpperIntervalType(-f._uniform_error,+f._uniform_error);
    }
    return r;
}


OutputStream& operator<<(OutputStream& os, const C1TaylorSeries& f) {
    os << "(" << f._coefficients[0] << "±" << f._zero_error << "/" << f._uniform_error << ")+(" << f._coefficients[1] << plusminus << f._derivative_error <<")*x";
    for(Nat i=2; i<f._coefficients.size(); ++i) {
        if (f._coefficients[i]>=0) { os << "+"; }
        os << f._coefficients[i] << "*x^" << i;
    }
    return os;
}



C1TaylorFunction::C1TaylorFunction()
    : C1TaylorFunction(0u) { }

C1TaylorFunction::C1TaylorFunction(SizeType as)
    : _expansion(as)
    , _zero_error(0)
    , _uniform_error(0)
    , _derivative_errors(as)
{
    MultiIndex ind(as);
    for(Nat i=0; i!=as; ++i) {
        ind[i]=1;
        _expansion.append(ind,0);
        ind[i]=0;
    }
    _expansion.append(ind,0);
    _expansion.reverse_lexicographic_sort();
}

C1TaylorFunction C1TaylorFunction::constant(SizeType as, FloatDP c) {
    C1TaylorFunction result(as);
    MultiIndex ind=MultiIndex::zero(as);
    //result._expansion[ind]=FloatDP(c);
    result._expansion.append(ind,c);
    return result;
}

C1TaylorFunction C1TaylorFunction::coordinate(SizeType as, SizeType j) {
    C1TaylorFunction result(as);
    MultiIndex ind=MultiIndex::unit(as,j);
    //result._expansion[ind]=1;
    result._expansion.append(ind,1.0);
    return result;
}

Nat C1TaylorFunction::argument_size() const {
    return this->_expansion.argument_size();
}

Void C1TaylorFunction::clear() {
    this->_expansion.clear();
    this->_zero_error=0;
    this->_uniform_error=0;
    for(Nat j=0; j!=this->_derivative_errors.size(); ++j) {
        this->_derivative_errors[j]=0;
    }
}

C1TaylorFunction& C1TaylorFunction::operator=(NumericType ic) {
    this->clear();
    FloatDP::set_rounding_upward();
    FloatDP e=(ic.upper().raw()-ic.lower().raw())/2;
    this->_zero_error=e;
    this->_uniform_error=e;
    FloatDP::set_rounding_to_nearest();
    FloatDP c=(ic.upper().raw()-ic.lower().raw())/2;
    this->_expansion.append(MultiIndex(this->argument_size()),c);
    return *this;
}

C1TaylorFunction& operator+=(C1TaylorFunction& f, FloatDP ec) {
    const FloatDP& c=ec;
    if(f._expansion.empty() || (--f._expansion.end())->index().degree()!=0) {
        f._expansion.append(MultiIndex(f.argument_size()),c);
        return f;
    }
    //ARIADNE_DEBUG_ASSERT(f._expansion.back().key().degree()==0);
    FloatDP::set_rounding_upward();
    //FloatDP& fv=f._expansion.back().data();
    ReferenceType<FloatDP> fv=(--f._expansion.end())->coefficient();
    FloatDP& fze=f._zero_error;
    FloatDP& fe=f._uniform_error;
    FloatDP::set_rounding_upward();
    VOLATILE FloatDP fvu=fv+c;
    VOLATILE FloatDP mfvl=(-fv)-c;
    fze+=(fvu+mfvl)/2;
    fe+=(fvu+mfvl)/2;
    FloatDP::set_rounding_to_nearest();
    fv+=c;
    ARIADNE_ASSERT_MSG(f._zero_error>=0,"f="<<f<<" c="<<c);
    return f;
}

C1TaylorFunction& operator*=(C1TaylorFunction& f, FloatDP ec) {
    FloatDP::set_rounding_upward();
    FloatDP& fze=f._zero_error;
    FloatDP& fue=f._uniform_error;
    Array<FloatDP>& fde=f._derivative_errors;
    const FloatDP& c=ec;
    const FloatDP ac=abs(c);

    FloatDP::set_rounding_upward();
    fze*=ac;
    fue*=ac;
    for(Nat j=0; j!=f.argument_size(); ++j) {
        fde[j]*=ac;
    }

    for(Expansion<MultiIndex,FloatDP>::Iterator iter=f._expansion.begin();
        iter!=f._expansion.end(); ++iter)
    {
        ConstReferenceType<MultiIndex> a=iter->index();
        ReferenceType<FloatDP> fv=iter->coefficient();
        VOLATILE FloatDP fvu=fv*c;
        VOLATILE FloatDP mfvl=(-fv)*c;
        const FloatDP e=(fvu+mfvl)/2;
        if(a.degree()==0) { fze+=e; }
        fue+=e;
        for(Nat j=0; j!=f.argument_size(); ++j) {
            fde[j]+=a[j]*e;
        }
    }

    FloatDP::set_rounding_to_nearest();
    for(Expansion<MultiIndex,FloatDP>::Iterator iter=f._expansion.begin();
        iter!=f._expansion.end(); ++iter)
    {
        ReferenceType<FloatDP> fv=iter->coefficient();
        fv*=c;
    }

    return f;

}



C1TaylorFunction operator+(C1TaylorFunction f1, C1TaylorFunction f2) {
    ARIADNE_PRECONDITION(f1.argument_size()==f2.argument_size());
    const Nat n=f1.argument_size();
    C1TaylorFunction f0(n);
    f0._expansion.clear();
    f0._expansion.reserve(f1._expansion.number_of_nonzeros()+f2._expansion.number_of_nonzeros());

    Expansion<MultiIndex,FloatDP>::ConstIterator i1=f1._expansion.begin();
    Expansion<MultiIndex,FloatDP>::ConstIterator i2=f2._expansion.begin();
    while(i1!=f1._expansion.end() && i2!=f2._expansion.end()) {
        if(i1->index()==i2->index()) {
            ConstReferenceType<MultiIndex> a = i1->index();
            FloatDP::set_rounding_upward();
            VOLATILE FloatDP fvu=i1->coefficient()+i2->coefficient();
            VOLATILE FloatDP mfvl=(-i1->coefficient())-i2->coefficient();
            const FloatDP e=(fvu+mfvl)/2;
            if(a.degree()==0) {
                f0._zero_error+=e;
            }
            f0._uniform_error+=e;
            for(Nat j=0; j!=n; ++j) {
                f0._derivative_errors[j]+=a[j]*e;
            }
            FloatDP::set_rounding_to_nearest();
            f0._expansion.append(a,i1->coefficient()+i2->coefficient());
            ++i1;
            ++i2;
        } else if(reverse_lexicographic_less(i1->index(),i2->index())) {
            f0._expansion.append(i1->index(),i1->coefficient());
            ++i1;
        } else {
            f0._expansion.append(i2->index(),i2->coefficient());
            ++i2;
        }
    }
    while(i1!=f1._expansion.end()) {
        f0._expansion.append(i1->index(),i1->coefficient());
        ++i1;
    }
    while(i2!=f2._expansion.end()) {
        f0._expansion.append(i2->index(),i2->coefficient());
        ++i2;
    }

    FloatDP::set_rounding_upward();
    f0._zero_error+=(f1._zero_error+f2._zero_error);
    f0._uniform_error+=(f1._uniform_error+f2._uniform_error);
    for(Nat j=0; j!=n; ++j) {
        f0._derivative_errors[j]+=(f1._derivative_errors[j]+f2._derivative_errors[j]);
    }
    FloatDP::set_rounding_to_nearest();

    return f0;
}


Void fma(C1TaylorFunction& f0, const C1TaylorFunction& f1, const C1TaylorFunction& f2, const FloatDP& c3, const MultiIndex& a3) {
    ARIADNE_PRECONDITION(f0.argument_size()==f1.argument_size());
    ARIADNE_PRECONDITION(f1.argument_size()==f2.argument_size());
    const Nat n=f1.argument_size();
    MultiIndex a(n);
    f0.clear();

    Expansion<MultiIndex,FloatDP>::ConstIterator i1=f1._expansion.begin();
    Expansion<MultiIndex,FloatDP>::ConstIterator i2=f2._expansion.begin();
    while(i1!=f1._expansion.end() && i2!=f2._expansion.end()) {
        if(i1->index()==i2->index()+a3) {
            ConstReferenceType<MultiIndex> a = i1->index();
            FloatDP::set_rounding_upward();
            VOLATILE FloatDP fvu=i1->coefficient()+i2->coefficient()*c3;
            VOLATILE FloatDP mfvl=(-i1->coefficient())+i2->coefficient()*(-c3);
            const FloatDP e=(fvu+mfvl)/2;
            if(a.degree()==0) {
                f0._zero_error+=e;
            }
            f0._uniform_error+=e;
            for(Nat j=0; j!=n; ++j) {
                f0._derivative_errors[j]+=a[j]*e;
            }
            FloatDP::set_rounding_to_nearest();
            f0._expansion.append(a,i1->coefficient()+i2->coefficient()*c3);
            ++i1;
            ++i2;
        } else if(reverse_lexicographic_less(i1->index(),i2->index()+a3)) {
            f0._expansion.append(i1->index(),i1->coefficient());
            ++i1;
        } else {
            a=i2->index()+a3;
            FloatDP::set_rounding_upward();
            VOLATILE FloatDP fvu=i2->coefficient()*c3;
            VOLATILE FloatDP mfvl=i2->coefficient()*(-c3);
            const FloatDP e=(fvu+mfvl)/2;
            if(a.degree()==0) {
                f0._zero_error+=e;
            }
            f0._uniform_error+=e;
            for(Nat j=0; j!=n; ++j) {
                f0._derivative_errors[j]+=(uchar)a[j]*e;
            }
            FloatDP::set_rounding_to_nearest();
            f0._expansion.append(a,i2->coefficient()*c3);
            ++i2;
        }
    }
    while(i1!=f1._expansion.end()) {
        f0._expansion.append(i1->index(),i1->coefficient());
        ++i1;
    }
    while(i2!=f2._expansion.end()) {
        FloatDP::set_rounding_upward();
        VOLATILE FloatDP fvu=i2->coefficient()*c3;
        VOLATILE FloatDP mfvl=i2->coefficient()*(-c3);
        const FloatDP e=(fvu+mfvl)/2;
        if(a.degree()==0) {
            f0._zero_error+=e;
        }
        f0._uniform_error+=e;
        for(Nat j=0; j!=n; ++j) {
            f0._derivative_errors[j]+=(uchar)a[j]*e;
        }
        FloatDP::set_rounding_to_nearest();
        f0._expansion.append(i2->index()+a3,i2->coefficient()*c3);
        ++i2;
    }

    FloatDP::set_rounding_upward();
    f0._zero_error+=(f1._zero_error+f2._zero_error*abs(c3));
    f0._uniform_error+=(f1._uniform_error+f2._uniform_error*abs(c3));
    for(Nat j=0; j!=n; ++j) {
        f0._derivative_errors[j]+=(f1._derivative_errors[j]+f2._derivative_errors[j]*abs(c3));
    }
    FloatDP::set_rounding_to_nearest();

}

C1TaylorFunction operator*(C1TaylorFunction f1, C1TaylorFunction f2) {
    C1TaylorFunction f0a(f2);
    f0a._expansion.clear();
    C1TaylorFunction f0b(f0a);
    f0b.clear();
    C1TaylorFunction* ftp=&f0a;
    C1TaylorFunction* frp=&f0b;
    for(Expansion<MultiIndex,FloatDP>::ConstIterator i2=f2._expansion.begin();
        i2!=f2._expansion.end(); ++i2)
    {
        fma(*frp,*ftp,f1,i2->coefficient(),i2->index());
        std::swap(ftp,frp);
    }
    return *frp;
}

UpperIntervalType evaluate(C1TaylorFunction f, Vector<UpperIntervalType> x) {
    UpperIntervalType r=horner_evaluate(reinterpret_cast<Expansion<MultiIndex,FloatDPValue>const&>(f._expansion),x);
    r += UpperIntervalType(-f._uniform_error,+f._uniform_error);
    return r;
}


C1TaylorFunction compose(C1TaylorSeries f, C1TaylorFunction g) {
    C1TaylorFunction r = g;
    r.clear();

    Nat i=f.degree();
    r+=FloatDP(f._coefficients[i]);

    while(i!=0) {
        r=r*g;
        --i;
        r+=FloatDP(f._coefficients[i]);
    }
    std::cerr<<"intermediate="<<r<<"\n";

    // TODO: How do first derivatives change?
    ARIADNE_NOT_IMPLEMENTED;
    return r;
}

C1TaylorFunction compose(C1TaylorFunction f, Vector<C1TaylorFunction> g) {
    C1TaylorFunction r=horner_evaluate(reinterpret_cast<Expansion<MultiIndex,FloatDPValue>const&>(f._expansion),g);
    std::cerr<<"intermediate="<<r<<"\n";
    r._uniform_error += f._uniform_error;
    r._zero_error += f._zero_error;
    // TODO: How do first derivatives change?
    ARIADNE_NOT_IMPLEMENTED;
    return r;
}

template<class T> struct ListForm {
    ListForm(const T& t) : value(t) { } const T& value;
};
template<class T> ListForm<T> list_form(const T& t) { return ListForm<T>(t); }

OutputStream& operator<<(OutputStream& os, const ListForm<Expansion<MultiIndex,FloatDP>>& lfe) {
    const Expansion<MultiIndex,FloatDP>& e=lfe.value;
    os << "{ ";
    for(Expansion<MultiIndex,FloatDP>::ConstIterator iter=e.begin();
        iter!=e.end(); ++iter)
    {
        if(iter!=e.begin()) { os << ", "; }
        for(Nat i=0; i!=iter->index().size(); ++i) {
            os << Nat(iter->index()[i]);
            if(i+1!=iter->index().size()) { os << ","; }
        }
        os << ":" << iter->coefficient();
    }
    return os << " }";
}

OutputStream& operator<<(OutputStream& os, const C1TaylorFunction& f) {
    os << "C1TaylorFunction( _expansion=" << list_form(f._expansion)
       << ", _zero_error=" << f._zero_error
       << ", _uniform_error=" << f._uniform_error
       << ", _derivative_errors=" << f._derivative_errors
       << ")";
    return os;

    for(Expansion<MultiIndex,FloatDP>::ConstIterator iter=f._expansion.begin();
        iter!=f._expansion.end(); ++iter)
    {
        os << *iter;
    }
}



} // namespace Ariadne
