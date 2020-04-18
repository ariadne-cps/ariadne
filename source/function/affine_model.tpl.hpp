/***************************************************************************
 *            affine_model.tcc
 *
 *  Copyright  2009-20  Pieter Collins
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

#include "../numeric/numeric.hpp"

#include "../algebra/vector.hpp"
#include "../function/function.hpp"
#include "../function/taylor_model.hpp"
#include "../function/affine_model.hpp"

#include "../function/affine.hpp"
#include "../function/taylor_function.hpp"
#include "../algebra/vector.hpp"

namespace Ariadne {

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

template<class X> decltype(auto) values(Vector<X> const& v) {
    //typedef typename std::remove_reference<decltype(a.zero_element().value())>::type C;
    return elementwise([&](X const& x){return x.value();},v);
}

/*
FloatDPBounds const& make_singleton(IntervalRangeType const& ivl, DoublePrecision) {
    return reinterpret_cast<FloatDPBounds const&>(ivl);
}

FloatMPBounds const& make_singleton(IntervalRangeType const& ivl, MultiplePrecision pr) {
    return FloatMPBounds(ivl.lower(),ivl.upper());
}
*/

/*
Bounds<FloatMP> cast_singleton(Interval<UpperBound<FloatMP>> const& ivl) {
    return Bounds<FloatMP>(ivl.lower(),ivl.upper());
}

Vector<Bounds<FloatMP>> cast_singleton(Vector<Interval<UpperBound<FloatMP>>> const& bx) {
    return elementwise(bx,[&](Interval<UpperBound<FloatMP>> const& ivl){return cast_singleton(ivl);});
}
*/

template<class F> struct AlgebraOperations<AffineModel<ApproximateTag,F>> {
    typedef typename F::PrecisionType PrecisionType;
    typedef AffineModel<ApproximateTag,F> AffineModelType;
    typedef typename AffineModelType::NumericType NumericType;
    typedef typename AffineModelType::CoefficientType CoefficientType;

    static AffineModelType apply(Neg, const AffineModelType& a) {
        SizeType n=a.argument_size();
        PrecisionType prec=a.precision();

        AffineModelType r(n,prec);
        r=CoefficientType(-a.value().raw());
        for(SizeType i=0; i!=n; ++i) {
            r[i]=CoefficientType(-a.gradient(i).raw());
        }
        return r;
    }

    static AffineModelType apply(Add, const AffineModelType& a1, const AffineModelType& a2) {
        ARIADNE_ASSERT_MSG(a1.argument_size()==a2.argument_size(),"a1="<<a1<<" a2="<<a2);
        SizeType n=a1.argument_size();
        PrecisionType prec=max(a1.precision(),a2.precision());

        AffineModelType r(n,prec);
        r=CoefficientType(a1.value().raw()+a2.value().raw());
        for(SizeType i=0; i!=n; ++i) {
            r[i]=CoefficientType(a1.gradient(i).raw()+a2.gradient(i).raw());
        }
        return r;
    }

    static AffineModelType apply(Mul, const AffineModelType& a1, const AffineModelType& a2) {
        ARIADNE_ASSERT_MSG(a1.argument_size()==a2.argument_size(),"a1="<<a1<<" a2="<<a2);
        SizeType n=a1.argument_size();
        PrecisionType prec=max(a1.precision(),a2.precision());

        AffineModelType r(n,prec);
        r=CoefficientType(a1.value().raw()+a2.value().raw());
        for(SizeType i=0; i!=n; ++i) {
            r[i]=CoefficientType(a1.value().raw()*a2.gradient(i).raw()+a1.gradient(i).raw()+a2.value().raw());
        }
        return r;
    }

    static AffineModelType apply(Add, const NumericType& c, const AffineModelType& a) {
        AffineModelType r=a;
        r._c=CoefficientType(c.raw() + a.value().raw());
        return r;
    }

    static AffineModelType apply(Add, const AffineModelType& a, const NumericType& c) {
        return apply(Add(), c,a);
    }

    static AffineModelType apply(Mul, const NumericType& c, const AffineModelType& a) {
        SizeType n=a.argument_size();
        PrecisionType prec=a.precision();

        AffineModelType r(n,prec);
        r=CoefficientType(a.value().raw()*c.raw());
        for(SizeType i=0; i!=n; ++i) {
            r[i]=CoefficientType(a.gradient(i).raw()*c.raw());
        }
        return r;
    }

    static AffineModelType apply(Mul, const AffineModelType& a, const NumericType& c) {
        return apply(Mul(), c,a);
    }
};


template<class F> struct AlgebraOperations<AffineModel<ValidatedTag,F>> {
    typedef typename F::PrecisionType PrecisionType;
    typedef AffineModel<ValidatedTag,F> AffineModelType;
    typedef typename AffineModelType::NumericType NumericType;
    typedef typename AffineModelType::CoefficientType CoefficientType;
    typedef typename AffineModelType::ErrorType ErrorType;

    static AffineModelType apply(Neg, const AffineModelType& a) {
        SizeType n=a.argument_size();
        AffineModelType r(n,a.precision());
        r=CoefficientType( -a.value().raw() );
        for(SizeType i=0; i!=n; ++i) {
            r[i]=CoefficientType(-a.gradient(i).raw());
        }
        r.set_error( ErrorType( a.error().raw() ) );

        return r;
    }

    static AffineModelType apply(Add, const AffineModelType& a1, const AffineModelType& a2) {
        ARIADNE_ASSERT_MSG(a1.argument_size()==a2.argument_size(),"a1="<<a1<<" a2="<<a2);
        SizeType n=a1.argument_size();
        PrecisionType prec=max(a1.precision(),a2.precision());

        AffineModelType r(n,prec);
        r=CoefficientType( a1.value().raw()+a2.value().raw() );
        for(SizeType i=0; i!=n; ++i) {
            r[i]=CoefficientType(a1.gradient(i).raw()+a2.gradient(i).raw());
        }

        F::set_rounding_upward();

        F te=0.0;
        for(SizeType j=0; j!=n; ++j) {
            F mrjl = (-a1.gradient(j).raw())-a2.gradient(j).raw();
            F  rju = ( a1.gradient(j).raw())+a2.gradient(j).raw();
            te+=(rju+mrjl);
        }
        F mrl = (-a1.value().raw())-a2.value().raw();
        F  ru = ( a1.value().raw())+a2.value().raw();
        te += (ru+mrl);

        r.set_error( ErrorType(te/2 + (a1.error().raw()+a2.error().raw()) ) );

        F::set_rounding_to_nearest();

        return r;
    }

    static AffineModelType apply(Mul, const AffineModelType& a1, const AffineModelType& a2) {
        ARIADNE_ASSERT_MSG(a1.argument_size()==a2.argument_size(),"a1="<<a1<<" a2="<<a2);
        SizeType n=a1.argument_size();
        PrecisionType prec=max(a1.precision(),a2.precision());

        AffineModelType r(n,prec);
        r=CoefficientType( a1.value().raw()+a2.value().raw() );
        for(SizeType i=0; i!=n; ++i) {
            r[i]=CoefficientType(a1.gradient(i).raw()*a2.value().raw()+a1.value().raw()*a2.gradient(i).raw());
        }

        F::set_rounding_upward();

        F mrvl = ((-a1.value().raw())*a2.value().raw());
        F  rvu = ((+a1.value().raw())*a2.value().raw());
        F te=(mrvl+rvu);
        for(SizeType j=0; j!=n; ++j) {
            F mrjl = ((-a1.value().raw())*a1.gradient(j).raw()+a2.gradient(j).raw()*(-a2.value().raw()));
            F  rju = ((+a1.value().raw())*a1.gradient(j).raw()+a2.gradient(j).raw()*(+a2.value().raw()));
            te+=(rju+mrjl);
        }
        F roe=hlf(te); // roundoff error

        F tre(prec); // truncation error
        for(SizeType j1=0; j1!=n; ++j1) {
            for(SizeType j2=0; j2!=n; ++j2) {
                tre = tre + abs(a1.gradient(j1).raw()) * abs(a2.gradient(j2).raw());
            }
        }

        F nrm1=abs(a1.value().raw());
        F nrm2=abs(a2.value().raw());
        for(SizeType j=0; j!=n; ++j) {
            nrm1=nrm1+abs(a1.value().raw());
            nrm2=nrm2+abs(a2.value().raw());
        }
        F ace = a1.error().raw()*nrm2 + nrm1 * a1.error().raw();  // accumulated error
        r.set_error( ErrorType(tre + roe + ace) );

        F::set_rounding_to_nearest();

        return r;
    }

    static AffineModelType apply(Add, const NumericType& c, const AffineModelType& a) {
        AffineModelType r=a;
        F cm=c.value().raw();
        r.set_value( CoefficientType( cm + a.value().raw() ) );

        F::set_rounding_upward();

        F mrl = (-a.value().raw())-cm;
        F  ru = ( a.value().raw())+cm;
        F te = (ru+mrl)/2;

        r.set_error( ErrorType( a.error().raw() + max(c.upper().raw()-cm,cm-c.lower().raw()) + te) );

        F::set_rounding_to_nearest();

        return r;
    }

    static AffineModelType apply(Add, const AffineModelType& a, const NumericType& c) {
        return apply(Add(), c,a);
    }

    static AffineModelType apply(Mul, const NumericType& c, const AffineModelType& a) {
        SizeType n=a.argument_size();
        PrecisionType prec=a.precision();

        AffineModelType r(n,prec);
        F cm=c.value().raw();
        r=CoefficientType(a.value().raw()*cm);
        for(SizeType i=0; i!=n; ++i) {
            r[i].raw()=a.gradient(i).raw()*cm;
        }

        F::set_rounding_upward();

        F te=0.0;
        for(SizeType j=0; j!=n; ++j) {
            F mca=(-cm)*a.gradient(j).raw();
            F ca= cm*a.gradient(j).raw();
            te+=(ca+mca);
        }
        F mca=(-cm)*a.value().raw();
        F ca= cm*a.value().raw();

        F re=0.0;
        if(c.lower_raw()!=c.upper_raw()) {
            F ce=max(c.upper().raw()-cm,cm-c.lower().raw());
            for(SizeType j=0; j!=n; ++j) {
                re+=abs(a.gradient(j).raw()*ce);
            }
        }

        r.set_error(ErrorType(abs(cm)*a.error().raw() + ((ca+mca) + te)/2 + re));

        F::set_rounding_to_nearest();

        return r;
    }

    static AffineModelType apply(Mul, const AffineModelType& a, const NumericType& c) {
        return apply(Mul(), c,a);
    }
};


template<class F> AffineModel<ValidatedTag,F>::AffineModel(const Affine<NumericType>& affine)
    : AffineModel(affine.argument_size(),affine.value().precision())
{
    AffineModel<ValidatedTag,F>& affine_model=*this;
    affine_model = affine.value().value();
    for(SizeType j=0; j!=affine.argument_size(); ++j) {
        affine_model[j] = affine[j].value();
    }
    F::set_rounding_upward();
    F e = 0.0;
    for(SizeType j=0; j!=affine.argument_size(); ++j) {
        e += max(affine.gradient(j).upper().raw()-affine_model.gradient(j).raw(),affine_model.gradient(j).raw()-affine.gradient(j).lower().raw());
    }
    e += max(affine.value().upper().raw()-affine_model.value().raw(),affine_model.value().raw()-affine.value().lower().raw());
    affine_model.set_error(ErrorType(e));
    F::set_rounding_to_nearest();
}

template<class F> AffineModel<ValidatedTag,F>::AffineModel(const TaylorModel<ValidatedTag,F>& taylor_model)
    : AffineModel(taylor_model.argument_size(),taylor_model.value().precision())
{
    AffineModel<ValidatedTag,F>& affine_model=*this;

    typename F::RoundingModeType rnd=F::get_rounding_mode();
    F::set_rounding_upward();
    for(auto iter=taylor_model.begin(); iter!=taylor_model.end(); ++iter) {
        if(iter->index().degree()>=2) {
            affine_model.set_error(mag(iter->coefficient())+affine_model.error());
        } else if(iter->index().degree()==1) {
            for(SizeType i=0; i!=taylor_model.argument_size(); ++i) {
                if(iter->index()[i]==1) {
                    affine_model.set_gradient(i,iter->coefficient());
                    break;
                }
            }
        } else {
            affine_model.set_value(iter->coefficient());
        }
    }
    affine_model.set_error(taylor_model.error()+affine_model.error());
    F::set_rounding_mode(rnd);
}

template<class P, class PR> AffineModel<P,RawFloatType<PR>> affine_model(const BoxDomainType& domain, const ScalarMultivariateFunction<P>& function, PR precision)
{
    ValidatedScalarMultivariateTaylorFunctionModelDP tf(domain,function,AffineSweeper<RawFloat<PR>>(precision));
    return affine_model(tf.model());
}

inline decltype(auto) operator<(FloatMPApproximation x,const FloatDPValue y) { return x<FloatMPApproximation(y.get_d(),x.precision()); }
inline decltype(auto) operator<(FloatMPValue x,const FloatDPValue y) { return x<FloatMPValue(Dyadic(y.get_d()),x.precision()); }

// FIXME: Clean up construction of interval with precision PrecisionType<F>
template<class F> AffineModel<ApproximateTag,F> AffineModel<ApproximateTag,F>::scaling(SizeType n, SizeType j, const IntervalDomainType& codom, PrecisionType pr)
{
    AffineModel<ApproximateTag,F> r(n,pr);
    FloatApproximation<PR> l(codom.lower().get_d(),pr);
    FloatApproximation<PR> u(codom.upper().get_d(),pr);
    r.set_gradient(j,hlf(u-l));
    r.set_value(hlf(l+u));
    return r;
}

template<class F> Vector<AffineModel<ApproximateTag,F>> AffineModel<ApproximateTag,F>::scalings(const BoxDomainType& codom, PrecisionType pr)
{
    SizeType n=codom.dimension();
    Vector<AffineModel<ApproximateTag,F>> r(n,AffineModel<ApproximateTag,F>(n,pr));
    for(SizeType i=0; i!=n; ++i) {
        r[i] = AffineModel<ApproximateTag,F>::scaling(n,i,codom[i],pr);
    }
    return r;
}

template<class F> AffineModel<ValidatedTag,F> AffineModel<ValidatedTag,F>::scaling(SizeType n, SizeType j, const IntervalDomainType& codom, PrecisionType pr)
{
    AffineModel<ValidatedTag,F> r(n,pr);
    FloatValue<PR> l(Dyadic(codom.lower()),pr);
    FloatValue<PR> u(Dyadic(codom.upper()),pr);
    Interval<FloatValue<PR>> ivl(l,u);
    r.set_gradient(j,1);
    r*=ivl.radius();
    r+=ivl.midpoint();
    return r;
}

template<class F> Vector<AffineModel<ValidatedTag,F>> AffineModel<ValidatedTag,F>::scalings(const BoxDomainType& codom, PrecisionType pr)
{
    SizeType n=codom.size();
    Vector<AffineModel<ValidatedTag,F>> r(n,AffineModel<ValidatedTag,F>(n,pr));
    for(SizeType i=0; i!=n; ++i) {
        r[i] = AffineModel<ValidatedTag,F>::scaling(n,i,codom[i],pr);
    }
    return r;
}

template<class F> auto AffineModel<ApproximateTag,F>::range() const -> RangeType
{
    auto v=this->value();
    auto e=this->gradient().zero_element();
    for(SizeType i=0; i!=this->argument_size(); ++i) {
        e+=abs(this->gradient(i));
    }
    return RangeType(v-e,v+e);
}

template<class F> auto AffineModel<ValidatedTag,F>::range() const -> RangeType
{
    auto v=this->value();
    ErrorType e=this->error();
    for(SizeType i=0; i!=this->argument_size(); ++i) {
        e+=abs(this->gradient(i));
    }
    return RangeType(v-e,v+e);
}

template<class F> OutputStream& AffineModel<ApproximateTag,F>::_write(OutputStream& os) const
{
    os << this->value().raw();
    for(SizeType j=0; j!=this->argument_size(); ++j) {
        if(decide(this->gradient(j)>0)) { os << "+"; }
        if(definitely(this->gradient(j)!=0)) { os << this->gradient(j) << "*x" << j; }
    }
    return os;
}

template<class F> OutputStream& AffineModel<ValidatedTag,F>::_write(OutputStream& os) const
{
    os << this->value().raw();
    for(SizeType j=0; j!=this->argument_size(); ++j) {
        if(decide(this->gradient(j)>0)) { os << "+"; }
        if(definitely(this->gradient(j)!=0)) { os << this->gradient(j) << "*x" << j; }
    }
    return os << "+/-" << this->error();

}

template<class F> auto
AffineModel<ApproximateTag,F>::_compose(const ScalarMultivariateFunction<P>& f, Vector<AffineModel<P,F>> const& g) -> AffineModel<P,F> {
    auto c=values(g);
    auto b=f(c);
    auto A=f.gradient(c);
    return b+A*(g-c);
}

template<class F> auto
AffineModel<ApproximateTag,F>::_compose(const VectorMultivariateFunction<P>& f, Vector<AffineModel<P,F>> const& g) -> Vector<AffineModel<P,F>> {
    auto c=values(g);
    auto b=f(c);
    auto A=f.jacobian(c);
    return b+A*(g-c);
}

template<class F> auto
AffineModel<ValidatedTag,F>::_compose(const ScalarMultivariateFunction<P>& f, Vector<AffineModel<P,F>> const& g) -> AffineModel<P,F> {
    auto d = ranges(g);
    auto r = cast_singleton(d);
    auto c=values(g);
    auto b=f(c);
    auto A=f.gradient(r);
    return b+A*(g-c);
}

template<class F> auto
AffineModel<ValidatedTag,F>::_compose(const VectorMultivariateFunction<P>& f, Vector<AffineModel<P,F>> const& g) -> Vector<AffineModel<P,F>> {
    auto d = ranges(g);
    auto r = cast_singleton(d);
    auto c=values(g);
    auto b=f(c);
    auto A=f.jacobian(r);
    return b+A*(g-c);
}







} //namespace Ariadne


