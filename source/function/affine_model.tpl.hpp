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

#include "numeric/numeric.hpp"

#include "algebra/vector.hpp"
#include "function/function.hpp"
#include "function/taylor_model.hpp"
#include "function/affine_model.hpp"

#include "model_utilities.hpp"

#include "function/affine.hpp"
#include "function/taylor_function.hpp"
#include "algebra/vector.hpp"

namespace Ariadne {

template<class X> decltype(auto) values(Vector<X> const& v) {
    //typedef typename std::remove_reference<decltype(a.zero_element().value())>::type C;
    return elementwise([&](X const& x){return x.value();},v);
}

/*
FloatDPBounds const& make_singleton(IntervalRangeType const& ivl, DoublePrecision) {
    return reinterpret_cast<FloatDPBounds const&>(ivl);
}

FloatMPBounds const& make_singleton(IntervalRangeType const& ivl, MultiplePrecision pr) {
    return FloatMPBounds(ivl.lower_bound(),ivl.upper_bound());
}
*/

/*
Bounds<FloatMP> cast_singleton(Interval<UpperBound<FloatMP>> const& ivl) {
    return Bounds<FloatMP>(ivl.lower_bound(),ivl.upper_bound());
}

Vector<Bounds<FloatMP>> cast_singleton(Vector<Interval<UpperBound<FloatMP>>> const& bx) {
    return elementwise(bx,[&](Interval<UpperBound<FloatMP>> const& ivl){return cast_singleton(ivl);});
}
*/

template<class FLT> struct AlgebraOperations<AffineModel<ApproximateTag,FLT>> {
    typedef typename FLT::PrecisionType PrecisionType;
    typedef AffineModel<ApproximateTag,FLT> AffineModelType;
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
        r=a1.value()+a2.value();
        for(SizeType i=0; i!=n; ++i) {
            r[i]=a1.gradient(i)+a2.gradient(i);
        }
        return r;
    }

    static AffineModelType apply(Mul, const AffineModelType& a1, const AffineModelType& a2) {
        ARIADNE_ASSERT_MSG(a1.argument_size()==a2.argument_size(),"a1="<<a1<<" a2="<<a2);
        SizeType n=a1.argument_size();
        PrecisionType prec=max(a1.precision(),a2.precision());

        AffineModelType r(n,prec);
        r=a1.value()+a2.value();
        for(SizeType i=0; i!=n; ++i) {
            r[i]=a1.value()*a2.gradient(i)+a1.gradient(i)+a2.value();
        }
        return r;
    }

    static AffineModelType apply(Add, const NumericType& c, const AffineModelType& a) {
        AffineModelType r=a;
        r._c=c+a.value();
        return r;
    }

    static AffineModelType apply(Add, const AffineModelType& a, const NumericType& c) {
        return apply(Add(), c,a);
    }

    static AffineModelType apply(Mul, const NumericType& c, const AffineModelType& a) {
        SizeType n=a.argument_size();
        PrecisionType prec=a.precision();

        AffineModelType r(n,prec);
        r=a.value()*c;
        for(SizeType i=0; i!=n; ++i) {
            r[i]=a.gradient(i)*c;
        }
        return r;
    }

    static AffineModelType apply(Mul, const AffineModelType& a, const NumericType& c) {
        return apply(Mul(), c,a);
    }
};



template<class FLT> PositiveUpperBound<FLT> norm(Covector<FLT> const& g) {
    PositiveUpperBound<FLT> r=mag(g.zero_element());
    for (SizeType i=0; i!=g.size(); ++i) {
        r += mag(g[i]);
    }
    return r;
}

template<class FLT> struct AlgebraOperations<AffineModel<ValidatedTag,FLT>> {
    typedef typename FLT::PrecisionType PrecisionType;
    typedef AffineModel<ValidatedTag,FLT> AffineModelType;
    typedef typename AffineModelType::NumericType NumericType;
    typedef typename AffineModelType::CoefficientType CoefficientType;
    typedef typename AffineModelType::ErrorValueType ErrorValueType;

    static AffineModelType apply(Neg, const AffineModelType& a) {
        SizeType n=a.argument_size();
        AffineModelType r(n,a.precision());
        r=neg(a.value());
        for(SizeType i=0; i!=n; ++i) {
            r[i]=neg(a.gradient(i));
        }
        r.set_error( a.error() );

        return r;
    }
    static AffineModelType apply(Add, const AffineModelType& a1, const AffineModelType& a2) {
        ARIADNE_ASSERT_MSG(a1.argument_size()==a2.argument_size(),"a1="<<a1<<" a2="<<a2);
        SizeType n=a1.argument_size();
        PrecisionType prec=max(a1.precision(),a2.precision());

        AffineModelType r(n,prec);
        r.value()=add_err(a1.value(),a2.value(),r.error());
        for(SizeType i=0; i!=n; ++i) {
            r.gradient(i)=add_err(a1.gradient(i),a2.gradient(i),r.error());
        }
        r.error()+=(a1.error()+a2.error());
        return r;
    }

    static AffineModelType apply(Mul, const AffineModelType& a1, const AffineModelType& a2) {
        ARIADNE_ASSERT_MSG(a1.argument_size()==a2.argument_size(),"a1="<<a1<<" a2="<<a2);
        SizeType n=a1.argument_size();
        PrecisionType prec=max(a1.precision(),a2.precision());

        AffineModelType r(n,prec);

        r.value()=mul_err(a1.value(),a2.value(),r.error());
        for(SizeType i=0; i!=n; ++i) {
            r.gradient(i)=lin_err(a1.value(),a1.gradient(i),a2.value(),a2.gradient(i),r.error());
        }

        r.error() = (r.error()+a1.error()*a2.error())
                        + ( (mag(a1.value())+norm(a1.gradient()) * a2.error()) + a1.error()*(mag(a2.value())+norm(a2.gradient())) )
                        + norm(a1.gradient()) * norm(a2.gradient());

        return r;
    }

    static AffineModelType apply(Add, const NumericType& c, const AffineModelType& a) {
        AffineModelType r=a;
        r.value()=add_err(a.value(),c,r.error());
        return r;
    }

    static AffineModelType apply(Add, const AffineModelType& a, const NumericType& c) {
        return apply(Add(), c,a);
    }

    static AffineModelType apply(Mul, const NumericType& c, const AffineModelType& a) {
        SizeType n=a.argument_size();
        PrecisionType prec=a.precision();

        AffineModelType r(n,prec);
        r.value()=mul_err(a.value(),c,r.error());
        for(SizeType i=0; i!=n; ++i) {
            r.gradient(i)=mul_err(a.gradient(i),c,r.error());
        }
        r.error()+=mag(c)*a.error();
        return r;
    }

    static AffineModelType apply(Mul, const AffineModelType& a, const NumericType& c) {
        return apply(Mul(), c,a);
    }
};


template<class FLT> AffineModel<ValidatedTag,FLT>::AffineModel(const Affine<NumericType>& affine)
    : AffineModel(affine.argument_size(),affine.value().precision())
{
    AffineModel<ValidatedTag,FLT>& affine_model=*this;
    affine_model = affine.value().value();
    for(SizeType j=0; j!=affine.argument_size(); ++j) {
        affine_model[j] = affine[j].value();
    }
    FLT::set_rounding_upward();
    FLT e(affine_model.error().precision());
    for(SizeType j=0; j!=affine.argument_size(); ++j) {
        e = add(rounded,e,max(sub(rounded,affine.gradient(j).upper().raw(),affine_model.gradient(j)),
                              sub(rounded,affine_model.gradient(j),affine.gradient(j).lower().raw())));
    }
    e = add(rounded,e,max(sub(rounded,affine.value().upper().raw(),affine_model.value()),
                          sub(rounded,affine_model.value(),affine.value().lower().raw())));
    affine_model.set_error(ErrorValueType(e));
    FLT::set_rounding_to_nearest();
}

template<class FLT> AffineModel<ValidatedTag,FLT>::AffineModel(const TaylorModel<ValidatedTag,FLT>& taylor_model)
    : AffineModel(taylor_model.argument_size(),taylor_model.value().precision())
{
    AffineModel<ValidatedTag,FLT>& affine_model=*this;

    typename FLT::RoundingModeType rnd=FLT::get_rounding_mode();
    FLT::set_rounding_upward();
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
    FLT::set_rounding_mode(rnd);
}

template<class P, class PR> AffineModel<P,RawFloatType<PR>> affine_model(const BoxDomainType& domain, const ScalarMultivariateFunction<P>& function, PR precision)
{
    ValidatedScalarMultivariateTaylorFunctionModelDP tf(domain,function,AffineSweeper<RawFloat<PR>>(precision));
    return affine_model(tf.model());
}

inline decltype(auto) operator<(FloatMPApproximation x,const FloatDP y) { return x<FloatMPApproximation(y.get_d(),x.precision()); }
inline decltype(auto) operator<(FloatMP x,const FloatDP y) { return x<FloatMP(Dyadic(y.get_d()),x.precision()); }

// FIXME: Clean up construction of interval with precision PrecisionType<FLT>
template<class FLT> AffineModel<ApproximateTag,FLT> AffineModel<ApproximateTag,FLT>::scaling(SizeType n, SizeType j, const IntervalDomainType& codom, PrecisionType pr)
{
    AffineModel<ApproximateTag,FLT> r(n,pr);
    FloatApproximation<PR> l(codom.lower_bound().get_d(),pr);
    FloatApproximation<PR> u(codom.upper_bound().get_d(),pr);
    r.set_gradient(j,hlf(u-l));
    r.set_value(hlf(l+u));
    return r;
}

template<class FLT> Vector<AffineModel<ApproximateTag,FLT>> AffineModel<ApproximateTag,FLT>::scalings(const BoxDomainType& codom, PrecisionType pr)
{
    SizeType n=codom.dimension();
    Vector<AffineModel<ApproximateTag,FLT>> r(n,AffineModel<ApproximateTag,FLT>(n,pr));
    for(SizeType i=0; i!=n; ++i) {
        r[i] = AffineModel<ApproximateTag,FLT>::scaling(n,i,codom[i],pr);
    }
    return r;
}

template<class FLT> AffineModel<ValidatedTag,FLT> AffineModel<ValidatedTag,FLT>::scaling(SizeType n, SizeType j, const IntervalDomainType& codom, PrecisionType pr)
{
    AffineModel<ValidatedTag,FLT> r(n,pr);
    Float<PR> l(Dyadic(codom.lower_bound()),pr);
    Float<PR> u(Dyadic(codom.upper_bound()),pr);
    Interval<Float<PR>> ivl(l,u);
    r.set_gradient(j,1);
    r*=ivl.radius();
    r+=ivl.midpoint();
    return r;
}

template<class FLT> Vector<AffineModel<ValidatedTag,FLT>> AffineModel<ValidatedTag,FLT>::scalings(const BoxDomainType& codom, PrecisionType pr)
{
    SizeType n=codom.size();
    Vector<AffineModel<ValidatedTag,FLT>> r(n,AffineModel<ValidatedTag,FLT>(n,pr));
    for(SizeType i=0; i!=n; ++i) {
        r[i] = AffineModel<ValidatedTag,FLT>::scaling(n,i,codom[i],pr);
    }
    return r;
}

template<class FLT> auto AffineModel<ApproximateTag,FLT>::range() const -> RangeType
{
    auto v=this->value();
    auto e=this->gradient().zero_element();
    for(SizeType i=0; i!=this->argument_size(); ++i) {
        e+=abs(this->gradient(i));
    }
    return RangeType(v-e,v+e);
}

template<class FLT> auto AffineModel<ValidatedTag,FLT>::range() const -> RangeType
{
    auto v=this->value();
    ErrorValueType e=this->error();
    for(SizeType i=0; i!=this->argument_size(); ++i) {
        e+=abs(this->gradient(i));
    }
    return RangeType(v-e,v+e);
}

template<class FLT> OutputStream& AffineModel<ApproximateTag,FLT>::_write(OutputStream& os) const
{
    os << this->value().raw();
    for(SizeType j=0; j!=this->argument_size(); ++j) {
        if(decide(this->gradient(j)>0)) { os << "+"; }
        if(definitely(this->gradient(j)!=0)) { os << this->gradient(j) << "*x" << j; }
    }
    return os;
}

template<class FLT> OutputStream& AffineModel<ValidatedTag,FLT>::_write(OutputStream& os) const
{
    os << this->value().raw();
    for(SizeType j=0; j!=this->argument_size(); ++j) {
        if(decide(this->gradient(j)>0)) { os << "+"; }
        if(definitely(this->gradient(j)!=0)) { os << this->gradient(j) << "*x" << j; }
    }
    return os << "+/-" << this->error();

}

template<class FLT> auto
AffineModel<ApproximateTag,FLT>::_compose(const ScalarMultivariateFunction<P>& f, Vector<AffineModel<P,FLT>> const& g) -> AffineModel<P,FLT> {
    auto c=values(g);
    auto b=f(c);
    auto A=f.gradient(c);
    return b+A*(g-c);
}

template<class FLT> auto
AffineModel<ApproximateTag,FLT>::_compose(const VectorMultivariateFunction<P>& f, Vector<AffineModel<P,FLT>> const& g) -> Vector<AffineModel<P,FLT>> {
    auto c=values(g);
    auto b=f(c);
    auto A=f.jacobian(c);
    return b+A*(g-c);
}

template<class FLT> auto
AffineModel<ValidatedTag,FLT>::_compose(const ScalarMultivariateFunction<P>& f, Vector<AffineModel<P,FLT>> const& g) -> AffineModel<P,FLT> {
    auto d = ranges(g);
    auto r = cast_singleton(d);
    auto c=values(g);
    auto b=f(c);
    auto A=f.gradient(r);
    return b+A*(g-c);
}

template<class FLT> auto
AffineModel<ValidatedTag,FLT>::_compose(const VectorMultivariateFunction<P>& f, Vector<AffineModel<P,FLT>> const& g) -> Vector<AffineModel<P,FLT>> {
    auto d = ranges(g);
    auto r = cast_singleton(d);
    auto c=values(g);
    auto b=f(c);
    auto A=f.jacobian(r);
    return b+A*(g-c);
}







} //namespace Ariadne


