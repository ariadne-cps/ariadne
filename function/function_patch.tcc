/***************************************************************************
 *            function_patch.tcc
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

#include "function/functional.h"
#include "config.h"

#include <iostream>
#include <iomanip>

#include "utility/macros.h"
#include "utility/exceptions.h"
#include "numeric/numeric.h"
#include "algebra/vector.h"
#include "algebra/matrix.h"
#include "algebra/multi_index.h"
#include "function/polynomial.h"
#include "algebra/differential.h"

#include "function/function.h"
#include "function/function_mixin.h"

#include "algebra/evaluate.h"

#define VOLATILE ;

namespace Ariadne {

template<class M> Void _set_scaling(FunctionPatch<M>& x, const ExactInterval& ivl, SizeType j)
{
    rounding_mode_t rounding_mode=get_rounding_mode();
    set_rounding_mode(upward);
    const Float& l=ivl.lower().raw();
    const Float& u=ivl.upper().raw();
    VOLATILE Float pc=u; pc+=l;
    VOLATILE Float nc=-u; nc-=l;
    VOLATILE Float pg=u; pg-=l;
    VOLATILE Float ng=l; ng-=u;
    x.error()=ErrorType((pc+nc+pg+ng)/4);
    set_rounding_mode(to_nearest);
    MultiIndex a(x.argument_size());
    x.expansion().raw().append(a,(l+u)/2);
    ++a[j];
    x.expansion().raw().append(a,(l+u)/2);
    set_rounding_mode(rounding_mode);
}



inline OutputStream& operator<<(OutputStream& os, const Representation<Float>& flt_repr)
{
    const Float& flt=*flt_repr.pointer;
    Int precision=os.precision(); std::ios_base::fmtflags flags = os.flags();
    os.precision(17); os.setf(std::ios_base::showpoint);
    os << "Float(" << flt << ")";
    os.precision(precision); os.flags(flags);
    return os;
}

inline OutputStream& operator<<(OutputStream& os, const Representation<UpperFloat>& flt_repr)
{
    return os << reinterpret_cast<Representation<Float>const&>(flt_repr);
}

inline OutputStream& operator<<(OutputStream& os, const Representation<ExactInterval>& ivl_repr)
{
    const ExactInterval& ivl=*ivl_repr.pointer;
    Int precision=os.precision(); std::ios_base::fmtflags flags = os.flags();
    os.precision(17); os.setf(std::ios_base::showpoint);
    os << "ExactInterval("<<ivl.lower()<<","<<ivl.upper()<<")";
    os.precision(precision); os.flags(flags);
    return os;
}


inline OutputStream& operator<<(OutputStream& os, const Representation< Expansion<Float> >& exp_repr)
{
    const Expansion<Float>& exp=*exp_repr.pointer;
    Int precision=os.precision(); std::ios_base::fmtflags flags = os.flags();
    os.precision(17); os.setf(std::ios_base::showpoint);
    os << "Expansion<Float>(" << exp.argument_size() << "," << exp.number_of_nonzeros();
    for(Expansion<Float>::ConstIterator iter=exp.begin(); iter!=exp.end(); ++iter) {
        for(SizeType j=0; j!=iter->key().size(); ++j) {
            os << "," << Nat(iter->key()[j]);
        }
        os << "," << iter->data();
    }
    os << ")";
    os.precision(precision); os.flags(flags);
    return os;
}

inline OutputStream& operator<<(OutputStream& os, const Representation< Expansion<ExactFloat> >& exp_repr) {
    return os << reinterpret_cast<Expansion<Float>const&>(exp_repr);
}

inline OutputStream& operator<<(OutputStream& os, const Representation< Expansion<ApproximateFloat> >& exp_repr) {
    return os << reinterpret_cast<Expansion<Float>const&>(exp_repr);
}

template<class X> inline OutputStream& operator<<(OutputStream& os, const Representation< Vector<X> >& vec_repr)
{
    const Vector<X>& vec=*vec_repr.pointer;
    ARIADNE_ASSERT(vec.size()!=0);
    os << "(";
    for(SizeType i=0; i!=vec.size(); ++i) {
        if(i!=0) { os << ","; }
        os << representation(vec[i]);
    }
    os << ")";
    return os;
}

inline OutputStream& operator<<(OutputStream& os, const Representation< ExactBox >& box_repr)
{
    const Vector<ExactInterval>& vec=*box_repr.pointer;
    ARIADNE_ASSERT(vec.size()!=0);
    os << "(";
    for(SizeType i=0; i!=vec.size(); ++i) {
        if(i!=0) { os << ","; }
        os << representation(vec[i]);
    }
    os << ")";
    return os;
}

inline OutputStream& operator<<(OutputStream& os, const Representation< Sweeper >& swp_repr)
{
    const Sweeper& swp=*swp_repr.pointer;
    Int precision=os.precision(); std::ios_base::fmtflags flags = os.flags();
    os.precision(17); os.setf(std::ios_base::showpoint);
    os << swp;
    os.precision(precision); os.flags(flags);
    return os;
}

template<class T> inline OutputStream& operator<<(OutputStream& os, const Representation< List<T> >& lst_repr)
{
    const List<T>& lst=*lst_repr.pointer;
    ARIADNE_ASSERT(lst.size()!=0);
    os << "(";
    for(SizeType i=0; i!=lst.size(); ++i) {
        if(i!=0) { os << ","; }
        os << representation(lst[i]);
    }
    os << ")";
    return os;
}



template<class M> FunctionPatch<M>::FunctionPatch()
    : _domain(), _model()
{ }

template<class M> FunctionPatch<M>::FunctionPatch(const ExactBox& d, Sweeper swp)
    : _domain(d), _model(d.size(),swp)
{
}

template<class M> FunctionPatch<M>::FunctionPatch(const ExactBox& d, const Expansion<RawFloat>& p, const RawFloat& e, const Sweeper& swp)
    : _domain(d), _model(p,e,swp)
{
}

template<class M> FunctionPatch<M>::FunctionPatch(const ExactBox& d, const Expansion<ExactFloat>& p, const ErrorFloat& e, const Sweeper& swp)
    : _domain(d), _model(p,e,swp)
{
}

template<class M> FunctionPatch<M>::FunctionPatch(const ExactBox& d, const ModelType& m)
    : _domain(d), _model(m)
{
}

template<class M> FunctionPatch<M>::FunctionPatch(const ExactBox& d, const ScalarFunctionType<M>& f, Sweeper swp)
    : _domain(d), _model(f.argument_size(),swp)
{
    ARIADNE_ASSERT_MSG(d.size()==f.argument_size(),"d="<<d<<" f="<<f);
    Vector<ModelType> x=ModelType::scalings(d,swp);
    this->_model=f.evaluate(x);
    this->_model.sweep();
}

template<class M> FunctionPatch<M>& FunctionPatch<M>::operator=(const ValidatedScalarFunctionModel& f)
{
    return (*this)=FunctionPatch<M>(this->domain(),f,this->sweeper());
}



template<class M> FunctionPatch<M> FunctionPatch<M>::zero(const ExactBox& d, Sweeper swp)
{
    return FunctionPatch<M>(d,ModelType::zero(d.size(),swp));
}

template<class M> FunctionPatch<M> FunctionPatch<M>::constant(const ExactBox& d, const NumericType& c, Sweeper swp)
{
    return FunctionPatch<M>(d,ModelType::constant(d.size(),c,swp));
}

template<class M> FunctionPatch<M> FunctionPatch<M>::coordinate(const ExactBox& d, SizeType j, Sweeper swp)
{
    ARIADNE_ASSERT(j<d.size());
    return FunctionPatch<M>(d,ModelType::scaling(d.size(),j,d[j],swp));
}

template<class M> VectorFunctionPatch<M> FunctionPatch<M>::identity(const ExactBox& d, Sweeper swp)
{
    return VectorFunctionPatch<M>(d,ModelType::scalings(d,swp));
}


template<class M> Vector<FunctionPatch<M>> FunctionPatch<M>::constants(const ExactBox& d, const Vector<NumericType>& c, Sweeper swp)
{
    ARIADNE_DEPRECATED("FunctionPatch<M>::constants","Use VectorFunctionPatch<M>::constant instead");
    Vector<FunctionPatch<M>> x(c.size(),FunctionPatch<M>(d,swp));
    for(SizeType i=0; i!=c.size(); ++i) {
        x[i]=c[i];
    }
    return x;
}

template<class M> Vector<FunctionPatch<M>> FunctionPatch<M>::coordinates(const ExactBox& d, Sweeper swp)
{
    ARIADNE_DEPRECATED("FunctionPatch<M>::coordinates","Use VectorFunctionPatch<M>::identity instead");
    Vector<FunctionPatch<M>> x(d.dimension(),FunctionPatch<M>(d,swp));
    for(SizeType i=0; i!=x.size(); ++i) {
        x[i]=FunctionPatch<M>::coordinate(d,i,swp);
    }
    return x;
}

template<class M> Vector<FunctionPatch<M>> FunctionPatch<M>::coordinates(const ExactBox& d, SizeType imin, SizeType imax, Sweeper swp)
{
    ARIADNE_DEPRECATED("FunctionPatch<M>::coordinates","Use VectorFunctionPatch<M>::projection instead");
    ARIADNE_ASSERT(imin<=imax);
    ARIADNE_ASSERT(imax<=d.size());

    Vector<FunctionPatch<M>> x(imax-imin);
    for(SizeType i=imin; i!=imax; ++i) {
        x[i-imin]=FunctionPatch<M>::coordinate(d,i,swp);
    }
    return x;
}


template<class M> FunctionPatch<M> FunctionPatch<M>::create_zero() const
{
    return FunctionPatch<M>(this->domain(),this->_model.sweeper());
}

template<class M> FunctionPatch<M>* FunctionPatch<M>::_clone() const
{
    return new FunctionPatch<M>(*this);
}

template<class M> FunctionPatch<M>* FunctionPatch<M>::_create() const
{
    return new FunctionPatch<M>(this->domain(),this->_model.sweeper());
}

template<class M> VectorFunctionModelInterface<typename M::Paradigm>* FunctionPatch<M>::_create_identity() const
{
    Sweeper sweeper=this->sweeper();
    VectorFunctionPatch<M>* result = new VectorFunctionPatch<M>(this->domain().size(), FunctionPatch<M>(this->domain(),sweeper));
    for(SizeType i=0; i!=result->size(); ++i) { (*result)[i]=FunctionPatch<M>::coordinate(this->domain(),i,sweeper); }
    return result;
}

template<class M> VectorFunctionModelInterface<typename M::Paradigm>* FunctionPatch<M>::_create_vector(SizeType i) const
{
    return new VectorFunctionPatch<M>(i,this->domain(),this->_model.sweeper());
}


// To scale from a model on [a,b] to a model on [c,d], use scale factor s=(d-c)/(b-a)
// and translation t=((c+d)-(a+b))/(b-a)
// Because we are scaling the model on [-1,+1], this is not the same as
// the mapping taking [a,b] to [c,d]
template<class M> FunctionPatch<M> partial_restriction(const FunctionPatch<M>& tv, SizeType k, const ExactInterval& new_ivl) {
    ARIADNE_ASSERT(k<tv.argument_size())
    const ExactInterval& old_ivl=tv.domain()[k];
    ARIADNE_ASSERT(subset(new_ivl,old_ivl));
    if(new_ivl==old_ivl) { return tv; }
    Float a=old_ivl.lower().raw(); Float b=old_ivl.upper().raw();
    Float c=new_ivl.lower().raw(); Float d=new_ivl.upper().raw();
    if(a==b) { ARIADNE_ASSERT( a<b || (a==b && c==d) ); return tv; }
    ValidatedNumber s=static_cast<ValidatedNumber>(sub_ivl(d,c)/sub_ivl(b,a));
    // ValidatedNumber t=(mul_ivl(b,c)-mul_ivl(a,d))/sub_ivl(b,a);  // WRONG!!
    ValidatedNumber t=static_cast<ValidatedNumber>((add_ivl(c,d)-add_ivl(a,b))/sub_ivl(b,a));
    ExactBox new_dom=tv.domain();
    new_dom[k]=new_ivl;
    return FunctionPatch<M>(new_dom,preaffine(tv.model(),k,s,t));
}

template<class M> FunctionPatch<M> restriction(const FunctionPatch<M>& tv, const ExactBox& d) {
    ARIADNE_ASSERT(subset(d,tv.domain()));
    const ExactBox& od=tv.domain();
    FunctionPatch<M> r=tv;
    for(SizeType j=0; j!=d.size(); ++j) {
        if(od[j]!=d[j]) { r=partial_restriction(r,j,d[j]); }
    }
    return r;
}

template<class M> Void FunctionPatch<M>::restrict(const ExactBox& d) {
    (*this)=restriction(*this,d);
}

template<class M> FunctionPatch<M> extension(const FunctionPatch<M>& tv, const ExactBox& d) {
    const ExactBox& domain=tv.domain();
    ARIADNE_ASSERT(subset(domain,d));
    for(SizeType i=0; i!=d.size(); ++i) {
        ARIADNE_ASSERT(domain[i]==d[i] || domain[i].lower()==domain[i].upper());
    }
    return FunctionPatch<M>(d,tv.model());
}



template<class M> Polynomial<ValidatedFloat> polynomial(const M& tm);

inline Bool operator==(ExactFloat x1, Int n2) { return x1.raw()==Float(n2); }
inline Bool operator==(ValidatedFloat x1, Int n2) { return x1.upper_raw()==Float(n2) && x1.lower_raw()==Float(n2); }
inline Bool operator==(ApproximateFloat x1, Int n2) { return x1.raw()==Float(n2); }

inline Bool operator!=(ExactFloat x1, Int n2) { return x1.raw()!=Float(n2); }
inline Bool operator!=(ValidatedFloat x1, Int n2) { return x1.upper_raw()!=Float(n2) || x1.lower_raw()!=Float(n2); }
inline Bool operator!=(ApproximateFloat x1, Int n2) { return x1.raw()!=Float(n2); }

inline Bool operator> (ExactFloat x1, Int n2) { return x1.raw()> Float(n2); }
inline Bool operator> (ValidatedFloat x1, Int n2) { return x1.lower_raw()> Float(n2); }
inline Bool operator> (ApproximateFloat x1, Int n2) { return x1.raw()> Float(n2); }

template<class M> Polynomial<ValidatedFloat> FunctionPatch<M>::polynomial() const
{
    Polynomial<ValidatedFloat> z(this->argument_size());

    Polynomial<ValidatedFloat> p=Ariadne::polynomial(this->model());

    Vector<Polynomial<ValidatedFloat> > s(this->argument_size(),z);
    for(SizeType j=0; j!=this->argument_size(); ++j) {
        ExactInterval const& domj=this->domain()[j];
        if(domj.width()<=0) {
            ARIADNE_ASSERT(this->domain()[j].width()==0);
            s[j]=Polynomial<ValidatedFloat>::constant(this->argument_size(),0);
        } else {
            //s[j]=Ariadne::polynomial(ModelType::unscaling(this->argument_size(),j,this->domain()[j],this->sweeper()));
            s[j]=(Polynomial<ValidatedFloat>::coordinate(this->argument_size(),j)-domj.midpoint())/domj.radius();
        }
    }

    return compose(p,s);
}

template<class M> ScalarFunctionType<M> FunctionPatch<M>::function() const
{
    return ValidatedScalarFunction(new FunctionPatch<M>(*this));
}


template<class M> Bool FunctionPatch<M>::operator==(const FunctionPatch<M>& tv) const
{
    return this->_domain==tv._domain && this->_model==tv._model;
}



template<class M> FunctionPatch<M>& operator+=(FunctionPatch<M>& f, const FunctionPatch<M>& g) {
    ARIADNE_ASSERT_MSG(subset(f.domain(),g.domain()),"f="<<f<<", g="<<g);
    if(f.domain()==g.domain()) { f.model()+=g.model(); }
    else { f.model()+=restriction(g,f.domain()).model(); }
    return f;
}

template<class M> FunctionPatch<M>& operator-=(FunctionPatch<M>& f, const FunctionPatch<M>& g) {
    ARIADNE_ASSERT_MSG(subset(f.domain(),g.domain()),"f="<<f<<", g="<<g);
    if(f.domain()==g.domain()) { f.model()-=g.model(); }
    else { f.model()-=restriction(g,f.domain()).model(); }
    return f;
}


template<class M> FunctionPatch<M>& operator+=(FunctionPatch<M>& f, const NumericType<M>& c) {
    f.model()+=c;
    return f;
}

template<class M> FunctionPatch<M>& operator-=(FunctionPatch<M>& f, const NumericType<M>& c) {
    f.model()-=c;
    return f;
}

template<class M> FunctionPatch<M>& operator*=(FunctionPatch<M>& f, const NumericType<M>& c) {
    f.model()*=c;
    return f;
}

template<class M> FunctionPatch<M>& operator/=(FunctionPatch<M>& f, const NumericType<M>& c) {
    f.model()/=c;
    return f;
}

template<class M> FunctionPatch<M> operator+(const FunctionPatch<M>& x) {
    return FunctionPatch<M>(x.domain(),x.model());
}

template<class M> FunctionPatch<M> operator-(const FunctionPatch<M>& x) {
    return FunctionPatch<M>(x.domain(),-x.model());
}


template<class M> FunctionPatch<M> operator+(const FunctionPatch<M>& x1, const FunctionPatch<M>& x2) {
    if(x1.domain()==x2.domain()) {
        return FunctionPatch<M>(x1.domain(),x1.model()+x2.model()); }
    else {
        ExactBox domain=intersection(x1.domain(),x2.domain());
        return FunctionPatch<M>(domain,restriction(x1,domain).model()+restriction(x2,domain).model());}
}

template<class M> FunctionPatch<M> operator-(const FunctionPatch<M>& x1, const FunctionPatch<M>& x2) {
    if(x1.domain()==x2.domain()) {
        return FunctionPatch<M>(x1.domain(),x1.model()-x2.model()); }
    else {
        ExactBox domain=intersection(x1.domain(),x2.domain());
        return FunctionPatch<M>(domain,restriction(x1,domain).model()-restriction(x2,domain).model());}
}

template<class M> FunctionPatch<M> operator*(const FunctionPatch<M>& x1, const FunctionPatch<M>& x2) {
    if(x1.domain()==x2.domain()) {
        return FunctionPatch<M>(x1.domain(),x1.model()*x2.model()); }
    else {
        ExactBox domain=intersection(x1.domain(),x2.domain());
        return FunctionPatch<M>(domain,restriction(x1,domain).model()*restriction(x2,domain).model());}
}

template<class M> FunctionPatch<M> operator/(const FunctionPatch<M>& x1, const FunctionPatch<M>& x2) {
    if(x1.domain()==x2.domain()) {
        return FunctionPatch<M>(x1.domain(),x1.model()/x2.model()); }
    else {
        ExactBox domain=intersection(x1.domain(),x2.domain());
        return FunctionPatch<M>(domain,restriction(x1,domain).model()/restriction(x2,domain).model());}
}

template<class M> FunctionPatch<M> operator+(const FunctionPatch<M>& x, const NumericType<M>& c) {
    FunctionPatch<M> r(x); r+=c; return r;
}

template<class M> FunctionPatch<M> operator-(const FunctionPatch<M>& x, const NumericType<M>& c) {
    FunctionPatch<M> r(x); r+=neg(c); return r;
}

template<class M> FunctionPatch<M> operator*(const FunctionPatch<M>& x, const NumericType<M>& c) {
    FunctionPatch<M> r(x); r*=c; return r;
}

template<class M> FunctionPatch<M> operator/(const FunctionPatch<M>& x, const NumericType<M>& c) {
    FunctionPatch<M> r(x); r*=rec(c); return r;
}

template<class M> FunctionPatch<M> operator+(const NumericType<M>& c, const FunctionPatch<M>& x) {
    FunctionPatch<M> r(x); r+=c; return r;
}

template<class M> FunctionPatch<M> operator-(const NumericType<M>& c, const FunctionPatch<M>& x) {
    FunctionPatch<M> r(neg(x)); r+=c; return r;
}

template<class M> FunctionPatch<M> operator*(const NumericType<M>& c, const FunctionPatch<M>& x) {
    FunctionPatch<M> r(x); r*=c; return r;
}

template<class M> FunctionPatch<M> operator/(const NumericType<M>& c, const FunctionPatch<M>& x) {
    FunctionPatch<M> r(rec(x)); r*=c; return r;
}

template<class M> FunctionPatch<M> operator+(const FunctionType<M>& f1, const FunctionPatch<M>& tf2) {
    return FunctionPatch<M>(tf2.domain(),f1,tf2.sweeper())+tf2; }
template<class M> FunctionPatch<M> operator-(const FunctionType<M>& f1, const FunctionPatch<M>& tf2) {
    return FunctionPatch<M>(tf2.domain(),f1,tf2.sweeper())-tf2; }
template<class M> FunctionPatch<M> operator*(const FunctionType<M>& f1, const FunctionPatch<M>& tf2) {
    return FunctionPatch<M>(tf2.domain(),f1,tf2.sweeper())*tf2; }
template<class M> FunctionPatch<M> operator/(const FunctionType<M>& f1, const FunctionPatch<M>& tf2) {
    return FunctionPatch<M>(tf2.domain(),f1,tf2.sweeper())/tf2; }
template<class M> FunctionPatch<M> operator+(const FunctionPatch<M>& tf1, const FunctionType<M>& f2) {
    return tf1+FunctionPatch<M>(tf1.domain(),f2,tf1.sweeper()); }
template<class M> FunctionPatch<M> operator-(const FunctionPatch<M>& tf1, const FunctionType<M>& f2) {
    return tf1-FunctionPatch<M>(tf1.domain(),f2,tf1.sweeper()); }
template<class M> FunctionPatch<M> operator*(const FunctionPatch<M>& tf1, const FunctionType<M>& f2) {
    return tf1*FunctionPatch<M>(tf1.domain(),f2,tf1.sweeper()); }
template<class M> FunctionPatch<M> operator/(const FunctionPatch<M>& tf1, const FunctionType<M>& f2) {
    return tf1/FunctionPatch<M>(tf1.domain(),f2,tf1.sweeper()); }






template<class M> FunctionPatch<M> max(const FunctionPatch<M>& x1, const FunctionPatch<M>& x2) {
    if(x1.domain()==x2.domain()) {
        return FunctionPatch<M>(x1.domain(),max(x1.model(),x2.model())); }
    else {
        ExactBox domain=intersection(x1.domain(),x2.domain());
        return FunctionPatch<M>(domain,max(restriction(x1,domain).model(),restriction(x2,domain).model()));}
}

template<class M> FunctionPatch<M> min(const FunctionPatch<M>& x1, const FunctionPatch<M>& x2) {
    if(x1.domain()==x2.domain()) {
        return FunctionPatch<M>(x1.domain(),min(x1.model(),x2.model())); }
    else {
        ExactBox domain=intersection(x1.domain(),x2.domain());
        return FunctionPatch<M>(domain,min(restriction(x1,domain).model(),restriction(x2,domain).model()));}
}

template<class M> FunctionPatch<M> abs(const FunctionPatch<M>& x) {
    return FunctionPatch<M>(x.domain(),abs(x.model())); }
template<class M> FunctionPatch<M> neg(const FunctionPatch<M>& x) {
    return FunctionPatch<M>(x.domain(),-x.model()); }
template<class M> FunctionPatch<M> rec(const FunctionPatch<M>& x) {
    return FunctionPatch<M>(x.domain(),rec(x.model())); }
template<class M> FunctionPatch<M> sqr(const FunctionPatch<M>& x) {
    return FunctionPatch<M>(x.domain(),sqr(x.model())); }
template<class M> FunctionPatch<M> pow(const FunctionPatch<M>& x, Int n) {
    return FunctionPatch<M>(x.domain(),pow(x.model(),n)); }
template<class M> FunctionPatch<M> sqrt(const FunctionPatch<M>& x) {
    return FunctionPatch<M>(x.domain(),sqrt(x.model())); }
template<class M> FunctionPatch<M> exp(const FunctionPatch<M>& x) {
    return FunctionPatch<M>(x.domain(),exp(x.model())); }
template<class M> FunctionPatch<M> log(const FunctionPatch<M>& x) {
    return FunctionPatch<M>(x.domain(),log(x.model())); }
template<class M> FunctionPatch<M> sin(const FunctionPatch<M>& x) {
    return FunctionPatch<M>(x.domain(),sin(x.model())); }
template<class M> FunctionPatch<M> cos(const FunctionPatch<M>& x) {
    return FunctionPatch<M>(x.domain(),cos(x.model())); }
template<class M> FunctionPatch<M> tan(const FunctionPatch<M>& x) {
    return FunctionPatch<M>(x.domain(),tan(x.model())); }
template<class M> FunctionPatch<M> asin(const FunctionPatch<M>& x) {
    ARIADNE_NOT_IMPLEMENTED; }
template<class M> FunctionPatch<M> acos(const FunctionPatch<M>& x) {
    ARIADNE_NOT_IMPLEMENTED; }
template<class M> FunctionPatch<M> atan(const FunctionPatch<M>& x) {
    ARIADNE_NOT_IMPLEMENTED; }


template<class M> FunctionPatch<M> antiderivative(const FunctionPatch<M>& x, SizeType k) {
    ValidatedNumber sf=rad_val(x.domain()[k]);
    return FunctionPatch<M>(x.domain(),antiderivative(x.model(),k)*sf); }

template<class M> FunctionPatch<M> derivative(const FunctionPatch<M>& x, SizeType k) {
    ValidatedNumber sf=rec(rad_val(x.domain()[k]));
    return FunctionPatch<M>(x.domain(),derivative(x.model(),k)*sf); }

template<class M> FunctionPatch<M> embed(const ExactBox& dom1, const FunctionPatch<M>& tv2,const ExactBox& dom3) {
    return FunctionPatch<M>(product(dom1,tv2.domain(),dom3),embed(dom1.size(),tv2.model(),dom3.size())); }

template<class M> FunctionPatch<M>* FunctionPatch<M>::_derivative(SizeType j) const
{
    return new FunctionPatch<M>(Ariadne::derivative(*this,j));
}


template<class M> ApproximateNumber FunctionPatch<M>::operator() (const Vector<ApproximateNumber>& x) const
{
    const FunctionPatch<M>& f=*this;
    if(!contains(f.domain(),make_exact(x))) {
        ARIADNE_THROW(DomainException,"tf.evaluate(ax) with tf="<<f<<", ax="<<x," ax is not an element of tf.domain()="<<f.domain());
    }
    Vector<ApproximateNumber> sx=Ariadne::unscale(x,f._domain);
    return Ariadne::evaluate(this->_model.expansion(),sx);
}

template<class M> ValidatedNumber FunctionPatch<M>::operator()(const Vector<ValidatedNumber>& x) const
{
    const FunctionPatch<M>& f=*this;
    if(!contains(f.domain(),x)) {
        ARIADNE_THROW(DomainException,"evaluate(f,x) with f="<<f<<", x="<<x,"x is not a subset of f.domain()="<<f.domain());
    }
    return unchecked_evaluate(f,x);
}

template<class M> ValidatedNumber FunctionPatch<M>::operator()(const Vector<ExactNumber>& x) const
{
    return Ariadne::evaluate(*this,Vector<ValidatedNumber>(x));
}

template<class M> Covector<NumericType<M>> FunctionPatch<M>::gradient(const Vector<NumericType>& x) const
{
    Covector<NumericType> g=Ariadne::gradient(this->_model,unscale(x,this->_domain));
    for(SizeType j=0; j!=g.size(); ++j) {
        NumericType rad=rad_val(this->_domain[j]);
        g[j]/=rad;
    }
    return g;
}



template<class M> NumericType<M> unchecked_evaluate(const FunctionPatch<M>& f, const Vector<NumericType<M>>& x) {
    return evaluate(f.model(),unscale(x,f.domain()));
}


template<class M> FunctionPatch<M> compose(const FunctionType<M>& g, const VectorFunctionPatch<M>& f) {
    return FunctionPatch<M>(f.domain(),g.evaluate(f.models()));
}

template<class M> FunctionPatch<M> compose(const FunctionPatch<M>& g, const VectorFunctionPatch<M>& f)
{
    if(!subset(f.codomain(),g.domain())) {
        ARIADNE_THROW(DomainException,"compose(g,f) with g="<<g<<", f="<<f,"f.codomain()="<<f.codomain()<<" is not a subset of g.domain()="<<g.domain());
    }
    return unchecked_compose(g,f);
}

template<class M> FunctionPatch<M> unchecked_compose(const FunctionPatch<M>& g, const VectorFunctionPatch<M>& f)
{
    return FunctionPatch<M>(f.domain(),compose(g.model(),unscale(f.models(),g.domain())));
}



template<class M> FunctionPatch<M> partial_evaluate(const FunctionPatch<M>& te, SizeType k, const NumericType<M>& c)
{
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



template<class M> FunctionPatch<M> antiderivative(const FunctionPatch<M>& f, SizeType k, ExactFloat c)
{
    ARIADNE_ASSERT(k<f.argument_size());
    ARIADNE_ASSERT(contains(f.domain()[k],c));

    FunctionPatch<M> g = antiderivative(f,k);
    VectorFunctionPatch<M> h = VectorFunctionPatch<M>::identity(g.domain(),g.sweeper());
    h[k] = FunctionPatch<M>::constant(g.domain(),c,g.sweeper());

    return g-compose(g,h);
}





template<class M> Pair<FunctionPatch<M>,FunctionPatch<M>> split(const FunctionPatch<M>& tv, SizeType j)
{
    typedef M ModelType;
    Pair<ModelType,ModelType> models=split(tv.model(),j);
    Pair<ExactBox,ExactBox> subdomains=split(tv.domain(),j);
    return make_pair(FunctionPatch<M>(subdomains.first,models.first),
                     FunctionPatch<M>(subdomains.second,models.second));

}

template<class M> Bool refines(const FunctionPatch<M>& tv1, const FunctionPatch<M>& tv2)
{
    if(tv1.domain()==tv2.domain()) { return refines(tv1.model(),tv2.model()); }
    if(subset(tv2.domain(),tv1.domain())) { return refines(restriction(tv1,tv2.domain()).model(),tv2.model()); }
    else { return false; }
}

template<class M> Bool inconsistent(const FunctionPatch<M>& tv1, const FunctionPatch<M>& tv2)
{
    if(tv1.domain()==tv2.domain()) {
        return inconsistent(tv1.model(),tv2.model());
    } else {
        ExactBox domain=intersection(tv1.domain(),tv2.domain());
        return inconsistent(restriction(tv1,domain).model(),restriction(tv2,domain).model());
    }
}

template<class M> FunctionPatch<M> refinement(const FunctionPatch<M>& tv1, const FunctionPatch<M>& tv2)
{
    ARIADNE_ASSERT(tv1.domain()==tv2.domain());
    return FunctionPatch<M>(tv1.domain(),intersection(tv1.model(),tv2.model()));
}

template<class M> ErrorFloat norm(const FunctionPatch<M>& f) {
    return norm(f.model());
}

template<class M> ErrorFloat distance(const FunctionPatch<M>& f1, const FunctionPatch<M>& f2) {
    return norm(f1-f2);
}

template<class M> ErrorFloat distance(const FunctionPatch<M>& f1, const ValidatedScalarFunction& f2) {
    return distance(f1,FunctionPatch<M>(f1.domain(),f2,f1.sweeper()));
}


template<class M> FunctionPatch<M> midpoint(const FunctionPatch<M>& f)
{
    M tm=f.model();
    tm.set_error(0u);
    return FunctionPatch<M>(f.domain(),tm);
}



template<class M> OutputStream& FunctionPatch<M>::write(OutputStream& os) const
{
    Polynomial<ValidatedFloat> p=this->polynomial();
    Polynomial<ApproximateFloat> ap=p;
    os << "FP" << this->domain();
    os << "(";
    os << ap;
    if(this->error()>0.0) { os << "+/-" << this->error(); }
    os << ")";
    return os;
}

template<class M> OutputStream& FunctionPatch<M>::repr(OutputStream& os) const
{
    return os << "FunctionPatch<M>(" << representation(this->domain()) << ", " << representation(this->model().expansion().raw())
              << "," << representation(this->error().raw())<<","<<this->sweeper()<<")";
}

template<class M> OutputStream& operator<<(OutputStream& os, const FunctionPatch<M>& tf)
{
    return tf.write(os);
}

template<class M> OutputStream& operator<<(OutputStream& os, const Representation<FunctionPatch<M>>& tf)
{
    return tf.pointer->repr(os);
}

/*
template<class M> OutputStream& operator<<(OutputStream& os, const Representation<FunctionPatch<M>>& frepr)
{
    FunctionPatch<M> const& function=*frepr.pointer;
    FunctionPatch<M> truncated_function=function;
    truncated_function.set_error(0.0);
    truncated_function.sweep(ThresholdSweeper(TAYLOR_FUNCTION_WRITING_ACCURACY));

    os << midpoint(truncated_function.polynomial());
    if(truncated_function.error()>0.0) { os << "+/-" << truncated_function.error(); }
    if(function.error()>0.0) { os << "+/-" << function.error(); }
    // TODO: Use Unicode +/- literal when this becomes avaialable in C++0x
    return os;
}
*/
/*
template<class M> OutputStream& operator<<(OutputStream& os, const ModelRepresentation<FunctionPatch<M>>& frepr)
{
    FunctionPatch<M> const& f=*frepr.pointer;
    Float truncatation_error = 0.0;
    os << "<"<<f.domain()<<"\n";
    for(ModelType::ConstIterator iter=f.begin(); iter!=f.end(); ++iter) {
        if(abs(iter->data())>frepr.threshold) { truncatation_error+=abs(iter->data()); }
        else { os << iter->key() << ":" << iter->data() << ","; }
    }
    os << "+/-" << truncatation_error << "+/-" << f.error();
    return os;
}
*/

template<class M> OutputStream& operator<<(OutputStream& os, const ModelRepresentation<FunctionPatch<M>>& frepr)
{
    FunctionPatch<M> const& f=*frepr.pointer;
    FunctionPatch<M> tf=f;
    tf.clobber();
    tf.sweep(ThresholdSweeper(frepr.threshold));
    os << "("<<tf.model()<<"+/-"<<f.error();
    return os;
}

template<class M> OutputStream& operator<<(OutputStream& os, const PolynomialRepresentation<FunctionPatch<M>>& frepr)
{
    FunctionPatch<M> const& function=*frepr.pointer;
    FunctionPatch<M> truncated_function = function;
    truncated_function.clobber();
    truncated_function.sweep(ThresholdSweeper(frepr.threshold));
    ErrorFloat truncatation_error = truncated_function.error();
    truncated_function.clobber();
    Polynomial<ValidatedFloat> validated_polynomial_function=polynomial(truncated_function);
    Polynomial<ExactFloat> polynomial_function = midpoint(validated_polynomial_function);
    if(frepr.names.empty()) { os << polynomial_function; }
    else { os << named_argument_repr(polynomial_function,frepr.names); }
    os << "+/-" << truncatation_error << "+/-" << function.error();
    return os;
}





template<class M> Bool check(const Vector<FunctionPatch<M>>& tv)
{
    for(SizeType i=0; i!=tv.size(); ++i) {
        if(tv.zero_element().domain()!=tv[i].domain()) { return false; }
    }
    return true;
}

template<class M> Vector<Expansion<ExactFloat>> expansion(const Vector<FunctionPatch<M>>& x)
{
    Vector< Expansion<ExactFloat> > r(x.size());
    for(SizeType i=0; i!=x.size(); ++i) {
        r[i]=x[i].expansion();
    }
    return r;
}

template<class M> Vector<ErrorFloat> error(const Vector<FunctionPatch<M>>& x)
{
    Vector<ErrorFloat> r(x.size());
    for(SizeType i=0; i!=x.size(); ++i) {
        r[i]=x[i].error();
    }
    return r;
}

template<class M> Vector<ExactFloat> value(const Vector<FunctionPatch<M>>& x)
{
    Vector<ExactFloat> r(x.size());
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







template<class M> VectorFunctionPatch<M>::VectorFunctionPatch()
    : _domain(), _models()
{
}

template<class M> VectorFunctionPatch<M>::VectorFunctionPatch(SizeType k)
    : _domain(), _models(k)
{
}

template<class M> VectorFunctionPatch<M>::VectorFunctionPatch(SizeType m, const ExactBox& d, Sweeper swp)
    : _domain(d), _models(m,ModelType(d.size(),swp))
{
}

template<class M> VectorFunctionPatch<M>::VectorFunctionPatch(SizeType k, const FunctionPatch<M>& f)
    : _domain(f.domain()), _models(k,f.model())
{
}

template<class M> VectorFunctionPatch<M>::VectorFunctionPatch(const ValidatedVectorFunctionModel& f)
    : _domain(), _models()
{
    ARIADNE_ASSERT(dynamic_cast<const VectorFunctionPatch<M>*>(&f.reference()));
    *this = dynamic_cast<const VectorFunctionPatch<M>&>(f.reference());
}

template<class M> VectorFunctionPatch<M>& VectorFunctionPatch<M>::operator=(const ValidatedVectorFunctionModel& f)
{
    ARIADNE_ASSERT(dynamic_cast<const VectorFunctionPatch<M>*>(&f.reference()));
    *this = dynamic_cast<const VectorFunctionPatch<M>&>(f.reference());
    return *this;
}




template<class M> VectorFunctionPatch<M>::VectorFunctionPatch(const ExactBox& d,
                                           const Vector<ModelType>& f)
    : _domain(d), _models(f)
{
    for(SizeType i=0; i!=f.size(); ++i) {
        ARIADNE_ASSERT_MSG(d.size()==f[i].argument_size(),"d="<<d<<", f="<<f);
    }
}

template<class M> VectorFunctionPatch<M>::VectorFunctionPatch(const ExactBox& d,
                                           const Vector<Expansion<ExactFloat>>& f,
                                           const Vector<ErrorFloat>& e,
                                           Sweeper swp)
    : _domain(d), _models(f.size(),ModelType(d.size(),swp))
{
    ARIADNE_ASSERT(f.size()==e.size());
    for(SizeType i=0; i!=f.size(); ++i) {
        ARIADNE_ASSERT(d.size()==f[i].argument_size());
        _models[i]=ModelType(f[i],e[i],swp);
    }
}

template<class M> VectorFunctionPatch<M>::VectorFunctionPatch(const ExactBox& d,
                                           const Vector<Expansion<ExactFloat>>& f,
                                           Sweeper swp)
    : VectorFunctionPatch<M>(d,f,Vector<ErrorFloat>(f.size()),swp)
{
}

template<class M> VectorFunctionPatch<M>::VectorFunctionPatch(const ExactBox& d,
                                           const Vector<Expansion<RawFloat>>& f,
                                           const Vector<RawFloat>& e,
                                           Sweeper swp)
    : VectorFunctionPatch<M>(d,reinterpret_cast<Vector<Expansion<ExactFloat>>const&>(f),
                           reinterpret_cast<Vector<ErrorFloat>const&>(e),swp)
{
}

template<class M> VectorFunctionPatch<M>::VectorFunctionPatch(const ExactBox& d,
                                           const Vector<Expansion<RawFloat>>& f,
                                           Sweeper swp)
    : VectorFunctionPatch<M>(d,reinterpret_cast<Vector<Expansion<ExactFloat>>const&>(f),Vector<ErrorFloat>(f.size()),swp)
{
}



template<class M> VectorFunctionPatch<M>::VectorFunctionPatch(const ExactBox& d,
                                           const VectorFunctionType<M>& f,
                                           const Sweeper& swp)
    : _domain(d), _models(f.result_size())
{
    //ARIADNE_ASSERT_MSG(f.result_size()>0, "d="<<d<<", f="<<f<<", swp="<<swp);
    ARIADNE_ASSERT(d.size()==f.argument_size());
    Vector<ModelType> x=ModelType::scalings(d,swp);
    this->_models=f.evaluate(x);
    this->sweep();
}

template<class M> VectorFunctionPatch<M>::VectorFunctionPatch(const Vector<FunctionPatch<M>>& v)
    : _domain(), _models(v.size())
{
    ARIADNE_ASSERT(v.size()>0);
    for(SizeType i=0; i!=v.size(); ++i) { ARIADNE_ASSERT(v[i].domain()==v.zero_element().domain()); }
    this->_domain=v.zero_element().domain();
    for(SizeType i=0; i!=v.size(); ++i) {
        this->_models[i]=v[i].model();
    }
}

template<class M> VectorFunctionPatch<M>::VectorFunctionPatch(const List<FunctionPatch<M>>& v)
    : _domain(), _models(v.size())
{
    ARIADNE_ASSERT(v.size()>0);
    for(SizeType i=1; i!=v.size(); ++i) { ARIADNE_ASSERT(v[i].domain()==v[0].domain()); }
    this->_domain=v[0].domain();
    for(SizeType i=0; i!=v.size(); ++i) {
        this->_models[i]=v[i].model();
    }
}

template<class M> VectorFunctionPatch<M>::VectorFunctionPatch(InitializerList<FunctionPatch<M>> lst)
    : _domain(), _models(lst.size())
{
    *this=VectorFunctionPatch<M>(List<FunctionPatch<M>>(lst));
}


template<class M> VectorFunctionPatch<M>* VectorFunctionPatch<M>::_clone() const
{
    return new VectorFunctionPatch<M>(*this);
}

template<class M> VectorFunctionPatch<M>* VectorFunctionPatch<M>::_create() const
{
    return new VectorFunctionPatch<M>(this->result_size(), FunctionPatch<M>(this->domain(),this->sweeper()));
}

template<class M> FunctionPatch<M>* VectorFunctionPatch<M>::_create_zero() const
{
    return new FunctionPatch<M>(this->domain(),this->sweeper());
}

template<class M> VectorFunctionPatch<M>* VectorFunctionPatch<M>::_create_identity() const
{
    Sweeper sweeper=this->sweeper();
    VectorFunctionPatch<M>* result = new VectorFunctionPatch<M>(this->domain().size(), FunctionPatch<M>(this->domain(),sweeper));
    for(SizeType i=0; i!=result->size(); ++i) { (*result)[i]=FunctionPatch<M>::coordinate(this->domain(),i,sweeper); }
    return result;
}

template<class M> Void VectorFunctionPatch<M>::adjoin(const FunctionPatch<M>& sf)
{
    ARIADNE_ASSERT_MSG(sf.domain()==this->domain(),"sf="<<sf);
    this->_models=join(this->_models,sf.model());
}







template<class M> VectorFunctionPatch<M> VectorFunctionPatch<M>::constant(const ExactBox& d, const Vector<NumericType>& c, Sweeper swp)
{
    return VectorFunctionPatch<M>(d,ModelType::constants(d.size(),c,swp));
}

template<class M> VectorFunctionPatch<M> VectorFunctionPatch<M>::identity(const ExactBox& d, Sweeper swp)
{
    return VectorFunctionPatch<M>(d,ModelType::scalings(d,swp));
}

template<class M> VectorFunctionPatch<M> VectorFunctionPatch<M>::projection(const ExactBox& d, SizeType imin, SizeType imax, Sweeper swp)
{
    return VectorFunctionPatch<M>(FunctionPatch<M>::coordinates(d,imin,imax,swp));
}


template<class M> Polynomial<ValidatedFloat> polynomial(const M& tm) {
    return Polynomial<ValidatedFloat>(tm.expansion())+ValidatedNumber(-tm.error(),+tm.error());
}

template<class M> Vector<Polynomial<ValidatedFloat>> VectorFunctionPatch<M>::polynomials() const
{
    Vector<Polynomial<ValidatedFloat> > p(this->result_size(),Polynomial<ValidatedFloat>(this->argument_size()));
    for(SizeType i=0; i!=this->result_size(); ++i) {
        p[i]=static_cast<FunctionPatch<M>>((*this)[i]).polynomial();
    }
    return p;
}

template<class M> Vector<Expansion<ExactFloat>> const VectorFunctionPatch<M>::expansions() const
{
    Vector<Expansion<ExactFloat>> e(this->result_size(),Expansion<ExactFloat>(this->argument_size()));
    for(SizeType i=0; i!=this->result_size(); ++i) {
        e[i]=this->models()[i].expansion();
    }
    return e;
}

template<class M> Vector<typename VectorFunctionPatch<M>::ErrorType> const VectorFunctionPatch<M>::errors() const
{
    Vector<ErrorFloat> e(this->result_size());
    for(SizeType i=0; i!=this->result_size(); ++i) {
        e[i]=this->models()[i].error();
    }
    return e;
}

template<class M> typename VectorFunctionPatch<M>::ErrorType const VectorFunctionPatch<M>::error() const
{
    ErrorFloat e=0;
    for(SizeType i=0; i!=this->result_size(); ++i) {
        e=max(e,this->models()[i].error());
    }
    return e;
}

template<class M> VectorFunctionType<M> VectorFunctionPatch<M>::function() const
{
    return ValidatedVectorFunction(new VectorFunctionPatch<M>(*this));
}

template<class M> Bool VectorFunctionPatch<M>::operator==(const VectorFunctionPatch<M>& tm) const
{
    return this->_models==tm._models;
}



template<class M> Bool VectorFunctionPatch<M>::operator!=(const VectorFunctionPatch<M>& p2) const
{
    return !(*this==p2);
}



template<class M> Sweeper VectorFunctionPatch<M>::sweeper() const
{
    ARIADNE_ASSERT(this->size()>0); return this->_models[0].sweeper();
}


template<class M> Void VectorFunctionPatch<M>::set_sweeper(Sweeper swp)
{
    for(SizeType i=0; i!=this->result_size(); ++i) {
        this->_models[i].set_sweeper(swp);
    }
}

template<class M> const ExactBox& VectorFunctionPatch<M>::domain() const
{
    return this->_domain;
}

template<class M> const ExactBox VectorFunctionPatch<M>::codomain() const
{
    ExactBox result(this->result_size());
    for(SizeType i=0; i!=result.size(); ++i) {
        result[i]=this->_models[i].codomain();
    }
    return result;
}


template<class M> const UpperBox VectorFunctionPatch<M>::range() const
{
    Vector<UpperInterval> result(this->result_size());
    for(SizeType i=0; i!=result.size(); ++i) {
        result[i]=this->_models[i].range();
    }
    return UpperBox(result);
}


template<class M> const Vector<typename VectorFunctionPatch<M>::CoefficientType> VectorFunctionPatch<M>::centre() const
{
    Vector<ExactFloat> result(this->result_size());
    for(SizeType i=0; i!=result.size(); ++i) {
        result[i]=this->_models[i].value();
    }
    return result;
}


template<class M> const Vector<typename VectorFunctionPatch<M>::ModelType>& VectorFunctionPatch<M>::models() const
{
    return this->_models;
}

template<class M> Vector<typename VectorFunctionPatch<M>::ModelType>& VectorFunctionPatch<M>::models()
{
    return this->_models;
}

template<class M> const typename VectorFunctionPatch<M>::ModelType& VectorFunctionPatch<M>::model(SizeType i) const
{
    return this->_models[i];
}

template<class M> typename VectorFunctionPatch<M>::ModelType& VectorFunctionPatch<M>::model(SizeType i)
{
    return this->_models[i];
}




template<class M> SizeType VectorFunctionPatch<M>::argument_size() const
{
    return this->_domain.size();
}


template<class M> SizeType VectorFunctionPatch<M>::result_size() const
{
    return this->_models.size();
}


template<class M> FunctionPatch<M> VectorFunctionPatch<M>::operator[](SizeType i) const
{
    return this->get(i);
}

template<class M> VectorFunctionPatchElementReference<M> VectorFunctionPatch<M>::operator[](SizeType i)
{
    return VectorFunctionPatchElementReference<M>(*this,i);
}

template<class M> FunctionPatch<M> VectorFunctionPatch<M>::get(SizeType i) const
{
    return FunctionPatch<M>(this->_domain,this->_models[i]);
}

template<class M> Void VectorFunctionPatch<M>::set(SizeType i, const FunctionPatch<M>& e)
{
    ARIADNE_ASSERT_MSG(this->size()>i,"Cannot set "<<i<<"th element of VectorFunctionPatch<M> "<<(*this));
    if(this->domain().size()!=0) {
        ARIADNE_ASSERT_MSG(e.domain()==this->domain(),"Domain of "<<e<<" conflicts with existing domain "<<this->domain());
    } else {
        this->_domain=e.domain();
    }
    this->_models[i]=e.model();
}














template<class M> template<class T> Void VectorFunctionPatch<M>::_compute(Vector<T>& r, const Vector<T>& a) const
{
    typedef typename T::NumericType X;
    const VectorFunctionPatch<M>& f=*this;
    Vector<T> sx=Ariadne::unscale(a,f._domain);
    for(SizeType i=0; i!=r.size(); ++i) {
        T ri=Ariadne::evaluate(this->_models[i].expansion(),sx);
        X e=convert_error<X>(this->_models[i].error());
        r[i]=ri+e;
    }
}




template<class M> VectorFunctionPatch<M>& VectorFunctionPatch<M>::sweep()
{
    for(SizeType i=0; i!=this->size(); ++i) {
        this->_models[i].sweep();
    }
    return *this;
}

template<class M> VectorFunctionPatch<M>& VectorFunctionPatch<M>::sweep(const SweeperInterface& sweeper)
{
    for(SizeType i=0; i!=this->size(); ++i) {
        this->_models[i].sweep(sweeper);
    }
    return *this;
}


template<class M> Void VectorFunctionPatch<M>::clobber()
{
    for(SizeType i=0; i!=this->size(); ++i) {
        this->_models[i].clobber();
    }
}





template<class M> Vector<ApproximateNumber> VectorFunctionPatch<M>::operator()(const Vector<ApproximateNumber>& x) const
{
    const VectorFunctionPatch<M>& f=*this;
    if(!decide(contains(f.domain(),x))) {
        ARIADNE_THROW(DomainException,"tf.evaluate(ax) with tf="<<f<<", ax="<<x,"ax is not an element of tf.domain()="<<f.domain());
    }
    Vector<ApproximateNumber> sx=Ariadne::unscale(x,f._domain);
    Vector<ApproximateNumber> r(this->result_size());
    for(SizeType i=0; i!=r.size(); ++i) {
        r[i]=Ariadne::evaluate(this->_models[i].expansion(),sx);
    }
    return r;
}

template<class M> Vector<ValidatedNumber> VectorFunctionPatch<M>::operator()(const Vector<ValidatedNumber>& x) const
{
    const VectorFunctionPatch<M>& f=*this;
    if(!contains(f.domain(),x)) {
        ARIADNE_THROW(DomainException,"tf.evaluate(vx) with tf="<<f<<", x="<<x,"vx is not a subset of tf.domain()="<<f.domain());
    }
    Vector<ValidatedNumber> sx=Ariadne::unscale(x,f._domain);
    return Ariadne::evaluate(f._models,sx);
}

template<class M> Matrix<typename VectorFunctionPatch<M>::NumericType> VectorFunctionPatch<M>::jacobian(const Vector<NumericType>& x) const
{
    Matrix<NumericType> J=Ariadne::jacobian(this->_models,unscale(x,this->_domain));
    for(SizeType j=0; j!=J.column_size(); ++j) {
        NumericType rad=rad_val(this->_domain[j]);
        for(SizeType i=0; i!=J.row_size(); ++i) {
            J[i][j]/=rad;
        }
    }
    return J;
}

template<class M> Void VectorFunctionPatch<M>::restrict(const ExactBox& x)
{
    *this=restriction(*this,x);
}


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
    return VectorFunctionPatch<M>(f1.domain(),join(f1.model(),f2.model()));
}

template<class M> VectorFunctionPatch<M> join(const FunctionPatch<M>& f1, const VectorFunctionPatch<M>& f2)
{
    ARIADNE_ASSERT(f1.domain()==f2.domain());
    return VectorFunctionPatch<M>(f1.domain(),join(f1.model(),f2.models()));
}

template<class M> VectorFunctionPatch<M> combine(const FunctionPatch<M>& f1, const FunctionPatch<M>& f2)
{
    return VectorFunctionPatch<M>(product(f1.domain(),f2.domain()),combine(f1.model(),f2.model()));
}

template<class M> VectorFunctionPatch<M> combine(const FunctionPatch<M>& f1, const VectorFunctionPatch<M>& f2)
{
    return VectorFunctionPatch<M>(product(f1.domain(),f2.domain()),combine(f1.model(),f2.models()));
}

template<class M> VectorFunctionPatch<M> combine(const VectorFunctionPatch<M>& f1, const FunctionPatch<M>& f2)
{
    return VectorFunctionPatch<M>(product(f1.domain(),f2.domain()),combine(f1.models(),f2.model()));
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

template<class M> Bool disjoint(const VectorFunctionPatch<M>& f1, const VectorFunctionPatch<M>& f2) {
    ARIADNE_ASSERT(f1.result_size()==f2.result_size());
    for(SizeType i=0; i!=f1.result_size(); ++i) {
        if(disjoint(f1[i],f2[i])) { return true; }
    }
    return false;
}

template<class M> VectorFunctionPatch<M> intersection(const VectorFunctionPatch<M>& f1, const VectorFunctionPatch<M>& f2) {
    ARIADNE_ASSERT(f1.result_size()==f2.result_size());
    VectorFunctionPatch<M> r(f1.result_size());
    for(SizeType i=0; i!=r.result_size(); ++i) {
        r[i]=intersection(f1[i],f2[i]);
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
    ARIADNE_ASSERT(!empty(intersection(f1.domain(),f2.domain())));
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
    ARIADNE_ASSERT(!empty(intersection(f1.domain(),f2.domain())));
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

template<class M> VectorFunctionPatch<M> operator*(const Matrix<Float>& A, const VectorFunctionPatch<M>& f)
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
    if(!contains(f.domain(),x)) {
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






template<class M> ErrorFloat norm(const VectorFunctionPatch<M>& f) {
    ErrorFloat res=0u;
    for(SizeType i=0; i!=f.result_size(); ++i) {
        res=max(res,norm(f[i]));
    }
    return res;
}

template<class M> ErrorFloat distance(const VectorFunctionPatch<M>& f1, const VectorFunctionPatch<M>& f2) {
    return norm(f1-f2);
}

template<class M> ErrorFloat distance(const VectorFunctionPatch<M>& f1, const ValidatedVectorFunction& f2) {
    return distance(f1,VectorFunctionPatch<M>(f1.domain(),f2,f1.sweeper()));
}


template<class M> OutputStream& VectorFunctionPatch<M>::write(OutputStream& os) const
{
    os << "[";
    for(SizeType i=0; i!=this->result_size(); ++i) {
        if(i!=0) { os << ", "; }
        FunctionPatch<M> tfi=(*this)[i];
        os << tfi;
    }
    os << "]";
    return os;
}

template<class M> OutputStream& VectorFunctionPatch<M>::repr(OutputStream& os) const
{
    return os << "VectorFunctionPatch<M>("
              << representation(this->domain()) << ", " << representation(this->expansions()) << ", "
              << representation(this->errors()) << ", " << representation(this->sweeper()) << ")";
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

template<class M> Polynomial<ValidatedFloat> polynomial(const FunctionPatch<M>& tfn) {
    return Polynomial<ValidatedFloat>(tfn.polynomial());
}

template<class M> Vector< Polynomial<ValidatedFloat> > polynomials(const VectorFunctionPatch<M>& tfn) {
    return tfn.polynomials();
}



} // namespace Ariadne
