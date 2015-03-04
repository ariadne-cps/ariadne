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

template<class M> Void _set_scaling(FunctionPatch<M>& x, const ExactIntervalType& ivl, SizeType j)
{
    Float64::RoundingModeType rounding_mode=Float64::get_rounding_mode();
    Float64::set_rounding_upward();
    const Float64& l=ivl.lower().raw();
    const Float64& u=ivl.upper().raw();
    VOLATILE Float64 pc=u; pc+=l;
    VOLATILE Float64 nc=-u; nc-=l;
    VOLATILE Float64 pg=u; pg-=l;
    VOLATILE Float64 ng=l; ng-=u;
    x.error()=ErrorType((pc+nc+pg+ng)/4);
    Float64::set_rounding_to_nearest();
    MultiIndex a(x.argument_size());
    x.expansion().raw().append(a,(l+u)/2);
    ++a[j];
    x.expansion().raw().append(a,(l+u)/2);
    Float64::set_rounding_mode(rounding_mode);
}



inline OutputStream& operator<<(OutputStream& os, const Representation<Float64>& flt_repr)
{
    const Float64& flt=*flt_repr.pointer;
    Int precision=os.precision(); std::ios_base::fmtflags flags = os.flags();
    os.precision(17); os.setf(std::ios_base::showpoint);
    os << "Float64(" << flt << ")";
    os.precision(precision); os.flags(flags);
    return os;
}

inline OutputStream& operator<<(OutputStream& os, const Representation<PositiveUpperFloat64>& flt_repr)
{
    return os << reinterpret_cast<Representation<Float64>const&>(flt_repr);
}

inline OutputStream& operator<<(OutputStream& os, const Representation<ExactIntervalType>& ivl_repr)
{
    const ExactIntervalType& ivl=*ivl_repr.pointer;
    Int precision=os.precision(); std::ios_base::fmtflags flags = os.flags();
    os.precision(17); os.setf(std::ios_base::showpoint);
    os << "ExactIntervalType("<<ivl.lower()<<","<<ivl.upper()<<")";
    os.precision(precision); os.flags(flags);
    return os;
}


inline OutputStream& operator<<(OutputStream& os, const Representation< Expansion<Float64> >& exp_repr)
{
    const Expansion<Float64>& exp=*exp_repr.pointer;
    Int precision=os.precision(); std::ios_base::fmtflags flags = os.flags();
    os.precision(17); os.setf(std::ios_base::showpoint);
    os << "Expansion<Float64>(" << exp.argument_size() << "," << exp.number_of_nonzeros();
    for(Expansion<Float64>::ConstIterator iter=exp.begin(); iter!=exp.end(); ++iter) {
        for(SizeType j=0; j!=iter->key().size(); ++j) {
            os << "," << Nat(iter->key()[j]);
        }
        os << "," << iter->data();
    }
    os << ")";
    os.precision(precision); os.flags(flags);
    return os;
}

inline OutputStream& operator<<(OutputStream& os, const Representation< Expansion<ExactFloat64> >& exp_repr) {
    return os << reinterpret_cast<Expansion<Float64>const&>(exp_repr);
}

inline OutputStream& operator<<(OutputStream& os, const Representation< Expansion<ApproximateFloat64> >& exp_repr) {
    return os << reinterpret_cast<Expansion<Float64>const&>(exp_repr);
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

inline OutputStream& operator<<(OutputStream& os, const Representation< ExactBoxType >& box_repr)
{
    const Vector<ExactIntervalType>& vec=*box_repr.pointer;
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

template<class M> FunctionPatch<M>::FunctionPatch(const ExactBoxType& d, Sweeper swp)
    : _domain(d), _model(d.size(),swp)
{
}

template<class M> FunctionPatch<M>::FunctionPatch(const ExactBoxType& d, const Expansion<RawFloat64>& p, const RawFloat64& e, const Sweeper& swp)
    : _domain(d), _model(p,e,swp)
{
}

template<class M> FunctionPatch<M>::FunctionPatch(const ExactBoxType& d, const Expansion<ExactFloat64>& p, const ErrorFloat64& e, const Sweeper& swp)
    : _domain(d), _model(p,e,swp)
{
}

template<class M> FunctionPatch<M>::FunctionPatch(const ExactBoxType& d, const ModelType& m)
    : _domain(d), _model(m)
{
}

template<class M> FunctionPatch<M>::FunctionPatch(const ExactBoxType& d, const ScalarFunctionType<M>& f, Sweeper swp)
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



template<class M> FunctionPatch<M> FunctionPatch<M>::zero(const ExactBoxType& d, Sweeper swp)
{
    return FunctionPatch<M>(d,ModelType::zero(d.size(),swp));
}

template<class M> FunctionPatch<M> FunctionPatch<M>::constant(const ExactBoxType& d, const NumericType& c, Sweeper swp)
{
    return FunctionPatch<M>(d,ModelType::constant(d.size(),c,swp));
}

template<class M> FunctionPatch<M> FunctionPatch<M>::coordinate(const ExactBoxType& d, SizeType j, Sweeper swp)
{
    ARIADNE_ASSERT(j<d.size());
    return FunctionPatch<M>(d,ModelType::scaling(d.size(),j,d[j],swp));
}

template<class M> VectorFunctionPatch<M> FunctionPatch<M>::identity(const ExactBoxType& d, Sweeper swp)
{
    return VectorFunctionPatch<M>(d,ModelType::scalings(d,swp));
}


template<class M> Vector<FunctionPatch<M>> FunctionPatch<M>::constants(const ExactBoxType& d, const Vector<NumericType>& c, Sweeper swp)
{
    ARIADNE_DEPRECATED("FunctionPatch<M>::constants","Use VectorFunctionPatch<M>::constant instead");
    Vector<FunctionPatch<M>> x(c.size(),FunctionPatch<M>(d,swp));
    for(SizeType i=0; i!=c.size(); ++i) {
        x[i]=c[i];
    }
    return x;
}

template<class M> Vector<FunctionPatch<M>> FunctionPatch<M>::coordinates(const ExactBoxType& d, Sweeper swp)
{
    ARIADNE_DEPRECATED("FunctionPatch<M>::coordinates","Use VectorFunctionPatch<M>::identity instead");
    Vector<FunctionPatch<M>> x(d.dimension(),FunctionPatch<M>(d,swp));
    for(SizeType i=0; i!=x.size(); ++i) {
        x[i]=FunctionPatch<M>::coordinate(d,i,swp);
    }
    return x;
}

template<class M> Vector<FunctionPatch<M>> FunctionPatch<M>::coordinates(const ExactBoxType& d, SizeType imin, SizeType imax, Sweeper swp)
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

template<class M> FunctionPatch<M> FunctionPatch<M>::create_constant(NumericType const& c) const
{
    return FunctionPatch<M>::constant(this->domain(),c,this->_model.sweeper());
}

template<class M> FunctionPatch<M> FunctionPatch<M>::create(GenericType const& f) const
{
    return FunctionPatch<M>(this->domain(),f,this->_model.sweeper());
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

template<class M> FunctionPatch<M>* FunctionPatch<M>::_create_zero(DomainType const& dom) const
{
    return new FunctionPatch<M>(dom,this->sweeper());
}

template<class M> FunctionPatch<M>* FunctionPatch<M>::_create_constant(DomainType const& dom, NumericType const& c) const
{
    return new FunctionPatch<M>(FunctionPatch<M>::constant(dom,c,this->sweeper()));
}

template<class M> FunctionPatch<M>* FunctionPatch<M>::_create_coordinate(DomainType const& dom, SizeType j) const
{
    return new FunctionPatch<M>(FunctionPatch<M>::coordinate(dom,j,this->sweeper()));
}

template<class M> Void VectorFunctionPatch<M>::adjoin(const FunctionPatch<M>& sf)
{
    ARIADNE_ASSERT_MSG(sf.domain()==this->domain(),"sf="<<sf);
    this->_models=join(this->_models,sf.model());
}


template<class M> Void FunctionPatch<M>::restrict(const ExactBoxType& dom) {
    (*this)=restriction(*this,dom);
}



inline Bool operator==(ExactFloat64 x1, Int n2) { return x1.raw()==Float64(n2); }
inline Bool operator==(BoundedFloat64 x1, Int n2) { return x1.upper_raw()==Float64(n2) && x1.lower_raw()==Float64(n2); }
inline Bool operator==(ApproximateFloat64 x1, Int n2) { return x1.raw()==Float64(n2); }

inline Bool operator!=(ExactFloat64 x1, Int n2) { return x1.raw()!=Float64(n2); }
inline Bool operator!=(BoundedFloat64 x1, Int n2) { return x1.upper_raw()!=Float64(n2) || x1.lower_raw()!=Float64(n2); }
inline Bool operator!=(ApproximateFloat64 x1, Int n2) { return x1.raw()!=Float64(n2); }

inline Bool operator> (ExactFloat64 x1, Int n2) { return x1.raw()> Float64(n2); }
inline Bool operator> (BoundedFloat64 x1, Int n2) { return x1.lower_raw()> Float64(n2); }
inline Bool operator> (ApproximateFloat64 x1, Int n2) { return x1.raw()> Float64(n2); }

template<class M> Polynomial<BoundedFloat64> FunctionPatch<M>::polynomial() const
{
    Vector<Polynomial<BoundedFloat64> > pid=Polynomial<NumericType>::coordinates(this->argument_size());
    return horner_evaluate(this->expansion(),unscale(pid,this->domain()))+BoundedFloat64(-this->error(),+this->error());

    Polynomial<BoundedFloat64> z(this->argument_size());
    Polynomial<BoundedFloat64> p;//=Ariadne::polynomial(this->model());

    Vector<Polynomial<BoundedFloat64> > s(this->argument_size(),z);
    for(SizeType j=0; j!=this->argument_size(); ++j) {
        ExactIntervalType const& domj=this->domain()[j];
        if(domj.lower()>=domj.upper()) {
            ARIADNE_ASSERT(this->domain()[j].width()==0);
            s[j]=Polynomial<BoundedFloat64>::constant(this->argument_size(),0);
        } else {
            //s[j]=Ariadne::polynomial(ModelType::unscaling(this->argument_size(),j,this->domain()[j],this->sweeper()));
            s[j]=(Polynomial<BoundedFloat64>::coordinate(this->argument_size(),j)-domj.midpoint())/domj.radius();
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



template<class M> FunctionPatch<M>* FunctionPatch<M>::_derivative(SizeType j) const
{
    return new FunctionPatch<M>(Ariadne::derivative(*this,j));
}


template<class M> ApproximateNumericType FunctionPatch<M>::operator() (const Vector<ApproximateNumericType>& x) const
{
    const FunctionPatch<M>& f=*this;
    if(!contains(f.domain(),cast_exact(x))) {
        ARIADNE_THROW(DomainException,"tf.evaluate(ax) with tf="<<f<<", ax="<<x," ax is not an element of tf.domain()="<<f.domain());
    }
    Vector<ApproximateNumericType> sx=Ariadne::unscale(x,f._domain);
    return Ariadne::evaluate(this->_model.expansion(),sx);
}

template<class M> ValidatedNumericType FunctionPatch<M>::operator()(const Vector<ValidatedNumericType>& x) const
{
    const FunctionPatch<M>& f=*this;
    if(!definitely(contains(f.domain(),x))) {
        ARIADNE_THROW(DomainException,"evaluate(f,x) with f="<<f<<", x="<<x,"x is not definitely and element of f.domain()="<<f.domain());
    }
    return unchecked_evaluate(f,x);
}

template<class M> ValidatedNumericType FunctionPatch<M>::operator()(const Vector<ExactNumericType>& x) const
{
    return Ariadne::evaluate(*this,Vector<ValidatedNumericType>(x));
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



template<class M> OutputStream& FunctionPatch<M>::write(OutputStream& os) const {
    Polynomial<BoundedFloat64> p=this->polynomial();
    Polynomial<ApproximateFloat64> ap=p;
    os << "FP" << this->domain();
    os << "(";
    os << ap;
    if(this->error().raw()>0.0) { os << "+/-" << this->error(); }
    os << ")";
    return os;
}

template<class M> OutputStream& FunctionPatch<M>::repr(OutputStream& os) const
{
    return os << "FunctionPatch<M>(" << representation(this->domain()) << ", " << representation(this->model().expansion().raw())
              << "," << representation(this->error().raw())<<","<<this->sweeper()<<")";
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
    Float64 truncatation_error = 0.0;
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
    ErrorFloat64 truncatation_error = truncated_function.error();
    truncated_function.clobber();
    Polynomial<BoundedFloat64> validated_polynomial_function=polynomial(truncated_function);
    Polynomial<ExactFloat64> polynomial_function = midpoint(validated_polynomial_function);
    if(frepr.names.empty()) { os << polynomial_function; }
    else { os << named_argument_repr(polynomial_function,frepr.names); }
    os << "+/-" << truncatation_error << "+/-" << function.error();
    return os;
}








template<class M> VectorFunctionPatch<M>::VectorFunctionPatch()
    : _domain(), _models()
{
}

template<class M> VectorFunctionPatch<M>::VectorFunctionPatch(SizeType k)
    : _domain(), _models(k)
{
}

template<class M> VectorFunctionPatch<M>::VectorFunctionPatch(SizeType m, const ExactBoxType& d, Sweeper swp)
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




template<class M> VectorFunctionPatch<M>::VectorFunctionPatch(const ExactBoxType& d,
                                           const Vector<ModelType>& f)
    : _domain(d), _models(f)
{
    for(SizeType i=0; i!=f.size(); ++i) {
        ARIADNE_ASSERT_MSG(d.size()==f[i].argument_size(),"d="<<d<<", f="<<f);
    }
}

template<class M> VectorFunctionPatch<M>::VectorFunctionPatch(const ExactBoxType& d,
                                           const Vector<Expansion<ExactFloat64>>& f,
                                           const Vector<ErrorFloat64>& e,
                                           Sweeper swp)
    : _domain(d), _models(f.size(),ModelType(d.size(),swp))
{
    ARIADNE_ASSERT(f.size()==e.size());
    for(SizeType i=0; i!=f.size(); ++i) {
        ARIADNE_ASSERT(d.size()==f[i].argument_size());
        _models[i]=ModelType(f[i],e[i],swp);
    }
}

template<class M> VectorFunctionPatch<M>::VectorFunctionPatch(const ExactBoxType& d,
                                           const Vector<Expansion<ExactFloat64>>& f,
                                           Sweeper swp)
    : VectorFunctionPatch<M>(d,f,Vector<ErrorFloat64>(f.size()),swp)
{
}

template<class M> VectorFunctionPatch<M>::VectorFunctionPatch(const ExactBoxType& d,
                                           const Vector<Expansion<RawFloat64>>& f,
                                           const Vector<RawFloat64>& e,
                                           Sweeper swp)
    : VectorFunctionPatch<M>(d,reinterpret_cast<Vector<Expansion<ExactFloat64>>const&>(f),
                           reinterpret_cast<Vector<ErrorFloat64>const&>(e),swp)
{
}

template<class M> VectorFunctionPatch<M>::VectorFunctionPatch(const ExactBoxType& d,
                                           const Vector<Expansion<RawFloat64>>& f,
                                           Sweeper swp)
    : VectorFunctionPatch<M>(d,reinterpret_cast<Vector<Expansion<ExactFloat64>>const&>(f),Vector<ErrorFloat64>(f.size()),swp)
{
}



template<class M> VectorFunctionPatch<M>::VectorFunctionPatch(const ExactBoxType& d,
                                           const VectorFunctionType<M>& f,
                                           const Sweeper& swp)
    : _domain(d), _models(f.result_size())
{
    //ARIADNE_ASSERT_MSG(f.result_size()>0, "d="<<d<<", f="<<f<<", swp="<<swp);
    ARIADNE_ASSERT(d.size()==f.argument_size());
    Vector<ModelType> x=ModelType::scalings(d,swp);
    this->_models=f(x);
    ARIADNE_DEBUG_ASSERT(this->argument_size()==f.argument_size());
    ARIADNE_DEBUG_ASSERT_MSG(this->result_size()==f.result_size(),"  f="<<f<<"\n  r="<<*this<<"\n");
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







template<class M> VectorFunctionPatch<M> VectorFunctionPatch<M>::constant(const ExactBoxType& d, const Vector<NumericType>& c, Sweeper swp)
{
    return VectorFunctionPatch<M>(d,ModelType::constants(d.size(),c,swp));
}

template<class M> VectorFunctionPatch<M> VectorFunctionPatch<M>::identity(const ExactBoxType& d, Sweeper swp)
{
    return VectorFunctionPatch<M>(d,ModelType::scalings(d,swp));
}

template<class M> VectorFunctionPatch<M> VectorFunctionPatch<M>::projection(const ExactBoxType& d, SizeType imin, SizeType imax, Sweeper swp)
{
    return VectorFunctionPatch<M>(FunctionPatch<M>::coordinates(d,imin,imax,swp));
}


template<class M> Vector<Polynomial<BoundedFloat64>> VectorFunctionPatch<M>::polynomials() const
{
    Vector<Polynomial<BoundedFloat64> > p(this->result_size(),Polynomial<BoundedFloat64>(this->argument_size()));
    for(SizeType i=0; i!=this->result_size(); ++i) {
        p[i]=static_cast<FunctionPatch<M>>((*this)[i]).polynomial();
    }
    return p;
}

template<class M> Vector<Expansion<ExactFloat64>> const VectorFunctionPatch<M>::expansions() const
{
    Vector<Expansion<ExactFloat64>> e(this->result_size(),Expansion<ExactFloat64>(this->argument_size()));
    for(SizeType i=0; i!=this->result_size(); ++i) {
        e[i]=this->models()[i].expansion();
    }
    return e;
}

template<class M> Vector<typename VectorFunctionPatch<M>::ErrorType> const VectorFunctionPatch<M>::errors() const
{
    Vector<ErrorFloat64> e(this->result_size());
    for(SizeType i=0; i!=this->result_size(); ++i) {
        e[i]=this->models()[i].error();
    }
    return e;
}

template<class M> typename VectorFunctionPatch<M>::ErrorType const VectorFunctionPatch<M>::error() const
{
    ErrorFloat64 e=0u;
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

template<class M> const ExactBoxType VectorFunctionPatch<M>::domain() const
{
    return this->_domain;
}

template<class M> const ExactBoxType VectorFunctionPatch<M>::codomain() const
{
    ExactBoxType result(this->result_size());
    for(SizeType i=0; i!=result.size(); ++i) {
        result[i]=this->_models[i].codomain();
    }
    return result;
}


template<class M> const UpperBoxType VectorFunctionPatch<M>::range() const
{
    Vector<UpperIntervalType> result(this->result_size());
    for(SizeType i=0; i!=result.size(); ++i) {
        result[i]=this->_models[i].range();
    }
    return UpperBoxType(result);
}


template<class M> const Vector<typename VectorFunctionPatch<M>::CoefficientType> VectorFunctionPatch<M>::centre() const
{
    Vector<ExactFloat64> result(this->result_size());
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


template<class M> FunctionPatch<M> const VectorFunctionPatch<M>::operator[](SizeType i) const
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
    ARIADNE_DEBUG_ASSERT_MSG(r.size()==f.result_size(),"\nr="<<r<<"\nf="<<f<<"\n");
    Vector<T> sx=Ariadne::unscale(a,f._domain);
    for(SizeType i=0; i!=r.size(); ++i) {
        T ri=Ariadne::evaluate(this->_models[i].expansion(),sx);
        X e=convert_error_to_bounds(this->_models[i].error());
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





template<class M> Vector<ApproximateNumericType> VectorFunctionPatch<M>::operator()(const Vector<ApproximateNumericType>& x) const
{
    const VectorFunctionPatch<M>& f=*this;
    if(!decide(contains(f.domain(),x))) {
        ARIADNE_THROW(DomainException,"tf.evaluate(ax) with tf="<<f<<", ax="<<x,"ax is not an element of tf.domain()="<<f.domain());
    }
    Vector<ApproximateNumericType> sx=Ariadne::unscale(x,f._domain);
    Vector<ApproximateNumericType> r(this->result_size());
    for(SizeType i=0; i!=r.size(); ++i) {
        r[i]=Ariadne::evaluate(this->_models[i].expansion(),sx);
    }
    return r;
}

template<class M> Vector<ValidatedNumericType> VectorFunctionPatch<M>::operator()(const Vector<ValidatedNumericType>& x) const
{
    const VectorFunctionPatch<M>& f=*this;
    if(!definitely(contains(f.domain(),x))) {
        ARIADNE_THROW(DomainException,"tf.evaluate(vx) with tf="<<f<<", x="<<x,"vx is not a definitely and element of tf.domain()="<<f.domain());
    }
    Vector<ValidatedNumericType> sx=Ariadne::unscale(x,f._domain);
    return Ariadne::evaluate(f._models,sx);
}

template<class M> Matrix<typename VectorFunctionPatch<M>::NumericType> VectorFunctionPatch<M>::jacobian(const Vector<NumericType>& x) const
{
    Vector<NumericType> y=unscale(x,this->_domain);
    Matrix<NumericType> J(this->size(),x.size());
    for(SizeType i=0; i!=J.row_size(); ++i) {
        J[i]=gradient(this->_models[i],y);
    }
    for(SizeType j=0; j!=J.column_size(); ++j) {
        NumericType rad=rad_val(this->_domain[j]);
        for(SizeType i=0; i!=J.row_size(); ++i) {
            J[i][j]/=rad;
        }
    }
    return J;
}

template<class M> Void VectorFunctionPatch<M>::restrict(const ExactBoxType& x)
{
    *this=restriction(*this,x);
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



} // namespace Ariadne
