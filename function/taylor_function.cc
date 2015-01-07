/***************************************************************************
 *            taylor_function.cc
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
#include "function/taylor_model.h"

#include "function/function.h"
#include "function/taylor_function.h"

#include "function_mixin.tcc"
#include "function/taylor_function.h"

#define VOLATILE ;

namespace Ariadne {

typedef unsigned int uint;

static double TAYLOR_FUNCTION_WRITING_ACCURACY = 1e-8;

void _set_scaling(ScalarTaylorFunction& x, const ExactInterval& ivl, uint j)
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



ScalarFunctionModel<ValidatedTag>& ScalarFunctionModel<ValidatedTag>::operator=(const ScalarTaylorFunction& f) {
    this->_ptr=clone_on_copy_ptr< ScalarFunctionModelInterface<ValidatedTag> >(new ScalarTaylorFunction(f)); return *this;
}


ScalarTaylorFunction::ScalarTaylorFunction()
    : _domain(), _model()
{ }

ScalarTaylorFunction::ScalarTaylorFunction(const ExactBox& d, Sweeper swp)
    : _domain(d), _model(d.size(),swp)
{
}

ScalarTaylorFunction::ScalarTaylorFunction(const ExactBox& d, const Expansion<Float>& p, const double& e, const Sweeper& swp)
    : _domain(d), _model(p,e,swp)
{
}

ScalarTaylorFunction::ScalarTaylorFunction(const ExactBox& d, const Expansion<ExactFloat>& p, const ErrorFloat& e, const Sweeper& swp)
    : _domain(d), _model(p,e,swp)
{
}

ScalarTaylorFunction::ScalarTaylorFunction(const ExactBox& d, const ValidatedTaylorModel& m)
    : _domain(d), _model(m)
{
}

ScalarTaylorFunction::ScalarTaylorFunction(const ExactBox& d, const ValidatedScalarFunction& f, Sweeper swp)
    : _domain(d), _model(f.argument_size(),swp)
{
    ARIADNE_ASSERT_MSG(d.size()==f.argument_size(),"d="<<d<<" f="<<f);
    Vector<ValidatedTaylorModel> x=ValidatedTaylorModel::scalings(d,swp);
    this->_model=f.evaluate(x);
    this->_model.sweep();
}

ScalarTaylorFunction::ScalarTaylorFunction(const ValidatedScalarFunctionModel& f) {
     ARIADNE_ASSERT_MSG(dynamic_cast<const ScalarTaylorFunction*>(f._ptr.operator->())," f="<<f);
     *this=dynamic_cast<const ScalarTaylorFunction&>(*f._ptr);
}

ScalarTaylorFunction& ScalarTaylorFunction::operator=(const ValidatedScalarFunctionModel& f)
{
    return (*this)=ScalarTaylorFunction(this->domain(),f,this->sweeper());
}



ScalarTaylorFunction ScalarTaylorFunction::zero(const ExactBox& d, Sweeper swp)
{
    return ScalarTaylorFunction(d,ValidatedTaylorModel::zero(d.size(),swp));
}

ScalarTaylorFunction ScalarTaylorFunction::constant(const ExactBox& d, const ValidatedNumber& c, Sweeper swp)
{
    return ScalarTaylorFunction(d,ValidatedTaylorModel::constant(d.size(),c,swp));
}

ScalarTaylorFunction ScalarTaylorFunction::identity(const ExactInterval& ivl, Sweeper swp)
{
    Vector<ExactInterval> d(1,ivl);
    return ScalarTaylorFunction(d,ValidatedTaylorModel::scaling(d.size(),0u,d[0u],swp));
}

ScalarTaylorFunction ScalarTaylorFunction::coordinate(const ExactBox& d, uint j, Sweeper swp)
{
    ARIADNE_ASSERT(j<d.size());
    return ScalarTaylorFunction(d,ValidatedTaylorModel::scaling(d.size(),j,d[j],swp));
}

ScalarTaylorFunction ScalarTaylorFunction::variable(const ExactBox& d, uint j, Sweeper swp)
{
    ARIADNE_DEPRECATED("ScalarTaylorFunction::variable","Use ScalarTaylorFunction::coordinate instead");
    ARIADNE_ASSERT(j<d.size());
    return ScalarTaylorFunction(d,ValidatedTaylorModel::scaling(d.size(),j,d[j],swp));
}


Vector<ScalarTaylorFunction> ScalarTaylorFunction::constants(const ExactBox& d, const Vector<ValidatedNumber>& c, Sweeper swp)
{
    ARIADNE_DEPRECATED("ScalarTaylorFunction::constants","Use VectorTaylorFunction::constant instead");
    Vector<ScalarTaylorFunction> x(c.size(),ScalarTaylorFunction(d,swp));
    for(uint i=0; i!=c.size(); ++i) {
        x[i]=c[i];
    }
    return x;
}

Vector<ScalarTaylorFunction> ScalarTaylorFunction::variables(const ExactBox& d, Sweeper swp)
{
    ARIADNE_DEPRECATED("ScalarTaylorFunction::variables","Use VectorTaylorFunction::identity instead");
    return variables(d,0u,d.size(),swp);
}

Vector<ScalarTaylorFunction> ScalarTaylorFunction::variables(const ExactBox& d, uint imin, uint imax, Sweeper swp)
{
    ARIADNE_DEPRECATED("ScalarTaylorFunction::variables","Use VectorTaylorFunction::projection instead");
    ARIADNE_ASSERT(imin<=imax);
    ARIADNE_ASSERT(imax<=d.size());

    Vector<ScalarTaylorFunction> x(imax-imin);
    for(uint i=imin; i!=imax; ++i) {
        x[i-imin]=ScalarTaylorFunction::variable(d,i,swp);
    }
    return x;
}


ScalarTaylorFunction ScalarTaylorFunction::create_zero() const
{
    return ScalarTaylorFunction(this->domain(),this->_model.sweeper());
}

ScalarTaylorFunction* ScalarTaylorFunction::_clone() const
{
    return new ScalarTaylorFunction(*this);
}

ScalarTaylorFunction* ScalarTaylorFunction::_create() const
{
    return new ScalarTaylorFunction(this->domain(),this->_model.sweeper());
}

VectorFunctionModelInterface<ValidatedTag>* ScalarTaylorFunction::_create_identity() const
{
    Sweeper sweeper=this->sweeper();
    VectorTaylorFunction* result = new VectorTaylorFunction(this->domain().size(), ScalarTaylorFunction(this->domain(),sweeper));
    for(uint i=0; i!=result->size(); ++i) { (*result)[i]=ScalarTaylorFunction::coordinate(this->domain(),i,sweeper); }
    return result;
}

VectorFunctionModelInterface<ValidatedTag>* ScalarTaylorFunction::_create_vector(Nat i) const
{
    return new VectorTaylorFunction(i,this->domain(),this->_model.sweeper());
}



Sweeper ScalarTaylorFunction::sweeper() const
{
    return this->_model.sweeper();
}


Polynomial<ValidatedNumber> polynomial(const ValidatedTaylorModel& tm);

inline bool operator==(ExactFloat x1, int n2) { return x1.raw()==Float(n2); }
inline bool operator==(ValidatedFloat x1, int n2) { return x1.upper_raw()==Float(n2) && x1.lower_raw()==Float(n2); }
inline bool operator==(ApproximateFloat x1, int n2) { return x1.raw()==Float(n2); }

inline bool operator!=(ExactFloat x1, int n2) { return x1.raw()!=Float(n2); }
inline bool operator!=(ValidatedFloat x1, int n2) { return x1.upper_raw()!=Float(n2) || x1.lower_raw()!=Float(n2); }
inline bool operator!=(ApproximateFloat x1, int n2) { return x1.raw()!=Float(n2); }

inline bool operator> (ExactFloat x1, int n2) { return x1.raw()> Float(n2); }
inline bool operator> (ValidatedFloat x1, int n2) { return x1.lower_raw()> Float(n2); }
inline bool operator> (ApproximateFloat x1, int n2) { return x1.raw()> Float(n2); }

Polynomial<ValidatedNumber>
ScalarTaylorFunction::polynomial() const
{
    Polynomial<ValidatedNumber> z(this->argument_size());

    Polynomial<ValidatedNumber> p=Ariadne::polynomial(this->model());

    Vector<Polynomial<ValidatedNumber> > s(this->argument_size(),z);
    for(uint j=0; j!=this->argument_size(); ++j) {
        ExactInterval const& domj=this->domain()[j];
        if(domj.width()<=0) {
            ARIADNE_ASSERT(this->domain()[j].width()==0);
            s[j]=Polynomial<ValidatedNumber>::constant(this->argument_size(),0);
        } else {
            //s[j]=Ariadne::polynomial(ValidatedTaylorModel::unscaling(this->argument_size(),j,this->domain()[j],this->sweeper()));
            s[j]=(Polynomial<ValidatedNumber>::coordinate(this->argument_size(),j)-domj.midpoint())/domj.radius();
        }
    }

    return compose(p,s);
}

ValidatedScalarFunction
ScalarTaylorFunction::function() const
{
    return ValidatedScalarFunction(new ScalarTaylorFunction(*this));
}


bool ScalarTaylorFunction::operator==(const ScalarTaylorFunction& tv) const
{
    return this->_domain==tv._domain && this->_model==tv._model;
}



ScalarTaylorFunction& operator+=(ScalarTaylorFunction& x, const ScalarTaylorFunction& y) {
    ARIADNE_ASSERT(subset(x.domain(),y.domain()));
    if(x.domain()==y.domain()) { x._model+=y._model; }
    else { x._model+=restrict(y,x.domain())._model; }
    return x;
}

ScalarTaylorFunction& operator-=(ScalarTaylorFunction& x, const ScalarTaylorFunction& y) {
    ARIADNE_ASSERT(subset(x.domain(),y.domain()));
    if(x.domain()==y.domain()) { x._model-=y._model; }
    else { x._model-=restrict(y,x.domain())._model; }
    return x;
}


ScalarTaylorFunction& operator+=(ScalarTaylorFunction& f, const ValidatedNumber& c) {
    f._model+=c;
    return f;
}

ScalarTaylorFunction& operator-=(ScalarTaylorFunction& f, const ValidatedNumber& c) {
    f._model-=c;
    return f;
}

ScalarTaylorFunction& operator*=(ScalarTaylorFunction& f, const ValidatedNumber& c) {
    f._model*=c;
    return f;
}

ScalarTaylorFunction& operator/=(ScalarTaylorFunction& f, const ValidatedNumber& c) {
    f._model/=c;
    return f;
}

ScalarTaylorFunction operator+(const ScalarTaylorFunction& x) {
    return ScalarTaylorFunction(x._domain,x._model);
}

ScalarTaylorFunction operator-(const ScalarTaylorFunction& x) {
    return ScalarTaylorFunction(x._domain,-x._model);
}


ScalarTaylorFunction operator+(const ScalarTaylorFunction& x1, const ScalarTaylorFunction& x2) {
    if(x1._domain==x2._domain) {
        return ScalarTaylorFunction(x1._domain,x1._model+x2._model); }
    else {
        ExactBox domain=intersection(x1._domain,x2._domain);
        return ScalarTaylorFunction(domain,restrict(x1,domain)._model+restrict(x2,domain)._model);}
}

ScalarTaylorFunction operator-(const ScalarTaylorFunction& x1, const ScalarTaylorFunction& x2) {
    if(x1._domain==x2._domain) {
        return ScalarTaylorFunction(x1._domain,x1._model-x2._model); }
    else {
        ExactBox domain=intersection(x1._domain,x2._domain);
        return ScalarTaylorFunction(domain,restrict(x1,domain)._model-restrict(x2,domain)._model);}
}

ScalarTaylorFunction operator*(const ScalarTaylorFunction& x1, const ScalarTaylorFunction& x2) {
    if(x1._domain==x2._domain) {
        return ScalarTaylorFunction(x1._domain,x1._model*x2._model); }
    else {
        ExactBox domain=intersection(x1._domain,x2._domain);
        return ScalarTaylorFunction(domain,restrict(x1,domain)._model*restrict(x2,domain)._model);}
}

ScalarTaylorFunction operator/(const ScalarTaylorFunction& x1, const ScalarTaylorFunction& x2) {
    if(x1._domain==x2._domain) {
        return ScalarTaylorFunction(x1._domain,x1._model/x2._model); }
    else {
        ExactBox domain=intersection(x1._domain,x2._domain);
        return ScalarTaylorFunction(domain,restrict(x1,domain)._model/restrict(x2,domain)._model);}
}





ScalarTaylorFunction max(const ScalarTaylorFunction& x1, const ScalarTaylorFunction& x2) {
    if(x1._domain==x2._domain) {
        return ScalarTaylorFunction(x1._domain,max(x1._model,x2._model)); }
    else {
        ExactBox domain=intersection(x1._domain,x2._domain);
        return ScalarTaylorFunction(domain,max(restrict(x1,domain)._model,restrict(x2,domain)._model));}
}

ScalarTaylorFunction min(const ScalarTaylorFunction& x1, const ScalarTaylorFunction& x2) {
    if(x1._domain==x2._domain) {
        return ScalarTaylorFunction(x1._domain,min(x1._model,x2._model)); }
    else {
        ExactBox domain=intersection(x1._domain,x2._domain);
        return ScalarTaylorFunction(domain,min(restrict(x1,domain)._model,restrict(x2,domain)._model));}
}



ScalarTaylorFunction* ScalarTaylorFunction::_derivative(uint j) const
{
    return new ScalarTaylorFunction(Ariadne::derivative(*this,j));
}


ApproximateNumber
ScalarTaylorFunction::evaluate(const Vector<ApproximateNumber>& x) const
{
    const ScalarTaylorFunction& f=*this;
    if(!contains(f.domain(),make_exact(x))) {
        ARIADNE_THROW(DomainException,"tf.evaluate(ax) with tf="<<f<<", ax="<<x," ax is not an element of tf.domain()="<<f.domain());
    }
    Vector<ApproximateNumber> sx=Ariadne::unscale(x,f._domain);
    return Ariadne::evaluate(this->_model.expansion(),sx);
}

ValidatedNumber
ScalarTaylorFunction::evaluate(const Vector<ValidatedNumber>& x) const
{
    return Ariadne::evaluate(*this,x);
}

ValidatedNumber
ScalarTaylorFunction::operator()(const Vector<ValidatedNumber>& x) const
{
    return Ariadne::evaluate(*this,x);
}


void ScalarTaylorFunction::restrict(const ExactBox& d) {
    ARIADNE_ASSERT(subset(d,this->domain()));
    const ExactBox& od=this->domain();
    for(uint j=0; j!=d.size(); ++j) {
        if(od[j]!=d[j]) { *this=Ariadne::restrict(*this,j,d[j]); }
    }
}

ScalarTaylorFunction restrict(const ScalarTaylorFunction& tv, const ExactBox& d) {
    ARIADNE_ASSERT(subset(d,tv.domain()));
    const ExactBox& od=tv.domain();
    ScalarTaylorFunction r=tv;
    for(uint j=0; j!=d.size(); ++j) {
        if(od[j]!=d[j]) { r=restrict(r,j,d[j]); }
    }
    return r;
}

ScalarTaylorFunction extend(const ScalarTaylorFunction& tv, const ExactBox& d) {
    const ExactBox& domain=tv.domain();
    ARIADNE_ASSERT(subset(domain,d));
    for(uint i=0; i!=d.size(); ++i) {
        ARIADNE_ASSERT(domain[i]==d[i] || domain[i].lower()==domain[i].upper());
    }
    return ScalarTaylorFunction(d,tv._model);
}

ValidatedNumber
evaluate(const ScalarTaylorFunction& f, const Vector<ValidatedNumber>& x) {
    if(!contains(f.domain(),x)) {
        ARIADNE_THROW(DomainException,std::setprecision(17)<<"evaluate(tf,vx) with tf="<<f<<", vx="<<x," vx is not a subset of tf.domain()="<<f.domain());
    }
    return unchecked_evaluate(f,x);
}

ValidatedNumber
unchecked_evaluate(const ScalarTaylorFunction& f, const Vector<ValidatedNumber>& x) {
    return evaluate(f.model(),unscale(x,f.domain()));
}


ScalarTaylorFunction
compose(const EffectiveScalarFunction& g, const VectorTaylorFunction& f)
{
    return ScalarTaylorFunction(f.domain(),g.evaluate(f.models()));
}

ScalarTaylorFunction
compose(const ValidatedScalarFunction& g, const VectorTaylorFunction& f)
{
    return ScalarTaylorFunction(f.domain(),g.evaluate(f.models()));
}

ScalarTaylorFunction
compose(const ValidatedScalarFunctionInterface& g, const VectorTaylorFunction& f)
{
    return ScalarTaylorFunction(f.domain(),g.evaluate(f.models()));
}

ScalarTaylorFunction
compose(const ScalarTaylorFunction& g, const VectorTaylorFunction& f)
{
    if(!subset(f.codomain(),g.domain())) {
        ARIADNE_THROW(DomainException,"compose(g,f) with g="<<g<<", f="<<f,"f.codomain()="<<f.codomain()<<" is not a subset of g.domain()="<<g.domain());
    }
    return unchecked_compose(g,f);
}

ScalarTaylorFunction
unchecked_compose(const ScalarTaylorFunction& g, const VectorTaylorFunction& f)
{
    return ScalarTaylorFunction(f.domain(),compose(g.model(),unscale(f.models(),g.domain())));
}



ScalarTaylorFunction
partial_evaluate(const ScalarTaylorFunction& tf, uint k, const ExactNumber& c)
{
    return partial_evaluate(tf,k,ValidatedNumber(c));
}

ScalarTaylorFunction
partial_evaluate(const ScalarTaylorFunction& te, uint k, const ValidatedNumber& c)
{
    // Scale c to domain
    const uint as=te.argument_size();
    ARIADNE_ASSERT(k<as);
    const ExactBox& domain=te.domain();
    const ExactInterval& dk=domain[k];
    ValidatedNumber sc=(c-med_val(dk))/rad_val(dk);

    ExactBox new_domain(as-1);
    for(uint i=0; i!=k; ++i) { new_domain[i]=domain[i]; }
    for(uint i=k; i!=as-1; ++i) { new_domain[i]=domain[i+1]; }

    ValidatedTaylorModel new_model=partial_evaluate(te.model(),k,sc);

    return ScalarTaylorFunction(new_domain,new_model);
}


// To scale from a model on [a,b] to a model on [c,d], use scale factor s=(d-c)/(b-a)
// and translation t=((c+d)-(a+b))/(b-a)
// Because we are scaling the model on [-1,+1], this is not the same as
// the mapping taking [a,b] to [c,d]
ScalarTaylorFunction partial_restrict(const ScalarTaylorFunction& tv, uint k, const ExactInterval& new_ivl) {
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
    return ScalarTaylorFunction(new_dom,preaffine(tv.model(),k,s,t));
}


ScalarTaylorFunction restrict(const ScalarTaylorFunction& tv, uint k, const ExactInterval& new_ivl) {
    return partial_restrict(tv,k,new_ivl);
}

ScalarTaylorFunction antiderivative(const ScalarTaylorFunction& f, uint k, ExactFloat c)
{
    ARIADNE_ASSERT(k<f.argument_size());
    ARIADNE_ASSERT(contains(f.domain()[k],c));

    ScalarTaylorFunction g = antiderivative(f,k);
    VectorTaylorFunction h = VectorTaylorFunction::identity(g.domain(),g.sweeper());
    h[k] = ScalarTaylorFunction::constant(g.domain(),c,g.sweeper());

    return g-compose(g,h);
}





Pair<ScalarTaylorFunction,ScalarTaylorFunction>
split(const ScalarTaylorFunction& tv, uint j)
{
    Pair<ValidatedTaylorModel,ValidatedTaylorModel> models=split(tv.model(),j);
    Pair<ExactBox,ExactBox> subdomains=split(tv.domain(),j);
    return make_pair(ScalarTaylorFunction(subdomains.first,models.first),
                     ScalarTaylorFunction(subdomains.second,models.second));

}

bool refines(const ScalarTaylorFunction& tv1, const ScalarTaylorFunction& tv2)
{
    if(tv1.domain()==tv2.domain()) { return refines(tv1.model(),tv2.model()); }
    if(subset(tv2.domain(),tv1.domain())) { return refines(restrict(tv1,tv2.domain()).model(),tv2.model()); }
    else { return false; }
}

bool disjoint(const ScalarTaylorFunction& tv1, const ScalarTaylorFunction& tv2)
{
    if(tv1.domain()==tv2.domain()) {
        return disjoint(tv1.model(),tv2.model());
    } else {
        ExactBox domain=intersection(tv1.domain(),tv2.domain());
        return disjoint(restrict(tv1,domain).model(),restrict(tv2,domain).model());
    }
}

ScalarTaylorFunction intersection(const ScalarTaylorFunction& tv1, const ScalarTaylorFunction& tv2)
{
    ARIADNE_ASSERT(tv1.domain()==tv2.domain());
    return ScalarTaylorFunction(tv1.domain(),intersection(tv1.model(),tv2.model()));
}

ErrorFloat norm(const ScalarTaylorFunction& f) {
    return norm(f.model());
}

ErrorFloat distance(const ScalarTaylorFunction& f1, const ScalarTaylorFunction& f2) {
    return norm(f1-f2);
}

ErrorFloat distance(const ScalarTaylorFunction& f1, const EffectiveScalarFunction& f2) {
    return distance(f1,ScalarTaylorFunction(f1.domain(),f2,f1.sweeper()));
}

Vector<ScalarTaylorFunction> compose(const Vector<ScalarTaylorFunction>& x, const Vector<ScalarTaylorFunction>& y) {
    ARIADNE_NOT_IMPLEMENTED; }
ScalarTaylorFunction compose(const ScalarTaylorFunction& x, const Vector<ScalarTaylorFunction>& y) {
    ARIADNE_NOT_IMPLEMENTED; }

Vector<ScalarTaylorFunction>
prod(const Matrix<ValidatedNumber>& A,
     const Vector<ScalarTaylorFunction>& x)
{
    ARIADNE_ASSERT(x.size()>0);
    ARIADNE_ASSERT(A.column_size()==x.size());
    for(uint i=0; i!=x.size(); ++i) { ARIADNE_ASSERT(x[i].argument_size()==x.zero_element().argument_size()); }

    Vector<ScalarTaylorFunction> r(A.row_size(),x.zero_element());
    for(uint i=0; i!=A.row_size(); ++i) {
        for(uint j=0; j!=A.column_size(); ++j) {
            //r[i]+=A[i][j]*x[j];
            const ValidatedNumber& Aij=A[i][j]; const ScalarTaylorFunction& xj=x[j]; ScalarTaylorFunction& ri=r[i]; ri+=Aij*xj;
        }
    }
    return r;
}

Matrix<ExactInterval>
jacobian(const Vector<ScalarTaylorFunction>& tv, const ExactBox& x);

ScalarTaylorFunction
midpoint(const ScalarTaylorFunction& f)
{
    ValidatedTaylorModel tm=f.model();
    tm.set_error(0u);
    return ScalarTaylorFunction(f.domain(),tm);
}


ScalarTaylorFunction
compose(const EffectiveScalarFunction& f, const Vector<ScalarTaylorFunction>& g)
{
    ARIADNE_ASSERT(f.argument_size()==g.size());
    for(uint i=0; i!=g.size(); ++i) {
        ARIADNE_ASSERT(g.zero_element().domain()==g[i].domain());
    }

    ExactBox gdomain=g.zero_element().domain();
    Vector<ValidatedTaylorModel> gmodels(g.size());
    for(uint i=0; i!=g.size(); ++i) { gmodels[i]=g[i].model(); }

    return ScalarTaylorFunction(gdomain,f.evaluate(gmodels));
}

std::ostream& operator<<(std::ostream& os, const Representation<Float>& flt_repr)
{
    const Float& flt=*flt_repr.pointer;
    int precision=os.precision(); std::ios_base::fmtflags flags = os.flags();
    os.precision(17); os.setf(std::ios_base::showpoint);
    os << "Float(" << flt << ")";
    os.precision(precision); os.flags(flags);
    return os;
}

std::ostream& operator<<(std::ostream& os, const Representation<UpperFloat>& flt_repr)
{
    return os << reinterpret_cast<Representation<Float>const&>(flt_repr);
}

std::ostream& operator<<(std::ostream& os, const Representation<ExactInterval>& ivl_repr)
{
    const ExactInterval& ivl=*ivl_repr.pointer;
    int precision=os.precision(); std::ios_base::fmtflags flags = os.flags();
    os.precision(17); os.setf(std::ios_base::showpoint);
    os << "ExactInterval("<<ivl.lower()<<","<<ivl.upper()<<")";
    os.precision(precision); os.flags(flags);
    return os;
}


std::ostream& operator<<(std::ostream& os, const Representation< Expansion<Float> >& exp_repr)
{
    const Expansion<Float>& exp=*exp_repr.pointer;
    int precision=os.precision(); std::ios_base::fmtflags flags = os.flags();
    os.precision(17); os.setf(std::ios_base::showpoint);
    os << "Expansion<Float>(" << exp.argument_size() << "," << exp.number_of_nonzeros();
    for(Expansion<Float>::const_iterator iter=exp.begin(); iter!=exp.end(); ++iter) {
        for(uint j=0; j!=iter->key().size(); ++j) {
            os << "," << uint(iter->key()[j]);
        }
        os << "," << iter->data();
    }
    os << ")";
    os.precision(precision); os.flags(flags);
    return os;
}

std::ostream& operator<<(std::ostream& os, const Representation< Expansion<ExactFloat> >& exp_repr) {
    return os << reinterpret_cast<Expansion<Float>const&>(exp_repr);
}

std::ostream& operator<<(std::ostream& os, const Representation< Expansion<ApproximateFloat> >& exp_repr) {
    return os << reinterpret_cast<Expansion<Float>const&>(exp_repr);
}

template<class X> std::ostream& operator<<(std::ostream& os, const Representation< Vector<X> >& vec_repr)
{
    const Vector<X>& vec=*vec_repr.pointer;
    ARIADNE_ASSERT(vec.size()!=0);
    os << "(";
    for(uint i=0; i!=vec.size(); ++i) {
        if(i!=0) { os << ","; }
        os << representation(vec[i]);
    }
    os << ")";
    return os;
}
template std::ostream& operator<< <Float> (std::ostream&, const Representation< Vector<Float> >&);
template std::ostream& operator<< <ExactInterval> (std::ostream&, const Representation< Vector<ExactInterval> >&);

std::ostream& operator<<(std::ostream& os, const Representation< ExactBox >& box_repr)
{
    const Vector<ExactInterval>& vec=*box_repr.pointer;
    ARIADNE_ASSERT(vec.size()!=0);
    os << "(";
    for(uint i=0; i!=vec.size(); ++i) {
        if(i!=0) { os << ","; }
        os << representation(vec[i]);
    }
    os << ")";
    return os;
}
std::ostream& operator<<(std::ostream& os, const Representation< Sweeper >& swp_repr)
{
    const Sweeper& swp=*swp_repr.pointer;
    int precision=os.precision(); std::ios_base::fmtflags flags = os.flags();
    os.precision(17); os.setf(std::ios_base::showpoint);
    os << swp;
    os.precision(precision); os.flags(flags);
    return os;
}

template<class T> std::ostream& operator<<(std::ostream& os, const Representation< List<T> >& lst_repr)
{
    const List<T>& lst=*lst_repr.pointer;
    ARIADNE_ASSERT(lst.size()!=0);
    os << "(";
    for(uint i=0; i!=lst.size(); ++i) {
        if(i!=0) { os << ","; }
        os << representation(lst[i]);
    }
    os << ")";
    return os;
}


std::ostream&
ScalarTaylorFunction::write(std::ostream& os) const
{
    Polynomial<ValidatedFloat> p=this->polynomial();
    Polynomial<ApproximateFloat> ap=p;
    os << ap;
    if(this->error()>0.0) { os << "+/-" << this->error(); }
    return os;
}

std::ostream&
ScalarTaylorFunction::repr(std::ostream& os) const
{
    return os << "ScalarTaylorFunction(" << representation(this->domain()) << ", " << representation(this->model().expansion().raw())
              << "," << representation(this->error().raw())<<","<<this->sweeper()<<")";
}

std::ostream&
operator<<(std::ostream& os, const ScalarTaylorFunction& tf)
{
    return tf.write(os);
}

std::ostream& operator<<(std::ostream& os, const Representation<ScalarTaylorFunction>& tf)
{
    return tf.pointer->repr(os);
}

/*
std::ostream& operator<<(std::ostream& os, const Representation<ScalarTaylorFunction>& frepr)
{
    ScalarTaylorFunction const& function=*frepr.pointer;
    ScalarTaylorFunction truncated_function=function;
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
std::ostream& operator<<(std::ostream& os, const ModelsRepresentation<ScalarTaylorFunction>& frepr)
{
    ScalarTaylorFunction const& f=*frepr.pointer;
    Float truncatation_error = 0.0;
    os << "<"<<f.domain()<<"\n";
    for(ValidatedTaylorModel::const_iterator iter=f.begin(); iter!=f.end(); ++iter) {
        if(abs(iter->data())>frepr.threshold) { truncatation_error+=abs(iter->data()); }
        else { os << iter->key() << ":" << iter->data() << ","; }
    }
    os << "+/-" << truncatation_error << "+/-" << f.error();
    return os;
}
*/

std::ostream& operator<<(std::ostream& os, const ModelsRepresentation<ScalarTaylorFunction>& frepr)
{
    ScalarTaylorFunction const& f=*frepr.pointer;
    ScalarTaylorFunction tf=f;
    tf.clobber();
    tf.sweep(ThresholdSweeper(frepr.threshold));
    os << "("<<tf.model()<<"+/-"<<f.error();
    return os;
}

std::ostream& operator<<(std::ostream& os, const PolynomialRepresentation<ScalarTaylorFunction>& frepr)
{
    ScalarTaylorFunction const& function=*frepr.pointer;
    ScalarTaylorFunction truncated_function = function;
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





bool
check(const Vector<ScalarTaylorFunction>& tv)
{
    for(uint i=0; i!=tv.size(); ++i) {
        if(tv.zero_element().domain()!=tv[i].domain()) { return false; }
    }
    return true;
}

Vector< Expansion<ExactFloat> >
expansion(const Vector<ScalarTaylorFunction>& x)
{
    Vector< Expansion<ExactFloat> > r(x.size());
    for(uint i=0; i!=x.size(); ++i) {
        r[i]=x[i].expansion();
    }
    return r;
}

Vector<ErrorFloat>
error(const Vector<ScalarTaylorFunction>& x)
{
    Vector<ErrorFloat> r(x.size());
    for(uint i=0; i!=x.size(); ++i) {
        r[i]=x[i].error();
    }
    return r;
}

Vector<CoefficientType>
value(const Vector<ScalarTaylorFunction>& x)
{
    Vector<CoefficientType> r(x.size());
    for(uint i=0; i!=x.size(); ++i) {
        r[i]=x[i].value();
    }
    return r;
}

Vector<UpperInterval>
ranges(const Vector<ScalarTaylorFunction>& x)
{
    Vector<UpperInterval> r(x.size());
    for(uint i=0; i!=x.size(); ++i) {
        r[i]=x[i].range();
    }
    return r;
}


Vector<ValidatedNumber>
evaluate(const Vector<ScalarTaylorFunction>& tv, const Vector<ValidatedNumber>& x)
{
    Vector<ValidatedNumber> r(tv.size());
    for(uint i=0; i!=tv.size(); ++i) {
        r[i]=evaluate(tv[i],x);
    }
    return r;
}

Matrix<ValidatedNumber>
jacobian(const Vector<ScalarTaylorFunction>& tv, const Vector<ValidatedNumber>& x)
{
    ARIADNE_ASSERT(check(tv));
    const Vector<ExactInterval> dom=tv.zero_element().domain();
    const uint n=dom.size();
    Vector< Differential<ValidatedNumber> > s(n,n,1u);
    for(uint j=0; j!=n; ++j) {
        ExactInterval dj=dom[j];
        s[j].set_value((x[j]-med_val(dj))/rad_val(dj));
        s[j].set_gradient(j,rec(rad_val(dj)));
    }
    Vector< Expansion<ExactFloat> > p=expansion(tv);
    Vector< Differential<ValidatedNumber> > d=evaluate(p,s);
    return d.jacobian();
}





VectorTaylorFunction::VectorTaylorFunction()
    : _domain(), _models()
{
}

VectorTaylorFunction::VectorTaylorFunction(uint k)
    : _domain(), _models(k)
{
}

VectorTaylorFunction::VectorTaylorFunction(uint m, const ExactBox& d, Sweeper swp)
    : _domain(d), _models(m,ValidatedTaylorModel(d.size(),swp))
{
}

VectorTaylorFunction::VectorTaylorFunction(uint k, const ScalarTaylorFunction& f)
    : _domain(f.domain()), _models(k,f.model())
{
}

VectorTaylorFunction::VectorTaylorFunction(const ValidatedVectorFunctionModel& f)
    : _domain(), _models()
{
    ARIADNE_ASSERT(dynamic_cast<const VectorTaylorFunction*>(&f.reference()));
    *this = dynamic_cast<const VectorTaylorFunction&>(f.reference());
}

VectorTaylorFunction& VectorTaylorFunction::operator=(const ValidatedVectorFunctionModel& f)
{
    ARIADNE_ASSERT(dynamic_cast<const VectorTaylorFunction*>(&f.reference()));
    *this = dynamic_cast<const VectorTaylorFunction&>(f.reference());
    return *this;
}


VectorTaylorFunction substitute(const VectorTaylorFunction& f, uint k, const ScalarTaylorFunction& h) {
    ARIADNE_ASSERT_MSG(f.argument_size()==h.argument_size()+1u,"f="<<f<<", k="<<k<<", h="<<h);
    const ExactInterval& dk=f.domain()[k];
    ValidatedTaylorModel unscaled_model=unscale(h.model(),dk);
    return VectorTaylorFunction(h.domain(),substitute(f.models(),k,unscaled_model));
}

ScalarTaylorFunction substitute(const ScalarTaylorFunction& f, uint k, const ScalarTaylorFunction& h) {
    const ExactInterval& dk=f.domain()[k];
    ValidatedTaylorModel unscaled_model=unscale(h.model(),dk);
    return ScalarTaylorFunction(h.domain(),substitute(f.model(),k,unscaled_model));
}


VectorTaylorFunction::VectorTaylorFunction(const ExactBox& d,
                                           const Vector<ValidatedTaylorModel>& f)
    : _domain(d), _models(f)
{
    for(uint i=0; i!=f.size(); ++i) {
        ARIADNE_ASSERT_MSG(d.size()==f[i].argument_size(),"d="<<d<<", f="<<f);
    }
}

VectorTaylorFunction::VectorTaylorFunction(const ExactBox& d,
                                           const Vector< Expansion<ExactFloat> >& f,
                                           const Vector<ErrorFloat>& e,
                                           Sweeper swp)
    : _domain(d), _models(f.size(),ValidatedTaylorModel(d.size(),swp))
{
    ARIADNE_ASSERT(f.size()==e.size());
    for(uint i=0; i!=f.size(); ++i) {
        ARIADNE_ASSERT(d.size()==f[i].argument_size());
        _models[i]=ValidatedTaylorModel(f[i],e[i],swp);
    }
}

VectorTaylorFunction::VectorTaylorFunction(const ExactBox& d,
                                           const Vector< Expansion<ExactFloat> >& f,
                                           Sweeper swp)
    : VectorTaylorFunction(d,f,Vector<ErrorFloat>(f.size()),swp)
{
}

VectorTaylorFunction::VectorTaylorFunction(const ExactBox& d,
                                           const Vector< Expansion<Float> >& f,
                                           const Vector<Float>& e,
                                           Sweeper swp)
    : VectorTaylorFunction(d,reinterpret_cast<Vector<Expansion<ExactFloat>>const&>(f),
                           reinterpret_cast<Vector<ErrorFloat>const&>(e),swp)
{
}

VectorTaylorFunction::VectorTaylorFunction(const ExactBox& d,
                                           const Vector< Expansion<Float> >& f,
                                           Sweeper swp)
    : VectorTaylorFunction(d,reinterpret_cast<Vector<Expansion<ExactFloat>>const&>(f),Vector<ErrorFloat>(f.size()),swp)
{
}



VectorTaylorFunction::VectorTaylorFunction(const ExactBox& d,
                                           const ValidatedVectorFunction& f,
                                           const Sweeper& swp)
    : _domain(d), _models(f.result_size())
{
    //ARIADNE_ASSERT_MSG(f.result_size()>0, "d="<<d<<", f="<<f<<", swp="<<swp);
    ARIADNE_ASSERT(d.size()==f.argument_size());
    Vector<ValidatedTaylorModel> x=ValidatedTaylorModel::scalings(d,swp);
    this->_models=f.evaluate(x);
    this->sweep();
}

/*
VectorTaylorFunction::VectorTaylorFunction(const ExactBox& d,
                                           const EffectiveVectorFunction& f,
                                           const Sweeper& swp)
    : _domain(d), _models(f.result_size())
{
    ARIADNE_ASSERT(f.result_size()>0);
    ARIADNE_ASSERT(d.size()==f.argument_size());
    Vector<ValidatedTaylorModel> x=ValidatedTaylorModel::scalings(d,swp);
    for(uint i=0; i!=x.size(); ++i) { x[i].set_sweeper(swp); }
    this->_models=f.evaluate(x);
    this->sweep();
}
*/

VectorTaylorFunction::VectorTaylorFunction(const Vector<ScalarTaylorFunction>& v)
    : _domain(), _models(v.size())
{
    ARIADNE_ASSERT(v.size()>0);
    for(uint i=0; i!=v.size(); ++i) { ARIADNE_ASSERT(v[i].domain()==v.zero_element().domain()); }
    this->_domain=v.zero_element().domain();
    for(uint i=0; i!=v.size(); ++i) {
        this->_models[i]=v[i].model();
    }
}

VectorTaylorFunction::VectorTaylorFunction(const List<ScalarTaylorFunction>& v)
    : _domain(), _models(v.size())
{
    ARIADNE_ASSERT(v.size()>0);
    for(uint i=1; i!=v.size(); ++i) { ARIADNE_ASSERT(v[i].domain()==v[0].domain()); }
    this->_domain=v[0].domain();
    for(uint i=0; i!=v.size(); ++i) {
        this->_models[i]=v[i].model();
    }
}

VectorTaylorFunction::VectorTaylorFunction(std::initializer_list<ScalarTaylorFunction> lst)
    : _domain(), _models(lst.size())
{
    *this=VectorTaylorFunction(List<ScalarTaylorFunction>(lst));
}


VectorTaylorFunction* VectorTaylorFunction::_clone() const
{
    return new VectorTaylorFunction(*this);
}

VectorTaylorFunction* VectorTaylorFunction::_create() const
{
    return new VectorTaylorFunction(this->result_size(), ScalarTaylorFunction(this->domain(),this->sweeper()));
}

ScalarTaylorFunction* VectorTaylorFunction::_create_zero() const
{
    return new ScalarTaylorFunction(this->domain(),this->sweeper());
}

VectorTaylorFunction* VectorTaylorFunction::_create_identity() const
{
    Sweeper sweeper=this->sweeper();
    VectorTaylorFunction* result = new VectorTaylorFunction(this->domain().size(), ScalarTaylorFunction(this->domain(),sweeper));
    for(uint i=0; i!=result->size(); ++i) { (*result)[i]=ScalarTaylorFunction::coordinate(this->domain(),i,sweeper); }
    return result;
}

Void VectorTaylorFunction::adjoin(const ScalarTaylorFunction& sf)
{
    ARIADNE_ASSERT_MSG(sf.domain()==this->domain(),"sf="<<sf);
    this->_models=join(this->_models,sf.model());
}







VectorTaylorFunction
VectorTaylorFunction::constant(const ExactBox& d, const Vector<ExactNumber>& c, Sweeper swp)
{
    return VectorTaylorFunction(d,ValidatedTaylorModel::constants(d.size(),c,swp));
}

VectorTaylorFunction
VectorTaylorFunction::constant(const ExactBox& d, const Vector<ValidatedNumber>& c, Sweeper swp)
{
    return VectorTaylorFunction(d,ValidatedTaylorModel::constants(d.size(),c,swp));
}

VectorTaylorFunction
VectorTaylorFunction::identity(const ExactBox& d, Sweeper swp)
{
    return VectorTaylorFunction(d,ValidatedTaylorModel::scalings(d,swp));
}

VectorTaylorFunction
VectorTaylorFunction::projection(const ExactBox& d, uint imin, uint imax, Sweeper swp)
{
    return VectorTaylorFunction(ScalarTaylorFunction::variables(d,imin,imax,swp));
}


Polynomial<ValidatedNumber> polynomial(const ValidatedTaylorModel& tm) {
    return Polynomial<ValidatedNumber>(tm.expansion())+ValidatedNumber(-tm.error(),+tm.error());
}

Vector< Polynomial<ValidatedNumber> >
VectorTaylorFunction::polynomials() const
{
    Vector<Polynomial<ValidatedNumber> > p(this->result_size(),Polynomial<ValidatedNumber>(this->argument_size()));
    for(uint i=0; i!=this->result_size(); ++i) {
        p[i]=static_cast<ScalarTaylorFunction>((*this)[i]).polynomial();
    }
    return p;
}

Vector< Expansion<ExactFloat> > const
VectorTaylorFunction::expansions() const
{
    Vector< Expansion<ExactFloat> > e(this->result_size(),Expansion<ExactFloat>(this->argument_size()));
    for(uint i=0; i!=this->result_size(); ++i) {
        e[i]=this->models()[i].expansion();
    }
    return e;
}

Vector<ErrorFloat> const
VectorTaylorFunction::errors() const
{
    Vector<ErrorFloat> e(this->result_size());
    for(uint i=0; i!=this->result_size(); ++i) {
        e[i]=this->models()[i].error();
    }
    return e;
}

ErrorFloat const
VectorTaylorFunction::error() const
{
    ErrorFloat e=0;
    for(uint i=0; i!=this->result_size(); ++i) {
        e=max(e,this->models()[i].error());
    }
    return e;
}

ValidatedVectorFunction
VectorTaylorFunction::function() const
{
    return ValidatedVectorFunction(new VectorTaylorFunction(*this));
}

bool
VectorTaylorFunction::operator==(const VectorTaylorFunction& tm) const
{
    return this->_models==tm._models;
}



bool
VectorTaylorFunction::operator!=(const VectorTaylorFunction& p2) const
{
    return !(*this==p2);
}



Sweeper
VectorTaylorFunction::sweeper() const
{
    ARIADNE_ASSERT(this->size()>0); return this->_models[0].sweeper();
}


void
VectorTaylorFunction::set_sweeper(Sweeper swp)
{
    for(uint i=0; i!=this->result_size(); ++i) {
        this->_models[i].set_sweeper(swp);
    }
}

const ExactBox&
VectorTaylorFunction::domain() const
{
    return this->_domain;
}

const ExactBox
VectorTaylorFunction::codomain() const
{
    ExactBox result(this->result_size());
    for(uint i=0; i!=result.size(); ++i) {
        result[i]=this->_models[i].codomain();
    }
    return result;
}


const UpperBox
VectorTaylorFunction::range() const
{
    Vector<UpperInterval> result(this->result_size());
    for(uint i=0; i!=result.size(); ++i) {
        result[i]=this->_models[i].range();
    }
    return UpperBox(result);
}


const Vector<ExactFloat>
VectorTaylorFunction::centre() const
{
    Vector<ExactFloat> result(this->result_size());
    for(uint i=0; i!=result.size(); ++i) {
        result[i]=this->_models[i].value();
    }
    return result;
}


const Vector<ValidatedTaylorModel>&
VectorTaylorFunction::models() const
{
    return this->_models;
}

const ValidatedTaylorModel&
VectorTaylorFunction::model(uint i) const
{
    return this->_models[i];
}

ValidatedTaylorModel&
VectorTaylorFunction::model(uint i)
{
    return this->_models[i];
}




uint
VectorTaylorFunction::argument_size() const
{
    return this->_domain.size();
}


uint
VectorTaylorFunction::result_size() const
{
    return this->_models.size();
}


ScalarTaylorFunction
VectorTaylorFunction::operator[](uint i) const
{
    return this->get(i);
}

VectorTaylorFunctionElementReference
VectorTaylorFunction::operator[](uint i)
{
    return VectorTaylorFunctionElementReference(*this,i);
}

ScalarTaylorFunction
VectorTaylorFunction::get(uint i) const
{
    return ScalarTaylorFunction(this->_domain,this->_models[i]);
}

void
VectorTaylorFunction::set(uint i, const ScalarTaylorFunction& e)
{
    ARIADNE_ASSERT_MSG(this->size()>i,"Cannot set "<<i<<"th element of VectorTaylorFunction "<<(*this));
    if(this->domain().size()!=0) {
        ARIADNE_ASSERT_MSG(e.domain()==this->domain(),"Domain of "<<e<<" conflicts with existing domain "<<this->domain());
    } else {
        this->_domain=e.domain();
    }
    this->_models[i]=e.model();
}














template<class T>
void
VectorTaylorFunction::_compute(Vector<T>& r, const Vector<T>& a) const
{
    typedef typename T::NumericType X;
    const VectorTaylorFunction& f=*this;
    Vector<T> sx=Ariadne::unscale(a,f._domain);
    for(uint i=0; i!=r.size(); ++i) {
        T ri=Ariadne::evaluate(this->_models[i].expansion(),sx);
        X e=convert_error<X>(this->_models[i].error());
        r[i]=ri+e;
    }
}




VectorTaylorFunction&
VectorTaylorFunction::sweep()
{
    for(uint i=0; i!=this->size(); ++i) {
        this->_models[i].sweep();
    }
    return *this;
}

VectorTaylorFunction&
VectorTaylorFunction::sweep(const SweeperInterface& sweeper)
{
    for(uint i=0; i!=this->size(); ++i) {
        this->_models[i].sweep(sweeper);
    }
    return *this;
}


Void
VectorTaylorFunction::clobber()
{
    for(uint i=0; i!=this->size(); ++i) {
        this->_models[i].clobber();
    }
}



Vector<ApproximateNumber>
VectorTaylorFunction::evaluate(const Vector<ApproximateNumber>& x) const
{
    const VectorTaylorFunction& f=*this;
    if(!decide(contains(f.domain(),x))) {
        ARIADNE_THROW(DomainException,"tf.evaluate(ax) with tf="<<f<<", ax="<<x,"ax is not an element of tf.domain()="<<f.domain());
    }
    Vector<ApproximateNumber> sx=Ariadne::unscale(x,f._domain);
    Vector<ApproximateNumber> r(this->result_size());
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=Ariadne::evaluate(this->_models[i].expansion(),sx);
    }
    return r;
}


Vector<ValidatedNumber>
VectorTaylorFunction::evaluate(const Vector<ValidatedNumber>& x) const
{
    const VectorTaylorFunction& f=*this;
    if(!contains(f.domain(),x)) {
        ARIADNE_THROW(DomainException,"tf.evaluate(vx) with tf="<<f<<", x="<<x,"vx is not a subset of tf.domain()="<<f.domain());
    }
    Vector<ValidatedNumber> sx=Ariadne::unscale(x,f._domain);
    return Ariadne::evaluate(f._models,sx);
}

Vector<ApproximateNumber>
VectorTaylorFunction::operator()(const Vector<ApproximateNumber>& x) const
{
    return this->evaluate(x);
}

Vector<ValidatedNumber>
VectorTaylorFunction::operator()(const Vector<ValidatedNumber>& x) const
{
    return this->evaluate(x);
}

Matrix<ValidatedNumber>
VectorTaylorFunction::jacobian(const Vector<ValidatedNumber>& x) const
{
    Matrix<ValidatedNumber> J=Ariadne::jacobian(this->_models,unscale(x,this->_domain));
    for(uint j=0; j!=J.column_size(); ++j) {
        ValidatedNumber rad=rad_val(this->_domain[j]);
        for(uint i=0; i!=J.row_size(); ++i) {
            J[i][j]/=rad;
        }
    }
    return J;
}

Void
VectorTaylorFunction::restrict(const ExactBox& x)
{
    *this=Ariadne::restrict(*this,x);
}


VectorTaylorFunction
join(const VectorTaylorFunction& f1, const ScalarTaylorFunction& f2)
{
    ARIADNE_ASSERT_MSG(f1.domain()==f2.domain(),"f1="<<f1<<", f2="<<f2);
    return VectorTaylorFunction(f1.domain(),join(f1.models(),f2.model()));
}

VectorTaylorFunction
join(const VectorTaylorFunction& f, const VectorTaylorFunction& g)
{
    ARIADNE_ASSERT(f.domain()==g.domain());
    return VectorTaylorFunction(f.domain(),join(f.models(),g.models()));
}

VectorTaylorFunction
join(const ScalarTaylorFunction& f1, const ScalarTaylorFunction& f2)
{
    ARIADNE_ASSERT(f1.domain()==f2.domain());
    return VectorTaylorFunction(f1.domain(),join(f1.model(),f2.model()));
}

VectorTaylorFunction
join(const ScalarTaylorFunction& f1, const VectorTaylorFunction& f2)
{
    ARIADNE_ASSERT(f1.domain()==f2.domain());
    return VectorTaylorFunction(f1.domain(),join(f1.model(),f2.models()));
}

VectorTaylorFunction
combine(const ScalarTaylorFunction& f1, const ScalarTaylorFunction& f2)
{
    return VectorTaylorFunction(product(f1.domain(),f2.domain()),combine(f1.model(),f2.model()));
}

VectorTaylorFunction
combine(const ScalarTaylorFunction& f1, const VectorTaylorFunction& f2)
{
    return VectorTaylorFunction(product(f1.domain(),f2.domain()),combine(f1.model(),f2.models()));
}

VectorTaylorFunction
combine(const VectorTaylorFunction& f1, const ScalarTaylorFunction& f2)
{
    return VectorTaylorFunction(product(f1.domain(),f2.domain()),combine(f1.models(),f2.model()));
}

VectorTaylorFunction
combine(const VectorTaylorFunction& f1, const VectorTaylorFunction& f2)
{
    return VectorTaylorFunction(product(f1.domain(),f2.domain()),combine(f1.models(),f2.models()));
}


VectorTaylorFunction
embed(const VectorTaylorFunction& f, const ExactInterval& d)
{
    return embed(ExactBox(),f,ExactBox(1u,d));
}

VectorTaylorFunction
embed(const VectorTaylorFunction& f, const ExactBox& d)
{
    return embed(ExactBox(),f,d);
}

VectorTaylorFunction
embed(const ExactBox& d, const VectorTaylorFunction& f)
{
    return embed(d,f,ExactBox());
}

VectorTaylorFunction
embed(const ExactBox& d1, const VectorTaylorFunction& f, const ExactBox& d2)
{
    return VectorTaylorFunction(join(d1,f.domain(),d2),embed(embed(d1.size(),f.models()),d2.size()));
}

VectorTaylorFunction
restrict(const VectorTaylorFunction& f, const ExactBox& d)
{
    ARIADNE_ASSERT_MSG(subset(d,f.domain()),"Cannot restrict "<<f<<" to non-sub-domain "<<d);
    if(d==f.domain()) { return f; }
    VectorTaylorFunction r(f.result_size(),d,f.sweeper());
    for(uint i=0; i!=r.result_size(); ++i) {
        r.set(i,restrict(f[i],d));
    }
    return r;
}

Pair<VectorTaylorFunction,VectorTaylorFunction>
split(const VectorTaylorFunction& tf, uint j)
{
    Pair< Vector<ValidatedTaylorModel>,Vector<ValidatedTaylorModel> > models=split(tf.models(),j);
    Pair<ExactBox,ExactBox> subdomains=split(tf.domain(),j);
    return make_pair(VectorTaylorFunction(subdomains.first,models.first),
                     VectorTaylorFunction(subdomains.second,models.second));

}

bool
refines(const VectorTaylorFunction& f1, const VectorTaylorFunction& f2) {
    ARIADNE_ASSERT(f1.result_size()==f2.result_size());
    for(uint i=0; i!=f1.result_size(); ++i) {
        if(!refines(f1[i],f2[i])) { return false; }
    }
    return true;
}

bool
disjoint(const VectorTaylorFunction& f1, const VectorTaylorFunction& f2) {
    ARIADNE_ASSERT(f1.result_size()==f2.result_size());
    for(uint i=0; i!=f1.result_size(); ++i) {
        if(disjoint(f1[i],f2[i])) { return true; }
    }
    return false;
}

VectorTaylorFunction
intersection(const VectorTaylorFunction& f1, const VectorTaylorFunction& f2) {
    ARIADNE_ASSERT(f1.result_size()==f2.result_size());
    VectorTaylorFunction r(f1.result_size());
    for(uint i=0; i!=r.result_size(); ++i) {
        r[i]=intersection(f1[i],f2[i]);
    }
    return r;
}

VectorTaylorFunction&
operator+=(VectorTaylorFunction& f, const VectorTaylorFunction& g)
{
    ARIADNE_ASSERT(f.result_size()==g.result_size());
    ARIADNE_ASSERT(subset(f.domain(),g.domain()));
    ARIADNE_ASSERT(f.domain()==g.domain());
    f._models+=g._models;
    return f;
}

VectorTaylorFunction&
operator-=(VectorTaylorFunction& f, const VectorTaylorFunction& g)
{
    ARIADNE_ASSERT(f.result_size()==g.result_size());
    ARIADNE_ASSERT(subset(f.domain(),g.domain()));
    ARIADNE_ASSERT(f.domain()==g.domain());
    f._models+=g._models;
    return f;
}

VectorTaylorFunction&
operator+=(VectorTaylorFunction& f, const Vector<ValidatedNumber>& c)
{
    ARIADNE_ASSERT(f.result_size()==c.size());
    f._models+=c;
    return f;
}

VectorTaylorFunction&
operator-=(VectorTaylorFunction& f, const Vector<ValidatedNumber>& c)
{
    ARIADNE_ASSERT(f.result_size()==c.size());
    f._models-=c;
    return f;
}

VectorTaylorFunction&
operator*=(VectorTaylorFunction& f, const ValidatedNumber& c)
{
    f._models*=c;
    return f;
}

VectorTaylorFunction&
operator/=(VectorTaylorFunction& f, const ValidatedNumber& c)
{
    f._models/=c;
    return f;
}


VectorTaylorFunction
operator+(const VectorTaylorFunction& f1, const VectorTaylorFunction& f2)
{
    ARIADNE_ASSERT_MSG(!empty(intersection(f1.domain(),f2.domain())),
                       "operator+(VectorTaylorFunction f1, VectorTaylorFunction f2) with f1="<<f1<<" f2="<<f2<<
                       ": domains are disjoint");
    if(f1.domain()==f2.domain()) {
        return VectorTaylorFunction(f1.domain(),Vector<ValidatedTaylorModel>(f1.models()+f2.models()));
    } else {
        ExactBox new_domain=intersection(f1.domain(),f2.domain());
        return operator+(restrict(f1,new_domain),restrict(f2,new_domain));
    }
}


VectorTaylorFunction
operator-(const VectorTaylorFunction& f1, const VectorTaylorFunction& f2)
{
    ARIADNE_ASSERT(!empty(intersection(f1.domain(),f2.domain())));
    if(f1.domain()==f2.domain()) {
        return VectorTaylorFunction(f1.domain(),Vector<ValidatedTaylorModel>(f1.models()-f2.models()));
    } else {
        ExactBox new_domain=intersection(f1.domain(),f2.domain());
        return operator-(restrict(f1,new_domain),restrict(f2,new_domain));
    }
}

VectorTaylorFunction
operator*(const VectorTaylorFunction& f1, const ScalarTaylorFunction& f2)
{
    ARIADNE_ASSERT(!empty(intersection(f1.domain(),f2.domain())));
    if(f1.domain()==f2.domain()) {
        return VectorTaylorFunction(f1.domain(),Vector<ValidatedTaylorModel>(f1.models()*f2.model()));
    } else {
        ExactBox new_domain=intersection(f1.domain(),f2.domain());
        return operator*(restrict(f1,new_domain),restrict(f2,new_domain));
    }
}

VectorTaylorFunction
operator*(const ScalarTaylorFunction& f1, const VectorTaylorFunction& f2)
{
    ARIADNE_ASSERT(!empty(intersection(f1.domain(),f2.domain())));
    if(f1.domain()==f2.domain()) {
        return VectorTaylorFunction(f1.domain(),Vector<ValidatedTaylorModel>(f1.model()*f2.models()));
    } else {
        ExactBox new_domain=intersection(f1.domain(),f2.domain());
        return operator*(restrict(f1,new_domain),restrict(f2,new_domain));
    }
}



VectorTaylorFunction
operator-(const VectorTaylorFunction& f)
{
    return VectorTaylorFunction(f.domain(),Vector<ValidatedTaylorModel>(-f.models()));
}

VectorTaylorFunction
operator*(const ValidatedNumber& c, const VectorTaylorFunction& f)
{
    return VectorTaylorFunction(f.domain(),Vector<ValidatedTaylorModel>(f.models()*c));
}

VectorTaylorFunction
operator*(const VectorTaylorFunction& f, const ValidatedNumber& c)
{
    return VectorTaylorFunction(f.domain(),Vector<ValidatedTaylorModel>(f.models()*c));
}

VectorTaylorFunction
operator/(const VectorTaylorFunction& f, const ValidatedNumber& c)
{
    return VectorTaylorFunction(f.domain(),Vector<ValidatedTaylorModel>(f.models()/c));
}

VectorTaylorFunction
operator+(const VectorTaylorFunction& f, const Vector<ValidatedNumber>& c)
{
    return VectorTaylorFunction(f.domain(),Vector<ValidatedTaylorModel>(f.models()+c));
}

VectorTaylorFunction
operator-(const VectorTaylorFunction& f, const Vector<ValidatedNumber>& c)
{
    return VectorTaylorFunction(f.domain(),Vector<ValidatedTaylorModel>(f.models()-c));
}

VectorTaylorFunction
operator*(const Matrix<Float>& A, const VectorTaylorFunction& f)
{
    ARIADNE_PRECONDITION(A.column_size()==f.size());
    Vector<ValidatedTaylorModel> models(A.row_size(),ValidatedTaylorModel(f.argument_size(),f.sweeper()));
    for(size_t i=0; i!=A.row_size(); ++i) {
        for(size_t j=0; j!=A.column_size(); ++j) {
            models[i] += A.get(i,j) * f.model(j);
        }
    }
    return VectorTaylorFunction(f.domain(),models);
}

VectorTaylorFunction
operator*(const Matrix<ValidatedNumber>& A, const VectorTaylorFunction& f)
{
    ARIADNE_PRECONDITION(A.column_size()==f.size());
    Vector<ValidatedTaylorModel> models(A.row_size(),ValidatedTaylorModel(f.argument_size(),f.sweeper()));
    for(size_t i=0; i!=A.row_size(); ++i) {
        for(size_t j=0; j!=A.column_size(); ++j) {
            models[i] += A.get(i,j) * f.model(j);
        }
    }
    return VectorTaylorFunction(f.domain(),models);
}


VectorTaylorFunction
operator-(const VectorTaylorFunction& f1, const EffectiveVectorFunction& f2) {
    return f1 - VectorTaylorFunction(f1.domain(),f2,f1.sweeper());
}





VectorTaylorFunction
partial_evaluate(const VectorTaylorFunction& tf, uint k, const ExactNumber& c)
{
    return partial_evaluate(tf,k,ValidatedNumber(c));
}

VectorTaylorFunction
partial_evaluate(const VectorTaylorFunction& tf, uint k, const ValidatedNumber& c)
{
    // Scale c to domain
    const uint as=tf.argument_size();
    ARIADNE_ASSERT(k<as);
    const Vector<ExactInterval>& domain=tf.domain();
    const ExactInterval& dk=domain[k];
    ValidatedNumber sc=(c-med_val(dk))/rad_val(dk);

    Vector<ExactInterval> new_domain(as-1);
    for(uint i=0; i!=k; ++i) { new_domain[i]=domain[i]; }
    for(uint i=k; i!=as-1; ++i) { new_domain[i]=domain[i+1]; }

    Vector<ValidatedTaylorModel> new_models=partial_evaluate(tf.models(),k,sc);

    return VectorTaylorFunction(new_domain,new_models);
}


VectorTaylorFunction
partial_restrict(const VectorTaylorFunction& tf, uint k, const ExactInterval& d)
{
    VectorTaylorFunction r(tf.result_size(),tf.domain(),tf.sweeper());
    for(uint i=0; i!=tf.result_size(); ++i) {
        r[i]=partial_restrict(tf[i],k,d);
    }
    return r;
}

VectorTaylorFunction
restrict(const VectorTaylorFunction& tf, uint k, const ExactInterval& d)
{
    return partial_restrict(tf,k,d);
}


Vector<ValidatedNumber>
evaluate(const VectorTaylorFunction& f, const Vector<ValidatedNumber>& x) {
    if(!contains(f.domain(),x)) {
        ARIADNE_THROW(DomainException,"evaluate(f,x) with f="<<f<<", x="<<x,"x is not a subset of f.domain()="<<f.domain());
    }
    return unchecked_evaluate(f,x);
}

Vector<ValidatedNumber>
unchecked_evaluate(const VectorTaylorFunction& f, const Vector<ValidatedNumber>& x) {
    return evaluate(f.models(),unscale(x,f.domain()));
}

VectorTaylorFunction
compose(const EffectiveVectorFunction& g, const VectorTaylorFunction& f)
{
    return VectorTaylorFunction(f.domain(),g.evaluate(f.models()));
}

VectorTaylorFunction
compose(const ValidatedVectorFunction& g, const VectorTaylorFunction& f)
{
    return VectorTaylorFunction(f.domain(),g.evaluate(f.models()));
}

VectorTaylorFunction
compose(const VectorTaylorFunction& g, const VectorTaylorFunction& f)
{
    if(!subset(f.codomain(),g.domain())) {
        ARIADNE_THROW(DomainException,"compose(g,f) with g="<<g<<", f="<<f,"f.codomain()="<<f.codomain()<<" is not a subset of g.domain()="<<g.domain());
    }
    return unchecked_compose(g,f);
}


VectorTaylorFunction
unchecked_compose(const VectorTaylorFunction& g, const VectorTaylorFunction& f)
{
    return VectorTaylorFunction(f.domain(),compose(g.models(),unscale(f.models(),g.domain())));
}



VectorTaylorFunction
derivative(const VectorTaylorFunction& f, uint k)
{
    ARIADNE_ASSERT_MSG(k<f.argument_size(),"f="<<f<<"\n f.argument_size()="<<f.argument_size()<<" k="<<k);
    ValidatedNumber fdomkrad=rad_val(f.domain()[k]);
    VectorTaylorFunction g=f;
    for(uint i=0; i!=g.size(); ++i) {
        g[i]=derivative(f[i],k);
    }
    return g;
}

VectorTaylorFunction
antiderivative(const VectorTaylorFunction& f, uint k)
{
    ARIADNE_ASSERT_MSG(k<f.argument_size(),"f="<<f<<"\n f.argument_size()="<<f.argument_size()<<" k="<<k);
    ValidatedNumber fdomkrad=rad_val(f.domain()[k]);
    VectorTaylorFunction g=f;
    for(uint i=0; i!=g.size(); ++i) {
        g._models[i].antidifferentiate(k);
        g._models[i]*=fdomkrad;
    }
    return g;
}

VectorTaylorFunction
antiderivative(const VectorTaylorFunction& f, uint k, ExactNumber c)
{
    ARIADNE_ASSERT_MSG(k<f.argument_size(),"f="<<f<<"\n f.argument_size()="<<f.argument_size()<<" k="<<k);
    ValidatedNumber fdomkrad=rad_val(f.domain()[k]);
    VectorTaylorFunction g=f;
    for(uint i=0; i!=g.size(); ++i) {
        g[i]=antiderivative(f[i],k,c);
    }
    return g;
}






ErrorFloat norm(const VectorTaylorFunction& f) {
    ErrorFloat res=0u;
    for(uint i=0; i!=f.result_size(); ++i) {
        res=max(res,norm(f[i]));
    }
    return res;
}

ErrorFloat distance(const VectorTaylorFunction& f1, const VectorTaylorFunction& f2) {
    return norm(f1-f2);
}

ErrorFloat distance(const VectorTaylorFunction& f1, const EffectiveVectorFunction& f2) {
    return distance(f1,VectorTaylorFunction(f1.domain(),f2,f1.sweeper()));
}


std::ostream&
VectorTaylorFunction::write(std::ostream& os) const
{
    os << "[";
    for(uint i=0; i!=this->result_size(); ++i) {
        if(i!=0) { os << ", "; }
        ScalarTaylorFunction tfi=(*this)[i];
        os << tfi;
    }
    os << "]";
    return os;
}

std::ostream&
VectorTaylorFunction::repr(std::ostream& os) const
{
    return os << "VectorTaylorFunction("
              << representation(this->domain()) << ", " << representation(this->expansions()) << ", "
              << representation(this->errors()) << ", " << representation(this->sweeper()) << ")";
}

std::ostream& operator<<(std::ostream& os, const Representation<VectorTaylorFunction>& repr)
{
    return repr.pointer->repr(os);
}

std::ostream& operator<<(std::ostream& os, const PolynomialRepresentation<VectorTaylorFunction>& repr)
{
    const VectorTaylorFunction& function = *repr.pointer;
    os << "[";
    for(uint i=0; i!=function.result_size(); ++i) {
        if(i!=0) { os << ","; }
        os << polynomial_repr(function[i],repr.threshold,repr.names);
    }
    return os << "]";
}

std::ostream& operator<<(std::ostream& os, const PolynomialRepresentation< List<ScalarTaylorFunction> >& repr)
{
    const List<ScalarTaylorFunction>& functions = *repr.pointer;
    os << "[";
    for(uint i=0; i!=functions.size(); ++i) {
        if(i!=0) { os << ","; }
        os << polynomial_repr(functions[i],repr.threshold);
    }
    return os << "]";
}



std::ostream&
operator<<(std::ostream& os, const VectorTaylorFunction& p)
{
    return p.write(os);
}

Polynomial<ValidatedNumber> polynomial(const ScalarTaylorFunction& tfn) {
    return Polynomial<ValidatedNumber>(tfn.polynomial());
}

Vector< Polynomial<ValidatedNumber> > polynomial(const VectorTaylorFunction& tfn) {
    return tfn.polynomials();
}

List< Polynomial<ValidatedNumber> > polynomials(const List<ScalarTaylorFunction>& tfns) {
    List< Polynomial<ValidatedNumber> > result;
    for(uint i=0; i!=tfns.size(); ++i) {
        result.append(polynomial(tfns[i]));
    }
    return result;
}


ScalarTaylorFunction
TaylorFunctionFactory::create(const ExactBox& domain, const ValidatedScalarFunctionInterface& function) const
{
    return ScalarTaylorFunction(domain,function,this->_sweeper);
}

ScalarTaylorFunction*
TaylorFunctionFactory::_create(const ExactBox& domain, const ValidatedScalarFunctionInterface& function) const
{
    return new ScalarTaylorFunction(domain,function,this->_sweeper);
}

VectorTaylorFunction
TaylorFunctionFactory::create(const ExactBox& domain, const ValidatedVectorFunctionInterface& function) const
{
    return VectorTaylorFunction(domain,function,this->_sweeper);
}

VectorTaylorFunction*
TaylorFunctionFactory::_create(const ExactBox& domain, const ValidatedVectorFunctionInterface& function) const
{
    return new VectorTaylorFunction(domain,function,this->_sweeper);
}

ScalarTaylorFunction
TaylorFunctionFactory::create_zero(const ExactBox& domain) const
{
    return ScalarTaylorFunction::zero(domain,this->_sweeper);
}

ScalarTaylorFunction
TaylorFunctionFactory::create_constant(const ExactBox& domain, ValidatedNumber value) const
{
    return ScalarTaylorFunction::constant(domain,value,this->_sweeper);
}

ScalarTaylorFunction
TaylorFunctionFactory::create_coordinate(const ExactBox& domain, uint k) const
{
    return ScalarTaylorFunction::coordinate(domain,k,this->_sweeper);
}

ScalarTaylorFunction
TaylorFunctionFactory::create_identity(const ExactInterval& domain) const
{
    return ScalarTaylorFunction::identity(domain,this->_sweeper);
}

VectorTaylorFunction
TaylorFunctionFactory::create_identity(const ExactBox& domain) const
{
    return VectorTaylorFunction::identity(domain,this->_sweeper);
}


FunctionModelFactoryInterface<ValidatedTag>* make_taylor_function_factory() {
    return new TaylorFunctionFactory(Sweeper());
}

FunctionModelFactoryInterface<ValidatedTag>* make_taylor_function_factory(double sweep_threshold) {
    return new TaylorFunctionFactory(ThresholdSweeper(sweep_threshold));
}

/*
latexstream&
operator<<(Output::latexstream& texs, const VectorTaylorFunction& p)
{
    using namespace Function;
    texs << "%VectorTaylorFunction\n";
    texs << "\\ensuremath{\n";
    texs << "\\left( \\begin{Array}{c}\n";
    char var='x';
    for(uint i=0; i!=p.result_size(); ++i) {
        bool first = true;
        if(i!=0) { texs << "\\\\"; }
        for(MultiIndex j(p.argument_size()); j.degree()<=p.order(); ++j) {
            const ValidatedNumber& a=p.centre_derivatives()[i][j];
            if(a!=0) {
                if(first) { first=false; }
                else { if(a>0) { texs << '+'; } }
                if(a==1) { if(j.degree()==0) { texs << a; } }
                else if(a==-1) { if(j.degree()==0) { texs << a; } else { texs << '-'; } }
                else { texs << a << ' '; }
                for(uint k=0; k!=p.argument_size(); ++k) {
                    if(j[k]!=0) {
                        texs << var << "_{ " << k+1 << "}";
                        if(j[k]!=1) {
                            texs << "^{" << j[k] << "}";
                        }
                    }
                    texs << " ";
                }
            }
        }
        texs << "\n";
    }
    texs << "\\end{Array}\\right)\n}\n";
    return texs;
}
*/


} // namespace Ariadne
