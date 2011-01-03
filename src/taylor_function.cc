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

#include <iostream>
#include <iomanip>

#include "macros.h"
#include "exceptions.h"
#include "numeric.h"
#include "vector.h"
#include "matrix.h"
#include "multi_index.h"
#include "polynomial.h"
#include "differential.h"

#include "function.h"
#include "taylor_function.h"

#include "../src/function_mixin.tcc"

#define VOLATILE ;

namespace Ariadne {

typedef unsigned int uint;

typedef Vector<Float> Point;
typedef Vector<Interval> Box;

static double TAYLOR_FUNCTION_WRITING_ACCURACY = 1e-8;

void _set_scaling(ScalarTaylorFunction& x, const Interval& ivl, uint j)
{
    rounding_mode_t rounding_mode=get_rounding_mode();
    set_rounding_mode(upward);
    const Float& l=ivl.lower();
    const Float& u=ivl.upper();
    VOLATILE Float pc=u; pc+=l;
    VOLATILE Float nc=-u; nc-=l;
    VOLATILE Float pg=u; pg-=l;
    VOLATILE Float ng=l; ng-=u;
    x.error()=(pc+nc+pg+ng)/4;
    set_rounding_mode(to_nearest);
    MultiIndex a(x.argument_size());
    x.expansion().append(a,(l+u)/2);
    ++a[j];
    x.expansion().append(a,(l+u)/2);
    set_rounding_mode(rounding_mode);
}

/*
Interval _unscale(const Interval& d, const Interval& x) {
    if(d.lower()==d.upper()) { return Interval(0.0); }
    return ( 2*x - add_ivl(d.lower(),d.upper()) ) / sub_ivl(d.upper(),d.lower());
}

IntervalTaylorModel _unscale(const Interval& d, const IntervalTaylorModel& x) {
    if(d.lower()==d.upper()) { return IntervalTaylorModel(x.argument_size()); }
    return ( 2*x - add_ivl(d.lower(),d.upper()) ) / sub_ivl(d.upper(),d.lower());
}
*/

ScalarTaylorFunction::ScalarTaylorFunction()
    : _domain(0), _model(0)
{ }

ScalarTaylorFunction::ScalarTaylorFunction(const DomainType& d)
    : _domain(d), _model(d.size())
{
}

ScalarTaylorFunction::ScalarTaylorFunction(const DomainType& d, shared_ptr<TaylorModelAccuracy> acc_ptr)
    : _domain(d), _model(d.size(),acc_ptr)
{
}

ScalarTaylorFunction::ScalarTaylorFunction(const DomainType& d, const IntervalTaylorModel& m)
    : _domain(d), _model(m)
{
}

ScalarTaylorFunction::ScalarTaylorFunction(const DomainType& d, const ExpansionType& f, const ErrorType& e)
    : _domain(d), _model(f,e)
{
    ARIADNE_ASSERT_MSG(d.size()==f.argument_size(),"d="<<d<<" f="<<f);
}

ScalarTaylorFunction::ScalarTaylorFunction(const DomainType& d, const RealScalarFunction& f)
    : _domain(d), _model(f.argument_size())
{
    ARIADNE_ASSERT_MSG(d.size()==f.argument_size(),"d="<<d<<" f="<<f);
    Vector<IntervalTaylorModel> x=IntervalTaylorModel::scalings(d);
    this->_model=f.evaluate(x);
    this->_model.clean();
}

ScalarTaylorFunction::ScalarTaylorFunction(const DomainType& d, const IntervalScalarFunction& f)
    : _domain(d), _model(f.argument_size())
{
    ARIADNE_ASSERT_MSG(d.size()==f.argument_size(),"d="<<d<<" f="<<f);
    Vector<IntervalTaylorModel> x=IntervalTaylorModel::scalings(d);
    this->_model=f.evaluate(x);
    this->_model.clean();
}

ScalarTaylorFunction::ScalarTaylorFunction(const DomainType& d, const Polynomial<Float>& p)
    : _domain(d), _model(p.argument_size())
{
    ARIADNE_ASSERT_MSG(d.size()==p.argument_size(),"d="<<d<<" p="<<p);
    Vector<IntervalTaylorModel> x=IntervalTaylorModel::scalings(d);
    this->_model=Ariadne::evaluate(p,x);
    this->_model.clean();
}

ScalarTaylorFunction::ScalarTaylorFunction(const DomainType& d, const Polynomial<Interval>& p)
    : _domain(d), _model(p.argument_size())
{
    ARIADNE_ASSERT_MSG(d.size()==p.argument_size(),"d="<<d<<" p="<<p);
    Vector<IntervalTaylorModel> x=IntervalTaylorModel::scalings(d);
    this->_model=Ariadne::evaluate(p,x);
}


ScalarTaylorFunction ScalarTaylorFunction::constant(const Vector<Interval>& d, const Interval& c)
{
    return ScalarTaylorFunction(d,IntervalTaylorModel::constant(d.size(),c));
}

ScalarTaylorFunction ScalarTaylorFunction::identity(const Interval& ivl)
{
    Vector<Interval> d(1,ivl);
    return ScalarTaylorFunction(d,IntervalTaylorModel::scaling(d.size(),0u,d[0u]));
}

ScalarTaylorFunction ScalarTaylorFunction::coordinate(const Vector<Interval>& d, uint j)
{
    ARIADNE_ASSERT(j<d.size());
    return ScalarTaylorFunction(d,IntervalTaylorModel::scaling(d.size(),j,d[j]));
}

ScalarTaylorFunction ScalarTaylorFunction::variable(const Vector<Interval>& d, uint j)
{
    ARIADNE_DEPRECATED("ScalarTaylorFunction::variable","Use ScalarTaylorFunction::coordinate instead");
    ARIADNE_ASSERT(j<d.size());
    return ScalarTaylorFunction(d,IntervalTaylorModel::scaling(d.size(),j,d[j]));
}


Vector<ScalarTaylorFunction> ScalarTaylorFunction::constants(const Vector<Interval>& d, const Vector<Interval>& c)
{
    ARIADNE_DEPRECATED("ScalarTaylorFunction::constants","Use VectorTaylorFunction::constant instead");
    Vector<ScalarTaylorFunction> x(c.size(),ScalarTaylorFunction(d));
    for(uint i=0; i!=c.size(); ++i) {
        x[i]=c[i];
    }
    return x;
}

Vector<ScalarTaylorFunction> ScalarTaylorFunction::variables(const Vector<Interval>& d)
{
    ARIADNE_DEPRECATED("ScalarTaylorFunction::variables","Use VectorTaylorFunction::identity instead");
    return variables(d,0u,d.size());
}

Vector<ScalarTaylorFunction> ScalarTaylorFunction::variables(const Vector<Interval>& d, uint imin, uint imax)
{
    ARIADNE_DEPRECATED("ScalarTaylorFunction::variables","Use VectorTaylorFunction::projection instead");
    ARIADNE_ASSERT(imin<=imax);
    ARIADNE_ASSERT(imax<=d.size());

    Vector<ScalarTaylorFunction> x(imax-imin);
    for(uint i=imin; i!=imax; ++i) {
        x[i-imin]=ScalarTaylorFunction::variable(d,i);
    }
    return x;
}


ScalarTaylorFunction* ScalarTaylorFunction::_clone() const
{
    return new ScalarTaylorFunction(*this);
}

ScalarTaylorFunction* ScalarTaylorFunction::_create() const
{
    return new ScalarTaylorFunction(this->domain(),this->_model.accuracy_ptr());
}



shared_ptr<TaylorModelAccuracy> ScalarTaylorFunction::accuracy_ptr() const
{
    return this->_model.accuracy_ptr();
}


Polynomial<Interval> polynomial(const IntervalTaylorModel& tm);

Polynomial<Interval>
ScalarTaylorFunction::polynomial() const
{
    Polynomial<Interval> p(this->argument_size());
    p=Ariadne::polynomial(this->model());

    Vector<Polynomial<Interval> > s(this->argument_size());
    for(uint j=0; j!=this->argument_size(); ++j) {
        if(this->domain()[j].radius()<=0) {
            ARIADNE_ASSERT(this->domain()[j].radius()==0);
            s[j]=Polynomial<Interval>::constant(this->argument_size(),this->domain()[j]);
        } else {
            s[j]=Ariadne::polynomial(IntervalTaylorModel::unscaling(this->argument_size(),j,this->domain()[j]));
        }
    }

    return compose(p,s);
}

IntervalScalarFunction
ScalarTaylorFunction::function() const
{
    return new ScalarTaylorFunction(*this);
}

RealScalarFunction
ScalarTaylorFunction::real_function() const
{
    return RealScalarFunction(Polynomial<Real>(this->polynomial()));
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


ScalarTaylorFunction& operator+=(ScalarTaylorFunction& f, const Interval& c) {
    f._model+=c;
    return f;
}

ScalarTaylorFunction& operator-=(ScalarTaylorFunction& f, const Interval& c) {
    f._model-=c;
    return f;
}

ScalarTaylorFunction& operator*=(ScalarTaylorFunction& f, const Interval& c) {
    f._model*=c;
    return f;
}

ScalarTaylorFunction& operator/=(ScalarTaylorFunction& f, const Interval& c) {
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
        ScalarTaylorFunction::DomainType domain=intersection(x1._domain,x2._domain);
        return ScalarTaylorFunction(domain,restrict(x1,domain)._model+restrict(x2,domain)._model);}
}

ScalarTaylorFunction operator-(const ScalarTaylorFunction& x1, const ScalarTaylorFunction& x2) {
    if(x1._domain==x2._domain) {
        return ScalarTaylorFunction(x1._domain,x1._model-x2._model); }
    else {
        ScalarTaylorFunction::DomainType domain=intersection(x1._domain,x2._domain);
        return ScalarTaylorFunction(domain,restrict(x1,domain)._model-restrict(x2,domain)._model);}
}

ScalarTaylorFunction operator*(const ScalarTaylorFunction& x1, const ScalarTaylorFunction& x2) {
    if(x1._domain==x2._domain) {
        return ScalarTaylorFunction(x1._domain,x1._model*x2._model); }
    else {
        ScalarTaylorFunction::DomainType domain=intersection(x1._domain,x2._domain);
        return ScalarTaylorFunction(domain,restrict(x1,domain)._model*restrict(x2,domain)._model);}
}

ScalarTaylorFunction operator/(const ScalarTaylorFunction& x1, const ScalarTaylorFunction& x2) {
    if(x1._domain==x2._domain) {
        return ScalarTaylorFunction(x1._domain,x1._model/x2._model); }
    else {
        ScalarTaylorFunction::DomainType domain=intersection(x1._domain,x2._domain);
        return ScalarTaylorFunction(domain,restrict(x1,domain)._model/restrict(x2,domain)._model);}
}





ScalarTaylorFunction max(const ScalarTaylorFunction& x1, const ScalarTaylorFunction& x2) {
    if(x1._domain==x2._domain) {
        return ScalarTaylorFunction(x1._domain,max(x1._model,x2._model)); }
    else {
        ScalarTaylorFunction::DomainType domain=intersection(x1._domain,x2._domain);
        return ScalarTaylorFunction(domain,max(restrict(x1,domain)._model,restrict(x2,domain)._model));}
}

ScalarTaylorFunction min(const ScalarTaylorFunction& x1, const ScalarTaylorFunction& x2) {
    if(x1._domain==x2._domain) {
        return ScalarTaylorFunction(x1._domain,min(x1._model,x2._model)); }
    else {
        ScalarTaylorFunction::DomainType domain=intersection(x1._domain,x2._domain);
        return ScalarTaylorFunction(domain,min(restrict(x1,domain)._model,restrict(x2,domain)._model));}
}



ScalarTaylorFunction* ScalarTaylorFunction::_derivative(uint j) const
{
    return new ScalarTaylorFunction(Ariadne::derivative(*this,j));
}


Float
ScalarTaylorFunction::evaluate(const Vector<Float>& x) const
{
    return Ariadne::evaluate(this->_model.expansion(),x);
}

Interval
ScalarTaylorFunction::evaluate(const Vector<Interval>& x) const
{
    return Ariadne::evaluate(*this,x);
}

Interval
ScalarTaylorFunction::operator()(const Vector<Interval>& x) const
{
    return Ariadne::evaluate(*this,x);
}


void ScalarTaylorFunction::restrict(const Vector<Interval>& d) {
    ARIADNE_ASSERT(subset(d,this->domain()));
    const Vector<Interval>& od=this->domain();
    for(uint j=0; j!=d.size(); ++j) {
        if(od[j]!=d[j]) { *this=Ariadne::restrict(*this,j,d[j]); }
    }
}

ScalarTaylorFunction restrict(const ScalarTaylorFunction& tv, const Vector<Interval>& d) {
    ARIADNE_ASSERT(subset(d,tv.domain()));
    const Vector<Interval>& od=tv.domain();
    ScalarTaylorFunction r=tv;
    for(uint j=0; j!=d.size(); ++j) {
        if(od[j]!=d[j]) { r=restrict(r,j,d[j]); }
    }
    return r;
}

ScalarTaylorFunction extend(const ScalarTaylorFunction& tv, const Vector<Interval>& d) {
    const Vector<Interval>& domain=tv.domain();
    ARIADNE_ASSERT(subset(domain,d));
    for(uint i=0; i!=d.size(); ++i) {
        ARIADNE_ASSERT(domain[i]==d[i] || domain[i].lower()==domain[i].upper());
    }
    return ScalarTaylorFunction(d,tv._model);
}

Interval
evaluate(const ScalarTaylorFunction& f, const Vector<Interval>& x) {
    if(!subset(x,f.domain())) {
        ARIADNE_THROW(DomainException,"evaluate(f,x) with f="<<f<<", x="<<x,"x is not a subset of f.domain()="<<f.domain());
    }
    return unchecked_evaluate(f,x);
}

Interval
unchecked_evaluate(const ScalarTaylorFunction& f, const Vector<Interval>& x) {
    return evaluate(f.model(),unscale(x,f.domain()));
}


ScalarTaylorFunction
compose(const RealScalarFunction& g, const VectorTaylorFunction& f)
{
    return ScalarTaylorFunction(f.domain(),g.evaluate(f.models()));
}

ScalarTaylorFunction
compose(const IntervalScalarFunctionInterface& g, const VectorTaylorFunction& f)
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
partial_evaluate(const ScalarTaylorFunction& tf, uint k, const Float& c)
{
    return partial_evaluate(tf,k,Interval(c));
}

ScalarTaylorFunction
partial_evaluate(const ScalarTaylorFunction& te, uint k, const Interval& c)
{
    // Scale c to domain
    const uint as=te.argument_size();
    ARIADNE_ASSERT(k<as);
    const Vector<Interval>& domain=te.domain();
    const Interval& dk=domain[k];
    Interval sc=(c-med_ivl(dk))/rad_ivl(dk);

    Vector<Interval> new_domain(as-1);
    for(uint i=0; i!=k; ++i) { new_domain[i]=domain[i]; }
    for(uint i=k; i!=as-1; ++i) { new_domain[i]=domain[i+1]; }

    IntervalTaylorModel new_model=partial_evaluate(te.model(),k,sc);

    return ScalarTaylorFunction(new_domain,new_model);
}


// To scale from a model on [a,b] to a model on [c,d], use scale factor s=(d-c)/(b-a)
// and translation t=((c+d)-(a+b))/(b-a)
// Because we are scaling the model on [-1,+1], this is not the same as
// the mapping taking [a,b] to [c,d]
ScalarTaylorFunction partial_restrict(const ScalarTaylorFunction& tv, uint k, const Interval& new_ivl) {
    ARIADNE_ASSERT(k<tv.argument_size())
    const Interval& old_ivl=tv.domain()[k];
    ARIADNE_ASSERT(subset(new_ivl,old_ivl));
    if(new_ivl==old_ivl) { return tv; }
    Float a=old_ivl.lower(); Float b=old_ivl.upper();
    Float c=new_ivl.lower(); Float d=new_ivl.upper();
    if(a==b) { ARIADNE_ASSERT( a<b || (a==b && c==d) ); return tv; }
    Interval s=sub_ivl(d,c)/sub_ivl(b,a);
    // Interval t=(mul_ivl(b,c)-mul_ivl(a,d))/sub_ivl(b,a);  // WRONG!!
    Interval t=(add_ivl(c,d)-add_ivl(a,b))/sub_ivl(b,a);
    Vector<Interval> new_dom=tv.domain();
    new_dom[k]=new_ivl;
    return ScalarTaylorFunction(new_dom,preaffine(tv.model(),k,s,t));
}


ScalarTaylorFunction restrict(const ScalarTaylorFunction& tv, uint k, const Interval& new_ivl) {
    return partial_restrict(tv,k,new_ivl);
}




pair<ScalarTaylorFunction,ScalarTaylorFunction>
split(const ScalarTaylorFunction& tv, uint j)
{
    pair<IntervalTaylorModel,IntervalTaylorModel> models=split(tv.model(),j);
    pair<Box,Box> subdomains=split(tv.domain(),j);
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
        Vector<Interval> domain=intersection(tv1.domain(),tv2.domain());
        return disjoint(restrict(tv1,domain).model(),restrict(tv2,domain).model());
    }
}

ScalarTaylorFunction intersection(const ScalarTaylorFunction& tv1, const ScalarTaylorFunction& tv2)
{
    ARIADNE_ASSERT(tv1.domain()==tv2.domain());
    return ScalarTaylorFunction(tv1.domain(),intersection(tv1.model(),tv2.model()));
}

Float norm(const ScalarTaylorFunction& f) {
    return norm(f.model());
}

Float distance(const ScalarTaylorFunction& f1, const ScalarTaylorFunction& f2) {
    return norm(f1-f2);
}

Float distance(const ScalarTaylorFunction& f1, const RealScalarFunction& f2) {
    return distance(f1,ScalarTaylorFunction(f1.domain(),f2));
}

Vector<ScalarTaylorFunction> compose(const Vector<ScalarTaylorFunction>& x, const Vector<ScalarTaylorFunction>& y) {
    ARIADNE_NOT_IMPLEMENTED; }
ScalarTaylorFunction compose(const ScalarTaylorFunction& x, const Vector<ScalarTaylorFunction>& y) {
    ARIADNE_NOT_IMPLEMENTED; }

Vector<ScalarTaylorFunction>
prod(const Matrix<Interval>& A,
     const Vector<ScalarTaylorFunction>& x)
{
    ARIADNE_ASSERT(x.size()>0);
    ARIADNE_ASSERT(A.column_size()==x.size());
    for(uint i=1; i!=x.size(); ++i) { ARIADNE_ASSERT(x[i].argument_size()==x[0].argument_size()); }

    Vector<ScalarTaylorFunction> r(A.row_size(),ScalarTaylorFunction(x[0].domain()));
    for(uint i=0; i!=A.row_size(); ++i) {
        for(uint j=0; j!=A.column_size(); ++j) {
            r[i]+=A[i][j]*x[j];
        }
    }
    return r;
}

Matrix<Interval>
jacobian(const Vector<ScalarTaylorFunction>& tv, const Vector<Interval>& x);

ScalarTaylorFunction
midpoint(const ScalarTaylorFunction& f)
{
    IntervalTaylorModel tm=f.model();
    tm.set_error(0.0);
    return ScalarTaylorFunction(f.domain(),tm);
}


ScalarTaylorFunction
compose(const RealScalarFunction& f, const Vector<ScalarTaylorFunction>& g)
{
    ARIADNE_ASSERT(f.argument_size()==g.size());
    for(uint i=0; i!=g.size(); ++i) {
        ARIADNE_ASSERT(g[0].domain()==g[i].domain());
    }

    Vector<Interval> gdomain=g[0].domain();
    Vector<IntervalTaylorModel> gmodels(g.size());
    for(uint i=0; i!=g.size(); ++i) { gmodels[i]=g[i].model(); }

    return ScalarTaylorFunction(gdomain,f.evaluate(gmodels));
}


std::ostream&
ScalarTaylorFunction::write(std::ostream& os) const
{
    os << "ScalarTaylorFunction"<< this->domain()
       << "(" << midpoint(this->polynomial());
    if(this->error()>0.0) { os << "+/-" << this->error(); }
    os <<  ")";
    return os;
}

std::ostream&
operator<<(std::ostream& os, const ScalarTaylorFunction& tf)
{
    return tf.write(os);
}

std::ostream& operator<<(std::ostream& os, const Representation<ScalarTaylorFunction>& frepr)
{
    ScalarTaylorFunction const& function=*frepr.pointer;
    ScalarTaylorFunction truncated_function=function;
    truncated_function.set_error(0.0);
    truncated_function.sweep(TAYLOR_FUNCTION_WRITING_ACCURACY);

    os << midpoint(truncated_function.polynomial());
    if(truncated_function.error()>0.0) { os << "+/-" << truncated_function.error(); }
    if(function.error()>0.0) { os << "+/-" << function.error(); }
    // TODO: Use Unicode +/- literal when this becomes avaialable in C++0x
    return os;
}

/*
std::ostream& operator<<(std::ostream& os, const ModelsRepresentation<ScalarTaylorFunction>& frepr)
{
    ScalarTaylorFunction const& f=*frepr.pointer;
    Float truncatation_error = 0.0;
    os << "<"<<f.domain()<<"\n";
    for(IntervalTaylorModel::const_iterator iter=f.begin(); iter!=f.end(); ++iter) {
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
    tf.sweep(frepr.threshold);
    os << "("<<tf.model()<<"+/-"<<f.error();
    return os;
}

std::ostream& operator<<(std::ostream& os, const PolynomialRepresentation<ScalarTaylorFunction>& frepr)
{
    ScalarTaylorFunction const& function=*frepr.pointer;
    ScalarTaylorFunction truncated_function = function;
    truncated_function.clobber();
    truncated_function.sweep(frepr.threshold);
    Float truncatation_error = truncated_function.error();
    truncated_function.clobber();
    Polynomial<Float> polynomial_function = midpoint(polynomial(truncated_function));
    if(frepr.names.empty()) { os << polynomial_function; }
    else { os << named_argument_repr(polynomial_function,frepr.names); }
    os << "+/-" << truncatation_error << "+/-" << function.error();
    return os;
}





bool
check(const Vector<ScalarTaylorFunction>& tv)
{
    for(uint i=0; i!=tv.size(); ++i) {
        if(tv[0].domain()!=tv[i].domain()) { return false; }
    }
    return true;
}

Vector< Expansion<Float> >
expansion(const Vector<ScalarTaylorFunction>& x)
{
    Vector< Expansion<Float> > r(x.size());
    for(uint i=0; i!=x.size(); ++i) {
        r[i]=x[i].expansion();
    }
    return r;
}

Vector<Float>
error(const Vector<ScalarTaylorFunction>& x)
{
    Vector<Float> r(x.size());
    for(uint i=0; i!=x.size(); ++i) {
        r[i]=x[i].error();
    }
    return r;
}

Vector<Float>
value(const Vector<ScalarTaylorFunction>& x)
{
    Vector<Float> r(x.size());
    for(uint i=0; i!=x.size(); ++i) {
        r[i]=x[i].value();
    }
    return r;
}

Vector<Interval>
ranges(const Vector<ScalarTaylorFunction>& x)
{
    Vector<Interval> r(x.size());
    for(uint i=0; i!=x.size(); ++i) {
        r[i]=x[i].range();
    }
    return r;
}


Vector<Interval>
evaluate(const Vector<ScalarTaylorFunction>& tv, const Vector<Interval>& x)
{
    Vector<Interval> r(tv.size());
    for(uint i=0; i!=tv.size(); ++i) {
        r[i]=evaluate(tv[i],x);
    }
    return r;
}

Matrix<Interval>
jacobian(const Vector<ScalarTaylorFunction>& tv, const Vector<Interval>& x)
{
    ARIADNE_ASSERT(check(tv));
    const Vector<Interval>& dom=tv[0].domain();
    const uint n=dom.size();
    Vector< Differential<Interval> > s(n,n,1u);
    for(uint j=0; j!=n; ++j) {
        Interval dj=dom[j];
        s[j].set_value((x[j]-med_ivl(dj))/rad_ivl(dj));
        s[j].set_gradient(j,rec(rad_ivl(dj)));
    }
    Vector< Expansion<Float> > p=expansion(tv);
    Vector< Differential<Interval> > d=evaluate(p,s);
    //std::cerr<<"  x="<<x<<"\n  p="<<p<<"\n"<<"  s="<<s<<"\n  p.s="<<d<<"\n  J="<<d.jacobian()<<"\n"<<std::endl;
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

VectorTaylorFunction::VectorTaylorFunction(uint m, const Vector<Interval>& d)
    : _domain(d), _models(m,d.size())
{
}

VectorTaylorFunction::VectorTaylorFunction(uint k, const ScalarTaylorFunction& f)
    : _domain(f.domain()), _models(k,f.model())
{
}


VectorTaylorFunction substitute(const VectorTaylorFunction& f, uint k, const ScalarTaylorFunction& h) {
    ARIADNE_ASSERT_MSG(f.argument_size()==h.argument_size()+1u,"f="<<f<<", k="<<k<<", h="<<h);
    const Interval& dk=f.domain()[k];
    IntervalTaylorModel unscaled_model=unscale(h.model(),dk);
    return VectorTaylorFunction(h.domain(),substitute(f.models(),k,unscaled_model));
}

ScalarTaylorFunction substitute(const ScalarTaylorFunction& f, uint k, const ScalarTaylorFunction& h) {
    const Interval& dk=f.domain()[k];
    IntervalTaylorModel unscaled_model=unscale(h.model(),dk);
    return ScalarTaylorFunction(h.domain(),substitute(f.model(),k,unscaled_model));
}


VectorTaylorFunction::VectorTaylorFunction(const Vector<Interval>& d,
                               const Vector<IntervalTaylorModel>& f)
    : _domain(d), _models(f)
{
    for(uint i=0; i!=f.size(); ++i) {
        ARIADNE_ASSERT_MSG(d.size()==f[i].argument_size(),"d="<<d<<"f="<<f);
    }
}

VectorTaylorFunction::VectorTaylorFunction(const Vector<Interval>& d,
                               const Vector< Expansion<Float> >& f)
    : _domain(d), _models(f.size())
{
    for(uint i=0; i!=f.size(); ++i) {
        ARIADNE_ASSERT(d.size()==f[i].argument_size());
        _models[i]=IntervalTaylorModel(f[i],0.0);
    }
}

VectorTaylorFunction::VectorTaylorFunction(const Vector<Interval>& d,
                               const Vector< Expansion<Float> >& f,
                               const Vector<Float>& e)
    : _domain(d), _models(f.size())
{
    ARIADNE_ASSERT(f.size()==e.size());
    for(uint i=0; i!=f.size(); ++i) {
        ARIADNE_ASSERT(d.size()==f[i].argument_size());
        _models[i]=IntervalTaylorModel(f[i],e[i]);
    }
}



VectorTaylorFunction::VectorTaylorFunction(const Vector<Interval>& d,
                                           const IntervalVectorFunction& f)
    : _domain(d), _models(f.result_size())
{
    ARIADNE_ASSERT(f.result_size()>0);
    ARIADNE_ASSERT(d.size()==f.argument_size());
    Vector<IntervalTaylorModel> x=IntervalTaylorModel::scalings(d);
    this->_models=f.evaluate(x);
    this->sweep();
}

VectorTaylorFunction::VectorTaylorFunction(const Vector<Interval>& d,
                               const RealVectorFunction& f)
    : _domain(d), _models(f.result_size())
{
    ARIADNE_ASSERT(f.result_size()>0);
    ARIADNE_ASSERT(d.size()==f.argument_size());
    Vector<IntervalTaylorModel> x=IntervalTaylorModel::scalings(d);
    this->_models=f.evaluate(x);
    this->sweep();
}

VectorTaylorFunction::VectorTaylorFunction(const Vector<Interval>& d,
                               const RealVectorFunction& f,
                               const shared_ptr<IntervalTaylorModel::Accuracy> accuracy_ptr)
    : _domain(d), _models(f.result_size())
{
    ARIADNE_ASSERT(f.result_size()>0);
    ARIADNE_ASSERT(d.size()==f.argument_size());
    Vector<IntervalTaylorModel> x=IntervalTaylorModel::scalings(d);
    for(uint i=0; i!=x.size(); ++i) { x[i].accuracy_ptr()=accuracy_ptr; }
    this->_models=f.evaluate(x);
    this->sweep();
}


VectorTaylorFunction::VectorTaylorFunction(const Vector<Interval>& d,
                               const Vector< Polynomial<Float> >& p)
    : _domain(d), _models(p.size())
{
    for(uint i=0; i!=p.size(); ++i) { ARIADNE_ASSERT(d.size()==p[i].argument_size()); }

    Vector<IntervalTaylorModel> x=IntervalTaylorModel::scalings(d);
    this->_models=Ariadne::evaluate(p,x);
    this->sweep();
}

VectorTaylorFunction::VectorTaylorFunction(const Vector<Interval>& d,
                               const Vector< Polynomial<Interval> >& p)
    : _domain(d), _models(p.size())
{
    for(uint i=0; i!=p.size(); ++i) { ARIADNE_ASSERT(d.size()==p[i].argument_size()); }

    Vector<IntervalTaylorModel> x=IntervalTaylorModel::scalings(d);
    this->_models=Ariadne::evaluate(p,x);
}

VectorTaylorFunction::VectorTaylorFunction(const Vector<ScalarTaylorFunction>& v)
    : _domain(), _models(v.size())
{
    ARIADNE_ASSERT(v.size()>0);
    for(uint i=1; i!=v.size(); ++i) { ARIADNE_ASSERT(v[i].domain()==v[0].domain()); }
    this->_domain=v[0].domain();
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


VectorTaylorFunction* VectorTaylorFunction::_clone() const
{
    return new VectorTaylorFunction(*this);
}

VectorTaylorFunction* VectorTaylorFunction::_create() const
{
    return new VectorTaylorFunction(this->result_size(), ScalarTaylorFunction(this->domain(),this->_models[0].accuracy_ptr()));
}







VectorTaylorFunction
VectorTaylorFunction::constant(const Vector<Interval>& d, const Vector<Float>& c)
{
    return VectorTaylorFunction(d,IntervalTaylorModel::constants(d.size(),c));
}

VectorTaylorFunction
VectorTaylorFunction::constant(const Vector<Interval>& d, const Vector<Interval>& c)
{
    return VectorTaylorFunction(d,IntervalTaylorModel::constants(d.size(),c));
}

VectorTaylorFunction
VectorTaylorFunction::identity(const Vector<Interval>& d)
{
    return VectorTaylorFunction(d,IntervalTaylorModel::scalings(d));
}

VectorTaylorFunction
VectorTaylorFunction::projection(const Vector<Interval>& d, uint imin, uint imax)
{
    return VectorTaylorFunction(ScalarTaylorFunction::variables(d,imin,imax));
}


Polynomial<Interval> polynomial(const IntervalTaylorModel& tm) {
    return Polynomial<Interval>(tm.expansion())+Interval(-tm.error(),+tm.error());
}

Vector< Polynomial<Interval> >
VectorTaylorFunction::polynomial() const
{
    Vector<Polynomial<Interval> > p(this->result_size());
    for(uint i=0; i!=this->result_size(); ++i) {
        p[i]=Ariadne::polynomial(this->models()[i]);
    }

    Vector<Polynomial<Interval> > s(this->argument_size());
    for(uint j=0; j!=this->argument_size(); ++j) {
        if(this->domain()[j].radius()<=0) { ARIADNE_WARN("zero radius in domain of VectorTaylorFunction"<<std::endl); }
        else { s[j]=Ariadne::polynomial(IntervalTaylorModel::unscaling(this->argument_size(),j,this->domain()[j])); }
    }

    return compose(p,s);
}

Vector<Float>
VectorTaylorFunction::errors() const
{
    Vector<Float> e(this->result_size());
    for(uint i=0; i!=this->result_size(); ++i) {
        e[i]=this->models()[i].error();
    }
    return e;
}

Float
VectorTaylorFunction::error() const
{
    Float e=0.0;
    for(uint i=0; i!=this->result_size(); ++i) {
        e=max(e,this->models()[i].error());
    }
    return e;
}

IntervalVectorFunction
VectorTaylorFunction::function() const
{
    return new VectorTaylorFunction(*this);
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



shared_ptr<IntervalTaylorModel::Accuracy>
VectorTaylorFunction::accuracy_ptr() const
{
    return this->_models[0].accuracy_ptr();
}


void
VectorTaylorFunction::set_accuracy(shared_ptr<IntervalTaylorModel::Accuracy> acc_ptr)
{
    for(uint i=0; i!=this->result_size(); ++i) {
        this->_models[i].set_accuracy(acc_ptr);
    }
}

const Vector<Interval>&
VectorTaylorFunction::domain() const
{
    return this->_domain;
}

const Vector<Interval>
VectorTaylorFunction::codomain() const
{
    Vector<Interval> result(this->result_size());
    for(uint i=0; i!=result.size(); ++i) {
        result[i]=this->_models[i].range();
    }
    return result;
}


const Vector<Float>
VectorTaylorFunction::centre() const
{
    Vector<Float> result(this->result_size());
    for(uint i=0; i!=result.size(); ++i) {
        result[i]=this->_models[i].value();
    }
    return result;
}

const Vector<Interval>
VectorTaylorFunction::range() const
{
    Vector<Interval> result(this->result_size());
    for(uint i=0; i!=result.size(); ++i) {
        result[i]=this->_models[i].range();
    }
    return result;
}



const Vector<IntervalTaylorModel>&
VectorTaylorFunction::models() const
{
    return this->_models;
}

const IntervalTaylorModel&
VectorTaylorFunction::model(uint i) const
{
    return this->_models[i];
}

IntervalTaylorModel&
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
VectorTaylorFunction::truncate()
{
    for(uint i=0; i!=this->size(); ++i) {
        this->_models[i].truncate();
    }
    return *this;
}

VectorTaylorFunction&
VectorTaylorFunction::sweep(double threshold)
{
    for(uint i=0; i!=this->size(); ++i) {
        this->_models[i].sweep(threshold);
    }
    return *this;
}

VectorTaylorFunction&
VectorTaylorFunction::truncate(uint degree)
{
    for(uint i=0; i!=this->size(); ++i) {
        this->_models[i].truncate(degree);
    }
    return *this;
}

void
VectorTaylorFunction::set_sweep_threshold(double threshold)
{
    for(uint i=0; i!=this->size(); ++i) {
        this->_models[i].set_sweep_threshold(threshold);
    }
}

void
VectorTaylorFunction::set_maximum_degree(uint degree)
{
    for(uint i=0; i!=this->size(); ++i) {
        this->_models[i].set_maximum_degree(degree);
    }
}

VectorTaylorFunction&
VectorTaylorFunction::clobber()
{
    for(uint i=0; i!=this->size(); ++i) {
        this->_models[i].clobber();
    }
    return *this;
}



Vector<Float>
VectorTaylorFunction::evaluate(const Vector<Float>& x) const
{
    const VectorTaylorFunction& f=*this;
    if(!contains(f.domain(),x)) {
        ARIADNE_THROW(DomainException,"f.evaluate(x) with f="<<f<<", x="<<x,"x is not an element of f.domain()="<<f.domain());
    }
    Vector<Float> sx=Ariadne::unscale(x,f._domain);
    Vector<Float> r(this->result_size());
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=Ariadne::evaluate(this->_models[i].expansion(),sx);
    }
    return r;
}


Vector<Interval>
VectorTaylorFunction::evaluate(const Vector<Interval>& x) const
{
    const VectorTaylorFunction& f=*this;
    if(!subset(x,f.domain())) {
        ARIADNE_THROW(DomainException,"f.evaluate(x) with f="<<f<<", x="<<x,"x is not a subset of f.domain()="<<f.domain());
    }
    Vector<Interval> sx=Ariadne::evaluate(IntervalTaylorModel::unscalings(f._domain),x);
    return Ariadne::evaluate(f._models,sx);
}

Vector<Interval>
VectorTaylorFunction::operator()(const Vector<Interval>& x) const
{
    return this->evaluate(x);
}

Matrix<Interval>
VectorTaylorFunction::jacobian(const Vector<Interval>& x) const
{
    Matrix<Interval> J=Ariadne::jacobian(this->_models,unscale(x,this->_domain));
    for(uint j=0; j!=J.column_size(); ++j) {
        Interval rad=rad_ivl(this->_domain[j]);
        for(uint i=0; i!=J.row_size(); ++i) {
            J[i][j]/=rad;
        }
    }
    return J;
}


VectorTaylorFunction
join(const VectorTaylorFunction& f1, const ScalarTaylorFunction& f2)
{
    ARIADNE_ASSERT(f1.domain()==f2.domain());
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
    return VectorTaylorFunction(join(f1.domain(),f2.domain()),combine(f1.model(),f2.model()));
}

VectorTaylorFunction
combine(const ScalarTaylorFunction& f1, const VectorTaylorFunction& f2)
{
    return VectorTaylorFunction(join(f1.domain(),f2.domain()),combine(f1.model(),f2.models()));
}

VectorTaylorFunction
combine(const VectorTaylorFunction& f1, const ScalarTaylorFunction& f2)
{
    return VectorTaylorFunction(join(f1.domain(),f2.domain()),combine(f1.models(),f2.model()));
}

VectorTaylorFunction
combine(const VectorTaylorFunction& f1, const VectorTaylorFunction& f2)
{
    return VectorTaylorFunction(join(f1.domain(),f2.domain()),combine(f1.models(),f2.models()));
}


VectorTaylorFunction
embed(const VectorTaylorFunction& f, const Interval& d)
{
    return VectorTaylorFunction(join(f.domain(),d),embed(f.models(),1u));
}

VectorTaylorFunction
embed(const VectorTaylorFunction& f, const Vector<Interval>& d)
{
    return VectorTaylorFunction(join(f.domain(),d),embed(f.models(),d.size()));
}

VectorTaylorFunction
embed(const Vector<Interval>& d, const VectorTaylorFunction& f)
{
    return VectorTaylorFunction(join(d,f.domain()),embed(d.size(),f.models()));
}

VectorTaylorFunction
restrict(const VectorTaylorFunction& f, const Vector<Interval>& d)
{
    ARIADNE_ASSERT_MSG(subset(d,f.domain()),"Cannot restrict "<<f<<" to non-sub-domain "<<d);
    if(d==f.domain()) { return f; }
    VectorTaylorFunction r(f.result_size(),d);
    for(uint i=0; i!=r.result_size(); ++i) {
        r.set(i,restrict(f[i],d));
    }
    r.set_accuracy(f.accuracy_ptr());
    return r;
}

std::pair<VectorTaylorFunction,VectorTaylorFunction>
split(const VectorTaylorFunction& tf, uint j)
{
    pair< Vector<IntervalTaylorModel>,Vector<IntervalTaylorModel> > models=split(tf.models(),j);
    pair<Box,Box> subdomains=split(tf.domain(),j);
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
operator+=(VectorTaylorFunction& f, const Vector<Interval>& e)
{
    ARIADNE_ASSERT(f.result_size()==e.size());
    f._models+=e;
    return f;
}

VectorTaylorFunction&
operator-=(VectorTaylorFunction& f, const Vector<Interval>& e)
{
    ARIADNE_ASSERT(f.result_size()==e.size());
    f._models-=e;
    return f;
}

VectorTaylorFunction&
operator*=(VectorTaylorFunction& f, const Float& c)
{
    f._models*=c;
    return f;
}

VectorTaylorFunction&
operator/=(VectorTaylorFunction& f, const Float& c)
{
    f._models/=c;
    return f;
}


VectorTaylorFunction
operator+(const VectorTaylorFunction& f1, const VectorTaylorFunction& f2)
{
    ARIADNE_ASSERT_MSG(!intersection(f1.domain(),f2.domain()).empty(),
                       "operator+(VectorTaylorFunction f1, VectorTaylorFunction f2) with f1="<<f1<<" f2="<<f2<<
                       ": domains are disjoint");
    if(f1.domain()==f2.domain()) {
        return VectorTaylorFunction(f1.domain(),Vector<IntervalTaylorModel>(f1.models()+f2.models()));
    } else {
        Box new_domain=intersection(f1.domain(),f2.domain());
        return operator+(restrict(f1,new_domain),restrict(f2,new_domain));
    }
}


VectorTaylorFunction
operator-(const VectorTaylorFunction& f1, const VectorTaylorFunction& f2)
{
    ARIADNE_ASSERT(!intersection(f1.domain(),f2.domain()).empty());
    if(f1.domain()==f2.domain()) {
        return VectorTaylorFunction(f1.domain(),Vector<IntervalTaylorModel>(f1.models()-f2.models()));
    } else {
        Box new_domain=intersection(f1.domain(),f2.domain());
        return operator-(restrict(f1,new_domain),restrict(f2,new_domain));
    }
}



VectorTaylorFunction
operator-(const VectorTaylorFunction& f)
{
    return VectorTaylorFunction(f.domain(),Vector<IntervalTaylorModel>(-f.models()));
}

VectorTaylorFunction
operator*(const Float& c, const VectorTaylorFunction& f)
{
    return VectorTaylorFunction(f.domain(),Vector<IntervalTaylorModel>(f.models()*c));
}

VectorTaylorFunction
operator*(const VectorTaylorFunction& f, const Float& c)
{
    return VectorTaylorFunction(f.domain(),Vector<IntervalTaylorModel>(f.models()*c));
}

VectorTaylorFunction
operator/(const VectorTaylorFunction& f, const Float& c)
{
    return VectorTaylorFunction(f.domain(),Vector<IntervalTaylorModel>(f.models()/c));
}

VectorTaylorFunction
operator+(const VectorTaylorFunction& f, const Vector<Float>& c)
{
    return VectorTaylorFunction(f.domain(),Vector<IntervalTaylorModel>(f.models()+c));
}

VectorTaylorFunction
operator+(const VectorTaylorFunction& f, const Vector<Interval>& c)
{
    return VectorTaylorFunction(f.domain(),Vector<IntervalTaylorModel>(f.models()+c));
}

VectorTaylorFunction
operator-(const VectorTaylorFunction& f, const Vector<Float>& c)
{
    return VectorTaylorFunction(f.domain(),Vector<IntervalTaylorModel>(f.models()-c));
}

VectorTaylorFunction
operator-(const VectorTaylorFunction& f, const Vector<Interval>& c)
{
    return VectorTaylorFunction(f.domain(),Vector<IntervalTaylorModel>(f.models()-c));
}

VectorTaylorFunction
operator*(const Matrix<Float>& A, const VectorTaylorFunction& f)
{
    return VectorTaylorFunction(f.domain(),Vector<IntervalTaylorModel>(prod(A,f.models())));
}

VectorTaylorFunction
operator*(const Matrix<Interval>& A, const VectorTaylorFunction& f)
{
    return VectorTaylorFunction(f.domain(),Vector<IntervalTaylorModel>(boost::numeric::ublas::prod(A,f.models())));
}





VectorTaylorFunction
partial_evaluate(const VectorTaylorFunction& tf, uint k, const Float& c)
{
    return partial_evaluate(tf,k,Interval(c));
}

VectorTaylorFunction
partial_evaluate(const VectorTaylorFunction& tf, uint k, const Interval& c)
{
    // Scale c to domain
    const uint as=tf.argument_size();
    ARIADNE_ASSERT(k<as);
    const Vector<Interval>& domain=tf.domain();
    const Interval& dk=domain[k];
    Interval sc=(c-med_ivl(dk))/rad_ivl(dk);

    Vector<Interval> new_domain(as-1);
    for(uint i=0; i!=k; ++i) { new_domain[i]=domain[i]; }
    for(uint i=k; i!=as-1; ++i) { new_domain[i]=domain[i+1]; }

    Vector<IntervalTaylorModel> new_models=partial_evaluate(tf.models(),k,sc);

    return VectorTaylorFunction(new_domain,new_models);
}


VectorTaylorFunction
partial_restrict(const VectorTaylorFunction& tf, uint k, const Interval& c)
{
    VectorTaylorFunction r(tf.result_size(),tf.domain());
    for(uint i=0; i!=tf.result_size(); ++i) {
        r[i]=partial_restrict(tf[i],k,c);
    }
    return r;
}

VectorTaylorFunction
restrict(const VectorTaylorFunction& tf, uint k, const Interval& c)
{
    return partial_restrict(tf,k,c);
}


Vector<Interval>
evaluate(const VectorTaylorFunction& f, const Vector<Interval>& x) {
    if(!subset(x,f.domain())) {
        ARIADNE_THROW(DomainException,"evaluate(f,x) with f="<<f<<", x="<<x,"x is not a subset of f.domain()="<<f.domain());
    }
    return unchecked_evaluate(f,x);
}

Vector<Interval>
unchecked_evaluate(const VectorTaylorFunction& f, const Vector<Interval>& x) {
    return evaluate(f.models(),unscale(x,f.domain()));
}

VectorTaylorFunction
compose(const RealVectorFunction& g, const VectorTaylorFunction& f)
{
    return VectorTaylorFunction(f.domain(),g.evaluate(f.models()));
}

VectorTaylorFunction
compose(const IntervalVectorFunctionInterface& g, const VectorTaylorFunction& f)
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
antiderivative(const VectorTaylorFunction& f, uint k)
{
    ARIADNE_ASSERT_MSG(k<f.argument_size(),"f="<<f<<"\n f.argument_size()="<<f.argument_size()<<" k="<<k);
    Interval fdomkrad=rad_ivl(f.domain()[k].lower(),f.domain()[k].upper());
    VectorTaylorFunction g=f;
    for(uint i=0; i!=g.size(); ++i) {
        g._models[i].antidifferentiate(k);
        g._models[i]*=fdomkrad;
    }
    return g;
}




template<class X> inline Vector<X> join(const Vector<X>& v1, const Vector<X>& v2, const X& s3) {
    Vector<X> r(v1.size()+v2.size()+1u);
    for(uint i=0; i!=v1.size(); ++i) { r[i]=v1[i]; }
    for(uint i=0; i!=v2.size(); ++i) { r[v1.size()+i]=v2[i]; }
    r[v1.size()+v2.size()]=s3;
    return r;
}



Float norm(const VectorTaylorFunction& f) {
    Float res=0.0;
    for(uint i=0; i!=f.result_size(); ++i) {
        res=max(res,norm(f[i]));
    }
    return res;
}

Float distance(const VectorTaylorFunction& f1, const VectorTaylorFunction& f2) {
    return norm(f1-f2);
}

Float distance(const VectorTaylorFunction& f1, const RealVectorFunction& f2) {
    return distance(f1,VectorTaylorFunction(f1.domain(),f2));
}


std::ostream&
VectorTaylorFunction::write(std::ostream& os) const
{
    //os << "VectorTaylorFunction( "<<this->domain()<<" , ";
    //for(uint i=0; i!=this->result_size(); ++i) {
    //    os << (i==0?'[':',')<<this->_models[i].expansion()<<","<<this->_models[i].error();
    //}
    //return os << "] )";
    return os << "VectorTaylorFunction(d=" << this->domain() << ", p~" << midpoint(this->polynomial()) << ", e=" << this->errors() << ", m=" << this->models() << ")";
}

std::ostream& operator<<(std::ostream& os, const Representation<VectorTaylorFunction>& repr)
{
    const VectorTaylorFunction& function = *repr.pointer;
    os << "[";
    for(uint i=0; i!=function.result_size(); ++i) {
        if(i!=0) { os << ","; }
        os << Ariadne::repr(function[i]);
    }
    return os << "]";
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

Polynomial<Interval> polynomial(const ScalarTaylorFunction& tfn) {
    return Polynomial<Interval>(tfn.polynomial());
}

Vector< Polynomial<Interval> > polynomial(const VectorTaylorFunction& tfn) {
    return Vector< Polynomial<Interval> >(tfn.polynomial());
}

List< Polynomial<Interval> > polynomials(const List<ScalarTaylorFunction>& tfns) {
    List< Polynomial<Interval> > result;
    for(uint i=0; i!=tfns.size(); ++i) {
        result.append(polynomial(tfns[i]));
    }
    return result;
}


ScalarTaylorFunction
TaylorFunctionFactory::create_zero(const IntervalVector& domain) const
{
    return ScalarTaylorFunction(domain,IntervalTaylorModel(domain.size(),this->_acc_ptr));
}

ScalarTaylorFunction*
TaylorFunctionFactory::_create_zero(const IntervalVector& domain) const
{
    return new ScalarTaylorFunction(domain,IntervalTaylorModel(domain.size(),this->_acc_ptr));
}

VectorTaylorFunction
TaylorFunctionFactory::create_identity(const IntervalVector& domain) const
{
    VectorTaylorFunction result(domain.size(),domain);
    for(uint i=0; i!=domain.size(); ++i) {
        result._models[i]=ScalarTaylorFunction::coordinate(domain,i)._model;
        result._models[i].set_accuracy(this->_acc_ptr);
    }
    return result;
}

VectorTaylorFunction*
TaylorFunctionFactory::_create_identity(const IntervalVector& domain) const
{
    VectorTaylorFunction* result=new VectorTaylorFunction(domain.size(),domain);
    for(uint i=0; i!=domain.size(); ++i) {
        result->_models[i]=ScalarTaylorFunction::coordinate(domain,i)._model;
        result->_models[i].set_accuracy(this->_acc_ptr);
    }
    return result;
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
            const Interval& a=p.centre_derivatives()[i][j];
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
