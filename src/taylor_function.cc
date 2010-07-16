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

namespace Ariadne {

typedef unsigned int uint;

typedef Vector<Float> Point;
typedef Vector<Interval> Box;


void _set_scaling(ScalarTaylorFunction& x, const Interval& ivl, uint j)
{
    rounding_mode_t rounding_mode=get_rounding_mode();
    set_rounding_mode(upward);
    const double& l=ivl.lower();
    const double& u=ivl.upper();
    volatile double pc=u; pc+=l;
    volatile double nc=-u; nc-=l;
    volatile double pg=u; pg-=l;
    volatile double ng=l; ng-=u;
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

TaylorModel _unscale(const Interval& d, const TaylorModel& x) {
    if(d.lower()==d.upper()) { return TaylorModel(x.argument_size()); }
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

ScalarTaylorFunction::ScalarTaylorFunction(const DomainType& d, const TaylorModel& m)
    : _domain(d), _model(m)
{
}

ScalarTaylorFunction::ScalarTaylorFunction(const DomainType& d, const ExpansionType& f, const ErrorType& e)
    : _domain(d), _model(f,e)
{
    ARIADNE_ASSERT_MSG(d.size()==f.argument_size(),"d="<<d<<" f="<<f);
}

ScalarTaylorFunction::ScalarTaylorFunction(const DomainType& d, const ScalarFunction& f)
    : _domain(d), _model(f.argument_size())
{
    ARIADNE_ASSERT_MSG(d.size()==f.argument_size(),"d="<<d<<" f="<<f);
    Vector<TaylorModel> x=TaylorModel::scalings(d);
    this->_model=f.evaluate(x);
}

ScalarTaylorFunction::ScalarTaylorFunction(const DomainType& d, const Polynomial<Float>& p)
    : _domain(d), _model(p.argument_size())
{
    ARIADNE_ASSERT_MSG(d.size()==p.argument_size(),"d="<<d<<" p="<<p);
    Vector<TaylorModel> x=TaylorModel::scalings(d);
    this->_model=Ariadne::evaluate(p,x);
}

ScalarTaylorFunction::ScalarTaylorFunction(const DomainType& d, const Polynomial<Interval>& p)
    : _domain(d), _model(p.argument_size())
{
    ARIADNE_ASSERT_MSG(d.size()==p.argument_size(),"d="<<d<<" p="<<p);
    Vector<TaylorModel> x=TaylorModel::scalings(d);
    this->_model=Ariadne::evaluate(p,x);
}


ScalarTaylorFunction ScalarTaylorFunction::constant(const Vector<Interval>& d, const Float& c)
{
    return ScalarTaylorFunction(d,TaylorModel::constant(d.size(),c));
}

ScalarTaylorFunction ScalarTaylorFunction::constant(const Vector<Interval>& d, const Interval& c)
{
    return ScalarTaylorFunction(d,TaylorModel::constant(d.size(),c));
}

ScalarTaylorFunction ScalarTaylorFunction::identity(const Interval& ivl)
{
    Vector<Interval> d(1,ivl);
    return ScalarTaylorFunction(d,TaylorModel::scaling(d.size(),0u,d[0u]));
}

ScalarTaylorFunction ScalarTaylorFunction::coordinate(const Vector<Interval>& d, uint j)
{
    ARIADNE_ASSERT(j<d.size());
    return ScalarTaylorFunction(d,TaylorModel::scaling(d.size(),j,d[j]));
}

ScalarTaylorFunction ScalarTaylorFunction::variable(const Vector<Interval>& d, uint j)
{
    ARIADNE_DEPRECATED("ScalarTaylorFunction::variable","Use ScalarTaylorFunction::coordinate instead");
    ARIADNE_ASSERT(j<d.size());
    return ScalarTaylorFunction(d,TaylorModel::scaling(d.size(),j,d[j]));
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


Polynomial<Interval> polynomial(const TaylorModel& tm);

Polynomial<Interval>
ScalarTaylorFunction::polynomial() const
{
    Polynomial<Interval> p(this->argument_size());
    p=Ariadne::polynomial(this->model());

    Vector<Polynomial<Interval> > s(this->argument_size());
    for(uint j=0; j!=this->argument_size(); ++j) {
        if(this->domain()[j].radius()<=0) {
            s[j]=Polynomial<Interval>::constant(this->argument_size(),this->domain()[j].lower());
        } else {
            s[j]=Ariadne::polynomial(TaylorModel::unscaling(this->argument_size(),j,this->domain()[j]));
        }
    }

    return compose(p,s);
}

ScalarFunction
ScalarTaylorFunction::function() const
{
    VectorFunction unscaling=VectorUnscalingFunction(this->domain());
    ScalarFunction polynomial=ScalarFunction(this->expansion());
    ScalarFunction error=ScalarFunction::constant(this->argument_size(),Real(Interval(-this->error(),+this->error())));
    return compose(polynomial+error,unscaling);
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





Interval
ScalarTaylorFunction::evaluate(const Vector<Float>& x) const
{
    return Ariadne::evaluate(*this,Vector<Interval>(x));
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
compose(const ScalarFunction& g, const VectorTaylorFunction& f)
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

    TaylorModel new_model=partial_evaluate(te.model(),k,sc);

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
    pair<TaylorModel,TaylorModel> models=split(tv.model(),j);
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

Float distance(const ScalarTaylorFunction& f1, const ScalarFunction& f2) {
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
    TaylorModel tm=f.model();
    tm.set_error(0.0);
    return ScalarTaylorFunction(f.domain(),tm);
}


ScalarTaylorFunction
compose(const ScalarFunction& f, const Vector<ScalarTaylorFunction>& g)
{
    ARIADNE_ASSERT(f.argument_size()==g.size());
    for(uint i=0; i!=g.size(); ++i) {
        ARIADNE_ASSERT(g[0].domain()==g[i].domain());
    }

    Vector<Interval> gdomain=g[0].domain();
    Vector<TaylorModel> gmodels(g.size());
    for(uint i=0; i!=g.size(); ++i) { gmodels[i]=g[i].model(); }

    return ScalarTaylorFunction(gdomain,f.evaluate(gmodels));
}


Vector<ScalarTaylorFunction>
implicit(const Vector<ScalarTaylorFunction>& f)
{
    ARIADNE_ASSERT(f.size()>0);
    ARIADNE_ASSERT(f.size()<=f[0].argument_size());
    // Solve the equation f(x,h(x))=0
    // Use D1f + D2f Dh = 0, so Dh=-D2f^-1 D1f
    uint rs=f.size();
    uint fas=f[0].argument_size();
    uint has=fas-rs;

    Vector<Interval> h_domain=project(f[0].domain(),range(0u,has));
    Vector<Interval> h_range=project(f[0].domain(),range(has,fas));
    Vector<ScalarTaylorFunction> id=ScalarTaylorFunction::variables(h_domain);
    Vector<ScalarTaylorFunction> h=ScalarTaylorFunction::constants(h_domain,h_range);

    for(uint k=0; k!=10; ++k) {
        Vector<Interval> ih_range=join(h_domain,h_range);
        Matrix<Interval> Df=jacobian(f,ih_range);
        //std::cerr<<"  Df="<<Df<<std::endl;
        Matrix<Interval> D2f=project(Df,range(0,rs),range(has,fas));
        //std::cerr<<"  D2f="<<J<<std::endl;
        Matrix<Interval> D2finv=inverse(D2f);

        for(uint i=0; i!=rs; ++i) {
            h[i].set_error(0);
        }
        Vector<ScalarTaylorFunction> idh=join(id,h);
        Vector<ScalarTaylorFunction> fidxhx=compose(f,idh);
        //std::cerr<<"  f(x,h(x))="<<fh<<std::endl;
        Vector<ScalarTaylorFunction> dh=prod(D2finv,fidxhx);
        //std::cerr<<"  dh="<<dh<<std::endl;
        h=h-dh;
    }
    //std::cerr<<"\n  f="<<f<<"\n  h[0]="<<h0<<"\n  h[1]="<<h1<<"\n\n";
    ARIADNE_ASSERT(h.size()==f.size());
    ARIADNE_ASSERT(h[0].argument_size()+h.size()==f[0].argument_size());
    return h;

}

ScalarTaylorFunction implicit(const ScalarTaylorFunction& f) {
    Vector<Interval> h_domain=project(f.domain(),range(0u,f.argument_size()-1u));
    Interval h_codomain=f.domain()[f.argument_size()-1u];
    TaylorModel h_model=implicit(f.model());
    ARIADNE_ASSERT(h_model.argument_size()+1==f.model().argument_size());
    TaylorModel hrs_model=h_model.rescale(Interval(-1,+1),h_codomain);
    ARIADNE_ASSERT(hrs_model.argument_size()+1==f.model().argument_size());
    return ScalarTaylorFunction(h_domain,hrs_model);
}


TaylorModel implicit(const ScalarFunction& f, const Vector<TaylorModel>& g);

ScalarTaylorFunction
implicit(const ScalarFunction& f, const Vector<ScalarTaylorFunction>& g)
{
    ARIADNE_ASSERT(f.argument_size()>g.size());
    ARIADNE_ASSERT(g.size()>0);
    for(uint i=1; i!=g.size(); ++i) {
        ARIADNE_ASSERT(g[i].domain()==g[0].domain());
    }

    Vector<Interval> domain=g[0].domain();
    Vector<TaylorModel> models(g.size());
    for(uint i=0; i!=g.size(); ++i) {
        models[i]=g[i].model();
    }

    return ScalarTaylorFunction(domain,implicit(f,models));
}

ScalarTaylorFunction
implicit(const ScalarFunction& f, const Vector<Interval>& d)
{
    ARIADNE_ASSERT(f.argument_size()>=1u);
    ARIADNE_ASSERT(d.size()+1u==f.argument_size());

    Vector<ScalarTaylorFunction> id=ScalarTaylorFunction::variables(d);
    //for(uint i=0; i!=d.size(); ++i) { id[i].model().set_accuracy(a); }
    return implicit(f,id);
}


ScalarTaylorFunction
crossing_time(const ScalarFunction& g, const VectorTaylorFunction& phi)
{
    Vector<Interval> d=project(phi.domain(),range(0,phi.result_size()));
    Interval h=phi.domain()[phi.result_size()];

    ScalarTaylorFunction gphi(phi.domain(),g.evaluate(phi.models()));

    return implicit(gphi);
}


std::ostream&
operator<<(std::ostream& os, const ScalarTaylorFunction& tv) {
    return os << "ScalarTaylorFunction( domain=" << tv.domain() << ", polynomial=" << midpoint(tv.polynomial()) << "+/-" << tv.error() << ", model=" << tv.model() << ")";
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
        s[j].set_value((x[j]-dj.midpoint())/dj.radius());
        s[j].set_gradient(j,1/dj.radius());
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
    TaylorModel unscaled_model=unscale(h.model(),dk);
    return VectorTaylorFunction(h.domain(),substitute(f.models(),k,unscaled_model));
}

ScalarTaylorFunction substitute(const ScalarTaylorFunction& f, uint k, const ScalarTaylorFunction& h) {
    const Interval& dk=f.domain()[k];
    TaylorModel unscaled_model=unscale(h.model(),dk);
    return ScalarTaylorFunction(h.domain(),substitute(f.model(),k,unscaled_model));
}


VectorTaylorFunction::VectorTaylorFunction(const Vector<Interval>& d,
                               const Vector<TaylorModel>& f)
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
        _models[i]=TaylorModel(f[i],0.0);
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
        _models[i]=TaylorModel(f[i],e[i]);
    }
}



VectorTaylorFunction::VectorTaylorFunction(const Vector<Interval>& d,
                               const VectorFunction& f)
    : _domain(d), _models(f.result_size())
{
    ARIADNE_ASSERT(f.result_size()>0);
    ARIADNE_ASSERT(d.size()==f.argument_size());
    Vector<TaylorModel> x=TaylorModel::scalings(d);
    this->_models=f.evaluate(x);
}

VectorTaylorFunction::VectorTaylorFunction(const Vector<Interval>& d,
                               const VectorFunction& f,
                               const shared_ptr<TaylorModel::Accuracy> accuracy_ptr)
    : _domain(d), _models(f.result_size())
{
    ARIADNE_ASSERT(f.result_size()>0);
    ARIADNE_ASSERT(d.size()==f.argument_size());
    Vector<TaylorModel> x=TaylorModel::scalings(d);
    for(uint i=0; i!=x.size(); ++i) { x[i].accuracy_ptr()=accuracy_ptr; }
    this->_models=f.evaluate(x);
}


VectorTaylorFunction::VectorTaylorFunction(const Vector<Interval>& d,
                               const Vector< Polynomial<Float> >& p)
    : _domain(d), _models(p.size())
{
    for(uint i=0; i!=p.size(); ++i) { ARIADNE_ASSERT(d.size()==p[i].argument_size()); }

    Vector<TaylorModel> x=TaylorModel::scalings(d);
    this->_models=Ariadne::evaluate(p,x);
}

VectorTaylorFunction::VectorTaylorFunction(const Vector<Interval>& d,
                               const Vector< Polynomial<Interval> >& p)
    : _domain(d), _models(p.size())
{
    for(uint i=0; i!=p.size(); ++i) { ARIADNE_ASSERT(d.size()==p[i].argument_size()); }

    Vector<TaylorModel> x=TaylorModel::scalings(d);
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










VectorTaylorFunction
VectorTaylorFunction::constant(const Vector<Interval>& d, const Vector<Float>& c)
{
    return VectorTaylorFunction(d,TaylorModel::constants(d.size(),c));
}

VectorTaylorFunction
VectorTaylorFunction::constant(const Vector<Interval>& d, const Vector<Interval>& c)
{
    return VectorTaylorFunction(d,TaylorModel::constants(d.size(),c));
}

VectorTaylorFunction
VectorTaylorFunction::identity(const Vector<Interval>& d)
{
    return VectorTaylorFunction(d,TaylorModel::scalings(d));
}

VectorTaylorFunction
VectorTaylorFunction::projection(const Vector<Interval>& d, uint imin, uint imax)
{
    return VectorTaylorFunction(ScalarTaylorFunction::variables(d,imin,imax));
}


Polynomial<Interval> polynomial(const TaylorModel& tm) {
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
        if(this->domain()[j].radius()<=0) { std::cerr<<"WARNING: zero radius in domain of VectorTaylorFunction"<<std::endl; }
        else { s[j]=Ariadne::polynomial(TaylorModel::unscaling(this->argument_size(),j,this->domain()[j])); }
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

VectorFunction
VectorTaylorFunction::function() const
{
    VectorFunction unscaling=VectorUnscalingFunction(this->domain());
    VectorFunction polynomials(this->result_size(),this->argument_size());
    for(uint i=0; i!=this->result_size(); ++i) {
        polynomials[i]=ScalarFunction(this->_models[i].expansion())+Real(Interval(-this->_models[i].error(),+this->_models[i].error()));
    }
    return compose(polynomials,unscaling);

    Vector< Polynomial<Interval> > polynomial=this->polynomial();
    VectorFunction result(this->result_size(),this->argument_size());
    for(uint i=0; i!=this->result_size(); ++i) {
        result.set(i,ScalarFunction(polynomial[i]));
    }
    return result;
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



shared_ptr<TaylorModel::Accuracy>
VectorTaylorFunction::accuracy_ptr() const
{
    return this->_models[0].accuracy_ptr();
}


void
VectorTaylorFunction::set_accuracy(shared_ptr<TaylorModel::Accuracy> acc_ptr)
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



const Vector<TaylorModel>&
VectorTaylorFunction::models() const
{
    return this->_models;
}

const TaylorModel&
VectorTaylorFunction::model(uint i) const
{
    return this->_models[i];
}

TaylorModel&
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



Vector<Interval>
VectorTaylorFunction::evaluate(const Vector<Float>& x) const
{
    return this->evaluate(Vector<Interval>(x));
}


Vector<Interval>
VectorTaylorFunction::evaluate(const Vector<Interval>& x) const
{
    const VectorTaylorFunction& f=*this;
    if(!subset(x,f.domain())) {
        ARIADNE_THROW(DomainException,"f.evaluate(x) with f="<<f<<", x="<<x,"x is not a subset of f.domain()="<<f.domain());
    }
    Vector<Interval> sx=Ariadne::evaluate(TaylorModel::unscalings(f._domain),x);
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
    pair< Vector<TaylorModel>,Vector<TaylorModel> > models=split(tf.models(),j);
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
        return VectorTaylorFunction(f1.domain(),Vector<TaylorModel>(f1.models()+f2.models()));
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
        return VectorTaylorFunction(f1.domain(),Vector<TaylorModel>(f1.models()-f2.models()));
    } else {
        Box new_domain=intersection(f1.domain(),f2.domain());
        return operator-(restrict(f1,new_domain),restrict(f2,new_domain));
    }
}



VectorTaylorFunction
operator-(const VectorTaylorFunction& f)
{
    return VectorTaylorFunction(f.domain(),Vector<TaylorModel>(-f.models()));
}

VectorTaylorFunction
operator*(const Float& c, const VectorTaylorFunction& f)
{
    return VectorTaylorFunction(f.domain(),Vector<TaylorModel>(f.models()*c));
}

VectorTaylorFunction
operator*(const VectorTaylorFunction& f, const Float& c)
{
    return VectorTaylorFunction(f.domain(),Vector<TaylorModel>(f.models()*c));
}

VectorTaylorFunction
operator/(const VectorTaylorFunction& f, const Float& c)
{
    return VectorTaylorFunction(f.domain(),Vector<TaylorModel>(f.models()/c));
}

VectorTaylorFunction
operator+(const VectorTaylorFunction& f, const Vector<Float>& c)
{
    return VectorTaylorFunction(f.domain(),Vector<TaylorModel>(f.models()+c));
}

VectorTaylorFunction
operator+(const VectorTaylorFunction& f, const Vector<Interval>& c)
{
    return VectorTaylorFunction(f.domain(),Vector<TaylorModel>(f.models()+c));
}

VectorTaylorFunction
operator-(const VectorTaylorFunction& f, const Vector<Float>& c)
{
    return VectorTaylorFunction(f.domain(),Vector<TaylorModel>(f.models()-c));
}

VectorTaylorFunction
operator-(const VectorTaylorFunction& f, const Vector<Interval>& c)
{
    return VectorTaylorFunction(f.domain(),Vector<TaylorModel>(f.models()-c));
}

VectorTaylorFunction
operator*(const Matrix<Float>& A, const VectorTaylorFunction& f)
{
    return VectorTaylorFunction(f.domain(),Vector<TaylorModel>(prod(A,f.models())));
}

VectorTaylorFunction
operator*(const Matrix<Interval>& A, const VectorTaylorFunction& f)
{
    return VectorTaylorFunction(f.domain(),Vector<TaylorModel>(boost::numeric::ublas::prod(A,f.models())));
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

    Vector<TaylorModel> new_models=partial_evaluate(tf.models(),k,sc);

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
compose(const VectorFunction& g, const VectorTaylorFunction& f)
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

ScalarTaylorFunction
implicit(const ScalarFunction& f, const VectorTaylorFunction& g)
{
    return ScalarTaylorFunction(g.domain(),implicit(f,g.models()));
}

VectorTaylorFunction
implicit(const VectorTaylorFunction& f)
{
    ARIADNE_ASSERT_MSG(f.argument_size()>f.result_size(),"f.argument_size()<=f.result_size() in implicit(f): f="<<f);
    uint fas=f.argument_size();
    uint has=f.argument_size()-f.result_size();
    Vector<Interval> hdom=project(f.domain(),range(0,has));
    Vector<Interval> hcodom=project(f.domain(),range(has,fas));
    return VectorTaylorFunction(hdom,scale(implicit(f.models()),hcodom));
}

VectorTaylorFunction
flow(const VectorFunction& vf, const Vector<Interval>& d, const Float& h, const uint o)
{
    VectorTaylorFunction phi0=embed(VectorTaylorFunction::identity(d),Vector<Interval>(1u,Interval(-h,+h)));
    Vector<Interval> hvfd=0*h*vf(d);
    VectorTaylorFunction phi=phi0+hvfd;

    for(uint i=0; i!=10; ++i) {
        phi=phi0+antiderivative(compose(vf,phi),vf.result_size())*h;
    }

    return phi;
}

VectorTaylorFunction
flow(const VectorTaylorFunction& vf, const Vector<Interval>& d, const Float& h, const uint o)
{
    return flow(vf,d,Interval(-h,+h),o);
}

VectorTaylorFunction
unchecked_flow(const VectorTaylorFunction& vf, const Vector<Interval>& d, const Float& h, const uint o)
{
    return unchecked_flow(vf,d,Interval(-h,+h),o);
}

VectorTaylorFunction
flow(const VectorTaylorFunction& vf, const Vector<Interval>& d, const Interval& h, const uint o)
{
    ARIADNE_ASSERT(subset(d,vf.domain()));
    Vector<Interval> vf_range=vf.range();
    Vector<Interval> euler_step=d+h*vf_range;
    if(!subset(euler_step,vf.domain())) {
        ARIADNE_THROW(FlowBoundsException,"flow(VectorTaylorFunction,Box,Interval,Nat)",
                      std::setprecision(20)<<"Euler step "<<euler_step<<" of vector field "<<vf<<
                      " with range "<<vf_range<<" starting in domain "<<d<<
                      " over time interval "<<h<<" does not remain in domain of vector field.");
    }

    return unchecked_flow(vf,d,h,o);
}


template<class X> inline Vector<X> join(const Vector<X>& v1, const Vector<X>& v2, const X& s3) {
    Vector<X> r(v1.size()+v2.size()+1u);
    for(uint i=0; i!=v1.size(); ++i) { r[i]=v1[i]; }
    for(uint i=0; i!=v2.size(); ++i) { r[v1.size()+i]=v2[i]; }
    r[v1.size()+v2.size()]=s3;
    return r;
}


VectorTaylorFunction
parameterised_flow(const VectorTaylorFunction& vf, const Vector<Interval>& d, const Float& h, const uint o)
{
    ARIADNE_ASSERT(vf.result_size()<=vf.argument_size());
    ARIADNE_ASSERT(d.size()==vf.result_size());
    ARIADNE_ASSERT(h>0.0);

    uint np=vf.argument_size()-vf.result_size();
    uint nx=vf.result_size();

    Vector<Interval> const bp=project(vf.domain(),range(0,np)); // The bounding box of the parameters
    Vector<Interval> const bx=project(vf.domain(),range(np,np+nx)); // The bounding box of the variables
    Vector<Interval> const& dp=bp;  // The domain of the parameters
    Vector<Interval> const& dx=d;  // The domain of the variables

    //std::cerr<<"dx="<<dx<<"\nbx="<<bx<<"\ndp="<<dp<<"\n";
    ARIADNE_ASSERT(subset(dx,bx));
    Vector<Interval> vf_range=vf.range();
    Vector<Interval> euler_step=dx+Interval(0,+h)*vf_range;
    if(!subset(euler_step,bx)) {
        ARIADNE_THROW(FlowBoundsException,"parameterised_flow(VectorTaylorFunction,Box,Float,Nat)",
                      std::setprecision(20)<<"Euler step "<<euler_step<<" of vector field "<<vf<<
                      " with range "<<vf_range<<" starting in domain "<<d<<
                      " over time interval "<<h<<" does not remain in domain of vector field.");
    }

    // Sanity check that vector field domain has nonempty interior
    for(uint i=0; i!=nx; ++i) { ARIADNE_ASSERT_MSG(bx[i].radius()>0.0,"Domain of vector field "<<bx<<" has non-empty interior."); }
    for(uint i=0; i!=nx; ++i) { if(dx[i].radius()<=0) { std::cerr<<"WARNING: Initial set "<<dx<<" has non-empty interior.\n";  } }
    for(uint i=0; i!=np; ++i) { if(dp[i].radius()<=0) { std::cerr<<"WARNING: Parameter set "<<dp<<" has non-empty interior.\n";  } }

    // Scale multiply models by inverse reciprocal of radius
    Vector<TaylorModel> unit_scaled_vf(nx);
    for(uint i=0; i!=nx; ++i) { unit_scaled_vf[i]=vf.models()[i]*(h/rad_ivl(bx[i])); }
    //std::cerr<<"unit_scaled_vf="<<unit_scaled_vf<<std::endl;

    Vector<TaylorModel> y0=embed(np,TaylorModel::scalings(d));
    //std::cerr<<"y0="<<y0<<std::endl;
    Vector<TaylorModel> unit_scaled_y0=unscale(y0,bx);
    //std::cerr<<"unit_scaled_y0"<<unit_scaled_y0<<std::endl;

    Vector<TaylorModel> unit_scaled_flow=parameterised_flow(unit_scaled_vf,unit_scaled_y0,o);
    //std::cerr<<"unit_scaled_flow="<<unit_scaled_flow<<std::endl;

    Vector<TaylorModel> model_flow=scale(unit_scaled_flow,bx);
    //std::cerr<<"model_flow="<<model_flow<<std::endl;

    // Check if flow is only positive
    if(0.0==0) { model_flow=split(model_flow,np+nx,true); }
    //std::cerr<<"split_model_flow="<<model_flow<<std::endl;

    VectorTaylorFunction flow(join(dp,dx,Interval(0,+h)),model_flow);
    //std::cerr<<"\nflow="<<flow<<"\n"<<std::endl;

    return flow;
}

// This method should be used if we know already that the flow over time h remains in
// the domain of the vector field approximation, for example, if this has been
// checked for the original flow
VectorTaylorFunction
unchecked_flow(const VectorTaylorFunction& vf, const Vector<Interval>& d, const Interval& h, const uint o)
{
    uint n=vf.size();
    Float hmag=mag(h);
    const Vector<Interval>& b=vf.domain();

    assert(h.lower()==-h.upper() || h.lower()==0);

    // Sanity check that vector field domain has nonempty interior
    for(uint i=0; i!=n; ++i) { ARIADNE_ASSERT_MSG(b[i].radius()>0.0,"WARNING: Domain of vector field "<<b<<" has non-empty interior."); }
    for(uint i=0; i!=n; ++i) { if(d[i].radius()<=0) { std::cerr<<"WARNING: Initial set "<<d<<" has non-empty interior.\n";  } }


    // Scale multiply models by inverse reciprocal of radius
    Vector<TaylorModel> unit_scaled_vf(vf.size());
    for(uint i=0; i!=vf.size(); ++i) { unit_scaled_vf[i]=vf.models()[i]*(hmag/rad_ivl(b[i])); }
    //std::cerr<<"unit_scaled_vf="<<unit_scaled_vf<<std::endl;

    Vector<TaylorModel> y0=TaylorModel::scalings(d);
    //std::cerr<<"y0="<<y0<<std::endl;
    Vector<TaylorModel> unit_scaled_y0=embed(unscale(y0,b),1u);
    //std::cerr<<"unit_scaled_y0"<<unit_scaled_y0<<std::endl;

    Vector<TaylorModel> unit_scaled_flow=unchecked_flow(unit_scaled_vf,unit_scaled_y0,o);
    //std::cerr<<"unit_scaled_flow="<<unit_scaled_flow<<std::endl;

    Vector<TaylorModel> model_flow=scale(unit_scaled_flow,b);
    //std::cerr<<"model_flow="<<model_flow<<std::endl;

    if(h.lower()==0) { model_flow=split(model_flow,d.size(),true); }

    VectorTaylorFunction flow(join(d,h),model_flow);
    //std::cerr<<"\nflow="<<flow<<"\n"<<std::endl;

    return flow;
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

Float distance(const VectorTaylorFunction& f1, const VectorFunction& f2) {
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



/*
latexstream&
operator<<(Output::latexstream& texs, const VectorTaylorFunction& p)
{
    using namespace Function;
    texs << "%VectorTaylorFunction\n";
    texs << "\\ensuremath{\n";
    texs << "\\left( \\begin{array}{c}\n";
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
    texs << "\\end{array}\\right)\n}\n";
    return texs;
}
*/


} // namespace Ariadne
