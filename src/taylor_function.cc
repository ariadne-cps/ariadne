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
#include "function_interface.h"
#include "taylor_variable.h"
#include "taylor_function.h"

namespace Ariadne {

typedef unsigned int uint;

typedef Vector<Float> Point;
typedef Vector<Interval> Box;


namespace {
Vector<Interval> unscale(const Vector<Interval>& x, const Vector<Interval>& d) {
    Vector<Interval> r(x);
    for(uint i=0; i!=r.size(); ++i) {
        (r[i]-=med_ivl(d[i]))/=rad_ivl(d[i]); }
    return r;
}
} // namespace


TaylorFunction::TaylorFunction()
    : _domain(), _models()
{
}




TaylorFunction::TaylorFunction(const Vector<Interval>& d,
                               const Vector<TaylorModel>& f)
    : _domain(d), _models(f)
{
    for(uint i=0; i!=f.size(); ++i) {
        ARIADNE_ASSERT(d.size()==f[i].argument_size());
    }
}

TaylorFunction::TaylorFunction(const Vector<Interval>& d,
                               const Vector< Expansion<Float> >& f)
    : _domain(d), _models(f.size())
{
    for(uint i=0; i!=f.size(); ++i) {
        ARIADNE_ASSERT(d.size()==f[i].argument_size());
        _models[i]=TaylorModel(f[i]);
    }
}

TaylorFunction::TaylorFunction(const Vector<Interval>& d,
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



TaylorFunction::TaylorFunction(const Vector<Interval>& d,
                               const FunctionInterface& f)
    : _domain(d), _models(f.result_size())
{
    ARIADNE_ASSERT(f.result_size()>0);
    ARIADNE_ASSERT(d.size()==f.argument_size());
    Vector<TaylorModel> x=TaylorModel::scalings(d);
    this->_models=f.evaluate(x);
}

TaylorFunction::TaylorFunction(const Vector<Interval>& d,
                               const FunctionInterface& f,
                               const shared_ptr<TaylorModel::Accuracy> accuracy_ptr)
    : _domain(d), _models(f.result_size())
{
    ARIADNE_ASSERT(f.result_size()>0);
    ARIADNE_ASSERT(d.size()==f.argument_size());
    Vector<TaylorModel> x=TaylorModel::scalings(d);
    for(uint i=0; i!=x.size(); ++i) { x[i].accuracy_ptr()=accuracy_ptr; }
    this->_models=f.evaluate(x);
}


TaylorFunction::TaylorFunction(const Vector<Interval>& d,
                               const Vector< Polynomial<Float> >& p)
    : _domain(d), _models(p.size())
{
    for(uint i=0; i!=p.size(); ++i) { ARIADNE_ASSERT(d.size()==p[i].argument_size()); }

    Vector<TaylorModel> x=TaylorModel::scalings(d);
    this->_models=Ariadne::evaluate(p,x);
}

TaylorFunction::TaylorFunction(const Vector<TaylorVariable>& v)
    : _domain(), _models(v.size())
{
    ARIADNE_ASSERT(v.size()>0);
    for(uint i=1; i!=v.size(); ++i) { ARIADNE_ASSERT(v[i].domain()==v[0].domain()); }
    this->_domain=v[0].domain();
    for(uint i=0; i!=v.size(); ++i) {
        this->_models[i]=v[i].model();
    }
}










TaylorFunction
TaylorFunction::constant(const Vector<Interval>& d, const Vector<Float>& c)
{
    return TaylorFunction(d,TaylorModel::constants(d.size(),c));
}

TaylorFunction
TaylorFunction::constant(const Vector<Interval>& d, const Vector<Interval>& c)
{
    return TaylorFunction(d,TaylorModel::constants(d.size(),c));
}

TaylorFunction
TaylorFunction::identity(const Vector<Interval>& d)
{
    return TaylorFunction(d,TaylorModel::scalings(d));
}


Polynomial<Interval> polynomial(const TaylorModel& tm) {
    return Polynomial<Interval>(tm.expansion())+Interval(-tm.error(),+tm.error());
}

Vector< Polynomial<Interval> > 
TaylorFunction::polynomial() const
{
    Vector<Polynomial<Interval> > p(this->result_size());
    for(uint i=0; i!=this->result_size(); ++i) {
        p[i]=Ariadne::polynomial(this->models()[i]); }
    Vector<Polynomial<Interval> > s(this->argument_size());
    for(uint j=0; j!=this->argument_size(); ++j) {
        s[j]=Ariadne::polynomial(TaylorModel::unscaling(this->argument_size(),j,this->domain()[j])); }
    return compose(p,s);
}

bool
TaylorFunction::operator==(const TaylorFunction& tm) const
{
    return this->_models==tm._models;
}



bool
TaylorFunction::operator!=(const TaylorFunction& p2) const
{
    return !(*this==p2);
}



const Vector<Interval>&
TaylorFunction::domain() const
{
    return this->_domain;
}


const Vector<Interval>
TaylorFunction::range() const
{
    Vector<Interval> result(this->result_size());
    for(uint i=0; i!=result.size(); ++i) {
        result[i]=this->_models[i].range();
    }
    return result;
}



const Vector<TaylorModel>&
TaylorFunction::models() const
{
    return this->_models;
}




uint
TaylorFunction::argument_size() const
{
    return this->_domain.size();
}


uint
TaylorFunction::result_size() const
{
    return this->_models.size();
}


TaylorVariable
TaylorFunction::operator[](uint i) const
{
    return TaylorVariable(this->_domain,this->_models[i]);
}





















TaylorFunction
TaylorFunction::truncate(ushort degree) const
{
    ARIADNE_NOT_IMPLEMENTED;
}



Vector<Interval>
TaylorFunction::evaluate(const Vector<Float>& x) const
{
    return this->evaluate(Vector<Interval>(x));
}


Vector<Interval>
TaylorFunction::evaluate(const Vector<Interval>& x) const
{
    Vector<Interval> sx=Ariadne::evaluate(TaylorModel::scalings(this->_domain),x);
    return Ariadne::evaluate(this->_models,sx);
}


Matrix<Interval>
TaylorFunction::jacobian(const Vector<Interval>& x) const
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


TaylorFunction
join(const TaylorFunction& f, const TaylorFunction& g)
{
    ARIADNE_ASSERT(f.domain()==g.domain());
    return TaylorFunction(f.domain(),join(f.models(),g.models()));
}

TaylorFunction
restrict(const TaylorFunction& f, const Vector<Interval>& d)
{
    ARIADNE_ASSERT(subset(d,f.domain()));
    if(d==f.domain()) { return f; }
    Vector<TaylorModel> s=TaylorModel::rescalings(f.domain(),d);
    return TaylorFunction(d,compose(f._models,s));
}



TaylorFunction&
operator+=(TaylorFunction& f, const Vector<Interval>& e)
{
    ARIADNE_ASSERT(f.result_size()==e.size());
    f._models+=e;
    return f;
}


TaylorFunction
operator+(const TaylorFunction& f1, const TaylorFunction& f2)
{
    ARIADNE_ASSERT_MSG(!intersection(f1.domain(),f2.domain()).empty(),
                       "operator+(TaylorFunction f1, TaylorFunction f2) with f1="<<f1<<" f2="<<f2<<
                       ": domains are disjoint");
    if(f1.domain()==f2.domain()) {
        return TaylorFunction(f1.domain(),Vector<TaylorModel>(f1.models()+f2.models()));
    } else {
        Box new_domain=intersection(f1.domain(),f2.domain());
        return operator+(restrict(f1,new_domain),restrict(f2,new_domain));
    }
}


TaylorFunction
operator-(const TaylorFunction& f1, const TaylorFunction& f2)
{
    ARIADNE_ASSERT(!intersection(f1.domain(),f2.domain()).empty());
    if(f1.domain()==f2.domain()) {
        return TaylorFunction(f1.domain(),Vector<TaylorModel>(f1.models()+f2.models()));
    } else {
        Box new_domain=intersection(f1.domain(),f2.domain());
        return operator-(restrict(f1,new_domain),restrict(f2,new_domain));
    }
}



TaylorFunction
operator+(const TaylorFunction& f, const Vector<Float>& c)
{
    return TaylorFunction(f.domain(),Vector<TaylorModel>(f.models()+c));
}

TaylorFunction
operator+(const TaylorFunction& f, const Vector<Interval>& c)
{
    return TaylorFunction(f.domain(),Vector<TaylorModel>(f.models()+c));
}

TaylorFunction
operator-(const TaylorFunction& f, const Vector<Float>& c)
{
    return TaylorFunction(f.domain(),Vector<TaylorModel>(f.models()-c));
}

TaylorFunction
operator-(const TaylorFunction& f, const Vector<Interval>& c)
{
    return TaylorFunction(f.domain(),Vector<TaylorModel>(f.models()-c));
}

TaylorFunction
operator*(const Matrix<Interval>& A, const TaylorFunction& f)
{
    return TaylorFunction(f.domain(),Vector<TaylorModel>(prod(A,f.models())));
}





TaylorFunction
combine(const TaylorFunction& f1, const TaylorFunction& f2)
{
    return TaylorFunction(join(f1.domain(),f2.domain()),combine(f1.models(),f2.models()));
}


TaylorFunction
compose(const FunctionInterface& g, const TaylorFunction& f)
{
    return TaylorFunction(f.domain(),g.evaluate(f.models()));
}

TaylorFunction
compose(const TaylorFunction& g, const TaylorFunction& f)
{
    ARIADNE_ASSERT_MSG(subset(f.range(),g.domain()),"f.range()="<<f.range()<<" is not a subset of g.domain()="<<g.domain());
    return TaylorFunction(f.domain(),Ariadne::compose(g.models(),unscale(f.models(),g.domain())));
}


TaylorVariable
compose(const TaylorVariable& g, const TaylorFunction& f)
{
    ARIADNE_ASSERT_MSG(subset(f.range(),g.domain()),"f.range()="<<f.range()<<" is not a subset of g.domain()="<<g.domain());
    return TaylorVariable(f.domain(),Ariadne::unchecked_compose(g.model(),unscale(f.models(),g.domain())));
}

TaylorVariable
unchecked_compose(const TaylorVariable& g, const TaylorFunction& f)
{
        return TaylorVariable(f.domain(),Ariadne::unchecked_compose(g.model(),unscale(f.models(),g.domain())));
}



TaylorFunction
antiderivative(const TaylorFunction& f, uint k)
{
    ARIADNE_ASSERT_MSG(k<f.argument_size(),"f="<<f<<"\n f.argument_size()="<<f.argument_size()<<" k="<<k);
    Interval fdomkrad=rad_ivl(f.domain()[k].lower(),f.domain()[k].upper());
    TaylorFunction g=f;
    for(uint i=0; i!=g.size(); ++i) {
        g._models[i].antidifferentiate(k);
        g._models[i]*=fdomkrad;
    }
    return g;
}

TaylorFunction
implicit(const TaylorFunction& f)
{
    ARIADNE_ASSERT_MSG(f.argument_size()>f.result_size(),"f.argument_size()<=f.result_size() in implicit(f): f="<<f);
    uint fas=f.argument_size();
    uint has=f.argument_size()-f.result_size();
    Vector<Interval> hdom=project(f.domain(),range(0,has));
    Vector<Interval> hcodom=project(f.domain(),range(has,fas));
    return TaylorFunction(hdom,scale(implicit(f.models()),hcodom));
}

TaylorFunction
flow(const TaylorFunction& vf, const Vector<Interval>& d, const Interval& h, const uint o)
{
    ARIADNE_ASSERT(subset(d,vf.domain()));
    Vector<Interval> vf_range=vf.range();
    Vector<Interval> euler_step=d+h*vf_range;
    if(!subset(euler_step,vf.domain())) {
        ARIADNE_THROW(FlowBoundsException,"flow(TaylorFunction,Box,Interval,Nat)",
                      std::setprecision(20)<<"Euler step "<<euler_step<<" of vector field "<<vf<<
                      " with range "<<vf_range<<" starting in domain "<<d<<
                      " over time interval "<<h<<" does not remain in domain of vector field.");
    }

    return unchecked_flow(vf,d,h,o);
}

// This method should be used if we know already that the flow over time h remains in
// the domain of the vector field approximation, for example, if this has been
// checked for the original flow
TaylorFunction
unchecked_flow(const TaylorFunction& vf, const Vector<Interval>& d, const Interval& h, const uint o)
{
    Vector<TaylorModel> unscaled_vf=unscale(vf.models(),vf.domain());
    Vector<Interval> unscaled_d=unscale(d,vf.domain());
    //std::cerr<<"\n"<<unscaled_vf<<"\n"<<unscaled_d<<"\n";
    //Vector<TaylorModel> unscaled_flw=flow(unscaled_vf,unscaled_d,h,o);
    //std::cerr<<unscaled_flw<<"\n";
    //std::cerr<<scale(unscaled_flw,vf.domain());
    return TaylorFunction(join(d,h),
        scale(unchecked_flow(unscale(vf.models(),vf.domain()),unscale(d,vf.domain()),h,o),vf.domain()));
    
}




bool
refines(const TaylorFunction& f, const TaylorFunction& g)
{
    return refines(f.models(),g.models());
}



std::ostream&
TaylorFunction::write(std::ostream& os) const
{
    os << "TaylorFunction( "<<this->domain()<<" , ";
    for(uint i=0; i!=this->result_size(); ++i) {
        os << (i==0?'[':',')<<this->_models[i].expansion()<<","<<this->_models[i].error();
    }
    return os << "] )";
}



std::ostream&
operator<<(std::ostream& os, const TaylorFunction& p)
{
    return p.write(os);
}



/*
latexstream&
operator<<(Output::latexstream& texs, const TaylorFunction& p)
{
    using namespace Function;
    texs << "%TaylorFunction\n";
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



TaylorExpression::TaylorExpression()
    : _domain(), _model() { }

TaylorExpression::TaylorExpression(const Vector<Interval>& domain)
    : _domain(domain), _model(domain.size()) { }

TaylorExpression::TaylorExpression(const Vector<Interval>& domain, const TaylorModel& model)
    : _domain(domain), _model(model) { }

TaylorExpression TaylorExpression::variable(const Vector<Interval>& domain, Nat i) {
    TaylorExpression t(domain); t._model=TaylorModel::scaling(domain.size(),i,domain[i]); return t; }

Vector<TaylorExpression> TaylorExpression::variables(const Vector<Interval>& domain) {
    Vector<TaylorExpression> r(domain.size(),TaylorExpression(domain));
    for(uint i=0; i!=r.size(); ++i) { r[i]=variable(domain,i); }
    return r; }

const Vector<Interval>& TaylorExpression::domain() const {
    return  this->_domain; }
const TaylorModel& TaylorExpression::model() const {
    return  this->_model; }

TaylorExpression operator+(const Float& c, const TaylorExpression& t) {
    return TaylorExpression(t._domain,c+t._model); }
TaylorExpression operator-(const Float& c, const TaylorExpression& t) {
    return TaylorExpression(t._domain,c-t._model); }
TaylorExpression operator*(const Float& c, const TaylorExpression& t) {
    return TaylorExpression(t._domain,c*t._model); }
TaylorExpression operator/(const Float& c, const TaylorExpression& t) {
    return TaylorExpression(t._domain,c/t._model); }

TaylorExpression operator+(const Interval& c, const TaylorExpression& t) {
    return TaylorExpression(t._domain,c+t._model); }
TaylorExpression operator-(const Interval& c, const TaylorExpression& t) {
    return TaylorExpression(t._domain,c-t._model); }
TaylorExpression operator*(const Interval& c, const TaylorExpression& t) {
    return TaylorExpression(t._domain,c*t._model); }
TaylorExpression operator/(const Interval& c, const TaylorExpression& t) {
    return TaylorExpression(t._domain,c/t._model); }

TaylorExpression operator+(const TaylorExpression& t, const Float& c) {
    return TaylorExpression(t._domain,t._model+c); }
TaylorExpression operator-(const TaylorExpression& t, const Float& c) {
    return TaylorExpression(t._domain,t._model-c); }
TaylorExpression operator*(const TaylorExpression& t, const Float& c) {
    return TaylorExpression(t._domain,t._model*c); }
TaylorExpression operator/(const TaylorExpression& t, const Float& c) {
    return TaylorExpression(t._domain,t._model/c); }

TaylorExpression operator+(const TaylorExpression& t, const Interval& c) {
    return TaylorExpression(t._domain,t._model+c); }
TaylorExpression operator-(const TaylorExpression& t, const Interval& c) {
    return TaylorExpression(t._domain,t._model-c); }
TaylorExpression operator*(const TaylorExpression& t, const Interval& c) {
    return TaylorExpression(t._domain,t._model*c); }
TaylorExpression operator/(const TaylorExpression& t, const Interval& c) {
    return TaylorExpression(t._domain,t._model/c); }

TaylorExpression operator+(const TaylorExpression& t) {
    return TaylorExpression(t._domain,+t._model);  }
TaylorExpression operator-(const TaylorExpression& t) {
    return TaylorExpression(t._domain,-t._model);  }
TaylorExpression operator+(const TaylorExpression& t1, const TaylorExpression& t2) {
    ARIADNE_ASSERT(t1._domain==t2._domain); return TaylorExpression(t1._domain,t1._model+t2._model); }
TaylorExpression operator-(const TaylorExpression& t1, const TaylorExpression& t2) {
    ARIADNE_ASSERT(t1._domain==t2._domain); return TaylorExpression(t1._domain,t1._model-t2._model); }
TaylorExpression operator*(const TaylorExpression& t1, const TaylorExpression& t2) {
    ARIADNE_ASSERT(t1._domain==t2._domain); return TaylorExpression(t1._domain,t1._model*t2._model); }
TaylorExpression operator/(const TaylorExpression& t1, const TaylorExpression& t2) {
    ARIADNE_ASSERT(t1._domain==t2._domain); return TaylorExpression(t1._domain,t1._model/t2._model); }

TaylorFunction operator*(const TaylorExpression& t, const Vector<Float>& v) {
    return TaylorFunction(t.domain(),Vector<TaylorModel>(t._model*v)); }
TaylorFunction operator*(const Vector<Float>& v, const TaylorExpression& t) {
    return t*v; };

std::ostream& operator<<(std::ostream& os, const TaylorExpression& t) {
    return os<<"("<<t._domain<<","<<t._model<<")"; }

} // namespace Ariadne
