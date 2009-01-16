/***************************************************************************
 *            taylor_model.cc
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
 
#include "macros.h"
#include "exceptions.h"
#include "numeric.h"
#include "vector.h"
#include "matrix.h"
#include "multi_index.h"
#include "sparse_differential.h"
#include "differential_vector.h"
#include "function_interface.h"
#include "taylor_variable.h"
#include "taylor_model.h"

namespace Ariadne {

typedef Vector<Float> Point;
typedef Vector<Interval> Box;

TaylorModel::TaylorModel() 
    : _domain(), 
      _expansion()
{
}


TaylorModel::TaylorModel(uint rs, uint as, ushort o, ushort s) 
    : _domain(as,I(-1,1)),
      _expansion(rs,as,o)
{
}


TaylorModel::TaylorModel(const Vector<Interval>& d,
                         const Vector<TaylorVariable>& e)
    : _domain(d),
      _expansion(e)
{
    for(uint i=0; i!=e.size(); ++i) {
        ARIADNE_ASSERT(d.size()==e[i].argument_size());
    }
}

TaylorModel::TaylorModel(const Vector<Interval>& d,
                         const Vector<SparseDifferential<Float> >& e)
    : _domain(d),
      _expansion(e.size())
{
    for(uint i=0; i!=e.size(); ++i) {
        ARIADNE_ASSERT(d.size()==e[i].argument_size());
        this->_expansion[i]=TaylorVariable(e[i],Interval(0));
    }
}




TaylorModel::TaylorModel(const Vector<Interval>& d, 
                         ushort o, ushort s,
                         const FunctionInterface& f)
    : _domain(d),
      _expansion(f.result_size())
{
    ARIADNE_NOT_IMPLEMENTED;
}




/*
TaylorModel
TaylorModel::zero(uint rs, uint as)  
{
    return TaylorModel(Vector<Interval>(as,Interval(-inf(),+inf()),TaylorVariable::constants(rs,as,Vector<Float>(rs,0.0))));
}
*/




TaylorModel
TaylorModel::constants(uint as, const Vector<Float>& c)  
{
    ARIADNE_NOT_IMPLEMENTED;
}

TaylorModel
TaylorModel::variables(uint as, const Vector<Float>& c)  
{
    ARIADNE_NOT_IMPLEMENTED;
}



bool
TaylorModel::operator==(const TaylorModel& tm) const
{
    return this->_domain==tm._domain
        && this->_expansion==tm._expansion;
}



bool
TaylorModel::operator!=(const TaylorModel& p2) const
{
    return !(*this==p2);
}



const Vector<Interval>&
TaylorModel::domain() const
{ 
    return this->_domain; 
}


const Vector<Interval>
TaylorModel::range() const
{ 
    Vector<Interval> result(this->result_size());
    for(uint i=0; i!=result.size(); ++i) {
        result[i]=this->_expansion[i].range();
    }
    return result;
}



const Vector<TaylorVariable>&
TaylorModel::expansion() const
{ 
    return this->_expansion;
}




uint 
TaylorModel::argument_size() const
{ 
    return this->_expansion[0].argument_size(); 
}


uint 
TaylorModel::result_size() const 
{ 
    return this->_expansion.size();
}


    
      


















TaylorModel 
TaylorModel::truncate(ushort degree) const
{
    ARIADNE_NOT_IMPLEMENTED;
}



Vector<Interval> 
TaylorModel::evaluate(const Vector<Float>& x) const
{
    return this->evaluate(Vector<Interval>(x));
}


Vector<Interval> 
TaylorModel::evaluate(const Vector<Interval>& x) const
{
    if(this->argument_size()!=x.size()) {
        ARIADNE_THROW(std::runtime_error,"TaylorModel::evaluate(Vector)","Incompatible argument size");
    }

    //TODO: Make this MUCH more efficient!

    // Scale x to domain
    Vector<Interval> scaled_x(x.size());
    for(uint i=0; i!=x.size(); ++i) {
        const Interval& d=this->_domain[i];
        Interval dm=add_ivl(d.l/2,d.u/2);
        Interval dr=sub_ivl(d.u/2,d.l/2);
        scaled_x[i]=(x[i]-dm)/dr;
    }

    // Compute result
    Vector<Interval> result(this->result_size());
    for(uint i=0; i!=result.size(); ++i) {
        result[i]=this->_expansion[i].evaluate(scaled_x);
    }
    return result;
}


Matrix<Float> 
TaylorModel::jacobian(const Vector<Float>& x) const
{
    DifferentialVector< SparseDifferential<Float> > y(this->argument_size(),SparseDifferential<Float>(this->argument_size(),1u));
    for(uint j=0; j!=this->argument_size(); ++j) {
        y[j].set_value(this->domain()[j].midpoint());
        y[j].set_gradient(j,this->domain()[j].radius());
    }
    DifferentialVector< SparseDifferential<Float> > t(this->result_size());
    for(uint i=0; i!=this->result_size(); ++i) {
        t[i]=this->_expansion[i].expansion();
    }
    return get_jacobian(compose(t-x,y));
}


TaylorModel
recentre(const TaylorModel& tm, const Box& bx, const Point& c)
{
    ARIADNE_NOT_IMPLEMENTED;
}



TaylorModel
operator+(const TaylorModel& p1, const TaylorModel& p2)
{
    ARIADNE_ASSERT(!intersection(p1.domain(),p2.domain()).empty());
    if(p1.domain()==p2.domain()) {
        return TaylorModel(p1._domain,Vector<TaylorVariable>(p1._expansion+p2._expansion));
    } else {
        ARIADNE_NOT_IMPLEMENTED;
    }
}


TaylorModel
operator-(const TaylorModel& p1, const TaylorModel& p2)
{
    ARIADNE_ASSERT(!intersection(p1.domain(),p2.domain()).empty());
    if(p1.domain()==p2.domain()) {
        return TaylorModel(p1._domain,Vector<TaylorVariable>(p1._expansion-p2._expansion));
    } else {
        ARIADNE_NOT_IMPLEMENTED;
        Box new_domain=intersection(p1.domain(),p2.domain());
        //Point new_centre=new_domain.centre();
        Point new_centre=midpoint(new_domain);
        return operator-(recentre(p1,new_domain,new_centre),recentre(p2,new_domain,new_centre));
    }
}



/*

TaylorModel
mul(const TaylorModel& tm, const Float& x)
{
    return TaylorModel(tm.domain(),tm.centre(),tm.smoothness(),tm.centre_derivatives()*x);
}

*/


/*

void
mul(TaylorModel& p0, const TaylorModel& p1, const TaylorModel& p2)
{
    if(p1.result_size()!=1u) {
        ARIADNE_THROW(IncompatibleSizes,"mul(TaylorModel,TaylorModel)","p1.result_size()="<<p1.result_size());
    }
    if(p2.result_size()!=1u) {
        ARIADNE_THROW(IncompatibleSizes,"mul(TaylorModel,TaylorModel)","p2.result_size()="<<p2.result_size());
    }
    if(p1.argument_size()!=p2.argument_size()) {
        ARIADNE_THROW(IncompatibleSizes,"add(TaylorModel p1,TaylorModel p2)","p1.result_size()="<<p1.result_size()<<", p2.result_size()="<<p2.result_size());
    }

    uint rs=1u;
    uint as=p1.argument_size();
    uint d1=p1.order();
    uint d2=p2.order();
    uint d=d1+d2;
    uint s1=p1.smoothness();
    uint s2=p2.smoothness();
    uint s=std::max(s1,s2);

    p0.resize(rs,as,d,s);
    MultiIndex j0(as);
    for(MultiIndex j1(as); j1.degree()<=d1; ++j1) {
        const Float& x1=p1.get(0u,j1);
        for(MultiIndex j2(as); j2.degree()<=d2; ++j2) {
            const R& x2=p2.get(0u,j2);
            j0=j1+j2;
            p0.at(0u,j0)+=x1*x2;
        }
    }
}


*/





/*
TaylorModel&
operator*=(TaylorModel& p0, const Interval& x1)
{
    for(uint i=0; i!=p0._data.size(); ++i) {
        p0._data[i]*=x1;
    }
    return p0;
}
*/

TaylorModel
compose(const TaylorModel& p1, const TaylorModel& p2)
{
    ARIADNE_ASSERT(subset(p2.range(),p1.domain()));
    return TaylorModel(p2.domain(),Ariadne::compose(p1.expansion(),p1.domain(),p2.expansion()));
}



TaylorModel
antiderivative(const TaylorModel& tm, uint k) 
{
    return TaylorModel(tm.domain(),antiderivative(tm.expansion(),k));
}



/*
void
compose(TaylorModel& p0, const TaylorModel& p1, const TaylorModel& p2)
{
    // TODO: Improve this algorithm as it's critical!!
    if(p1.argument_size()!=p2.result_size()) {
        ARIADNE_THROW(IncompatibleSizes,"compose(TaylorModel p1,TaylorModel p2)","p1.argument_size()="<<p1.argument_size()<<", p2.result_size()="<<p2.result_size());
    }

    if(p1.order()==0) { 
        p0=static_cast< TaylorModel >(p1); 
        return;
    }

    p0.resize(p1.result_size(),p2.argument_size(),p1.order()*p2.order(),std::max(p1.smoothness(),p2.smoothness()));
  
    TaylorModel* all_powers=new TaylorModel[p2.result_size()*(p1.order()+1)];
    TaylorModel* powers[p2.result_size()];
    for(uint i=0; i!=p2.result_size(); ++i) {
        powers[i]=all_powers+i*(p1.order()+1);
    }

    for(uint i=0; i!=p2.result_size(); ++i) {
        powers[i][0]=TaylorModel::one(p2.argument_size());
        powers[i][1]=p2.component(i);
        if(p1.order()>=2) {
            powers[i][2]=pow(powers[i][1],2);
        }
        for(uint j=3; j<=p1.order(); ++j) {
            powers[i][j]=powers[i][2]*powers[i][j-2];
        }
    }
  
    TaylorModel* results=new TaylorModel[p1.result_size()];
    for(uint i=0; i!=p1.result_size(); ++i) {
        results[i]=TaylorModel::zero(1u,p2.argument_size());
    }

    for(uint i=0; i!=p1.result_size(); ++i) {
        for(MultiIndex j(p1.argument_size()); j.degree()<=p1.order(); ++j) {
            TaylorModel t=TaylorModel::constant(p2.argument_size(),p1.get(i,j));
            for(uint k=0; k!=p1.argument_size(); ++k) {
                t=t*powers[k][j[k]];
            }
            results[i]=results[i]+t;
        }
    }
  
    for(uint i=0; i!=p0.result_size(); ++i) {
        for(uint j=0; j!=p0.data().size()/p0.result_size(); ++j) {
            p0._data[i+j*p0.result_size()]=results[i].data()[j];
        }
    }
  
    delete[] results;
    delete[] all_powers;
}

*/

/*

void
derivative(TaylorModel& p0, const TaylorModel& p1, uint k)
{
    if(p1.smoothness()==0) {
        ARIADNE_THROW(std::runtime_error,"derivative(TaylorModel,uint)"," model has smoothness 0");
    }

    p0.resize(p1.result_size(),p1.argument_size(),p1.order()-1,p1.smoothness()-1);
  
    MultiIndex dj(p1.argument_size());

    for(uint i=0; i!=p0.result_size(); ++i) {
        for(MultiIndex j(p1.argument_size()); j.degree()<=p1.order(); ++j) {
            if(j[k]!=0) {
                dj=j;
                dj.decrement_index(k);
                p0.at(i,dj)+=static_cast<int>(j[k])*p1.get(i,j);
            }
        }
    }
}

*/

 /*
Matrix<Interval> 
TaylorModel::jacobian(const Vector<Float>& s) const
{
    return this->jacobian(Vector<Interval>(s));
}


Matrix<Interval> 
TaylorModel::jacobian(const Vector<Interval>& x) const
{
    Matrix<Interval> J(this->result_size(),this->argument_size());
    Box w=x-this->centre();
    array< array<I> > powers=this->_powers(w);

    for(uint j=0; j!=this->argument_size(); ++j) {
        for(MultiIndex m(this->argument_size()); m.degree()<this->order(); ++m) {
            MultiIndex n=m;
            int c=n[j];
            if(c!=0) {
                n.decrement_index(j);
                I a=c;
                for(uint k=0; k!=this->argument_size(); ++k) {
                    a*=powers[k][n[k]];
                }
                for(uint i=0; i!=this->result_size(); ++i) {
                    I amim=a*this->_centre_derivatives[i][m];
                    J[i][j]+=amim;
                }
            }
        }
        for(MultiIndex m=MultiIndex::first(this->argument_size(),this->order()); m.degree()<=this->order(); ++m) {
            MultiIndex n=m;
            int c=n[j];
            if(c!=0) {
                n.decrement_index(j);
                I a=c;
                for(uint k=0; k!=this->argument_size(); ++k) {
                    a*=powers[k][n[k]];
                }
                for(uint i=0; i!=this->result_size(); ++i) {
                    J[i][j]+=a*this->_domain_derivatives[i][m];
                }
            }
        }
    }
  
    return J;
}
 */



 
TaylorModel 
inverse(const TaylorModel& p, const Point& v)
{
    assert(p.result_size()==p.argument_size());
    assert(p.argument_size()==v.size());

    // The following are only to simplfy testing.
    assert(v==Point(v.size(),0));
    assert(p.evaluate(v)==Point(p.result_size(),0));
    typedef Interval I;

    Point c=midpoint(Box(p.evaluate(v)));
    Matrix<Float> J=p.jacobian(v);
  
    Matrix<Float> invJ=inverse(J);

    // FIXME: Need to re-solve for image of centre. What should initial set be? Different code needed for Rational?
    Box invf=v;

    // FIXME: Give correct initial conditions
    TaylorModel result;

    for(MultiIndex m(p.result_size()); m.degree()<=2; ++m) {
        if(m.degree()==0) {
            for(uint i=0; i!=p.argument_size(); ++i) {
                result._expansion[i][m]=v[i];
            }
        } else if(m.degree()==1) {
            for(uint i=0; i!=p.argument_size(); ++i) {
                result._expansion[i][m]=invJ[i][m.position()-1];
            }
        } else {
            // FIXME: Add code for higher indices
        }
    }
    return result;
}


 
TaylorModel 
implicit(const TaylorModel& p, const Point& v)
{
    ARIADNE_NOT_IMPLEMENTED;
}


array< array<Interval> >
TaylorModel::_powers(const Vector<Interval>& v) const
{
    uint order=21;
    array< array<I> > powers(this->argument_size(), array<I>(order));
    for(uint i=0; i!=this->argument_size(); ++i) {
        powers[i][0]=1;
        if(order>=1) {
            powers[i][1]=v(i);
            if(order>=2) {
                powers[i][2]=pow(v(i),2);
                for(uint j=3; j<=order; ++j) {
                    powers[i][j]=powers[i][2]*powers[i][j-2];
                }
            }
        }
    }
    return powers;
}



std::ostream&
TaylorModel::write(std::ostream& os) const 
{
    os << "TaylorModel(\n";
    for(uint i=0; i!=this->result_size(); ++i) {
        os << "  domain=" << this->domain() << ",\n" << std::flush;
        os << "  range=" << this->range() << ",\n" << std::flush;
        os << "  expansion=" << this->_expansion << ",\n" << std::flush;
    }
    os << ")\n";
    return os;
}



std::ostream&
operator<<(std::ostream& os, const TaylorModel& p)
{
    return p.write(os);
}



/*
latexstream&
operator<<(Output::latexstream& texs, const TaylorModel& p)
{
    using namespace Function;
    texs << "%TaylorModel\n";
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
