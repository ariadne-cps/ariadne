/***************************************************************************
 *            taylor_set.cc
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
#include "function.h"
#include "approximate_taylor_model.h"

#include "geometry.h"
#include "taylor_variable.h"
#include "taylor_set.h"

#include "grid_set.h"

#include "zonotope.h"
#include "polytope.h"
#include "graphics.h"

namespace Ariadne {

TaylorSet::TaylorSet(uint d) 
    : _variables(d)
{
}


TaylorSet::TaylorSet(const Vector<TaylorVariable>& tvs) 
    : _variables(tvs.size())
{
    for(size_t i=0; i!=tvs.size(); ++i) {
        this->_variables[i]=TaylorVariable(tvs[i]);
    }
}








tribool 
TaylorSet::disjoint(const Box& bx) const
{
    return this->bounding_box().disjoint(bx) || indeterminate;
}


tribool 
TaylorSet::overlaps(const Box& bx) const
{
    ARIADNE_NOT_IMPLEMENTED;
}


tribool
TaylorSet::inside(const Box& bx) const
{
    Vector<Interval> bb=this->bounding_box();
    return Ariadne::inside(bb,bx) || indeterminate;
}


Box
TaylorSet::bounding_box() const
{
    Box r(this->dimension());
    for(uint i=0; i!=this->dimension(); ++i) {
        r[i]=(*this)[i].range();
    }
    return r;
}


Zonotope 
zonotope(const TaylorSet& ts)
{
    uint d=ts.dimension();
    uint ng=ts.number_of_generators();
    Vector<Float> c(d);
    Matrix<Float> G(d,ng);
    Vector<Float> e(d);
    for(uint i=0; i!=d; ++i) {
        if(ts[i].error()==0) {
            c[i]=ts[i].expansion().value();
            e[i]=0;
        } else {
            Interval ce=ts[i].expansion().value()+ts[i].error();
            c[i]=midpoint(ce);
            e[i]=radius(ce);
        }
        for(uint j=0; j!=ng; ++j) {
            G[i][j]=ts[i].expansion().gradient(j);
        }
    }

    for(uint i=0; i!=d; ++i) {
        const SparseDifferential<Float>& sd=ts[i].expansion();
        for(SparseDifferential<Float>::const_iterator iter=sd.begin(); iter!=sd.end(); ++iter) {
            if(iter->first.degree()>=2) {
                e[i]=add_up(e[i],abs(iter->second));
            }
        }
    }
    return Zonotope(c,G,e);
}




pair<TaylorSet,TaylorSet>
TaylorSet::split(uint j) const
{
    const TaylorSet& ts=*this;
    pair<TaylorSet,TaylorSet> result(ts.dimension(),ts.dimension());
    for(uint i=0; i!=ts.dimension(); ++i) {
        make_lpair(result.first[i],result.second[i])=Ariadne::split(ts[i],j);
    }
    return result;
}


void
_adjoin_outer_approximation(const TaylorSet& set, const Box& domain, Float eps, GridTreeSet& grid_set, uint depth)
{
    //std::cerr<<"adjoin_outer_approximation(TaylorSet,Box,Float,GridTreeSet,Nat)\n  domain="<<domain<<"\n";
    uint d=set.dimension();
    Box range(set.dimension());
    for(uint i=0; i!=set.dimension(); ++i) {
        range[i]=evaluate(set[i].expansion(),domain);
    }
    if(range.radius()<eps) {
        for(uint i=0; i!=set.dimension(); ++i) {
            range[i]+=set[i].error();
        }
        grid_set.adjoin_over_approximation(range,depth);
    } else {
        Box domain1(d),domain2(d);
        make_lpair(domain1,domain2)=split(domain);
        _adjoin_outer_approximation(set,domain1,eps,grid_set,depth);
        _adjoin_outer_approximation(set,domain2,eps,grid_set,depth);
    }    
}

 
void
adjoin_outer_approximation(GridTreeSet& grid_set, const TaylorSet& set, uint depth) 
{
    Box domain(set.number_of_generators(),Interval(-1,+1));
    Float eps=1.0/(1<<(depth/set.dimension()+1));
    return _adjoin_outer_approximation(set,domain,eps,grid_set,depth);
}


GridTreeSet 
outer_approximation(const TaylorSet& set, const Grid& grid, uint depth)
{
    ARIADNE_ASSERT(set.dimension()==grid.dimension());
    GridTreeSet grid_set(grid);
    adjoin_outer_approximation(grid_set,set,depth);
    return grid_set;
}

std::ostream&
TaylorSet::write(std::ostream& os) const
{
    os << "TaylorSet(\n";
    os << "  dimension=" << this->dimension() << ",\n" << std::flush;
    os << "  variables=" << this->_variables << ",\n" << std::flush;
    os << ")\n";
    return os;
}

Vector<Interval> evaluate(const TaylorSet& ts, const Vector<Interval>& iv) {
    Vector<Interval> r(ts.dimension());
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=ts[i].evaluate(iv);
    }
    return r;
}

Vector<Interval> evaluate(const TaylorSet& ts, const Vector<Float>& v) {
    Vector<Interval> r(ts.dimension());
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=ts[i].evaluate(v);
    }
    return r;
}

Vector<Interval> error(const TaylorSet& ts) {
    Vector<Interval> r(ts.dimension());
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=ts[i].error();
    }
    return r;
}


void draw(Figure& fig, const TaylorSet& ts) {
    static const double MAX_NEGLIGABLE_NORM=1e-10;
    if(ts.dimension()==2 && ts.number_of_generators()==2 && norm(error(ts))<MAX_NEGLIGABLE_NORM) {
        std::vector<Point> pts;
        Vector<Float> s(2); 
        s[0]=-1; s[1]=-1;
        for(uint i=0; i!=16; ++i) {
            pts.push_back(midpoint(evaluate(ts,s)));
            s[0]+=1./8;
        }
        for(uint i=0; i!=16; ++i) {
            pts.push_back(midpoint(evaluate(ts,s)));
            s[1]+=1./8;
        }
        for(uint i=0; i!=16; ++i) {
            pts.push_back(midpoint(evaluate(ts,s)));
            s[0]-=1./8;
        }
        for(uint i=0; i!=16; ++i) {
            pts.push_back(midpoint(evaluate(ts,s)));
            s[1]-=1./8;
        }
       fig.draw(pts);
    }
    else {
        draw(fig,zonotope(ts));
    }

}

} // namespace Ariadne
 
