/****************************************************************************
 *            polytope.cc
 *
 *  Copyright  2007-8  Alberto Casagrande, Pieter Collins
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

#include "stlio.h"
#include "stlio.h"
#include "vector.h"
#include "matrix.h"
#include "point.h"
#include "box.h"
#include "zonotope.h"
#include "curve.h"
#include "function.h"

#include "taylor_set.h"

#include "polytope.h"

namespace Ariadne { 

Zonotope zonotope(const TaylorSet& ts);

tribool 
Polytope::disjoint(const Box& bx) const {
    return this->bounding_box().disjoint(bx) || indeterminate;
}

tribool 
Polytope::overlaps(const Box& bx) const {
    return bx.covers(baricentre(*this)) || indeterminate;
}

tribool 
Polytope::inside(const Box& bx) const {
    return this->bounding_box().inside(bx) || indeterminate;
}

Box 
Polytope::bounding_box() const
{
    const Polytope& p=*this;
    Box res(p._vertices[0]);
    for(uint i=1; i!=p._vertices.size(); ++i) {
        res=hull(res,p._vertices[i]);
    }
    return res;
}


std::ostream& 
operator<<(std::ostream& os, const Polytope& p)
{
    return os << "Polytope( vertices=" << p.vertices() << " )";
}






Float 
slope2d(const Point& pt1, const Point& pt2) 
{
    ARIADNE_ASSERT(pt1.dimension()==2 && pt2.dimension()==2);
    return (pt2[1]-pt1[1])/(pt2[0]-pt1[0]);
}



Polytope polytope(const Box& bx) 
{
    ARIADNE_ASSERT(bx.dimension()==2);
    Polytope result(2);
    Point pt(2);
    pt[0]=bx[0].lower();
    pt[1]=bx[1].lower();
    result.new_vertex(pt);
    pt[0]=bx[0].upper();
    result.new_vertex(pt);
    pt[1]=bx[1].upper();
    result.new_vertex(pt);
    pt[0]=bx[0].lower();
    result.new_vertex(pt);
    return result;
}


Polytope polytope(const Zonotope& z) 
{
    //std::cerr<<ARIADNE_PRETTY_FUNCTION<<std::endl;
    Polytope res;
    size_t ng=z.generators().column_size();
    size_t nv=1<<ng;
    for(size_t i=0; i!=nv; ++i) {
        Point pt=z.centre();
        for(size_t j=0; j!=ng; ++j) {
            if(i & 1<<j) {
                pt+=column(z.generators(),j);
            } else {
                pt-=column(z.generators(),j);
            }
        }
        res.new_vertex(pt);
    }
    if(res.dimension()==2) {
        return reduce2d(res);
    } else {
        return res;
    }
}


Polytope polytope(const Polytope& p) 
{
    Polytope r=p;
    return reduce2d(r);
}


Polytope polytope(const TaylorSet& ts) 
{
    return polytope(zonotope(ts)); 
}


Polytope& reduce2d(Polytope& p) 
{
    ARIADNE_ASSERT(p.dimension()==2);
    //std::cerr << this->_vertices << std::endl;
    std::sort(p.vertices().begin(),p.vertices().end());
    //std::cerr << this->_vertices << std::endl;

    const std::vector<Point>& old_vertices=p.vertices();
    std::vector<Point> new_vertices;
  
    // Sweep lower boundary from bottom-left to top right
    size_t min_size=1;
    for(std::vector<Point>::const_iterator vertex_iter=old_vertices.begin();
        vertex_iter!=old_vertices.end(); ++vertex_iter) 
        {
            const Point& vertex=*vertex_iter;
            while(new_vertices.size()>min_size) {
                const Point& penultimate=new_vertices[new_vertices.size()-2];
                const Point& last=new_vertices[new_vertices.size()-1];
                if(penultimate[0]==last[0]) {
                    new_vertices.pop_back();
                } else if(last[0]==vertex[0]) {
                    break;
                } else if(slope2d(penultimate,vertex) <= slope2d(penultimate,last)) {
                    new_vertices.pop_back();
                } else {
                    break;
                }
            }
            new_vertices.push_back(vertex);
        }
    new_vertices.pop_back();
  
    // Upper sweep
    //std::cerr << "Forward pass\n";
    min_size=new_vertices.size()+1;
    for(std::vector<Point>::const_reverse_iterator vertex_iter=old_vertices.rbegin();
        vertex_iter!=old_vertices.rend(); ++vertex_iter) 
        {
            const Point& vertex=*vertex_iter;
            while(new_vertices.size()>min_size) {
                const Point& penultimate=new_vertices[new_vertices.size()-2];
                const Point& last=new_vertices[new_vertices.size()-1];
                if(penultimate[0]==last[0]) {
                    new_vertices.pop_back();
                } else if(last[0]==vertex[0]) {
                    break;
                } else if(slope2d(penultimate,vertex) <= slope2d(penultimate,last)) {
                    new_vertices.pop_back();
                } else {
                    break;
                }
            }
            new_vertices.push_back(vertex);
        }
    new_vertices.pop_back();
  
    new_vertices.swap(p.vertices());
    return p;
}


Point 
baricentre(const Polytope& p) 
{
    const std::vector<Point>& vertices=p.vertices();
    Point baricentre(p.dimension());
  
    for (size_t j=0; j!=vertices.size(); ++j) {
        for (size_t i=0; i<2; i++) {
            baricentre[i]=baricentre[i]+vertices[j][i];
        }
    }
    for (size_t i=0; i!=2; ++i) {
        baricentre[i]/=vertices.size();
    }
    return baricentre;
}








Vector<Float>
apply(const ProjectionFunction& map, const Vector<Float>& v)  
{
    ARIADNE_ASSERT(v.size()==map.argument_size());
    Vector<Float> result(map.result_size());
    for(uint i=0; i!=map.result_size(); ++i) {
        result[i]=v[map[i]];
    }
    return result;
}

Point
apply(const ProjectionFunction& map, const Point& pt)  
{
    ARIADNE_ASSERT(pt.dimension()==map.argument_size());
    Point result(map.result_size());
    for(uint i=0; i!=map.result_size(); ++i) {
        result[i]=pt[map[i]];
    }
    return result;
}


Box
apply(const ProjectionFunction& map, const Box& bx)  
{
    ARIADNE_ASSERT(bx.dimension()==map.argument_size());
    Box result(map.result_size()); 
    for(uint i=0; i!=map.result_size(); ++i) {
        result[i].l = bx[map[i]].l;
        result[i].u = bx[map[i]].u;
    }
    return result;
}


Zonotope 
apply(const ProjectionFunction& map, const Zonotope& z)  
{
    Point new_centre=apply(map,z.centre());
    Matrix<Float> new_generators(2u,z.number_of_generators());
    for(size_t i=0; i!=z.number_of_generators(); ++i) {
        column(new_generators,i)=apply(map,Vector<Float>(column(z.generators(),i)));
    }
    return Zonotope(new_centre,new_generators);
}


Polytope
apply(const ProjectionFunction& map, const Polytope& p)  
{
    Polytope result(map.result_size());
    for(size_t i=0; i!=p.number_of_vertices(); ++i) {
        Point const& v=p.vertex(i);
        Point pt=apply(map,v);
        result.new_vertex(pt);
    }
    if(result.dimension()==2) {
        return reduce2d(result); 
    } else {
        return result;
    }
}


InterpolatedCurve
apply(const ProjectionFunction& map, const InterpolatedCurve& curve) 
{
    InterpolatedCurve::const_iterator iter=curve.begin();
    InterpolatedCurve::const_iterator end=curve.end();
    InterpolatedCurve result(apply(map,iter->second));
    while(iter!=curve.end()) {
        ++iter;
        result.insert(iter->first,apply(map,iter->second));
    }
    return result;
}




} // namespace Ariadne
