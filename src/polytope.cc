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

#include "functional.h"
#include "config.h"

#include "stlio.h"
#include "vector.h"
#include "matrix.h"
#include "point.h"
#include "box.h"
#include "zonotope.h"
#include "curve.h"
#include "function.h"

#include "polytope.h"

namespace Ariadne {

tribool
Polytope::separated(const Box& bx) const {
    return this->bounding_box().separated(bx) || indeterminate;
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
Polytope::write(std::ostream& os) const
{
    return os << "Polytope( vertices=" << this->vertices() << " )";
}






Float
slope2d(const ExactPoint& pt1, const ExactPoint& pt2)
{
    ARIADNE_ASSERT(pt1.dimension()==2 && pt2.dimension()==2);
    return (pt2[1]-pt1[1])/(pt2[0]-pt1[0]);
}



Polytope polytope(const Box& bx)
{
    return Polytope(bx.vertices());
}


Polytope polytope(const Zonotope& z)
{
    //std::cerr<<ARIADNE_PRETTY_FUNCTION<<std::endl;
    Polytope res;
    size_t ng=z.generators().column_size();
    size_t nv=1<<ng;
    for(size_t i=0; i!=nv; ++i) {
        ExactPoint pt=z.centre();
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
    if(r.dimension()==2) {
        return reduce2d(r);
    } else {
        return r;
    }
}


Polytope& reduce2d(Polytope& p)
{
    ARIADNE_ASSERT(p.dimension()==2);
    //std::cerr << this->_vertices << std::endl;
    std::sort(p.vertices().begin(),p.vertices().end());
    //std::cerr << this->_vertices << std::endl;

    const std::vector<ExactPoint>& old_vertices=p.vertices();
    std::vector<ExactPoint> new_vertices;

    // Sweep lower boundary from bottom-left to top right
    size_t min_size=1;
    for(std::vector<ExactPoint>::const_iterator vertex_iter=old_vertices.begin();
        vertex_iter!=old_vertices.end(); ++vertex_iter)
        {
            const ExactPoint& vertex=*vertex_iter;
            while(new_vertices.size()>min_size) {
                const ExactPoint& penultimate=new_vertices[new_vertices.size()-2];
                const ExactPoint& last=new_vertices[new_vertices.size()-1];
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
    for(std::vector<ExactPoint>::const_reverse_iterator vertex_iter=old_vertices.rbegin();
        vertex_iter!=old_vertices.rend(); ++vertex_iter)
        {
            const ExactPoint& vertex=*vertex_iter;
            while(new_vertices.size()>min_size) {
                const ExactPoint& penultimate=new_vertices[new_vertices.size()-2];
                const ExactPoint& last=new_vertices[new_vertices.size()-1];
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


ExactPoint
baricentre(const Polytope& p)
{
    const std::vector<ExactPoint>& vertices=p.vertices();
    ExactPoint baricentre(p.dimension());

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









void Polytope::draw(CanvasInterface& c, const Projection2d& p) const {
    uint xi=p.x_coordinate(); uint yi=p.y_coordinate();
    ARIADNE_ASSERT(max(xi,yi)<this->dimension());

    Polytope pr(2);
    ExactPoint prv(2);

    // Project polytope onto canvas coordinates
    for(uint i=0; i!=this->number_of_vertices(); ++i) {
        const ExactPoint& v=this->vertex(i);
        prv[0]=v[xi]; prv[1]=v[yi];
        pr.new_vertex(prv);
    }

    // Reduce boundary of projected polytope
    reduce2d(pr);

    // Trace boundary
    prv=pr.vertex(pr.number_of_vertices()-1);
    c.move_to(prv[0],prv[1]);
    for(uint i=0; i!=pr.number_of_vertices(); ++i) {
        prv=pr.vertex(i);
        c.line_to(prv[0],prv[1]);
    }

    c.fill();
}


} // namespace Ariadne
