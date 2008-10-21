/****************************************************************************
 *            geometry2d.cc
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
#include "list_set.h"
#include "grid_set.h"

#include "geometry2d.h"

namespace Ariadne { 

Box Polytope::bounding_box() const
{
  const Polytope& p=*this;
  Box res(p[0]);
  for(uint i=1; i!=p.size(); ++i) {
    res=hull(res,p[i]);
  }
  return res;
}


bool 
operator<(const Point& pt1, const Point& pt2)
{
  ARIADNE_ASSERT(pt1.dimension()==pt2.dimension());
  for(uint i=0; i!=pt1.dimension(); ++i) {
    if(pt1[i]<pt2[i]) {
      return true; 
    } else if(pt1[i]>pt2[i]) {
      return false;
    }
  }
  return false;
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
  //std::cerr<<__PRETTY_FUNCTION__<<std::endl;
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
  return reduce2d(res);
}


Polytope polytope(const Polytope& p) 
{
  Polytope r=p;
  return reduce2d(r);
}


Polytope& reduce2d(Polytope& p) 
{
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
  
  new_vertices.swap(p);
  return p;
}


Point 
baricentre(const Polytope& p) 
{
  const std::vector<Point>& vertices=p;
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






PlanarProjectionMap::PlanarProjectionMap()
  : _d(2), _i(0), _j(1) 
{ 
}


PlanarProjectionMap::PlanarProjectionMap(uint d, uint i, uint j)
  : _d(d), _i(i), _j(j) 
{ 
  if(i>=d || j>=d) { 
    ARIADNE_THROW(InvalidCoordinate,"PlanarProjectionMap::PlanarProjectionMap(uint d, uint i, uint j)",
                  "d="<<d<<", i="<<i<<", j="<<j); 
  }
}
  

Point
PlanarProjectionMap::operator()(const Point& pt) const 
{
  Point result(2);
  ARIADNE_ASSERT(pt.dimension()==this->_d);
  result[0]=pt[this->_i]; 
  result[1]=pt[this->_j]; 
  return result;
}


Box
PlanarProjectionMap::operator()(const Box& bx) const 
{
  Box result(2); 
  ARIADNE_ASSERT(bx.dimension()==this->_d);
  result[0].l = bx[this->_i].l;
  result[0].u = bx[this->_i].u;
  result[1].l = bx[this->_j].l;
  result[1].u = bx[this->_j].u;
  return result;
}


Zonotope 
PlanarProjectionMap::operator() (const Zonotope& z) const 
{
  const PlanarProjectionMap& map=*this;
  Point new_centre=map(z.centre());
  Matrix<Float> new_generators(2u,z.number_of_generators());
  for(size_t i=0; i!=z.number_of_generators(); ++i) {
    column(new_generators,i)=map(Vector<Float>(column(z.generators(),i)));
  }
  return Zonotope(new_centre,new_generators);
}


Polytope
PlanarProjectionMap::operator() (const Polytope& p) const 
{
  Polytope result(2);
  for(size_t i=0; i!=p.size(); ++i) {
    Point v=p[i];
    Point pt=(*this)(v);
    result.new_vertex(pt);
  }
  return reduce2d(result);
}


/*
InterpolatedCurve
PlanarProjectionMap::operator()(const InterpolatedCurve& curve) const
{
  const PlanarProjectionMap& self=*this;
  InterpolatedCurve::const_iterator iter=curve.begin();
  InterpolatedCurve::const_iterator end=curve.end();
  InterpolatedCurve result(self(iter->second));
  while(iter!=curve.end()) {
    ++iter;
    result.push_back(self(iter->second));
  }
  return result;
}
*/

std::ostream& 
operator<<(std::ostream& os, const PlanarProjectionMap& ppm) 
{
  return os << "PlanarProjectionMap( argument_dimension=" << ppm._d
            << ", x_variable=" << ppm._i << ", y_variable=" << ppm._j << " )";
}



} // namespace Ariadne
