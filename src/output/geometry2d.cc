/****************************************************************************
 *            geometry2d.cc
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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

#include "base/stlio.h"
#include "output/geometry2d.h"

namespace Ariadne { 

Output::Polygon2d& 
Output::Polygon2d::reduce() 
{
  //std::cerr << this->_vertices << std::endl;
  std::sort(this->_vertices.begin(),this->_vertices.end());
  //std::cerr << this->_vertices << std::endl;

  const std::vector<Point2d>& old_vertices=this->_vertices;
  std::vector<Point2d> new_vertices;
  
  // Sweep lower boundary from bottom-left to top right
  std::size_t min_size=1;
  for(std::vector<Point2d>::const_iterator vertex_iter=old_vertices.begin();
      vertex_iter!=old_vertices.end(); ++vertex_iter) 
  {
    const Point2d& vertex=*vertex_iter;
    while(new_vertices.size()>min_size) {
      const Point2d& penultimate=new_vertices[new_vertices.size()-2];
      const Point2d& last=new_vertices[new_vertices.size()-1];
      if(penultimate[0]==last[0]) {
        new_vertices.pop_back();
      } else if(last[0]==vertex[0]) {
        break;
      } else if(slope(penultimate,vertex) <= slope(penultimate,last)) {
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
  for(std::vector<Point2d>::const_reverse_iterator vertex_iter=old_vertices.rbegin();
      vertex_iter!=old_vertices.rend(); ++vertex_iter) 
  {
    const Point2d& vertex=*vertex_iter;
    while(new_vertices.size()>min_size) {
      const Point2d& penultimate=new_vertices[new_vertices.size()-2];
      const Point2d& last=new_vertices[new_vertices.size()-1];
      if(penultimate[0]==last[0]) {
        new_vertices.pop_back();
      } else if(last[0]==vertex[0]) {
        break;
      } else if(slope(penultimate,vertex) <= slope(penultimate,last)) {
        new_vertices.pop_back();
      } else {
        break;
      }
    }
    new_vertices.push_back(vertex);
  }
  new_vertices.pop_back();
  new_vertices.swap(this->_vertices);
  //std::cerr << this->_vertices << std::endl;
  return *this;
}




std::ostream& 
Output::operator<<(std::ostream& os, const Point2d& pt) 
{
  return os << "(" << pt[0] << "," << pt[1] << ")";
}

std::ostream& 
Output::operator<<(std::ostream& os, const Rectangle2d& r) 
{
  return os << "[" << r.lower_bound(0) << "," << r.upper_bound(0) << "]x["
            << r.lower_bound(1) << "," << r.upper_bound(1) << "]";
}

std::ostream& 
Output::operator<<(std::ostream& os, const PlanarProjectionMap& ppm) 
{
  return os << "PlanarProjectionMap( argument_dimension=" << ppm._d
            << ", x_variable=" << ppm._i << ", y_variable=" << ppm._j << " )";
}


Output::Point2d 
Output::Polygon2d::baricentre() const
{
  const Polygon2d& vertices=*this;
  Point2d baricentre(2);
  
  for (size_type j=0; j!=vertices.size(); ++j) {
    for (size_type i=0; i<2; i++) {
      baricentre[i]=baricentre[i]+vertices[j][i];
    }
  }
  for (size_type i=0; i!=2; ++i) {
    baricentre[i]/=vertices.size();
  }
  return baricentre;
}



Output::PlanarProjectionMap::PlanarProjectionMap()
  : _d(2), _i(0), _j(1) 
{ 
}


Output::PlanarProjectionMap::PlanarProjectionMap(dimension_type d, dimension_type i, dimension_type j)
  : _d(d), _i(i), _j(j) 
{ 
  if(i>=d || j>=d) { 
    using Geometry::InvalidCoordinate; 
    ARIADNE_THROW(InvalidCoordinate,"PlanarProjectionMap::PlanarProjectionMap(dimension_type d, dimension_type i, dimension_type j)",
                  "d="<<d<<", i="<<i<<", j="<<j); 
  }
}
  


} // namespace Ariadne
