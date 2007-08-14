/****************************************************************************
 *            epsfstream.cc
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
#include "output/epsfstream.h"

namespace Ariadne { 


const uint Output::epsfstream::xBBoxSide=300;
const uint Output::epsfstream::yBBoxSide=300;
const double Output::epsfstream::linewidth=0.0000001;
const double Output::epsfstream::scale_dimension=3.5;
    

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
  





Output::epsfstream::epsfstream()
  : std::ofstream(), line_colour(black), fill_colour(green), line_style(true), fill_style(true)
{
}



Output::epsfstream::~epsfstream() {
  this->close();
}

     
 
void
Output::epsfstream::open(const char* fn, const Rectangle2d& bbox, 
                         const PlanarProjectionMap& p_map)
{
  this->p_map=p_map;
  this->bbox=bbox;

  this->std::ofstream::open(fn);
  this->write_header();
  static_cast<std::ostream&>(*this) << "black\n";
  this->trace(bbox);
  static_cast<std::ostream&>(*this) << " stroke\n\n";
}



void 
Output::epsfstream::close() 
{
  this->write_trailer();
  this->std::ofstream::close();
}



void 
Output::epsfstream::write_header() 
{
  std::ostream& os=*this;
  os      << "%!PS-Adobe-2.0\n"
          << "%%Title: Ariadne\n"
          << "%%Creator: Ariadne\n"
          << "%%CreationDate: Unknown\n"
          << "%%For: Test\n"
          << "%%Pages: 1\n"
          << "%%DocumentFonts: /Times-Roman\n"
          << "%%BoundingBox: -10 -10 " << xBBoxSide+10 << " " << yBBoxSide+10 << "\n"
          << "%%EndComments\n"
          << "\n";
  
  os      << "/linewidth " << epsfstream::linewidth << " def\n";
  os      << "/black {0.0 0.0 0.0 setrgbcolor} def\n"
          << "/gray {0.5 0.5 0.5 setrgbcolor} def\n"
          << "/white {1.0 1.0 1.0 setrgbcolor} def\n"
          << "/red {1.0 0.0 0.0 setrgbcolor} def\n"
          << "/green {0.0 1.0 0.0 setrgbcolor} def\n"
          << "/blue {0.0 0.0 1.0 setrgbcolor} def\n"
          << "/yellow {1.0 1.0 0.0 setrgbcolor} def\n"
          << "/magenta {1.0 0.0 1.0 setrgbcolor} def\n"
          << "/cyan {0.0 1.0 1.0 setrgbcolor} def\n";
  
  os      << "/bordercolour { black } def\n"
          << "/fillcolour { green } def\n";
  
  os      << "/xl " << this->bbox.lower_bound(0) << " def\n"
          << "/xu " << this->bbox.upper_bound(0) << " def\n"
          << "/yl " << this->bbox.lower_bound(1) << " def\n"
          << "/yu " << this->bbox.upper_bound(1) << " def\n";
  os      << "/xBBoxSide " << xBBoxSide << " def\n"
          << "/yBBoxSide " << yBBoxSide << " def\n"
          << "/xscale xBBoxSide xu xl sub div def  % horizontal scaling factor\n"
          << "/yscale yBBoxSide yu yl sub div def  % vertical scaling factor\n"
          << "/xscalei 1 xscale div def     % horizontal scaling factor inverse\n"
          << "/yscalei 1 yscale div def     % vertical scaling factor inverse\n";
  
  os      << "linewidth setlinewidth\n";
  os      << "green\n";
  
  os      << "gsave\n"
          << "xscale yscale scale xl neg yl neg translate\n";
}



void 
Output::epsfstream::write_trailer() 
{
  std::ostream& os=*this;
  os << "grestore showpage\n"
    "%%Trailer\n";
}



void
Output::epsfstream::trace_scale(const char* x_name, const char* y_name,
                        const int& x_step, const int& y_step) 
{
  double lx=this->bbox.lower_bound(0);
  double ly=this->bbox.lower_bound(1);
  double ux=this->bbox.upper_bound(0);
  double uy=this->bbox.upper_bound(1);
  
  double scale_x=((ux-lx)/(scale_dimension + 1));
  double scale_y=((uy-ly)/(scale_dimension + 1));
  
  lx=lx+0.8*scale_x;
  ly=ly+0.8*scale_y;
  
  double setlinewidth = ((ux-lx)+(uy-ly))/(scale_dimension*300);
  double fontsize = ((ux-lx)+(uy-ly))/(scale_dimension*20);
  
  *this << setlinewidth << " setlinewidth"<< std::endl
        << "fillcolour black" << std::endl
        << "/Times-Roman findfont" << std::endl
        << fontsize  << " scalefont" << std::endl
        << "setfont" << std::endl;
  
  *this << "newpath" << std::endl
        << lx << ' ' << ly << " moveto"<< std::endl
        << ux << ' ' << ly << " lineto"<< std::endl
        << "stroke"<< std::endl;
  
  size_t prec=this->precision();
  
  if (x_step>0) {
    for (double i=lx; i< ux; i=i+((ux-lx)/(x_step))) {
      *this << "newpath" << std::endl
            << i << ' ' << ly - 0.05*scale_y << " moveto"<< std::endl
            << i << ' ' << ly + 0.05*scale_y << " lineto"<< std::endl
            << "stroke"<< std::endl;
      
      *this << "newpath" << std::endl
            << i-0.05*scale_x << ' ' << ly - 0.2*scale_y 
            << " moveto"<< std::endl;
      
      this->precision(1);
      *this << "("<< i << ") show"<< std::endl;
      this->precision(prec);
    }
    
    *this << "newpath" << std::endl
          << ux << ' ' << ly - 0.05*scale_y << " moveto"<< std::endl
          << ux << ' ' << ly + 0.05*scale_y << " lineto"<< std::endl
          << "stroke"<< std::endl;
    
    *this << "newpath" << std::endl
          << ux-0.05*scale_x << ' ' << ly - 0.2*scale_y 
          << " moveto"<< std::endl;
    
    this->precision(1);
    *this << "("<< ux << ") show"<< std::endl;
    this->precision(prec);
    
    double x_name_len=((double)strlen(x_name))/28;
    
    *this << "/Times-Roman findfont" << std::endl
          << 1.5* fontsize  << " scalefont" << std::endl
          << "setfont" << std::endl
          << "newpath" << std::endl
          << (ux+lx-x_name_len)/2 << ' ' << ly - 0.5*scale_y
          << " moveto"<< std::endl
          << "("<< x_name << ") show"<< std::endl
          << "/Times-Roman findfont" << std::endl
          << fontsize  << " scalefont" << std::endl
          << "setfont" << std::endl;
  }
  
  *this << "newpath" << std::endl
        << lx << ' ' << ly << " moveto"<< std::endl
        << lx << ' ' << uy << " lineto"<< std::endl
        << "stroke"<< std::endl;
  
  if (y_step>0) {
    for (double i=ly; i<= uy; i=i+((uy-ly)/(y_step))) {
      
      *this << "newpath" << std::endl
            << lx - 0.05*scale_x << ' ' << i << " moveto"<< std::endl
            << lx + 0.05*scale_x << ' ' << i << " lineto"<< std::endl
            << "stroke"<< std::endl;
      
      *this << "newpath" << std::endl
            << lx - 0.30*scale_x << ' ' << i-0.05*scale_y 
            << " moveto"<< std::endl;
      
      this->precision(1);
      *this << "("<< i << ") show"<< std::endl;
      this->precision(prec);
    }
    
    *this << "newpath" << std::endl
          << lx - 0.05*scale_x << ' ' << uy << " moveto"<< std::endl
          << lx + 0.05*scale_x << ' ' << uy << " lineto"<< std::endl
          << "stroke"<< std::endl;
    
    *this << "newpath" << std::endl
          << lx - 0.30*scale_x << ' ' << uy-0.05*scale_y 
          << " moveto"<< std::endl;
    
    this->precision(1);
    *this << "("<< uy << ") show"<< std::endl;
    this->precision(prec);
    
    double y_name_len=((double)strlen(y_name))/28;
    
    *this << "/Times-Roman findfont" << std::endl
          << 1.5* fontsize  << " scalefont" << std::endl
          << "setfont" << std::endl
          << "newpath" << std::endl
          << lx - 0.5*scale_x << ' ' << (uy+ly-y_name_len)/2 
          << " moveto"<< std::endl
          << "90 rotate" << std::endl
          << "("<< y_name << ") show"<< std::endl
          << "/Times-Roman findfont" << std::endl
          << fontsize  << " scalefont" << std::endl
          << "setfont" << std::endl
          << "stroke" << std::endl
          << "-90 rotate" << std::endl;
  }
  
  *this << "linewidth setlinewidth" << std::endl
        << "fillcolour green" << std::endl;
}



void
Output::epsfstream::trace(const Point2d& pt)
{        
  std::ostream& os=*this;
  os << pt[0] << " " << pt[1] << " 0.001 0 360 arc closepath\n";
}

void
Output::epsfstream::trace(const Rectangle2d& r)
{
  double lx=r.lower_bound(0);
  double ux=r.upper_bound(0);
  double ly=r.lower_bound(1);
  double uy=r.upper_bound(1);
  
  std::ostream& os=*this;
  os << lx << ' ' << ly << " moveto\n"
     << ux << ' ' << ly << " lineto\n"
     << ux << ' ' << uy << " lineto\n"
     << lx << ' ' << uy << " lineto\n"
     << lx << ' ' << ly << " lineto\n";
}

void
Output::epsfstream::trace(const Polygon2d& vertices)
{
  std::ostream& os=*this;
  os << vertices[0][0] << ' ' << vertices[0][1] 
     << " moveto\n";
  for (size_type i=1; i!=vertices.size(); ++i) {
    os << vertices[i][0] << ' ' << vertices[i][1] 
       << " lineto\n";
  }
  os << vertices[0][0] << ' ' << vertices[0][1] 
     << " lineto\n";
}



void
Output::epsfstream::draw(const Point2d& pt)
{
  std::ostream& os=*this;
  epsfstream& eps=*this;
  if(eps.fill_style) {
    os << pt[0] << " " << pt[1] << " 0.001 0 360 arc closepath\n";
    os << eps.fill_colour.name() << " fill\n";
  }
  if(eps.line_style) {
    os << pt[0] << " " << pt[1] << " 0.001 0 360 arc closepath\n";
    os << eps.line_colour.name() << " stroke\n\n";
  }
}

void
Output::epsfstream::draw(const Rectangle2d& r) 
{
  epsfstream& eps=*this;
  if(eps.fill_style) {
    eps.trace(r);
    eps.fill();
  }
  if(eps.line_style) {
    eps.trace(r);
    eps.stroke();
  }
}

void
Output::epsfstream::draw(const Polygon2d& vertices)
{
  epsfstream& eps=*this;
  if(eps.fill_style) {
    eps.trace(vertices);
    eps.fill();
  }
  if(eps.line_style) {
    eps.trace(vertices);
    eps.stroke();
  }
}


void
Output::epsfstream::draw(std::vector<Rectangle2d>& rl)
{
  epsfstream& eps=*this;
  std::sort(rl.begin(),rl.end());
  rl.resize(std::unique(rl.begin(),rl.end())-rl.begin());
  if(eps.fill_style) {
    for(uint i=0; i!=rl.size(); ++i) {
      eps.trace(rl[i]);
      eps.fill();
    }
  }
  if(eps.line_style) {
    for(uint i=0; i!=rl.size(); ++i) {
      eps.trace(rl[i]);
      eps.stroke();
    }
  }
}


void 
Output::epsfstream::fill()
{
  std::ostream& os=*this;
  os << this->fill_colour.name() << " fill\n";
}


void 
Output::epsfstream::stroke()
{
  std::ostream& os=*this;
  os << this->line_colour.name() << " stroke\n";
}


}
