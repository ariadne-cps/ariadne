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

#include "output/epsfstream.h"

namespace Ariadne { namespace Output {

const uint epsfstream::xBBoxSide=300;
const uint epsfstream::yBBoxSide=300;
const double epsfstream::linewidth=0.0000001;
    

Polygon2d& Polygon2d::reduce() 
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
operator<<(std::ostream& os, const Point2d& pt) 
{
  return os << "(" << pt[0] << "," << pt[1] << ")";
}

std::ostream& 
operator<<(std::ostream& os, const PlanarProjectionMap& ppm) 
{
  return os << "PlanarProjectionMap( argument_dimension=" << ppm._d
            << ", x_variable=" << ppm._i << ", y_variable=" << ppm._j << " )";
}


Point2d 
baricentre(const Polygon2d& vertices)
{
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



Polygon2d
order_around_a_point(const Polygon2d& vertices, 
                     const Point2d &centre)
{
  radiant_pointer_type pointers[vertices.size()];
  
  double tangent=0;
  double tangent_R=0;
  
  for (size_t i=0; i< vertices.size(); i++) {
    pointers[i].pos=i;
    if (vertices[i][1]==0) {
      if (vertices[i][0] >0)
        pointers[i].radiant=std::acos(0.0);
      else
        pointers[i].radiant=-std::acos(0.0);
    } else {
      tangent_R=vertices[i][0]/vertices[i][1];
      tangent=tangent_R;
      pointers[i].radiant=std::atan(tangent);
      
      if (vertices[i][1] <0)
        pointers[i].radiant+=std::acos(-1.0);
    }
  }
  
  std::sort(pointers, pointers + vertices.size(), is_smaller_than);
  
  Polygon2d new_vector(2,vertices.size()); 
  
  for (size_t i=0; i< vertices.size(); i++) {
    Point2d pt(vertices[pointers[i].pos]);
    new_vector[i]=pt;
  }
  
  return new_vector;
}
    





epsfstream::epsfstream()
  : std::ofstream(), line_colour("black"), fill_colour("green"), line_style(true), fill_style(true)
{
}



epsfstream::~epsfstream() {
  this->close();
}

     
 
void
epsfstream::open(const char* fn, const Rectangle2d& bbox, 
                 const PlanarProjectionMap& p_map)
{
  this->p_map=p_map;
  this->bbox=bbox;

  this->std::ofstream::open(fn);
  this->write_header();
}



void 
epsfstream::close() 
{
  this->write_trailer();
  this->std::ofstream::close();
}



void 
epsfstream::write_header() 
{
  (*this) << "%!PS-Adobe-2.0\n"
          << "%%Title: Ariadne\n"
          << "%%Creator: Ariadne\n"
          << "%%CreationDate: Unknown\n"
          << "%%For: Test\n"
          << "%%Pages: 1\n"
          << "%%DocumentFonts: /Times-Roman\n"
          << "%%BoundingBox: 0 0 " << xBBoxSide+10 << " " << yBBoxSide+10 << "\n"
          << "%%EndComments\n"
          << "\n";
  
  (*this) << "/linewidth " << epsfstream::linewidth << " def\n";
  (*this) << "/black {0.0 0.0 0.0 setrgbcolor} def\n"
          << "/gray {0.5 0.5 0.5 setrgbcolor} def\n"
          << "/white {1.0 1.0 1.0 setrgbcolor} def\n"
          << "/red {1.0 0.0 0.0 setrgbcolor} def\n"
          << "/green {0.0 1.0 0.0 setrgbcolor} def\n"
          << "/blue {0.0 0.0 1.0 setrgbcolor} def\n"
          << "/yellow {1.0 1.0 0.0 setrgbcolor} def\n"
          << "/magenta {1.0 0.0 1.0 setrgbcolor} def\n"
          << "/cyan {0.0 1.0 1.0 setrgbcolor} def\n";
  
  (*this) << "/bordercolour { black } def\n"
          << "/fillcolour { green } def\n";
  
  (*this) << "/xl " << this->bbox.lower_bound(0) << " def\n"
          << "/xu " << this->bbox.upper_bound(0) << " def\n"
          << "/yl " << this->bbox.lower_bound(1) << " def\n"
          << "/yu " << this->bbox.upper_bound(1) << " def\n";
  (*this) << "/xBBoxSide " << xBBoxSide << " def\n"
          << "/yBBoxSide " << yBBoxSide << " def\n"
          << "/xscale xBBoxSide xu xl sub div def  % horizontal scaling factor\n"
          << "/yscale yBBoxSide yu yl sub div def  % vertical scaling factor\n"
          << "/xscalei 1 xscale div def     % horizontal scaling factor inverse\n"
          << "/yscalei 1 yscale div def     % vertical scaling factor inverse\n";
  
  (*this) << "linewidth setlinewidth\n";
  (*this) << "green\n";
  
  (*this) << "gsave\n"
          << "xscale yscale scale xl neg yl neg translate\n";
}



void 
epsfstream::write_trailer() 
{
  (*this) << "grestore showpage\n"
    "%%Trailer\n";
}



void
epsfstream::trace_scale(const char* x_name, const char* y_name,
                        const int& x_step, const int& y_step) 
{
  double lx=this->bbox.lower_bound(0);
  double ly=this->bbox.lower_bound(1);
  double ux=this->bbox.upper_bound(0);
  double uy=this->bbox.upper_bound(1);
  
  double scale_x=((ux-lx)/(SCALE_DIMENSION + 1));
  double scale_y=((uy-ly)/(SCALE_DIMENSION + 1));
  
  lx=lx+0.8*scale_x;
  ly=ly+0.8*scale_y;
  
  double setlinewidth = ((ux-lx)+(uy-ly))/(SCALE_DIMENSION*300);
  double fontsize = ((ux-lx)+(uy-ly))/(SCALE_DIMENSION*20);
  
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



epsfstream&
trace(epsfstream& eps, const Point2d& pt)
{        
  eps << pt[0] << " " << pt[1] << " 0.001 0 360 arc closepath\n";
  return eps;
}

epsfstream&
trace(epsfstream& eps, const Rectangle2d& r)
{
  double lx=r.lower_bound(0);
  double ux=r.upper_bound(0);
  double ly=r.lower_bound(1);
  double uy=r.upper_bound(1);
  
  eps << lx << ' ' << ly << " moveto\n"
      << ux << ' ' << ly << " lineto\n"
      << ux << ' ' << uy << " lineto\n"
      << lx << ' ' << uy << " lineto\n"
      << lx << ' ' << ly << " lineto\n";
  return eps;
}

epsfstream&
trace(epsfstream& eps, const Polygon2d& vertices)
{
  eps << vertices[0][0] << ' ' << vertices[0][1] 
      << " moveto\n";
  for (size_type i=1; i!=vertices.size(); ++i) {
    eps << vertices[i][0] << ' ' << vertices[i][1] 
        << " lineto\n";
  }
  eps << vertices[0][0] << ' ' << vertices[0][1] 
      << " lineto\n";
  return eps;
}



epsfstream&
draw(epsfstream& eps, const Point2d& pt)
{
  if(eps.fill_style) {
    eps << pt[0] << " " << pt[1] << " 0.001 0 360 arc closepath\n";
    eps << eps.fill_colour << " fill\n";
  }
  if(eps.line_style) {
    eps << pt[0] << " " << pt[1] << " 0.001 0 360 arc closepath\n";
    eps << eps.line_colour << " stroke\n\n";
  }
  return eps;
}

epsfstream&
draw(epsfstream& eps, const Rectangle2d& r) 
{
  if(eps.fill_style) {
    trace(eps,r) << eps.fill_colour << " fill\n";
  }
  if(eps.line_style) {
    trace(eps,r) << eps.line_colour << " stroke\n\n";
  }
  return eps;
}

epsfstream&
draw(epsfstream& eps, const Polygon2d& vertices)
{
  if(eps.fill_style) {
    trace(eps,vertices) << eps.fill_colour << " fill\n";
  }
  if(eps.line_style) {
    trace(eps,vertices) << eps.line_colour << " stroke\n\n";
  }
  return eps;
}



}}
