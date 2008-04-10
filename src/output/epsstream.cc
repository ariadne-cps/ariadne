/****************************************************************************
 *            epsstream.cc
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
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

#include "base/stlio.h"
#include "output/epsstream.h"

namespace Ariadne { 


const uint Output::epsstream::xBBoxSide=300;
const uint Output::epsstream::yBBoxSide=300;
const double Output::epsstream::linewidth=0.0000001;
const double Output::epsstream::scale_dimension=3.5;
    


void 
Output::epsfstream::open(const char* fn, 
                         const Rectangle2d& bbox, 
                         const PlanarProjectionMap& p_map)
{
  this->_ofs_ptr->open(fn); 

  this->set_bounding_box(bbox);
  this->set_projection_map(p_map);

  this->write_header();
  this->ostream() << "black\n";
  this->trace(bbox);
  this->ostream() << "clip stroke\n\n";
}




Output::epsstream::epsstream()
  : _os_ptr(&std::cout), line_colour(black), fill_colour(green), line_style(true), fill_style(true)
{
}

Output::epsstream::epsstream(std::ostream& os)
  : _os_ptr(&os), line_colour(black), fill_colour(green), line_style(true), fill_style(true)
{
}

void
Output::epsstream::redirect(std::ostream& os)
{
  this->_os_ptr=&os;
}



Output::epsstream::~epsstream() {
}

     
 






void 
Output::epsstream::write_header() 
{
  std::ostream& os=this->ostream();
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
  
  os      << "/linewidth " << epsstream::linewidth << " def\n";
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
Output::epsstream::write_trailer() 
{
  std::ostream& os=this->ostream();

  // redraw the bounding box
  os << "black\n";
  this->trace(this->bounding_box());
  os << "stroke\n";

  os << "grestore showpage\n"
    "%%Trailer\n";
}



void
Output::epsstream::trace_scale(const char* x_name, const char* y_name,
                        const int& x_step, const int& y_step) 
{
  std::ostream& os=this->ostream();

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
  
  os << setlinewidth << " setlinewidth"<< std::endl
     << "fillcolour black" << std::endl
     << "/Times-Roman findfont" << std::endl
     << fontsize  << " scalefont" << std::endl
     << "setfont" << std::endl;
  
  os << "newpath" << std::endl
     << lx << ' ' << ly << " moveto"<< std::endl
     << ux << ' ' << ly << " lineto"<< std::endl
     << "stroke"<< std::endl;
  
  size_t prec=os.precision();
  
  if (x_step>0) {
    for (double i=lx; i< ux; i=i+((ux-lx)/(x_step))) {
      os << "newpath" << std::endl
         << i << ' ' << ly - 0.05*scale_y << " moveto"<< std::endl
         << i << ' ' << ly + 0.05*scale_y << " lineto"<< std::endl
         << "stroke"<< std::endl;
      
      os << "newpath" << std::endl
         << i-0.05*scale_x << ' ' << ly - 0.2*scale_y 
         << " moveto"<< std::endl;
      
      os.precision(1);
      os << "("<< i << ") show"<< std::endl;
      os.precision(prec);
    }
    
    os << "newpath" << std::endl
       << ux << ' ' << ly - 0.05*scale_y << " moveto"<< std::endl
       << ux << ' ' << ly + 0.05*scale_y << " lineto"<< std::endl
       << "stroke"<< std::endl;
    
    os << "newpath" << std::endl
       << ux-0.05*scale_x << ' ' << ly - 0.2*scale_y 
       << " moveto"<< std::endl;
    
    os.precision(1);
    os << "("<< ux << ") show"<< std::endl;
    os.precision(prec);
    
    double x_name_len=((double)strlen(x_name))/28;
    
    os << "/Times-Roman findfont" << std::endl
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
  
  os << "newpath" << std::endl
     << lx << ' ' << ly << " moveto"<< std::endl
     << lx << ' ' << uy << " lineto"<< std::endl
     << "stroke"<< std::endl;
  
  if (y_step>0) {
    for (double i=ly; i<= uy; i=i+((uy-ly)/(y_step))) {
      
      os << "newpath" << std::endl
         << lx - 0.05*scale_x << ' ' << i << " moveto"<< std::endl
         << lx + 0.05*scale_x << ' ' << i << " lineto"<< std::endl
         << "stroke"<< std::endl;
      
      os << "newpath" << std::endl
         << lx - 0.30*scale_x << ' ' << i-0.05*scale_y 
         << " moveto"<< std::endl;
      
      os.precision(1);
      os << "("<< i << ") show"<< std::endl;
      os.precision(prec);
    }
    
    os << "newpath" << std::endl
       << lx - 0.05*scale_x << ' ' << uy << " moveto"<< std::endl
       << lx + 0.05*scale_x << ' ' << uy << " lineto"<< std::endl
       << "stroke"<< std::endl;
    
    os << "newpath" << std::endl
       << lx - 0.30*scale_x << ' ' << uy-0.05*scale_y 
       << " moveto"<< std::endl;
    
    os.precision(1);
    os << "("<< uy << ") show"<< std::endl;
    os.precision(prec);
    
    double y_name_len=((double)strlen(y_name))/28;
    
    os << "/Times-Roman findfont" << std::endl
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
  
  os << "linewidth setlinewidth" << std::endl
     << "fillcolour green" << std::endl;
}



void
Output::epsstream::trace(const Point2d& pt)
{        
  std::ostream& os=this->ostream();
  os << pt[0] << " " << pt[1] << " 0.001 0 360 arc closepath\n";
}

void
Output::epsstream::trace(const Rectangle2d& r)
{
  double lx=r.lower_bound(0);
  double ux=r.upper_bound(0);
  double ly=r.lower_bound(1);
  double uy=r.upper_bound(1);
  
  std::ostream& os=this->ostream();
  os << lx << ' ' << ly << " moveto\n"
     << ux << ' ' << ly << " lineto\n"
     << ux << ' ' << uy << " lineto\n"
     << lx << ' ' << uy << " lineto\n"
    //<< lx << ' ' << ly << " lineto\n";
     << "closepath\n";
}

void
Output::epsstream::trace(const Zonotope2d& z)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

void
Output::epsstream::trace(const Polygon2d& vertices)
{
  std::ostream& os=this->ostream();
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
Output::epsstream::draw(const Point2d& pt)
{
  std::ostream& os=this->ostream();
  epsstream& eps=*this;
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
Output::epsstream::draw(const Rectangle2d& r) 
{
  epsstream& eps=*this;
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
Output::epsstream::draw(const Zonotope2d& z)
{
  draw(Polygon2d(z));
}

void
Output::epsstream::draw(const Polygon2d& p)
{
  epsstream& eps=*this;
  if(eps.fill_style) {
    eps.trace(p);
    eps.fill();
  }
  if(eps.line_style) {
    eps.trace(p);
    eps.stroke();
  }
}


void
Output::epsstream::draw(std::vector<Rectangle2d>& rl)
{
  //std::cerr<<__PRETTY_FUNCTION__<<std::endl;
  epsstream& eps=*this;
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
Output::epsstream::fill()
{
  std::ostream& os=this->ostream();
  os << this->fill_colour.name() << " fill\n";
}


void 
Output::epsstream::stroke()
{
  std::ostream& os=this->ostream();
  os << this->line_colour.name() << " stroke\n";
}

void 
Output::epsstream::_instantiate() {
}

}
