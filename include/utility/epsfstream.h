/****************************************************************************
 * epsstream.h
 *
 *  22 June 2005
 *  Copyright  2005  Alberto Casagrande, Pieter Collins
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

/*! \file epsfstream.h
 *  \brief Encapsulated postscript output.
 */
 
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string.h>

#include "../geometry/point.h"
#include "../geometry/rectangle.h"
#include "../geometry/list_set.h"
#include "../geometry/grid_set.h"
#include "../geometry/parallelotope.h"
#include "../geometry/zonotope.h"
#include "../geometry/polyhedron.h"
#include "../geometry/partition_tree_set.h"
#include "../linear_algebra/matrix.h"
#include "../evaluation/affine_map.h"

#define SCALE_DIMENSION 3.5 

namespace Ariadne {
  namespace Postscript {
    
    typedef struct{
      size_t pos;
      double radiant;
    } radiant_pointer_type;

    inline bool
    is_smaller_than(const radiant_pointer_type &a, 
                    const radiant_pointer_type &b) {
      return (a.radiant<b.radiant);
    }
   
    template <typename R>
    std::vector< Geometry::Point<R> >
    order_around_a_point(const std::vector< Geometry::Point<R> > &vertices, 
                         const Geometry::Point<R> &centre)
    {
      std::vector< LinearAlgebra::Vector<R> > vert_pos(vertices.size());
      
      for (size_t i=0; i< vertices.size(); i++) 
        vert_pos[i]=vertices[i].position_vector()-
                         centre.position_vector();
      
      radiant_pointer_type pointers[vertices.size()];
      
      double tangent=0;
      R tangent_R=0;
      
      for (size_t i=0; i< vert_pos.size(); i++) {
        pointers[i].pos=i;
        if (vert_pos[i](1)==0) {
          if (vert_pos[i](0) >0)
            pointers[i].radiant=acos(0.0);
          else
            pointers[i].radiant=-acos(0.0);
        } else {
          tangent_R=vert_pos[i](0)/vert_pos[i](1);
          tangent=convert_to<double>(tangent_R);
          pointers[i].radiant=atan(tangent);
	
          if (vert_pos[i](1) <0)
            pointers[i].radiant+=acos(-1.0);
        }
      }
      
      std::sort(pointers, pointers + vert_pos.size(), is_smaller_than);
      
      std::vector< Geometry::Point<R> > new_vector(vertices.size()); 
   
      for (size_t i=0; i< vertices.size(); i++) {
        new_vector[i]=vertices[pointers[i].pos];
      }

      return new_vector;
    }
    
    template<typename R>
    Ariadne::Geometry::Rectangle<double>
    convert_to_double_rectangle(const Ariadne::Geometry::Rectangle<R>& r) {
      Ariadne::Geometry::Rectangle<double> result(r.dimension());
      for(dimension_type i=0; i!=result.dimension(); ++i) {
        result.set_lower_bound(i,convert_to<double>(r.lower_bound(i)));
        result.set_upper_bound(i,convert_to<double>(r.upper_bound(i)));
      }
      return result;
    }
    
    template<typename R>
    class epsfstream : public std::ofstream {
      
      typedef Ariadne::Evaluation::AffineMap<R> ProjectionMap;

     private:
      static const uint xBBoxSide=300;
      static const uint yBBoxSide=300;
      static const double linewidth=0.0000001;
      
      ProjectionMap p_map;
      
     public:
      //FIXME: only needed for Boost Python interface.
      epsfstream(const epsfstream<R>& ofs) : std::ofstream("fnord") { };

      std::string line_colour;
      std::string fill_colour;
      bool line_style;
      bool fill_style;
      
      //template<typename R>
      epsfstream(const char* fn, const Ariadne::Geometry::Rectangle<R>& bbox)
       : std::ofstream(fn), line_colour("black"), fill_colour("green"), line_style(true), fill_style(true)
      {
	if (bbox.dimension()!=2)
	  throw std::runtime_error("epsfstream: the bounding box hasn't dimension 2."); 

	Ariadne::LinearAlgebra::identity_matrix<R> p_matrix(2);
        Ariadne::LinearAlgebra::Vector<R> p_vector(2);
	
	this->p_map=ProjectionMap(p_matrix,p_vector);
	
        this->open(bbox);
      }

      //template<typename R>
      epsfstream(const char* fn, const Ariadne::Geometry::Rectangle<R>& bbox, const unsigned int &x,  const unsigned int& y)
       : std::ofstream(fn), line_colour("black"), fill_colour("green"), line_style(true), fill_style(true)
      {
	if (bbox.dimension()<=x)
	  throw std::runtime_error("epsfstream: the given x dimension is greater than the space dimension."); 
	if (bbox.dimension()<=y)
	  throw std::runtime_error("epsfstream: the given y dimension is greater than the space dimension."); 
	
        Ariadne::LinearAlgebra::Matrix<R> p_matrix(2,bbox.dimension());
        Ariadne::LinearAlgebra::Vector<R> p_vector(2);
	
        p_matrix(0,x)=1.0;
        p_matrix(1,y)=1.0;

	this->p_map=ProjectionMap(p_matrix,p_vector);

        Ariadne::Geometry::Rectangle<R> proj_bbox=this->p_map(bbox);
	
	this->open(proj_bbox);
      }

      epsfstream(const char* fn, const Ariadne::Geometry::Rectangle<R>& bbox, const unsigned int &x,  const unsigned int& y, const char* x_name, const char* y_name)
       : std::ofstream(fn), line_colour("black"), fill_colour("green"), line_style(true), fill_style(true)
      {
	if (bbox.dimension()<=x)
	  throw std::runtime_error("epsfstream: the given x dimension is greater than the space dimension."); 
	if (bbox.dimension()<=y)
	  throw std::runtime_error("epsfstream: the given y dimension is greater than the space dimension."); 
	
        Ariadne::LinearAlgebra::Matrix<R> p_matrix(2,bbox.dimension());
        Ariadne::LinearAlgebra::Vector<R> p_vector(2);
	
        p_matrix(0,x)=1.0;
        p_matrix(1,y)=1.0;

	this->p_map=ProjectionMap(p_matrix,p_vector);

        Ariadne::Geometry::Rectangle<R> proj_bbox=this->p_map(bbox);
	
	Ariadne::Geometry::Point<R> l(proj_bbox.lower_corner());
	Ariadne::Geometry::Point<R> u(proj_bbox.upper_corner());

	l[0]-=((u[0]-l[0])/SCALE_DIMENSION);
	l[1]-=((u[1]-l[1])/SCALE_DIMENSION);
	//u[0]+=SCALE_DIMENSION;u[1]+=SCALE_DIMENSION;
	
	proj_bbox=Ariadne::Geometry::Rectangle<R>(l,u);
	this->open(proj_bbox);

	this->trace_scale(proj_bbox,6,6,x_name,y_name);
      }
      
      ~epsfstream() {
        this->close();
      }

      inline const ProjectionMap &projection() const { return this->p_map; }

      //template<typename R> 
      inline void open(const Ariadne::Geometry::Rectangle<R>& bbox) {
        open(convert_to_double_rectangle(bbox));
      }

      void
      trace_scale(const Ariadne::Geometry::Rectangle<R>& bb, const int &x_step, 
		      const int &y_step, const char* x_name, const char* y_name) 
      {
        R lx=bb.lower_bound(0);
        R ly=bb.lower_bound(1);
        R ux=bb.upper_bound(0);
        R uy=bb.upper_bound(1);

	R scale_x=((ux-lx)/(SCALE_DIMENSION + 1));
	R scale_y=((uy-ly)/(SCALE_DIMENSION + 1));

        lx=lx+0.8*scale_x;
        ly=ly+0.8*scale_y;

        double rlx=convert_to<double>(lx);
        double rly=convert_to<double>(ly);
        double rux=convert_to<double>(ux);
        double ruy=convert_to<double>(ux);

        double lscale_x=convert_to<double>(scale_x);
        double lscale_y=convert_to<double>(scale_y);

	double setlinewidth = ((rux-rlx)+(ruy-rly))/(SCALE_DIMENSION*300);
		
	double fontsize = ((rux-rlx)+(ruy-rly))/(SCALE_DIMENSION*20);
		
        *this << setlinewidth << " setlinewidth"<< std::endl
	      << "fillcolour black" << std::endl
	      << "/Times-Roman findfont" << std::endl
	      << fontsize  << " scalefont" << std::endl
	      << "setfont" << std::endl;
	
        *this << "newpath" << std::endl
	    << rlx << ' ' << rly << " moveto"<< std::endl
            << rux << ' ' << rly << " lineto"<< std::endl
            << "stroke"<< std::endl;
	    
	size_t prec=this->precision();

	if (x_step>0) {
	  for (double i=rlx; i< rux; i=i+((rux-rlx)/(x_step))) {
	  
	    *this << "newpath" << std::endl
	        << i << ' ' << rly - 0.05*lscale_y << " moveto"<< std::endl
                << i << ' ' << rly + 0.05*lscale_y << " lineto"<< std::endl
                << "stroke"<< std::endl;

	    *this << "newpath" << std::endl
  	        << i-0.05*lscale_x << ' ' << rly - 0.2*lscale_y 
		<< " moveto"<< std::endl;

	    this->precision(1);
            *this << "("<< i << ") show"<< std::endl;
	    this->precision(prec);
	  }

	  *this << "newpath" << std::endl
	        << rux << ' ' << rly - 0.05*lscale_y << " moveto"<< std::endl
                << rux << ' ' << rly + 0.05*lscale_y << " lineto"<< std::endl
                << "stroke"<< std::endl;

	  *this << "newpath" << std::endl
  	        << rux-0.05*lscale_x << ' ' << rly - 0.2*lscale_y 
		<< " moveto"<< std::endl;

	  this->precision(1);
          *this << "("<< rux << ") show"<< std::endl;
	  this->precision(prec);
	 
          double x_name_len=((double)strlen(x_name))/28;

          *this << "/Times-Roman findfont" << std::endl
	      << 1.5* fontsize  << " scalefont" << std::endl
	      << "setfont" << std::endl
	      << "newpath" << std::endl
	      << (rux+rlx-x_name_len)/2 << ' ' << rly - 0.5*lscale_y
	      << " moveto"<< std::endl
	      << "("<< x_name << ") show"<< std::endl
	      << "/Times-Roman findfont" << std::endl
	      << fontsize  << " scalefont" << std::endl
	      << "setfont" << std::endl;
	}
	
	*this << "newpath" << std::endl
	    << rlx << ' ' << rly << " moveto"<< std::endl
            << rlx << ' ' << ruy << " lineto"<< std::endl
            << "stroke"<< std::endl;

        if (y_step>0) {
	  for (double i=rly; i<= ruy; i=i+((ruy-rly)/(y_step))) {
	  
	    *this << "newpath" << std::endl
	        << rlx - 0.05*lscale_x << ' ' << i << " moveto"<< std::endl
                << rlx + 0.05*lscale_x << ' ' << i << " lineto"<< std::endl
                << "stroke"<< std::endl;

	    *this << "newpath" << std::endl
	        << rlx - 0.30*lscale_x << ' ' << i-0.05*lscale_y 
		<< " moveto"<< std::endl;

            this->precision(1);
            *this << "("<< i << ") show"<< std::endl;
	    this->precision(prec);
	  }
          
	  *this << "newpath" << std::endl
	        << rlx - 0.05*lscale_x << ' ' << ruy << " moveto"<< std::endl
                << rlx + 0.05*lscale_x << ' ' << ruy << " lineto"<< std::endl
                << "stroke"<< std::endl;

	  *this << "newpath" << std::endl
	        << rlx - 0.30*lscale_x << ' ' << ruy-0.05*lscale_y 
		<< " moveto"<< std::endl;

          this->precision(1);
          *this << "("<< ruy << ") show"<< std::endl;
	  this->precision(prec);

          double y_name_len=((double)strlen(y_name))/28;

          *this << "/Times-Roman findfont" << std::endl
	      << 1.5* fontsize  << " scalefont" << std::endl
	      << "setfont" << std::endl
	      << "newpath" << std::endl
	      << rlx - 0.5*lscale_x << ' ' << (ruy+rly-y_name_len)/2 
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
      
      void open(const Ariadne::Geometry::Rectangle<double>& bbox) {
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
       
        (*this) << "/xl " << bbox.lower_bound(0) << " def\n"
                << "/xu " << bbox.upper_bound(0) << " def\n"
                << "/yl " << bbox.lower_bound(1) << " def\n"
                << "/yu " << bbox.upper_bound(1) << " def\n";
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

      void close() {
        (*this) << "grestore showpage\n"
                   "%%Trailer\n";
        std::ofstream::close();
      }
      
      void set_line_style(bool ls) {
        this->line_style=ls;
      }
    
      void set_fill_style(bool fs) {
        this->fill_style=fs;
      }
    
      void set_pen_colour(const char* pc) {
        this->line_colour=pc;
      }
    
      void set_fill_colour(const char* fc) {
        this->fill_colour=fc;
      }

    };

        
    template<typename R>
    epsfstream<R>&
    trace(epsfstream<R>& eps, const Ariadne::Geometry::Rectangle<R>& r)
    {
      Ariadne::Geometry::Rectangle<R> proj_r=eps.projection()(r);

      double rlx=convert_to<double>(proj_r.lower_bound(0));
      double rux=convert_to<double>(proj_r.upper_bound(0));
      double rly=convert_to<double>(proj_r.lower_bound(1));
      double ruy=convert_to<double>(proj_r.upper_bound(1));

      eps << rlx << ' ' << rly << " moveto\n"
          << rux << ' ' << rly << " lineto\n"
          << rux << ' ' << ruy << " lineto\n"
          << rlx << ' ' << ruy << " lineto\n"
          << rlx << ' ' << rly << " lineto\n";
      return eps;
    }
   
    template<typename R>
    epsfstream<R>&
    trace(epsfstream<R>& eps, const std::vector< Geometry::Point<R> > &vertices) 
    {
      assert(vertices.size()>0);
     
      if (vertices[0].dimension()>2) {
        std::vector< Geometry::Point<R> > proj_vert(vertices.size());

	for (size_t j=0; j<vertices.size(); j++)
          proj_vert[j]=eps.projection()(vertices[j]);
	
	return trace(eps, proj_vert);
      }
      
      Geometry::Point<R> baricentre=vertices[0];

      for (size_t j=1; j<vertices.size(); j++) 
        for (size_t i=0; i<2; i++) 
          baricentre[i]=baricentre[i]+vertices[j][i];

      for (size_t i=0; i<2; i++) 
        baricentre[i]/=vertices.size();

      return trace(eps, vertices, baricentre);
    }
    
    template<typename R>
    epsfstream<R>&
    trace(epsfstream<R>& eps, const std::vector< Geometry::Point<R> > &vertices,
                   const Geometry::Point<R> &baricentre) {
      
      assert(vertices.size()>0);
     
      if (baricentre.dimension()>2) {
        std::vector< Geometry::Point<R> > proj_vert(vertices.size());
        Geometry::Point<R> proj_baric=eps.projection()(baricentre);
		
	for (size_t j=0; j<vertices.size(); j++) 
          proj_vert[j]=eps.projection()(vertices[j]);
	
	return trace(eps, proj_vert,proj_baric);
      }
	    	    
      std::vector< Geometry::Point<R> > ordered_vertices;
   
      ordered_vertices=order_around_a_point(vertices,baricentre);

      eps << ordered_vertices[0][0] << ' ' << ordered_vertices[0][1] 
          << " moveto\n";
      for (size_t i=1; i< ordered_vertices.size(); i++) {
        eps << ordered_vertices[i][0] << ' ' << ordered_vertices[i][1] 
            << " lineto\n";
      }
      return eps;
    }


    template<typename R>
    epsfstream<R>&
    operator<<(epsfstream<R>& eps, const Ariadne::Geometry::Point<R>& p) 
    {
      assert(p.dimension()>=2);

      if (p.dimension()>2) {
        Geometry::Point<R> proj_point=eps.projection()(p);
		
	return trace(eps, proj_point);
      }

      if(eps.fill_style) {
        eps << double(p[0]) << " " << double(p[0]) << "0.001 0 360 closepath\n";
        eps << eps.fill_colour << " fill\n";
      }
      if(eps.line_style) {
        eps << double(p[0]) << " " << double(p[0]) << "0.001 0 360 closepath\n";
        eps << eps.line_colour << " stroke\n";
      }
      return eps;
    }

    template<typename R>
    epsfstream<R>&
    operator<<(epsfstream<R>& eps, const Ariadne::Geometry::Rectangle<R>& r) 
    {
      assert(r.dimension()>=2);

      if(eps.fill_style) {
        trace(eps,r);
        eps << eps.fill_colour << " fill\n";
      }
      if(eps.line_style) {
        trace(eps,r);
        eps << eps.line_colour << " stroke\n";
      }
      return eps;
    }

    template<typename R>
    epsfstream<R>&
    operator<<(epsfstream<R>& eps, const Ariadne::Geometry::Parallelotope<R>& p)
    {
      
      if(eps.fill_style) {
        trace(eps,p.vertices(),p.centre());
        eps << eps.fill_colour << " fill\n";
      }
      if(eps.line_style) {
        trace(eps,p.vertices(),p.centre());
        eps << eps.line_colour << " stroke\n";
      }
      return eps;
    }
    
    template<typename R>
    epsfstream<R>&
    operator<<(epsfstream<R>& eps, const Ariadne::Geometry::Zonotope<R>& z)
    {
      if(eps.fill_style) {
        trace(eps,z.vertices(),z.centre());
        eps << eps.fill_colour << " fill\n";
      }
      if(eps.line_style) {
        trace(eps,z.vertices(),z.centre());
        eps << eps.line_colour << " stroke\n";
      }
      return eps;
    }
   
    template<typename R>
    epsfstream<R>&
    operator<<(epsfstream<R>& eps, const Ariadne::Geometry::Polyhedron<R>& p)
    {
      
      if(eps.fill_style) {
        trace(eps,p.vertices());
        eps << eps.fill_colour << " fill\n";
      }
      if(eps.line_style) {
        trace(eps,p.vertices());
        eps << eps.line_colour << " stroke\n";
      }
      return eps;
    }
    
    template<typename R, template<typename> class BS>
    epsfstream<R>&
    operator<<(epsfstream<R>& eps, const Ariadne::Geometry::ListSet<R,BS>& ds)
    {
      typedef typename Ariadne::Geometry::ListSet<R,BS>::const_iterator const_iterator;
        for(const_iterator set_iter=ds.begin(); set_iter!=ds.end(); ++set_iter) {
          //trace(eps,*set_iter) << eps.fill_colour << " fill\n";
          eps << *set_iter << eps.fill_colour << " fill\n";
        }
      if(eps.line_style) {
        for(const_iterator set_iter=ds.begin(); set_iter!=ds.end(); ++set_iter) {
          //trace(eps,*set_iter) << eps.line_colour << " stroke\n";
          eps<< *set_iter << eps.line_colour << " stroke\n";
        }
      }
      return eps;
    }

    template<typename R>
    epsfstream<R>&
    operator<<(epsfstream<R>& eps, const Ariadne::Geometry::GridMaskSet<R>& ds)
    {
      return eps << Ariadne::Geometry::ListSet<R,Ariadne::Geometry::Rectangle>(ds);
    }
    
    template<typename R>
    epsfstream<R>&
    operator<<(epsfstream<R>& eps, const Ariadne::Geometry::GridRectangleListSet<R>& ds)
    {
      return eps << Ariadne::Geometry::ListSet<R,Ariadne::Geometry::Rectangle>(ds);
    }
    
    template<typename R>
    epsfstream<R>&
    operator<<(epsfstream<R>& eps, const Ariadne::Geometry::GridCellListSet<R>& ds)
    {
      return eps << Ariadne::Geometry::ListSet<R,Ariadne::Geometry::Rectangle>(ds);
    }
    

    template<typename R>
    epsfstream<R>&
    operator<<(epsfstream<R>& eps, const Ariadne::Geometry::PartitionTreeSet<R>& ds)
    {
      return eps << Ariadne::Geometry::ListSet<R,Ariadne::Geometry::Rectangle>(ds);
    }

    template<typename R>
    epsfstream<R>&
    operator<<(epsfstream<R>& eps, const Ariadne::Geometry::PartitionTree<R>& pt)
    {
      for(typename Ariadne::Geometry::PartitionTree<R>::const_iterator iter = pt.begin(); iter!=pt.end(); ++iter) {
        eps << Ariadne::Geometry::Rectangle<R>(*iter);
      }
      return eps;
    }
    
  }
}

