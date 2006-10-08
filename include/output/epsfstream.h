/****************************************************************************
 *            epsfstream.h
 *
 *  Copyright  2005-6  Alberto Casagrande, Pieter Collins
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

#include "../numeric/conversion.h"
#include "../geometry/point.h"
#include "../geometry/rectangle.h"
#include "../geometry/list_set.h"
#include "../geometry/grid_set.h"
#include "../geometry/parallelotope.h"
#include "../geometry/zonotope.h"
#include "../geometry/polytope.h"
#include "../geometry/partition_tree_set.h"
#include "../linear_algebra/matrix.h"
#include "../system/affine_map.h"

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
   
    template <typename R1, typename R2>
    Geometry::Point<R1> 
    approximate_point(const Geometry::Point<R2>& pt) 
    {
      Geometry::Point<R1> result(pt.dimension());
      for(size_type i=0; i!= result.dimension(); ++i) {
        result[i]=conv_approx<R1>(pt[i]);
      }
      return result;
    }
    
    template <typename R1, typename R2>
    Geometry::PointList<R1> 
    approximate_point_list(const Geometry::PointList<R2>& ptl) 
    {
      Geometry::PointList<R1> result(ptl.dimension(),ptl.size());
      for(size_type i=0; i!= ptl.size(); ++i) {
        result.push_back(approximate_point<R1>(ptl[i]));
      }
      return result;
    }
    
    class PlanarProjectionMap
    {
     public:
      PlanarProjectionMap() : _d(2), _i(0), _j(1) { }
      PlanarProjectionMap(dimension_type d, dimension_type i, dimension_type j)
        : _d(d), _i(i), _j(j) { assert(i<d && j<d); }
      template<typename R> Geometry::Point<double> operator() (const Geometry::Point<R>& pt) const {
        Geometry::Point<double> result(2); 
        assert(pt.dimension()==_d);
        result[0]=conv_approx<double>(pt[_i]); 
        result[1]=conv_approx<double>(pt[_j]); 
        return result;
      }
      template<typename R> Geometry::PointList<double> operator() (const Geometry::PointList<R>& ptl) const {
        Geometry::PointList<double> result(2,ptl.size());
        for(size_type i=0; i!=ptl.size(); ++i) { result[i]=this->operator()(ptl[i]); }
        return result;
      }
      template<typename R> Geometry::Rectangle<double> operator() (const Geometry::Rectangle<R>& r) const {
        Geometry::Rectangle<double> result(2); 
        assert(r.dimension()==_d);
        result[0]=Interval<double>(conv_approx<double>(r[_i].lower()),conv_approx<double>(r[_i].upper())); 
        result[1]=Interval<double>(conv_approx<double>(r[_j].lower()),conv_approx<double>(r[_j].upper())); 
        return result;
      }
      template<typename R> Geometry::PointList<double> operator() (const Geometry::Zonotope<R>& z) const {
        return this->operator()(z.vertices());
      }
      template<typename R> Geometry::PointList<double> operator() (const Geometry::Polytope<R>& p) const {
        return this->operator()(p.vertices());
      }
     private:
      dimension_type _d;
      dimension_type _i;
      dimension_type _j;
    };
    
    Geometry::Point<double> baricentre_of_points(const Geometry::PointList<double>& vertices)
    {
      Geometry::Point<double> baricentre(2);

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


    Geometry::PointList<double>
    order_around_a_point(const Geometry::PointList<double>& vertices, 
                         const Geometry::Point<double> &centre)
    {
      std::vector< LinearAlgebra::Vector<double> > vert_pos(vertices.size());
      
      for (size_t i=0; i< vertices.size(); i++) 
        vert_pos[i]=vertices[i].position_vector()-
                         centre.position_vector();
      
      radiant_pointer_type pointers[vertices.size()];
      
      double tangent=0;
      double tangent_R=0;
      
      for (size_t i=0; i< vert_pos.size(); i++) {
        pointers[i].pos=i;
        if (vert_pos[i](1)==0) {
          if (vert_pos[i](0) >0)
            pointers[i].radiant=acos(0.0);
          else
            pointers[i].radiant=-acos(0.0);
        } else {
          tangent_R=vert_pos[i](0)/vert_pos[i](1);
          tangent=conv_approx<double>(tangent_R);
          pointers[i].radiant=atan(tangent);
        
          if (vert_pos[i](1) <0)
            pointers[i].radiant+=acos(-1.0);
        }
      }
      
      std::sort(pointers, pointers + vert_pos.size(), is_smaller_than);
      
      Geometry::PointList<double> new_vector(2,vertices.size()); 
   
      for (size_t i=0; i< vertices.size(); i++) {
        Geometry::PointListReference<double> plref=new_vector[i];
        Geometry::Point<double> pt(vertices[pointers[i].pos]);
        plref=pt;
      }

      return new_vector;
    }
    
    template<typename R>
    Geometry::Rectangle<double>
    approximate_rectangle(Geometry::Rectangle<R>& r) 
    {
      Geometry::Rectangle<double> result(r.dimension());
      for(dimension_type i=0; i!=result.dimension(); ++i) {
        result.set_lower_bound(i,conv_approx<double>(r.lower_bound(i)));
        result.set_upper_bound(i,conv_approx<double>(r.upper_bound(i)));
      }
      return result;
    }
    
    class epsfstream : public std::ofstream {
     private:
      static const uint xBBoxSide=300;
      static const uint yBBoxSide=300;
      static const double linewidth=0.0000001;
      
      PlanarProjectionMap p_map;
     public:
      //FIXME: only needed for Boost Python interface.
      epsfstream(const epsfstream& ofs) : std::ofstream("fnord") { };

      std::string line_colour;
      std::string fill_colour;
      bool line_style;
      bool fill_style;
      
      epsfstream()
        : std::ofstream(), line_colour("black"), fill_colour("green"), line_style(true), fill_style(true)
      {
        this->p_map=PlanarProjectionMap(2,0,1);
      }
      
      template<typename R>
      epsfstream(const char* fn, const Ariadne::Geometry::Rectangle<R>& bbox)
        : std::ofstream(fn), line_colour("black"), fill_colour("green"), line_style(true), fill_style(true)
      {
        if (bbox.dimension()!=2) {
           throw std::runtime_error("epsfstream: the bounding box hasn't dimension 2."); 
        }
                
        this->p_map=PlanarProjectionMap(2,0,1);
        this->open(this->p_map(bbox));
      }

      template<typename R>
      epsfstream(const char* fn, const Ariadne::Geometry::Rectangle<R>& bbox, const unsigned int & ix,  const unsigned int& iy)
       : std::ofstream(fn), line_colour("black"), fill_colour("green"), line_style(true), fill_style(true)
      {
        if (bbox.dimension()<=ix)
          throw std::runtime_error("epsfstream: the given x dimension is greater than the space dimension."); 
        if (bbox.dimension()<=iy)
          throw std::runtime_error("epsfstream: the given y dimension is greater than the space dimension."); 
        
        this->p_map=PlanarProjectionMap(bbox.dimension(),ix,iy);
        this->open(this->p_map(bbox));
      }

      template<typename R>
      epsfstream(const char* fn, const Ariadne::Geometry::Rectangle<R>& bbox, 
                 const unsigned int &ix,  const unsigned int& iy, 
                 const char* x_name, const char* y_name)
        : std::ofstream(fn), line_colour("black"), fill_colour("green"), line_style(true), fill_style(true)
      {
        if (bbox.dimension()<=ix)
          throw std::runtime_error("epsfstream: the given x dimension is greater than the space dimension."); 
        if (bbox.dimension()<=iy)
          throw std::runtime_error("epsfstream: the given y dimension is greater than the space dimension."); 
        
        this->p_map=PlanarProjectionMap(bbox.dimension(),ix,iy);

        Geometry::Rectangle<double> proj_bbox=this->p_map(bbox);

        Geometry::Point<double> l(proj_bbox.lower_corner());
        Geometry::Point<double> u(proj_bbox.upper_corner());

        l[0]-=((u[0]-l[0])/SCALE_DIMENSION);
        l[1]-=((u[1]-l[1])/SCALE_DIMENSION);
        //u[0]+=SCALE_DIMENSION;u[1]+=SCALE_DIMENSION;
        
        proj_bbox=Geometry::Rectangle<double>(l,u);
        this->open(proj_bbox);

        this->trace_scale(proj_bbox,6,6,x_name,y_name);
      }
      
      ~epsfstream() {
        this->close();
      }

      inline const PlanarProjectionMap& projection_map() const { return this->p_map; }

      template<typename R> 
      inline void open(const char* fn, const Ariadne::Geometry::Rectangle<R>& bbox)
      {
        if (bbox.dimension()!=2) {
           throw std::runtime_error("epsfstream: the bounding box does not have dimension 2."); 
        }
        
        std::ofstream::open(fn);
        this->line_colour="black";
        this->fill_colour="green";
        this->line_style=true;
        this->fill_style=true;
        
        this->p_map=PlanarProjectionMap(2,0,1);
        
        this->open(p_map(bbox));
      }

      void
      trace_scale(const Ariadne::Geometry::Rectangle<double>& bb, const int &x_step, 
                  const int &y_step, const char* x_name, const char* y_name) 
      {
        double lx=bb.lower_bound(0);
        double ly=bb.lower_bound(1);
        double ux=bb.upper_bound(0);
        double uy=bb.upper_bound(1);

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
    
    
    
    epsfstream&
    trace_polygon(epsfstream& eps, const Geometry::PointList<double>& vertices)
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
    trace(epsfstream& eps, const Ariadne::Geometry::Point<double>& pt)
    {        
      eps << pt[0] << " " << pt[1] << "0.001 0 360 closepath\n";
      return eps;
    }
    
    
    epsfstream&
    trace(epsfstream& eps, const Ariadne::Geometry::Rectangle<double>& r)
    {
      Geometry::Rectangle<double> proj_r=eps.projection_map()(r);
      
      double lx=proj_r.lower_bound(0);
      double ux=proj_r.upper_bound(0);
      double ly=proj_r.lower_bound(1);
      double uy=proj_r.upper_bound(1);

      eps << lx << ' ' << ly << " moveto\n"
          << ux << ' ' << ly << " lineto\n"
          << ux << ' ' << uy << " lineto\n"
          << lx << ' ' << uy << " lineto\n"
          << lx << ' ' << ly << " lineto\n";
      return eps;
    }
    
    
    epsfstream&
    trace(epsfstream& eps, const Geometry::PointList<double>& vertices)
    {
      Geometry::Point<double> baricentre=baricentre_of_points(vertices);
      Geometry::PointList<double> ordered_vertices=order_around_a_point(vertices,baricentre);
      return trace_polygon(eps,ordered_vertices);
    }
      
    epsfstream&
    draw_polygon(epsfstream& eps, const Geometry::PointList<double>& vertices)
    {
      if(eps.fill_style) {
        trace_polygon(eps,vertices) << eps.fill_colour << " fill\n";
      }
      if(eps.line_style) {
        trace_polygon(eps,vertices) << eps.line_colour << " stroke\n";
      }
      return eps;
    }
    

    epsfstream&
    draw(epsfstream& eps, const Geometry::PointList<double>& vertices) 
    {
      Geometry::Point<double> baricentre=baricentre_of_points(vertices);
      Geometry::PointList<double> ordered_vertices=order_around_a_point(vertices,baricentre);
      return draw_polygon(eps,ordered_vertices);
    }
  
    epsfstream&
    draw(epsfstream& eps, const Geometry::PointList<double>& vertices,const Geometry::Point<double> baricentre) 
    {
      Geometry::PointList<double> ordered_vertices=order_around_a_point(vertices,baricentre);
      return draw_polygon(eps,ordered_vertices);
    }
     
     
    
    epsfstream&
    draw(epsfstream& eps, 
         const Geometry::Point<double>& pt)
    {
      if(eps.fill_style) {
        eps << pt[0] << " " << pt[1] << "0.001 0 360 closepath\n";
        eps << eps.fill_colour << " fill\n";
      }
      if(eps.line_style) {
        eps << pt[0] << " " << pt[1] << "0.001 0 360 closepath\n";
        eps << eps.line_colour << " stroke\n";
      }
      return eps;
    }

    epsfstream&
    draw(epsfstream& eps, const Geometry::Rectangle<double>& r) 
    {
      if(eps.fill_style) {
        trace(eps,r) << eps.fill_colour << " fill\n";
      }
      if(eps.line_style) {
        trace(eps,r) << eps.line_colour << " stroke\n";
      }
      return eps;
    }
    


    template<typename R>
    epsfstream&
    operator<<(epsfstream& eps, const Ariadne::Geometry::Point<R>& pt) 
    {
      return draw(eps, eps.projection_map()(pt));
    }

    template<typename R>
    epsfstream&
    operator<<(epsfstream& eps, const Ariadne::Geometry::Rectangle<R>& r) 
    {
      Geometry::Rectangle<double> dr=eps.projection_map()(r);
      return draw(eps,dr);
    }
    
    template<typename R>
    epsfstream&
    operator<<(epsfstream& eps, const Geometry::Zonotope<R>& z)
    {
      Geometry::Point<double> centre=eps.projection_map()(z.centre());
      Geometry::PointList<double> vertices=eps.projection_map()(z.vertices());      
      return draw(eps,vertices,centre);
    }
       
    template<typename R>
    epsfstream&
    operator<<(epsfstream& eps, const Geometry::Parallelotope<R>& p)
    {
      const Geometry::Zonotope<R>& z=p;
      return eps << z;
    }
       
    template<typename R>
    epsfstream&
    operator<<(epsfstream& eps, const Ariadne::Geometry::Polytope<R>& p)
    {
      return draw(eps,eps.projection_map()(p.vertices()));      
    }
    
    template<typename R, template<typename> class BS>
    epsfstream&
    operator<<(epsfstream& eps, const Ariadne::Geometry::ListSet<R,BS>& ds)
    {
      typedef typename Ariadne::Geometry::ListSet<R,BS>::const_iterator const_iterator;
      if(eps.fill_style) {
        for(const_iterator set_iter=ds.begin(); set_iter!=ds.end(); ++set_iter) {
          trace(eps,eps.projection_map()(*set_iter)) << eps.fill_colour << " fill\n";
        }
      }
      if(eps.line_style) {
        for(const_iterator set_iter=ds.begin(); set_iter!=ds.end(); ++set_iter) {
          trace(eps,eps.projection_map()(*set_iter)) << eps.line_colour << " stroke\n";
        }
      }
      return eps;
    }

    template<typename R>
    epsfstream&
    operator<<(epsfstream& eps, const Ariadne::Geometry::GridMaskSet<R>& ds)
    {
      return eps << Ariadne::Geometry::ListSet<R,Ariadne::Geometry::Rectangle>(ds);
    }
    
    template<typename R>
    epsfstream&
    operator<<(epsfstream& eps, const Ariadne::Geometry::GridBlockListSet<R>& ds)
    {
      return eps << Ariadne::Geometry::ListSet<R,Ariadne::Geometry::Rectangle>(ds);
    }
    
    template<typename R>
    epsfstream&
    operator<<(epsfstream& eps, const Ariadne::Geometry::GridCellListSet<R>& ds)
    {
      return eps << Ariadne::Geometry::ListSet<R,Ariadne::Geometry::Rectangle>(ds);
    }
    

    template<typename R>
    epsfstream&
    operator<<(epsfstream& eps, const Ariadne::Geometry::PartitionTreeSet<R>& ds)
    {
      return eps << Ariadne::Geometry::ListSet<R,Ariadne::Geometry::Rectangle>(ds);
    }

    template<typename R>
    epsfstream&
    operator<<(epsfstream& eps, const Ariadne::Geometry::PartitionTree<R>& pt)
    {
      for(typename Ariadne::Geometry::PartitionTree<R>::const_iterator iter = pt.begin(); iter!=pt.end(); ++iter) {
        eps << Ariadne::Geometry::Rectangle<R>(*iter);
      }
      return eps;
    }
    
  }
}
