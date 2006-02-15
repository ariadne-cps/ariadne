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

#include <iostream>
#include <fstream>

#include "geometry/rectangle.h"
#include "geometry/list_set.h"
#include "geometry/grid_set.h"
#include "geometry/parallelopiped.h"
#include "geometry/partition_tree_set.h"

namespace Ariadne {
  namespace Postscript {

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
    
    class epsfstream : public std::ofstream {
     private:
      static const uint xBBoxSide=300;
      static const uint yBBoxSide=300;
      static const double linewidth=0.0000001;

     public:
      //FIXME: only needed for Boost Python interface.
      epsfstream(const epsfstream& ofs) : std::ofstream("fnord") { };

      std::string line_colour;
      std::string fill_colour;
      bool line_style;
      
      template<typename R>
      epsfstream(const char* fn, const Ariadne::Geometry::Rectangle<R>& bbox)
       : std::ofstream(fn), line_colour("black"), fill_colour("green"), line_style(true)
      {
        this->open(bbox);
      }

      ~epsfstream() {
        this->close();
      }

      template<typename R> 
      inline void open(const Ariadne::Geometry::Rectangle<R>& bbox) {
        open(convert_to_double_rectangle(bbox));
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
    
      void set_pen_colour(const char* pc) {
        this->line_colour=pc;
      }
    
      void set_fill_colour(const char* fc) {
        this->fill_colour=fc;
      }

    };

    
    
    template<typename R>
    epsfstream&
    trace(epsfstream& eps, const Ariadne::Geometry::Rectangle<R>& r)
    {
      double rlx=convert_to<double>(r.lower_bound(0));
      double rux=convert_to<double>(r.upper_bound(0));
      double rly=convert_to<double>(r.lower_bound(1));
      double ruy=convert_to<double>(r.upper_bound(1));

      eps << rlx << ' ' << rly << " moveto\n"
          << rux << ' ' << rly << " lineto\n"
          << rux << ' ' << ruy << " lineto\n"
          << rlx << ' ' << ruy << " lineto\n"
          << rlx << ' ' << rly << " lineto\n";
      return eps;
    }
    
    template<typename R>
    epsfstream&
    trace(epsfstream& eps, const Ariadne::Geometry::Parallelopiped<R>& p)
    {
      Ariadne::Geometry::Point<R> bl=p.centre()-p.generator(0)-p.generator(1);
      Ariadne::Geometry::Point<R> br=p.centre()+p.generator(0)-p.generator(1);
      Ariadne::Geometry::Point<R> tr=p.centre()+p.generator(0)+p.generator(1);
      Ariadne::Geometry::Point<R> tl=p.centre()-p.generator(0)+p.generator(1);

      eps << bl[0] << ' ' << bl[1] << " moveto\n"
          << br[0] << ' ' << br[1] << " lineto\n"
          << tr[0] << ' ' << tr[1] << " lineto\n"
          << tl[0] << ' ' << tl[1] << " lineto\n"
          << bl[0] << ' ' << bl[1] << " lineto\n";
      return eps;
    }

    template<typename R>
    epsfstream&
    operator<<(epsfstream& eps, const Ariadne::Geometry::Rectangle<R>& r) 
    {
      assert(r.dimension()==2);

      trace(eps,r);
      eps << eps.fill_colour << " fill\n";
      if(eps.line_style) {
        trace(eps,r);
        eps << eps.line_colour << " stroke\n";
      }
      return eps;
    }

    template<typename R>
    epsfstream&
    operator<<(epsfstream& eps, const Ariadne::Geometry::Parallelopiped<R>& p)
    {
      assert(p.dimension()==2);

      trace(eps,p);
      eps << eps.fill_colour << " fill\n";
      if(eps.line_style) {
        trace(eps,p);
        eps << eps.line_colour << " stroke\n";
      }
      return eps;
    }
    
    template<typename R, template<typename> class BS>
    epsfstream&
    operator<<(epsfstream& eps, const Ariadne::Geometry::ListSet<R,BS>& ds)
    {
      typedef typename Ariadne::Geometry::ListSet<R,BS>::const_iterator const_iterator;
        for(const_iterator set_iter=ds.begin(); set_iter!=ds.end(); ++set_iter) {
          trace(eps,*set_iter) << eps.fill_colour << " fill\n";
        }
      if(eps.line_style) {
        for(const_iterator set_iter=ds.begin(); set_iter!=ds.end(); ++set_iter) {
          trace(eps,*set_iter) << eps.line_colour << " stroke\n";
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
    operator<<(epsfstream& eps, const Ariadne::Geometry::GridRectangleListSet<R>& ds)
    {
      return eps << Ariadne::Geometry::ListSet<R,Ariadne::Geometry::Rectangle>(ds);
    }
    
    template<typename R>
    epsfstream&
    operator<<(epsfstream& eps, const Ariadne::Geometry::GridCellListSet<R>& ds)
    {
      return eps << Ariadne::Geometry::ListSet<R,Ariadne::Geometry::Rectangle>(ds);
    }
    
/*
    template<typename R>
    epsfstream&
    operator<<(epsfstream& eps, const Ariadne::Geometry::PartitionTreeSet<R>& ds)
    {
      return eps << Ariadne::Geometry::ListSet<R,Ariadne::Geometry::Rectangle>(ds);
    }
*/
    
  }
}
