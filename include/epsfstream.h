/****************************************************************************            epsstream.h
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
#include "rectangle.h"
#include "list_set.h"

namespace Ariadne {
  namespace Postscript {
    // TODO: Move these to numerical_types.h
    template<typename R> double convert_to_double(const R& x) {
      return double(x); }

    template<> double convert_to_double(const Rational& x) {
      return x.get_d(); }

    template<> double convert_to_double(const Dyadic& x) {
      return x.get_d(); }


    class epsfstream : public std::ofstream {
     private:
      static const uint xBBoxSide=300;
      static const uint yBBoxSide=300;
      static const double linewidth=0.001;
    public:
      epsfstream(const char* fn, const Ariadne::Geometry::Rectangle<double>& bbox)
       : std::ofstream(fn)
      {
        open(bbox);
      }

      ~epsfstream() {
        close();
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

        (*this) << "/linewidth " << linewidth << " def\n";
        (*this) << "/black {0.0 0.0 0.0 setrgbcolor} def\n"
                << "/gray {0.5 0.5 0.5 setrgbcolor} def\n"
                << "/red {1.0 0.0 0.0 setrgbcolor} def\n"
                << "/green {0.0 1.0 0.0 setrgbcolor} def\n"
                << "/blue {0.0 0.0 1.0 setrgbcolor} def\n";

        (*this) << "/xl " << bbox.lower(0) << " def\n"
                << "/xu " << bbox.upper(0) << " def\n"
                << "/yl " << bbox.lower(1) << " def\n"
                << "/yu " << bbox.upper(1) << " def\n";
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
    };

    template<typename R>
    epsfstream&
    operator<<(epsfstream& eps, const Ariadne::Geometry::Rectangle<R>& r)
    {
      assert(r.dimension()==2);

      double rlx=convert_to_double(r.lower(0));
      double rux=convert_to_double(r.upper(0));
      double rly=convert_to_double(r.lower(1));
      double ruy=convert_to_double(r.upper(1));

      eps << rlx << ' ' << rly << " moveto\n"
          << rux << ' ' << rly << " lineto\n"
          << rux << ' ' << ruy << " lineto\n"
          << rlx << ' ' << ruy << " lineto\n"
          << rlx << ' ' << rly << " lineto\n"
          << "fill\n";
      return eps;
    }

    template<typename R, template<typename> class BS>
    epsfstream&
    operator<<(epsfstream& eps, const Ariadne::Geometry::ListSet<R,BS>& ds)
    {
      typedef typename Ariadne::Geometry::ListSet<R,BS>::const_iterator const_iterator;
      for(const_iterator set_iter=ds.begin(); set_iter!=ds.end(); ++set_iter) {
        eps << (*set_iter);
      }
    }

  }
}
