/***************************************************************************
 *            cairo.hpp
 *
 *  Copyright 2011-17 Pieter Collins
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

#include "function/functional.hpp"
#include "config.h"

#include "output/graphics.hpp"

#ifdef HAVE_CAIRO_H
#include <cairo/cairo.h>
#endif

namespace Ariadne {

struct ImageSize2d {
    Nat nx,ny;
    ImageSize2d(Nat _nx,Nat _ny) : nx(_nx), ny(_ny) { }
};


#ifdef HAVE_CAIRO_H

class CairoCanvas
    : public CanvasInterface
{
    friend class Figure;
  private:
    cairo_t *cr;
    double lw; // The line width in pixels
    double dr; // The dot radius in pixels
    Colour lc,fc; // The line and fill colours
  public:
    ~CairoCanvas();
    CairoCanvas(const ImageSize2d& size);
    CairoCanvas(const ImageSize2d& size, const Box2d& bounds);
    CairoCanvas(cairo_t *c);
    Void initialise(StringType x, StringType y, double xl, double xu, double yl, double yu);
    Void write(const char* filename) const;
    Void finalise();
    Void move_to(double x, double y);
    Void line_to(double x, double y);
    Void circle(double x, double y, double r);
    Void dot(double x, double y);
    Void stroke();
    Void fill();
    Void set_dot_radius(double dr);
    Void set_line_width(double lw);
    Void set_line_colour(double r, double g, double b);
    Void set_fill_opacity(double o);
    Void set_fill_colour(double r, double g, double b);

    Vector2d scaling() const;
    Box2d bounds() const;
  public:
    ImageSize2d size_in_pixels() const;
};

#endif // HAVE_CAIRO_H

} // namespace Ariadne


