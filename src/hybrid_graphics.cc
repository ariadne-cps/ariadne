/***************************************************************************
 *            hybrid_graphics.cc
 *
 *  Copyright 2011  Pieter Collins
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

#include "config.h"
#undef HAVE_GMPXX_H

#include "macros.h"
#include "stlio.h"
#include "numeric.h"
#include "space.h"
#include "point.h"
#include "box.h"
#include "geometry2d.h"
#include "discrete_location.h"
#include "hybrid_graphics.h"

#ifdef HAVE_GTK_H
#include <gtk/gtk.h>
#endif

#ifdef HAVE_CAIRO_H
#include <cairo/cairo.h>
#endif

namespace Ariadne {

static const int DEFAULT_WIDTH = 800;
static const int DEFAULT_HEIGHT = 800;

static const int LEFT_MARGIN = 160;
static const int BOTTOM_MARGIN = 40;
static const int TOP_MARGIN = 10;
static const int RIGHT_MARGIN = 10;


struct ImageSize2d {
    uint nx,ny;
    ImageSize2d(uint _nx,uint _ny) : nx(_nx), ny(_ny) { }
};

Projection2d projection(const RealSpace& spc, const Variables2d& axes) {
    return Projection2d(spc.dimension(),spc.index(axes.x_variable()),spc.index(axes.y_variable()));
}


void set_properties(CanvasInterface& canvas, const GraphicsProperties& properties);

void draw(CanvasInterface& canvas, const Set<DiscreteLocation>& locations, const Variables2d& axes, const HybridDrawableInterface& shape) {
    shape.draw(canvas,locations,axes);
}

void paint(CanvasInterface& canvas, const Set<DiscreteLocation>& locations, const Variables2d& axes, const List<HybridGraphicsObject>& objects) {
    for(uint i=0; i!=objects.size(); ++i) {
        const HybridDrawableInterface* shape_ptr=objects[i].shape_ptr.operator->();
        set_properties(canvas, objects[i].properties);
        shape_ptr->draw(canvas,locations,axes);
    }
}

HybridFigure::~HybridFigure() {
}


HybridFigure::HybridFigure()
    : axes(RealVariable("x"),RealVariable("y"))
{
}





class CairoCanvas
    : public CanvasInterface
{
    friend class Figure;
  private:
    cairo_t *cr;
    double lw; // The line width in pixels
    Colour lc,fc; // The line and fill colours
    double fo; // The fill opacity
  public:
    ~CairoCanvas();
    CairoCanvas(const ImageSize2d& size, const Box2d& bounds);
    CairoCanvas(cairo_t *c);
    void initialise(std::string x, std::string y, double xl, double xu, double yl, double yu);
    void finalise();
    void move_to(double x, double y) { cairo_move_to (cr, x, y); }
    void line_to(double x, double y) { cairo_line_to (cr, x, y); }
    void circle(double x, double y, double r) { cairo_arc (cr, x, y, r, 0, 2*M_PI); }
    void dot(double x, double y) { static const double RADIUS=0.01; cairo_arc (cr, x, y, RADIUS, 0, 2*M_PI); }
    void stroke();
    void fill() { cairo_set_source_rgba(cr,fc.red,fc.green,fc.blue,fo); cairo_fill_preserve (cr); this->stroke(); }
    void set_line_width(double lw) { this->lw=lw; }
    void set_line_colour(double r, double g, double b) { lc=Colour(r,g,b); }
    void set_fill_opacity(double o) { fo=o; }
    void set_fill_colour(double r, double g, double b) { fc=Colour(r,g,b); }

    Vector2d scaling() const;
    Box2d bounds() const;
  public:
    ImageSize2d size_in_pixels() const {
        return ImageSize2d(cairo_image_surface_get_width(cairo_get_target(cr))-(LEFT_MARGIN+RIGHT_MARGIN),
                           cairo_image_surface_get_height(cairo_get_target(cr))-(BOTTOM_MARGIN+TOP_MARGIN)); }
};


void
HybridFigure::write(const char* cfilename) const
{
    this->write(cfilename, DEFAULT_WIDTH, DEFAULT_HEIGHT);
}


void
HybridFigure::write(const char* cfilename, uint drawing_width, uint drawing_height) const
{
    cairo_surface_t *surface;
    cairo_t *cr;

    const int canvas_width = drawing_width+LEFT_MARGIN+RIGHT_MARGIN;
    const int canvas_height = drawing_height+BOTTOM_MARGIN+TOP_MARGIN;;

    surface = cairo_image_surface_create (CAIRO_FORMAT_ARGB32, canvas_width, canvas_height);
    cr = cairo_create (surface);
    CairoCanvas canvas(cr);

    this->_paint_all(canvas);

    std::string filename(cfilename);
    if(filename.rfind(".") != std::string::npos) {
    } else {
        filename=filename+".png";
    }

    cairo_surface_write_to_png (surface, filename.c_str());
    //cairo_surface_destroy (surface);
}

void HybridFigure::_paint_all(CanvasInterface& canvas) const
{
    // Project the bounding box onto the canvas
    double xl=numeric_cast<double>(bounds[axes.x_variable()].lower());
    double xu=numeric_cast<double>(bounds[axes.x_variable()].upper());
    double yl=numeric_cast<double>(bounds[axes.y_variable()].lower());
    double yu=numeric_cast<double>(bounds[axes.y_variable()].upper());

    canvas.initialise(axes.x_variable().name(),axes.y_variable().name(),xl,xu,yl,yu);

    // Draw shapes
    for(uint i=0; i!=objects.size(); ++i) {
        const HybridDrawableInterface& shape=objects[i].shape_ptr.operator*();
        set_properties(canvas, objects[i].properties);
        shape.draw(canvas,this->locations,this->axes);
    }

    canvas.finalise();
}



} // namespace Ariadne


