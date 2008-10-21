/***************************************************************************
 *            graphics.cc
 *
 *  Copyright 2008  Pieter Collins
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

#include "macros.h"
#include "stlio.h"
#include "vector.h"
#include "matrix.h"
#include "point.h"
#include "box.h"
#include "geometry2d.h"
#include "graphics.h"

#ifdef HAVE_GTK_H
#include <gtk/gtk.h>
#endif

#ifdef HAVE_CAIRO_H
#include <cairo/cairo.h>
#endif

namespace Ariadne {

static const int DEFAULT_WIDTH = 800;
static const int DEFAULT_HEIGHT = 800;

static void 
trace(cairo_t *cr, const Box& bx);

static void 
draw(cairo_t *cr, const std::vector<Box>& boxes, 
     int canvas_width, int canvas_height);


struct GraphicsObject {
  Colour fill_colour;
  Polytope shape;
};

struct Graphic::Impl 
{
  Box bounding_box;
  PlanarProjectionMap projection_map;
  Colour fill_colour;
  std::vector<Polytope> polytopes;
};


Graphic::~Graphic()
{
    delete this->_impl;
}

 
Graphic::Graphic()
    : _impl(new Impl()) 
{ 
}


void Graphic::set_line_style(bool ls) 
{
}

void Graphic::set_line_width(double lw) 
{
}

void Graphic::set_line_colour(Colour lc)
{ 
}

void Graphic::set_fill_style(bool fs) 
{
}

void Graphic::set_fill_colour(Colour fc)
{ 
  this->_impl->fill_colour=fc;
}


void Graphic::plot(const Box& bx) {
    ARIADNE_ASSERT(bx.size()==2);
    _impl->polytopes.push_back(polytope(bx));
}

void Graphic::clear() {
    _impl->polytopes.clear();
}



void trace(cairo_t *cr, const Box& bx) 
{
    cairo_move_to (cr, bx[0].lower(), bx[1].lower());
    cairo_line_to (cr, bx[0].upper(), bx[1].lower());
    cairo_line_to (cr, bx[0].upper(), bx[1].upper());
    cairo_line_to (cr, bx[0].lower(), bx[1].upper());
    cairo_line_to (cr, bx[0].lower(), bx[1].lower());
}

void trace(cairo_t *cr, const Polytope& p) 
{
    ARIADNE_ASSERT(p.dimension()==2);
    ARIADNE_ASSERT(p.size()>=3);
    cairo_move_to (cr, p[0][0], p[0][1]);
    for(uint i=1; i!=p.size(); ++i) {
      cairo_line_to (cr, p[i][0], p[i][1]); 
    }
    cairo_line_to (cr, p[0][0], p[0][1]); 
}

void draw(cairo_t *cr, const std::vector<Polytope>& polytopes, const Colour& fill_colour, int canvas_width, int canvas_height) 
{
    //std::cerr << "draw(...)\n  polytopes=" << polytopes << std::endl;

    // Compute extreme values
    if(polytopes.empty()) {
      cairo_destroy (cr);
      return; 
    }

    Box bbox=polytopes[0].bounding_box();
    for(uint i=1; i!=polytopes.size(); ++i) {
        bbox=hull(bbox,polytopes[i].bounding_box()); 
    }

    // The bounding box for the actual used area
    Box lbbox=bbox+Vector<Interval>(2,Interval(-0.075,0.075));

    // The bounding box for the entire figure
    Box gbbox=bbox+Vector<Interval>(2,Interval(-0.1,0.1));
  
    // std::cerr << "  bbox="<<bbox<<std::endl;
    bbox[0]+=Interval(-0.1,0.1);
    bbox[1]+=Interval(-0.1,0.1);

    // clear background 
    cairo_set_source_rgb (cr, 1,1,1);
    cairo_paint (cr);

    // compute user to canvas coordinate transformation
    cairo_set_line_width (cr,0.001);

    // compute user to canvas coordinate transformation
    cairo_scale (cr, canvas_width/gbbox[0].width(), 
                    -canvas_height/gbbox[1].width());
    cairo_translate(cr, 0.0, 0.05-gbbox[1].width());
    
    cairo_set_source_rgb (cr, fill_colour.red, fill_colour.green, fill_colour.blue);
    for(uint i=0; i!=polytopes.size(); ++i) {
        const Polytope& p=polytopes[i];
        trace (cr,p);
        cairo_fill (cr);
    }
    
    cairo_set_source_rgb (cr, 0,0,0);
    for(uint i=0; i!=polytopes.size(); ++i) {
        const Polytope& p=polytopes[i];
        trace (cr,p);
        cairo_stroke (cr);
    }

    trace (cr,lbbox);
    cairo_stroke (cr);

    cairo_destroy (cr);
}


void 
Graphic::write(const char* filename) 
{
    cairo_surface_t *surface;
    cairo_t *cr;

    const int canvas_width = DEFAULT_WIDTH;
    const int canvas_height = DEFAULT_HEIGHT;

    std::vector<Polytope>& polytopes=this->_impl->polytopes;

    surface = cairo_image_surface_create (CAIRO_FORMAT_ARGB32, canvas_width, canvas_height);
    cr = cairo_create (surface);
 
    draw(cr, polytopes, this->_impl->fill_colour, canvas_width, canvas_height);

   cairo_destroy (cr);
   
   cairo_surface_write_to_png (surface, (std::string(filename)+".png").c_str());
   
   cairo_surface_destroy (surface);


}



#ifdef HAVE_GTK_H

void
paint (GtkWidget      *widget,
       GdkEventExpose *eev,
       gpointer        data)
{
    cairo_t *cr;
  
    Graphic::Impl* impl=static_cast<Graphic::Impl*>(data);
    //Graphic::Impl* impl=(Graphic::Impl*)data;
    std::vector<Polytope>& polytopes=impl->polytopes;

    gint canvas_width  = widget->allocation.width;
    gint canvas_height = widget->allocation.height;
    
    // Get Cairo drawing context
    cr = gdk_cairo_create (widget->window);

    // Draw Cairo objects
    draw(cr, polytopes, impl->fill_colour, canvas_width, canvas_height);
}

void Graphic::display() 
{
    
    GtkWidget *window;
    GtkWidget *canvas;

    int argc=0;
    char **argv;

    // initialize gtk
    gtk_init (&argc,&argv);

    // create a new top level window
    window   = gtk_window_new (GTK_WINDOW_TOPLEVEL);
    // make the gtk terminate the process the close button is pressed
    g_signal_connect (G_OBJECT (window), "delete-event",
                      G_CALLBACK (gtk_main_quit), NULL);

    // create a new drawing area widget
    canvas = gtk_drawing_area_new ();

    // set a requested (minimum size) for the canvas
    gtk_widget_set_size_request (canvas, DEFAULT_WIDTH, DEFAULT_HEIGHT);

    // connect our drawing method to the "expose" signal
    g_signal_connect (G_OBJECT (canvas), "expose-event",
                      G_CALLBACK (paint),
                      const_cast<Graphic::Impl*>(this->_impl));  //  here we can pass a pointer to a custom data structure 

    // pack canvas widget into window
    gtk_container_add (GTK_CONTAINER (window), canvas);

    // show window and all it's children (just the canvas widget)
    gtk_widget_show_all (window);

    // enter main loop
    gtk_main ();

}

#else

void Graphic::display() 
{
    throw std::runtime_error("No facilities for displaying graphics are available.");
}

#endif



Colour::Colour()
    : name("transparant"), red(1.0), green(1.0), blue(1.0), transparant(true) { }
Colour::Colour(double rd, double gr, double bl, bool tr) 
    : name(), red(rd), green(gr), blue(bl), transparant(tr) { }
Colour::Colour(const char* nm, double rd, double gr, double bl, bool tr) 
    : name(nm), red(rd), green(gr), blue(bl), transparant(tr) { }
std::ostream& operator<<(std::ostream& os, const Colour& c) {
  return os << "Colour( name=" << c.name << ", r=" << c.red << ", g=" << c.green << ", b=" << c.blue << " )"; }


const Colour transparant=Colour();

const Colour white=Colour("white",1.0,1.0,1.0);
const Colour black=Colour("black",0.0,0.0,0.0);
const Colour red=Colour("red",1.0,0.0,0.0);
const Colour green=Colour("green",0.0,1.0,0.0);
const Colour blue=Colour("blue",0.0,0.0,1.0);
const Colour yellow=Colour("yellow",1.0,1.0,0.0);
const Colour cyan=Colour("cyan",0.0,1.0,1.0);
const Colour magenta=Colour("magenta",1.0,0.0,1.0);




} // namespace Ariadne

 
