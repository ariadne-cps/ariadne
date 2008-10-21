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
trace(cairo_t *cr, const Polytope& p);

struct GraphicsObject {
  enum ShapeKind { BOX, POLYTOPE, CURVE };
  GraphicsObject(Colour fc, Polytope sh) : fill_colour(fc), shape(sh) { }
  Colour fill_colour;
  Polytope shape;
};

static void 
draw(cairo_t *cr, const std::vector<GraphicsObject>& objects, 
     int canvas_width, int canvas_height);



struct Graphic::Impl 
{
  Box bounding_box;
  PlanarProjectionMap projection_map;
  Colour fill_colour;
  std::vector<GraphicsObject> objects;
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
    ARIADNE_ASSERT(bx.dimension()==2);
    this->_impl->objects.push_back(GraphicsObject(this->_impl->fill_colour,polytope(bx)));
}

void Graphic::plot(const Polytope& p) {
    ARIADNE_ASSERT(p.dimension()==2);
    this->_impl->objects.push_back(GraphicsObject(this->_impl->fill_colour,polytope(p)));
}

void Graphic::clear() {
    this->_impl->objects.clear();
}



void trace(cairo_t *cr, const Box& bx) 
{
    // std::cerr << "trace(cairo_t *cr, Box bx) bx="<<bx<<std::endl;
    cairo_move_to (cr, bx[0].lower(), bx[1].lower());
    cairo_line_to (cr, bx[0].upper(), bx[1].lower());
    cairo_line_to (cr, bx[0].upper(), bx[1].upper());
    cairo_line_to (cr, bx[0].lower(), bx[1].upper());
    cairo_line_to (cr, bx[0].lower(), bx[1].lower());
}

void trace(cairo_t *cr, const Polytope& p) 
{
    // std::cerr << "trace(cairo_t *cr, Polytope p) p="<<p<<std::endl;
    ARIADNE_ASSERT(p.dimension()==2);
    ARIADNE_ASSERT(p.size()>=3);
    cairo_move_to (cr, p[0][0], p[0][1]);
    for(uint i=1; i!=p.size(); ++i) {
      cairo_line_to (cr, p[i][0], p[i][1]); 
    }
    cairo_line_to (cr, p[0][0], p[0][1]); 
}

void draw(cairo_t *cr, const std::vector<GraphicsObject>& objects, int canvas_width, int canvas_height) 
{
    //std::cerr<<"draw(...)\n  polytopes=("<<polytopes.size()<<")"<<polytopes<<std::endl;

    // Compute extreme values
    if(objects.empty()) {
      cairo_destroy (cr);
      return; 
    }

    Box bbox=objects[0].shape.bounding_box();
    for(uint i=1; i!=objects.size(); ++i) {
        bbox=hull(bbox,objects[i].shape.bounding_box()); 
    }

    // The bounding box for the actual used area
    Box lbbox=bbox+Vector<Interval>(2,Interval(-0.075,0.075));

    // The bounding box for the entire figure
    Box gbbox=bbox+Vector<Interval>(2,Interval(-0.1,0.1));
  
    //std::cerr<<"  bbox="<<bbox<<"\n  lbbox="<<lbbox<<"\n  gbbox="<<gbbox<<std::endl;

    // clear background 
    cairo_set_source_rgb (cr, 1,1,1);
    cairo_paint (cr);

    // compute user to canvas coordinate transformation
    cairo_set_line_width (cr,0.001);

    // compute user to canvas coordinate transformation
    double sc0=canvas_width/gbbox[0].width();
    double sc1=-canvas_height/gbbox[1].width();
    double tr0=-gbbox[0].lower();
    double tr1=-gbbox[1].upper();
    //std::cerr << "  sc0="<<sc0<<", sc1="<<sc1<<", tr0="<<tr0<<", tr1="<<tr1<<std::endl;
    cairo_scale (cr, sc0,sc1);
    cairo_translate(cr, tr0, tr1);
    
    for(uint i=0; i!=objects.size(); ++i) {
        const Colour& fc=objects[i].fill_colour;
        const Polytope& p=objects[i].shape;
        cairo_set_source_rgb (cr, fc.red, fc.green, fc.blue);
        trace (cr,p);
        cairo_fill (cr);
    }
    
    cairo_set_source_rgb (cr, 0,0,0);
    for(uint i=0; i!=objects.size(); ++i) {
        const Polytope& p=objects[i].shape;
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

    std::vector<GraphicsObject>& objects=this->_impl->objects;

    surface = cairo_image_surface_create (CAIRO_FORMAT_ARGB32, canvas_width, canvas_height);
    cr = cairo_create (surface);
 
    draw(cr, objects, canvas_width, canvas_height);

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
    std::vector<GraphicsObject>& objects=impl->objects;

    gint canvas_width  = widget->allocation.width;
    gint canvas_height = widget->allocation.height;
    
    // Get Cairo drawing context
    cr = gdk_cairo_create (widget->window);

    // Draw Cairo objects
    draw(cr, objects, canvas_width, canvas_height);
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

 
