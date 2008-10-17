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
#include "box.h"
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


struct Graphic::Impl 
{
    std::vector<Box> boxes;
};

Colour::Colour()
    : name("transparant"), red(255), green(255), blue(255), transparant(true) { }
Colour::Colour(unsigned char rd, unsigned char gr, unsigned char bl, bool tr) 
    : name(), red(rd), green(gr), blue(bl), transparant(tr) { }
Colour::Colour(const char* nm, unsigned char rd, unsigned char gr, unsigned char bl, bool tr) 
    : name(nm), red(rd), green(gr), blue(bl), transparant(tr) { }
std::ostream& operator<<(std::ostream& os, const Colour& c) {
  return os << "Colour( name=" << c.name << ", r=" << c.red << ", g=" << c.green << ", b=" << c.blue << " )"; }


const Colour transparant=Colour();

const Colour white=Colour("white",255,255,255);
const Colour black=Colour("black",0,0,0);
const Colour red=Colour("red",255,0,0);
const Colour green=Colour("green",0,255,0);
const Colour blue=Colour("blue",0,0,255);
const Colour yellow=Colour("yellow",255,255,0);
const Colour cyan=Colour("cyan",0,255,255);
const Colour magenta=Colour("magenta",255,0,255);

Graphic::~Graphic()
{
    delete this->_impl;
}

 
Graphic::Graphic()
    : _impl(new Impl()) 
{ 
}

void Graphic::plot(const Box& bx) {
    ARIADNE_ASSERT(bx.size()==2);
    _impl->boxes.push_back(bx);
}

void Graphic::clear() {
    _impl->boxes.clear();
}



void trace(cairo_t *cr, const Box& bx) 
{
    cairo_move_to (cr, bx[0].lower(), bx[1].lower());
    cairo_line_to (cr, bx[0].upper(), bx[1].lower());
    cairo_line_to (cr, bx[0].upper(), bx[1].upper());
    cairo_line_to (cr, bx[0].lower(), bx[1].upper());
    cairo_line_to (cr, bx[0].lower(), bx[1].lower());
}

void draw(cairo_t *cr, const std::vector<Box>& boxes, int canvas_width, int canvas_height) 
{
    //std::cerr << "draw boxes: " << boxes << std::endl;

    // Compute extreme values
    assert(!boxes.empty());
    Box bbox=boxes[0];
    for(uint i=1; i!=boxes.size(); ++i) {
        bbox=hull(bbox,boxes[i]); 
    }

    // The bounding box for the actual used area
    Box lbbox=bbox+Vector<Interval>(2,Interval(-0.075,0.075));

    // The bounding box for the entire figure
    Box gbbox=bbox+Vector<Interval>(2,Interval(-0.1,0.1));
   
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
    
    //std::cerr<<"gbbox="<<gbbox<<std::endl;
    //std::cerr<<"bbox="<<bbox<<std::endl;
    cairo_set_source_rgb (cr, 0.5, 0.5, 1.0);
    for(uint i=0; i!=boxes.size(); ++i) {
        const Box& bx=boxes[i];
        //std::cerr<<"boxes[i]="<<bx<<"\n" << std::endl;
        trace (cr,bx);
        cairo_fill (cr);
    }
    
    cairo_set_source_rgb (cr, 0,0,0);
    for(uint i=0; i!=boxes.size(); ++i) {
        const Box& bx=boxes[i];
        trace (cr,bx);
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

    std::vector<Box>& boxes=this->_impl->boxes;

    surface = cairo_image_surface_create (CAIRO_FORMAT_ARGB32, canvas_width, canvas_height);
    cr = cairo_create (surface);
 
   draw(cr, boxes, canvas_width, canvas_height);

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
    std::vector<Box>& boxes=impl->boxes;

    gint canvas_width  = widget->allocation.width;
    gint canvas_height = widget->allocation.height;
    
    // Get Cairo drawing context
    cr = gdk_cairo_create (widget->window);

    // Draw Cairo objects
    draw(cr, boxes, canvas_width, canvas_height);
}

void Graphic::display() 
{
    const std::vector<Box>& boxes=this->_impl->boxes;
    std::cerr << "Graphic::boxes=" << boxes << std::endl;
    
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
                      const_cast<std::vector<Box>*>(&boxes));  //  here we can pass a pointer to a custom data structure 

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

} // namespace Ariadne

 
