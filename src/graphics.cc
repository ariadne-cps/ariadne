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
#undef HAVE_GMPXX_H

#include "macros.h"
#include "stlio.h"
#include "numeric.h"
#include "vector.h"
#include "matrix.h"
#include "function.h"
#include "point.h"
#include "box.h"
#include "curve.h"
#include "polytope.h"
#include "graphics.h"

#ifdef HAVE_GTK_H
#include <gtk/gtk.h>
#endif

#ifdef HAVE_CAIRO_H
#include <cairo/cairo.h>
#endif

namespace Ariadne {

// static const int DEFAULT_WIDTH = 800;
// static const int DEFAULT_HEIGHT = 800;
static const int DEFAULT_WIDTH = 1600;
static const int DEFAULT_HEIGHT = 1600;




std::vector<Point> extremal(const std::vector<Point> & points) {
    Polytope polytope(points);
    return reduce2d(polytope).vertices();
}
  

struct GraphicsObject {
    GraphicsObject(double lw, Colour lc, Colour fc, const DrawableInterface& sh)
        : line_width(lw), line_colour(lc), fill_colour(fc), shape_ptr(sh.clone()) { }
    double line_width;
    Colour line_colour;
    Colour fill_colour;
    shared_ptr<const DrawableInterface> shape_ptr;
};



struct Figure::Data 
{
    Data() : bounding_box(0), projection(2,0,1), current_line_width(1.0), current_line_colour(0,0,0), current_fill_colour(0.75,0.75,0.75) { }
    Box bounding_box;
    PlanarProjectionMap projection;
    bool current_line_style;
    double current_line_width;
    Colour current_line_colour;
    bool current_fill_style;
    Colour current_fill_colour;
    std::vector<GraphicsObject> objects;
};


Figure::~Figure()
{
    delete this->_data;
}

 
Figure::Figure()
    : _data(new Data()) 
{ 
    this->_data->bounding_box=Box(0);
    this->_data->projection=PlanarProjectionMap(2,0,1);
}

void Figure::draw(const DrawableInterface& shape)
{ 
    this->_data->objects.push_back(
        GraphicsObject(this->_data->current_line_width,this->_data->current_line_colour,this->_data->current_fill_colour,shape));
}


void Figure::set_projection(uint as, uint ix, uint iy)
{
    this->_data->projection=PlanarProjectionMap(as,ix,iy);
}

void Figure::set_projection_map(const ProjectionFunction& pf) 
{
    this->_data->projection=PlanarProjectionMap(pf.argument_size(),pf.p(0),pf.p(1));
}

void Figure::set_projection_map(const PlanarProjectionMap& p) 
{
    this->_data->projection=p;
}

void Figure::set_bounding_box(const Box& bx) 
{
    this->_data->bounding_box=bx;
}

PlanarProjectionMap Figure::get_projection_map() const
{
    return this->_data->projection;
}

Box Figure::get_bounding_box() const
{
    return this->_data->bounding_box;
}

void Figure::set_line_style(bool ls)
{
}

void Figure::set_line_width(double lw) 
{
    this->_data->current_line_width=lw;
}

void Figure::set_line_colour(Colour lc)
{ 
    this->_data->current_line_colour=lc;
}

void Figure::set_line_colour(double r, double g, double b)
{ 
    this->set_line_colour(Colour(r,g,b));
}

void Figure::set_fill_style(bool fs) 
{
}

void Figure::set_fill_colour(Colour fc)
{ 
    this->_data->current_fill_colour=fc;
}

void Figure::set_fill_colour(double r, double g, double b)
{ 
    this->set_fill_colour(Colour(r,g,b));
}

bool Figure::get_line_style() const
{ 
    return this->_data->current_line_style;
}

double Figure::get_line_width() const
{ 
    return this->_data->current_line_width;
}

Colour Figure::get_line_colour() const
{ 
    return this->_data->current_line_colour;
}


bool Figure::get_fill_style() const
{ 
    return this->_data->current_fill_style;
}

Colour Figure::get_fill_colour() const
{ 
    return this->_data->current_fill_colour;
}



void Figure::clear() {
    this->_data->objects.clear();
}


#ifdef HAVE_CAIRO_H

class CairoCanvas
    : public CanvasInterface
{
    friend class Figure;
  private:
    cairo_t *cr; uint ix; uint iy;
    Colour lc,fc;
  public:
    CairoCanvas(cairo_t *c, uint i, uint j) : cr(c), ix(i), iy(j), lc(0,0,0), fc(1,1,1) { }
    uint x_coordinate() const { return ix; }
    uint y_coordinate() const { return iy; }
    void move_to(double x, double y) { cairo_move_to (this->cr, x, y); }
    void line_to(double x, double y) { cairo_line_to (this->cr, x, y); }
    void circle(double x, double y, double r) { cairo_arc (cr, x, y, r, 0, 2*M_PI); }
    void dot(double x, double y) { static const double RADIUS=0.01; cairo_arc (cr, x, y, RADIUS, 0, 2*M_PI); }
    void stroke() { cairo_set_source_rgb(cr, lc.red,lc.green,lc.blue); cairo_stroke (this->cr); }
    void fill() { cairo_set_source_rgb(cr, fc.red,fc.green,fc.blue); cairo_fill_preserve (this->cr); this->stroke(); }
    void clip() { cairo_clip (this->cr); }
    void set_line_width(double lw) { cairo_set_line_width (cr,lw); }
    void set_line_colour(double r, double g, double b) { lc=Colour(r,g,b); }
    void set_fill_colour(double r, double g, double b) { fc=Colour(r,g,b); }
    void set_bounding_box(double x0, double x1, double y0, double y1) { ARIADNE_NOT_IMPLEMENTED; }
    double get_line_width() const { return cairo_get_line_width (cr); }

};






std::string str(double x) {
    std::stringstream ss;
    ss << x;
    return ss.str();
}





void Figure::_paint_all(CanvasInterface& canvas)
{
    //std::cerr<<"Figure::_paint_all(canvas)\n";

    Box bounding_box=this->_data->bounding_box;
    const PlanarProjectionMap proj=this->_data->projection;
    const std::vector<GraphicsObject>& objects=this->_data->objects;

   // Test if there are no objects to be drawn
    if(objects.empty()) {
        return;
    }

    // Don't attempt to compute a bounding box, as this relies on
    // a drawable object having one. Instead, the bounding box must be
    // specified explicitly
    if(bounding_box.dimension()==0) {
        bounding_box=Box(proj.n,Interval(-1,1));
/*
        if(objects.empty()) {
            bounding_box=Box(proj.n,Interval(-1,1));
        } else {
            bounding_box=objects[0].shape_ptr->bounding_box();
            for(uint i=1; i!=objects.size(); ++i) {
                bounding_box=hull(bounding_box,objects[i].shape_ptr->bounding_box());
            }
        }
*/
    }

    //std::cerr << "bounding_box="<<bounding_box<<"\n";

    // Project the bounding box onto the canvas
    Box bbox(2);
    if(bounding_box.dimension()==2) {
        bbox=bounding_box;
    } else {
        bbox[0]=bounding_box[proj.i];
        bbox[1]=bounding_box[proj.j];
    }

    //std::cerr << "bbox="<<bbox<<"\n";

    // The bounding box for the actual used area
    Box lbbox=bbox;
    lbbox[0]+=Interval(-1,1)*(bbox[0].radius()/25);
    lbbox[1]+=Interval(-1,1)*(bbox[1].radius()/25);

    //std::cerr << "lbbox="<<bbox<<"\n";


    cairo_t *cr=static_cast<CairoCanvas&>(canvas).cr;

    const int canvas_height=DEFAULT_HEIGHT;
    const int canvas_width=DEFAULT_WIDTH;

    const int left_margin=160;
    const int bottom_margin=40;
    const int top_margin=10;
    const int right_margin=10;

    // clear background
    cairo_set_source_rgb (cr, 1,1,1);
    cairo_paint (cr);

    // Save unclipped state and canvas coordinates
    cairo_save (cr);

    // Set clipping region
    cairo_move_to (cr, left_margin, canvas_height-bottom_margin);
    cairo_line_to (cr, canvas_width-right_margin, canvas_height-bottom_margin);
    cairo_line_to (cr, canvas_width-right_margin, top_margin);
    cairo_line_to (cr, left_margin, top_margin);
    cairo_line_to (cr, left_margin, canvas_height-bottom_margin);
    cairo_clip (cr);
    cairo_new_path (cr);
    
    cairo_set_line_width (cr,0.002*min(bbox[0].width(),bbox[1].width()));  
    
    // compute user to canvas coordinate transformation
    double ctr0=left_margin;
    double ctr1=top_margin;
    double sc0=(canvas_width-left_margin-right_margin)/lbbox[0].width();
    double sc1=-(canvas_height-top_margin-bottom_margin)/lbbox[1].width();
    double utr0=-lbbox[0].lower();
    double utr1=-lbbox[1].upper();

    // Scale to user coordinates
    cairo_translate(cr, ctr0, ctr1);
    cairo_scale (cr, sc0,sc1);
    cairo_translate(cr, utr0, utr1);


    

    // Check projection and bounding box have same values;
    ARIADNE_ASSERT(bounding_box.dimension()==proj.argument_size());
    ARIADNE_ASSERT(bounding_box.dimension()>proj.i);
    ARIADNE_ASSERT(bounding_box.dimension()>proj.j);

    // Draw shapes
    for(uint i=0; i!=objects.size(); ++i) {
        const DrawableInterface* shape_ptr=objects[i].shape_ptr.operator->();
        ARIADNE_ASSERT_MSG(bounding_box.dimension()==shape_ptr->dimension(),
                           "Shape dimension="<<shape_ptr->dimension()<<", bounding_box="<<bounding_box);
        const Colour& fc=objects[i].fill_colour;
        const Colour& lc=objects[i].line_colour;
        canvas.set_fill_colour(fc.red, fc.green, fc.blue);
        canvas.set_line_colour(lc.red, lc.green, lc.blue);
        shape_ptr->draw(canvas);
    }
    

    // Restore canvas coordinates and unclipped state
    cairo_restore (cr);

    // Set text font
    cairo_select_font_face (cr, "roman",CAIRO_FONT_SLANT_NORMAL,CAIRO_FONT_WEIGHT_NORMAL);
    cairo_set_font_size (cr, 30);

    // Set text colour
    cairo_set_source_rgb (cr, 0., 0., 0.);

    // Get axis label text
    std::string text_xl=str(bbox[0].lower());
    std::string text_xu=str(bbox[0].upper());
    std::string text_yl=str(bbox[1].lower());
    std::string text_yu=str(bbox[1].upper());

    // Write axis labels
    cairo_text_extents_t te;
    cairo_text_extents (cr, text_xl.c_str(), &te);
    cairo_move_to(cr, left_margin-2, canvas_height-bottom_margin+4+te.height);
    cairo_show_text (cr, text_xl.c_str());
    cairo_text_extents (cr, text_xu.c_str(), &te);
    cairo_move_to(cr, canvas_width-te.width-4, canvas_height-bottom_margin+4+te.height);
    cairo_show_text (cr, text_xu.c_str());

    cairo_text_extents (cr, text_yl.c_str(), &te);
    cairo_move_to(cr, left_margin-te.width-6, canvas_width-bottom_margin+2);
    cairo_show_text (cr, text_yl.c_str());
    cairo_text_extents (cr, text_yu.c_str(), &te);
    cairo_move_to(cr, left_margin-te.width-6, top_margin+te.height+2);
    cairo_show_text (cr, text_yu.c_str());
    
    // Draw bounding box
    cairo_set_line_width (cr, 2.0);
    cairo_set_source_rgb (cr, 0.0, 0.0, 0.0);
    cairo_move_to (cr, left_margin, canvas_height-bottom_margin);
    cairo_line_to (cr, canvas_width-right_margin, canvas_height-bottom_margin);
    cairo_line_to (cr, canvas_width-right_margin, top_margin);
    cairo_line_to (cr, left_margin, top_margin);
    cairo_line_to (cr, left_margin, canvas_height-bottom_margin);
    cairo_stroke (cr);

    cairo_destroy (cr);
}


void 
Figure::write(const char* cfilename) 
{
    //std::cerr<<"Figure::write(filename="<<cfilename<<")\n";
    cairo_surface_t *surface;
    cairo_t *cr;

    const int canvas_width = DEFAULT_WIDTH;
    const int canvas_height = DEFAULT_HEIGHT;

    const Box& bounding_box=this->_data->bounding_box;
    const PlanarProjectionMap& projection=this->_data->projection;
    const std::vector<GraphicsObject>& objects=this->_data->objects;

    surface = cairo_image_surface_create (CAIRO_FORMAT_ARGB32, canvas_width, canvas_height);
    cr = cairo_create (surface);
    CairoCanvas canvas(cr,projection.i,projection.j);

    this->_paint_all(canvas);
    
    std::string filename(cfilename);
    if(filename.rfind(".") != std::string::npos) {
    } else {
        filename=filename+".png";
    }

    cairo_surface_write_to_png (surface, filename.c_str());
    cairo_surface_destroy (surface);
}


#ifdef HAVE_GTK_H

void
paint (GtkWidget      *widget,
       GdkEventExpose *eev,
       gpointer        gdata)
{
    cairo_t *cr;
  
    Figure* figure=static_cast<Figure*>(gdata);
    PlanarProjectionMap projection=figure->get_projection_map();

    gint canvas_width  = widget->allocation.width;
    gint canvas_height = widget->allocation.height;
    
    // Get Cairo drawing context
    cr = gdk_cairo_create (widget->window);

    // Draw Cairo objects
    CairoCanvas canvas(cr,projection.i,projection.j);
    figure->_paint_all(canvas);
}

void Figure::display() 
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
                      const_cast<Figure*>(this));  //  here we can pass a pointer to a custom data structure

    // pack canvas widget into window
    gtk_container_add (GTK_CONTAINER (window), canvas);

    // show window and all it's children (just the canvas widget)
    gtk_widget_show_all (window);

    // enter main loop
    gtk_main ();

}

#else // NO GTK_H

void Figure::display() 
{
    throw std::runtime_error("No facilities for displaying graphics are available.");
}

#endif // HAVE_GTK_H

#else // NO CAIRO_H

void 
Figure::write(const char* filename) 
{
    throw std::runtime_error("No facilities for drawing graphics are available.");
}

void Figure::display() 
{
    throw std::runtime_error("No facilities for displaying graphics are available.");
}

#endif // HAVE_CAIRO_H



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

#include "box.h"
#include "zonotope.h"

namespace Ariadne {







} // namespace Ariadne

 
