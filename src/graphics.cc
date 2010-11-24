/***************************************************************************
 *      graphics.cc
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
#include "function.h"
#include "point.h"
#include "box.h"
#include "graphics.h"

#ifdef HAVE_GTK_H
#include <gtk/gtk.h>
#endif

#ifdef HAVE_CAIRO_H
#include <cairo/cairo.h>
#endif

namespace Ariadne {

//static const int DEFAULT_WIDTH = 1200;
//static const int DEFAULT_HEIGHT = 800;

static const int LEFT_MARGIN = 160;
static const int BOTTOM_MARGIN = 40;
static const int TOP_MARGIN = 10;
static const int RIGHT_MARGIN = 10;


struct GraphicsProperties {
    GraphicsProperties()
     : line_style(true), line_width(1.0), line_colour(black), fill_style(true), fill_opacity(1.0), fill_colour(white) { }
    GraphicsProperties(bool ls, double lw, Colour lc, bool fs, double fo, Colour fc)
     : line_style(ls), line_width(lw), line_colour(lc), fill_style(fs), fill_opacity(fo), fill_colour(fc) { }
    bool line_style;
    double line_width;
    Colour line_colour;
    bool fill_style;
    double fill_opacity;
    Colour fill_colour;
};


struct GraphicsObject {
    GraphicsObject(const GraphicsProperties& gp, const DrawableInterface& sh)
     : properties(gp), shape_ptr(sh.clone()) { }
    GraphicsProperties properties;
    shared_ptr<const DrawableInterface> shape_ptr;
};



struct Figure::Data
{
    Data() : bounding_box(0), projection(2,0,1), properties() { }
    Box bounding_box;
    PlanarProjectionMap projection;
    GraphicsProperties properties;
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
    this->_data->objects.push_back(GraphicsObject(this->_data->properties,shape));
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
    this->_data->properties.line_style=ls;
}

void Figure::set_line_width(double lw)
{
    this->_data->properties.line_width=lw;
}

void Figure::set_line_colour(Colour lc)
{
    this->_data->properties.line_colour=lc;
}

void Figure::set_line_colour(double r, double g, double b)
{
    this->set_line_colour(Colour(r,g,b));
}

void Figure::set_fill_style(bool fs)
{
    this->_data->properties.fill_style=fs;
}

void Figure::set_fill_opacity(double fo)
{
    this->_data->properties.fill_opacity=fo;
}

void Figure::set_fill_colour(Colour fc)
{
    this->_data->properties.fill_colour=fc;
}

void Figure::set_fill_colour(double r, double g, double b)
{
    this->set_fill_colour(Colour(r,g,b));
}

bool Figure::get_line_style() const
{
    return this->_data->properties.line_style;
}

double Figure::get_line_width() const
{
    return this->_data->properties.line_width;
}

Colour Figure::get_line_colour() const
{
    return this->_data->properties.line_colour;
}


bool Figure::get_fill_style() const
{
    return this->_data->properties.fill_style;
}

double Figure::get_fill_opacity() const
{
    return this->_data->properties.fill_opacity;
}

Colour Figure::get_fill_colour() const
{
    return this->_data->properties.fill_colour;
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
    cairo_t *cr;
    uint ix; uint iy; // The indices of the variables plotted on the x and y axes
    Colour lc,fc; // The line and fill colours
    double fo; // The fill opacity
  public:
    CairoCanvas(cairo_t *c, uint i, uint j) : cr(c), ix(i), iy(j), lc(0,0,0), fc(1,1,1), fo(1.0) { }
    uint x_coordinate() const { return ix; }
    uint y_coordinate() const { return iy; }
    uint x_size_in_pixels() const { return cairo_image_surface_get_width(cairo_get_target(cr))-(LEFT_MARGIN+RIGHT_MARGIN); }
    uint y_size_in_pixels() const { return cairo_image_surface_get_height(cairo_get_target(cr))-(BOTTOM_MARGIN+TOP_MARGIN); }
    void move_to(double x, double y) { cairo_move_to (this->cr, x, y); }
    void line_to(double x, double y) { cairo_line_to (this->cr, x, y); }
    void circle(double x, double y, double r) { cairo_arc (cr, x, y, r, 0, 2*M_PI); }
    void dot(double x, double y) { static const double RADIUS=0.01; cairo_arc (cr, x, y, RADIUS, 0, 2*M_PI); }
    void stroke() { cairo_set_source_rgb(cr, lc.red,lc.green,lc.blue); cairo_stroke (this->cr); }
    void fill() { cairo_set_source_rgba(cr,fc.red,fc.green,fc.blue,fo); cairo_fill_preserve (this->cr); this->stroke(); }
    void clip() { cairo_clip (this->cr); }
    void set_line_width(double lw) { cairo_set_line_width (cr,0.01*lw); }
    void set_line_colour(double r, double g, double b) { lc=Colour(r,g,b); }
    void set_fill_opacity(double o) { fo=o; }
    void set_fill_colour(double r, double g, double b) { fc=Colour(r,g,b); }
    void get_bounding_box(double& xl, double& xu, double& yl, double& yu) const {
     xl=LEFT_MARGIN; yu=TOP_MARGIN;
     xu=x_size_in_pixels()+LEFT_MARGIN; yl=y_size_in_pixels()+TOP_MARGIN;
     std::cerr<<"Device: "<<xl<<" "<<xu<<" "<<yl<<" "<<yu<<"\n";
     cairo_device_to_user(cr,&xl,&yu);
     cairo_device_to_user(cr,&xu,&yl);
     std::cerr<<"User: "<<xl<<" "<<xu<<" "<<yl<<" "<<yu<<"\n";
    }
    double get_line_width() const { return cairo_get_line_width (cr); }

};






std::string str(Float x) {
    std::stringstream ss;
    ss << x;
    return ss.str();
}



std::ostream& operator<<(std::ostream& os, const DrawableInterface& drawable) {
    return drawable.write(os);
}


void Figure::_paint_all(CanvasInterface& canvas)
{
    //std::cerr<<"Figure::_paint_all(canvas)\n";

    Box bounding_box=this->_data->bounding_box;
    const PlanarProjectionMap proj=this->_data->projection;
    const std::vector<GraphicsObject>& objects=this->_data->objects;

    uint dimension=proj.argument_size();

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

    // Check projection and bounding box have same values.
    // If the bounding box dimension is 2, assume these are the bounds on the projected variables.
    // Project the bounding box onto the canvas
    Box2d bbox;
    //if(bounding_box.dimension()==2) {
    //    proj=PlanarProjectionMap(2,0,1);
    //}

    ARIADNE_ASSERT_MSG(bounding_box.dimension()==proj.argument_size(),"bounding_box="<<static_cast<const DrawableInterface&>(bounding_box)<<", projection="<<proj);
    ARIADNE_ASSERT(bounding_box.dimension()>proj.i);
    ARIADNE_ASSERT(bounding_box.dimension()>proj.j);
    bbox.xl=numeric_cast<double>(bounding_box[proj.i].lower());
    bbox.xu=numeric_cast<double>(bounding_box[proj.i].upper());
    bbox.yl=numeric_cast<double>(bounding_box[proj.j].lower());
    bbox.yu=numeric_cast<double>(bounding_box[proj.j].upper());

    cairo_t *cr=static_cast<CairoCanvas&>(canvas).cr;

    const int drawing_width = canvas.x_size_in_pixels();
    const int drawing_height = canvas.y_size_in_pixels();

    const int canvas_width = cairo_image_surface_get_width(cairo_get_target(cr));
    const int canvas_height = cairo_image_surface_get_height(cairo_get_target(cr));

    const int left_margin = LEFT_MARGIN;
    const int right_margin = RIGHT_MARGIN;
    const int bottom_margin = BOTTOM_MARGIN;
    const int top_margin = TOP_MARGIN;

    // clear background
    cairo_set_source_rgb (cr, 1,1,1);
    cairo_paint (cr);

    // Save unclipped state and canvas coordinates
    cairo_save (cr);

    cairo_set_line_width (cr,0.002*min(bbox.xu-bbox.xl,bbox.yu-bbox.yl));

    // Set clipping region
    cairo_move_to (cr, left_margin, top_margin+drawing_height);
    cairo_line_to (cr, left_margin+drawing_width, top_margin+drawing_height);
    cairo_line_to (cr, left_margin+drawing_width, top_margin);
    cairo_line_to (cr, left_margin, top_margin);
    cairo_line_to (cr, left_margin, top_margin+drawing_height);

    // Fill clipping region with a very light colour to indicate where figure
    // should be drawn
    cairo_set_source_rgb(cr, 0.95,1.00,0.95);
    cairo_fill_preserve (cr);
    cairo_clip (cr);
    cairo_new_path (cr);

    //std::cerr<<"cw="<<canvas_width<<" lm="<<left_margin<<" dw="<<drawing_width<<" rm="<<right_margin<<" xl="<<bbox.xl<<" xu="<<bbox.xu<<"\n";
    //std::cerr<<"ch="<<canvas_height<<"tm="<<top_margin<<" dw="<<drawing_height<<" bm="<<bottom_margin<<" yl="<<bbox.yl<<" yu="<<bbox.yu<<"\n";

    // compute device to user coordinate transformation
    double ctr0=left_margin;
    double ctr1=top_margin;
    double sc0=(drawing_width)/(bbox.xu-bbox.xl);
    double sc1=-(drawing_height)/(bbox.yu-bbox.yl);
    double utr0=(-bbox.xl);
    double utr1=(-bbox.yu);

    // Scale to user coordinates
    cairo_translate(cr, ctr0, ctr1);
    cairo_scale (cr, sc0,sc1);
    cairo_translate(cr, utr0, utr1);


    // Draw shapes
    for(uint i=0; i!=objects.size(); ++i) {
     const DrawableInterface* shape_ptr=objects[i].shape_ptr.operator->();
     if(shape_ptr->dimension()==0) { break; } // The dimension may be equal to two for certain empty sets.
     ARIADNE_ASSERT_MSG(dimension==shape_ptr->dimension(),
      "Shape "<<*shape_ptr<<", dimension="<<shape_ptr->dimension()<<", bounding_box="<<static_cast<const DrawableInterface&>(bounding_box));
     //const double& lw=objects[i].line_width;
     //canvas.set_line_width(lw);
     //canvas.set_line_width(objects[i].properties.line_width);
     canvas.set_fill_opacity(objects[i].properties.fill_opacity);
     const Colour& lc=objects[i].properties.line_colour;
     const Colour& fc=objects[i].properties.fill_colour;
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
    std::string text_xl=str(bbox.xl);
    std::string text_xu=str(bbox.xu);
    std::string text_yl=str(bbox.yl);
    std::string text_yu=str(bbox.yu);

    // Write axis labels
    cairo_text_extents_t te;
    cairo_text_extents (cr, text_xl.c_str(), &te);
    cairo_move_to(cr, left_margin-2, top_margin+drawing_height+4+te.height);
    cairo_show_text (cr, text_xl.c_str());
    cairo_text_extents (cr, text_xu.c_str(), &te);
    cairo_move_to(cr, left_margin+drawing_width-te.width-4, top_margin+drawing_height+4+te.height);
    cairo_show_text (cr, text_xu.c_str());

    cairo_text_extents (cr, text_yl.c_str(), &te);
    cairo_move_to(cr, left_margin-te.width-6, top_margin+drawing_height+2);
    cairo_show_text (cr, text_yl.c_str());
    cairo_text_extents (cr, text_yu.c_str(), &te);
    cairo_move_to(cr, left_margin-te.width-6, top_margin+te.height+2);
    cairo_show_text (cr, text_yu.c_str());

    // Draw bounding box
    cairo_set_line_width (cr, 2.0);
    cairo_set_source_rgb (cr, 0.0, 0.0, 0.0);
    cairo_move_to (cr, left_margin, top_margin+drawing_height);
    cairo_line_to (cr, left_margin+drawing_width, top_margin+drawing_height);
    cairo_line_to (cr, left_margin+drawing_width, top_margin);
    cairo_line_to (cr, left_margin, top_margin);
    cairo_line_to (cr, left_margin, top_margin+drawing_height);
    cairo_stroke (cr);

    cairo_destroy (cr);
}


void
Figure::write(const char* cfilename, uint drawing_width, uint drawing_height)
{
    //std::cerr<<"Figure::write(filename="<<cfilename<<")\n";
    cairo_surface_t *surface;
    cairo_t *cr;

    const int canvas_width = drawing_width+LEFT_MARGIN+RIGHT_MARGIN;
    const int canvas_height = drawing_height+BOTTOM_MARGIN+TOP_MARGIN;;

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
    gpointer     gdata)
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


