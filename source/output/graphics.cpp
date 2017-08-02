/***************************************************************************
 *            graphics.cpp
 *
 *  Copyright 2008--17  Pieter Collins
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

#include "utility/standard.hpp"
#include "config.h"

#include "utility/macros.hpp"
#include "utility/stlio.hpp"
#include "numeric/numeric.hpp"
#include "function/function.hpp"
#include "geometry/point.hpp"
#include "geometry/box.hpp"
#include "output/geometry2d.hpp"
#include "output/graphics.hpp"
#include "output/cairo.hpp"

#ifdef HAVE_GTK_H
#include <gtk/gtk.h>
#endif

namespace Ariadne {

static const Int DEFAULT_WIDTH = 800;
static const Int DEFAULT_HEIGHT = 800;

static const Int LEFT_MARGIN = 160;
static const Int BOTTOM_MARGIN = 40;
static const Int TOP_MARGIN = 10;
static const Int RIGHT_MARGIN = 10;


StringType str(FloatDP x) {
    StringStream ss;
    ss << x;
    return ss.str();
}

Void draw(Figure& fig, const DrawableInterface& shape) {
    fig.draw(shape);
}

Void draw(Figure& fig, FloatDPApproximateBox const& box) {
    fig.draw(box);
}

OutputStream& operator<<(OutputStream& os, const DrawableInterface& drawable) {
    return drawable.write(os);
}




struct GraphicsObject {
    GraphicsObject(const GraphicsProperties& gp, const DrawableInterface& sh)
        : properties(gp), shape_ptr(sh.clone()) { }
    GraphicsProperties properties;
    std::shared_ptr<const DrawableInterface> shape_ptr;
};



struct Figure::Data
{
    Data() : bounding_box(0), projection(2,0,1), properties() { }
    ApproximateBoxType bounding_box;
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
    this->_data->bounding_box=ApproximateBoxType(0);
    this->_data->projection=PlanarProjectionMap(2,0,1);
}

Void Figure::draw(const DrawableInterface& shape)
{
    this->_data->objects.push_back(GraphicsObject(this->_data->properties,shape));
}


Void Figure::set_projection(Nat as, Nat ix, Nat iy)
{
    this->_data->projection=PlanarProjectionMap(as,ix,iy);
}

Void Figure::set_projection_map(const PlanarProjectionMap& p)
{
    this->_data->projection=p;
}

Void Figure::set_bounding_box(const ApproximateBoxType& bx)
{
    this->_data->bounding_box=bx;
}

PlanarProjectionMap Figure::get_projection_map() const
{
    return this->_data->projection;
}

ApproximateBoxType Figure::get_bounding_box() const
{
    return this->_data->bounding_box;
}

Void Figure::set_dot_radius(double dr)
{
    this->_data->properties.dot_radius=dr;
}

Void Figure::set_line_style(Bool ls)
{
    this->_data->properties.line_style=ls;
}

Void Figure::set_line_width(double lw)
{
    this->_data->properties.line_width=lw;
}

Void Figure::set_line_colour(Colour lc)
{
    this->_data->properties.line_colour=lc;
}

Void Figure::set_line_colour(double r, double g, double b)
{
    this->set_line_colour(Colour(r,g,b));
}

Void Figure::set_fill_style(Bool fs)
{
    this->_data->properties.fill_style=fs;
}

Void Figure::set_fill_opacity(double fo)
{
    this->_data->properties.fill_colour.opacity=fo;
}

Void Figure::set_fill_colour(Colour fc)
{
    this->_data->properties.fill_colour=fc;
}

Void Figure::set_fill_colour(double r, double g, double b)
{
    this->set_fill_colour(Colour(r,g,b,this->_data->properties.fill_colour.opacity));
}

Bool Figure::get_line_style() const
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


Bool Figure::get_fill_style() const
{
    return this->_data->properties.fill_style;
}

double Figure::get_fill_opacity() const
{
    return this->_data->properties.fill_colour.opacity;
}

Colour Figure::get_fill_colour() const
{
    return this->_data->properties.fill_colour;
}

Void Figure::draw(ApproximateBoxType const& box)
{
    ApproximateBoxSet box_set(box);
    DrawableInterface const& shape=box_set;
    this->draw(shape);
}



Void Figure::clear() {
    this->_data->objects.clear();
}


#ifdef HAVE_CAIRO_H

CairoCanvas::~CairoCanvas()
{
    cairo_surface_destroy(cairo_get_target(cr));
    cairo_destroy(cr);
}

CairoCanvas::CairoCanvas(cairo_t *c)
    : cr(c), lw(1.0), dr(1.0), lc(0,0,0), fc(1,1,1, 1.0)
{
}

CairoCanvas::CairoCanvas(const ImageSize2d& size)
    : cr(0), lw(1.0), lc(0.0,0.0,0.0), fc(1.0,1.0,1.0, 1.0)
{
    const Int canvas_width = size.nx+LEFT_MARGIN+RIGHT_MARGIN;
    const Int canvas_height = size.ny+BOTTOM_MARGIN+TOP_MARGIN;;

    cairo_surface_t* surface = cairo_image_surface_create (CAIRO_FORMAT_ARGB32, canvas_width, canvas_height);
    cr = cairo_create (surface);
}

CairoCanvas::CairoCanvas(const ImageSize2d& size, const Box2d& bounds)
    : CairoCanvas(size)
{
}

Vector2d CairoCanvas::scaling() const
{
    ImageSize2d sz = this->size_in_pixels();
    Box2d bb=this->bounds();
    return Vector2d((bb.xu-bb.xl)/sz.nx,(bb.yu-bb.yl)/sz.ny);
}

Box2d CairoCanvas::bounds() const
{
    double xl=LEFT_MARGIN;
    double yu=TOP_MARGIN;
    double xu=cairo_image_surface_get_width(cairo_get_target(cr))-RIGHT_MARGIN;
    double yl=cairo_image_surface_get_height(cairo_get_target(cr))-BOTTOM_MARGIN;
    cairo_device_to_user(cr,&xl,&yu);
    cairo_device_to_user(cr,&xu,&yl);
    return Box2d(xl,xu,yl,yu);
}

Void CairoCanvas::stroke()
{
    cairo_save(cr);

    // Set user and device space identical so that the line width is interpreted as pixels
    cairo_identity_matrix(cr);

    cairo_set_source_rgb(cr, lc.red,lc.green,lc.blue);
    cairo_set_line_width(cr, lw);
    cairo_stroke (this->cr);

    cairo_restore(cr);
}

Void CairoCanvas::fill() { 
    cairo_set_source_rgba(cr,fc.red,fc.green,fc.blue,fc.opacity); 
    cairo_fill_preserve (cr); 
    this->stroke(); 
}

ImageSize2d CairoCanvas::size_in_pixels() const {
    return ImageSize2d(cairo_image_surface_get_width(cairo_get_target(cr))-(LEFT_MARGIN+RIGHT_MARGIN),
                       cairo_image_surface_get_height(cairo_get_target(cr))-(BOTTOM_MARGIN+TOP_MARGIN)); 
}

Void CairoCanvas::move_to(double x, double y) { cairo_move_to (cr, x, y); }
Void CairoCanvas::line_to(double x, double y) { cairo_line_to (cr, x, y); }
Void CairoCanvas::circle(double x, double y, double r) { cairo_arc (cr, x, y, r, 0, 2*M_PI); }
Void CairoCanvas::dot(double x, double y) { cairo_arc (cr, x, y, dr/1000, 0, 2*M_PI); }
Void CairoCanvas::set_dot_radius(double dr) { this->dr=dr; }
Void CairoCanvas::set_line_width(double lw) { this->lw=lw; }
Void CairoCanvas::set_line_colour(double r, double g, double b) { lc.red=r; lc.green=g; lc.blue=b; }
Void CairoCanvas::set_fill_opacity(double o) { fc.opacity=o; }
Void CairoCanvas::set_fill_colour(double r, double g, double b) { fc.red=r; fc.green=g; fc.blue=b; }






// TODO: Use generic canvas routines; move cairo-specific functionality
// into CairoCanvas class.
Void CairoCanvas::initialise(StringType text_x, StringType text_y, double xl, double xu, double yl, double yu)
{


    CairoCanvas& cairo_canvas=*this;
    cairo_t *cr=cairo_canvas.cr;

    const ImageSize2d drawing_size = cairo_canvas.size_in_pixels();
    const Int drawing_width = drawing_size.nx;
    const Int drawing_height = drawing_size.ny;

    //const Int canvas_width = cairo_image_surface_get_width(cairo_get_target(cr));
    //const Int canvas_height = cairo_image_surface_get_height(cairo_get_target(cr));

    const Int left_margin = LEFT_MARGIN;
    //const Int right_margin = RIGHT_MARGIN;
    //const Int bottom_margin = BOTTOM_MARGIN;
    const Int top_margin = TOP_MARGIN;

    // clear background
    cairo_set_source_rgb (cr, 1,1,1);
    cairo_paint (cr);

    // Set text font
    cairo_select_font_face (cr, "roman",CAIRO_FONT_SLANT_NORMAL,CAIRO_FONT_WEIGHT_NORMAL);
    cairo_set_font_size (cr, 30);

    // Set text colour
    cairo_set_source_rgb (cr, 0., 0., 0.);

    // Get axis label text
    StringType text_xl=str(xl);
    StringType text_xu=str(xu);
    StringType text_yl=str(yl);
    StringType text_yu=str(yu);

    // Write axis labels
    cairo_text_extents_t te;
    cairo_text_extents (cr, text_xl.c_str(), &te);
    cairo_move_to(cr, left_margin-2, top_margin+drawing_height+4+te.height);
    cairo_show_text (cr, text_xl.c_str());
    cairo_text_extents (cr, text_xu.c_str(), &te);
    cairo_move_to(cr, left_margin+drawing_width-te.width-4, top_margin+drawing_height+4+te.height);
    cairo_show_text (cr, text_xu.c_str());
    cairo_text_extents (cr, text_x.c_str(), &te);
    cairo_move_to(cr, left_margin+drawing_width/2-te.width/2-3, top_margin+drawing_height+4+te.height);
    cairo_show_text (cr, text_x.c_str());

    cairo_text_extents (cr, text_yl.c_str(), &te);
    cairo_move_to(cr, left_margin-te.width-6, top_margin+drawing_height+2);
    cairo_show_text (cr, text_yl.c_str());
    cairo_text_extents (cr, text_yu.c_str(), &te);
    cairo_move_to(cr, left_margin-te.width-6, top_margin+te.height+2);
    cairo_show_text (cr, text_yu.c_str());
    cairo_text_extents (cr, text_y.c_str(), &te);
    cairo_move_to(cr, left_margin-te.width-6, top_margin+drawing_height/2+te.height+2);
    cairo_show_text (cr, text_y.c_str());


    // Save unclipped state and canvas coordinates
    cairo_save (cr);

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

    //std::cerr<<"cw="<<canvas_width<<" lm="<<left_margin<<" dw="<<drawing_width<<" rm="<<right_margin<<" xl="<<xl<<" xu="<<xu<<"\n";
    //std::cerr<<"ch="<<canvas_height<<"tm="<<top_margin<<" dw="<<drawing_height<<" bm="<<bottom_margin<<" yl="<<yl<<" yu="<<yu<<"\n";

    // compute device to user coordinate transformation
    double ctr0=left_margin;
    double ctr1=top_margin;
    double sc0=(drawing_width)/(xu-xl);
    double sc1=-(drawing_height)/(yu-yl);
    double utr0=(-xl);
    double utr1=(-yu);

    // Scale to user coordinates
    cairo_translate(cr, ctr0, ctr1);
    cairo_scale (cr, sc0,sc1);
    cairo_translate(cr, utr0, utr1);
}

Void CairoCanvas::write(const char* filename) const {
    cairo_surface_t* surface = cairo_get_target (cr);
    cairo_surface_write_to_png (surface, filename);
}

Void CairoCanvas::finalise()
{
    cairo_t *cr=this->cr;

    // Restore canvas coordinates and unclipped state
    cairo_restore (cr);

    const ImageSize2d drawing_size = this->size_in_pixels();
    const Int drawing_width = drawing_size.nx;
    const Int drawing_height = drawing_size.ny;

    const Int left_margin = LEFT_MARGIN;
    const Int top_margin = TOP_MARGIN;

    cairo_set_line_width (cr, 2.0);
    cairo_set_source_rgb (cr, 0.0, 0.0, 0.0);
    cairo_move_to (cr, left_margin, top_margin+drawing_height);
    cairo_line_to (cr, left_margin+drawing_width, top_margin+drawing_height);
    cairo_line_to (cr, left_margin+drawing_width, top_margin);
    cairo_line_to (cr, left_margin, top_margin);
    cairo_line_to (cr, left_margin, top_margin+drawing_height);
    cairo_stroke (cr);
}


Void set_properties(CanvasInterface& canvas, const GraphicsProperties& properties) {
    const Colour& line_colour=properties.line_colour;
    const Colour& fill_colour=properties.fill_colour;
    canvas.set_fill_opacity(properties.fill_colour.opacity);
    canvas.set_fill_colour(fill_colour.red, fill_colour.green, fill_colour.blue);
    canvas.set_line_colour(line_colour.red, line_colour.green, line_colour.blue);
}

inline OutputStream& operator<<(OutputStream& os, const ExactBoxType& bx) { return os << static_cast<const ExactIntervalVectorType&>(bx); }

Void Figure::_paint_all(CanvasInterface& canvas) const
{
    ApproximateBoxType bounding_box=this->_data->bounding_box;
    const PlanarProjectionMap projection=this->_data->projection;
    const std::vector<GraphicsObject>& objects=this->_data->objects;

    Nat dimension=projection.argument_size();

    // Don't attempt to compute a bounding box, as this relies on
    // a drawable object having one. Instead, the bounding box must be
    // specified explicitly
    if(bounding_box.dimension()==0) {
        bounding_box=ExactBoxType(dimension,ExactIntervalType(-1,1));
    }

    // Check projection and bounding box have same values.
    ARIADNE_ASSERT_MSG(bounding_box.dimension()==projection.argument_size(),"bounding_box="<<bounding_box<<", projection="<<projection);
    ARIADNE_ASSERT(bounding_box.dimension()>projection.x_coordinate());
    ARIADNE_ASSERT(bounding_box.dimension()>projection.y_coordinate());

    // Project the bounding box onto the canvas
    double xl=numeric_cast<double>(bounding_box[projection.x_coordinate()].lower());
    double xu=numeric_cast<double>(bounding_box[projection.x_coordinate()].upper());
    double yl=numeric_cast<double>(bounding_box[projection.y_coordinate()].lower());
    double yu=numeric_cast<double>(bounding_box[projection.y_coordinate()].upper());

    StringType tx=StringType("x")+str(projection.x_coordinate());
    StringType ty=StringType("x")+str(projection.y_coordinate());
    canvas.initialise(tx,ty,xl,xu,yl,yu);

    // Draw shapes
    for(Nat i=0; i!=objects.size(); ++i) {
        const DrawableInterface& shape=objects[i].shape_ptr.operator*();
        if(shape.dimension()==0) { break; } // The dimension may be equal to two for certain empty sets.
        ARIADNE_ASSERT_MSG(dimension==shape.dimension(),
                           "Shape "<<shape<<", dimension="<<shape.dimension()<<", bounding_box="<<bounding_box);
        set_properties(canvas, objects[i].properties);
        shape.draw(canvas,this->_data->projection);
    }

    canvas.finalise();
}


Void
Figure::write(const char* cfilename) const
{
    this->write(cfilename, DEFAULT_WIDTH, DEFAULT_HEIGHT);
}


Void
Figure::write(const char* cfilename, Nat drawing_width, Nat drawing_height) const
{
    CairoCanvas canvas(ImageSize2d(drawing_width,drawing_height));

    this->_paint_all(canvas);

    StringType filename(cfilename);
    if(filename.rfind(".") != StringType::npos) {
    } else {
        filename=filename+".png";
    }
    
    canvas.write(filename.c_str());
}


#ifdef HAVE_GTK_H

Void
paint (GtkWidget      *widget,
       GdkEventExpose *eev,
       gpointer        gdata)
{
    cairo_t *cr;

    Figure* figure=static_cast<Figure*>(gdata);

    //gint canvas_width  = widget->allocation.width;
    //gint canvas_height = widget->allocation.height;

    // Get Cairo drawing context
    cr = gdk_cairo_create (widget->window);

    // Draw Cairo objects
    CairoCanvas canvas(cr);
    figure->_paint_all(canvas);
}

Void Figure::display() const
{

    GtkWidget *window;
    GtkWidget *canvas;

    Int argc=0;
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

Void Figure::display() const
{
    throw std::runtime_error("No facilities for displaying graphics are available.");
}

#endif // HAVE_GTK_H

#else // NO CAIRO_H

Void
Figure::write(const char* filename) const
{
    throw std::runtime_error("No facilities for drawing graphics are available.");
}

Void Figure::display() const
{
    throw std::runtime_error("No facilities for displaying graphics are available.");
}

#endif // HAVE_CAIRO_H



Colour::Colour()
    : Colour("transparant", 1.0, 1.0, 1.0, 0.0) { }
Colour::Colour(double rd, double gr, double bl, Bool tr)
    : Colour("",rd,gr,bl,tr?0.0:1.0) { }
Colour::Colour(double rd, double gr, double bl, double op)
    : Colour("",rd,gr,bl,op) { }
Colour::Colour(const char* nm, double rd, double gr, double bl, Bool tr)
    : Colour("",rd,gr,bl,tr?0.0:1.0) { }
Colour::Colour(const char* nm, double rd, double gr, double bl, double op)
    : name(nm), red(rd), green(gr), blue(bl), opacity(op) { }
OutputStream& operator<<(OutputStream& os, const Colour& c) {
    return os << "Colour( name=" << c.name << ", r=" << c.red << ", g=" << c.green << ", b=" << c.blue << ", op=" << c.opacity << " )"; }


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


