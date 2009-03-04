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


std::vector<Point> interpolation_points(const InterpolatedCurve& c) {
    std::vector<Point> result;
    for(InterpolatedCurve::const_iterator iter=c.begin(); iter!=c.end(); ++iter) {
        result.push_back(iter->second);
    }
    return result;
}

std::vector<Point> vertices(const Polytope& pl) {
    return pl.vertices();
}

std::vector<Point> corner_points(const Box& bx) {
    std::vector<Point> result(2,Point(bx.dimension()));
    for(uint i=0; i!=bx.dimension(); ++i) {
        result[0][i]=bx[i].lower();
        result[1][i]=bx[i].upper();
    }
    //std::cerr<<"corner_points("<<bx<<")="<<result<<std::endl;
    return result;
}


Box bounding_box(const std::vector<Point>& points) {
    ARIADNE_ASSERT(!points.empty());
    Box result(points[0]);
    for(uint i=1; i!=points.size(); ++i) {
        const Point& pt=points[i];
        result=hull(result,pt);
    }
    return result;
}


std::vector<Point> apply(const ProjectionFunction& map, const std::vector<Point> & points) {
    //std::cerr<<"apply(Projection map, list<Point> pts): map="<<map<<", pts="<<points<<std::endl;
    std::vector<Point> result;
    for(size_t i=0; i!=points.size(); ++i) {
        const Point& pt=points[i];
        Point ppt(map.result_size());
        for(size_t j=0; j!=map.result_size(); ++j) {
            ppt[j]=pt[map[j]];
        }
        result.push_back(ppt);
    }
    return result;
}


std::vector<Point> extremal(const std::vector<Point> & points) {
    Polytope polytope(points);
    return reduce2d(polytope).vertices();
}
  

struct GraphicsObject {
    enum ShapeKind { POINT, BOX, POLYTOPE, CURVE, SHAPE };
    GraphicsObject(Colour lc, Point pt) : kind(POINT), line_colour(lc), fill_colour(), shape(std::vector<Point>(1,pt)) { }
    GraphicsObject(Colour lc, Colour fc, Box bx) : kind(BOX), line_colour(lc), fill_colour(fc), shape(corner_points(bx)) { }
    GraphicsObject(Colour lc, Colour fc, Polytope pl) : kind(POLYTOPE), line_colour(lc), fill_colour(fc), shape(vertices(pl)) { }
    GraphicsObject(Colour lc, InterpolatedCurve cv) : kind(CURVE), line_colour(lc), fill_colour(), shape(interpolation_points(cv)) { }
    GraphicsObject(Colour lc, Colour fc, const std::vector<Point>& pts) : kind(SHAPE), line_colour(lc), fill_colour(fc), shape(pts) { }
    ShapeKind kind;
    Colour line_colour;
    Colour fill_colour;
    std::vector<Point> shape;
};



struct Figure::Data 
{
    Data() : bounding_box(0), projection(2), current_line_colour(0,0,0), current_fill_colour(1,1,1) { }
    Box bounding_box;
    ProjectionFunction projection;
    Colour current_line_colour;
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
    this->_data->projection=ProjectionFunction(2);
}


void Figure::set_projection_map(const ProjectionFunction& p) 
{
    this->_data->projection=p;
}

void Figure::set_projection_map(const PlanarProjectionMap& p) 
{
    array<uint> ary(2);
    ary[0]=p.i; ary[1]=p.j;
    this->_data->projection=ProjectionFunction(ary,p.n);
}

void Figure::set_bounding_box(const Box& bx) 
{
    ARIADNE_ASSERT(bx.dimension()==0 || bx.dimension()==2 || bx.dimension()==this->_data->projection.argument_size());
    if(bx.dimension()>2) {
        this->_data->bounding_box=this->_data->projection(bx);
    } else {
        this->_data->bounding_box=bx;
    }
}

void Figure::set_line_style(bool ls) 
{
}

void Figure::set_line_width(double lw) 
{
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


void Figure::draw(const Point& pt) {
    if(this->_data->objects.empty() && this->_data->projection.argument_size() != pt.dimension()) {
        this->_data->projection=ProjectionFunction(2,pt.dimension(),0); }
    ARIADNE_ASSERT(pt.dimension()==this->_data->projection.argument_size());
    this->_data->objects.push_back(GraphicsObject(this->_data->current_line_colour,pt));
}

void Figure::draw(const std::vector<Point>& pts) {
    ARIADNE_ASSERT(!pts.empty());
    ARIADNE_ASSERT(pts[0].dimension()==this->_data->projection.argument_size());
    this->_data->objects.push_back(GraphicsObject(this->_data->current_line_colour,this->_data->current_fill_colour,pts));
}

void Figure::draw(const Box& bx) {
    if(this->_data->objects.empty() && this->_data->projection.argument_size() != bx.dimension()) {
        this->_data->projection=ProjectionFunction(2,bx.dimension(),0); }
    ARIADNE_ASSERT(bx.dimension()==this->_data->projection.argument_size());
    this->_data->objects.push_back(GraphicsObject(this->_data->current_line_colour,this->_data->current_fill_colour,bx));
}

void Figure::draw(const Polytope& p) {
    if(this->_data->objects.empty() && this->_data->projection.argument_size() != p.dimension()) {
        this->_data->projection=ProjectionFunction(2,p.dimension(),0); }
    ARIADNE_ASSERT(p.dimension()==this->_data->projection.argument_size());
    this->_data->objects.push_back(GraphicsObject(this->_data->current_line_colour,this->_data->current_fill_colour,p));
}

void Figure::draw(const InterpolatedCurve& c) {
    if(this->_data->objects.empty() && this->_data->projection.argument_size() != c.dimension()) {
        this->_data->projection=ProjectionFunction(2,c.dimension(),0); }
    ARIADNE_ASSERT(c.dimension()==this->_data->projection.argument_size());
    this->_data->objects.push_back(GraphicsObject(this->_data->current_line_colour,c));
}

void Figure::clear() {
    this->_data->objects.clear();
}


#ifdef HAVE_CAIRO_H

static void 
plot(cairo_t *cr, const Box& bounding_box, const ProjectionFunction& projection,
     const std::vector<GraphicsObject>& objects, 
     int canvas_width, int canvas_height);

void trace_point(cairo_t *cr, const std::vector<Point>& pts) 
{
    static const double RADIUS=0.01;
    cairo_arc (cr, pts[0][0], pts[0][1], RADIUS, 0, 2*M_PI);
}

void trace_box(cairo_t *cr, const std::vector<Point>& pts) 
{
    //std::cerr << "trace_box(cairo_t *cr, std::vector<Point> pts) pts="<<pts<<std::endl;
    cairo_move_to (cr, pts[0][0], pts[0][1]);
    cairo_line_to (cr, pts[1][0], pts[0][1]);
    cairo_line_to (cr, pts[1][0], pts[1][1]);
    cairo_line_to (cr, pts[0][0], pts[1][1]);
    cairo_line_to (cr, pts[0][0], pts[0][1]);
}


void trace_polytope(cairo_t *cr, const std::vector<Point>& p) 
{
    //std::cerr << "trace_polytope(cairo_t *cr, std::vector<Point> p) p="<<p<<std::endl;
    ARIADNE_ASSERT(p.size()>=3);
    ARIADNE_ASSERT(p[0].size()==2);
    cairo_move_to (cr, p[0][0], p[0][1]);
    for(uint i=1; i!=p.size(); ++i) {
        cairo_line_to (cr, p[i][0], p[i][1]); 
    }
    cairo_line_to (cr, p[0][0], p[0][1]); 
}

void trace_curve(cairo_t *cr, const std::vector<Point>& cv) 
{
    //std::cerr << "trace_curve(cairo_t *cr, InterpolatedCurve cv) cv="<<cv<<std::endl;
    ARIADNE_ASSERT(cv.size()>=2);
    ARIADNE_ASSERT(cv[0].size()==2);
    cairo_move_to (cr, cv[0][0], cv[0][1]);
    for(uint i=1; i!=cv.size(); ++i) {
        cairo_line_to (cr, cv[i][0], cv[i][1]); 
    }
}

std::string str(double x) {
    std::stringstream ss;
    ss << x;
    return ss.str();
}


void trace(cairo_t *cr, const GraphicsObject::ShapeKind kind, const std::vector<Point>& pts) 
{
    switch(kind) {
    case GraphicsObject::POINT: trace_point(cr,pts); return;
    case GraphicsObject::BOX: trace_box(cr,pts); return;
    case GraphicsObject::POLYTOPE: trace_polytope(cr,extremal(pts)); return;
    case GraphicsObject::CURVE: trace_curve(cr,pts); return;
    case GraphicsObject::SHAPE: trace_polytope(cr,pts); return;
    }
}




void plot(cairo_t *cr, const Box& bounding_box, const ProjectionFunction& projection, 
          const std::vector<GraphicsObject>& objects, int canvas_width, int canvas_height) 
{
    const int left_margin=80;
    const int bottom_margin=40;
    const int top_margin=10;
    const int right_margin=10;

    //std::cerr<<"Plot(...)\n"<<std::flush;

    // Test if there are no objects to be drawn
    if(objects.empty()) {
        cairo_destroy (cr);
        return; 
    }

    // Test if a bounding box has been given explicitly; if not, compute one 
    // Gives the bounding box for the actual used area
    Box bbox;
    if(bounding_box.dimension()==0) {
        if(objects.empty()) { 
            bbox=Box(2,Interval(-1,1));
        } else {
            bbox=Ariadne::bounding_box(apply(projection,objects[0].shape));
            for(uint i=1; i!=objects.size(); ++i) {
                bbox=hull(bbox,Ariadne::bounding_box(apply(projection,objects[i].shape))); 
            }
        }
    } else {
        bbox=bounding_box;
    }


    // The bounding box for the actual used area
    Box lbbox=bbox;
    lbbox[0]+=Interval(-1,1)*(bbox[0].radius()/25);
    lbbox[1]+=Interval(-1,1)*(bbox[1].radius()/25);


    //std::cerr<<"  bbox="<<bbox<<"\n  lbbox="<<lbbox<<"\n  gbbox="<<gbbox<<std::endl;

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
    
    // Fill shapes
    for(uint i=0; i!=objects.size(); ++i) {
        if(objects[i].kind==GraphicsObject::BOX || objects[i].kind==GraphicsObject::POLYTOPE 
            || objects[i].kind==GraphicsObject::SHAPE) 
        {
            trace(cr,objects[i].kind,apply(projection,objects[i].shape));
            const Colour& fc=objects[i].fill_colour;
            cairo_set_source_rgb (cr, fc.red, fc.green, fc.blue);
            cairo_fill (cr);
        }
    }
    
    // Trace curves and points
    for(uint i=0; i!=objects.size(); ++i) {
        const Colour& lc=objects[i].line_colour;
        cairo_set_source_rgb (cr, lc.red, lc.green, lc.blue);
        trace(cr,objects[i].kind,apply(projection,objects[i].shape));
        if(objects[i].kind==GraphicsObject::POINT) {
            cairo_fill (cr);
        }
        cairo_stroke (cr);
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
    //std::cerr<<"  Done draw.\n"<<std::flush;
}


void 
Figure::write(const char* cfilename) 
{
    cairo_surface_t *surface;
    cairo_t *cr;

    const int canvas_width = DEFAULT_WIDTH;
    const int canvas_height = DEFAULT_HEIGHT;

    const Box& bounding_box=this->_data->bounding_box;
    const ProjectionFunction& projection=this->_data->projection;
    const std::vector<GraphicsObject>& objects=this->_data->objects;

    surface = cairo_image_surface_create (CAIRO_FORMAT_ARGB32, canvas_width, canvas_height);
    cr = cairo_create (surface);
 
    Ariadne::plot(cr, bounding_box, projection, objects, canvas_width, canvas_height);
    
    
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
  
    const Figure::Data* data=static_cast<Figure::Data*>(gdata);
    const Box& bounding_box =data->bounding_box;
    const ProjectionFunction& projection =data->projection;
    const std::vector<GraphicsObject>& objects=data->objects;

    gint canvas_width  = widget->allocation.width;
    gint canvas_height = widget->allocation.height;
    
    // Get Cairo drawing context
    cr = gdk_cairo_create (widget->window);

    // Draw Cairo objects
    Ariadne::plot(cr, bounding_box, projection, objects, canvas_width, canvas_height);
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
                      const_cast<Figure::Data*>(this->_data));  //  here we can pass a pointer to a custom data structure 

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

 
