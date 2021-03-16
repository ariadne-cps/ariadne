/***************************************************************************
 *            output/graphics.cpp
 *
 *  Copyright  2008-20  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "utility/standard.hpp"
#include "config.hpp"

#include "utility/macros.hpp"
#include "utility/stlio.hpp"
#include "numeric/numeric.hpp"
#include "function/function.hpp"
#include "geometry/point.hpp"
#include "geometry/box.hpp"
#include "symbolic/variable.hpp"
#include "symbolic/space.hpp"
#include "symbolic/expression_set.hpp"
#include "output/geometry2d.hpp"
#include "output/graphics.hpp"
#include "output/cairo.hpp"
#include "output/progress_indicator.hpp"
#include "output/logging.hpp"

#include "../output/gnuplot.hpp"
//#include "algebra/tensor.hpp"

namespace Ariadne {

static const Int DEFAULT_WIDTH = 800;
static const Int DEFAULT_HEIGHT = 800;

#ifdef HAVE_CAIRO_H

static const Int LEFT_MARGIN = 160;
static const Int BOTTOM_MARGIN = 40;
static const Int TOP_MARGIN = 10;
static const Int RIGHT_MARGIN = 10;

#endif

OutputStream& operator<<(OutputStream& os, const DrawableInterface& drawable) {
    if(auto writable=dynamic_cast<const WritableInterface*>(&drawable)) {
        os << *writable;
    } else {
        os << "Drawable";
    }
    return os;
}

Colour::Colour()
    : Colour("transparant", 1.0, 1.0, 1.0, 0.0) { }
Colour::Colour(Dbl rd, Dbl gr, Dbl bl, Bool tr)
    : Colour("",rd,gr,bl,tr?0.0:1.0) { }
Colour::Colour(Dbl rd, Dbl gr, Dbl bl, Dbl op)
    : Colour("",rd,gr,bl,op) { }
Colour::Colour(const Char* nm, Dbl rd, Dbl gr, Dbl bl, Bool tr)
    : Colour("",rd,gr,bl,tr?0.0:1.0) { }
Colour::Colour(const char* nm, Dbl rd, Dbl gr, Dbl bl, Dbl op)
    : name(nm), red(rd), green(gr), blue(bl), opacity(op) { }
OutputStream& operator<<(OutputStream& os, const Colour& c) {
    return os << "Colour( name=" << c.name << ", r=" << c.red << ", g=" << c.green << ", b=" << c.blue << ", op=" << c.opacity << " )"; }


const Colour transparent=Colour();

const Colour white=Colour("white",1.0,1.0,1.0);
const Colour black=Colour("black",0.0,0.0,0.0);
const Colour red=Colour("red",1.0,0.0,0.0);
const Colour green=Colour("green",0.0,1.0,0.0);
const Colour blue=Colour("blue",0.0,0.0,1.0);
const Colour yellow=Colour("yellow",1.0,1.0,0.0);
const Colour cyan=Colour("cyan",0.0,1.0,1.0);
const Colour magenta=Colour("magenta",1.0,0.0,1.0);
const Colour orange=Colour("orange",1.0,0.75,0.5);
const Colour grey=Colour("grey",0.5,0.5,0.5);
const Colour lightgrey=Colour("lightgrey",0.75,0.75,0.75);
const Colour darkgrey=Colour("darkgrey",0.25,0.25,0.25);

GraphicsProperties& GraphicsProperties::set_dot_radius(Dbl dr) { this->dot_radius=dr; return *this; }
GraphicsProperties& GraphicsProperties::set_line_style(Bool ls) { this->line_style=ls; return *this; }
GraphicsProperties& GraphicsProperties::set_line_width(Dbl lw) { this->line_width=lw; return *this; }
GraphicsProperties& GraphicsProperties::set_line_colour(Colour lc) { this->line_colour=lc; return *this; }
GraphicsProperties& GraphicsProperties::set_line_colour(Dbl r, Dbl g, Dbl b) {
    this->set_line_colour(Colour(r,g,b)); return *this; }
GraphicsProperties& GraphicsProperties::set_fill_style(Bool fs) { this->fill_style=fs; return *this; }
GraphicsProperties& GraphicsProperties::set_fill_opacity(Dbl fo) { this->fill_colour.opacity=fo; return *this; }
GraphicsProperties& GraphicsProperties::set_fill_colour(Colour fc) { this->fill_colour=fc; return *this; }
GraphicsProperties& GraphicsProperties::set_fill_colour(Dbl r, Dbl g, Dbl b) {
    this->set_fill_colour(Colour(r,g,b,this->fill_colour.opacity)); return *this; }
GraphicsProperties& GraphicsProperties::set_3d(Bool dim) { this->is3D = true; return *this; }
GraphicsProperties& GraphicsProperties::set_proj_xy() {this->isProj=true; this->isXY=true; return *this; }
GraphicsProperties& GraphicsProperties::set_proj_xz() {this->isProj=true; this->isXZ=true; return *this; }
GraphicsProperties& GraphicsProperties::set_proj_yz() {this->isProj=true; this->isYZ=true; return *this; }

OutputStream& operator<<(OutputStream& os, GraphicsProperties const& gp) {
    return os << "GraphicsProperties(" << "dot_radius=" << gp.dot_radius
        << ", line_style=" << gp.line_style<<", line_width=" << gp.line_width << ", line_colour=" << gp.line_colour
        << ", fill_style=" << gp.fill_style << ", fill_colour=" << gp.fill_colour << ")"; }

Variables2d::Variables2d(const RealVariable& x, const RealVariable& y) : _x(x.name()), _y(y.name()) { }

Variables3d::Variables3d(const RealVariable& x, const RealVariable& y, const RealVariable& z) : _x(x.name()), _y(y.name()), _z(z.name()) { }

RealVariable Variables2d::x() const { return RealVariable(_x); }

RealVariable Variables2d::y() const { return RealVariable(_y); }

RealVariable Variables2d::x_variable() const { return RealVariable(_x); }

RealVariable Variables2d::y_variable() const { return RealVariable(_y); }

RealVariable Variables3d::x() const { return RealVariable(_x); }

RealVariable Variables3d::y() const { return RealVariable(_y); }

RealVariable Variables3d::z() const { return RealVariable(_z); }

RealVariable Variables3d::x_variable() const { return RealVariable(_x); }

RealVariable Variables3d::y_variable() const { return RealVariable(_y); }

RealVariable Variables3d::z_variable() const { return RealVariable(_z); }

Axes2d::Axes2d(const ApproximateDoubleVariableInterval x, const ApproximateDoubleVariableInterval& y)
        : variables(x.variable(),y.variable()), bounds() {
    bounds.insert(x.variable(),x.interval());
    bounds.insert(y.variable(),y.interval());
}

Axes2d::Axes2d(ApproximateDouble xl, const RealVariable& x, ApproximateDouble xu, ApproximateDouble yl, const RealVariable& y, ApproximateDouble yu)
        : variables(x,y), bounds() {
    bounds.insert(x,ApproximateDoubleInterval(xl,xu));
    bounds.insert(y,ApproximateDoubleInterval(yl,yu));
}

Axes3d::Axes3d(const ApproximateDoubleVariableInterval x, const ApproximateDoubleVariableInterval y, const ApproximateDoubleVariableInterval z)
        : variables3d(x.variable(), y.variable(), z.variable()), bounds() {
    bounds.insert(x.variable(), x.interval());
    bounds.insert(y.variable(), y.interval());
    bounds.insert(z.variable(), z.interval());
}

Axes3d::Axes3d(ApproximateDouble xl, const RealVariable x, ApproximateDouble xu, ApproximateDouble yl, const RealVariable y, ApproximateDouble yu, ApproximateDouble zl, const RealVariable z, ApproximateDouble zu)
        : variables3d(x, y, z), bounds() {
    bounds.insert(x, ApproximateDoubleInterval(xl, xu));
    bounds.insert(y, ApproximateDoubleInterval(yl, yu));
    bounds.insert(z, ApproximateDoubleInterval(zl, zu));
}

Bool valid_axis_variables(const RealSpace& space, const Variables2d& variables) {
    return ( (variables.x_variable().name()==TimeVariable().name()) || space.contains(variables.x_variable()) ) && space.contains(variables.y_variable());
}

Projection2d projection(const RealSpace& space, const Variables2d& variables) {
    ARIADNE_ASSERT(valid_axis_variables(space,variables));
    Nat x_index = (variables.x_variable()==TimeVariable() && !space.contains(variables.x_variable())) ? space.dimension() : space.index(variables.x_variable());
    Nat y_index = space.index(variables.y_variable());
    return Projection2d(space.dimension(),x_index,y_index);
}


Void draw(Figure& fig, const DrawableInterface& shape) {
    fig.draw(shape);
}

Void draw(Figure& fig, FloatDPApproximateBox const& box) {
    fig.draw(box);
}

struct GraphicsObject {
    GraphicsObject(const GraphicsProperties& gp, const DrawableInterface& sh)
        : properties(gp), shape_ptr(sh.clone()) { }
    GraphicsProperties properties;
    std::shared_ptr<const DrawableInterface> shape_ptr;
};

struct LabelledGraphicsObject {
    LabelledGraphicsObject(const GraphicsProperties& gp, const LabelledDrawableInterface& sh)
        : properties(gp), shape_ptr(sh.clone()) { }
    GraphicsProperties properties;
    std::shared_ptr<const LabelledDrawableInterface> shape_ptr;
};

struct Figure::Data
{
    Data() : properties(), projection(2,0,1), projection3d(3,0,1,2),bounding_box(0), variables(RealVariable(""),RealVariable("")), variables3d(RealVariable(""), RealVariable(""), RealVariable("")), tensor2d({0,0},0), tensor3d({0,0,0},0), arrayBound(0), array(0) { }
    Data(ApproximateBoxType bbx, Projection2d prj) : properties(), projection(prj), projection3d(3, prj.x_coordinate(), prj.x_coordinate(), 2),bounding_box(bbx),  variables(RealVariable(""),RealVariable("")), variables3d(RealVariable(""), RealVariable(""), RealVariable("")), tensor2d({0,0},0), tensor3d({0,0,0},0), arrayBound(0), array(0) { }
    Data(ApproximateBoxType bbx, Projection3d prj) : properties(), projection(2, prj.x_coordinate(), prj.y_coordinate()), projection3d(prj), bounding_box(bbx), variables(RealVariable(""), RealVariable("")), variables3d(RealVariable(""), RealVariable(""), RealVariable("")), tensor2d({0,0},0), tensor3d({0,0,0},0), arrayBound(0), array(0) { }

    GraphicsProperties properties;

    Projection2d projection;
    Projection3d projection3d;
    ApproximateBoxType bounding_box;
    List<GraphicsObject> objects;

    Variables2d variables;
    Variables3d variables3d;
    Map<RealVariable,ApproximateDoubleInterval> bounds;
    List<LabelledGraphicsObject> labelled_objects;

    Tensor<2, double> tensor2d;
    Tensor<3, double> tensor3d;
    Array<Array<double>> arrayBound;
    Array<double> array;
};

Figure::~Figure()
{
    delete this->_data;
}

Figure::Figure(const GraphicsBoundingBoxType& bbx, const Projection2d& proj)
    : _data(new Data(bbx,proj))
{
    ARIADNE_ASSERT_MSG(proj.argument_size() == bbx.dimension(), "Coordinate projection "<<proj<<" must take same number of arguments as the dimension of the bounding box "<<bbx);
}

Figure::Figure(const GraphicsBoundingBoxType& bbx, DimensionType ix, DimensionType iy)
    : Figure(bbx,Projection2d(bbx.dimension(),ix,iy))
{
}

Figure::Figure(const GraphicsBoundingBoxType& bbx, Pair<DimensionType,DimensionType> ixy)
    : Figure(bbx,Projection2d(bbx.dimension(),ixy.first,ixy.second))
{
}

Figure::Figure(const GraphicsBoundingBoxType& bbx, const Projection3d& proj)
    : _data(new Data(bbx, proj))
{
    ARIADNE_ASSERT_MSG(proj.argument_size() == bbx.dimension(), "Coordinate projection "<<proj<<" must take same number of arguments as the dimension of the bounding box "<<bbx);
}

Figure::Figure(const GraphicsBoundingBoxType& bbx, DimensionType ix, DimensionType iy, DimensionType iz)
    : Figure(bbx, Projection3d(bbx.dimension(), ix, iy, iz))
{
}

Figure& Figure::draw(const DrawableInterface& shape)
{
    this->_data->objects.push_back(GraphicsObject(this->_data->properties,shape)); return *this;
}

Figure& Figure::set_projection(DimensionType as, DimensionType ix, DimensionType iy)
{
    this->_data->projection=Projection2d(as,ix,iy); return *this;
}

Figure& Figure::set_projection(DimensionType as, DimensionType ix, DimensionType iy, DimensionType iz)
{
    this->_data->projection3d=Projection3d(as, ix,iy,iz); return *this;
}

Figure& Figure::set_projection_map(const Projection2d& p)
{
    this->_data->projection=p; return *this;
}

Figure& Figure::set_projection_map(const Projection3d& p)
{
    this->_data->projection3d=p; return *this;
}

Figure& Figure::set_bounding_box(const ApproximateBoxType& bx)
{
    this->_data->bounding_box=bx; return *this;
}

Projection2d Figure::get_projection_map() const
{
    return this->_data->projection;
}

Projection3d Figure::get_3dprojection_map() const
{
    return this->_data->projection3d;
}

ApproximateBoxType Figure::get_bounding_box() const
{
    return this->_data->bounding_box;
}

GraphicsProperties& Figure::properties() {
    return this->_data->properties;
}

const GraphicsProperties& Figure::properties() const {
    return this->_data->properties;
}

Figure& Figure::set_dot_radius(Dbl dr)
{
    this->_data->properties.dot_radius=dr; return *this;
}

Figure& Figure::set_line_style(Bool ls)
{
    this->_data->properties.line_style=ls; return *this;
}

Figure& Figure::set_line_width(Dbl lw)
{
    this->_data->properties.line_width=lw; return *this;
}

Figure& Figure::set_line_colour(Colour lc)
{
    this->_data->properties.line_colour=lc; return *this;
}

Figure& Figure::set_line_colour(Dbl r, Dbl g, Dbl b)
{
    this->set_line_colour(Colour(r,g,b)); return *this;
}

Figure& Figure::set_fill_style(Bool fs)
{
    this->_data->properties.fill_style=fs; return *this;
}

Figure& Figure::set_fill_opacity(Dbl fo)
{
    this->_data->properties.fill_colour.opacity=fo; return *this;
}

Figure& Figure::set_fill_colour(Colour fc)
{
    this->_data->properties.fill_colour=fc; return *this;
}

Figure& Figure::set_fill_colour(Dbl r, Dbl g, Dbl b)
{
    this->set_fill_colour(Colour(r,g,b,this->_data->properties.fill_colour.opacity)); return *this;
}

Figure& Figure::set_proj_xy()
{
    this->_data->properties.isProj = true;
    this->_data->properties.isXY = true;
    return *this;
}
Figure& Figure::set_proj_xz()
{
    this->_data->properties.isProj = true;
    this->_data->properties.isXZ = true;
    return *this;
}
Figure& Figure::set_proj_yz()
{
    this->_data->properties.isProj = true;
    this->_data->properties.isYZ = true;
    return *this;
}
Bool Figure::get_line_style() const
{
    return this->_data->properties.line_style;
}

Dbl Figure::get_line_width() const
{
    return this->_data->properties.line_width;
}

Dbl Figure::get_dot_radius() const
{
    return this->_data->properties.dot_radius;
}

Colour Figure::get_line_colour() const
{
    return this->_data->properties.line_colour;
}


Bool Figure::get_fill_style() const
{
    return this->_data->properties.fill_style;
}

Dbl Figure::get_fill_opacity() const
{
    return this->_data->properties.fill_colour.opacity;
}

Colour Figure::get_fill_colour() const
{
    return this->_data->properties.fill_colour;
}

Figure& Figure::draw(ApproximateBoxType const& box)
{
    ApproximateBoxSetType box_set(box);
    DrawableInterface const& shape=box_set;
    this->draw(shape); return *this;
}

Figure& Figure::draw(RealBox const& box)
{
    return this->draw(ApproximateBoxType(box,DoublePrecision()));
}

Void Figure::function_to_draw(Tensor<2, double> tensor)
{
    this->_data->tensor2d = tensor;
}

Void Figure::function_to_draw(Tensor<3, double> tensor)
{
    this->_data->tensor3d = tensor;
}

Void Figure::function_to_draw(Array<double> arrayIn)
{
    this->_data->array = arrayIn;
}

Void Figure::function_to_draw(Array<Array<double>> vector)
{
    this->_data->arrayBound = vector;
}

Figure& Figure::clear() {
    this->_data->objects.clear(); return *this;
}

Figure& operator<<(Figure& fig, const DrawableInterface& shape) {
    fig.draw(shape); return fig;
}


Void set_properties(CanvasInterface& canvas, const GraphicsProperties& properties) {
    const Colour& line_colour=properties.line_colour;
    const Colour& fill_colour=properties.fill_colour;
    if (properties.line_style==false) { canvas.set_line_width(0); }
    else { canvas.set_line_width(properties.line_width); }
    if (properties.fill_style==false) { canvas.set_fill_opacity(0); }
    else { canvas.set_fill_opacity(properties.fill_colour.opacity); }
    canvas.set_fill_colour(fill_colour.red, fill_colour.green, fill_colour.blue);
    canvas.set_line_colour(line_colour.red, line_colour.green, line_colour.blue);
}

Void Figure::_paint2d(CanvasInterface& canvas, GnuplotFileType fileType) const
{
    DimensionType dimension=this->_data->projection.argument_size();

    ApproximateBoxType bounding_box=this->_data->bounding_box;
    const Projection2d projection=this->_data->projection;

    if(bounding_box.dimension()==0) {
        bounding_box=ExactBoxType(dimension,ExactIntervalType(-1,1));
    }

    ARIADNE_ASSERT_MSG(bounding_box.dimension()==projection.argument_size(),"bounding_box="<<bounding_box<<", projection="<<projection);
    ARIADNE_ASSERT(bounding_box.dimension()>projection.x_coordinate());
    ARIADNE_ASSERT(bounding_box.dimension()>projection.y_coordinate());

    // Project the bounding box onto the canvas
    Dbl xl=numeric_cast<Dbl>(bounding_box[projection.x_coordinate()].lower_bound());
    Dbl xu=numeric_cast<Dbl>(bounding_box[projection.x_coordinate()].upper_bound());
    Dbl yl=numeric_cast<Dbl>(bounding_box[projection.y_coordinate()].lower_bound());
    Dbl yu=numeric_cast<Dbl>(bounding_box[projection.y_coordinate()].upper_bound());

    String tx=String("x")+to_str(projection.x_coordinate());
    String ty=String("x")+to_str(projection.y_coordinate());

    canvas.initialise(tx,ty,xl,xu,yl,yu);
    set_properties(canvas, this->_data->properties);
}

Void Figure::_paint3d(CanvasInterface& canvas, GnuplotFileType fileType) const
{
    DimensionType dimension=this->_data->projection.argument_size();

    ApproximateBoxType bounding_box=this->_data->bounding_box;
    const Projection3d projection=this->_data->projection3d;

    if(bounding_box.dimension()==0) {
        bounding_box=ExactBoxType(dimension,ExactIntervalType(-1,1));
    }

    ARIADNE_ASSERT_MSG(bounding_box.dimension()==projection.argument_size(),"bounding_box="<<bounding_box<<", projection="<<projection);
    ARIADNE_ASSERT(bounding_box.dimension()>projection.x_coordinate());
    ARIADNE_ASSERT(bounding_box.dimension()>projection.y_coordinate());
    ARIADNE_ASSERT(bounding_box.dimension()>projection.z_coordinate());

    Dbl xl=numeric_cast<Dbl>(bounding_box[projection.x_coordinate()].lower_bound());
    Dbl xu=numeric_cast<Dbl>(bounding_box[projection.x_coordinate()].upper_bound());
    Dbl yl=numeric_cast<Dbl>(bounding_box[projection.y_coordinate()].lower_bound());
    Dbl yu=numeric_cast<Dbl>(bounding_box[projection.y_coordinate()].upper_bound());
    Dbl zl=numeric_cast<Dbl>(bounding_box[projection.z_coordinate()].lower_bound());
    Dbl zu=numeric_cast<Dbl>(bounding_box[projection.z_coordinate()].upper_bound());

    String tx=String("x")+to_str(projection.x_coordinate());
    String ty=String("x")+to_str(projection.y_coordinate());
    String tz=String("x")+to_str(projection.z_coordinate());

    canvas.initialise(tx, ty, tz, xl, xu, yl, yu, zl, zu);
}

Void Figure::_paint_all(CanvasInterface& canvas, GnuplotFileType filetype) const
{
    if (this->_data->properties.is3D == false){
        this->_paint2d(canvas, filetype);
    }   
    else{
        this->_paint3d(canvas, filetype);
    }
}

Void Figure::_paint_all(CanvasInterface& canvas) const
{
    DimensionType dimension=this->_data->projection.argument_size();

    ApproximateBoxType bounding_box=this->_data->bounding_box;
    const Projection2d projection=this->_data->projection;

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
    Dbl xl=numeric_cast<Dbl>(bounding_box[projection.x_coordinate()].lower_bound());
    Dbl xu=numeric_cast<Dbl>(bounding_box[projection.x_coordinate()].upper_bound());
    Dbl yl=numeric_cast<Dbl>(bounding_box[projection.y_coordinate()].lower_bound());
    Dbl yu=numeric_cast<Dbl>(bounding_box[projection.y_coordinate()].upper_bound());

    String tx=String("x")+to_str(projection.x_coordinate());
    String ty=String("x")+to_str(projection.y_coordinate());

    canvas.initialise(tx,ty,xl,xu,yl,yu);

    // Draw shapes
    for(const GraphicsObject& object : this->_data->objects) {
        const DrawableInterface& shape=object.shape_ptr.operator*();
        if(shape.dimension()==0) { break; } // The dimension may be equal to two for certain empty sets.
        ARIADNE_ASSERT_MSG(dimension==shape.dimension(),
                           "Shape "<<shape<<", dimension="<<shape.dimension()<<", bounding_box="<<this->_data->bounding_box);
        set_properties(canvas, object.properties);
        shape.draw(canvas,this->_data->projection);
    }
}

Void
Figure::write(const Char* cfilename) const  //Default PNG
{

    #ifdef HAVE_CAIRO_H
        this->write(cfilename, CairoFileType::PNG);
    #else
    #ifdef HAVE_GNUPLOT_H
        this->write(cfilename, GnuplotFileType::PNG);
    #else
        ARIADNE_ERROR("No facilities for displaying graphics are available.");
    #endif
    #endif

}

Void
Figure::write(const Char* cfilename, CairoFileType fileType) const
{
    this->write(cfilename, DEFAULT_WIDTH, DEFAULT_HEIGHT, fileType);
}

Void
Figure::write(const Char* cfilename, GnuplotFileType fileType) const
{

    this->write(cfilename, DEFAULT_WIDTH, DEFAULT_HEIGHT, fileType);
}

Void
Figure::write(const Char* cfilename, Nat drawing_width, Nat drawing_height, CairoFileType fileType) const
{
    SharedPointer<CanvasInterface> canvas=make_canvas(cfilename, drawing_width,drawing_height, fileType);

    this->_paint_all(*canvas);

    canvas->write(cfilename);
}
  
Void
Figure::write(const Char* cfilename, Nat drawing_width, Nat drawing_height, GnuplotFileType fileType) const
{
    if(this->_data->tensor3d.size(2) == 0){ //If not Tensor 3D
        if(this->_data->tensor2d.size(1) != 0){ //If is Tensor 2D
            if(this->_data->tensor2d.size(1) > 1){  //If is animation GIF
                this->write(cfilename, drawing_width, drawing_height, GnuplotFileType::GIF, this->_data->tensor2d);
            }else if(this->_data->tensor2d.size(1) == 1){   //if is static PNG
                this->write(cfilename, drawing_width, drawing_height, GnuplotFileType::PNG, this->_data->tensor2d);
            }
        }else if (!this->_data->arrayBound.empty()){ //If is array of bounds
            this->write(cfilename, drawing_width,drawing_height, GnuplotFileType::PNG, this->_data->arrayBound);
        }else if (!this->_data->array.empty()){ //If is array
            this->write(cfilename, drawing_width,drawing_height,GnuplotFileType::PNG, this->_data->array);
        }
        else{ 
            SharedPointer<CanvasInterface> canvas=make_canvas(cfilename, drawing_width,drawing_height, fileType);
            this->_paint_all(*canvas); 
            canvas->write(cfilename);
        }
    }else{
        if(this->_data->tensor3d.size(2) > 1){
            this->write(cfilename, drawing_width, drawing_height, GnuplotFileType::GIF, this->_data->tensor3d);
        }
        else if(this->_data->tensor3d.size(2) == 1){
            this->write(cfilename, drawing_width, drawing_height, GnuplotFileType::PNG, this->_data->tensor3d);
        }
    }
}

//2D
Void
Figure::write(const Char* cfilename, Nat drawing_width, Nat drawing_height, GnuplotFileType fileType, Tensor<2, double> tensor) const
{
    SharedPointer<CanvasInterface> canvas=make_canvas(cfilename, drawing_width,drawing_height, fileType);

    this->_paint_all(*canvas, fileType);
    if (this->_data->properties.is3D == false){
        canvas->plot_tensor_2d_image(tensor);
    }
    else{
        ARIADNE_ERROR("Error for trying to plot 3D img in a 2D tensor");
    }

    canvas->write(cfilename);
}

//3D
Void
Figure::write(const Char* cfilename, Nat drawing_width, Nat drawing_height, GnuplotFileType fileType, Tensor<3, double> tensor) const
{
    SharedPointer<CanvasInterface> canvas=make_canvas(cfilename, drawing_width,drawing_height, fileType);
    
    this->_data->properties.set_3d(true);
    this->_paint_all(*canvas, fileType);

    if (this->_data->properties.isProj){
        if (this->_data->properties.isXY){
            canvas->plot_xy_projection(tensor);
        }
        else if (this->_data->properties.isXZ){
            canvas->plot_xz_projection(tensor);
        }
        else{
            canvas->plot_yz_projection(tensor);
        } 
    }
    else if (this->_data->properties.is3D == true){           
        canvas->plot_tensor_3d_image(tensor);
    }
    else{
        ARIADNE_ERROR("Error trying to plot 2D img in a 3D tensor");
    }

    canvas->write(cfilename);
}

Void 
Figure::write(const Char* filename, Nat drawing_width, Nat drawing_height, GnuplotFileType fileType, Array<double> data) const
{
    SharedPointer<CanvasInterface> canvas=make_canvas(filename, drawing_width,drawing_height, fileType);

    this->_paint_all(*canvas, fileType);
    canvas->plot_data(data);

    canvas->write(filename);
}


Void
Figure::write(const Char* filename, Nat drawing_width, Nat drawing_height, GnuplotFileType fileType, Array<Array<double>> bounds) const
{
    SharedPointer<CanvasInterface> canvas=make_canvas(filename, drawing_width,drawing_height, fileType);

    this->_paint_all(*canvas, fileType);
    canvas->plot_bounds(bounds);

    canvas->write(filename);

}

struct LabelledFigure::Data
{
    Data(Axes2d axes)
        : properties(), variables(axes.variables), variables3d(axes.variables.x_variable(),axes.variables.y_variable(), RealVariable(Identifier("z"))), bounds(axes.bounds), tensor2d({0,0},0), tensor3d({0,0,0},0), arrayBound(0), array(0) { }
    Data(const Variables2d& vars, const Map<RealVariable,ApproximateDoubleInterval>& bnds)
        : properties(), variables(vars), variables3d(vars.x_variable(), vars.y_variable(), RealVariable(Identifier("z"))) ,bounds(bnds), tensor2d({0,0},0), tensor3d({0,0,0},0), arrayBound(0), array(0) { }
    Data(Axes3d axes)
        : properties(), variables(axes.variables3d.x_variable(), axes.variables3d.y_variable()) ,variables3d(axes.variables3d), bounds(axes.bounds), tensor2d({0,0},0), tensor3d({0,0,0},0), arrayBound(0), array(0) { }
    Data(const Variables3d& vars, const Map<RealVariable,ApproximateDoubleInterval>& bnds)
        : properties(), variables(vars.x_variable(), vars.y_variable()) ,variables3d(vars), bounds(bnds), tensor2d({0,0},0), tensor3d({0,0,0},0), arrayBound(0), array(0) { }

    GraphicsProperties properties;

    Variables2d variables;
    Variables3d variables3d;
    Map<RealVariable,ApproximateDoubleInterval> bounds;
    List<LabelledGraphicsObject> objects;

    Tensor<2, double> tensor2d;
    Tensor<3, double> tensor3d;
    Array<Array<double>> arrayBound;
    Array<double> array;
};

LabelledFigure::~LabelledFigure() {
    delete this->_data;
}

LabelledFigure::LabelledFigure(const Axes2d& axes)
    : _data(new Data(axes)) {
}

LabelledFigure::LabelledFigure(const Axes3d& axes) 
    : _data(new Data(axes)) {
}

LabelledFigure::LabelledFigure(const Variables2d& vars, const VariablesBox<ApproximateIntervalType>& bbx)
    : _data(new Data(vars,VariablesBox<ApproximateDoubleInterval>(bbx).bounds())) {
}

LabelledFigure::LabelledFigure(const Variables3d& vars, const VariablesBox<ApproximateIntervalType>& bnds)
    : _data(new Data(vars, VariablesBox<ApproximateDoubleInterval>(bnds).bounds())){
}

Void LabelledFigure::set_axes(const Axes2d& axes) {
    _data->bounds=axes.bounds; _data->variables=axes.variables;
}

Void LabelledFigure::set_axes(const Axes3d& axes) {
    _data->bounds=axes.bounds; _data->variables3d=axes.variables3d;
}

Void LabelledFigure::set_bounds(const RealVariable& x, const ApproximateDoubleInterval& ivl) {
    _data->bounds.insert(x,ivl);
}
Void LabelledFigure::set_bounds(const Map<RealVariable,ApproximateDoubleInterval>& b) {
    _data->bounds=b;
}

Void LabelledFigure::set_bounding_box(const VariablesBox<ApproximateIntervalType>& bx) {
    _data->bounds=bx.bounds();
}

GraphicsProperties& LabelledFigure::properties() {
    return this->_data->properties;
}

GraphicsProperties const& LabelledFigure::properties() const {
    return this->_data->properties;
}

LabelledFigure& LabelledFigure::draw(const LabelledDrawableInterface& shape)
{
    this->_data->objects.push_back(LabelledGraphicsObject(this->_data->properties,shape)); return *this;
}

LabelledFigure& LabelledFigure::clear() {
    this->_data->objects.clear(); return *this;
}

LabelledFigure& operator<<(LabelledFigure& fig, const LabelledDrawableInterface& shape) {
    fig.draw(shape); return fig;
}

Void LabelledFigure::function_to_draw(Tensor<2, double> tensor)
{
    this->_data->tensor2d = tensor;
}

Void LabelledFigure::function_to_draw(Tensor<3, double> tensor)
{
    this->_data->tensor3d = tensor;
}

Void LabelledFigure::function_to_draw(Array<double> arrayIn)
{
    this->_data->array = arrayIn;
}

Void LabelledFigure::function_to_draw(Array<Array<double>> vector)
{
    this->_data->arrayBound = vector;
}

Void LabelledFigure::_paint2d(CanvasInterface& canvas, GnuplotFileType filetype) const
{
    auto const& bounds = this->_data->bounds;
    RealVariable const& x=this->_data->variables.x_variable();
    RealVariable const& y=this->_data->variables.y_variable();

    Dbl xl=numeric_cast<Dbl>(bounds[x].lower_bound());
    Dbl xu=numeric_cast<Dbl>(bounds[x].upper_bound());
    Dbl yl=numeric_cast<Dbl>(bounds[y].lower_bound());
    Dbl yu=numeric_cast<Dbl>(bounds[y].upper_bound());

    ARIADNE_ASSERT(x.name()!="" && y.name()!="");

    String tx=x.name();
    String ty=y.name();

    canvas.initialise(tx,ty,xl,xu,yl,yu);

    set_properties(canvas, this->_data->properties);
}

Void LabelledFigure::_paint3d(CanvasInterface& canvas, GnuplotFileType fileType) const
{
    ARIADNE_LOG_SCOPE_CREATE;
    auto const& bounds = this->_data->bounds;
    RealVariable const& x=this->_data->variables.x_variable();
    RealVariable const& y=this->_data->variables.y_variable();
    RealVariable const& z=this->_data->variables3d.z_variable();

    String tx=x.name();
    String ty=y.name();
    String tz=z.name();

    Dbl xl=numeric_cast<Dbl>(bounds[x].lower_bound());
    Dbl xu=numeric_cast<Dbl>(bounds[x].upper_bound());
    Dbl yl=numeric_cast<Dbl>(bounds[y].lower_bound());
    Dbl yu=numeric_cast<Dbl>(bounds[y].upper_bound());
    Dbl zl=numeric_cast<Dbl>(bounds[z].lower_bound());
    Dbl zu=numeric_cast<Dbl>(bounds[z].upper_bound());   

    canvas.initialise(tx, ty, tz, xl, xu, yl, yu, zl, zu);
}

Void LabelledFigure::_paint_all(CanvasInterface& canvas) const
{
    ARIADNE_LOG_SCOPE_CREATE;
    auto const& bounds = this->_data->bounds;
    RealVariable const& x=this->_data->variables.x_variable();
    RealVariable const& y=this->_data->variables.y_variable();
    
    ARIADNE_ASSERT(x.name()!="" && y.name()!="");

    String tx=x.name();
    String ty=y.name();

    Dbl xl=numeric_cast<Dbl>(bounds[x].lower_bound());
    Dbl xu=numeric_cast<Dbl>(bounds[x].upper_bound());
    Dbl yl=numeric_cast<Dbl>(bounds[y].lower_bound());
    Dbl yu=numeric_cast<Dbl>(bounds[y].upper_bound());

    canvas.initialise(tx,ty,xl,xu,yl,yu);
    // Draw shapes
    SizeType total_objects = this->_data->objects.size();
    SizeType processed_objects = 0;
    ProgressIndicator indicator(total_objects);
    ARIADNE_LOG_PRINTLN("Writing " << total_objects << " object" << (total_objects > 1 ? "s..." : "..."));
    for(const LabelledGraphicsObject& object : this->_data->objects) {
        const LabelledDrawableInterface& shape=object.shape_ptr.operator*();
        set_properties(canvas, object.properties);
        shape.draw(canvas,this->_data->variables);
        indicator.update_current(processed_objects++);
        ARIADNE_LOG_SCOPE_PRINTHOLD("[" << indicator.symbol() << "] " << indicator.percentage() << "% ");
    }
    canvas.finalise();
}

Void LabelledFigure::_paint_all(CanvasInterface& canvas, GnuplotFileType fileType) const
{
    if (this->_data->properties.is3D == false){
        this->_paint2d(canvas, fileType);
    }   
    else{
        this->_paint3d(canvas, fileType);
    }
}

Void
LabelledFigure::write(const Char* cfilename) const  //Default PNG on static image
{

    #ifdef HAVE_CAIRO_H
        this->write(cfilename, CairoFileType::PNG);
    #else
        #ifdef HAVE_GNUPLOT_H
            this->write(cfilename, GnuplotFileType::PNG);
        #else
            ARIADNE_ERROR("No facilities for displaying graphics are available.");
        #endif
    #endif
}

Void
LabelledFigure::write(const Char* cfilename, CairoFileType fileType) const  //Cairo PNG on static image
{
        this->write(cfilename, DEFAULT_WIDTH, DEFAULT_HEIGHT, fileType);
}

Void
LabelledFigure::write(const Char* cfilename, GnuplotFileType fileType) const    //Gnuplot PNG or GIF on static image
{
    this->write(cfilename, DEFAULT_WIDTH, DEFAULT_HEIGHT, fileType);
}

Void
LabelledFigure::write(const Char* cfilename, Nat drawing_width, Nat drawing_height, CairoFileType fileType) const
{
    SharedPointer<CanvasInterface> canvas=make_canvas(cfilename, drawing_width,drawing_height, fileType);

    this->_paint_all(*canvas);

    canvas->write(cfilename);
}

Void
LabelledFigure::write(const Char* cfilename, Nat drawing_width, Nat drawing_height, GnuplotFileType fileType) const
{
    if(this->_data->tensor3d.size(2) == 0){ //If not Tensor 3D
        if(this->_data->tensor2d.size(1) != 0){ //If is Tensor 2D
            if(this->_data->tensor2d.size(1) > 1){  //If is animation GIF
                this->write(cfilename, drawing_width, drawing_height, GnuplotFileType::GIF, this->_data->tensor2d);
            }else if(this->_data->tensor2d.size(1) == 1){   //if is static PNG
                this->write(cfilename, drawing_width, drawing_height, GnuplotFileType::PNG, this->_data->tensor2d);
            }
        }else if (!this->_data->arrayBound.empty()){ //If is array of bounds
            this->write(cfilename, drawing_width,drawing_height, GnuplotFileType::PNG, this->_data->arrayBound);
        }else if (!this->_data->array.empty()){ //If is array
            this->write(cfilename, drawing_width,drawing_height,GnuplotFileType::PNG, this->_data->array);
        }
        else{ 
            SharedPointer<CanvasInterface> canvas=make_canvas(cfilename, drawing_width,drawing_height, fileType);
            this->_paint_all(*canvas); 
            canvas->write(cfilename);
        }
    }else{
        if(this->_data->tensor3d.size(2) > 1){
            this->write(cfilename, drawing_width, drawing_height, GnuplotFileType::GIF, this->_data->tensor3d);
        }
        else if(this->_data->tensor3d.size(2) == 1){
            this->write(cfilename, drawing_width, drawing_height, GnuplotFileType::PNG, this->_data->tensor3d);
        }
    }
}

Void
LabelledFigure::write(const Char* cfilename, Nat drawing_width, Nat drawing_height, GnuplotFileType fileType, Tensor<2, double> tensor) const
{
    SharedPointer<CanvasInterface> canvas=make_canvas(cfilename, drawing_width,drawing_height, fileType);


    this->_paint_all(*canvas, fileType);
    if (this->_data->properties.is3D == false)
        canvas->plot_tensor_2d_image(tensor);
    else{
        ARIADNE_ERROR("Error for trying to plot 3D img in a 2D tensor");
    }

    canvas->write(cfilename);

}

Void
LabelledFigure::write(const Char* cfilename, Nat drawing_width, Nat drawing_height, GnuplotFileType fileType, Tensor<3, double> tensor) const
{
    SharedPointer<CanvasInterface> canvas=make_canvas(cfilename, drawing_width,drawing_height, fileType);

    this->_data->properties.set_3d(true);
    this->_paint_all(*canvas, fileType);

    if (this->_data->properties.isProj){
        if (this->_data->properties.isXY){
            canvas->plot_xy_projection(tensor);
        }
        else if (this->_data->properties.isXZ){
            canvas->plot_xz_projection(tensor);
        }
        else{
            canvas->plot_yz_projection(tensor);
        } 
    }
    else if (this->_data->properties.is3D == true){  
        canvas->plot_tensor_3d_image(tensor);
    }
    else{
        ARIADNE_ERROR("Error for trying to plot 2D img in a 3D tensor");
    }

    canvas->write(cfilename);
}

Void 
LabelledFigure::write(const Char* filename, Nat drawing_width, Nat drawing_height, GnuplotFileType fileType, Array<double> data) const
{
    SharedPointer<CanvasInterface> canvas=make_canvas(filename, drawing_width,drawing_height, fileType);

    this->_paint_all(*canvas, fileType);
    canvas->plot_data(data);

    canvas->write(filename);
}


Void
LabelledFigure::write(const Char* filename, Nat drawing_width, Nat drawing_height, GnuplotFileType fileType, Array<Array<double>> bounds) const
{
    SharedPointer<CanvasInterface> canvas=make_canvas(filename, drawing_width,drawing_height, fileType);

    this->_paint_all(*canvas, fileType);
    canvas->plot_bounds(bounds);

    canvas->write(filename);
}

#ifdef HAVE_GNUPLOT_H

SharedPointer<CanvasInterface> make_canvas(const char* cfilename, Nat drawing_width, Nat drawing_height, GnuplotFileType fileType) {
    return std::make_shared<GnuplotCanvas>(cfilename, fileType, drawing_width, drawing_height);
}
#else
SharedPointer<CanvasInterface> make_canvas(const char* filename, Nat drawing_width, Nat drawing_height, GnuplotFileType fileType) {
    ARIADNE_WARN_ONCE("No facilities for displaying graphics are available.");
    return std::make_shared<NullCanvas>();
}
#endif

#ifdef HAVE_CAIRO_H

SharedPointer<CanvasInterface> make_canvas(const char* cfilename, Nat drawing_width, Nat drawing_height, CairoFileType filetype) {
    return std::make_shared<CairoCanvas>(ImageSize2d(drawing_width,drawing_height));
}

#else

SharedPointer<CanvasInterface> make_canvas(const char* filename, Nat drawing_width, Nat drawing_height, CairoFileType fileType) {
    ARIADNE_WARN_ONCE("No facilities for displaying graphics are available.");
    return std::make_shared<NullCanvas>();
}
#endif


#ifdef HAVE_CAIRO_H

CairoCanvas::~CairoCanvas()
{
    cairo_surface_destroy(cairo_get_target(cr));
    cairo_destroy(cr);
}

CairoCanvas::CairoCanvas(cairo_t *c)
    : cr(c), lw(1.0), dr(1.0), lc(0.0,0.0,0.0), fc(1.0,1.0,1.0, 1.0)
{
}


CairoCanvas::CairoCanvas(const ImageSize2d& size)
    : cr(0), lw(1.0), dr(1.0), lc(0.0,0.0,0.0), fc(1.0,1.0,1.0, 1.0)
{
    const Int canvas_width = static_cast<Int>(size.nx)+LEFT_MARGIN+RIGHT_MARGIN;
    const Int canvas_height = static_cast<Int>(size.ny)+BOTTOM_MARGIN+TOP_MARGIN;

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
Void CairoCanvas::set_dot_radius(double radius) { this->dr=radius; }
Void CairoCanvas::set_line_width(double width) { this->lw=width; }
Void CairoCanvas::set_line_colour(double r, double g, double b) { lc.red=r; lc.green=g; lc.blue=b; }
Void CairoCanvas::set_fill_opacity(double o) { fc.opacity=o; }
Void CairoCanvas::set_fill_colour(double r, double g, double b) { fc.red=r; fc.green=g; fc.blue=b; }


// TODO: Use generic canvas routines; move cairo-specific functionality
// into CairoCanvas class.
Void CairoCanvas::initialise(StringType x, StringType y, StringType z, double xl, double xu, double yl, double yu, double lz, double uz) {
    ARIADNE_NOT_IMPLEMENTED;
}
Void CairoCanvas::initialise(StringType text_x, StringType text_y, double xl, double xu, double yl, double yu)
{


    CairoCanvas& cairo_canvas=*this;
    cairo_t *crp=cairo_canvas.cr;

    const ImageSize2d drawing_size = cairo_canvas.size_in_pixels();
    const Int drawing_width = static_cast<Int>(drawing_size.nx);
    const Int drawing_height = static_cast<Int>(drawing_size.ny);

    //const Int canvas_width = cairo_image_surface_get_width(cairo_get_target(cr));
    //const Int canvas_height = cairo_image_surface_get_height(cairo_get_target(cr));

    const Int left_margin = LEFT_MARGIN;
    //const Int right_margin = RIGHT_MARGIN;
    //const Int bottom_margin = BOTTOM_MARGIN;
    const Int top_margin = TOP_MARGIN;

    // clear background
    cairo_set_source_rgb (crp, 1,1,1);
    cairo_paint (crp);

    // Set text font
    cairo_select_font_face (crp, "roman",CAIRO_FONT_SLANT_NORMAL,CAIRO_FONT_WEIGHT_NORMAL);
    cairo_set_font_size (crp, 30);

    // Set text colour
    cairo_set_source_rgb (crp, 0., 0., 0.);

    // Get axis label text
    StringType text_xl=to_str(xl);
    StringType text_xu=to_str(xu);
    StringType text_yl=to_str(yl);
    StringType text_yu=to_str(yu);

    // Write axis labels
    cairo_text_extents_t te;
    cairo_text_extents (crp, text_xl.c_str(), &te);
    cairo_move_to(crp, left_margin-2, top_margin+drawing_height+4+te.height);
    cairo_show_text (crp, text_xl.c_str());
    cairo_text_extents (crp, text_xu.c_str(), &te);
    cairo_move_to(crp, left_margin+drawing_width-te.width-4, top_margin+drawing_height+4+te.height);
    cairo_show_text (crp, text_xu.c_str());
    cairo_text_extents (crp, text_x.c_str(), &te);
    cairo_move_to(crp, left_margin+drawing_width/2-te.width/2-3, top_margin+drawing_height+4+te.height);
    cairo_show_text (crp, text_x.c_str());

    cairo_text_extents (crp, text_yl.c_str(), &te);
    cairo_move_to(crp, left_margin-te.width-6, top_margin+drawing_height+2);
    cairo_show_text (crp, text_yl.c_str());
    cairo_text_extents (crp, text_yu.c_str(), &te);
    cairo_move_to(crp, left_margin-te.width-6, top_margin+te.height+2);
    cairo_show_text (crp, text_yu.c_str());
    cairo_text_extents (crp, text_y.c_str(), &te);
    cairo_move_to(crp, left_margin-te.width-6, top_margin+drawing_height/2+te.height+2);
    cairo_show_text (crp, text_y.c_str());


    // Save unclipped state and canvas coordinates
    cairo_save (crp);

    // Set clipping region
    cairo_move_to (crp, left_margin, top_margin+drawing_height);
    cairo_line_to (crp, left_margin+drawing_width, top_margin+drawing_height);
    cairo_line_to (crp, left_margin+drawing_width, top_margin);
    cairo_line_to (crp, left_margin, top_margin);
    cairo_line_to (crp, left_margin, top_margin+drawing_height);

    // Fill clipping region with a very light colour to indicate where figure
    // should be drawn
    cairo_set_source_rgb(crp, 0.97,0.97,0.97);
    cairo_fill_preserve (crp);
    cairo_clip (crp);
    cairo_new_path (crp);

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
    cairo_translate(crp, ctr0, ctr1);
    cairo_scale (crp, sc0,sc1);
    cairo_translate(crp, utr0, utr1);
}

Void CairoCanvas::write(const char* cfilename) const {
    StringType filename(cfilename);
    if(filename.rfind(".") != StringType::npos) {
    } else {
        filename=filename+".png";
    }
    cairo_surface_t* surface = cairo_get_target (cr);
    cairo_surface_write_to_png (surface, filename.c_str());
}

Void CairoCanvas::finalise()
{
    cairo_t *crp=this->cr;

    // Restore canvas coordinates and unclipped state
    cairo_restore (crp);

    const ImageSize2d drawing_size = this->size_in_pixels();
    const Int drawing_width = static_cast<Int>(drawing_size.nx);
    const Int drawing_height = static_cast<Int>(drawing_size.ny);

    const Int left_margin = LEFT_MARGIN;
    const Int top_margin = TOP_MARGIN;

    cairo_set_line_width (crp, 2.0);
    cairo_set_source_rgb (crp, 0.0, 0.0, 0.0);
    cairo_move_to (crp, left_margin, top_margin+drawing_height);
    cairo_line_to (crp, left_margin+drawing_width, top_margin+drawing_height);
    cairo_line_to (crp, left_margin+drawing_width, top_margin);
    cairo_line_to (crp, left_margin, top_margin);
    cairo_line_to (crp, left_margin, top_margin+drawing_height);
    cairo_stroke (crp);
}

Void CairoCanvas::plot_data(Array<double> data)                  { ARIADNE_NOT_IMPLEMENTED; }
Void CairoCanvas::plot_tensor_2d_image(Tensor<2, double> tensor)   { ARIADNE_NOT_IMPLEMENTED; }
Void CairoCanvas::plot_tensor_3d_image(Tensor<3, double> tensor)   { ARIADNE_NOT_IMPLEMENTED; }
Void CairoCanvas::plot_xy_projection(Tensor<3, double> tensor)    { ARIADNE_NOT_IMPLEMENTED; } 
Void CairoCanvas::plot_xz_projection(Tensor<3, double> tensor)    { ARIADNE_NOT_IMPLEMENTED; } 
Void CairoCanvas::plot_yz_projection(Tensor<3, double> tensor)    { ARIADNE_NOT_IMPLEMENTED; } 
Void CairoCanvas::plot_bounds(Array<Array<double>> bounds)       { ARIADNE_NOT_IMPLEMENTED; }

#endif

#ifdef HAVE_GNUPLOT_H

//CANVAS
GnuplotCanvas::GnuplotCanvas(String cfilename, GnuplotFileType typeFile, Nat X, Nat Y): lc(0.0, 0.0, 0.0, 0.0),
                                            fc(1.0, 1.0, 1.0, 1.0),
                                            lw(1.0),
                                            dr(1.0),
                                            isdot(false),
                                            sizeX(X),
                                            sizeY(Y),
                                            isMultiplot(false),
                                            is2DPalette(false),
                                            is3DPalette(false)

{
    gnuplot = new Gnuplot("tee "+cfilename+".gnu | gnuplot > /dev/null 2>&1");
    
    if(typeFile == GnuplotFileType::PNG){ 
        *gnuplot << "set terminal png ";
    }
    else if (typeFile == GnuplotFileType::GIF){ 
        *gnuplot << "set terminal gif animate ";
    }

    *gnuplot << "size " << to_string(this->sizeX) << ", " <<
        to_string(this->sizeY);

    *gnuplot << "\n";

    if(typeFile == GnuplotFileType::PNG){
        *gnuplot << "set output \"" << cfilename << ".png\"\n";
        this->set_multiplot(true);
    }
    else if(typeFile == GnuplotFileType::GIF){
        *gnuplot << "set output \"" << cfilename << ".gif\"\n";
        *gnuplot << "unset multiplot\n";
        this->set_multiplot(false);
    }
    
    this->geom.resize(1024);
    this->dim = 0;

}
//Set Label and Range
void GnuplotCanvas::initialise(StringType x, StringType y, StringType z, double xl, double xu, double yl, double yu, double zl, double zu)
{
    this->set_x_label(x);
    this->set_y_label(y);
    this->set_z_label(z);
    this->set_range_3d(xl, xu, yl, yu, zl, zu);
}
//Set Label and Range
void GnuplotCanvas::initialise(StringType x, StringType y, double xl, double xu, double yl, double yu)
{
    this->set_x_label(x);
    this->set_y_label(y);
    this->set_range_2d(xl, xu, yl, yu);
}

void GnuplotCanvas::finalise() {}
void GnuplotCanvas::circle(double x, double y, double r) {}
void GnuplotCanvas::stroke() {}

void GnuplotCanvas::move_to(double x, double y)
{   
    this->Cpoint.x = this->geom[0].x = x;
    this->Cpoint.y = this->geom[0].y = y;
}

void GnuplotCanvas::line_to(double x, double y)
{
    this->dim++;

    if (this->dim >= this->geom.size())
    {
        this->geom.resize(this->geom.size()+this->geom.size());
        this->geom.resize(this->geom.size()+this->geom.size());
    }
  
    this->Cpoint.x = this->geom[this->dim].x = x;
    this->Cpoint.y = this->geom[this->dim].y = y;

}

void GnuplotCanvas::dot(double x, double y)
{
    this->isdot = true;
    this->Cpoint.x = x;
    this->Cpoint.y = y;
}

void GnuplotCanvas::fill()
{   
    //this->setMultiplot(true);
    char hex_string[20];
    if (this->isdot)
    {
        *gnuplot << "plot \"<echo '" << to_string(this->Cpoint.x) << " " << to_string(this->Cpoint.y) << "'\" w p ls 7 ps " << to_string(this->dr) << "\n";
    }
    else
    {
        *gnuplot << "plot '-' w filledcurves ";
        *gnuplot << "fc rgb \"#";
        if (this->fc.red < 9) { *gnuplot << "0" << this->fc.red;}
        else if (this->fc.red > 255){ *gnuplot << "FF";}
        else{   
            sprintf(hex_string, "%X", std::make_unsigned<int>::type(this->fc.red)); 
            *gnuplot << hex_string;
            }
        if (this->fc.green < 9) { *gnuplot << "0" << this->fc.green;}
        else if (this->fc.green > 255){ *gnuplot << "FF";}
        else{
            sprintf(hex_string, "%X", std::make_unsigned<int>::type(this->fc.green)); 
            *gnuplot << hex_string;
            }
        if (this->fc.blue < 9) { *gnuplot << "0" << this->fc.blue;}
        else if (this->fc.blue > 255){ *gnuplot << "FF";}
        else{sprintf(hex_string, "%X", std::make_unsigned<int>::type(this->fc.blue)); 
            *gnuplot <<hex_string;
            }
        *gnuplot << "\"\n";

        for (SizeType i = 0; i < this->dim; i++)
        {
            *gnuplot << to_string(this->geom[i].x) << " " << to_string(this->geom[i].y) << "\n";
        }
        *gnuplot << "e\n"; 

        if(this->lw == 0){}
        else{
            *gnuplot << "plot '-' w lines";
            *gnuplot << " lc rgb\"#";
            if (this->lc.red < 9) { *gnuplot << "0" << this->lc.red;}
            else if (this->lc.red > 255){ *gnuplot << "FF";}
            else{   
                sprintf(hex_string, "%X", std::make_unsigned<int>::type(this->lc.red)); 
                *gnuplot << hex_string;
                }
            if (this->lc.green < 9) { *gnuplot << "0" << this->lc.green;}
            else if (this->lc.green > 255){ *gnuplot << "FF";}
            else{
                sprintf(hex_string, "%X", std::make_unsigned<int>::type(this->lc.green)); 
                *gnuplot << hex_string;
                }
            if (this->lc.blue < 9) { *gnuplot << "0" << this->lc.blue;}
            else if (this->lc.blue > 255){ *gnuplot << "FF";}
            else{sprintf(hex_string, "%X", std::make_unsigned<int>::type(this->lc.blue)); 
                *gnuplot <<hex_string;
                }
            *gnuplot << "\" ";

            *gnuplot << "lw " << to_string(this->lw);

            *gnuplot << "\n";
            for (SizeType i = 0; i < this->dim; i++)
            {
                *gnuplot << to_string(this->geom[i].x) << " " << to_string(this->geom[i].y) << "\n";
            }
            *gnuplot << to_string(this->geom[0].x) << " " << to_string(this->geom[0].y) << "\n";
            *gnuplot << "e\n"; 
        }

        this->dim = 0;   
    }   
}

void GnuplotCanvas::write(const char* filename) const
{
    *gnuplot << "quit\n";
}

void GnuplotCanvas::set_dot_radius(double _dr)
{
    this->dr = _dr;
}

void GnuplotCanvas::set_line_width(double _lw)
{
    this->lw = _lw;
}

void GnuplotCanvas::set_line_colour(double _r, double _g, double _b) 
{
    this->lc.red = std::round(_r*255);
    this->lc.green = std::round(_g*255);
    this->lc.blue = std::round(_b*255); 
}

void GnuplotCanvas::set_fill_opacity(double _fo)
{
    this->fc.opacity = _fo;

    *gnuplot << "set style fill transparent solid " << to_string(_fo) << "\n";
}

void GnuplotCanvas::set_fill_colour(double _r, double _g, double _b) 
{
    this->fc.red = std::round(_r*255);
    this->fc.green = std::round(_g*255);
    this->fc.blue = std::round(_b*255);
}

Void GnuplotCanvas::set_3d_palette()
{
    is3DPalette = true;
    *gnuplot << "set cbrange [" << to_string(-0.5) << ":" << to_string(1) << "]\n";
    *gnuplot << "set cbtics " << to_string(0.2) << "\n";
    *gnuplot << "set palette defined\n";
}

Void GnuplotCanvas::plot_tensor_2d_image(Tensor<2, double> tensor)
{
    Array<double> data(tensor.size(0), 0);
    for (SizeType step = 0; step < tensor.size(1); step++)
    {
        for (SizeType x = 0; x < data.size(); x++)
        {
            data[x] = (tensor[{x, step}]);
        }
        plot_2d(data);
    }
}

Void GnuplotCanvas::plot_tensor_3d_image(Tensor<3, double> tensor)
{
    this->set_3d_palette();
    SizeType dimX = tensor.size(0);
    SizeType dimY = tensor.size(1);
    SizeType dimTime = tensor.size(2);
    Array<Array<double>> data(dimY);
    for(SizeType step = 0; step < dimTime; step++)
    {
        for (SizeType i = 0; i < dimY; i++)
        {
            data[i].resize(dimX);
            for (SizeType j = 0; j < dimX; j++)
            {
                data[i].at(j) = tensor[{j, i, step}];
            }    
        }
        plot_3d(data);
    }  
}

Void GnuplotCanvas::plot_xy_projection(Tensor<3, double> tensor)
{
    this->set_map();
    this->initialise(this->labels.xLabel, this->labels.yLabel, this->rng.Xmin, this->rng.Xmax, this->rng.Ymin, this->rng.Ymax);
    Array<Array<double>> data(tensor.size(1));
    this->plot_tensor_3d_image(tensor);
}

Void GnuplotCanvas::plot_yz_projection(Tensor<3, double> tensor)
{
    this->initialise(this->labels.yLabel, this->labels.zLabel, this->rng.Ymin, this->rng.Ymax, this->rng.Zmin, this->rng.Zmax);
    set_2d_palette(this->rng.Ymin, this->rng.Ymax, 0.2);
    Array<double> data(tensor.size(1), 0);
    for (SizeType step = 0; step < tensor.size(2); step++)
    { 
        set_multiplot(true);
        for (SizeType x = 0; x < tensor.size(0); x++)
        {
            for (SizeType y = 0; y < tensor.size(1); y++)
            {
                data[y] = tensor[{x, y, step}];
            }
            plot_2d(data);    
        }
    }
}

Void GnuplotCanvas::plot_xz_projection(Tensor<3, double> tensor)
{
    this->initialise(this->labels.xLabel, this->labels.zLabel, this->rng.Xmin, this->rng.Xmax, this->rng.Zmin, this->rng.Zmax);
    set_2d_palette(this->rng.Zmin, this->rng.Zmax, 0.2);
    Array<double> data(tensor.size(0), 0);
    for (SizeType step = 0; step < tensor.size(2); step++)
    { 
        set_multiplot(true);
        for (SizeType y = 0; y < tensor.size(1); y++)
        {
            for (SizeType x = 0; x < tensor.size(0); x++)
            {
                data[x] = tensor[{x, y, step}];
            }
            plot_2d(data);    
        }
    }
}

Void GnuplotCanvas::plot_bounds(Array<Array<double>> bounds)
{
    this->set_multiplot(true);
    this->plot_2d(bounds);
}

Void GnuplotCanvas::plot_data(Array<double> data)
{
    this->plot_2d(data);
}

Vector2d GnuplotCanvas::scaling() const { return Vector2d(0, 0); }
Box2d GnuplotCanvas::bounds() const { return Box2d(0, 0, 0, 0); }

void GnuplotCanvas::set_multiplot(bool s)
{

    if (this->isMultiplot == s){ }
    else{
        this->isMultiplot = s;
        if (s){ *gnuplot << "set multiplot\n"; }
        else { *gnuplot << "unset multiplot\n"; }
    }
}

void GnuplotCanvas::set_multiplot_layout(int nRow, int nCol, String title)
{
    *gnuplot << "set multiplot layout " << to_string(nRow) << "," << to_string(nCol) << " title \"" << title << "\"\n";
}

Void GnuplotCanvas::plot_2d(Array<double> data)
{
    char hex_string[20];
    // START PLOT SINTAX
    *gnuplot << "plot ";
    // Gnuplot wait for input
    *gnuplot << "'-' ";
    // set colour
    if (is2DPalette)
    {
        *gnuplot << " u ::1 " << "with lines lw " << to_string(this->lw) << " linecolor palette";
    }
    else
    {
        *gnuplot << "with lines lw " << to_string(this->lw) << " lc rgb\"#";
        if (this->lc.red < 9) { *gnuplot << "0" << this->lc.red;}
        else if (this->lc.red > 255){ *gnuplot << "FF";}
        else{   
            sprintf(hex_string, "%X", std::make_unsigned<int>::type(this->lc.red)); 
            *gnuplot << hex_string;
            }
        if (this->lc.green < 9) { *gnuplot << "0" << this->lc.green;}
        else if (this->lc.green > 255){ *gnuplot << "FF";}
        else{
            sprintf(hex_string, "%X", std::make_unsigned<int>::type(this->lc.green)); 
            *gnuplot << hex_string;
            }
        if (this->lc.blue < 9) { *gnuplot << "0" << this->lc.blue;}
        else if (this->lc.blue > 255){ *gnuplot << "FF";}
        else{sprintf(hex_string, "%X", std::make_unsigned<int>::type(this->lc.blue)); 
            *gnuplot <<hex_string;
            }
        *gnuplot << "\" ";
    }  

    *gnuplot << "\n";

    //Send data through pipe gp
    gnuplot->send1d(data);
}

void GnuplotCanvas::plot_2d(Array<Array<double>> dataBound)
{
    char hex_string[20];

    for (SizeType i = 0; i < 2; i++)    //Lower = 0 and Upper = 1 value
    {
        // START PLOT SINTAX
        *gnuplot<< "plot ";

        // Gnuplot wait for input
        *gnuplot << "'-' ";
        if (is2DPalette){
            *gnuplot << " u ::1 " << "with lines lw " << to_string(this->lw) << "linecolor palette";
        }
        else
        {
            *gnuplot << "with lines lw " << to_string(this->lw) << " lc rgb\"#";
            if (this->lc.red < 9) { *gnuplot << "0" << this->lc.red;}
            else if (this->lc.red > 255){ *gnuplot << "FF";}
            else{   
                sprintf(hex_string, "%X", std::make_unsigned<int>::type(this->lc.red)); 
                *gnuplot << hex_string;
                }
            if (this->lc.green < 9) { *gnuplot << "0" << this->lc.green;}
            else if (this->lc.green > 255){ *gnuplot << "FF";}
            else{
                sprintf(hex_string, "%X", std::make_unsigned<int>::type(this->lc.green)); 
                *gnuplot << hex_string;
                }
            if (this->lc.blue < 9) { *gnuplot << "0" << this->lc.blue;}
            else if (this->lc.blue > 255){ *gnuplot << "FF";}
            else{sprintf(hex_string, "%X", std::make_unsigned<int>::type(this->lc.blue)); 
                *gnuplot <<hex_string;
                }
            *gnuplot << "\" ";
        }
        *gnuplot<< "\n";
        //Send data through pipe gp
        gnuplot->send1d(dataBound[i]);
    }
}

Void GnuplotCanvas::plot_3d(Array<Array<double>> data)
{  
    // START PLOT SINTAX
    *gnuplot << "splot ";
    // Gnuplot wait for input 
    *gnuplot << "'-' ";
    // get linestyle
    *gnuplot << "with pm3d ";

    *gnuplot<< "\n";
    
    //Send data through pipe gp
    gnuplot->send2d(data); 
}

void GnuplotCanvas::set_x_label(String _xLabel)
{
    *gnuplot << "set xlabel '" << _xLabel << "'\n";
    this->labels.xLabel = _xLabel;
}

void GnuplotCanvas::set_y_label(String _yLabel)
{
    *gnuplot << "set ylabel '" << _yLabel << "'\n";
    this->labels.yLabel = _yLabel;
}

void GnuplotCanvas::set_z_label(String _zLabel)
{
    *gnuplot << "set zlabel '" << _zLabel << "'\n";
    this->labels.zLabel = _zLabel;
}

void GnuplotCanvas::set_title(String title)
{
    *gnuplot << "set title '" << title << "'\n";
}

void GnuplotCanvas::set_xyz_label(String xLabel, String yLabel, String zLabel = "")
{
    *gnuplot << "set xlabel '" << xLabel << "'\n";
    *gnuplot << "set ylabel '" << yLabel << "'\n";
    
    if (zLabel != "")
    {
        *gnuplot << "set zlabel '" << zLabel << "'\n";
    }

}

void GnuplotCanvas::set_labels(String xLabel, String yLabel, String zLabel, String title)
{
    *gnuplot << "set xlabel '" << xLabel << "'\n";
    *gnuplot << "set ylabel '" << yLabel << "'\n";
    *gnuplot << "set zlabel '" << zLabel << "'\n";
    *gnuplot << "set title '" << title << "'\n";
}

void GnuplotCanvas::set_range_2d(double minX, double maxX, double minY, double maxY)
{
    *gnuplot << "set xrange [" << to_string(minX) <<":" << to_string(maxX) << "] \n";
    this->rng.Xmin = minX;
    this->rng.Xmax = maxX;
    *gnuplot << "set yrange [" << to_string(minY) <<":" << to_string(maxY) << "] \n";
    this->rng.Ymin = minY;
    this->rng.Ymax = maxY;
}

void GnuplotCanvas::set_range_3d(double minX, double maxX, double minY, double maxY, double minZ, double maxZ)
{
    *gnuplot << "set xrange [" << to_string(minX) <<":" << to_string(maxX) << "] \n";
    this->rng.Xmin = minX;
    this->rng.Xmax = maxX;
    *gnuplot << "set yrange [" << to_string(minY) <<":" << to_string(maxY) << "] \n";
    this->rng.Ymin = minY;
    this->rng.Ymax = maxY;
    *gnuplot << "set zrange [" << to_string(minZ) <<":" << to_string(maxZ) << "] \n";
    this->rng.Zmin = minZ;
    this->rng.Zmax = maxZ;
}

void GnuplotCanvas::set_x_log_axis()
{
    *gnuplot << "set logscale x\n";
}

void GnuplotCanvas::set_y_log_axis()
{
    *gnuplot << "set logscale y\n";
}

void GnuplotCanvas::set_xy_log_axis()
{
    *gnuplot << "set logscale xy\n";
}

void GnuplotCanvas::set_xz_log_axis()
{
    *gnuplot << "set logscale xz\n";
}

void GnuplotCanvas::set_yz_log_axis()
{
    *gnuplot << "set logscale yz\n";
}

void GnuplotCanvas::set_xyz_log_axis()
{
    *gnuplot << "set logscale xyz\n";
}

void GnuplotCanvas::set_legend()
{
    *gnuplot << "set key default\n";
}
void GnuplotCanvas::set_map()
{
    //*gnuplot << "set pm3d map\n";
    *gnuplot << "set view map\n";
}

void GnuplotCanvas::set_2d_palette(double min, double max, double step)
{
    is2DPalette = true;
    *gnuplot << "set cbrange [" << to_string(min) << ":" << to_string(max) << "]\n";
    *gnuplot << "set cbtics " << to_string(step) << "\n";
    *gnuplot << "set palette defined\n";    
}

void GnuplotCanvas::unset_color_box()
{
    *gnuplot << "unset colorbox\n";
}

Void plot(const char* filename, const Projection2d& pr, const ApproximateBoxType& bbox, List<Pair<Colour,DrawableInterface const&>> const& csets) {
    Figure fig(bbox,pr);
    for (auto cset : csets) {
        fig.set_fill_colour(cset.first); fig << cset.second;
    }
    fig.write(filename);
}

#endif // DEBUG

} // namespace Ariadne


