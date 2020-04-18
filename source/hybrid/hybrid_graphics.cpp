/***************************************************************************
 *            hybrid/hybrid_graphics.cpp
 *
 *  Copyright  2011-20  Pieter Collins
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

#include "../function/functional.hpp"
#include "../config.hpp"

#include "../utility/macros.hpp"
#include "../utility/stlio.hpp"
#include "../numeric/numeric.hpp"
#include "../symbolic/space.hpp"
#include "../geometry/point.hpp"
#include "../geometry/box.hpp"
#include "../output/geometry2d.hpp"
#include "../output/graphics.hpp"
#include "../hybrid/discrete_location.hpp"
#include "../geometry/function_set.hpp"
#include "../symbolic/expression_set.hpp"
#include "../hybrid/hybrid_graphics.hpp"

namespace Ariadne {

static const Nat DEFAULT_WIDTH = 800;
static const Nat DEFAULT_HEIGHT = 800;

static const Nat LEFT_MARGIN = 160;
static const Nat BOTTOM_MARGIN = 40;
static const Nat TOP_MARGIN = 10;
static const Nat RIGHT_MARGIN = 10;

Bool valid_axis_variables(const RealSpace& space, const Variables2d& variables) {
    return ( (variables.x_variable().name()==TimeVariable().name()) || space.contains(variables.x_variable()) ) && space.contains(variables.y_variable());
}

Projection2d projection(const RealSpace& space, const Variables2d& variables) {
    ARIADNE_ASSERT(valid_axis_variables(space,variables));
    Nat x_index = (variables.x_variable()==TimeVariable() && !space.contains(variables.x_variable())) ? space.dimension() : space.index(variables.x_variable());
    Nat y_index = space.index(variables.y_variable());
    return Projection2d(space.dimension(),x_index,y_index);
}


Void paint(CanvasInterface& canvas, const Set<DiscreteLocation>& locations, const Variables2d& variables, const List<HybridGraphicsObject>& objects) {
    for(Nat i=0; i!=objects.size(); ++i) {
        const HybridDrawableInterface& shape=*objects[i].shape_ptr;
        set_properties(canvas, objects[i].properties);
        shape.draw(canvas,locations,variables);
    }
}

HybridFigure::~HybridFigure() {
}


HybridFigure::HybridFigure()
    : variables(RealVariable("x"),RealVariable("y"))
{
}




Void
HybridFigure::write(const char* cfilename) const
{
    this->write(cfilename, DEFAULT_WIDTH, DEFAULT_HEIGHT);
}


Void
HybridFigure::write(const char* cfilename, Nat drawing_width, Nat drawing_height) const
{
    const Nat canvas_width = drawing_width+LEFT_MARGIN+RIGHT_MARGIN;
    const Nat canvas_height = drawing_height+BOTTOM_MARGIN+TOP_MARGIN;

    SharedPointer<CanvasInterface> canvas=make_canvas(canvas_width, canvas_height);

    this->_paint_all(*canvas);

    StringType filename(cfilename);
    if(filename.rfind(".") != StringType::npos) {
    } else {
        filename=filename+".png";
    }

    canvas->write(filename.c_str());
}


Void HybridFigure::_paint_all(CanvasInterface& canvas) const
{
    // Project the bounding box onto the canvas
    double xl=numeric_cast<double>(bounds[variables.x_variable()].lower());
    double xu=numeric_cast<double>(bounds[variables.x_variable()].upper());
    double yl=numeric_cast<double>(bounds[variables.y_variable()].lower());
    double yu=numeric_cast<double>(bounds[variables.y_variable()].upper());

    canvas.initialise(variables.x_variable().name(),variables.y_variable().name(),xl,xu,yl,yu);

    // Draw shapes
    for(Nat i=0; i!=objects.size(); ++i) {
        const HybridDrawableInterface& shape=*objects[i].shape_ptr;
        set_properties(canvas, objects[i].properties);
        shape.draw(canvas,this->locations,this->variables);
    }

    canvas.finalise();
}



} // namespace Ariadne


