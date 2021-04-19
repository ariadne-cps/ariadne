/***************************************************************************
 *            graphics_submodule.cpp
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

#include "pybind11.hpp"
#include "utilities.hpp"

#include "config.hpp"

#include "output/graphics_interface.hpp"
#include "output/graphics.hpp"
#include "output/geometry2d.hpp"
#include "geometry/point.hpp"
#include "geometry/box.hpp"
#include "function/function.hpp"

using namespace Ariadne;

Void export_figure(pybind11::module& module)
{
    pybind11::class_<Projection2d> planar_projection_map_class(module,"Projection2d");
    planar_projection_map_class.def(pybind11::init<DimensionType,DimensionType,DimensionType>());

    pybind11::class_<Axes2d> axes2d_class(module,"Axes2d");
    axes2d_class.def(pybind11::init<double,RealVariable,double,double,RealVariable,double>());

    static constexpr auto reference_internal = pybind11::return_value_policy::reference_internal ;

    pybind11::class_<Colour> colour_class(module,"Colour");
    colour_class.def(pybind11::init<double,double,double>());

    pybind11::class_<Figure> figure_class(module,"Figure");
    figure_class.def(pybind11::init<GraphicsBoundingBoxType,Projection2d>());
    figure_class.def("set_projection_map",(Figure&(Figure::*)(const Projection2d&)) &Figure::set_projection_map, reference_internal);
    figure_class.def("set_projection",(Figure&(Figure::*)(DimensionType,DimensionType,DimensionType)) &Figure::set_projection, reference_internal);
    figure_class.def("set_bounding_box",&Figure::set_bounding_box, reference_internal);
    figure_class.def("set_dot_radius", (Figure&(Figure::*)(double)) &Figure::set_dot_radius, reference_internal);
    figure_class.def("set_line_style", (Figure&(Figure::*)(Bool)) &Figure::set_line_style, reference_internal);
    figure_class.def("set_line_width", (Figure&(Figure::*)(double)) &Figure::set_line_width, reference_internal);
    figure_class.def("set_line_colour", (Figure&(Figure::*)(double,double,double)) &Figure::set_line_colour, reference_internal);
    figure_class.def("set_fill_style", (Figure&(Figure::*)(Bool)) &Figure::set_fill_style, reference_internal);
    figure_class.def("set_fill_colour", (Figure&(Figure::*)(double,double,double)) &Figure::set_fill_colour, reference_internal);
    figure_class.def("set_fill_opacity", (Figure&(Figure::*)(double)) &Figure::set_fill_opacity, reference_internal);
    figure_class.def("draw",(Figure&(Figure::*)(const DrawableInterface&))&Figure::draw, reference_internal);
    figure_class.def("draw",(Figure&(Figure::*)(const RealBox&))&Figure::draw, reference_internal);
    figure_class.def("draw",(Figure&(Figure::*)(const ApproximateBoxType&))&Figure::draw, reference_internal);
    figure_class.def("clear",&Figure::clear, reference_internal);
    figure_class.def("write",(Void(Figure::*)(const Char*)const)&Figure::write);
    figure_class.def("write",(Void(Figure::*)(const Char*,Nat,Nat)const)&Figure::write);
}

Void export_plot(pybind11::module& module)
{
    module.def("plot",(Void(*)(const char*,Projection2d const&,ApproximateBoxType const&,List<Pair<Colour,DrawableInterface const&>> const&)) &plot);
}

Void graphics_submodule(pybind11::module& module) {
    export_figure(module);
    export_plot(module);
}

