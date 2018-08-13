/***************************************************************************
 *            graphics_submodule.cpp
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

#include "boost_python.hpp"
#include "utilities.hpp"

#include "config.hpp"

#include <boost/python.hpp>

#include "output/graphics_interface.hpp"
#include "output/graphics.hpp"
#include "output/geometry2d.hpp"
#include "geometry/point.hpp"
    #include "geometry/box.hpp"
#include "function/function.hpp"

using namespace boost::python;

using namespace Ariadne;


Void export_figure()
{
    class_<PlanarProjectionMap>("PlanarProjectionMap",init<Nat,Nat,Nat>());

    return_value_policy<reference_existing_object> ref_existing;

    class_<FigureInterface,boost::noncopyable>("FigureInterface",no_init);
    class_<Figure, bases<FigureInterface> > figure_class("Figure",init<>());

    //class_<Figure, bases<FigureInterface> > figure_class("Figure",init<>());
    figure_class.def("set_projection_map",(Figure&(Figure::*)(const PlanarProjectionMap&)) &Figure::set_projection_map, ref_existing);
    figure_class.def("set_projection",(Figure&(Figure::*)(Nat,Nat,Nat)) &Figure::set_projection, ref_existing);
    figure_class.def("set_bounding_box",&Figure::set_bounding_box, ref_existing);
    figure_class.def("set_dot_radius", (Figure&(Figure::*)(double)) &Figure::set_dot_radius, ref_existing);
    figure_class.def("set_line_style", (Figure&(Figure::*)(Bool)) &Figure::set_line_style, ref_existing);
    figure_class.def("set_line_width", (Figure&(Figure::*)(double)) &Figure::set_line_width, ref_existing);
    figure_class.def("set_line_colour", (Figure&(Figure::*)(double,double,double)) &Figure::set_line_colour, ref_existing);
    figure_class.def("set_fill_style", (Figure&(Figure::*)(Bool)) &Figure::set_fill_style, ref_existing);
    figure_class.def("set_fill_colour", (Figure&(Figure::*)(double,double,double)) &Figure::set_fill_colour, ref_existing);
    figure_class.def("set_fill_opacity", (Figure&(Figure::*)(double)) &Figure::set_fill_opacity, ref_existing);
    figure_class.def("draw",(Figure&(Figure::*)(const DrawableInterface&))&Figure::draw, ref_existing);
    figure_class.def("draw",(Figure&(Figure::*)(const RealBox&))&Figure::draw, ref_existing);
    figure_class.def("draw",(Figure&(Figure::*)(const ApproximateBoxType&))&Figure::draw, ref_existing);
    figure_class.def("clear",&Figure::clear, ref_existing);
    figure_class.def("display",&Figure::display);
    figure_class.def("write",(Void(Figure::*)(const char*)const)&Figure::write);
    figure_class.def("write",(Void(Figure::*)(const char*,Nat,Nat)const)&Figure::write);
}


Void graphics_submodule() {
    export_figure();
}

