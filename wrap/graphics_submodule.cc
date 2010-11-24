/***************************************************************************
 *      graphics_submodule.cc
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

#include <boost/python.hpp>

#include "graphics_interface.h"
#include "graphics.h"
#include "point.h"
#include "box.h"
#include "function.h"

using namespace boost::python;


namespace Ariadne {
    template<class S> void draw(FigureInterface& fig, const S& sh) { fig.draw(sh); }
}


using namespace Ariadne;



void export_figure()
{
    class_<FigureInterface,boost::noncopyable>("FigureInterface",no_init);
    class_<Figure, bases<FigureInterface> > figure_class("Figure",init<>());

    //class_<Figure, bases<FigureInterface> > figure_class("Figure",init<>());
    figure_class.def("set_projection_map",(void(Figure::*)(const PlanarProjectionMap&)) &Figure::set_projection_map);
    figure_class.def("set_projection_map",(void(Figure::*)(const ProjectionFunction&)) &Figure::set_projection_map);
    figure_class.def("set_projection",(void(Figure::*)(uint,uint,uint)) &Figure::set_projection);
    figure_class.def("set_bounding_box",&Figure::set_bounding_box);
    figure_class.def("set_line_style", (void(Figure::*)(bool)) &Figure::set_line_style);
    figure_class.def("set_line_width", (void(Figure::*)(double)) &Figure::set_line_width);
    figure_class.def("set_line_colour", (void(Figure::*)(double,double,double)) &Figure::set_line_colour);
    figure_class.def("set_fill_style", (void(Figure::*)(bool)) &Figure::set_fill_style);
    figure_class.def("set_fill_colour", (void(Figure::*)(double,double,double)) &Figure::set_fill_colour);
    figure_class.def("draw",(void(FigureInterface::*)(const DrawableInterface&))&FigureInterface::draw);
    figure_class.def("clear",&Figure::clear);
    figure_class.def("display",&Figure::display);
    figure_class.def("write",&Figure::write);
}


void graphics_submodule() {
    export_figure();
}

