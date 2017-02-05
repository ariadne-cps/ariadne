/***************************************************************************
 *            ariadne_module.cpp
 *
 *  Copyright  2007-8  Pieter Collins
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

void numeric_submodule();
void linear_algebra_submodule();
void optimization_submodule();
void differentiation_submodule();
void function_submodule();
void calculus_submodule();
void geometry_submodule();
void solver_submodule();
void storage_submodule();
void system_submodule();
void evolution_submodule();
void graphics_submodule();

BOOST_PYTHON_MODULE(ariadne)
{
    numeric_submodule();
    linear_algebra_submodule();
    optimization_submodule();
    differentiation_submodule();
    function_submodule();
    calculus_submodule();
    geometry_submodule();
    solver_submodule();
    storage_submodule();
    system_submodule();
    evolution_submodule();
    graphics_submodule();
}
