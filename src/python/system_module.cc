/***************************************************************************
 *            python/system_module.cc
 *
 *  21 October 2005
 *  Copyright  2005-6  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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

#include <boost/python.hpp>

#include "real_typedef.h"

template<typename R> void export_map();
template<typename R> void export_affine_map();
template<typename R> void export_affine_multimap();
template<typename R> void export_polynomial_map();
template<typename R> void export_henon_map();
 
template<typename R> void export_vector_field();
template<typename R> void export_affine_vector_field();
template<typename R> void export_lorenz_system();

BOOST_PYTHON_MODULE(system)
{
  export_map<Ariadne::Real>();
  export_affine_map<Ariadne::Real>();
  export_affine_multimap<Ariadne::Real>();
  export_polynomial_map<Ariadne::Real>();
  export_henon_map<Ariadne::Real>();
  
  export_vector_field<Ariadne::Real>();
  export_affine_vector_field<Ariadne::Real>();
  export_lorenz_system<Ariadne::Real>();
}
