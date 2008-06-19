/***************************************************************************
 *            build_map.h
 *
 *  Copyright 2007  Pieter Collins
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
 
/*! \file build_map.h
 *  \brief Macro to build a map from a raw function template.
 */

#ifndef ARIADNE_BUILD_MAP_H
#define ARIADNE_BUILD_MAP_H

#include "geometry/point.h"

#include "function/build_function.h"
#include "function/function_interface.h"
#include "system/map.h"
 


#define ARIADNE_BUILD_MAP(Nm,f,rd,ad,np,sm)   \
  ARIADNE_BUILD_FUNCTION(Nm##Function,f,rd,ad,np,sm) \
   \
  template<class R>                        \
  class Nm##Map \
    : public Map<R> \
  { \
   public: \
    template<class P> explicit Nm##Map(const P& p) \
      : Map<R>(Nm##Function<R>(p.position_vector())), _parameters(p) { } \
    const Point<R>& parameters() const { return this->_parameters; } \
    const R& parameter(size_type k) const { return this->_parameters[k]; }    \
   private: \
    Point<R> _parameters;                  \
  }; \


#endif /* ARIADNE_BUILD_MAP_H */
