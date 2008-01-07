/***************************************************************************
 *            name.h
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
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
 
#ifndef ARIADNE_GEOMETRY_NAME_H
#define ARIADNE_GEOMETRY_NAME_H

/*! \file name.h
 *  \brief Names of basic sets.
 */

#include <iostream>

#include "numeric/declarations.h"
#include "geometry/declarations.h"

namespace Ariadne {

  namespace Geometry {

    template<class BS> struct Name;
  
    template<class R> struct Name< Point<R> > {
      static std::string string() { return "Point"; }
    };

    template<class R> struct Name< Box<R> > {
      static std::string string() { return "Box"; }
    };

    template<class R> struct Name< Rectangle<R> > {
      static std::string string() { return "Rectangle"; }
    };

    template<class R> struct Name< Zonotope<R> > {
      static std::string string() { return "Zonotope"; }
    };
  
    template<class R> struct Name< Polytope<R> > {
      static std::string string() { return "Polytope"; }
    };
  
    template<class R> struct Name< Polyhedron<R> > {
      static std::string string() { return "Polyhedron"; }
    };
  

    template<class T, class BS> struct Name< TimedSet<T,BS> > {
      static std::string string() { return std::string("Timed")+Name<BS>().string(); }
    };

    template<class BS> 
    std::string name() {
      return Name<BS>::string();
    }


    }
}

#endif /* ARIADNE_GEOMETRY_NAME_H */
