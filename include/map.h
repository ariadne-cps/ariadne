/***************************************************************************
 *            map.h
 *
 *  Copyright  2004-8  Alberto Casagrande, Pieter Collins
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
 
/*! \file map.h
 *  \brief Main continuous dynamics system class.
 */

#ifndef ARIADNE_MAP_H
#define ARIADNE_MAP_H

#include <boost/shared_ptr.hpp>

#include "function.h"
#include "set_interface.h"
#include "grid.h"

namespace Ariadne {  

class VectorFunction;

/*! \brief An iterated function system in Euclidean space.
 */
class IteratedMap
{
  public:
    //! \brief The type used to represent time. 
    typedef Float TimeType;
    //! \brief The type used to represent real numbers. 
    typedef Float RealType ;
    //! \brief The type used to describe the state space. 
    typedef EuclideanSpace StateSpaceType;
  public:
    IteratedMap(const VectorFunction& f) : _function(f) { }
    const VectorFunction& function() const { return _function; }
    Grid grid() const { return Grid(_function.argument_size()); }
  private:
    VectorFunction _function;
};

inline std::ostream& operator<<(std::ostream& os, const IteratedMap& vf) {
    return os << "IteratedMap( " << vf.function() << " )";
}


} // namespace Ariadne

#endif // ARIADNE_MAP_H 
