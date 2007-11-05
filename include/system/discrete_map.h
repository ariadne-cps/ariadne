/***************************************************************************
 *            discrete_map.h
 *
 *  Copyright  2007  Pieter Collins
 *  pieter.collins@cwi.nl
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
 
/*! \file discrete_map.h
 *  \brief A class under defined by orbits of systems.
 */

#ifndef ARIADNE_DISCRETE_MAP_H
#define ARIADNE_DISCRETE_MAP_H

#include <boost/smart_ptr.hpp>

#include "../base/types.h"
#include "../base/declarations.h"
#include "../geometry/declarations.h"
#include "../system/declarations.h"
#include "../evaluation/declarations.h"

namespace Ariadne {
  namespace System {

    template<class R>
    class DiscreteMapInterface 
    {
     public:
      virtual ~DiscreteMapInterface<R>() { }
      virtual DiscreteMapInterface<R>* clone() const = 0;
      virtual size_type argument_dimension() const = 0;
      virtual size_type result_dimension() const = 0;
      virtual Geometry::Rectangle<R> apply(const Geometry::Rectangle<R>&) const = 0;
      virtual Geometry::GridCellListSet<R> apply(const Geometry::GridCell<R>&, const Geometry::Grid<R>& g) const = 0;
      virtual Geometry::DiscreteTimeOrbit<Numeric::Integer,Geometry::Rectangle<R> > orbit(const Geometry::Rectangle<R>&, const Numeric::Integer&, const R& s) const = 0;
      virtual Geometry::DiscreteTimeOrbit<Numeric::Integer,Geometry::GridCellListSet<R> > orbit(const Geometry::GridCell<R>&, const Numeric::Integer&) const = 0;
    };

    template<class R>
    class DiscreteMap
      : public DiscreteMapInterface<R>
    {
     public:
      DiscreteMap(const System::MapInterface<R>& map, const Evaluation::MapOrbiterInterface<R>& orbiter);
      const MapInterface<R>& map() const;
      virtual DiscreteMap<R>* clone() const;
      virtual size_type argument_dimension() const;
      virtual size_type result_dimension() const;
      virtual Geometry::Rectangle<R> apply(const Geometry::Rectangle<R>&) const;
      virtual Geometry::GridCellListSet<R> apply(const Geometry::GridCell<R>&, const Geometry::Grid<R>& g) const;
      virtual Geometry::DiscreteTimeOrbit<Numeric::Integer,Geometry::Rectangle<R> > orbit(const Geometry::Rectangle<R>&, const Numeric::Integer&, const R& s) const;
      virtual Geometry::DiscreteTimeOrbit<Numeric::Integer,Geometry::GridCellListSet<R> > orbit(const Geometry::GridCell<R>&, const Numeric::Integer&) const;
     private:
      System::MapInterface<R>* _map;
      Evaluation::MapOrbiterInterface<R>* _orbiter;
    };

    template<class R>
    class DiscretizedVectorField
      : public DiscreteMapInterface<R>
    {
     public:
      DiscretizedVectorField(
                             const System::VectorFieldInterface<R>& vf, 
                             const Evaluation::VectorFieldOrbiterInterface<R>& orbiter);
      const MapInterface<R>& map() const;
      virtual DiscretizedVectorField<R>* clone() const;
      virtual size_type argument_dimension() const;
      virtual size_type result_dimension() const;
      virtual Geometry::DiscreteTimeOrbit<Numeric::Integer,Geometry::Rectangle<R> > orbit(const Geometry::Rectangle<R>&, const Numeric::Integer&, const R& s) const;
      virtual Geometry::DiscreteTimeOrbit<Numeric::Integer,Geometry::GridCellListSet<R> > orbit(const Geometry::GridCell<R>&, const Numeric::Integer&) const;
     private:
      System::VectorFieldInterface<R>* _map;
      Evaluation::VectorFieldOrbiterInterface<R>* _orbiter;
    };

  }
}


#include "geometry/orbit.h"
#include "geometry/rectangle.h"
#include "geometry/grid_cell.h"
#include "geometry/grid_cell_list_set.h"
#include "system/map_interface.h"
#include "evaluation/orbiter_interface.h"



namespace Ariadne {

template<class R>
System::DiscreteMap<R>::DiscreteMap(const System::MapInterface<R>& map, const Evaluation::MapOrbiterInterface<R>& orbiter)
  : _map(map.clone()), _orbiter(orbiter.clone()) 
{
}

template<class R>
const System::MapInterface<R>&
System::DiscreteMap<R>::map() const 
{
  return *this->_map;
}

template<class R>
System::DiscreteMap<R>*
System::DiscreteMap<R>::clone() const 
{
  return new DiscreteMap<R>(*this->_map, *this->_orbiter);
}


template<class R>
size_type
System::DiscreteMap<R>::argument_dimension() const 
{
  return this->_map->argument_dimension();
}

template<class R>
size_type
System::DiscreteMap<R>::result_dimension() const 
{
  return this->_map->result_dimension();
}



template<class R>
Geometry::Rectangle<R> 
System::DiscreteMap<R>::apply(const Geometry::Rectangle<R>& r) const  
{
  return this->_orbiter->apply(this->map(),r); 
}

template<class R>
Geometry::GridCellListSet<R> 
System::DiscreteMap<R>::apply(const Geometry::GridCell<R>& gc, const Geometry::Grid<R>& g) const 
{
  return this->_orbiter->apply(this->map(),gc,g); 
}


template<class R>
Geometry::DiscreteTimeOrbit<Numeric::Integer,Geometry::Rectangle<R> > 
System::DiscreteMap<R>::orbit(const Geometry::Rectangle<R>& r, const Numeric::Integer& n, const R& s) const  
{
  return this->_orbiter->orbit(this->map(),r,n,s); 
}

template<class R>
Geometry::DiscreteTimeOrbit<Numeric::Integer,Geometry::GridCellListSet<R> > 
System::DiscreteMap<R>::orbit(const Geometry::GridCell<R>& gc, const Numeric::Integer& n) const 
{
  return this->_orbiter->orbit(this->map(),gc,n); 
}



}


#endif /* ARIADNE_DISCRETE_MAP_H */
