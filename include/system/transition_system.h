/***************************************************************************
 *            transition_system.h
 *
 *  Copyright  2007  Pieter Collins
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
 
/*! \file transition_system.h
 *  \brief A class under defined by orbits of systems.
 */

#ifndef ARIADNE_TRANSITION_SYSTEM_H
#define ARIADNE_TRANSITION_SYSTEM_H

#include <boost/smart_ptr.hpp>

#include "base/types.h"
#include "base/declarations.h"
#include "geometry/declarations.h"
#include "system/declarations.h"
#include "evaluation/declarations.h"

#include "system/transition_system_interface.h"

namespace Ariadne {
  namespace System {


    template<class R>
    class TransitionSystem
      : public TransitionSystemInterface<R>
    {
     public:
      TransitionSystem(const System::Map<R>& map, const Evaluation::MapEvolverInterface<R>& evolver, const Numeric::Integer& steps);
      const Map<R>& map() const;
      virtual TransitionSystem<R>* clone() const;
      virtual size_type dimension() const;
      virtual std::pair< Numeric::Integer, Geometry::Box<R> > lower_evolve(const Geometry::Box<R>&) const;
      virtual Geometry::BoxListSet<R> lower_reach(const Geometry::Box<R>&) const;
      virtual Geometry::GridCellListSet<R> upper_evolve(const Geometry::GridCell<R>&) const;
      virtual Geometry::GridCellListSet<R> upper_reach(const Geometry::GridCell<R>&) const;
     private:
      System::Map<R>* _map;
      Evaluation::MapOrbiterInterface<R>* _orbiter;
      Numeric::Integer _steps;
    };

    template<class R>
    class DiscretizedVectorField
      : public TransitionSystemInterface<R>
    {
     public:
      DiscretizedVectorField(const System::VectorField<R>& vf, 
                             const Evaluation::VectorFieldOrbiterInterface<R>& orbiter);
      const Map<R>& map() const;
      virtual DiscretizedVectorField<R>* clone() const;
      virtual size_type argument_dimension() const;
      virtual size_type result_dimension() const;
     private:
      System::VectorField<R>* _vf;
      Evaluation::VectorFieldOrbiterInterface<R>* _orbiter;
    };

  }
}


#include "geometry/box_list_set.h"
#include "geometry/orbit.h"
#include "geometry/rectangle.h"
#include "geometry/grid_cell.h"
#include "geometry/grid_cell_list_set.h"
#include "system/map.h"
#include "evaluation/map_orbiter_interface.h"



namespace Ariadne {

template<class R>
System::TransitionSystem<R>::TransitionSystem(const System::Map<R>& map, 
                                              const Evaluation::MapOrbiterInterface<R>& orbiter,
                                              const Numeric::Integer& steps)
  : _map(map.clone()), _orbiter(orbiter.clone()) , _steps(steps)
{
}

template<class R>
const System::Map<R>&
System::TransitionSystem<R>::map() const 
{
  return *this->_map;
}

template<class R>
System::TransitionSystem<R>*
System::TransitionSystem<R>::clone() const 
{
  return new TransitionSystem<R>(*this);
}


template<class R>
size_type
System::TransitionSystem<R>::dimension() const 
{
  return this->_map->argument_dimension();
}



template<class R>
std::pair< Numeric::Integer, Geometry::Box<R> >
System::TransitionSystem<R>::lower_evolve(const Geometry::Box<R>& bx) const 
{
  return this->_orbiter->lower_evolve(this->map(),bx,this->_steps); 
}

template<class R>
Geometry::BoxListSet<R>
System::TransitionSystem<R>::lower_reach(const Geometry::Box<R>& bx) const 
{
  return this->_orbiter->lower_reach(this->map(),bx,this->_steps); 
}

template<class R>
Geometry::GridCellListSet<R> 
System::TransitionSystem<R>::upper_evolve(const Geometry::GridCell<R>& gc) const 
{
  return this->_orbiter->upper_evolve(this->map(),gc,this->_steps); 
}

template<class R>
Geometry::GridCellListSet<R> 
System::TransitionSystem<R>::upper_reach(const Geometry::GridCell<R>& gc) const 
{
  return this->_orbiter->upper_reach(this->map(),gc,this->_steps); 
}


}


#endif /* ARIADNE_TRANSITION_SYSTEM_H */
