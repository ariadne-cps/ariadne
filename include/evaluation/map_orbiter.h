/***************************************************************************
 *            map_orbiter.h
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
 
/*! \file map_orbiter.h
 *  \brief Methods for computing orbits of boxes.
 */

#ifndef ARIADNE_MAP_ORBITER_H
#define ARIADNE_MAP_ORBITER_H

#include <boost/smart_ptr.hpp>

#include "base/types.h"
#include "base/declarations.h"
#include "geometry/declarations.h"
#include "system/declarations.h"
#include "evaluation/declarations.h"

#include "evaluation/evolution_parameters.h"
#include "evaluation/applicator_interface.h"
#include "evaluation/approximator_interface.h"
#include "evaluation/map_orbiter_interface.h"

namespace Ariadne {
  namespace Evaluation {


    /*! \brief A class for computing the evolution of a discrete-time autonomous system.
     *  \ingroup Applicators
     */
    template<class BS>
    class MapOrbiter 
      : public MapOrbiterInterface<typename BS::real_type>
    {
      typedef Numeric::Integer T;
      typedef typename BS::real_type R;
     public:
      /*! \brief Default construnctor. */
      //MapOrbiter();

      /*! \brief Construct from the evolution parameters and an applicator. */
      MapOrbiter(const EvolutionParameters<R>&, 
                 const ApplicatorInterface<BS>&, 
                 const ApproximatorInterface<BS>&);

      /*! \brief Copy constructor. */
      MapOrbiter(const MapOrbiter<BS>&);

      /*! \brief Make a dynamically-allocated copy. */
      MapOrbiter<BS>* clone() const;

      /*! \brief The maximum radius of a basic set. */
      virtual R maximum_basic_set_radius() const;

      /*! \brief Compute the evolved set under a map. */
      virtual
      Geometry::GridCellListSet<R> 
      upper_evolve(const System::Map<R>& f, const Geometry::GridCell<R>& gc, const Numeric::Integer& n) const;

      /*! \brief Compute the reach set under a map. */
      virtual
      Geometry::GridCellListSet<R> 
      upper_reach(const System::Map<R>& f, const Geometry::GridCell<R>& gc, const Numeric::Integer& n) const;

      /*! \brief Compute the evolved set under a map for at most \a n steps and a maximum radius of \a mr. Returns the number of steps as well as the evolved set. */
      virtual
      std::pair<Numeric::Integer, Geometry::Box<R> >
      lower_evolve(const System::Map<R>& f, const Geometry::Box<R>& s, const Numeric::Integer& n) const;

      /*! \brief Compute the reach set under a map. */
      virtual
      Geometry::BoxListSet<R> 
      lower_reach(const System::Map<R>& f, const Geometry::Box<R>& s, const Numeric::Integer& n) const;

      /*! \brief Compute the orbit of a basic set under \a n steps of continuous function. */
      virtual 
      Geometry::Orbit<Numeric::Integer,BS>*
      orbit(const System::Map<R>& f, const Geometry::Box<R>& r, const Numeric::Integer& n) const;

      /*! \brief Compute the orbit of a basic set under \a n steps of continuous function. */
      Geometry::Orbit<T,BS>
      orbit(const System::Map<R>& f, const BS& bs, const T& n) const;

      virtual std::ostream& write(std::ostream& os) const;

     private:
      // Intermediate computations
      std::vector< Geometry::ListSet<BS> > _orbit(const System::Map<R>& f, const BS& bs, const Numeric::Integer& n, Semantics s) const;
     private:
      // Functions used by the orbiter to compute the result
      typedef Geometry::ListSet<BS> BSL;
      typedef Geometry::Box<R> Bx;
      typedef Geometry::Grid<R> G;
      typedef Geometry::GridCellListSet<R> GCLS;
      typedef System::Map<R> Mp;
      BS apply(const Mp& f, const BS& bs) const {
        return this->_applicator->apply(f,bs); }
      BSL subdivide(const BS& bs) const {
        if(bs.radius()>this->maximum_basic_set_radius()) {
          return this->_approximator->subdivide(bs); }
        else { return bs; } }
      BS over_approximation(const Bx& bx) const {
        return this->_approximator->over_approximation(bx); }
      Bx bounding_box(const BS& bs) const {
        return this->_approximator->bounding_box(bs); }
      GCLS outer_approximation(const BS& bs, const G& g) const {
        return this->_approximator->outer_approximation(bs,g); }
      Bx bounding_box(const BSL& bsl) const {
        return this->_approximator->bounding_box(bsl); }
      GCLS outer_approximation(const BSL& bsl, const G& g) const {
        return this->_approximator->outer_approximation(bsl,g); }
     private:
      R _maximum_basic_set_radius;
      boost::shared_ptr< ApplicatorInterface<BS> > _applicator;
      boost::shared_ptr< ApproximatorInterface<BS> > _approximator;
    };

  }

}

#endif /* ARIADNE_MAP_ORBITER_H */
