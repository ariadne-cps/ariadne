/***************************************************************************
 *            vector_field_orbiter.h
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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
 
/*! \file vector_field_orbiter.h
 *  \brief Methods for computing orbits of boxes.
 */

#ifndef ARIADNE_VECTOR_FIELD_ORBITER_H
#define ARIADNE_VECTOR_FIELD_ORBITER_H

#include <boost/smart_ptr.hpp>

#include "base/types.h"
#include "base/declarations.h"
#include "geometry/declarations.h"
#include "system/declarations.h"

#include "evaluation/evolution_parameters.h"
#include "evaluation/bounder_interface.h"
#include "evaluation/integrator_interface.h"
#include "evaluation/approximator_interface.h"
#include "evaluation/vector_field_orbiter_interface.h"

namespace Ariadne {
  namespace Evaluation {

      
    /*! \brief A class for computing the evolution of a discrete-time autonomous system.
     *  \ingroup Applicators
     */
    template<class BS>
    class VectorFieldOrbiter 
      : public VectorFieldOrbiterInterface<typename BS::real_type>
    {
      typedef Numeric::Rational T;
      typedef typename BS::real_type R;
     public:
      /*! \brief Construct from evolution parameters and an integrator. */
      VectorFieldOrbiter<BS>(const EvolutionParameters<R>&, 
                             const IntegratorInterface<BS>&, 
                             const ApproximatorInterface<BS>&);

      VectorFieldOrbiter<BS>* clone() const;

      virtual 
      Geometry::GridCellListSet<R>
      upper_evolve(const System::VectorField<R>& vf, const Geometry::GridCell<R>& bx, const Numeric::Rational& t) const;

      virtual 
      Geometry::GridCellListSet<R>
      upper_reach(const System::VectorField<R>& vf, const Geometry::GridCell<R>& bx, const Numeric::Rational& t) const;

      virtual 
      Geometry::Box<R>
      lower_evolve(const System::VectorField<R>& vf, const Geometry::Box<R>& bx, const Numeric::Rational& t) const;

      virtual 
      Geometry::BoxListSet<R>
      lower_reach(const System::VectorField<R>& vf, const Geometry::Box<R>& bx, const Numeric::Rational& t) const;


      /*
      virtual 
      Geometry::Orbit<T,BS,BS>*
      orbit(const System::VectorField<R>& f, const Geometry::Box<R>& r, const Numeric::Rational& t) const;
      */

     private:
      // Convenience wrapper for integrator services
      typedef Geometry::ListSet<BS> BSL;
      typedef Geometry::Box<R> Bx;
      typedef Geometry::Grid<R> G;
      typedef Geometry::GridCellListSet<R> GCLS;
      typedef System::VectorField<R> VF;
      std::pair<T, Bx > flow_bounds(const VF& vf, const Bx& bx) {
        return this->_bounder->flow_bounds(vf,bx); }
      BS integration_step(const VF& vf, const BS& bs, const T& t, const Bx& bb) const {
        return this->_integrator->integration_step(vf,bs,t,bb); }
      BS reachability_step(const VF& vf, const BS& bs, const T& t, const Bx& bb) const {
        return this->_integrator->reachability_step(vf,bs,t,bb); }
      BSL subdivide(const BS& bs) const {
        return this->_approximator->subdivide(bs); }
      BS over_approximation(const Bx& bx) const {
        return this->_approximator->over_approximation(bx); }
      Bx bounding_box(const BS& bs) const {
        return this->_approximator->bounding_box(bs); }
      GCLS outer_approximation(const BS& bs, const G& g) const {
        return this->_approximator->outer_approximation(bs,g); }
     private:
      boost::shared_ptr< BounderInterface<R> > _bounder;
      boost::shared_ptr< IntegratorInterface<BS> > _integrator;
      boost::shared_ptr< ApproximatorInterface<BS> > _approximator;
      boost::shared_ptr< EvolutionParameters<R> > _parameters;
    };

  }

}

#endif /* ARIADNE_VECTOR_FIELD_ORBITER_H */
