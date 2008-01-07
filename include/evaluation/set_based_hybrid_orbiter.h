/***************************************************************************
 *            set_based_hybrid_orbiter.h
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
 
/*! \file set_based_hybrid_orbiter.h
 *  \brief Methods for computing orbits of boxes under a hybrid system.
 */

#ifndef ARIADNE_SET_BASED_HYBRID_ORBITER_H
#define ARIADNE_SET_BASED_HYBRID_ORBITER_H

#include <boost/smart_ptr.hpp>

#include "base/types.h"
#include "base/declarations.h"
#include "geometry/declarations.h"
#include "system/declarations.h"
#include "evaluation/declarations.h"

#include "evaluation/evolution_parameters.h"
#include "evaluation/bounder_interface.h"
#include "evaluation/integrator_interface.h"
#include "evaluation/applicator_interface.h"
#include "evaluation/approximator_interface.h"
#include "evaluation/satisfier_interface.h"
#include "evaluation/set_based_hybrid_orbiter_interface.h"

namespace Ariadne {
  namespace Evaluation {

      
    /*! \brief A class for computing the evolution of a hybrid-time system.
     *  \ingroup Evaluation
     */
    template<class BS>
    class SetBasedHybridOrbiter 
      : public SetBasedHybridOrbiterInterface<typename BS::real_type>
    {
      typedef Numeric::Rational Q;
      typedef Numeric::Integer Z;
      typedef typename BS::real_type R;
     public:
      /*! \brief Construct from evolution parameters and an integrator. */
      SetBasedHybridOrbiter<BS>(const EvolutionParameters<R>&, 
                                const IntegratorInterface<BS>&, 
                                const ApplicatorInterface<BS>&, 
                                const ApproximatorInterface<BS>&, 
                                const SatisfierInterface<BS>&);

      SetBasedHybridOrbiter<BS>* clone() const;

      virtual 
      Geometry::HybridListSet<BS>
      evolution(const System::SetBasedHybridAutomaton<R>& vf, 
                const Geometry::HybridListSet<BS>& bx, 
                const Numeric::Rational& t,
                const Semantics semantics,
                const EvolutionType evolution_type) const;

      virtual 
      Geometry::HybridGridCellListSet<R>
      upper_evolve(const System::SetBasedHybridAutomaton<R>& vf, 
                   const Geometry::HybridGridCell<R>& bx, 
                   const Numeric::Rational& t) const;

      virtual 
      Geometry::HybridGridCellListSet<R>
      upper_reach(const System::SetBasedHybridAutomaton<R>& vf, 
                  const Geometry::HybridGridCell<R>& bx, 
                  const Numeric::Rational& t) const;

      virtual 
      Geometry::HybridBox<R>
      lower_evolve(const System::SetBasedHybridAutomaton<R>& vf, 
                   const Geometry::HybridBox<R>& bx, 
                   const Numeric::Rational& t) const;

      virtual 
      Geometry::HybridBoxListSet<R>
      lower_reach(const System::SetBasedHybridAutomaton<R>& vf, 
                  const Geometry::HybridBox<R>& bx, 
                  const Numeric::Rational& t) const;


     private:     
      // Functions which should be provided elsewhere
      Geometry::HybridGrid<R> grid(const Geometry::HybridSpace& loc) const;
      Geometry::HybridGridCellListSet<R> outer_approximation(const Geometry::HybridListSet<BS>& hls, const Geometry::HybridGrid<R>& hg) const;
     private:
      // Convenience wrapper for services
      typedef Geometry::ListSet<BS> BSL;
      typedef Geometry::Box<R> Bx;
      typedef System::Map<R> Mp;
      typedef System::VectorField<R> VF;
      typedef System::SetBasedDiscreteMode<R> DM;
      typedef System::SetBasedDiscreteTransition<R> DT;
      typedef Geometry::ConstraintSet<R> CS;
      typedef Geometry::Grid<R> G;
      typedef Geometry::GridCell<R> GC;
      typedef Geometry::GridCellListSet<R> GCLS;
      typedef Geometry::HybridBox<R> HBx;
      std::pair<Q, Bx > flow_bounds(const VF& vf, const Bx& bx) const {
        return this->_bounder->flow_bounds(vf,bx); }
      BS integration_step(const VF& vf, const BS& bs, const Q& t, const Bx& bb) const {
        return this->_integrator->integration_step(vf,bs,t,bb); }
      BS reachability_step(const VF& vf, const BS& bs, const Q& t, const Bx& bb) const {
        return this->_integrator->reachability_step(vf,bs,t,bb); }
      BS apply(const Mp& f, const BS& bs) const {
        return this->_applicator->apply(f,bs); }
      BSL subdivide(const BS& bs) const {
        return this->_approximator->subdivide(bs); }
      BS over_approximation(const Bx& bx) const {
        return this->_approximator->over_approximation(bx); }
      Bx bounding_box(const BS& bs) const {
        return this->_approximator->bounding_box(bs); }
      R radius(const BS& bs) const {
        return bs.radius(); }
      GCLS outer_approximation(const BS& bs, const G& g) const {
        return this->_approximator->outer_approximation(bs,g); }
      tribool subset(const BS& bs, const CS& cs) const {
        return this->_satisfier->subset(bs,cs); }
      tribool intersects(const BS& bs, const CS& cs) const {
        return this->_satisfier->intersects(bs,cs); }
      tribool subset(const BS& bs, const Geometry::SetInterface<R>& set) const {
        return set.superset(this->bounding_box(bs)); }
      tribool intersects(const BS& bs, const Geometry::SetInterface<R>& set) const {
        return set.intersects(this->bounding_box(bs)); }
     private:
      boost::shared_ptr< BounderInterface<R> > _bounder;
      boost::shared_ptr< IntegratorInterface<BS> > _integrator;
      boost::shared_ptr< ApplicatorInterface<BS> > _applicator;
      boost::shared_ptr< ApproximatorInterface<BS> > _approximator;
      boost::shared_ptr< SatisfierInterface<BS> > _satisfier;
      boost::shared_ptr< EvolutionParameters<R> > _parameters;
    };

  }

}

#endif /* ARIADNE_SET_BASED_HYBRID_ORBITER_H */
