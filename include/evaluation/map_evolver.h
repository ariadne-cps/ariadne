/***************************************************************************
 *            map_evolver.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
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
 
/*! \file map_evolver.h
 *  \brief Methods for computing the images of sets under maps.
 *
 * This class works by approximating sets based on grids. 
 * A System::DiscreteMap class is made using a MapOrbiter which is fed into a ModelChecker which runs the compuations.
 */

#ifndef ARIADNE_MAP_EVOLVER_H
#define ARIADNE_MAP_EVOLVER_H

#include <boost/smart_ptr.hpp>

#include "base/types.h"
#include "base/declarations.h"
#include "geometry/declarations.h"
#include "system/declarations.h"
#include "evaluation/declarations.h"


#include "evaluation/map_evolver_interface.h"

// For approximation
#include "geometry/grid_approximation.h"
#include "geometry/timed_list_set.h"


namespace Ariadne {
  namespace Evaluation {

    /*! \brief A class for computing the evolution of a discrete-time autonomous system.
     *  \ingroup Applicators
     */
    template<class BS>
    class MapEvolver 
      : public MapEvolverInterface<typename BS::real_type>
    {
      typedef typename BS::real_type R;
      typedef Numeric::Integer T;
      typedef Numeric::Interval<R> I;
     private:
      boost::shared_ptr< EvolutionParameters<R> > _parameters;
      boost::shared_ptr< ApplicatorInterface<BS> > _applicator;
      boost::shared_ptr< ApproximatorInterface<BS> > _approximator;
      boost::shared_ptr< SubdividerInterface<BS> > _subdivider;
     public:
      //@{
      //! \name Constructors and destructors

      /*! \brief Default constructor chooses appropriate parameter values for maximum basic set radius and grid size. */
      MapEvolver();
      
      /*! \brief Construct from evolution parameters and a method for iterating basic sets. */
      MapEvolver(const EvolutionParameters<R>& parameters);
      
      /*! \brief Construct from evolution parameters and a method for iterating basic sets. */
      MapEvolver(const EvolutionParameters<R>& parameters, 
                 const ApplicatorInterface<BS>& applicator);
      
      /*! \brief Construct from evolution parameters, a method for iterating basic sets, and a scheme for approximating sets. */
      MapEvolver(const EvolutionParameters<R>& parameters, 
                 const ApplicatorInterface<BS>& applicator, 
                 const ApproximatorInterface<BS>& approximator,
                 const SubdividerInterface<BS>& subdivider);
      
      //@}


      //@{ 
      //! \name Methods to set and get the parameters controlling the accuracy.

      /*! \brief The parameters controlling the accuracy. */
      virtual const EvolutionParameters<R>& parameters() const;

      /*! \brief A reference to the parameters controlling the accuracy. */
      virtual EvolutionParameters<R>& parameters();

      //@}


      //@{
      //! \name Evaluation of maps on abstract sets

    
      /*! \brief Compute an approximation to the set obtained by iterating \a steps times \a map starting in \a initial_set. */
      virtual
      Geometry::SetInterface<R>*
      lower_evolve(const System::Map<R>& map, 
                   const Geometry::SetInterface<R>& initial_set,
                   const Numeric::Interval<Numeric::Integer>& steps) const;
    
      /*! \brief Compute an approximation to the set obtained by iterating \a steps times \a map starting in \a initial_set. */
      virtual
      Geometry::SetInterface<R>*
      upper_evolve(const System::Map<R>& map, 
                   const Geometry::SetInterface<R>& initial_set,
                   const Numeric::Interval<Numeric::Integer>& steps) const;
    


      /*! \brief Compute an approximation to the set obtained by iterating \a steps times \a map starting in \a initial_set. */
      virtual
      Geometry::SetInterface<R>*
      lower_evolve(const System::Map<R>& map, 
                   const Geometry::SetInterface<R>& initial_set,
                   const Numeric::Integer& steps) const;
    
      /*! \brief Compute an approximation to the reachable set of \a map starting in \a initial_set iterating at most \a steps times. */
      virtual
      Geometry::SetInterface<R>*
      lower_reach(const System::Map<R>& map, 
                  const Geometry::SetInterface<R>& initial_set,
                  const Numeric::Integer& steps) const;
    
      /*! \brief Compute an approximation to the set obtained by iterating \a steps times \a map starting in \a initial_set. */
      virtual
      Geometry::SetInterface<R>*
      upper_evolve(const System::Map<R>& map, 
                   const Geometry::SetInterface<R>& initial_set,
                   const Numeric::Integer& steps) const;
    
      /*! \brief Compute an approximation to the reachable set of \a map starting in \a initial_set iterating at most \a steps times. */
      virtual
      Geometry::SetInterface<R>*
      upper_reach(const System::Map<R>& map, 
                  const Geometry::SetInterface<R>& initial_set,
                  const Numeric::Integer& steps) const;
    
      /*! \brief Compute an outer-approximation to the chain-reachable set of \a map starting in \a initial_set. */
      virtual
      Geometry::SetInterface<R>*
      chainreach(const System::Map<R>& map, 
                 const Geometry::SetInterface<R>& initial_set) const;
    
      /*! \brief Compute an outer-approximation to the viability kernel of \a map within \a bounding_set. */
      virtual
      Geometry::SetInterface<R>* 
      viable(const System::Map<R>& map, 
             const Geometry::SetInterface<R>& bounding_set) const;
    
      /*! \brief Attempt to verify that the reachable set of \a map starting in \a initial_set remains in \a safe_set. */
      virtual
      tribool
      verify(const System::Map<R>& map, 
             const Geometry::SetInterface<R>& initial_set, 
             const Geometry::SetInterface<R>& safe_set) const;
      //@}


     private:
      // Simplifying typedefs
      typedef Geometry::SetInterface<R> SI;
      typedef Geometry::ListSet<BS> BSL;
      typedef Geometry::Box<R> Bx;
      typedef Geometry::BoxListSet<R> BxLS;
      typedef Geometry::Grid<R> Gr;
      typedef Geometry::GridBlock<R> GB;
      typedef Geometry::GridCellListSet<R> GCLS;
      typedef Geometry::GridMaskSet<R> GMS;
      typedef System::Map<R> Mp;
      typedef Geometry::TimedSet<T,BS> TBS;
      typedef Geometry::ListSet<TBS> TBSL;
     private:
      // Services provided by other classes
      BS apply(const Mp& f, const BS& bs) const {
        return this->_applicator->apply(f,bs); }
      R radius(const BS& bs) const {
        return this->_approximator->radius(bs); }
      BS basic_set(const Bx& bx) const {
        return this->_approximator->basic_set(bx); }
      BSL subdivide(const BS& bs) const {
        return this->_subdivider->subdivide(bs,this->maximum_basic_set_radius()); }
      Bx bounding_box(const BS& bs) const {
        return this->_approximator->bounding_box(bs); }
      Bx bounding_box(const BSL& bsl) const {
        return this->_approximator->bounding_box(bsl); }
      GCLS outer_approximation(const BS& bs, const Gr& g) const {
        return this->_approximator->outer_approximation(bs,g); }
      GCLS outer_approximation(const BSL& bsl, const Gr& g) const {
        return this->_approximator->outer_approximation(bsl,g); }
      GCLS lower_approximation(const BSL& bsl, const Gr& g) const {
        return this->_approximator->outer_approximation(bsl,g); }
     private:
      // Helper functions for timed sets
      BSL basic_set_list(const GCLS& gcls) const {
        BSL result; 
        for(size_type i=0; i!=gcls.size(); ++i) { 
          result.adjoin(this->basic_set(gcls[i])); }
        return result; }
      BSL basic_set_list(const BxLS& bxls) const {
        BSL result; 
        for(size_type i=0; i!=bxls.size(); ++i) { 
          result.adjoin(this->basic_set(bxls[i])); }
        return result; }
      TBSL timed_basic_set_list(const GCLS& gcls) const {
        TBSL result; 
        for(size_type i=0; i!=gcls.size(); ++i) { 
          result.adjoin(TBS(0,this->basic_set(gcls[i]))); }
        return result; }
      TBSL timed_basic_set_list(const BxLS& bxls) const {
        TBSL result; 
        for(size_type i=0; i!=bxls.size(); ++i) { 
          result.adjoin(TBS(0,this->basic_set(bxls[i]))); }
        return result; }
      R radius(const TBS& tbs) const {
        return this->_approximator->radius(tbs.set()); }
      TBSL subdivide(const TBS& tbs) const {
        BSL sets=this->subdivide(tbs.set());
        TBSL result; 
        for(size_type i=0; i!=sets.size(); ++i) { 
          result.adjoin(TBS(tbs.time(),sets[i])); } 
        return result; }
      TBS apply(const Mp& map, const TBS& tbs) const {
        return TBS(tbs.time()+1,this->apply(map,tbs.set())); }
     private:
      // Helper functions for computing orbits sets
      void _step(BSL& reach, BSL& evolve, TBSL& working, const Mp& map, const T& time, Semantics semantics) const;
     public:
      GCLS _upper_reach(const Mp& map, const GCLS& initial_set, const T& time) const;
      GCLS _upper_evolve(const Mp& map, const GCLS& initial_set, const T& time) const;
     private:
      // Helper functions for approximating sets
      Gr grid(dimension_type d) const { 
        return this->_parameters->grid(d); }
      GCLS outer_approximation(const SI& s) const {
        return Geometry::outer_approximation(s,this->grid(s.dimension())); }
      BxLS lower_approximation(const SI& s) const {
        return Geometry::lower_approximation(s,this->grid(s.dimension())); }
      GCLS inner_approximation(const SI& s) const {
        return Geometry::inner_approximation(s,this->grid(s.dimension())); }
     private:
      // Helper functions for accessing parameters
      T lock_to_grid_steps() const { return this->_parameters->lock_to_grid_steps(); }
      R maximum_basic_set_radius() const { return this->_parameters->maximum_basic_set_radius(); }
      Bx bounding_domain(const Mp& map) const { return this->_parameters->bounding_domain(map.dimension()); }
     private:
      // Helper functions for output parameters
    };


  }
}



#endif /* ARIADNE_MAP_EVOLVER_H */
