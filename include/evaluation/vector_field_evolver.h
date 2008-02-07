/***************************************************************************
 *            vector_field_evolver.h
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
 
/*! \file vector_field_evolver.h
 *  \brief Methods for integrating points and sets under a vector field.
 */

#ifndef ARIADNE_VECTOR_FIELD_EVOLVER_H
#define ARIADNE_VECTOR_FIELD_EVOLVER_H

#include "base/types.h"
#include "base/declarations.h"
#include "numeric/declarations.h"
#include "linear_algebra/declarations.h"
#include "geometry/declarations.h"
#include "system/declarations.h"
#include "evaluation/declarations.h"

#include "geometry/grid_approximation.h"
#include "geometry/timed_list_set.h"
#include "evaluation/vector_field_evolver_interface.h"
#include "evaluation/subdivider_interface.h"

namespace Ariadne {
  namespace Evaluation {

    /*! \brief Class for computing the evolution of a continuous-time autonomous system.
     *  \ingroup Integrate
     */
    template<class BS>
    class VectorFieldEvolver 
      : public VectorFieldEvolverInterface<typename BS::real_type>
    {
      typedef Numeric::Rational T;
      typedef typename BS::real_type R;
     private:
      boost::shared_ptr< EvolutionParameters<R> >  _parameters;
      boost::shared_ptr< IntegratorInterface<BS> >  _integrator;
      boost::shared_ptr< ApproximatorInterface<BS> >  _approximator;
      boost::shared_ptr< SubdividerInterface<BS> >  _subdivider;
      
      boost::shared_ptr< EvolutionProfiler >  _profiler;
     public:
      typedef R real_type;

      //@{
      //! \name Constructors and destructors

      /*! \brief Virtual destructor. */
      virtual ~VectorFieldEvolver();

      /*! \brief Construct from evolution parameters and an integration method, using the standard approximation scheme. */
      VectorFieldEvolver(const EvolutionParameters<R>& parameters, 
			 const IntegratorInterface<BS>& integrator);

      /*! \brief Construct from evolution parameters, an integration method, and an approximation scheme. */
      VectorFieldEvolver(const EvolutionParameters<R>& parameters, 
			 const IntegratorInterface<BS>& integrator,
			 const ApproximatorInterface<BS>& approximator);

      /*! \brief Cloning operator. */
      virtual VectorFieldEvolver<BS>* clone() const;

      //@}

      //@{
      //! \name Parameters controlling the accuracy

      /*! \brief The parameters controlling the accuracy. */
      const EvolutionParameters<R>& parameters() const;
      /*! \brief A reference to the parameters controlling the accuracy. */
      EvolutionParameters<R>& parameters();
      //@}

      //@{
      //! \name Integration routines

      /*! \brief Integrate \a intial_set for time \a time under \a vector_field using lower semantics. Returns a dynamically allocated set. */
      virtual Geometry::SetInterface<R>* lower_evolve(const System::VectorField<R>& vector_field,
                                                      const Geometry::SetInterface<R>& initial_set,
                                                      const Numeric::Rational& time) const;

      /*! \brief Integrate \a intial_set for times up to \a time under \a vector_field. Returns a dynamically allocated set. */
      virtual Geometry::SetInterface<R>* lower_reach(const System::VectorField<R>& vector_field,
                                                     const Geometry::SetInterface<R>& initial_set,
                                                     const Numeric::Rational& time) const;

      /*! \brief Integrate \a intial_set for time \a time under \a vector_field. Returns a dynamically allocated set. */
      virtual Geometry::SetInterface<R>* upper_evolve(const System::VectorField<R>& vector_field,
                                                      const Geometry::SetInterface<R>& initial_set,
                                                      const Numeric::Rational& time) const;

      /*! \brief Integrate \a intial_set for times up to \a time under \a vector_field. Returns a dynamically allocated set. */
      virtual Geometry::SetInterface<R>* upper_reach(const System::VectorField<R>& vector_field,
                                                     const Geometry::SetInterface<R>& initial_set,
                                                     const Numeric::Rational& time) const;


      /*! \brief Compute a lower-approximation to the infinite-time reachable set. */
      virtual 
      Geometry::SetInterface<R>*
      lower_reach(const System::VectorField<R>& vector_field,
                  const Geometry::SetInterface<R>& initial_set) const;

      /*! \brief Integrate \a intial_set for all times under \a vector_field, while remaining in \a bounding_box. 
       *
       * Implemented by repeated calls to integrate(...) followed by a single call to reach(...).
       */
      virtual Geometry::SetInterface<R>* chainreach(const System::VectorField<R>& vector_field,
                                                    const Geometry::SetInterface<R>& initial_set) const;

      /*! \brief Computes the set of points which remain in \a bounding_set under evolution of \a vector_field.
       */
      virtual Geometry::SetInterface<R>* viable(const System::VectorField<R>& vector_field,
                                                const Geometry::SetInterface<R>& bounding_set) const;

      /*! \brief  Verifies that the flow of \a vector_field starting in \a initial_set remains in \a safe_set all times.
       */
      virtual tribool verify(const System::VectorField<R>& vector_field,
                             const Geometry::SetInterface<R>& initial_set,
                             const Geometry::SetInterface<R>& safe_set) const;

      //@}

     private:
      //  Generate a grid
      Geometry::Grid<R> _grid(const dimension_type& d) const;

      //  Integrate a box.
      Geometry::GridCellListSet<R>
      _upper_evolve(const System::VectorField<R>& vector_field, 
                    const Geometry::GridCell<R>& initial_set, 
                    const Numeric::Rational& time) const;

      //  Compute the reachable set from a basic set.
      Geometry::GridCellListSet<R>
      _upper_reach(const System::VectorField<R>& vector_field, 
                   const Geometry::GridCell<R>& initial_set, 
                   const Numeric::Rational& time) const;

     private:
      // Simplifying typedefs
      typedef Geometry::SetInterface<R> SI;
      typedef Geometry::ListSet<BS> BSL;
      typedef Geometry::Box<R> Bx;
      typedef Geometry::BoxListSet<R> BxLS;
      typedef Geometry::Grid<R> Gr;
      typedef Geometry::GridCellListSet<R> GCLS;
      typedef Geometry::GridMaskSet<R> GMS;
      typedef System::VectorField<R> VF;
      typedef Geometry::TimedSet<T,BS> TBS;
      typedef Geometry::ListSet<TBS> TBSL;
     private:
      // Services provided by other classes
      std::pair<T,Bx> flow_bounds(const VF& vf, const Bx& bx) const {
        return this->_integrator->flow_bounds(vf,bx,this->maximum_step_size()); }
      std::pair<T,Bx> flow_bounds(const VF& vf, const Bx& bx, const T& h) const {
        return this->_integrator->flow_bounds(vf,bx,h); }
      TBS integration_step(const VF& vf, const TBS& tbs, const T& h, const Bx& bb) const {
        return TBS(tbs.time()+h,this->_integrator->integration_step(vf,tbs.set(),h,bb)); }
      BS reachability_step(const VF& vf, const BS& bs, const T& h, const Bx& bb) const {
        return this->_integrator->reachability_step(vf,bs,h,bb); }
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
     private:
      // Helper functions for computing orbits sets
      void _step(BSL& reach, BSL& evolve, TBSL& working, const VF& vf, const T& time, Semantics semantics) const;
      GCLS _upper_reach(const VF& vf, const GCLS& initial_set, const T& time) const;
      GCLS _upper_evolve(const VF& vf, const GCLS& initial_set, const T& time) const;
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
      T lock_to_grid_time() const { return this->_parameters->lock_to_grid_steps(); }
      T maximum_step_size() const { return this->_parameters->maximum_step_size(); }
      R maximum_basic_set_radius() const { return this->_parameters->maximum_basic_set_radius(); }
      Bx bounding_domain(const VF& vf) const { return this->_parameters->bounding_domain(vf.dimension()); }
    };

 

    
  }
}

#include "vector_field_evolver.inline.h"

#endif /* ARIADNE_VECTOR_FIELD_EVOLVER_H */
