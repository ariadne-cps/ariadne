/***************************************************************************
 *            applicator.h
 *
 *  17 January 2006
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
 
/*! \file applicator.h
 *  \brief Methods for computing the images of sets under maps.
 */

#ifndef ARIADNE_APPLICATOR_H
#define ARIADNE_APPLICATOR_H

#include <boost/smart_ptr.hpp>

#include "../base/types.h"
#include "../base/declarations.h"
#include "../geometry/declarations.h"
#include "../system/declarations.h"

namespace Ariadne {
  namespace Evaluation {

    /*! \brief A class for computing the image of a set under a map. 
     *  \ingroup Applicators
     */
    template<class R>
    class Applicator {
     private:
      R _maximum_basic_set_radius;
     public:
      /*! \brief Default constructor. */
      Applicator();
      
      /*! \brief Compute the image of a basic set under a continuous function. */
      virtual ~Applicator();
      
      /*! \brief Set the maximum allowable radius of a basic set during iteration. */
      void set_maximum_basic_set_radius(const R&);

      /*! \brief The maximum allowable radius of a basic set during iteration. */
      virtual R maximum_basic_set_radius() const;

      /*! \brief Compute the image of a rectangle under a continuous function. */
      virtual Geometry::Rectangle<R> image(const System::Map<R>& f, const Geometry::Rectangle<R>& s) const;

      /*! \brief Compute the image of a zonotope under a differentiable function. */
      virtual Geometry::Zonotope<R> image(const System::Map<R>& f, const Geometry::Zonotope<R>& s) const;

      /*! \brief Compute the image of an interval zonotope under a differentiable function. */
      virtual Geometry::Zonotope< Numeric::Interval<R> > 
      image(const System::Map<R>& f, const Geometry::Zonotope< Numeric::Interval<R> >& s) const;
     
     protected:
      /*! \brief Template for computing the image of a list set. */
      template<class BS>
      Geometry::ListSet<BS> 
      image_list_set(const System::Map<R>& f, 
                     const Geometry::ListSet<BS>& initial_set) const;

      
      /*! \brief Template for computing the image of a basic set. */
      template<class BS>
      BS
      image_basic_set(const System::Map<R>& f, 
                      const BS& initial_set) const;
     public:
      /*! \brief Compute the image of a list set under a map. */
      virtual 
      Geometry::ListSet< Geometry::Rectangle<R> > 
      image(const System::Map<R>& f, const Geometry::ListSet< Geometry::Rectangle<R> >& ds) const;
       
      /*! \brief Compute the image of a list set under a map. */
      virtual 
      Geometry::ListSet< Geometry::Zonotope<R> > 
      image(const System::Map<R>& f, const Geometry::ListSet< Geometry::Zonotope<R> >& ds) const;
      
      
      /*! \brief Compute the image of a list set under a map. */
      virtual 
      Geometry::ListSet< Geometry::Zonotope<Numeric::Interval<R> > >
      image(const System::Map<R>& f, const Geometry::ListSet< Geometry::Zonotope<Numeric::Interval<R> > >& ds) const;
      
      
      /*! \brief Compute the image of \a map starting in \a initial_set computing the result on \a grid. */
      virtual
      Geometry::GridCellListSet<R> 
      image(const System::Map<R>& map, 
            const Geometry::GridCellListSet<R>& initial_set,
            const Geometry::Grid<R>& grid) const;

      /*! \brief Compute the image of \a map starting in \a initial_set computing the result on \a grid. */
      virtual
      Geometry::GridMaskSet<R> 
      image(const System::Map<R>& map, 
            const Geometry::GridMaskSet<R>& initial_set,
            const Geometry::FiniteGrid<R>& grid) const;

      /*! \brief Compute the image of \a map starting in \a initial_set while remaining in \a bounding_set. */
      virtual
      Geometry::GridMaskSet<R> 
      image(const System::Map<R>& map, 
            const Geometry::GridMaskSet<R>& initial_set,
            const Geometry::GridMaskSet<R>& bounding_set) const;


      /*! \brief Compute the reachable set of \a map starting in \a initial_set. */
      virtual
      Geometry::GridMaskSet<R> 
      reach(const System::Map<R>& map, 
            const Geometry::GridMaskSet<R>& initial_set) const;

      
      /*! \brief Compute the chain-reachable set of \a map starting in \a initial_set while staying within \a bounding_set. */
      virtual
      Geometry::GridMaskSet<R> 
      chainreach(const System::Map<R>& map, 
                 const Geometry::GridMaskSet<R>& initial_set, 
                 const Geometry::GridMaskSet<R>& bounding_set) const;
    
      /*! \brief Attempt to verify that the reachable set of \a map starting in \a initial_set remains in \a safe_set. */
      virtual
      tribool
      verify(const System::Map<R>& map, 
             const Geometry::GridMaskSet<R>& initial_set, 
             const Geometry::GridMaskSet<R>& safe_set) const;

      /*! \brief Discretize a system on a grid. */ 
      virtual 
      System::GridMultiMap<R> 
      discretize(const System::Map<R>& f, 
                 const Geometry::GridMaskSet<R>& dom,
                 const Geometry::Grid<R>& range_grid) const;
     
      /*! \brief Compute a controller for a control-to-target problem. */ 
      virtual 
      System::GridMultiMap<R> 
      control_synthesis(const System::DiscreteTimeSystem<R>& f, 
                        const Geometry::SetInterface<R>& initial_set,
                        const Geometry::SetInterface<R>& target_set,
                        const Geometry::GridMaskSet<R>& state_bounding_set,
                        const Geometry::GridMaskSet<R>& control_set,
                        const Geometry::GridMaskSet<R>& noise_set) const;
     


       
    };



  }
}

#endif /* ARIADNE_APPLY_H */
