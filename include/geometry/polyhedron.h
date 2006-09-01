/***************************************************************************
 *            polyhedron.h
 *
 *  Thu Jan 27 10:26:36 2005
 *  Copyright  2005  Alberto Casagrande
 *  casagrande@dimi.uniud.it
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

/*! \file polyhedron.h
 *  \brief Polyhedra.
 */
 
#ifndef _ARIADNE_POLYHEDRON_H
#define _ARIADNE_POLYHEDRON_H

#include <iosfwd>
#include <ppl.hh>
#include <vector>

#include "../declarations.h"

#include "../linear_algebra/constraint_system.h"
#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"
#include "../linear_algebra/interval_vector.h"

#include "../linear_algebra/generator_system.h"
#include "../geometry/point.h"
#include "../geometry/rectangle.h"

namespace Parma_Polyhedra_Library {
  // Import all the output operators into the main PPL namespace.
  using IO_Operators::operator<<;
}
 
namespace Ariadne {  
  namespace Geometry {

    template<> 
    inline bool is_a<Polyhedron,Polyhedron>() { return true; }

    template<typename R> bool equal(const Polyhedron<R>& A, 
                                       const Polyhedron<R>& B);
    template<typename R> bool disjoint(const Polyhedron<R>& A, 
                                       const Polyhedron<R>& B);
    template<typename R> bool interiors_intersect(const Polyhedron<R>& A, 
                                                  const Polyhedron<R>& B);
    template<typename R> bool inner_subset(const Polyhedron<R>& A, 
                                           const Polyhedron<R>& B);
    template<typename R> bool subset(const Polyhedron<R>& A, 
                                     const Polyhedron<R>& B);
    template<typename R> Polyhedron<R> intersection(const Polyhedron<R>& A, 
                                                    const Polyhedron<R>& B);
    template<typename R> Polyhedron<R> regular_intersection(const Polyhedron<R>& A, 
                                                            const Polyhedron<R>& B);
    template<typename R> Polyhedron<R> convex_hull(const Polyhedron<R>& A, 
                                                   const Polyhedron<R>& B);
    template<typename R> Polyhedron<R> minkowski_sum(const Polyhedron<R>& A,
                                                     const Polyhedron<R>& B);
    template<typename R> Polyhedron<R> minkowski_difference(const Polyhedron<R>& A,
                                                            const Polyhedron<R>& B);
    template<typename R> std::ostream& operator<<(std::ostream& os, 
                                                  const Polyhedron<R>& P);
                                                   
                                                   
    /* Transforms an Point into a Parma_Polyhedra_Library::C_Polyhedron */
    template <typename R>
    Parma_Polyhedra_Library::NNC_Polyhedron 
    _from_Point_to_PPL_Polyhedron(const Point<R>& s); 
    
    /* Transforms an Rectangle into a Parma_Polyhedra_Library::C_Polyhedron */
    template <typename R>
    Parma_Polyhedra_Library::NNC_Polyhedron 
    _from_Rectangle_to_closed_PPL_Polyhedron(const Rectangle<R>& r);

    /*! \ingroup BasicSet
     *  \brief The polyhedron class.
     *
     * This is a wrapper from Parma Polyhedra Library Polyhedron class
     * to Ariadne representation.
     */ 
    template <typename R>
    class Polyhedron {
     public:
      /*! \brief The type of denotable real numbers used to describe the polyhedron. */
      typedef R real_type;
      /*! \brief The type of denotable point contained by the polyhedron. */
      typedef Point<R> state_type;
      /*! \brief The type of vector. */
      typedef Ariadne::LinearAlgebra::Vector<R> vector_type;
      /*! \brief The type of matrix. */
      typedef Ariadne::LinearAlgebra::Matrix<R> matrix_type;
      /*! \brief The type of a list of points. */
      typedef std::vector< Point<R> > state_list_type;
     public:
      friend class Parallelotope<R>;
    
      //@{ 
      //! \name Constructors
      /*! \brief Default constructor creates empty set. 
       */
      Polyhedron();
      
      /*! \brief Construct an empty polyhedron with dimension \a dim.
       *
       * Builds a polyhedron with dimension \a dim.
       * \param dim is the new polyhedron's number of dimensions.
       */
      Polyhedron(size_type dim);
        
      /*! \brief Construct a polyhedron from a constraint system.
       *
       * \param cs is the constraint system.
       */
      Polyhedron(LinearAlgebra::ConstraintSystem<R>& cs);
      
      /*! \brief Construct a polyhedron from a generator system.
       *
       * \param gen is the generator system.
       */
      Polyhedron(LinearAlgebra::GeneratorSystem<R>& gen);
      
      /*! \brief Construct the polyhedron defined by the matrix equations \f$Ax\leq b\f$.
       */
      Polyhedron(const LinearAlgebra::Matrix<R>& A, const LinearAlgebra::Vector<R>& b);
      
      /*! \brief Construct the polyhedron defined by the rectilinear generating system \f$c+Ae\f$.
       */
      Polyhedron(const Point<R>& c, const LinearAlgebra::Matrix<R>& A);
      
      /*! \brief Construct the polyhedron defined as the convex hull of a list of points.
       */
      Polyhedron(const std::vector< Point<R> >& points);
      
      /*! \brief Construct a polyhedron from an interval vector. */
      Polyhedron(const LinearAlgebra::IntervalVector<R>& iv);
            
      /*! \brief Construct a polyhedron from a rectangle. */
      Polyhedron(const Rectangle<R>& rect);
            
      /*! \brief Copy constructor. 
       */
      Polyhedron(const Polyhedron<R>& original);
          
      /*! \brief Copy assignment operator. 
       *
       * \param original is the original polyhedron.
       * \return A reference to the current object.
       */
      Polyhedron<R>& operator=(const Polyhedron<R>& original);
       //@}
      
      
      //@{
      //! \name Modifying operations
      /*! \brief Makes the polyhedron empty.
       */
      void clear();
      //@}
      
      
      //@{
      //! \name Geometric operations
      /*! \brief Returns the polyhedron's space dimension.
       */
      size_type dimension() const;
      
      /*! \brief Checks for emptyness.
       */
      bool empty() const;
      
   
      /*! \brief The vertices of the polyhedron. */
      state_list_type vertices() const; 
           
      /*! \brief Tests if a point is an element of the polyhedron.
       */
      bool contains(const state_type& point) const;

      /*! \brief Tests if the polyhedron contains a \a rectangle. */
      bool contains(const Rectangle<R>& rect) const; 
      
      /*! \brief Tests if a point is an element of the interior of the polyhedron.
       */
      bool interior_contains(const state_type& point) const;
      
      /*! \brief A rectangle containing the polyhedron. */
      Rectangle<R> bounding_box() const;
      
      /*! \brief An over-approximation of the polyhedron governed by the parameter \a delta.
       *
       * WARNING: The metric error of the approximation may be larger than \a delta.
       */
      Polyhedron<R> over_approximation(const real_type& delta) const;
      //@}
      

      //@{
      //! \name Geometric binary predicates
      /*! \brief Tests equality. */
      friend bool Geometry::equal <>(const Polyhedron<R>& A, 
                                        const Polyhedron<R>& B);
      /*! \brief Tests disjointness. */
      friend bool Geometry::disjoint <>(const Polyhedron<R>& A, 
                                        const Polyhedron<R>& B);
    
      /*! \brief Tests intersection of interiors. */
      friend bool Geometry::interiors_intersect<>(const Polyhedron<R>& A, 
                                                  const Polyhedron<R>& B);
      
      /*! \brief Tests inclusion of \a A in interior of \a B. */
      friend bool Geometry::inner_subset<>(const Polyhedron<R>& A, 
                                           const Polyhedron<R>& B);
    
      /*! \brief Tests inclusion of \a A in \a B. */
      friend bool Geometry::subset<>(const Polyhedron<R>& A, 
                                     const Polyhedron<R>& B);
    
      //@}
      

      //@{
      //! \name Geometric binary operations
      /*! \brief The intersection of two polyhedra. */
      friend Polyhedron<R> Geometry::intersection<>(const Polyhedron<R>& A, 
                                                    const Polyhedron<R>& B);
    
      /*! \brief The closure of the intersection of the interiors of two polyhedra. */
      friend Polyhedron<R> Geometry::regular_intersection<>(const Polyhedron<R>& A, 
                                                            const Polyhedron<R>& B);
    
      /*! \brief The convex hull of two polyhedra. */
      friend Polyhedron<R> Geometry::convex_hull<>(const Polyhedron<R>& A, 
                                                   const Polyhedron<R>& B);
      
      /*! \brief The Minkowski (pointwise) sum of two polyhedra. */
      friend Polyhedron<R> Geometry::minkowski_sum<>(const Polyhedron<R>& A, 
                                                     const Polyhedron<R>& B);
   
      /*! \brief The Minkowski (pointwise) difference of two polyhedra. */
      friend Polyhedron<R> Geometry::minkowski_difference<>(const Polyhedron<R>& A, 
                                                            const Polyhedron<R>& B);
       //@}
   
      friend std::ostream& operator<< <>(std::ostream& os, 
                                         const Polyhedron<R>& P);
    
     
     private:
      /* Construct a polyhedron from a PPL constraint system. */
      Polyhedron(Parma_Polyhedra_Library::Constraint_System& cs);
      
      const Parma_Polyhedra_Library::NNC_Polyhedron& _ppl_interior() const;
     private:
      /* The polyhedron data structure. */
      Parma_Polyhedra_Library::NNC_Polyhedron _ppl_poly; 
    
      /* This is mutable since we do not compute it unless needed */
      mutable Parma_Polyhedra_Library::NNC_Polyhedron _interior_poly;
    };


    template <typename R>
    bool 
    equal(const Polyhedron<R>& A, const Polyhedron<R>& B);
    
    template <typename R>
    bool 
    disjoint(const Polyhedron<R>& A, const Polyhedron<R>& B);
    
    template <typename R>
    bool
    interiors_intersect(const Polyhedron<R>& A, const Polyhedron<R>& B);
    
    template <typename R>
    bool 
    inner_subset(const Polyhedron<R>& A, const Polyhedron<R>& B);
    
    template <typename R>
    bool 
    subset(const Polyhedron<R>& A, const Polyhedron<R>& B);

    template <typename R>
    Polyhedron<R> 
    regular_intersection(const Polyhedron<R>& A, const Polyhedron<R>& B);

    template <typename R>
    Polyhedron<R> 
    intersection(const Polyhedron<R>& A, const Polyhedron<R>& B);

    template <typename R>
    Polyhedron<R> 
    convex_hull(const Polyhedron<R>& A, const Polyhedron<R>& B);

    template <typename R>
    Polyhedron<R> 
    minkowski_sum(const Polyhedron<R>& A, const Polyhedron<R>& B);

    template <typename R>
    Polyhedron<R> 
    minkowski_difference(const Polyhedron<R>& A, const Polyhedron<R>& B);

    template<typename R, template <typename> class BS> 
    inline Polyhedron<R> 
    minkowski_sum(const Polyhedron<R>& A, const BS<R>& B) {
      return minkowski_sum(A, Polyhedron<R>(B));
    }

    template <typename R>
    inline 
    bool 
    disjoint(const Polyhedron<R>& A, 
             const Rectangle<R>& B)
    {
      return disjoint(A,Polyhedron<R>(B));
    }
    
    template <typename R>
    inline 
    bool 
    disjoint(const Rectangle<R>& A, 
             const Polyhedron<R>& B)
    {
      return disjoint(B,A);
    }

    template <typename R>
    inline 
    bool 
    interiors_intersect(const Polyhedron<R>& A,
                        const Rectangle<R>& B) 
    {
       return interiors_intersect(A,Polyhedron<R>(B));
    }
    
    template <typename R>
    inline 
    bool 
    interiors_intersect(const Rectangle<R>& A, 
                        const Polyhedron<R>& B) 
    {
      return interiors_intersect(B,A);
    }    
    
    template <typename R>
    inline 
    bool 
    inner_subset(const Polyhedron<R>& A,
                 const Rectangle<R>& B) 
    {
      return inner_subset(A,Polyhedron<R>(B));
    }
    
    template <typename R>
    inline 
    bool 
    inner_subset(const Rectangle<R>& A, 
                 const Polyhedron<R>& B) 
    {
      return inner_subset(B,A);
    }
    
    template <typename R>
    inline 
    bool 
    subset(const Polyhedron<R>& A,
           const Rectangle<R>& B) 
    {
      return subset(A,Polyhedron<R>(B));
    }
    
    template <typename R>
    inline 
    bool 
    subset(const Rectangle<R>& A,
           const Polyhedron<R>& B) 
    {
      return subset(Polyhedron<R>(A),B);
    }

    template <typename R>
    inline 
    Polyhedron<R> 
    regular_intersection(const Polyhedron<R>& A, 
                         const Rectangle<R>& B) 
    {
      return regular_intersection(A, Polyhedron<R>(B));
    }

    template <typename R>
    inline 
    Polyhedron<R> 
    regular_intersection(const Rectangle<R>& A, 
                         const Polyhedron<R>& B) 
    {
      return regular_intersection(B,A);
    }

}
}

#endif /* _ARIADNE_POLYHEDRON_H */
