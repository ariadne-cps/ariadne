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

#include "../base/array.h"

#include "../linear_algebra/constraint.h"
#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

#include "../geometry/generator.h"
#include "../geometry/point.h"
#include "../geometry/rectangle.h"

namespace Parma_Polyhedra_Library {
  // Import all the output operators into the main PPL namespace.
  using IO_Operators::operator<<;
}
 
namespace Ariadne {  
namespace Geometry {

template<typename R> class Point;
template<typename R> class Polyhedron;
template<typename R, template<typename> class BS> class ListSet;

template<typename R> bool disjoint(const Polyhedron<R>& A, 
                                   const Polyhedron<R>& B);
template<typename R> bool disjoint(const Polyhedron<R>& A, 
                                   const Rectangle<R>& R);
template<typename R> bool disjoint(const Rectangle<R>& R, 
                                   const Polyhedron<R>& B);
template<typename R> bool interiors_intersect(const Polyhedron<R>& A, 
                                              const Polyhedron<R>& B);
template<typename R> bool interiors_intersect(const Polyhedron<R>& A,
                                              const Rectangle<R>& R);
template<typename R> bool interiors_intersect(const Rectangle<R>& R, 
                                              const Polyhedron<R>& B);
template<typename R> bool inner_subset(const Polyhedron<R>& A, 
                                       const Polyhedron<R>& B);
template<typename R> bool inner_subset(const Polyhedron<R>& A,
                                       const Rectangle<R>& R);
template<typename R> bool inner_subset(const Rectangle<R>& R, 
                                       const Polyhedron<R>& B);
template<typename R> bool inner_subset(const Polyhedron<R>& A, 
                                       const ListSet<R,Polyhedron>& B);
template<typename R> bool subset(const Polyhedron<R>& A, 
                                 const Polyhedron<R>& B);
template<typename R> bool subset(const Polyhedron<R>& A, 
                                 const ListSet<R,Polyhedron>& B);
template<typename R> Polyhedron<R> intersection(const Polyhedron<R>& A, 
                                                const Polyhedron<R>& B);
template<typename R> Polyhedron<R> regular_intersection(const Polyhedron<R>& A, 
                                                        const Polyhedron<R>& B);
template<typename R> Polyhedron<R> convex_hull(const Polyhedron<R>& A, 
                                               const Polyhedron<R>& B);
template<typename R> Polyhedron<R> minkowski_sum(const Polyhedron<R>& A,
                                                 const Polyhedron<R>& B);
                                               
template<typename R> std::ostream& operator<<(std::ostream& os, 
                                              const Polyhedron<R>& P);
                                               

template<typename R> class Rectangle;
template<typename R> class Simplex;
template<typename R> class Parallelopiped;

/* Transforms an Point into a Parma_Polyhedra_Library::C_Polyhedron */
template <typename R>
Parma_Polyhedra_Library::NNC_Polyhedron 
_from_Point_to_PPL_Polyhedron(const Point<R>& s); 

/* Transforms an Rectangle into a Parma_Polyhedra_Library::C_Polyhedron */
template <typename R>
Parma_Polyhedra_Library::NNC_Polyhedron 
_from_Rectangle_to_closed_PPL_Polyhedron(const Rectangle<R>& r);

/*! \brief The polyhedron class.
 *
 * This is a wrapper from Parma Polyhedra Library Polyhedron class
 * to Ariadne representation.
 */ 
template <typename R>
class Polyhedron {
 public:
  typedef R real_type;
  typedef Point<R> state_type;
  typedef std::vector< Point<R> > state_list_type;
  typedef Ariadne::LinearAlgebra::vector<R> vector_type;
  typedef Ariadne::LinearAlgebra::matrix<R> matrix_type;
 public:
    /*! \brief Default constructor creates empty set. 
     */
    Polyhedron()
      : _ppl_poly(Parma_Polyhedra_Library::EMPTY)
    { }      
    
    /*! \brief Construct an empty polyhedron with dimension \a dim.
     *
     * Builds a polyhedron with dimension \a dim.
     * \param dim is the new polyhedron's number of dimensions.
     */
    Polyhedron(size_type dim)
      : _ppl_poly(dim, Parma_Polyhedra_Library::EMPTY)
    { }
      
    /*! \brief Construct a polyhedron from a constraint system.
     *
     * \param cs is the constraint system.
     */
    Polyhedron(Ariadne::LinearAlgebra::ConstraintSystem<R>& cs)
      : _ppl_poly(Parma_Polyhedra_Library::Constraint_System(cs))
    {
    }
    
    /*! \brief Construct a polyhedron from a generator system.
     *
     * \param gen is the generator system.
     */
    Polyhedron(Ariadne::LinearAlgebra::GeneratorSystem<R>& gen)
      : _ppl_poly(Parma_Polyhedra_Library::Generator_System(gen)) 
    {
    }    
    
    /*! \brief Construct the polyhedron defined by the matrix equations \f$Ax\leq b\f$.
     */
    Polyhedron(const matrix_type& A, const vector_type& b);
    
    /*! \brief Construct the polyhedron defined as the convex hull of a list of points.
     */
    Polyhedron(const state_list_type& points);
    
    /*! \brief Construct a polyhedron from a rectangle. */
    Polyhedron(const Rectangle<R>& rect);
    
    /*! \brief Copy constructor. 
     */
    Polyhedron(const Polyhedron<R>& original)
      : _ppl_poly(original._ppl_poly)
    {
    }
        
    /*! \brief Copy assignment operator. 
     *
     * \param original is the original polyhedron.
     * \return A reference to the current object.
     */
    Polyhedron<R>& operator=(const Polyhedron<R>& original) {        
      if(this != &original) { this->_ppl_poly=original._ppl_poly; }
      return *this;
    }
    
    /*! \brief Returns the polyhedron's space dimension.
     */
    size_type dimension() const {
      return ((size_type)this->_ppl_poly.space_dimension());
    }
    
    /*! \brief Checks for emptyness.
     */
    bool empty() const {
      return (this->_ppl_poly.is_empty());
    }
    
    /*! \brief Makes the polyhedron empty.
     */
    void clear() {
      this->_ppl_poly=Parma_Polyhedra_Library::NNC_Polyhedron(this->dimension(), Parma_Polyhedra_Library::EMPTY);;
    }

    /*! \brief The vertices of the Polyhedron. */
    state_list_type vertices() const; 
    
    /*! \brief Tests if a Point is an element of the Polyhedron.
     */
    bool contains(const state_type& point) const {
      if (point.dimension()!=this->dimension()) {
        throw std::domain_error("This object and parameter have different space dimensions");
      }
      Parma_Polyhedra_Library::NNC_Polyhedron p(this->_ppl_poly);
      p.topological_closure_assign();
      
      return p.contains(_from_Point_to_PPL_Polyhedron(point));
      
    }
    
    /*! \brief Tests if a point is an element of the interior of the polyhedron.
     */
    bool interior_contains(const state_type& point) const {
      if (point.dimension()!=this->dimension()) {
        throw std::domain_error("This object and parameter have different space dimensions");
      }
      return this->_ppl_poly.contains(_from_Point_to_PPL_Polyhedron(point));
    }
    
    /*! \brief A rectangle containing the polyhedron. */
    Rectangle<R> bounding_box() const {
      throw std::runtime_error("Polyhedron::bounding_box() const not implemented.");
    }
     
  public:
    friend class Parallelopiped<R>;
      
    /*! \brief Tests disjointness. */
    friend bool disjoint <>(const Polyhedron<R>& A, 
                            const Polyhedron<R>& B);

    /*! \brief Tests disjointness. */
    friend bool disjoint <>(const Polyhedron<R>& A, 
                            const Rectangle<R>& R);

    /*! \brief Tests disjointness. */
    friend bool disjoint <>(const Rectangle<R>& R, 
                            const Polyhedron<R>& B);

    /*! \brief Tests intersection of interiors. */
    friend bool interiors_intersect<>(const Polyhedron<R>& A, 
                                      const Polyhedron<R>& B);
  
    /*! \brief Tests intersection of interiors. */  
    friend bool interiors_intersect<>(const Polyhedron<R>& A,
                                      const Rectangle<R>& R);
    
    /*! \brief Tests intersection of interiors. */  
    friend bool interiors_intersect<>(const Rectangle<R>& R,
                                      const Polyhedron<R>& B);
    
    /*! \brief Tests inclusion of \a A in interior of \a B. */
    friend bool inner_subset<>(const Polyhedron<R>& A, 
                               const Polyhedron<R>& B);
  
    /*! \brief Tests inclusion of \a R in interior of \a B. */
    friend bool inner_subset<>(const Rectangle<R>& R, 
                               const Polyhedron<R>& B);

    /*! \brief Tests inclusion of \a A in interior of \a R. */
    friend bool inner_subset<>(const Polyhedron<R>& A, 
                               const Rectangle<R>& R);

    /*! \brief Tests inclusion of \a A in the interior of \a LS. */
    friend bool inner_subset<>(const Polyhedron<R>& A, 
                               const ListSet<R,Ariadne::Geometry::Polyhedron>& LS);

    /*! \brief Tests inclusion. */
    friend bool subset<>(const Polyhedron<R>& A, 
                         const Polyhedron<R>& B);

    /*! \brief Tests inclusion. */
    friend bool subset<>(const Polyhedron<R>& A, 
                         const ListSet<R,Ariadne::Geometry::Polyhedron>& C);

    /*! \brief Makes intersection. */
    friend Polyhedron<R> intersection<>(const Polyhedron<R>& A, 
                                        const Polyhedron<R>& B);

    /*! \brief Makes closure of intersection of interiors. */
    friend Polyhedron<R> regular_intersection<>(const Polyhedron<R>& A, 
                                                const Polyhedron<R>& B);

    /*! \brief The convex hull of two polyhedra. */
    friend Polyhedron<R> convex_hull<>(const Polyhedron<R>& A, 
                                       const Polyhedron<R>& B);
    
    /*! \brief The Minkowski (pointwise) sum. */
    friend Polyhedron<R> minkowski_sum<>(const Polyhedron<R>& A, 
                                         const Polyhedron<R>& B);

    /*! \brief Prints polyhedron. */
    friend std::ostream& operator<< <>(std::ostream& os, 
                                       const Polyhedron<R>& P);

  private:
    /* Construct a polyhedron from a PPL constraint system. */
    Polyhedron(Parma_Polyhedra_Library::Constraint_System& cs)
      : _ppl_poly(cs)
    { }
    
    const Parma_Polyhedra_Library::NNC_Polyhedron& _ppl_interior() const;
 
    /* The polyhedron data structure. */
    Parma_Polyhedra_Library::NNC_Polyhedron _ppl_poly; //(thanks Parma, I LOVE you! :-) )
  
    /* This is mutable since we do not compute it unless needed */
    mutable Parma_Polyhedra_Library::NNC_Polyhedron _interior_poly;
};


template <typename R>
inline 
bool 
disjoint(const Polyhedron<R>& A, 
         const Polyhedron<R>& B)
{
  return (A._ppl_poly).is_disjoint_from(B._ppl_poly);
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
                    const Polyhedron<R>& B)
{
  B._ppl_interior();
  return !( (A._ppl_poly).is_disjoint_from(B._interior_poly) );        
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
             const Polyhedron<R>& B)
{        
  return (B._ppl_interior()).contains(A._ppl_poly);
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
inner_subset(const Polyhedron<R>& A, 
             const ListSet<R,Polyhedron>& B)
{
  typename ListSet<R,Polyhedron>::const_iterator iter;

  /* FIXME! */
  for (iter = B.begin(); iter != B.end(); ++iter) {
    if (inner_subset(A,*iter)) {
      return true;
    }
  }
  throw("Cannot determine whether Polyhedron is a subset of ListSet<R,Polyhedron>");
}


template <typename R>
inline 
bool 
subset(const Polyhedron<R>& A,
       const Polyhedron<R>& B) 
{
  return (B._ppl_poly).contains(A._ppl_poly);
}


template <typename R>
inline 
bool 
subset(const Polyhedron<R>& A, 
       const ListSet<R,Polyhedron>& B)
{
  typename ListSet<R,Polyhedron>::const_iterator iter;

  /* FIXME! */
  for (iter = B.begin(); iter != B.end(); ++iter) {
    if (subset(A,*iter)) {
      return true;
    }
  }
  throw("Cannot determine whether Polyhedron is a subset of ListSet<R,Polyhedron>");
}


template <typename R>
inline 
Polyhedron<R> 
regular_intersection(const Polyhedron<R>& A, 
                     const Polyhedron<R>& B) 
{
  Polyhedron<R> intersctn(A);        
  intersctn._ppl_poly.intersection_assign(B._ppl_interior());
  return intersctn;
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


template <typename R>
inline 
Polyhedron<R> 
intersection(const Polyhedron<R>& A, 
             const Polyhedron<R>& B) 
{
  Polyhedron<R> intersctn(A);        
  intersctn._ppl_poly.intersection_assign(B._ppl_poly);
  return intersctn;
}


template <typename R>
inline 
Polyhedron<R> 
convex_hull(const Polyhedron<R>& A, 
            const Polyhedron<R>& B)
{
  Polyhedron<R> chull(A);
  chull._ppl_poly.poly_hull_assign_and_minimize(B._ppl_poly);
  return chull;
}


/* FIXME: TO REIMPLEMENT */
template <typename R>
inline 
Polyhedron<R> 
minkowski_sum(const Polyhedron<R>& A,
              const Polyhedron<R>& B) 
{
  #ifdef DEBUG
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
  #endif
  
  if ((!(A._ppl_poly).is_bounded()) || (!(B._ppl_poly).is_bounded())) {
    throw std::domain_error("Minkosky sum is implemented for bounded polyhedra only");
  }  
        
  Ariadne::LinearAlgebra::GeneratorSystem<R> gen_A(A._ppl_poly.generators());
        
  Polyhedron<R> sum=B;
  Polyhedron<R> new_poly;
        
  for (size_t k=0; k<gen_A.number_of_points(); k++) {
    Ariadne::LinearAlgebra::GeneratorSystem<R> gen_B(B._ppl_poly.generators());
    for (size_t i=0; i<gen_B.space_dimension(); i++) {
      for (size_t j=0; j<gen_B.number_of_points() ; j++) {
        gen_B._points(i,j)+=gen_A._points(i,k);
      }
    }
    new_poly._ppl_poly=Parma_Polyhedra_Library::NNC_Polyhedron(gen_B.ppl_generator_system());
    sum=convex_hull(sum,new_poly);
  }
  
  #ifdef DEBUG
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
  #endif
  
  return sum;
}

}
}

#endif /* _ARIADNE_POLYHEDRON_H */
