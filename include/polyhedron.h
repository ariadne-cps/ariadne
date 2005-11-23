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
 
#ifndef _ARIADNE_POLYHEDRON_H
#define _ARIADNE_POLYHEDRON_H

#include<ppl.hh>
#include<iostream>
#include<vector>

#include <constraint.h>
#include <generator.h>

#include <rectangle.h>

namespace Ariadne {  
namespace Geometry {

template<typename R> class State;
template<typename R> class Polyhedron;

template<typename R> inline bool disjoint(const Polyhedron<R>& A, 
                                   const Polyhedron<R>& B);
template<typename R> inline bool interiors_intersect(const Polyhedron<R>& A, 
                                              const Polyhedron<R>& B);
template<typename R> inline bool interiors_intersect(const Rectangle<R>& A, 
                                              const Polyhedron<R>& B);
template<typename R> inline bool interiors_intersect(const Polyhedron<R>& A,
                                              const Rectangle<R>& B);
template<typename R> inline bool inner_subset(const Polyhedron<R>& A, 
                                       const Polyhedron<R>& B);
template<typename R> inline bool inner_subset(const Rectangle<R>& A, 
                                       const Polyhedron<R>& B);
template<typename R> inline bool inner_subset(const Polyhedron<R>& A,
                                       const Rectangle<R>& B);
template<typename R> inline bool subset(const Polyhedron<R>& A, 
                                 const Polyhedron<R>& B);
template<typename R> inline bool subset(const Polyhedron<R>& A, 
                                 const std::vector< Polyhedron<R> >& B);
template<typename R> inline Polyhedron<R> regular_intersection(const Polyhedron<R>& A, 
                                                        const Polyhedron<R>& B);
template<typename R> inline Polyhedron<R> convex_hull(const Polyhedron<R>& A, 
                                               const Polyhedron<R>& B);
template<typename R> inline Polyhedron<R> minkowski_sum(const Polyhedron<R>& A,
                                                 const Polyhedron<R>& B);
                                               
template<typename R> inline std::ostream& operator<<(std::ostream& os, 
                                              const Polyhedron<R>& P);
                                               

/* Transforms an State into a Parma_Polyhedra_Library::C_Polyhedron */
template <typename R>
inline 
Parma_Polyhedra_Library::NNC_Polyhedron 
_from_State_to_PPL_Polyhedron(const State<R>& s);

/* Transforms an Rectangle into a Parma_Polyhedra_Library::C_Polyhedron */
template <typename R>
inline Parma_Polyhedra_Library::NNC_Polyhedron 
_from_Rectangle_to_closed_PPL_Polyhedron(const Rectangle<R>& r);




  

/*! \brief The polyhedron class.
 *
 * This is a wrapper from Parma Polyhedra Library Polyhedron class
 * to Ariadne representation.
 */ 
template <typename R>
class Polyhedron {
 public:
  typedef size_t size_type;
  typedef R Real;
  typedef State<R> State;
 public:
    /*! \brief Default constructor creates empty set. 
     */
    inline Polyhedron()
      : _ppl_poly(Parma_Polyhedra_Library::Polyhedron::EMPTY),
        _interior_poly(Parma_Polyhedra_Library::Polyhedron::EMPTY)
    { }      
    
    /*! \brief Construct an empty polyhedron with dimension \a dim.
     *
     * Builds a polyhedron with dimension \a dim.
     * \param dim is the new polyhedron's number of dimensions.
     */
    inline Polyhedron(size_t dim)
      : _ppl_poly(dim, Parma_Polyhedra_Library::Polyhedron::EMPTY), 
        _interior_poly(dim, Parma_Polyhedra_Library::Polyhedron::EMPTY) 
    { }
      
    /*! \brief A polyhedron constructor. 
     *
     * Builds a hypercube of dimension \f$dim\f$ centered in the 
     * origin.
     * \param dim is the new polyhedron's space dimension.
     * \param size is the dimension per space dimension.
     */
    inline Polyhedron(const size_type& dim, const Real& size) {
      
      #ifdef DEBUG
        std::cout << __FILE__ << ":" << __LINE__ << std::endl;
      #endif
      
      Parma_Polyhedra_Library::Constraint_System cs;
      Real num=numerator(abs(size/2));
      Real den=denominator(abs(size/2));
      
      for (size_t i=0; i< dim; i++) {
        cs.insert( den* Parma_Polyhedra_Library::Variable(i) >=  -num );
        cs.insert( den* Parma_Polyhedra_Library::Variable(i) <=  num );
      }
      
      Parma_Polyhedra_Library::NNC_Polyhedron new_poly(cs);
      
      this->_ppl_poly=new_poly;  
      this->_evaluate_interior();

      #ifdef DEBUG
        std::cout << __FILE__ << ":" << __LINE__ << std::endl;
      #endif

    }
  
    /*! \brief Copy constructor. 
     */
    inline Polyhedron(const Polyhedron<R>& original)
      : _ppl_poly(original._ppl_poly), _interior_poly(original._interior_poly){}
        
    /*! \brief Construct a polyhedron from a constraint system.
     *
     * \param cs is the constraint system.
     */
    inline Polyhedron(Ariadne::LinearAlgebra::ConstraintSystem<R>& cs)
      : _ppl_poly(Parma_Polyhedra_Library::Constraint_System(cs))
    {
      this->_evaluate_interior();
    }
    
    /*! \brief Construct a polyhedron from a generator system.
     *
     * \param gen is the generator system.
     */
    inline Polyhedron(Ariadne::LinearAlgebra::GeneratorSystem<R>& gen)
      : _ppl_poly(Parma_Polyhedra_Library::Generator_System(gen)) 
    {
      this->_evaluate_interior();
    }    
    
    /*! \brief A polyhedron constructor. 
     *
     * Builds a polyhedron from a rectangle.
     */
    Polyhedron(const Rectangle<R>& rect) {
      _ppl_poly=_from_Rectangle_to_closed_PPL_Polyhedron(rect);
      this->_evaluate_interior();
    }
    
    /*! \brief Returns polyhedron's space dimension.
     *
     * \return The polyhedron's space dimension.
     */
    inline size_t dimension() const {
      return ((size_t)this->_ppl_poly.space_dimension());
    }
    
    /*! \brief Checks the emptyness.
     *
     * \return \a true if the polyhedron is empty,
     * \a false otherwise.    
     */
    inline bool empty() const {
      return (this->_ppl_poly.is_empty());
    }
    
    inline void clear() {
      this->_ppl_poly=Parma_Polyhedra_Library::NNC_Polyhedron(this->dimension(), Parma_Polyhedra_Library::Polyhedron::EMPTY);;
      this->_interior_poly=_ppl_poly;
    }

    /*! \brief Tests if a state is included into a polyhedron.
     *
     * \param state is a state in the polyhedron's space.
     * \return \a true if the state is contained into 
     * the current polyhedron, \a false otherwise.
     */
    inline bool contains(const State& state) const {
      
      if (state.dimension()!=this->dimension()) 
        throw std::domain_error("This object and parameter have different space dimensions");
      
      Parma_Polyhedra_Library::NNC_Polyhedron p(this->_ppl_poly);
      
      p.topological_closure_assign();
      
      return p.contains(_from_State_to_PPL_Polyhedron(state));
      
    }
    
    /*! \brief Tests if a state is included into the interior of a polyhedron.
     *
     * \param state is a state in the polyhedron's space.
     * \return \a true if the state is contained into the interior of 
     * the current polyhedron, \a false otherwise.
     */
    inline bool interior_contains(const State& state) const {
      
      if (state.dimension()!=this->dimension()) 
        throw std::domain_error("This object and parameter have different space dimensions");
      
      return this->_ppl_poly.contains(_from_State_to_PPL_Polyhedron(state));
      
    }
    
    
    
    inline Polyhedron<R> operator+(const Polyhedron<R>& A) const{
      
      #ifdef DEBUG
        std::cout << __FILE__ << ":" << __LINE__ << std::endl;
      #endif    
      
      return minkowski_sum(*this,A);
  
    }
    
    inline Polyhedron<R> over_approximation(const Real& delta) {
      Ariadne::LinearAlgebra::ConstraintSystem<R> cs(this->_ppl_poly.constraints());
      cs.expand_by(delta);
      return Polyhedron(cs);
    }
    
    /*! \brief Copy assignment operator. 
     *
     * \param original is the original polyhedron.
     * \return A reference to the current object.
     */
    inline Polyhedron<R>& operator=(const Polyhedron<R>& original) {        
      this->_ppl_poly=original._ppl_poly;
      this->_interior_poly=original._interior_poly;
      return *this;
    }
    
    inline Polyhedron<R>& set_precision_to_upperapproximating(const Real& delta) {
      Real denum=denominator(delta);
      Ariadne::LinearAlgebra::ConstraintSystem<R> cs(this->_ppl_poly.constraints());
      if (cs.already_at_precision(denum)) { 
        return *this;
      }
      cs.reduce_precision_to_expanding(denum);
      
      Parma_Polyhedra_Library::NNC_Polyhedron new_poly(cs.open_ppl_constraint_system());
      if (!new_poly.is_empty()) {
        this->_ppl_poly=Parma_Polyhedra_Library::NNC_Polyhedron(cs.ppl_constraint_system());
        this->_interior_poly=new_poly;
      }
      
      return *this;
    }
  
    inline Polyhedron<R>& set_precision_to_upperapproximating_for_output(const Real& delta) {
      Real denum=denominator(delta);
      Ariadne::LinearAlgebra::ConstraintSystem<R> cs(this->_ppl_poly.constraints());
      if (cs.already_at_precision(denum)) {
        return *this;
      }
      cs.reduce_precision_to_expanding(denum);
      this->_ppl_poly=Parma_Polyhedra_Library::NNC_Polyhedron(cs.ppl_constraint_system());
      this->_interior_poly=Parma_Polyhedra_Library::NNC_Polyhedron(cs.open_ppl_constraint_system());
      return *this;
    }

    inline Polyhedron<R> project_on_dimensions(const std::vector<uint>& dims) const {
      #ifdef DEBUG
        std::cout << __FILE__ << ":" << __LINE__ << std::endl;
      #endif    
      if (dims.size()==0) {
        throw "I can not project on zero dimensions.";  
      }
      if (dims.size()>this->dimension()) {
        throw "I can not project on more dimensions than the polyhedron ones.";
      }
      boost::numeric::ublas::matrix<Real> projection_map=
          Ariadne::LinearAlgebra::zero_matrix<Real>(this->dimension());
      for (size_t i=0; i< dims.size(); i++) {
        projection_map(dims[i],dims[i])=1.0;
      }
      
      #ifdef DEBUG
        std::cout << __FILE__ << ":" << __LINE__ << std::endl;
      #endif  
      /* FIXME! */
      return *this;
    }
  
  public:
    /*! \brief Tests disjointness. */
    friend bool disjoint <>(const Polyhedron<R>& A, 
                            const Polyhedron<R>& B);

    /*! \brief Tests intersection of interiors. */
    friend bool interiors_intersect<>(const Polyhedron<R>& A, 
                                      const Polyhedron<R>& B);
  
    /*! \brief Tests intersection of interiors. */  
    friend bool interiors_intersect<>(const Polyhedron<R>& P,
                                      const Rectangle<R>& R);
    
    /*! \brief Tests inclusion of \a A in interior of \a B. */
    friend bool inner_subset<>(const Polyhedron<R>& A, 
                               const Polyhedron<R>& B);
  
    /*! \brief Tests inclusion of \a R in interior of \a P. */
    friend bool inner_subset<>(const Rectangle<R>& R, 
                               const Polyhedron<R>& P);

    /*! \brief Tests inclusion of \a P in interior of \a R. */
    friend bool inner_subset<>(const Polyhedron<R>& P,
                               const Rectangle<R>& R);
        
    /*! \brief Tests inclusion. */
    friend bool subset<>(const Polyhedron<R>& A, 
                         const Polyhedron<R>& B);

    /*! \brief Tests inclusion. */
    friend bool subset<>(const Polyhedron<R>& A, 
                         const std::vector< Polyhedron<R> >& C);

    /*! \brief Makes closure of intersection of interiors. */
    friend Polyhedron<R> regular_intersection<>(const Polyhedron<R>& A, 
                                                const Polyhedron<R>& B);

    /*! \brief The convex hull of two polyhedra. */
    friend Polyhedron<R> convex_hull<>(const Polyhedron<R>& A, 
                                       const Polyhedron<R>& B);
    
    /*! \brief The Monkowski (pointwise) sum. */
    friend Polyhedron<R> minkowski_sum<>(const Polyhedron<R>& A, 
                                         const Polyhedron<R>& B);

    /*! \brief Prints polyhedron. */
    friend std::ostream& operator<< <>(std::ostream& os, 
                                       const Polyhedron<R>& P);

   
  private:
    inline void _evaluate_interior();
 
    /* The polyhedron data structure. */
    Parma_Polyhedra_Library::NNC_Polyhedron _ppl_poly; //(thanks Parma, I LOVE you! :-) )
  
    /* This is mutable since we do not compute it unless needed */
    mutable Parma_Polyhedra_Library::NNC_Polyhedron _interior_poly;
  private:
    template <typename STATE>
    friend class PolyhedronMatlabExporter;
    
    template <typename STATE>
    friend class PolyhedronOneDimMatlabExporter;
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
interiors_intersect(const Polyhedron<R>& A, 
                    const Polyhedron<R>& B)
{
  #ifdef DEBUG
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
  #endif
  return !((A._ppl_poly).is_disjoint_from(B._interior_poly));        
}

template <typename R>
inline 
bool 
interiors_intersect(const Rectangle<R>& rect, 
                    const Polyhedron<R>& A) 
{
  return interiors_intersect(A,rect);
}

template <typename R>
inline 
bool 
interiors_intersect(const Polyhedron<R>& A,
                    const Rectangle<R>& rect) 
{
  Parma_Polyhedra_Library::NNC_Polyhedron P_rect=_from_Rectangle_to_open_PPL_Polihedron(rect);
  return !((A._ppl_poly).is_disjoint_from(P_rect));
}

template <typename R>
inline 
bool 
inner_subset(const Polyhedron<R>& A, 
             const Polyhedron<R>& B)
{        
  return (B._interior_poly).contains(A._ppl_poly);
}


template <typename R>
inline 
bool 
inner_subset(const Rectangle<R>& rect, 
             const Polyhedron<R>& A) 
{
  Parma_Polyhedra_Library::NNC_Polyhedron P_rect=_from_Rectangle_to_open_PPL_Polihedron(rect);          
  return (A._interior_poly).contains(P_rect);
}

template <typename R>
inline 
bool 
inner_subset(const Polyhedron<R>& A,
             const Rectangle<R>& rect) 
{
  /* FIXME! */
  Parma_Polyhedra_Library::NNC_Polyhedron P_rect=_from_Rectangle_to_open_PPL_Polihedron(rect);          
  return (P_rect).contains(A);
}

template <typename R>
inline 
bool 
subset(const Polyhedron<R>& A,
       const Polyhedron<R>& B) 
{
  return B._ppl_poly.contains(A._ppl_poly);
}


template <typename R>
inline 
bool 
subset(const Polyhedron<R>& A, 
       const std::vector< Polyhedron<R> >& vector)
{
  typename std::vector< Polyhedron <R> >::const_iterator i;

  /* FIXME! */
  for (i = vector.begin(); i != vector.end(); i++) {
    if (inner_subset(A,*i)) {
      return true;
    }
  }
  return false;
}


template <typename R>
inline 
Polyhedron<R> 
regular_intersection(const Polyhedron<R>& A, 
                     const Polyhedron<R>& B) 
{
  Polyhedron<R> new_poly(A);        
    
  new_poly._ppl_poly.intersection_assign(B._ppl_poly);
  new_poly._interior_poly.intersection_assign(B._interior_poly);
  return new_poly;
}


template <typename R>
inline 
Polyhedron<R> 
convex_hull(const Polyhedron<R>& A, 
            const Polyhedron<R>& B)
{
  Polyhedron<R> chull(A);
  chull._ppl_poly.poly_hull_assign_and_minimize(B._ppl_poly);
  chull._evaluate_interior();
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
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
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
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
  #endif
  
  return sum;
}


template <typename R>
inline
std::ostream& 
operator<<(std::ostream &os, const Polyhedron<R>& p) {
  Parma_Polyhedra_Library::IO_Operators::operator<<(os,p._ppl_poly);
  return os;
}


template<typename R> 
inline
void 
Polyhedron<R>::_evaluate_interior()
{      
  #ifdef DEBUG
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
  #endif
  
  Ariadne::LinearAlgebra::ConstraintSystem<Real> cs(this->_ppl_poly.constraints());
  this->_interior_poly=Parma_Polyhedra_Library::NNC_Polyhedron(cs.open_ppl_constraint_system());
  
  #ifdef DEBUG
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
  #endif
}      





/* Transforms an State into a Parma_Polyhedra_Library::C_Polyhedron */
template <typename R>
inline 
Parma_Polyhedra_Library::NNC_Polyhedron 
_from_State_to_PPL_Polyhedron(const State<R>& s)
{
  Parma_Polyhedra_Library::Constraint_System cs;
  Rational num;
  
  for (size_t i=0; i<s.dimension(); i++) {
      num=convert_to<Rational>(s[i]);
      cs.insert(Parma_Polyhedra_Library::Variable(i) * 
                Ariadne::denominator(num) == Ariadne::numerator(num));
  }
  
  return Parma_Polyhedra_Library::NNC_Polyhedron(cs);
}


/* Transforms an Rectangle into a Parma_Polyhedra_Library::C_Polyhedron */
template <typename R>
inline Parma_Polyhedra_Library::NNC_Polyhedron 
_from_Rectangle_to_closed_PPL_Polyhedron(const Rectangle<R>& r)
{
  State<R> u_corner(r.upper_corner()), l_corner(r.lower_corner());
  
  Parma_Polyhedra_Library::Constraint_System cs;
  std::vector<Parma_Polyhedra_Library::Variable *> x;
  Rational num;
  
  for (size_t i=0; i<r.dimension(); i++) {
    x.push_back(new Parma_Polyhedra_Library::Variable(i));
  }
  
  for (size_t i=0; i<r.dimension(); i++) {
    num=convert_to<Rational>(u_corner[i]);
    cs.insert(
      (Parma_Polyhedra_Library::Variable(i))*denominator(num) <= numerator(num)
    );
  
    num=convert_to<Rational>(l_corner[i]);
    cs.insert(
      (Parma_Polyhedra_Library::Variable(i))*denominator(num) >= numerator(num)
    );
  }
  
  return Parma_Polyhedra_Library::NNC_Polyhedron(cs);
}


}
}

#endif /* _ARIADNE_POLYHEDRON_H */
