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

#include <ppl.hh>
#include <iostream>
#include <vector>

#include "constraint.h"
#include "generator.h"

#include "state.h"
#include "rectangle.h"
#include "linear_algebra.h"

namespace Parma_Polyhedra_Library {
  // Import all the output operators into the main PPL namespace.
  using IO_Operators::operator<<;
}
 
namespace Ariadne {  
namespace Geometry {

template<typename R> class State;
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
  typedef ::Ariadne::Geometry::State<R> State;
  typedef std::vector<State> StateList;
  typedef ::Ariadne::LinearAlgebra::vector<R> Vector;
  typedef ::Ariadne::LinearAlgebra::matrix<R> Matrix;
 public:
    /*! \brief Default constructor creates empty set. 
     */
    inline Polyhedron()
      : _ppl_poly(Parma_Polyhedra_Library::Polyhedron::EMPTY)
    { }      
    
    /*! \brief Construct an empty polyhedron with dimension \a dim.
     *
     * Builds a polyhedron with dimension \a dim.
     * \param dim is the new polyhedron's number of dimensions.
     */
    inline Polyhedron(size_type dim)
      : _ppl_poly(dim, Parma_Polyhedra_Library::Polyhedron::EMPTY)
    { }
      
    /*! \brief Construct the polyhedron defined by the matrix equations \f$Ax\leq b\f$.
     */
    inline Polyhedron(const Matrix& A, const Vector& v);
    
    /*! \brief Construct a polyhedron from a constraint system.
     *
     * \param cs is the constraint system.
     */
    inline Polyhedron(Ariadne::LinearAlgebra::ConstraintSystem<R>& cs)
      : _ppl_poly(Parma_Polyhedra_Library::Constraint_System(cs))
    {
    }
    
    /*! \brief Construct a polyhedron from a generator system.
     *
     * \param gen is the generator system.
     */
    inline Polyhedron(Ariadne::LinearAlgebra::GeneratorSystem<R>& gen)
      : _ppl_poly(Parma_Polyhedra_Library::Generator_System(gen)) 
    {
    }    
    
    /*! \brief Construct the polyhedron defined as the convex hull of a list of states.
     */
    inline Polyhedron(const StateList& states);
    
    /*! \brief Construct a polyhedron from a rectangle. */
    Polyhedron(const Rectangle<R>& rect);
    
    /*! \brief Copy constructor. 
     */
    inline Polyhedron(const Polyhedron<R>& original)
      : _ppl_poly(original._ppl_poly)
    {
    }
        
    /*! \brief Copy assignment operator. 
     *
     * \param original is the original polyhedron.
     * \return A reference to the current object.
     */
    inline Polyhedron<R>& operator=(const Polyhedron<R>& original) {        
      if(this != &original) { this->_ppl_poly=original._ppl_poly; }
      return *this;
    }
    
    /*! \brief Returns the polyhedron's space dimension.
     */
    inline size_type dimension() const {
      return ((size_type)this->_ppl_poly.space_dimension());
    }
    
    /*! \brief Checks for emptyness.
     */
    inline bool empty() const {
      return (this->_ppl_poly.is_empty());
    }
    
    /*! \brief Makes the polyhedron empty.
     */
    inline void clear() {
      this->_ppl_poly=Parma_Polyhedra_Library::NNC_Polyhedron(this->dimension(), Parma_Polyhedra_Library::Polyhedron::EMPTY);;
    }

    /*! \brief The vertices of the Polyhedron. */
    StateList vertices() const; 
    
    /*! \brief Tests if a State is an element of the Polyhedron.
     */
    inline bool contains(const State& state) const {
      if (state.dimension()!=this->dimension()) {
        throw std::domain_error("This object and parameter have different space dimensions");
      }
      Parma_Polyhedra_Library::NNC_Polyhedron p(this->_ppl_poly);
      p.topological_closure_assign();
      
      return p.contains(_from_State_to_PPL_Polyhedron(state));
      
    }
    
    /*! \brief Tests if a state is an element of the interior of the polyhedron.
     */
    inline bool interior_contains(const State& state) const {
      if (state.dimension()!=this->dimension()) {
        throw std::domain_error("This object and parameter have different space dimensions");
      }
      return this->_ppl_poly.contains(_from_State_to_PPL_Polyhedron(state));
    }
    
    
    
    /*! \brief An over-approximation of the Polyhdron governed by the parameter \a delta.
     *
     * WARNING: The metric error of the approximation may be larger than \a delta.
     */
    inline Polyhedron<R> over_approximation(const Real& delta) {
      Ariadne::LinearAlgebra::ConstraintSystem<R> cs(this->_ppl_poly.constraints());
      cs.expand_by(delta);
      return Polyhedron(cs);
    }
    
    /*! \brief WHAT DOES THIS DO??. */
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
  
    /*! \brief WHAT DOES THIS DO??. */
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

    /*! \brief Project onto the given coordinates. */
    inline Polyhedron<R> project_on_dimensions(const std::vector<uint>& dims) const {
      #ifdef DEBUG
        std::cout << __FILE__ << ":" << __LINE__ << std::endl;
      #endif    
      if (dims.size()==0) {
        throw "Cannot project on zero dimensions.";  
      }
      if (dims.size()>this->dimension()) {
        throw "Cannot project on more dimensions than the polyhedron ones.";
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
                               const ListSet<R,::Ariadne::Geometry::Polyhedron>& LS);

    /*! \brief Tests inclusion. */
    friend bool subset<>(const Polyhedron<R>& A, 
                         const Polyhedron<R>& B);

    /*! \brief Tests inclusion. */
    friend bool subset<>(const Polyhedron<R>& A, 
                         const ListSet<R,::Ariadne::Geometry::Polyhedron>& C);

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
    inline Polyhedron(Parma_Polyhedra_Library::Constraint_System& cs)
      : _ppl_poly(cs)
    { }
    
    inline const Parma_Polyhedra_Library::NNC_Polyhedron& _ppl_interior() const;
 
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
const Parma_Polyhedra_Library::NNC_Polyhedron& 
Polyhedron<R>::_ppl_interior() const
{      
  #ifdef DEBUG
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
  #endif
  
  Ariadne::LinearAlgebra::ConstraintSystem<Real> cs(this->_ppl_poly.constraints());
  this->_interior_poly=Parma_Polyhedra_Library::NNC_Polyhedron(cs.open_ppl_constraint_system());
  
  return this->_interior_poly;
  #ifdef DEBUG
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
  #endif
}


/*! Convert the vector \a v and constant \a c into the constraint $v\cdot x\leq c$. */
template<typename R>
inline
Parma_Polyhedra_Library::Constraint
_convert_to_PPL_constraint(const ::Ariadne::LinearAlgebra::vector<R> & v, 
                           const R& s) 
{
  Parma_Polyhedra_Library::Linear_Expression e;
  Parma_Polyhedra_Library::Linear_Expression c;
  Integer multiplier = denominator(convert_to<Rational>(s));
  
  for (uint i = 0; i!=v.size(); ++i) {
    multiplier=::Ariadne::lcm(multiplier,denominator(convert_to<Rational>(v[i])));
  }
  for (uint i = 0; i!=v.size(); ++i) {
    Rational coefficient = convert_to<Rational>(v[i]);
    Integer scaled_coefficient = numerator(coefficient)*(multiplier/denominator(coefficient));
    e += Parma_Polyhedra_Library::Variable(i) * convert_to<int>(scaled_coefficient);
  }
  Rational coefficient = convert_to<Rational>(s);
  Integer scaled_coefficient = numerator(coefficient)*(multiplier/denominator(coefficient));
  c = Parma_Polyhedra_Library::Linear_Expression(convert_to<int>(scaled_coefficient));
  
  return e<=c;
}

template <typename R>
Polyhedron<R>::Polyhedron(const Matrix& A, const Vector& b)
  : _ppl_poly(A.size1()) 
{
  Vector rw;
  Real cnst;
  Parma_Polyhedra_Library::Constraint_System cs;
  
  for(size_type i=0; i!=A.size1(); ++i) {
    rw=row(A,i);
    cnst=b[i];
    cs.insert(_convert_to_PPL_constraint(rw,cnst));
  }
  _ppl_poly=Parma_Polyhedra_Library::NNC_Polyhedron(cs);
}
  
template <typename R>
Polyhedron<R>::Polyhedron(const std::vector<State>& states)
  : _ppl_poly() 
{
  Parma_Polyhedra_Library::Generator_System ppl_gen;
  Parma_Polyhedra_Library::Linear_Expression ppl_lin_expr;
  
  size_type n=states[0].dimension();
  ::Ariadne::Geometry::State<Rational> s(n);
  for(size_type i=0; i!=states.size(); ++i) {
    for(size_type j=0; j!=n; ++j) {
      s[j]=convert_to<Rational>(states[i][j]);
    }
  
    Integer den=1;
    for (size_type j=0; j!=n; ++j) {
      den=lcm(numerator(den),denominator(s[j]));
    }

    Parma_Polyhedra_Library::Linear_Expression ppl_lin_expr;
    for (size_type j=0; j!=n; ++j) {
      ppl_lin_expr+=( (numerator(s[j])*(den/denominator(s[j]))) * Parma_Polyhedra_Library::Variable(j) );
    }
    std::cerr << ppl_lin_expr << " " << den << std::endl;
    ppl_gen.insert(Parma_Polyhedra_Library::Generator::point(ppl_lin_expr,den));
  }
  _ppl_poly=Parma_Polyhedra_Library::NNC_Polyhedron(ppl_gen);
}
  
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


template <typename R>
Polyhedron<R>::Polyhedron(const Rectangle<R>& r)
  : _ppl_poly(r.dimension())
{
  State u_corner(r.upper_corner());
  State l_corner(r.lower_corner());
  
  Parma_Polyhedra_Library::Constraint_System cs;
  Rational num;
    
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
  
  _ppl_poly.add_constraints_and_minimize(cs);
}

   


}
}

#endif /* _ARIADNE_POLYHEDRON_H */
