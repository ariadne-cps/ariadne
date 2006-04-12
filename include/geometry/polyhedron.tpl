/***************************************************************************
 *            polyhedron.tpl
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
 *  along with this program; if not, write to bouthe Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include <iostream>
#include <string>
#include <sstream>

#include <list>
#include <set>
#include <vector>
#include <valarray>

#include "../utility/stlio.h"
#include "../numeric/interval.h"

#include "../numeric/arithmetic.h"

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"
#include "../linear_algebra/constraint.h"

#include "../geometry/point.h"
#include "../geometry/rectangle.h"
#include "../geometry/list_set.h"
#include "../geometry/parallelotope.h"
#include "../geometry/polyhedron.h"

namespace Ariadne {
  namespace Geometry {
    
    template<typename R>
    Polyhedron<R>::Polyhedron()
      : _ppl_poly(Parma_Polyhedra_Library::EMPTY)
    { 
    }      
    
    template<typename R>
    Polyhedron<R>::Polyhedron(size_type dim)
      : _ppl_poly(dim, Parma_Polyhedra_Library::EMPTY)
    { 
    }
    
    template<typename R>
    Polyhedron<R>::Polyhedron(Parma_Polyhedra_Library::Constraint_System& cs)
      : _ppl_poly(cs)
    { 
    }
    
    template<typename R>
    Polyhedron<R>::Polyhedron(Ariadne::LinearAlgebra::ConstraintSystem<R>& cs)
      : _ppl_poly(Parma_Polyhedra_Library::Constraint_System(cs))
    {
    }
    
    template<typename R>
    Polyhedron<R>::Polyhedron(Ariadne::LinearAlgebra::GeneratorSystem<R>& gen)
      : _ppl_poly(Parma_Polyhedra_Library::Generator_System(gen)) 
    {
    }    
    
    template <typename R>
    Polyhedron<R>::Polyhedron(const matrix_type& A, const vector_type& b)
    //: _ppl_poly(A.size1()) 
    {
      vector_type rw;
      real_type cnst;
      Parma_Polyhedra_Library::Constraint_System cs;
      
      for(size_type i=0; i!=A.size1(); ++i) {
        //    rw=vector_type(row(A,i));
        rw=row(A,i);
        cnst=b[i];
        cs.insert(_convert_to_PPL_constraint(rw,cnst));
      }
      _ppl_poly=Parma_Polyhedra_Library::NNC_Polyhedron(cs);
    }
    

    template <typename R>
    Polyhedron<R>::Polyhedron(const state_list_type& points)
      : _ppl_poly() 
    {
      Parma_Polyhedra_Library::Generator_System ppl_gen;
      Parma_Polyhedra_Library::Linear_Expression ppl_lin_expr;
      
      size_type n=points[0].dimension();
      Ariadne::Geometry::Point<Rational> s(n);
      for(size_type i=0; i!=points.size(); ++i) {
        for(size_type j=0; j!=n; ++j) {
          s[j]=convert_to<Rational>(points[i][j]);
        }
        
        Integer den=1;
        for (size_type j=0; j!=n; ++j) {
          den=lcm(den,denominator(s[j]));
        }
        
        Parma_Polyhedra_Library::Linear_Expression ppl_lin_expr;
        for (size_type j=0; j!=n; ++j) {
          ppl_lin_expr+=( (numerator(s[j])*(den/denominator(s[j]))) * Parma_Polyhedra_Library::Variable(j) );
        }
        ppl_gen.insert(Parma_Polyhedra_Library::Generator::point(ppl_lin_expr,den));
      }
      _ppl_poly=Parma_Polyhedra_Library::NNC_Polyhedron(ppl_gen);
    }
    

    template <typename R>
    Polyhedron<R>::Polyhedron(const Rectangle<R>& r)
      : _ppl_poly(r.dimension())
    {
      state_type u_corner(r.upper_corner());
      state_type l_corner(r.lower_corner());
      
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
    
    template<typename R>
    Polyhedron<R>::Polyhedron(const Polyhedron<R>& original)
      : _ppl_poly(original._ppl_poly)
    {
    }
        
    template<typename R>
    Polyhedron<R>& 
    Polyhedron<R>::operator=(const Polyhedron<R>& original) 
    {        
      if(this != &original) { 
        this->_ppl_poly=original._ppl_poly; 
      }
      return *this;
    }
    


    template <typename R>
    size_type 
    Polyhedron<R>::dimension() const 
    {
      return ((size_type)this->_ppl_poly.space_dimension());
    }
    
    template <typename R>
    typename Polyhedron<R>::state_list_type
    Polyhedron<R>::vertices() const 
    {
      state_list_type result;
      const Parma_Polyhedra_Library::Generator_System& gs = this->_ppl_poly.minimized_generators();
      state_type point(this->dimension());
     
      for(Parma_Polyhedra_Library::Generator_System::const_iterator iter=gs.begin(); iter!=gs.end(); ++iter) {
        Parma_Polyhedra_Library::Generator gen = *iter;
        
        Parma_Polyhedra_Library::Coefficient_traits::const_reference denom = gen.divisor();
        for(size_type i=0; i!=this->dimension(); ++i) {
          Parma_Polyhedra_Library::Coefficient_traits::const_reference num = gen.coefficient(Parma_Polyhedra_Library::Variable(i));
          point[i] = convert_to<R>(num)/denom;
        }
        result.push_back(point);
      }
      
      return result;
    }
    
    /*! \brief Checks for emptyness.
     */
    template <typename R>
    bool 
    Polyhedron<R>::empty() const 
    {
      return (this->_ppl_poly.is_empty());
    }
    
    /*! \brief Makes the polyhedron empty.
     */
    template <typename R>
    void 
    Polyhedron<R>::clear() {
      this->_ppl_poly=Parma_Polyhedra_Library::NNC_Polyhedron(this->dimension(), Parma_Polyhedra_Library::EMPTY);;
    }

    /*! \brief Tests if a Point is an element of the Polyhedron.
     */
    template <typename R>
    bool 
    Polyhedron<R>::contains(const state_type& point) const {
      if (point.dimension()!=this->dimension()) {
        throw std::domain_error("This object and parameter have different space dimensions");
      }
      Parma_Polyhedra_Library::NNC_Polyhedron p(this->_ppl_poly);
      p.topological_closure_assign();
      
      return p.contains(_from_Point_to_PPL_Polyhedron(point));
      
    }
    
    /*! \brief Tests if a point is an element of the interior of the polyhedron.
     */
    template <typename R>
    bool 
    Polyhedron<R>::interior_contains(const state_type& point) const 
    {
      if (point.dimension()!=this->dimension()) {
        throw std::domain_error("This object and parameter have different space dimensions");
      }
      return this->_ppl_poly.contains(_from_Point_to_PPL_Polyhedron(point));
    }
    
    /*! \brief A rectangle containing the polyhedron. */
    template <typename R>
    Rectangle<R> 
    Polyhedron<R>::bounding_box() const 
    {
      throw std::runtime_error("Polyhedron::bounding_box() const not implemented.");
    }
    
    template <typename R>
    Polyhedron<R> 
    Polyhedron<R>::over_approximation(const real_type& delta) const
    {
      Ariadne::LinearAlgebra::ConstraintSystem<R> cs(this->_ppl_poly.constraints());
      cs.expand_by(delta);
      return Polyhedron(cs);
    }

    template <typename R>
    bool 
    disjoint(const Polyhedron<R>& A, const Polyhedron<R>& B)
    {
      return (A._ppl_poly).is_disjoint_from(B._ppl_poly);
    }
    
    template <typename R>
    bool
    interiors_intersect(const Polyhedron<R>& A, const Polyhedron<R>& B)
    {
      return !( (A._ppl_poly).is_disjoint_from(B._interior_poly) );        
    }
    
    template <typename R>
    inline 
    bool 
    inner_subset(const Polyhedron<R>& A, const Polyhedron<R>& B)
    {        
      return (B._ppl_interior()).contains(A._ppl_poly);
    }
    
    template <typename R>
    bool 
    subset(const Polyhedron<R>& A, const Polyhedron<R>& B)
    {
      return (B._ppl_poly).contains(A._ppl_poly);
    }

    template <typename R>
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
    Polyhedron<R> 
    convex_hull(const Polyhedron<R>& A, const Polyhedron<R>& B)
    {   
      Polyhedron<R> chull(A);
      chull._ppl_poly.poly_hull_assign_and_minimize(B._ppl_poly);
      return chull;
    }

    /* FIXME: TO REIMPLEMENT */
    template <typename R>
    Polyhedron<R> 
    minkowski_sum(const Polyhedron<R>& A, const Polyhedron<R>& B)
    {
      if ((!(A._ppl_poly).is_bounded()) || (!(B._ppl_poly).is_bounded())) {
        throw std::domain_error("Minkosky sum is implemented for bounded polyhedra only");
      }  
            
      Ariadne::LinearAlgebra::GeneratorSystem<R> gen_A(A._ppl_poly.generators());
            
      Polyhedron<R> sum=B;
      Polyhedron<R> new_poly;
            
      for (size_type k=0; k<gen_A.number_of_points(); k++) {
        Ariadne::LinearAlgebra::GeneratorSystem<R> gen_B(B._ppl_poly.generators());
        for (size_type i=0; i<gen_B.space_dimension(); i++) {
          for (size_type j=0; j<gen_B.number_of_points() ; j++) {
            gen_B._points(i,j)+=gen_A._points(i,k);
          }
        }
        new_poly._ppl_poly=Parma_Polyhedra_Library::NNC_Polyhedron(gen_B.ppl_generator_system());
        sum=Geometry::convex_hull(sum,new_poly);
      }      
      return sum;
    }
      

    template <typename R>
    std::ostream& 
    operator<<(std::ostream& os, const Polyhedron<R>& p) 
    {
      os << "Polyhedron(equations=[";
      Parma_Polyhedra_Library::IO_Operators::operator<<(os,p._ppl_poly);
      os << "], vertices=" << p.vertices() << ") ";
      return os;
    }
    
    
    
    template<typename R> 
    const Parma_Polyhedra_Library::NNC_Polyhedron& 
    Polyhedron<R>::_ppl_interior() const
    {      
      Ariadne::LinearAlgebra::ConstraintSystem<real_type> cs(this->_ppl_poly.constraints());
      this->_interior_poly=Parma_Polyhedra_Library::NNC_Polyhedron(cs.open_ppl_constraint_system());
      return this->_interior_poly;
    }
    

    template <typename R>
    Parma_Polyhedra_Library::NNC_Polyhedron 
    _from_Point_to_PPL_Polyhedron(const Point<R>& s)
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
    
    
    /* Convert the vector \a v and constant \a c into the constraint $v\cdot x\leq c$. */
    template<typename R>
    Parma_Polyhedra_Library::Constraint
    _convert_to_PPL_constraint(const Ariadne::LinearAlgebra::vector<R>&  v, 
                               const R& s) 
    {
      Parma_Polyhedra_Library::Linear_Expression e;
      Parma_Polyhedra_Library::Linear_Expression c;
      Integer multiplier = denominator(convert_to<Rational>(s));
      
      for (uint i = 0; i!=v.size(); ++i) {
        multiplier=Ariadne::lcm(multiplier,denominator(convert_to<Rational>(v[i])));
      }
      for (uint i = 0; i!=v.size(); ++i) {
        Rational coefficient = convert_to<Rational>(v[i]);
        Integer scaled_coefficient = numerator(coefficient)*(multiplier/denominator(coefficient));
        e += Parma_Polyhedra_Library::Variable(i) * convert_to<Integer>(scaled_coefficient);
      }
      Rational coefficient = convert_to<Rational>(s);
      Integer scaled_coefficient = numerator(coefficient)*(multiplier/denominator(coefficient));
      c = Parma_Polyhedra_Library::Linear_Expression(convert_to<Integer>(scaled_coefficient));
      
      return e<=c;
    }
    
  }
}
