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

#include <iosfwd>

#include <string>
#include <sstream>

#include <list>
#include <set>
#include <vector>
#include <valarray>

#include "../base/utility.h"
#include "../base/interval.h"

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
    
    template <typename R>
    Polyhedron<R>::Polyhedron(const Matrix& A, const Vector& b)
    //: _ppl_poly(A.size1()) 
    {
      Vector rw;
      Real cnst;
      Parma_Polyhedra_Library::Constraint_System cs;
      
      for(size_type i=0; i!=A.size1(); ++i) {
        //    rw=Vector(row(A,i));
        rw=row(A,i);
        cnst=b[i];
        cs.insert(_convert_to_PPL_constraint(rw,cnst));
      }
      _ppl_poly=Parma_Polyhedra_Library::NNC_Polyhedron(cs);
    }
    

    template <typename R>
    Polyhedron<R>::Polyhedron(const StateList& points)
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
    


    template <typename R>
    typename Polyhedron<R>::StateList
    Polyhedron<R>::vertices() const 
    {
      StateList result;
      const Parma_Polyhedra_Library::Generator_System& gs = this->_ppl_poly.minimized_generators();
      State point(this->dimension());
      
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
    
    

    template <typename R>
    std::ostream& 
    operator<<(std::ostream& os, const Polyhedron<R>& p) {
      os << "Polyhedron(equations=[";
      Parma_Polyhedra_Library::IO_Operators::operator<<(os,p._ppl_poly);
      os << "], vertices=";
      os << p.vertices();
      os << ") ";
      return os;
    }
    
    
    
    template<typename R> 
    const Parma_Polyhedra_Library::NNC_Polyhedron& 
    Polyhedron<R>::_ppl_interior() const
    {      
      Ariadne::LinearAlgebra::ConstraintSystem<Real> cs(this->_ppl_poly.constraints());
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
