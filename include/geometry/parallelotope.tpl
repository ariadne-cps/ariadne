/***************************************************************************
 *            parallelotope.tpl
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

#include <vector>

#include "../numeric/rational.h"
#include "../numeric/interval.h"

#include "../base/stlio.h"

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"
#include "../linear_algebra/linear_program.h"

#include "../combinatoric/lattice_set.h" 

#include "../geometry/point.h"
#include "../geometry/point_list.h"
#include "../geometry/rectangle.h"
#include "../geometry/list_set.h"
#include "../geometry/zonotope.h"
#include "../geometry/polyhedron.h"

#include "parallelotope.h"

namespace Ariadne {
  namespace Geometry {

    /*! \brief Tests if the parallelotope contains \a point. */
    template <typename R>
    bool Parallelotope<R>::contains(const Point<R>& point) const {
      LinearAlgebra::Vector<Rational> p=point.position_vector();
      LinearAlgebra::Vector<Rational> c=this->centre().position_vector();
      LinearAlgebra::Matrix<Rational> G=this->generators();

      LinearAlgebra::Vector<Rational> e=G.solve(c-p);

      for (size_t i=0; i<e.size(); i++) 
        if (abs(e(i))>1) return false;
      
      return true;
    }
    
    /*! \brief Tests if the interior of the parallelotope contains \a point. */
    template<typename R>
    bool Parallelotope<R>::interior_contains(const Point<R>& point) const 
    {
      LinearAlgebra::Vector<Rational> p=point.position_vector();
      LinearAlgebra::Vector<Rational> c=this->centre().position_vector();
      LinearAlgebra::Matrix<Rational> G=this->generators();

      LinearAlgebra::Vector<Rational> e=G.solve(c-p);

      for (size_t i=0; i<e.size(); i++) 
        if (abs(e(i))>=1) return false;
      
      return true;  
    }
      
    template <typename R>
    ListSet<R,Parallelotope>
    Parallelotope<R>::divide() const 
    {
      size_type n=this->dimension();
      ListSet<R,Geometry::Parallelotope> result(this->dimension());
      
      matrix_type new_generators=this->generators();
      
      R two=2;
      R max_norm=0;
      size_type max_column=0;
      for(size_type j=0; j!=n; ++j) {
        R norm = LinearAlgebra::norm(LinearAlgebra::Vector<R>(column(new_generators,j)));
        if(norm>max_norm) {
          max_norm=norm;
          max_column=j;
        }
      }
      
      size_type j=max_column;
      for(size_type i=0; i!=n; ++i) {
        new_generators(i,j)=div_up(new_generators(i,j),two);
      }
      
      state_type new_centre=sub_approx(this->centre(),div_approx(LinearAlgebra::Vector<R>(column(new_generators,j)),two));
      result.adjoin(Parallelotope<R>(new_centre,new_generators));
      new_centre=add_approx(new_centre,LinearAlgebra::Vector<R>(column(new_generators,j)));
      result.adjoin(Parallelotope(new_centre,new_generators));

      return result;
    }
    
    template <typename R>
    ListSet<R,Parallelotope>
    Parallelotope<R>::subdivide() const 
    {
      ListSet<R,Geometry::Parallelotope> result(this->dimension());
      
      R two=2;
      size_type n=this->dimension();
      
      matrix_type new_generators(n,n);
      for(size_type i=0; i!=n; ++i) {
        for(size_type j=0; j!=n; ++j) {
          new_generators(i,j)=div_up(this->generators()(i,j),two);
        }
      }
      
      state_type first_centre=this->centre();
      for(size_type i=0; i!=n; ++i) {
        first_centre=sub_approx(first_centre,div_approx(LinearAlgebra::Vector<R>(this->generator(i)),two));
      }
      
      Point<R> new_centre;
      for(size_type k=0; k!=1<<n; ++k) {
        new_centre=first_centre;
        for(size_type i=0; i!=n; ++i) {
          if(k&(1<<i)) {
            new_centre=add_approx(new_centre,this->generator(i));
          }
        }
        result.adjoin(Parallelotope(new_centre,new_generators));
      }
      return result;
    }
    
    
    template<typename R>
    LinearAlgebra::Vector<typename numerical_traits<R>::field_extension_type>
    Parallelotope<R>::coordinates(const state_type& s) const {
      typedef typename numerical_traits<R>::field_extension_type F;
      LinearAlgebra::Vector<F> p=s.position_vector();
      LinearAlgebra::Vector<F> c=this->centre().position_vector();
      LinearAlgebra::Vector<F> d=p-c;
      LinearAlgebra::Matrix<F> G=this->generators();
      return G.solve(d);
    }


    template<typename R>
    PointList<Rational>
    Parallelotope<R>::vertices() const
    {
      PointList<Rational> result;

      dimension_type d=this->dimension();
      assert(d<32);      
      Point<Rational> c(this->centre());
      LinearAlgebra::Matrix<Rational> g(this->generators());
      LinearAlgebra::Vector<Rational> e(d);

      size_type nv=(1<<d);
      result.reserve(nv);

      for (size_type i=0; i<nv; ++i) {
        for(size_type j=0; j!=d; ++j) {
          e[j]=(i&(1<<d) ? 1 : -1);
        }
        result.push_back(c+g*e);
      }
      return result;
    }
    
    

    
    template<typename R>
    Parallelotope<R>
    Parallelotope<R>::over_approximation(const Rectangle<R> &c, const LinearAlgebra::Matrix< Interval<R> >& A)
    {
#ifdef DEBUG
      std::cerr << "IntervalParallelotope<R>::over_approximating_parallelotope() const" << std::endl;
#endif
      size_type n=c.dimension();
      
      Point<R> cmid=c.centre();
      LinearAlgebra::Matrix<R> Amid=LinearAlgebra::centre(A);
      
      LinearAlgebra::Matrix<R> D(n,n);
      for(size_type i=0; i!=n; ++i) {
        D(i,i)=c[i].radius();
      }
      
      LinearAlgebra::Matrix< Interval<R> > Ainv=LinearAlgebra::inverse(LinearAlgebra::Matrix< Interval<R> >(Amid));
      R err = ((Ainv*D).norm()+(Ainv*A).norm()).upper();

      for(size_type i=0; i!=n; ++i) {
        for(size_type j=0; j!=n; ++j) {
          Amid(i,j)=mul_up(Amid(i,j),err);
        }
      }
      
      return Geometry::Parallelotope<R>(cmid,Amid);
    }
    
/*
    template<typename R>
    Parallelotope<R> 
    Parallelotope<R>::scale(const Parallelotope<R>& p, const R& scale_factor) {

      const Point<R>& centre=p.centre();
      const LinearAlgebra::Matrix<R>& generators=p.generators();
      
      Point<R> new_centre(p.dimension());

      for(size_type i=0; i!=p.dimension(); ++i) {
        new_centre[i]=mul_approx(scale_factor,centre[i]);
      }

      return Parallelotope<R>(new_centre, scale_factor*generators);
    }
*/
    
    template <typename R>
    std::ostream&
    Parallelotope<R>::write(std::ostream& os) const
    {
      const Parallelotope<R>& p=*this;
      if(p.dimension() > 0) {
        os << "Parallelotope(\n  centre=" << p.centre();
        os << "\n  directions=" << p.generators();
        os << "\n) ";
      }

      return os;
    }
    
    template <typename R>
    std::istream& 
    Parallelotope<R>::read(std::istream& is)
    {
      throw std::domain_error("Not implemented");
    }
      
  }
}
