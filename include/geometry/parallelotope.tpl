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

#include "parallelotope.h"

#include "../utility/stlio.h"
#include "../numeric/interval.h"

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"
#include "../linear_algebra/linear_program.h"

#include "../geometry/lattice_set.h" 
#include "../geometry/point.h"
#include "../geometry/rectangle.h"
#include "../geometry/list_set.h"
#include "../geometry/zonotope.h"
#include "../geometry/polyhedron.h"

namespace Ariadne {
  namespace Geometry {

    /*! \brief Tests if the parallelotope contains \a point. */
    template <typename R>
    bool Parallelotope<R>::contains(const state_type& point) const {
      matrix_type inv_gen(LinearAlgebra::inverse(this->generators()));
      state_type centre=this->centre();
      vector_type v(point.position_vector()-centre.position_vector());

      vector_type e=prod(inv_gen,v);

      for (size_t i=0; i<e.size(); i++) 
        if (abs(e(i))>1) return false;
      
      return true;
    }
    
    /*! \brief Tests if the interior of the parallelotope contains \a point. */
    template<typename R>
    bool Parallelotope<R>::interior_contains(const state_type& point) const {
      matrix_type inv_gen=LinearAlgebra::inverse(this->generators());
      state_type centre=this->centre();
      vector_type v(point.position_vector()-centre.position_vector());

      vector_type e=inv_gen*v;

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
        new_generators(i,j)/=2;
      }
      
      state_type new_centre=this->centre()-LinearAlgebra::Vector<R>(column(new_generators,j)/2);
      result.adjoin(Parallelotope<R>(new_centre,new_generators));
      new_centre=new_centre+LinearAlgebra::Vector<R>(column(new_generators,j));
      result.adjoin(Parallelotope(new_centre,new_generators));

      return result;
    }
    
    template <typename R>
    ListSet<R,Parallelotope>
    Parallelotope<R>::subdivide() const 
    {
      size_type n=this->dimension();
      ListSet<R,Geometry::Parallelotope> result(this->dimension());
      matrix_type new_generators=this->generators()/2;
      
      state_type first_centre=this->centre();
      for(size_type i=0; i!=n; ++i) {
        first_centre=first_centre-LinearAlgebra::Vector<R>((this->generator(i))/2);
      }
      
      array<index_type> lower(n,0);
      array<index_type> upper(n,2);
      array<index_type> finish(n,0);
      finish[n-1]=2;
      lattice_iterator end(finish,lower,upper);

      for(lattice_iterator iter(lower,lower,upper); iter!=end; ++iter) {
        array<index_type> ary=*iter;
        state_type new_centre=first_centre;
        for(size_type i=0; i!=n; ++i) {
          if(ary[i]==1) {
            new_centre=new_centre+this->generator(i);
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
      LinearAlgebra::Vector<F> diff=s-this->centre();
      LinearAlgebra::Matrix<F> inv = LinearAlgebra::inverse(this->generators());
      return prod(inv,diff);
    }


    template<typename R>
    std::vector< Point<Rational> > 
    Parallelotope<R>::vertices() const
    {
      std::vector< Point<Rational> > result;

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
    std::vector< Point<R> > 
    Parallelotope<R>::approximate_vertices() const
    {
      std::vector< Point<R> > result;
      
      dimension_type d=this->dimension();
      assert(d<32);      
      const Point<R>& c(this->centre());
      const LinearAlgebra::Matrix<R>& g(this->generators());
      LinearAlgebra::Vector<R> e(d);

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
    Parallelotope<R>::over_approximation(const Rectangle<R> &c, const LinearAlgebra::IntervalMatrix<R>& A)
    {
#ifdef DEBUG
      std::cerr << "IntervalParallelotope<R>::over_approximating_parallelotope() const" << std::endl;
#endif
      size_type n=c.dimension();
      
      Point<R> cmid=c.centre();
      LinearAlgebra::Matrix<R> Amid=A.centre();
      
      LinearAlgebra::Matrix<R> D(n,n);
      for(size_type i=0; i!=n; ++i) {
        D(i,i)=c[i].radius();
      }
      
      LinearAlgebra::IntervalMatrix<R> Ainv=LinearAlgebra::inverse(LinearAlgebra::IntervalMatrix<R>(Amid));
      R err = upper_norm(Ainv*D)+upper_norm(Ainv*A);

      return Geometry::Parallelotope<R>(cmid,err*Amid);
    }
    


    template <typename R>
    std::ostream&
    operator<<(std::ostream& os, const Parallelotope<R>& p) 
    {
      if(p.dimension() > 0) {
        os << "Parallelotope(\n  centre=" << p.centre();
        os << "\n  directions=" << p.generators();
        os << "\n) ";
      }

      return os;
    }
    
    template <typename R>
    std::istream& 
    operator>>(std::istream& is, Parallelotope<R>& p)
    {
      throw std::domain_error("Not implemented");
    }
      
  }
}
