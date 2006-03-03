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

#include "../geometry/lattice_set.h" 
#include "../geometry/point.h"
#include "../geometry/rectangle.h"
#include "../geometry/list_set.h"
#include "../geometry/parallelotope.h"
#include "../geometry/polyhedron.h"

namespace Ariadne {
  namespace Geometry {

    template<typename R>
    Rectangle<R> 
    Parallelotope<R>::bounding_box() const 
    {
      Vector offset(this->dimension());
      for(size_type i=0; i!=this->dimension(); ++i) {
        for(size_type j=0; j!=this->dimension(); ++j) {
          offset[i] += abs(this->_generators(i,j));
        }
      }
      Rectangle<R> result(this->centre()+offset, this->centre()-offset);
      return result;
    }
      
    
    template <typename R>
    Parallelotope<R>::operator Polyhedron<R>() const 
    {
      using namespace ::Ariadne::LinearAlgebra;
      
      typedef typename Parallelotope<R>::Real Real;
      typedef typename Parallelotope<R>::State State;
      
      size_type n = this->dimension();
      
      /* Express in form invs * x - offst in [-bnds,+bnds] */
      matrix<R> invs;
      vector<R> offst;
      vector<R> bnds;
      this->compute_linear_inequalities(invs,offst,bnds);
      
     
      matrix<R> A(2*n,n);
      vector<R> b(2*n);
      
      for(uint i=0; i!=n; ++i) {
        for(uint j=0; j!=n; ++j) {
          A(i,j) = -invs(i,j);
          A(i+n,j) = invs(i,j);
        }
        b(i) = bnds(i)-offst(i);
        b(i+n) = bnds(i)+offst(i);
      }
      return Polyhedron<R>(A,b);
    }


    template <typename R>
    ListSet<R,Parallelotope>
    Parallelotope<R>::subdivide() const 
    {
      size_type n=this->dimension();
      ListSet<R,Geometry::Parallelotope> result(this->dimension());
      Matrix new_generators=this->generators()/2;
      
      State first_centre=this->centre();
      for(size_type i=0; i!=n; ++i) {
        first_centre=first_centre-(this->generator(i))/2;
      }
      
      array<index_type> lower(n,0);
      array<index_type> upper(n,2);
      array<index_type> finish(n,0);
      finish[n-1]=2;
      lattice_iterator end(finish,lower,upper);

      for(lattice_iterator iter(lower,lower,upper); iter!=end; ++iter) {
        array<index_type> ary=*iter;
        State new_centre=first_centre;
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
    void 
    Parallelotope<R>::compute_linear_inequalities(Matrix& A, Vector& o, Vector& b) const
    {
      using namespace ::Ariadne::LinearAlgebra;
      size_type n=this->dimension();
      
      Vector c=this->centre() - State(n,0);
      A=inverse(this->generators());
      o=A*c;
      b=vector<R>(n);
      for(size_type i=0; i!=n; ++i) {
        b[i]=1;
      }
    }
      
    template<>
    void 
    Parallelotope<Dyadic>::compute_linear_inequalities(Matrix& A, Vector& o, Vector& b) const
    {
      using namespace ::Ariadne::LinearAlgebra;
      typedef Dyadic Real;
      size_type n=this->dimension();
      
      matrix<Rational> M(n,n);
      for(size_type i=0; i!=n; ++i) {
        for(size_type j=0; j!=n; ++j) {
          M(i,j) = convert_to<Rational>(this->_generators(i,j));
        }
      }
      M=inverse(M);
      
      vector<Integer> multipliers = row_common_denominators(M);
      
      A.resize(n,n);
      for(size_type i=0; i!=n; ++i) {
        for(size_type j=0; j!=n; ++j) {
          A(i,j) = numerator(M(i,j)) * (multipliers(i)/denominator(M(i,j)));
        }
      }

      Vector c = this->centre() - State(n,0);
      o = A * c;
      
      b.resize(n);
      for(size_type i=0; i!=n; ++i) {
        b(i)=multipliers(i);
      }
    }
      
  
  
    template<typename R>
    typename Parallelotope<R>::Vector 
    Parallelotope<R>::coordinates(const State& s) const {
      Vector diff = s-_centre;
      Matrix inv = LinearAlgebra::inverse(_generators);
      return prod(inv,diff);
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

