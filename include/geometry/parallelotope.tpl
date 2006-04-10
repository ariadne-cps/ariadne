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

#include "../declarations.h"

#include "../utility/stlio.h"
#include "../numeric/interval.h"

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"
#include "../linear_algebra/linear_program.h"

#include "../geometry/lattice_set.h" 
#include "../geometry/point.h"
#include "../geometry/rectangle.h"
#include "../geometry/list_set.h"
#include "../geometry/parallelotope.h"
#include "../geometry/zonotope.h"
#include "../geometry/polyhedron.h"

namespace Ariadne {
  namespace Geometry {

    template<typename R>
    Rectangle<R> 
    Parallelotope<R>::bounding_box() const 
    {
      vector_type offset(this->dimension());
      for(size_t i=0; i!=this->dimension(); ++i) {
        for(size_t j=0; j!=this->dimension(); ++j) {
          offset[i] += abs(this->_generators(i,j));
        }
      }
      return Rectangle<R>(this->centre()-offset, this->centre()+offset);
    }
   
    template <typename R>
    Parallelotope<R>::operator Zonotope<R>() const 
    {
      return Zonotope<R>(*this); 
    }

    
    template <typename R>
    Parallelotope<R>::operator Polyhedron<R>() const 
    {
      using namespace ::Ariadne::LinearAlgebra;
      
      size_type n = this->dimension();
      
      /* Express in form invs * x - offst in [-bnds,+bnds] */
      matrix_type invs;
      vector_type offst;
      vector_type bnds;
      this->compute_linear_inequalities(invs,offst,bnds);
      
      matrix_type A(2*n,n);
      vector_type b(2*n);
      
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

    /*! \brief The equality operator */
    template <typename R>
    bool Parallelotope<R>::operator==(const Parallelotope<R>& A) const {
      Zonotope<R> z_this(*this), z_A(A);

      return z_this==z_A;
    }
      
    /*! \brief The inequality operator */
    template <typename R>
    bool Parallelotope<R>::operator!=(const Parallelotope<R>& A) const {
      return !(*this == A);
    }  
     
    /*! \brief Tests if the parallelotope contains \a point. */
    template <typename R>
    bool Parallelotope<R>::contains(const state_type& point) const {
      matrix_type inv_gen(LinearAlgebra::inverse(this->_generators));
      const state_type &centre=this->_centre;
      vector_type v(point.position_vector()-centre.position_vector());

      vector_type e=prod(inv_gen,v);

      for (size_t i=0; i<e.size(); i++) 
        if (abs(e(i))>1) return false;
      
      return true;
    }
      
    /*! \brief Tests if the interior of the parallelotope contains \a point. */
    template<typename R>
    bool Parallelotope<R>::interior_contains(const state_type& point) const {
      matrix_type inv_gen=LinearAlgebra::inverse(this->_generators);
      const state_type &centre=this->_centre;
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
      
      matrix_type new_generators=generators();
      
      R max_norm=0;
      size_type max_column=0;
      for(size_type j=0; j!=n; ++j) {
        R norm = LinearAlgebra::norm(LinearAlgebra::vector<R>(column(new_generators,j)));
        if(norm>max_norm) {
          max_norm=norm;
          max_column=j;
        }
      }
      
      size_type j=max_column;
      for(size_type i=0; i!=n; ++i) {
        new_generators(i,j)/=2;
      }
      
      state_type new_centre=this->centre()-column(new_generators,j)/2;
      result.adjoin(Parallelotope<R>(new_centre,new_generators));
      new_centre=new_centre+column(new_generators,j);
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
        first_centre=first_centre-(this->generator(i))/2;
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
    void 
    Parallelotope<R>::compute_linear_inequalities(matrix_type& A, vector_type& o, vector_type& b) const
    {
      _compute_linear_inequalities(A,o,b,this->centre().position_vector(),this->generators(),typename numerical_traits<R>::algebraic_category());
    }
    
    template<typename R>
    inline
    void
    _compute_linear_inequalities(LinearAlgebra::matrix<R>& A,   
                                 LinearAlgebra::vector<R>& o, 
                                 LinearAlgebra::vector<R>& b, 
                                 const LinearAlgebra::vector<R>& c, 
                                 const LinearAlgebra::matrix<R>& G, 
                                 const field_tag&)
    {
      size_type n=c.size();
      
      A=LinearAlgebra::inverse(G);
      o=A*c;
      
      b=LinearAlgebra::vector<R>(n);
      for(size_type i=0; i!=n; ++i) {
        b[i]=1;
      }
    }
      
    template<typename R>
    inline
    void
    _compute_linear_inequalities(LinearAlgebra::matrix<R>& A,   
                                 LinearAlgebra::vector<R>& o, 
                                 LinearAlgebra::vector<R>& b, 
                                 const LinearAlgebra::vector<R>& c, 
                                 const LinearAlgebra::matrix<R>& G, 
                                 const ring_tag&)
    {
      typedef typename numerical_traits<R>::field_extension_type F;
      
      size_type n=c.size();
      
      LinearAlgebra::matrix<F> M(n,n);
      for(size_type i=0; i!=n; ++i) {
        for(size_type j=0; j!=n; ++j) {
          M(i,j) = convert_to<F>(G(i,j));
        }
      }
      M=LinearAlgebra::inverse(M);
      
      LinearAlgebra::vector<Integer> multipliers = row_common_denominators(M);
      
      A.resize(n,n);
      for(size_type i=0; i!=n; ++i) {
        for(size_type j=0; j!=n; ++j) {
          A(i,j) = Integer(numerator(M(i,j)) * (multipliers(i)/denominator(M(i,j))));
        }
      }

      o = A * c;
      
      b.resize(n);
      for(size_type i=0; i!=n; ++i) {
        b(i)=multipliers(i);
      }
    }
      
  
  
    template<typename R>
    LinearAlgebra::vector<typename numerical_traits<R>::field_extension_type>
    Parallelotope<R>::coordinates(const state_type& s) const {
      typedef typename numerical_traits<R>::field_extension_type F;
      LinearAlgebra::vector<F> diff(this->dimension());
      for(size_type i=0; i!=diff.size(); ++i) {
        diff(i)=s[i]-this->_centre[i];
      }
      LinearAlgebra::matrix<F> inv = LinearAlgebra::inverse(this->_generators);
      return prod(inv,diff);
    }

    template<typename R>
    bool
    Parallelotope<R>::disjoint(const Rectangle<R>& r) const
    {
      assert(this->dimension()==r.dimension());
      dimension_type n=this->dimension();
      
      // Construct tableau for testing intersection of rectangle and point
      // Rectangle  a<=x<=b
      // Parallelotope  x==c+Ae,  -1<=e<=1
      // 
      // Translate x'=x-a,  e'=e+1
      //   0<=x'<=b-a
      //   0<=e'<=2
      //   x'+a==c+A(e'-1) ->  x'-Ae' == c-a-A1
      //  
      // Introduce slack variables for first two inequalities
      // Introduce auxiliary variables for last equality, changing sign of RHS if necessary
      // 
      // Need to minimise sum of auxiliary variables -> add sum of last rows 
      // to get value function.
      LinearAlgebra::matrix<Rational> T(3*n+1,2*n+1);

      const Geometry::Point<R>& a=r.lower_corner();
      const Geometry::Point<R>& b=r.upper_corner();
      const Geometry::Point<R>& c=this->centre();
      const LinearAlgebra::matrix<R>& A=this->generators();

      for(size_type i=0; i!=n; ++i) {
        T(i,i)=1;
        T(i,2*n)=Rational(b[i])-Rational(a[i]);
        T(n+i,n+i)=1;
        T(n+i,2*n)=2;

        // Compute rhs = c[i]-a[i]-(A*1)[i]
        Rational rhs=Rational(c[i]) - Rational(a[i]);
        for(size_type j=0; j!=n; ++j) {
          rhs-=Rational(A(i,j));
        }
        
        if(rhs>=0) {
          T(2*n+i,i)=1;
          for(size_type j=0; j!=n; ++j) {
            T(2*n+i,n+j)=-A(i,j);
          }
          T(2*n+i,2*n)=rhs;
        }
        else {
          T(2*n+i,i)=-1;
          for(size_type j=0; j!=n; ++j) {
            T(2*n+i,n+j)=A(i,j);
          }
          T(2*n+i,2*n)=-rhs;
        }
        for(size_type j=0; j!=2u*n; ++j) {
          T(3*n,j)-=T(2*n+i,j);
        }
        T(3*n,2*n)-=T(2*n+i,2*n);
      }
      
      LinearAlgebra::LinearProgram<Rational> lp(T);
      
      bool result=(lp.optimal_value()!=0);
      /*
      if(result!=Geometry::disjoint(Polyhedron<R>(*this), Polyhedron<R>(r))) {
        std::cerr << "Incorrect result for \n  " << r << "\nand\n" << *this << "\n";
        std::cerr << T << "\n" << lp.tableau() << "\n";
        std::cerr << convert_to<double>(lp.tableau()(3*n,2*n)) << "\n";
        assert(false);
      }
      */
      return result;
    }
    
    template<typename R>
    bool
    Parallelotope<R>::disjoint(const Parallelotope<R>& p) const
    {
     return Geometry::disjoint(Polyhedron<R>(*this),Polyhedron<R>(p));
    }

    template <typename R>
    Point<R> 
    Parallelotope<R>::vertices(const size_t &i) const
    {
      Point<R> vertex(this->centre());
      const LinearAlgebra::matrix<R> gen=this->_generators;
      
      for (size_t j=0; j<gen.size1(); j++) {
	for (size_t k=0; k<gen.size2(); k++) {
          if ((1<<k)&(i)) {
	    vertex[j]+=gen(j,k);
	  } else {
	    vertex[j]-=gen(j,k);
	  }
	}
      }
 
      return vertex;

    }
    
    template<typename R>
    std::vector< Point<R> > 
    Parallelotope<R>::vertices() const
    {
      size_t vert_num=(1<<(this->_generators).size2());
      std::vector< Point<R> > vert(1<<(this->_generators).size2());

      assert((this->_generators).size2()<32);
      for (size_t i=0; i<vert_num; i++) {
	vert[i]=this->vertices(i);
      }

      return vert;
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
