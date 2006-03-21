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
//#include "../linear_algebra/constraint.h"
#include "../linear_algebra/linear_program.h"

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
      for(size_t i=0; i!=this->dimension(); ++i) {
        for(size_t j=0; j!=this->dimension(); ++j) {
          offset[i] += abs(this->_generators(i,j));
        }
      }
      return Rectangle<R>(this->centre()-offset, this->centre()+offset);
    }
      
    
    template <typename R>
    Parallelotope<R>::operator Polyhedron<R>() const 
    {
      std::cerr << "Parallelotope<" << name<R>() << "<::operator Polyhedron<" << name<R>() << ">() const" << std::endl;
      using namespace ::Ariadne::LinearAlgebra;
      
      typedef typename Parallelotope<R>::Real Real;
      typedef typename Parallelotope<R>::Point Point;
      
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
    bool 
    Parallelotope<R>::contains(const Point& point) const 
    {
      if (point.dimension()!=this->dimension()) {
        throw std::domain_error("This object and parameter have different space dimensions");
      }  
        
      Vector c=coordinates(point);
      for(dimension_type i=0; i!=this->dimension(); ++i) {
        if(c[i]<-1 || c[i]>1) {
          return false;
        }
      }
      return true;
    }
    
    template <typename R>
    bool 
    Parallelotope<R>::interior_contains(const Point& point) const {
      if (point.dimension()!=this->dimension()) {
        throw std::domain_error("This object and parameter have different space dimensions");
      }  
        
      Vector c=coordinates(point);
      for(dimension_type i=0; i!=this->dimension(); ++i) {
        if(c[i]<=-1 || c[i]>=1) {
          return false;
        }
      }
      return true;
    }
      
      
      
    template <typename R>
    ListSet<R,Parallelotope>
    Parallelotope<R>::subdivide() const 
    {
      size_type n=this->dimension();
      ListSet<R,Geometry::Parallelotope> result(this->dimension());
      Matrix new_generators=this->generators()/2;
      
      Point first_centre=this->centre();
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
        Point new_centre=first_centre;
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
      Vector c=(this->centre()).position_vector();
      
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

      Vector c = this->centre() - Point(n,0);
      o = A * c;
      
      b.resize(n);
      for(size_type i=0; i!=n; ++i) {
        b(i)=multipliers(i);
      }
    }
      
  
  
    template<typename R>
    typename Parallelotope<R>::Vector 
    Parallelotope<R>::coordinates(const Point& s) const {
      Vector diff = s-this->_centre;
      Matrix inv = LinearAlgebra::inverse(this->_generators);
      return prod(inv,diff);
    }

    template<typename R>
    bool
    Parallelotope<R>::disjoint(const Rectangle<R>& r) const
    {
      //std::cerr << "Parallelotope<" << name<R>() << ">::disjoint(const Rectangle<" << name<R>() << ">& r) const" << std::endl;
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
          rhs-=A(i,j);
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
