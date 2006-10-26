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

#include "../geometry/point.h"
#include "../geometry/point_list.h"
#include "../geometry/rectangle.h"
#include "../geometry/list_set.h"
#include "../geometry/zonotope.h"
#include "../geometry/polyhedron.h"
#include "../geometry/polytope.h"

#include "parallelotope.h"

namespace Ariadne {
  namespace Geometry {

    template<class R>
    tribool
    Parallelotope<R>::_instantiate_geometry_operators() 
    {
      Rectangle<R> r;
      Parallelotope<R> p;
      tribool b=Geometry::subset(r,p);
      b = b&&b;
      return b;
    }
    
    template <class R>
    tribool 
    _parallelotope_contains_coordinates(const LinearAlgebra::Vector<R>& e) 
    {  
      R one=1;
      tribool result=true;
      for (size_t i=0; i<e.size(); i++) {
        R av=abs(e(i));
        if (av>one) { return false; }
        if (!(av<one)) { result=indeterminate; }
      }
      return true;
    }

    template <class R>
    void 
    Parallelotope<R>::_compute_generators_inverse() const 
    {  
      this->_generators_inverse=this->generators().inverse();
    }
    
      
    template <class R>
    void 
    Parallelotope< Interval<R> >::_compute_generators_inverse() const 
    {  
      this->_generators_inverse=this->generators().inverse();
    }
    
    template <class R>
    tribool 
    Parallelotope<R>::contains(const Point<R>& pt) const {
      if(this->_generators_inverse.number_of_rows()==0) {
        this->_compute_generators_inverse();
      }
      return _parallelotope_contains_coordinates(this->_generators_inverse*(pt-this->centre()));
    }
      
    template <class R>
    tribool 
    Parallelotope< Interval<R> >::contains(const Point<I>& pt) const {
      if(this->_generators_inverse.number_of_rows()==0) {
        this->_compute_generators_inverse();
      }
      return _parallelotope_contains_coordinates(this->_generators_inverse*(pt-this->centre()));
    }

    template <class R>
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
    
    template <class R>
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
      for(unsigned long k=0; k!=1u<<n; ++k) {
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
    
    
    template<class R>
    LinearAlgebra::Vector<typename traits<R>::arithmetic_type>
    Parallelotope<R>::coordinates(const Point<R>& s) const {
      LinearAlgebra::Vector<F> p=s.position_vector();
      LinearAlgebra::Vector<F> c=this->centre().position_vector();
      LinearAlgebra::Vector<F> d=p-c;
      LinearAlgebra::Matrix<F> G=this->generators();
      return G.solve(d);
    }


    template<class R>
    PointList<typename Parallelotope<R>::F>
    Parallelotope<R>::vertices() const
    {
      PointList<F> result;

      dimension_type d=this->dimension();
      assert(d<32);      
      Point<F> c(this->centre());
      LinearAlgebra::Matrix<F> g(this->generators());
      LinearAlgebra::Vector<F> e(d);

      size_type nv=(1<<d);
      result.reserve(nv);

      for (size_type i=0; i<nv; ++i) {
        for(size_type j=0; j!=d; ++j) {
          e[j]=(i&(1<<d) ? 1 : -1);
        }
        result.push_back(c+LinearAlgebra::Vector<F>(g*e));
      }
      return result;
    }
    
    

    
    template<class R>
    Parallelotope<R>
    Parallelotope<R>::over_approximation(const Point<I> &c, const LinearAlgebra::Matrix<I>& A)
    {
#ifdef DEBUG
      std::cerr << "IntervalParallelotope<R>::over_approximating_parallelotope() const" << std::endl;
#endif
      size_type n=c.dimension();
      
      Point<R> cmid=approximation(c);
      LinearAlgebra::Matrix<R> Amid=LinearAlgebra::approximate_value(A);
      
      LinearAlgebra::Matrix<R> D(n,n);
      for(size_type i=0; i!=n; ++i) {
        D(i,i)=c[i].radius();
      }
      
      LinearAlgebra::Matrix<I> Ainv=LinearAlgebra::inverse(LinearAlgebra::Matrix<I>(Amid));
      R err = ((Ainv*D).norm()+(Ainv*A).norm()).upper();

      for(size_type i=0; i!=n; ++i) {
        for(size_type j=0; j!=n; ++j) {
          Amid(i,j)=mul_up(Amid(i,j),err);
        }
      }
      
      return Geometry::Parallelotope<R>(cmid,Amid);
    }
    
    template<class R>
    Parallelotope<R>
    Parallelotope<R>::over_approximation(const Zonotope<R>& z)
    {
      dimension_type n=z.dimension();
      LinearAlgebra::Matrix<R> A(n,n);
      for(dimension_type i=0; i!=n; ++i) {
        for(dimension_type j=0; j!=n; ++j) {
          A(i,j)=z.generators()(i,j);
        }
      }
      LinearAlgebra::Matrix<F> Aq=A;
      LinearAlgebra::Matrix<F> B=LinearAlgebra::inverse(Aq);
      const Point<R>& c=z.centre();
      
      throw std::runtime_error("Zonotope<R>::over_approximating_parallelotope() const not implemented");
      
      /*
      typedef Rational FF;
      Parallelotope<R> p(c,A);
      while(!Geometry::subset(Polytope<FF>(z),Polyhedron<FF>(p))) {
        Aq*=F(2);
        p=Parallelotope<R>(c,A);
      }
      return p;
      */
    }        
    
    template<class R>
    Parallelotope<R>
    Parallelotope< Interval<R> >::over_approximation() const
    { 
      throw std::runtime_error("Parallelotope<Interval<R>>::over_approximation() const not implemented");
    }
    
/*
    template<class R>
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
    
    /*! \brief Tests if the parallelotope contains \a pt. */
    template <class R>
    tribool 
    subset(const Rectangle<R>& r, const Parallelotope<R>& p)
    {
      //std::cerr << "subset(const Rectangle<R>& r, const Parallelotope<R>& p)" << std::endl;
      tribool result=true;
      for(typename Rectangle<R>::vertices_const_iterator v_iter=r.vertices_begin();
          v_iter!=r.vertices_end(); ++v_iter) 
      {
        result=result && p.contains(*v_iter);
        if(result==false) { 
          //std::cerr << "p.contains(" << *v_iter << ")==false" << std::endl;
          return result; 
        }
      }
      return result;
    }



    template <class R>
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
    
    template <class R>
    std::istream& 
    Parallelotope<R>::read(std::istream& is)
    {
      throw std::domain_error("Not implemented");
    }
      
  }
}
