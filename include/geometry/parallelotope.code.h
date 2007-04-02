/***************************************************************************
 *            parallelotope.code.h
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
#include "../linear_algebra/qr_matrix.h"
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
      //Parallelotope<R>(*goa)(const Parallelotope< Interval<R> >&) = &Geometry::over_approximation<R>;
      Rectangle<R> r;
      Parallelotope<R> p;
      Geometry::subset(r,p);
      return false;
    }
    
    template<class R>
    tribool
    Parallelotope< Numeric::Interval<R> >::_instantiate_geometry_operators() 
    {
      Parallelotope<R> p;
      Parallelotope< Numeric::Interval<R> > ip;
      Zonotope< Numeric::Interval<R> > iz;
      p=Geometry::over_approximation(ip);
      p=Geometry::orthogonal_over_approximation(iz);
      return p.empty();
    }
    
    template<class R>
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

    template<class R>
    R
    Parallelotope<R>::volume() const
    {
      return approximate_value(this->generators().determinant());
    }
   
    
    
    template<class R>
    void 
    Parallelotope<R>::_compute_generators_inverse() const 
    {  
      this->_generators_inverse=this->generators().inverse();
    }
    
      
    template<class R>
    void 
    Parallelotope< Numeric::Interval<R> >::_compute_generators_inverse() const 
    {  
      this->_generators_inverse=this->generators().inverse();
    }
    
    template<class R>
    tribool 
    Parallelotope<R>::contains(const Point<R>& pt) const {
      if(this->_generators_inverse.number_of_rows()==0) {
        this->_compute_generators_inverse();
      }
      return _parallelotope_contains_coordinates(this->_generators_inverse*(pt-this->centre()));
    }
      
    template<class R>
    tribool 
    Parallelotope< Numeric::Interval<R> >::contains(const Point<I>& pt) const {
      if(this->_generators_inverse.number_of_rows()==0) {
        this->_compute_generators_inverse();
      }
      return _parallelotope_contains_coordinates(this->_generators_inverse*(pt-this->centre()));
    }

    template<class R>
    ListSet< Parallelotope<R> >
    Parallelotope<R>::divide() const 
    {
      size_type n=this->dimension();
      ListSet< Parallelotope<R> > result(this->dimension());
      
      LinearAlgebra::Matrix<R> new_generators=this->generators();
      
      R two=2;
      R max_norm=0;
      size_type max_column=0;
      for(size_type j=0; j!=n; ++j) {
        R norm = LinearAlgebra::norm(LinearAlgebra::Vector<R>(new_generators.column(j)));
        if(norm>max_norm) {
          max_norm=norm;
          max_column=j;
        }
      }
      
      size_type j=max_column;
      for(size_type i=0; i!=n; ++i) {
        new_generators(i,j)=div_up(new_generators(i,j),two);
      }
      
      Point<R> new_centre=sub_approx(this->centre(),div_approx(LinearAlgebra::Vector<R>(new_generators.column(j)),two));
      result.adjoin(Parallelotope<R>(new_centre,new_generators));
      new_centre=add_approx(new_centre,LinearAlgebra::Vector<R>(new_generators.column(j)));
      result.adjoin(Parallelotope(new_centre,new_generators));

      return result;
    }
    
    template<class R>
    ListSet< Parallelotope<R> >
    Parallelotope<R>::subdivide() const 
    {
      ListSet< Parallelotope<R> > result(this->dimension());
      
      R two=2;
      size_type n=this->dimension();
      
      LinearAlgebra::Matrix<R> new_generators(n,n);
      for(size_type i=0; i!=n; ++i) {
        for(size_type j=0; j!=n; ++j) {
          new_generators(i,j)=div_up(this->generators()(i,j),two);
        }
      }
      
      Point<R> first_centre=this->centre();
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
    LinearAlgebra::Vector<typename Numeric::traits<R>::arithmetic_type>
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
          e(j)=(i&(1<<d) ? 1 : -1);
        }
        result.push_back(c+LinearAlgebra::Vector<F>(g*e));
      }
      return result;
    }
    
    
    template<class R>
    Parallelotope<R>
    over_approximation(const Parallelotope<R>& p)
    {
      return p;
    }
    
    template<class R>
    Parallelotope<R>
    over_approximation(const Parallelotope< Numeric::Interval<R> >& p)
    {
      //std::cerr << "Numeric::IntervalParallelotope<R>::over_approximating_parallelotope() const" << std::endl;
      typedef Numeric::Interval<R> I;
      const Point<I>& c=p.centre();
      const LinearAlgebra::Matrix<I> G=p.generators();
      
      size_type n=c.dimension();
      
      Point<R> cmid=approximate_value(c);
      LinearAlgebra::Matrix<R> Gmid=LinearAlgebra::approximate_value(G);
      
      LinearAlgebra::Matrix<R> D(n,n);
      for(size_type i=0; i!=n; ++i) {
        D(i,i)=c[i].radius();
      }
      
      LinearAlgebra::Matrix<I> Ginv=LinearAlgebra::inverse(LinearAlgebra::Matrix<I>(Gmid));
      R err = ((Ginv*D).norm()+(Ginv*G).norm()).upper();

      for(size_type i=0; i!=n; ++i) {
        for(size_type j=0; j!=n; ++j) {
          Gmid(i,j)=mul_up(Gmid(i,j),err);
        }
      }
      
      // Check to make such result is valid.
      assert((bool)(LinearAlgebra::norm(LinearAlgebra::row_norms(Gmid.inverse()*G)+(cmid-c))<=R(1)));
      return Geometry::Parallelotope<R>(cmid,Gmid);
    }
    
 
    template<class R>
    Parallelotope<R>
    orthogonal_over_approximation(const Zonotope<R>& z)
    {
      return orthogonal_over_approximation(Zonotope< Numeric::Interval<R> >(z));
    }
    
    template<class R>
    Parallelotope<R>
    orthogonal_over_approximation(const Zonotope< Numeric::Interval<R> >& z)
    {
      //std::cerr << "Parallelotope<R>::orthogonal_over_approximation(const Zonotope<I>&) const" << std::endl;
      typedef Numeric::Interval<R> I;
      typedef typename Numeric::traits<R>::approximate_arithmetic_type A;
      const Point<I>& c=z.centre();
      const LinearAlgebra::Matrix<I>& G=z.generators();
      
      size_type d=z.dimension();
      size_type ng=z.number_of_generators();
      
      Point<R> cmid=approximate_value(c);
      LinearAlgebra::Vector<I> cerr=c-cmid;
      
      LinearAlgebra::Matrix<A> Gapprx(d,ng);
      for(size_type i=0; i!=d; ++i) {
        for(size_type j=0; j!=ng; ++j) {
          Gapprx(i,j)=Numeric::conv_approx<A>(approximate_value(G(i,j)));
        }
      }
      
      LinearAlgebra::QRMatrix<I> QR(G);
      LinearAlgebra::Matrix<I> Q=QR.Q();
      LinearAlgebra::Matrix<R> Qmid=approximate_value(Q);
      LinearAlgebra::Vector<I> Rrwnrm=LinearAlgebra::row_norms(Qmid.transpose()*G);
      for(size_type i=0; i!=d; ++i) {
        R scale=(Rrwnrm(i)+cerr(i)).upper();
        for(size_type j=0; j!=d; ++j) {
          Qmid(i,j)=mul_up(Qmid(i,j),scale);
        }
      }
      
      
      // Check to make such result is valid.
      //assert(LinearAlgebra::norm(LinearAlgebra::row_norms(Qmid.inverse()*A)+(cmid-c))<=R(1));
      return Geometry::Parallelotope<R>(cmid,Qmid);
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
    template<class R>
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



    template<class R>
    std::ostream&
    Parallelotope<R>::write(std::ostream& os) const
    {
      const Parallelotope<R>& p=*this;
      if(p.dimension() > 0) {
        os << "Parallelotope( centre=" << p.centre()
           << ", directions=" << p.generators()
           << " ) ";
      }

      return os;
    }
    
    template<class R>
    std::ostream&
    Parallelotope< Numeric::Interval<R> >::write(std::ostream& os) const
    {
      const Parallelotope< Numeric::Interval<R> >& p=*this;
      if(p.dimension() > 0) {
        os << "Parallelotope( centre=" << p.centre()
           << ", directions=" << p.generators()
           << " ) ";
      }

      return os;
    }
    
    template<class R>
    std::istream& 
    Parallelotope<R>::read(std::istream& is)
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }
      
  }
}
