/***************************************************************************
 *            rectangle.tpl
 *
 *  Mon 2 May 2005
 *  Copyright 2005  Alberto Casagrande, Pieter Collins
 *  Email casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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
 
#include <iostream>
#include <string>
#include <sstream>
#include <exception>

#include "rectangle.h"

#include "../base/binary_word.h" 
#include "../base/array.h" 
#include "../linear_algebra/interval_vector.h" 
#include "../geometry/polyhedron.h" 
#include "../geometry/lattice_set.h" 
#include "../geometry/list_set.h" 
#include "../geometry/grid_set.h" 

namespace Ariadne {
  namespace Geometry {

    template <typename R>
    Rectangle<R>::Rectangle(const std::string& s)
      : _lower_corner(1), _upper_corner(1)
    {
      std::stringstream ss(s);
      ss >> *this;
    }

    template<typename R>
    Rectangle<R>::operator Polyhedron<Rational> () const {
      LinearAlgebra::IntervalVector<Rational> rpv(this->dimension());
      for(dimension_type i=0; i!=rpv.size(); ++i) {
        rpv[i]=Interval<Rational>(Rational(this->lower_bound(i)),Rational(this->upper_bound(i)));
      }
      return Polyhedron<Rational>(rpv);
    }
    
    template<typename R> 
    Rectangle<R> 
    minkowski_sum(const Rectangle<R>& A, const Rectangle<R>& B)
    {
      if (A.dimension()!=B.dimension()) {
        throw std::domain_error(
              "minkowski_sum: the two rectangles have different dimension.");
      }
      
      Rectangle<R> result(A.dimension());
      for(dimension_type i=0; i!=A.dimension(); ++i) {
        result.set_interval(i,A[i]+B[i]);
      }

      return result;

    }

    template<typename R> 
    Rectangle<R> 
    minkowski_difference(const Rectangle<R>& A, const Rectangle<R>& B)
    {
      if (A.dimension()!=B.dimension()) {
        throw std::domain_error(
              "minkowski_difference: the two rectangles have different dimension.");
      }

      Rectangle<R> result(A.dimension());
      for(dimension_type i=0; i!=A.dimension(); i++) {
          result.set_interval(i,A[i]-B[i]);
      }
      return result;
    }
    

    template <typename R>
    LinearAlgebra::IntervalVector<R> 
    Rectangle<R>::position_vectors() const 
    {
      LinearAlgebra::IntervalVector<R> result(this->dimension());
      for(size_type i=0; i!=result.size(); ++i) {
        result(i)=this->operator[](i);
      }
      return result;
    }
      
    template <typename R>
    Rectangle<R>& 
    Rectangle<R>::expand_by(const real_type& delta) 
    {
      for (size_type j=0; j< this->dimension(); ++j) {
      this->_upper_corner[j]+=delta;
        this->_lower_corner[j]-=delta;
      }
      return *this;
    }
      
    template <typename R>
    Rectangle<R>
    Rectangle<R>::quadrant(const BinaryWord& q) const 
    {
      assert(q.size() == this->dimension());
      Rectangle<R> quadrant(this->dimension());
      
      for (size_type j=0; j< this->dimension(); j++) {
        if (q[j]) {
          quadrant._lower_corner[j]=(this->_upper_corner[j]+
                                      this->_lower_corner[j])/2;
          quadrant._upper_corner[j]=this->_upper_corner[j];
        } 
        else {
          quadrant._upper_corner[j]=(this->_upper_corner[j]+
                                      this->_lower_corner[j])/2;
          quadrant._lower_corner[j]=this->_lower_corner[j];
        }
      }
      return quadrant;
    }
      
    template <typename R>
    ListSet<R,Rectangle>
    Rectangle<R>::subdivide() const 
    {
      size_type n=this->dimension();
      ListSet<R,Geometry::Rectangle> result(this->dimension());
      
      array<index_type> lower(n,0);
      array<index_type> upper(n,2);
      array<index_type> finish(n,0);
      finish[n-1]=2;
      lattice_iterator end(finish,lower,upper);

      Point<R> lwr_crnr=this->lower_corner();
      LinearAlgebra::Vector<R> offst=(this->upper_corner()-this->lower_corner())/2;
      for(lattice_iterator iter(lower,lower,upper); iter!=end; ++iter) {
        array<index_type> ary=*iter;
        state_type new_lwr_crnr=lwr_crnr;
        for(size_type i=0; i!=n; ++i) {
          if(ary[i]==1) {
            new_lwr_crnr[i]+=offst(i);
          }
        }
        result.adjoin(Rectangle(new_lwr_crnr,new_lwr_crnr+offst));
      }
      return result;
    }


    template<typename R>
    std::vector< Point<R> >
    Rectangle<R>::vertices() const
    {
      size_type number_of_vertices=(1<<this->dimension());
      std::vector< Point<R> > result(number_of_vertices);
      Point<R> vertex(this->dimension());
      
      for (size_type i=0; i<number_of_vertices; ++i) {
        for (size_type j=0; j<this->dimension(); ++j) {
          if ((1<<j)&(i)) {
            vertex[j]=this->upper_bound(j);
          } 
          else {
            vertex[j]=this->lower_bound(j);
          }
        }
        result[i]=vertex;
      }
      return result;   
    }
    
    
    /*! \brief Tests inclusion in an open cover.
     */
    template <typename R>
    bool 
    Rectangle<R>::subset_of_open_cover(const ListSet<R,Geometry::Rectangle>& cover) const
    {
      throw std::domain_error("subset_of_open_cover(Simplex, std::vector<Simplex>) not implemented");
    }
    
    /*! \brief Tests inclusion in the interior of a denotable set.
     * FIXME: Use some kind of GridSet for this operation.
     * WARNING: Maybe this is the wrong routing...
     */
    template <typename R>
    bool 
    Rectangle<R>::inner_subset(const ListSet<R,Geometry::Rectangle>& B) const
    {
      return Geometry::inner_subset(*this, GridMaskSet<R>(B));
    }
    
    /*! \brief Tests inclusion. 
     *  FIXME: Convert B to GridMaskSet<R> .
     */
    template <typename R>
    bool 
    Rectangle<R>::subset(const ListSet<R,Geometry::Rectangle>& B) const
    {
      return Geometry::subset(*this, GridMaskSet<R>(B));
    }
    
    template <typename R>
    std::ostream&
    operator<<(std::ostream& os, const Rectangle<R>& r) 
    {
      if(r.dimension() > 0) {
        os << "["<<r._lower_corner[0]<<","<<r._upper_corner[0]<<"]";
        for(size_type i=1; i!=r.dimension(); ++i) {
          os << "x["<<r._lower_corner[i]<<","<<r._upper_corner[i]<<"]";
        }
      }
      else {
        os << "ZeroDimensionalRectangle";
      }
      return os;
    }
    
    template <typename R>
    std::istream& 
    operator>>(std::istream& is, Rectangle<R>& r)
    {
     
      char c;
      is >> c;
      is.putback(c);
      if(c=='[') {
        /* Representation as a literal [a1,b1]x[a2,b2]x...x[an,bn] */
        std::vector< Interval<R> > v;
        Interval<R> i;
        c='x';
        while(c=='x') {
          is >> i;
          v.push_back(i);
          c=' ';
          while( is && c==' ') {
            is >> c;
          }
        }
        if(is) {
          is.putback(c);
        }
        r=Rectangle<R>(v.begin(),v.end());
      }
      else {
        /* representation as lower and upper corners */
        /* FIXME */
        // throw invalid_input("Not implemented");
      }
      return is;
    }
    
  }
}
