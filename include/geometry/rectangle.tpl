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
/*
#include <vector>
#include <list>
#include <set>
#include <valarray>
*/

#include "rectangle.h"

#include "../base/binary_word.h" 
#include "../base/array.h" 
#include "../linear_algebra/interval_vector.h" 
#include "../geometry/lattice_set.h" 
#include "../geometry/list_set.h" 
#include "../geometry/grid_set.h" 
#include "../geometry/parallelotope.h" 
#include "../geometry/zonotope.h" 
#include "../geometry/polyhedron.h" 

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
    Rectangle<R> 
    minkowski_sum(const Rectangle<R>& A, const Rectangle<R>& B)
    {
      using namespace Ariadne::LinearAlgebra;

      if (A.dimension()!=B.dimension())
        throw std::domain_error(
              "minkowski_sum: the two rectangles have different dimension.");

      Point<R> lower(A.dimension()), upper(A.dimension());

      for(dimension_type i=0; i<A.dimension(); i++) {
        lower[i]=A.lower_bound(i)+B.lower_bound(i);
        upper[i]=A.upper_bound(i)+B.upper_bound(i);
      }

      return Rectangle<R>(lower,upper);

    }

    template<typename R> 
    Rectangle<R> 
    minkowski_difference(const Rectangle<R>& A, const Rectangle<R>& B)
    {
      using namespace Ariadne::LinearAlgebra;

      if (A.dimension()!=B.dimension())
        throw std::domain_error(
              "minkowski_difference: the two rectangles have different dimension.");


      Point<R> lower(A.dimension()), upper(A.dimension());

      for(dimension_type i=0; i<A.dimension(); i++) {
        lower[i]=A.lower_bound(i)-B.lower_bound(i);
        upper[i]=A.upper_bound(i)-B.upper_bound(i);

      if ( lower[i]>upper[i])
        return Rectangle<R>(A.dimension());
      }

      return Rectangle<R>(lower,upper);
    }

    template <typename R>
    Rectangle<R>::operator Parallelotope<R>() const 
    {
       return Parallelotope<R>(*this); 	
    }
    
    template <typename R>
    Rectangle<R>::operator Zonotope<R>() const 
    {
       return Zonotope<R>(*this);
    }
    
    template <typename R>
    Rectangle<R>::operator Polyhedron<R>() const 
    {
       return Polyhedron<R>(*this);
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
    array< Point<R> >
    Rectangle<R>::vertices() const
    {
      size_type number_of_vertices=(1<<this->dimension());
      array< Point<R> > result(number_of_vertices);
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
    
    /* Compute all points in A on the grid of vertices of rectangles in the cover */
    template <typename R>
    void
    compute_gridpoints(std::vector< std::set<R> >& gridpoints,
                       const Rectangle<R>& A, 
                       const ListSet< R,Rectangle >& cover)
    {
      typedef typename ListSet< R, Rectangle >::const_iterator list_iterator;

      size_type dimension = A.dimension();
      
      for(size_type i=0; i!=dimension; ++i) {
        R lower=A.lower_bound(i);
        R upper=A.upper_bound(i);
        gridpoints[i].insert(lower);
        gridpoints[i].insert(upper);
        for(list_iterator rect=cover.begin(); rect!=cover.end(); ++rect) {
          R bound = rect->lower_bound(i);
          if(lower<bound && bound<upper) {
            gridpoints[i].insert(bound);
          }
          bound = rect->upper_bound(i);
          if(lower<bound && bound<upper) {
            gridpoints[i].insert(bound);
          }
        }
      }
      
#ifdef DEBUG
      std::cerr << "Gridpoints: " << gridpoints << '\n';
#endif
    }
    
    /*! \brief Tests inclusion in an open cover.
     * FIXME: Use some kind of GridSet for this operation.
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
    
      typedef typename ListSet<R,Geometry::Rectangle>::const_iterator list_iterator;
      typedef std::valarray<size_type> index_type;

      const Rectangle<R>& A=*this;
      size_type dimension = A.dimension();
      
      std::vector< std::set<R> > gridpoints(dimension);
      compute_gridpoints(gridpoints, A, B);
      
      
      /* Whether the jth gridpoint (in some ordering) is covered */
      std::vector<bool> cover_flags;
      
      /* Strides for indexing */
      index_type strides(dimension);
      size_type stride=1;
      for(size_type i=0; i!=dimension; ++i) {
        strides[i]=stride;
        stride*=gridpoints[i].size();
      }
      cover_flags.resize(stride);
      
      std::vector<size_type> lower_indices(dimension+1);
      std::vector<size_type> upper_indices(dimension+1);
      lower_indices[dimension] = 0;
      upper_indices[dimension] = 2;
      
      for(list_iterator rect=B.begin(); rect!=B.end(); ++rect) {
        for(size_type i=0; i!=dimension; ++i) {
          R lower_bnd=rect->lower_bound(i);
          size_type lower_indx=0;
          typename std::set<R>::const_iterator iter=gridpoints[i].begin();
          while(iter!=gridpoints[i].end() && (*iter) <= lower_bnd) {
            ++iter;
            ++lower_indx;
          }
          lower_indices[i]=lower_indx;
          
          R upper_bnd=rect->upper_bound(i);
          size_type upper_indx=0;
          iter=gridpoints[i].begin();
          while(iter!=gridpoints[i].end() && (*iter) < upper_bnd) {
            ++iter;
            ++upper_indx;
          }
          upper_indices[i]=upper_indx;
        }
        
#ifdef DEBUG
        std::cerr << "Rectangle: " <<  (*rect)
                  << " lower_indices: " << lower_indices
                  << " upper_indices: " << upper_indices << '\n';
#endif
        
        index_type index(dimension+1);
        for(size_type i=0; i!=dimension; ++i) {
          index[i] = lower_indices[i];
        }
        index[dimension]=0;
        
        while(index[dimension] != 1) {
#ifdef DEBUG
          std::cerr << index << " " << upper_indices << "\n";
#endif
          
          size_type entry = 0;
          for(size_type j=0; j!=dimension; ++j) {
            entry += index[j]*strides[j];
          }
          cover_flags[entry] = true;
          
          size_type inc=0;
          ++(index[inc]);
          while(index[inc] == upper_indices[inc]) {
            index[inc]=0;
            ++inc;
            ++(index[inc]);
          }
        }
        
#ifdef DEBUG
        std::cerr << index << "\n\n";
        std::cerr << cover_flags << '\n';
#endif
      }
      
      for( std::vector<bool>::const_iterator flag = cover_flags.begin();
           flag != cover_flags.end(); ++flag) {
        if(*flag == false) {
          return false;
        }
      }
      
      return true;
    }
    
    /*! \brief Tests inclusion. 
     *  FIXME: Convert B to GridMaskSet<R> .
     */
    template <typename R>
    bool 
    Rectangle<R>::subset(const ListSet<R,Geometry::Rectangle>& B) const
    {
      return Geometry::subset(*this, GridMaskSet<R>(B));
      typedef typename ListSet<R,Geometry::Rectangle>::const_iterator list_iterator;
      
      typedef std::valarray<size_type> index_type;

      const Rectangle<R>& A=*this;
      size_type dimension = A.dimension();
      
      std::vector< std::set<R> > gridpoints(dimension);
      compute_gridpoints(gridpoints, A, B);
      
      
      /* Whether the jth gridpoint (in some ordering) is covered */
      std::vector<bool> cover_flags;
      
      /* Strides for indexing */
      index_type strides(dimension);
      size_type stride=1;
      for(size_type i=0; i!=dimension; ++i) {
        strides[i]=stride;
        stride*=gridpoints[i].size();
      }
      cover_flags.resize(stride);
      
      std::vector<size_type> lower_indices(dimension+1);
      std::vector<size_type> upper_indices(dimension+1);
      lower_indices[dimension] = 0;
      upper_indices[dimension] = 2;
      
      for(list_iterator rect=B.begin(); rect!=B.end(); ++rect) {
        for(size_type i=0; i!=dimension; ++i) {
          R lower_bound=rect->lower_bound(i);
          size_type lower_index=0;
          typename std::set<R>::const_iterator iter=gridpoints[i].begin();
          while(iter!=gridpoints[i].end() && (*iter) < lower_bound) {
            ++iter;
            ++lower_index;
          }
          lower_indices[i]=lower_index;
          
          R upper_bound=rect->upper_bound(i);
          size_type upper_index=0;
          iter=gridpoints[i].begin();
          while(iter!=gridpoints[i].end() && (*iter) <= upper_bound) {
            ++iter;
            ++upper_index;
          }
          upper_indices[i]=upper_index;
        }
        
#ifdef DEBUG
        std::cerr << "Rectangle: " <<  (*rect)
                  << " lower_indices: " << lower_indices
                  << " upper_indices: " << upper_indices << '\n';
#endif
        
        index_type index(dimension+1);
        for(size_type i=0; i!=dimension; ++i) {
          index[i] = lower_indices[i];
        }
        index[dimension]=0;
        
        while(index[dimension] != 1) {
#ifdef DEBUG
          std::cerr << index << " " << upper_indices << "\n";
#endif
          
          size_type entry = 0;
          for(size_type j=0; j!=dimension; ++j) {
            entry += index[j]*strides[j];
          }
          cover_flags[entry] = true;
          
          size_type inc=0;
          ++(index[inc]);
          while(index[inc] == upper_indices[inc]) {
            index[inc]=0;
            ++inc;
            ++(index[inc]);
          }
        }
        
#ifdef DEBUG
        std::cerr << index << "\n\n";
        std::cerr << cover_flags << '\n';
#endif
        
      }
      
      for( std::vector<bool>::const_iterator flag = cover_flags.begin(); 
           flag != cover_flags.end(); ++flag) {
        if(*flag == false) {
          return false;
        }
      }
      
      return true;
    }
    

    template <typename R>
    std::ostream&
    operator<<(std::ostream& os, const Rectangle<R>& r) 
    {
      
      /*
        os << "{ lower=" << (r._lower_corner) << ", " ;
        os << "upper=" << (r._upper_corner) << " }" ;
      */
      
      /*if(r.empty()) {
        os << "Empty";
      }
      else if(r.dimension() > 0) {
        os << r[0];
        for(size_type i=1; i!=r.dimension(); ++i) {
          os << "x" << r[i];
        }
      }
      */
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
        /* Representation as list of intervals (deprecated) */ 
      /*
        std::vector< Interval > v;
        is >> v;
        r=Rectangle<R>(v.size(),&v[0]);
      */
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
