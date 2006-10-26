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

#include "../combinatoric/binary_word.h" 
#include "../base/array.h" 
#include "../geometry/point.h" 
#include "../geometry/point_list.h" 
#include "../geometry/list_set.h" 
#include "../geometry/grid_set.h" 

namespace Ariadne {
  namespace Geometry {

    template <class R>
    Rectangle<R>::Rectangle(const std::string& s)
      : _bounds()
    {
      std::stringstream ss(s);
      ss >> *this;
    }

    template <class R>
    Rectangle<R>& 
    Rectangle<R>::expand_by(const real_type& delta) 
    {
      Interval<R> expand(-delta,delta);
      for (size_type j=0; j< this->dimension(); ++j) {
        (*this)[j]+=expand;
      }
        
      return *this;
    }
      
    
    template <class R>
    Rectangle<R>
    Rectangle<R>::quadrant(const Combinatoric::BinaryWord& w) const 
    {
      assert(w.size() == this->dimension());
      Rectangle<R> quadrant(this->dimension());
      
      for (size_type i=0; i!=this->dimension(); ++i) {
        if(w[i]) {
          quadrant[i]=Interval<R>(this->interval(i).centre(),this->upper_bound(i));
        } 
        else {
          quadrant[i]=Interval<R>(this->lower_bound(i),this->interval(i).centre());
        }
      }
      return quadrant;
    }
      
    template <class R>
    ListSet<R,Rectangle>
    Rectangle<R>::subdivide() const 
    {
      size_type n=this->dimension();
      ListSet<R,Geometry::Rectangle> result(this->dimension());
      
      Point<R> lwr_crnr=this->lower_corner();
      Point<R> cntr=this->centre();
      Point<R> upr_crnr=this->upper_corner();
      Point<R> new_lwr_crnr;
      Point<R> new_upr_crnr;
      for(size_type i=0; i!=1u<<n; ++i) {
        for(size_type j=0; j!=n; ++j) {
          if(i&(1<<j)) {
            new_lwr_crnr[j]=lwr_crnr[j];
            new_upr_crnr[j]=cntr[j];
          }
          else {
            new_lwr_crnr[j]=cntr[j];
            new_upr_crnr[j]=upr_crnr[j];
          }
        }
        Rectangle<R> new_rect(new_lwr_crnr,new_upr_crnr);
        result.adjoin(new_rect);
      }
      return result;
    }


    template<class R>
    PointList<R>
    Rectangle<R>::vertices() const
    {
      size_type number_of_vertices=(1<<this->dimension());
      PointList<R> result(this->dimension(),number_of_vertices);
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
    
    template <class R>
    Point<R> 
    Rectangle<R>::vertex(size_type i) const 
    {
      size_type d=this->dimension();
      state_type result(d); 
            
      if (i >= (size_type)(1<<d)) {
         throw std::domain_error("Rectangle::vertex(i): Wrong vertex index.");
      }
      
      for (size_type j=0; j<d; ++j) {
        if (i%2) {
          result[j]=this->lower_bound(j);
        }
        else {
          result[j]=this->upper_bound(j);
        }
        i=i/2;
      }

      return result;
    }
    
    template<class R>
    typename Rectangle<R>::vertices_const_iterator
    Rectangle<R>::vertices_begin() const
    {
      return RectangleVerticesIterator<R>(*this,false);
    }
    
    template<class R>
    typename Rectangle<R>::vertices_const_iterator
    Rectangle<R>::vertices_end() const
    {
      return RectangleVerticesIterator<R>(*this,true);
    }
    
    
    template <class R>
    tribool 
    subset(const Rectangle<R>& A, 
           const ListSet<R,Geometry::Rectangle>& B)
    {
      return Geometry::subset(A, GridMaskSet<R>(B));
    }
    
    
    template <class R>
    std::ostream&
    Rectangle<R>::write(std::ostream& os) const 
    {
      const Rectangle<R>& self=*this;
      if(self.dimension()==0) {
        os << "EmptyRectangle";
      }
      else {
        os << "[" << self.lower_bound(0) << "," << self.upper_bound(0) << "]";
        for(dimension_type i=1; i!=self.dimension(); ++i) {
          os << "x[" << self.lower_bound(i) << "," << self.upper_bound(i) << "]";
        }
      }
      return os;
    }
    
    template <class R>
    std::istream& 
    Rectangle<R>::read(std::istream& is)
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
        (*this)=Rectangle<R>(v.begin(),v.end());
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
