/***************************************************************************
 *            grid.code.h
 *
 *  Copyright  2005-6  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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
 *  Foundation, Inc., 59 Templece Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "grid.h"

#include <ostream>
#include <iomanip>
#include <algorithm>

#include "exceptions.h"
#include "../base/stlio.h"

#include "../numeric/arithmetic.h"

#include "../geometry/rectangle.h"
#include "../geometry/list_set.h"

#include "../geometry/grid_set.h" // for GridBlock<R>

namespace Ariadne {
  namespace Geometry {


    template<class R> 
    Grid<R>::~Grid()
    {
    }
  
  

    template<class R> 
    Grid<R>::Grid(const dimension_type& d, const R& l)
      : _origin(d,static_cast<R>(0)),
        _lengths(d,l)
    {
      this->create();
    }
  

    template<class R> 
    Grid<R>::Grid(const LinearAlgebra::Vector<R>& v)
      : _origin(v.size(),static_cast<R>(0)),
        _lengths(v.data())
    {
      this->create();
    }
  

    template<class R> 
    Grid<R>::Grid(const Point<R>& pt, const LinearAlgebra::Vector<R>& v)
      : _origin(pt.data()),
        _lengths(v.data())
    {
      if(pt.dimension() != v.size()) {
        throw IncompatibleDimensions(__PRETTY_FUNCTION__);
      }
      this->create();
    }
  

    template<class R> 
    Grid<R>::Grid(const Rectangle<R>& r, const Combinatoric::LatticeBlock& lb)
      : _origin(r.dimension()),
        _lengths(r.dimension())
    {
      if(r.dimension() != lb.dimension()) {
        throw IncompatibleDimensions(__PRETTY_FUNCTION__);
      }
      for(dimension_type i=0; i!=r.dimension(); ++i) {
        const R& l=r.lower_bound(i);
        const R& u=r.upper_bound(i);
        const index_type& a=lb.lower_bound(i);
        const index_type& b=lb.upper_bound(i);
        R& o=this->_origin[i];
        R& s=this->_lengths[i];

        //std::cerr << "a=" << a << ", b=" << b << ";  l=" << l << ", u=" << u << ";  "  << std::flush;

        if(a==b) {
           throw std::runtime_error("Invalid grid block");
        }
        if(a>0) { // o < l
          s=div_up(sub_up(u,l),b-a);
          o=sub_down(l,mul_up(s,a));
          while(add_approx(o,mul_approx(s,b))<u) {
            s=next_up(s);
            o=sub_down(l,mul_up(s,a));
          }
        } else if(b<0) { // o > u
          s=div_up(sub_up(u,l),b-a);
          o=sub_up(u,mul_down(s,b));
          while(add_approx(o,mul_approx(s,a))>l) {
            s=next_up(s);
            o=sub_up(u,mul_down(s,b));
          }
        } else { // l <= o <= u
          s=div_up(sub_up(u,l),b-a);
          o=div_approx(sub_approx(mul_approx(l,b),mul_approx(u,a)),b-a);
          while(add_approx(o,mul_approx(s,a))>l || add_approx(o,mul_approx(s,b))<u ) {
            s=next_up(s);
            o=div_approx(sub_approx(mul_approx(l,b),mul_approx(u,a)),b-a);
          }
        }

        //std::cerr << "o=" << o << "; s=" << s << std::endl;

        // check  o+sa <= l <= u <= o+sb
        //std::cerr << "add_approx(o,mul_approx(s,a))=" << add_approx(o,mul_approx(s,a)) << std::endl;
        assert(add_approx(o,mul_approx(s,a)) <= l);
        //std::cerr << "add_approx(o,mul_approx(s,b))=" << add_approx(o,mul_approx(s,b)) << std::endl;
        assert(add_approx(o,mul_approx(s,b)) >= u);
      }

      this->create();
    }
  


     



    template<class R>
    void
    Grid<R>::create() 
    {
    }



    template<class R>
    dimension_type
    Grid<R>::dimension() const
    {
      return this->_lengths.size();
    }



    template<class R>
    R
    Grid<R>::subdivision_coordinate(dimension_type d, index_type n) const 
    {
        return add_approx(this->_origin[d],mul_approx(this->_lengths[d],n));
    }


    template<class R> 
    index_type 
    Grid<R>::subdivision_index(dimension_type d, const real_type& x) const 
    {
      using namespace Numeric;
      
      R half=0.5;
      index_type n=int_down<index_type>(add_approx(div_approx(sub_approx(x,this->_origin[d]),this->_lengths[d]),half));
      R sc=add_approx(this->_origin[d],mul_approx(this->_lengths[d],n));
      if(sc == x) { 
        return n; 
      } else {
        throw std::runtime_error("Value is not a grid coordinate");
      }
    }
    

    template<class R> 
    index_type 
    Grid<R>::subdivision_lower_index(dimension_type d, const real_type& x) const 
    {
      using namespace Numeric;
      
      index_type n=int_down<index_type>(div_down(sub_down(x,this->_origin[d]),this->_lengths[d]));
      if(x>=add_approx(this->_origin[d],mul_approx(this->_lengths[d],(n+1)))) {
        return n+1;
      } else {
        return n;
      }
    }


    template<class R> 
    index_type 
    Grid<R>::subdivision_upper_index(dimension_type d, const real_type& x) const 
    {
      using namespace Numeric;
      
      index_type n=int_up<index_type>(div_up(sub_up(x,this->_origin[d]),this->_lengths[d]));
      if(x<=add_approx(this->_origin[d],mul_approx(this->_lengths[d],(n-1)))) {
        return n-1;
      } else {
        return n;
      }
    }



    template<class R> 
    bool 
    Grid<R>::operator==(const Grid<R>& g) const
    {
      if(this==&g) { 
        return true; 
      } else {
        return this->_origin==g._origin && this->_lengths==g._lengths;
      }
    }
  

    template<class R> 
    bool 
    Grid<R>::operator!=(const Grid<R>& g) const
    {
      return !(*this==g);
    }

                   
    template<class R>
    IndexArray 
    Grid<R>::index(const Point<R>& s) const
    {
      IndexArray res(s.dimension());
      for(size_type i=0; i!=res.size(); ++i) {
        res[i]=subdivision_index(i,s[i]);
      }
      return res;
    }


    template<class R>
    IndexArray  
    Grid<R>::lower_index(const Rectangle<R>& r) const {
      IndexArray res(r.dimension());
      for(size_type i=0; i!=res.size(); ++i) {
        res[i]=subdivision_lower_index(i,r.lower_bound(i));
      }
      return res;
    }


    template<class R>
    IndexArray  
    Grid<R>::upper_index(const Rectangle<R>& r) const {
      IndexArray res(r.dimension());
      for(size_type i=0; i!=res.size(); ++i) {
        res[i]=subdivision_upper_index(i,r.upper_bound(i));
      }
      return res;
    }


    template<class R>
    Combinatoric::LatticeBlock  
    Grid<R>::index_block(const Rectangle<R>& r) const {
      Combinatoric::LatticeBlock res(r.dimension());
      for(size_type i=0; i!=res.dimension(); ++i) {
        res.set_lower_bound(i,this->subdivision_lower_index(i,r.lower_bound(i)));
        res.set_upper_bound(i,this->subdivision_upper_index(i,r.upper_bound(i)));
      }
      return res;
    }


    template<class R>
    Point<R> 
    Grid<R>::point(const IndexArray& a) const
    {
      Point<R> res(a.size());
      for(size_type i=0; i!=res.dimension(); ++i) {
        res[i]=subdivision_coordinate(i,a[i]);
      }
      return res;
    }


    template<class R>
    Rectangle<R> 
    Grid<R>::rectangle(const Combinatoric::LatticeBlock& lb) const
    {
      Rectangle<R> res(lb.dimension());
      for(size_type i=0; i!=res.dimension(); ++i) {
        res.set_lower_bound(i,this->subdivision_coordinate(i,lb.lower_bound(i)));
        res.set_upper_bound(i,this->subdivision_coordinate(i,lb.upper_bound(i)));
      }
      return res;
    }


    template<class R>
    std::ostream&
    Grid<R>::write(std::ostream& os) const
    {
      os << "Grid( origin=" << this->_origin
         << ", lengths=" << this->_lengths
         << " )";
      return os;
    }



    template<class R>
    std::istream&
    Grid<R>::read(std::istream& is) 
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }







    template<class R>
    FiniteGrid<R>::~FiniteGrid() 
    {
      if(this->_own_ptr) {
        delete this->_grid_ptr;
      }
    }


    template<class R> inline
    FiniteGrid<R>::FiniteGrid(const Grid<R>& g, const Combinatoric::LatticeBlock& b) 
      : _own_ptr(0), _grid_ptr(&g), _lattice_block(b)
    { 
    }
      
    template<class R> inline
    FiniteGrid<R>::FiniteGrid(const Grid<R>& g, const Rectangle<R>& bb) 
      : _own_ptr(0), _grid_ptr(&g), _lattice_block(over_approximation(bb,g).lattice_set())
    { 
    }
      
    template<class R>
    FiniteGrid<R>::FiniteGrid(const Rectangle<R>& bb, const size_type& s)
      : _own_ptr(1), _grid_ptr(0), _lattice_block(bb.dimension())
    {
      dimension_type d=bb.dimension();
      array<R> subdivision_lengths(bb.dimension());
      for(dimension_type i=0; i!=bb.dimension(); ++i) {
        subdivision_lengths[i]=div_up(bb[i].length(),R(s));
      }
      this->_grid_ptr=new Grid<R>(LinearAlgebra::Vector<R>(d,subdivision_lengths.begin()));
      this->_lattice_block=this->grid().index_block(bb);
    }


    template<class R>
    dimension_type 
    FiniteGrid<R>::dimension() const 
    {
      return this->_grid_ptr->dimension();
    }


    template<class R>
    Rectangle<R>
    FiniteGrid<R>::extent() const
    {
      return this->_grid_ptr->rectangle(this->_lattice_block);
    }

    
    template<class R>
    std::ostream&
    FiniteGrid<R>::write(std::ostream& os) const
    {
      return os << "FiniteGrid(grid=" << this->grid() << ", lattice_block=" << this->lattice_block() << ")";
    }


  }
}
